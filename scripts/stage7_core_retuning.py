"""Prespecified 20 gene resamples with fresh 3-fold tuning for all four core targets."""
from pathlib import Path
import argparse,json,time
import numpy as np,cupy as cp
from stage4_gpu import proximal,PATH,CFG,save_npz
R=Path(__file__).resolve().parents[1];O=R/'results/stage7/core_retuning';O.mkdir(parents=True,exist_ok=True)
def fit(X,Y,tr,va,ratios=PATH):
 x=X[tr];xm=x.mean(0);xs=x.std(0);xs=cp.where(xs>0,xs,1);x=(x-xm)/xs;y=Y[tr];ym=y.mean(0);ys=y.std(0);y=(y-ym)/ys
 Q=x.T@x/len(tr);Q=(Q+Q.T)/2;c=x.T@y/len(tr);s=cp.mean(y*y,axis=0);L=float(cp.linalg.eigvalsh(Q)[-1].get())*1.000001
 b,diag,amax=proximal(Q,c,s,L,ratios)
 mse=np.empty((4,0))
 if len(va):
  pred=cp.einsum('np,ptk->ntk',(X[va]-xm)/xs,b)*ys[None,:,None]+ym[None,:,None];mse=cp.mean((Y[va,:,None]-pred)**2,axis=0).get()
 return {'coef':b.get(),'diagnostics':diag,'alpha_max':amax,'x_mean':xm.get(),'x_scale':xs.get(),'y_mean':ym.get(),'y_scale':ys.get(),'MSE':mse}
def main(device,shard,shards):
 cp.cuda.Device(device).use();X=cp.asarray(np.load(R/'.cache/stage4/X.npy'));Y=cp.asarray(np.load(R/'.cache/stage4/Y_raw.npy'));samples=np.load(R/'data/splits/stage4_stability_indices.npz')['train_indices'][:20]
 for rep in range(shard,20,shards):
  out=O/f'repeat_{rep+1:02d}.npz';record=out.with_suffix('.json')
  if out.exists() and record.exists():continue
  t=time.time();tr=samples[rep];rng=np.random.default_rng(20268107+rep);folds=np.array_split(rng.permutation(tr),3);arrays={'train_indices':tr};errs=[]
  for fold,va in enumerate(folds):
   inner=np.setdiff1d(tr,va);result=fit(X,Y,inner,va);errs.append(result['MSE']);arrays[f'inner{fold+1}_train']=inner;arrays[f'inner{fold+1}_validation']=va
   arrays[f'inner{fold+1}_coef_path']=result['coef'];arrays[f'inner{fold+1}_diagnostics']=result['diagnostics'];arrays[f'inner{fold+1}_alpha_max']=result['alpha_max']
  err=np.stack(errs);mean=err.mean(0);se=err.std(0,ddof=1)/np.sqrt(3);minimum=np.argmin(mean,axis=1);simple=np.array([np.flatnonzero(mean[i]<=mean[i,minimum[i]]+se[i,minimum[i]])[0] for i in range(4)]);selected=np.stack([minimum,simple],axis=1)
  # Two independent full-subsample fits from zero for min and 1-SE parameter vectors.
  full_min=fit(X,Y,tr,[],[PATH[minimum]]);full_simple=fit(X,Y,tr,[],[PATH[simple]])
  arrays.update(cv_MSE=err,relative_path=PATH,selected_indices=selected,coef_min=full_min['coef'][:,:,0],coef_1se=full_simple['coef'][:,:,0],diagnostics_min=full_min['diagnostics'],diagnostics_1se=full_simple['diagnostics'],alpha_max=full_min['alpha_max'],x_mean=full_min['x_mean'],x_scale=full_min['x_scale'],y_mean=full_min['y_mean'],y_scale=full_min['y_scale'])
  save_npz(out,**arrays);meta={'repeat':rep+1,'seed':20268107+rep,'targets':4,'inner_folds':3,'all_converged':True,'initialization':'zero per inner path and per full selected fit','seconds':time.time()-t,'device':device};record.write_text(json.dumps(meta,indent=2)+'\n');print(meta,flush=True);cp.get_default_memory_pool().free_all_blocks()
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--device',type=int,required=True);p.add_argument('--shard',type=int,default=0);p.add_argument('--shards',type=int,default=1);a=p.parse_args();main(a.device,a.shard,a.shards)
