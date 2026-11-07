"""Read-only validation of the exploratory, equal-size non-landmark controls."""
from pathlib import Path
import json,time
import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1]
def close(a,b):np.testing.assert_allclose(a,b,atol=2e-9,rtol=2e-8,equal_nan=True)
def main():
    start=time.monotonic();cfg=json.loads((ROOT/'config/stage4_config.json').read_text());focus=cfg['focus_targets']
    with h5py.File(ROOT/'data/processed/analysis_matrices.h5') as f:
        X=f['model/X'][:].astype(float);targets=f['model/target_cells'].asstr()[:].tolist()
        Y=f['model/Y'][:][:,[targets.index(t) for t in focus]];lm=set(np.flatnonzero(f['model/is_landmark'][:]).tolist())
    split=json.loads((ROOT/'data/splits/gene_cv_splits.json').read_text());samples=json.loads((ROOT/'data/splits/stage4_size_control_samples.json').read_text())
    jobs=json.loads((ROOT/'data/splits/stage4_size_control_jobs.json').read_text())['jobs'];assert len(jobs)==60
    expected={}
    for dataset,sample in samples.items():
        rng=np.random.default_rng(sample['seed']);chosen=[]
        for e in split['outer']:
            count=sum(i in lm for i in e['test']);pool=[i for i in e['test'] if i not in lm]
            chosen.extend(rng.choice(pool,count,replace=False).tolist())
        assert sorted(chosen)==sample['gene_indices'] and len(chosen)==965 and not set(chosen)&lm
        keep=set(chosen)
        for e in split['outer']:
            prefix=f'outer_{e["fold"]:02d}';expected[dataset+'/'+prefix+'_fit']=([i for i in e['train'] if i in keep],[])
            for v in e['inner']:expected[dataset+'/'+prefix+f'_inner_{v["fold"]:02d}']=([i for i in v['train'] if i in keep],[i for i in v['validation'] if i in keep])
    assert set(expected)=={j['id'] for j in jobs}
    mse_error=0.;gap_error=0.;direct=0;saved=0
    for job in jobs:
        assert (job['train'],job['validation'])==expected[job['id']]
        path=ROOT/'results/stage4'/f'{job["id"]}.npz'
        with np.load(path) as z:a={k:z[k] for k in z.files}
        tr=a['train'];va=a['validation'];assert np.array_equal(tr,job['train']) and np.array_equal(va,job['validation'])
        assert not set(tr)&set(va)
        xm=X[tr].mean(axis=0);xs=X[tr].std(axis=0);xs[xs==0]=1;ym=Y[tr].mean(axis=0);ys=Y[tr].std(axis=0)
        close(xm,a['x_mean']);close(xs,a['x_scale']);close(ym,a['y_mean']);close(ys,a['y_scale'])
        xx=(X[tr]-xm)/xs;yz=(Y[tr]-ym)/ys;coef=a['Lasso_coef'];d=a['Lasso_diagnostics']
        ratios=np.geomspace(1,.001,50);amax=np.max(np.abs(xx.T@yz),axis=0)/len(tr)
        close(d[:,0,:],ratios[:,None]*amax[None,:]);close(d[:,1,:],d[:,2,:]-d[:,3,:])
        assert np.all(d[:,1,:]>=-1e-8) and np.all(d[:,1,:]<=1.0001e-4);saved+=d[:,1,:].size
        if len(va):
            p=(((X[va]-xm)/xs)@coef.reshape(X.shape[1],-1)).reshape(len(va),4,50)*ys[None,:,None]+ym[None,:,None]
            mse=np.mean((Y[va,:,None]-p)**2,axis=0);mse_error=max(mse_error,float(np.abs(mse-a['Lasso_validation_MSE']).max()));close(mse,a['Lasso_validation_MSE'])
        for j in [0,25,49]:
            beta=coef[:,:,j];res=yz-xx@beta;alpha=d[j,0];rr=np.mean(res**2,axis=0)
            scale=np.maximum(1,np.max(np.abs(xx.T@res),axis=0)/(len(tr)*alpha))
            primal=.5*rr+alpha*np.abs(beta).sum(axis=0);dual=np.mean(yz*res,axis=0)/scale-.5*rr/scale**2
            gap=primal-dual;gap_error=max(gap_error,float(np.abs(gap-d[j,1]).max()));close(gap,d[j,1]);direct+=4
    curves=pd.read_csv(ROOT/'tables/04_size_control_tuning.csv',float_precision='round_trip');choices=pd.read_csv(ROOT/'tables/04_size_control_choices.csv',float_precision='round_trip')
    pred=pd.read_csv(ROOT/'tables/04_size_control_predictions.csv',float_precision='round_trip');metrics=pd.read_csv(ROOT/'tables/04_size_control_metrics.csv',float_precision='round_trip')
    pairs=0
    for keys,f in curves.groupby(['dataset','target','scope']):
        f=f.sort_values('parameter_index');foldmse=f[[f'fold_{i}_MSE' for i in range(1,4)]].to_numpy()
        means=foldmse.mean(axis=1);se=foldmse.std(axis=1,ddof=1)/np.sqrt(3);best=int(np.argmin(means));one=int(np.flatnonzero(means<=means[best]+se[best])[0])
        close(means,f.mean_MSE);close(se,f.SE_MSE)
        with np.load(ROOT/'results/stage4'/keys[0]/(keys[2]+'_fit.npz')) as z:
            for name,j in [('min',best),('1se',one)]:
                row=choices[(choices.dataset==keys[0])&(choices.target==keys[1])&(choices.scope==keys[2])&(choices.selection==name)]
                assert len(row)==1 and int(row.iloc[0].parameter_index)==j
                t=focus.index(keys[1]);v=pred[(pred.dataset==keys[0])&(pred.target==keys[1])&(pred.outer_fold==int(keys[2][-2:]))]
                pp=((X[v.gene_index.to_numpy()]-z['x_mean'])/z['x_scale'])@z['Lasso_coef'][:,t,j]*z['y_scale'][t]+z['y_mean'][t]
                close(pp,v['Lasso_'+name])
        pairs+=1
    for _,r in metrics.iterrows():
        f=pred[(pred.dataset==r.dataset)&(pred.target==r.target)];obs=f.observed.to_numpy();p=f['Mean_baseline' if r.model=='Mean_baseline' else 'Lasso_'+r.selection].to_numpy()
        calculated=[1-np.sum((p-obs)**2)/np.sum((obs-obs.mean())**2),np.sqrt(np.mean((p-obs)**2)),np.corrcoef(obs,p)[0,1] if np.ptp(p)>0 else np.nan,spearmanr(obs,p).statistic if np.ptp(p)>0 else np.nan]
        close(r[['R2','RMSE','Pearson','Spearman']].to_numpy(float),calculated)
    out={'passed':True,'seconds':time.monotonic()-start,'matched_subsets':3,'contexts':60,'saved_certificates':saved,
        'independent_direct_certificates':direct,'maximum_validation_MSE_error':mse_error,'maximum_gap_error':gap_error,
        'selection_pairs':pairs,'metric_rows':len(metrics),'sampling_and_outer_counts_reconstructed':True}
    (ROOT/'logs/stage4_size_control_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=4):main()
