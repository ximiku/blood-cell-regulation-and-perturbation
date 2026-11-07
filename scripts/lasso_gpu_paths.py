"""Float64 batched Lasso paths with per-target primal/dual certificates on CUDA.

Only computation is batched: each column remains an independent Lasso objective.
Q=X_standardized.T@X_standardized/n and c=X_standardized.T@y_standardized/n
are constructed exclusively from that context's training genes.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse,json,time
import numpy as np
import cupy as cp

ROOT=Path(__file__).resolve().parents[1]
CACHE=ROOT/'.cache/lasso';OUT=ROOT/'results/lasso_gpu_paths'
CFG=json.loads((ROOT/'config/analysis_config.json').read_text())
RATIOS=np.geomspace(CFG['lasso_relative_alpha']['maximum'],CFG['lasso_relative_alpha']['minimum'],CFG['lasso_relative_alpha']['count'])


def certificates(Q,c,s,b,alpha):
    qb=Q@b
    cb=cp.sum(c*b,axis=0)
    residual_variance=cp.maximum(s-2*cb+cp.sum(b*qb,axis=0),0)
    primal=.5*residual_variance+alpha*cp.sum(cp.abs(b),axis=0)
    scale=cp.maximum(1,cp.max(cp.abs(c-qb),axis=0)/alpha)
    dual=(s-cb)/scale-.5*residual_variance/(scale*scale)
    return primal-dual,primal,dual


def solve_context(name,device):
    OUT.mkdir(parents=True,exist_ok=True)
    dest=OUT/f'{name}.npz'
    if dest.exists() and (OUT/f'{name}.json').exists():return
    start=time.monotonic()
    with np.load(CACHE/f'{name}_moments.npz') as z:m={k:z[k] for k in z.files}
    n=len(m['train'])
    cp.cuda.Device(device).use()
    Q=cp.asarray(np.load(CACHE/f'{name}_gram.npy')/n,dtype=cp.float64)
    Q=(Q+Q.T)*.5
    c=cp.asarray(m['xy']/n,dtype=cp.float64)
    Y=np.load(CACHE/'Y.npy',mmap_mode='r')[m['train']]
    s=cp.asarray(np.mean(((Y-m['y_mean'])/m['y_scale'])**2,axis=0))
    # Largest eigenvalue gives a valid smooth-loss step size, with a numerical margin.
    L=float(cp.linalg.eigvalsh(Q)[-1].get())*1.000001
    amax=cp.max(cp.abs(c),axis=0)
    tolerance=CFG['solver']['tol']*s
    b=cp.zeros(c.shape,dtype=cp.float64)
    paths=[];gaps=[];iterations=[];primals=[];duals=[];seconds=[]
    for j,ratio in enumerate(RATIOS):
        t0=time.monotonic();alpha=amax*ratio
        v=b.copy();momentum=cp.ones(c.shape[1],dtype=cp.float64)
        gap,primal,dual=certificates(Q,c,s,b,alpha)
        done=gap<=tolerance
        all_done=bool(cp.all(done).get())
        finished_at=cp.zeros(c.shape[1],dtype=cp.int64)
        for it in range(1,CFG['solver']['retry_max_iter']+1):
            if all_done:break
            gradient=Q@v-c
            step=v-gradient/L
            update=cp.sign(step)*cp.maximum(cp.abs(step)-alpha/L,0)
            update=cp.where(done[None,:],b,update)
            next_momentum=(1+cp.sqrt(1+4*momentum*momentum))/2
            restart=cp.sum((v-update)*(update-b),axis=0)>0
            next_momentum=cp.where(restart,1,next_momentum)
            v=update+((momentum-1)/next_momentum)[None,:]*(update-b)
            v=cp.where(restart[None,:],update,v)
            b=update;momentum=next_momentum
            if it==1 or it%25==0:
                gap,primal,dual=certificates(Q,c,s,b,alpha)
                newly=(~done)&(gap<=tolerance)
                finished_at=cp.where(newly,it,finished_at)
                done=done|newly
                all_done=bool(cp.all(done).get())
        cp.cuda.Stream.null.synchronize()
        gap,primal,dual=certificates(Q,c,s,b,alpha)
        if not bool(cp.all((gap<=tolerance*1.00001)&(gap>=-1e-8)).get()):
            raise RuntimeError(f'{name} ratio={ratio}: unconverged gaps {gap.get().tolist()}')
        paths.append(b.get());gaps.append(gap.get());iterations.append(finished_at.get())
        primals.append(primal.get());duals.append(dual.get());seconds.append(time.monotonic()-t0)
        print(json.dumps({'context':name,'device':device,'alpha_index':j,'ratio':ratio,'max_gap':float(gap.max().get()),
                          'max_iterations':int(finished_at.max().get()),'seconds':seconds[-1]}),flush=True)
    # feature x target x alpha, all paths use identical gene and perturbation order.
    temporary=dest.with_suffix('.npz.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream,coef=np.stack(paths,axis=2),gap=np.stack(gaps,axis=1),
            n_iter=np.stack(iterations,axis=1),primal=np.stack(primals,axis=1),dual=np.stack(duals,axis=1),
            alpha_max=amax.get(),ratios=RATIOS,dual_tolerance=tolerance.get(),Lipschitz=L,
            seconds_by_alpha=np.asarray(seconds),n_train=n)
    temporary.replace(dest)
    record={'context':name,'device':device,'solver':'float64_adaptive_restart_FISTA',
            'seconds':time.monotonic()-start,'cupy':cp.__version__,'numpy':np.__version__,
            'all_converged':True,'max_iterations':int(np.max(iterations)),
            'completed_utc':datetime.now(timezone.utc).isoformat()}
    (OUT/f'{name}.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'event':'context_done',**record}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--device',type=int,default=0);p.add_argument('--contexts',nargs='+',required=True)
    args=p.parse_args()
    for context in args.contexts:solve_context(context,args.device)
