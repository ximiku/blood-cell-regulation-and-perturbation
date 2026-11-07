"""Real-fold CPU references for all five newly implemented model families."""
from pathlib import Path
import json,time
import numpy as np
from sklearn.linear_model import ElasticNet,Lasso,Ridge
from sklearn.cross_decomposition import PLSRegression
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.monotonic();base=ROOT/'results/stage4/raw/outer_01_inner_01'
    with np.load(base.with_suffix('.npz')) as z:a={k:z[k] for k in z.files}
    meta=json.loads(base.with_suffix('.json').read_text());tr=a['train'];va=a['validation']
    X=np.load(ROOT/'.cache/stage4/X.npy',mmap_mode='r');Y=np.load(ROOT/'.cache/stage4/Y_raw.npy',mmap_mode='r')
    x=(X[tr]-X[tr].mean(axis=0))/X[tr].std(axis=0);y=(Y[tr,0]-Y[tr,0].mean())/Y[tr,0].std()
    np.testing.assert_allclose(a['x_mean'],X[tr].mean(axis=0),atol=2e-13)
    np.testing.assert_allclose(a['x_scale'],X[tr].std(axis=0),atol=2e-13)
    records=[]
    for model,index in [('ElasticNet',25),('ElasticNet',75),('ElasticNet',125),('NonnegativeLasso',25)]:
        par=meta['parameters'][model][index];mix=par['l1_ratio'];alpha=a[model+'_diagnostics'][index,0,0]
        ref=(Lasso(alpha=alpha,fit_intercept=False,positive=True,max_iter=100000,tol=1e-8)
             if model=='NonnegativeLasso' else ElasticNet(alpha=alpha,l1_ratio=mix,fit_intercept=False,max_iter=100000,tol=1e-8))
        t=time.monotonic();ref.fit(x,y);b=a[model+'_coef'][:,0,index]
        def primal(beta):return .5*np.mean((y-x@beta)**2)+alpha*mix*np.abs(beta).sum()+.5*alpha*(1-mix)*np.sum(beta**2)
        item={'model':model,'parameter_index':index,**par,'alpha':float(alpha),
          'GPU_primal':float(primal(b)),'reference_primal':float(primal(ref.coef_)),
          'reference_dual_gap':float(ref.dual_gap_),'reference_iterations':int(ref.n_iter_),
          'prediction_RMS_difference_response_SD':float(np.sqrt(np.mean((x@(b-ref.coef_))**2))),
          'max_coefficient_difference':float(np.max(np.abs(b-ref.coef_))),'seconds':time.monotonic()-t}
        assert abs(item['GPU_primal']-item['reference_primal'])<=1.0001e-4
        assert item['reference_dual_gap']<1.1e-8
        if model=='NonnegativeLasso':assert b.min()>=0
        records.append(item);print(json.dumps(item),flush=True)
    for j in [0,8,16]:
        lam=meta['parameters']['Ridge'][j]['lambda']
        ref=Ridge(alpha=len(tr)*lam,fit_intercept=False,solver='cholesky').fit(x,y)
        b=a['Ridge_coef'][:,0,j];error=float(np.max(np.abs(x@(b-ref.coef_))))
        assert error<1e-7
        records.append({'model':'Ridge','lambda':lam,'max_prediction_difference_response_SD':error})
    # An independent CPU symmetric eigensolver reconstructs PCR from raw training rows.
    vals,vec=eigh(x.T@x/len(tr),subset_by_index=[x.shape[1]-100,x.shape[1]-1])
    c=x.T@y/len(tr)
    for j,par in enumerate(meta['parameters']['PCR']):
        k=par['components'];v=vec[:,-k:];ref=v@((v.T@c)/vals[-k:])
        error=float(np.max(np.abs(x@(a['PCR_coef'][:,0,j]-ref))))
        assert error<1e-8
        records.append({'model':'PCR','components':k,'max_prediction_difference_response_SD':error})
    for j,par in enumerate(meta['parameters']['PLS']):
        ref=PLSRegression(n_components=par['components'],scale=False).fit(x,y)
        error=float(np.max(np.abs(ref.predict(x).ravel()-x@a['PLS_coef'][:,0,j])))
        assert error<1e-8
        records.append({'model':'PLS','components':par['components'],'max_prediction_difference_response_SD':error})
    out={'passed':True,'context':'raw/outer_01_inner_01','target':'B cell','seconds':time.monotonic()-start,'records':records}
    (ROOT/'logs/stage4_solver_crosscheck.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'n_reference_comparisons':len(records),'seconds':out['seconds']}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=4):main()
