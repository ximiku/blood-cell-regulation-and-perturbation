"""Compare two selected B-cell GPU solutions with sklearn on the same objective."""
from pathlib import Path
import json,time
import numpy as np
from sklearn.linear_model import Lasso
from threadpoolctl import threadpool_limits
import lasso_analysis as analysis

ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.monotonic();idx=4
    X=np.load(ROOT/'.cache/lasso/X.npy',mmap_mode='r');y=np.load(ROOT/'.cache/lasso/Y.npy',mmap_mode='r')[:,idx]
    split=json.loads((ROOT/'data/splits/gene_cv_splits.json').read_text())
    cv=analysis.tune(X,y,idx,split['full_inner'],'full')
    selected=[cv['min_index'],cv['one_se_index']]
    fit=analysis.fit_context(X,y,idx,'full_fit',analysis.RATIOS)
    xx=np.asfortranarray((X-fit['x_mean'])/fit['x_scale']);yy=(y-fit['y_mean'])/fit['y_scale']
    gram=xx.T@xx
    records=[]
    for name,j in zip(['minimum_MSE','one_SE'],selected):
        alpha=float(analysis.RATIOS[j]*fit['alpha_max'])
        t=time.monotonic()
        model=Lasso(alpha=alpha,fit_intercept=False,precompute=gram,tol=1e-8,max_iter=50000).fit(xx,yy)
        b=fit['coef'][:,j]
        def objective(beta):return float(.5*np.mean((yy-xx@beta)**2)+alpha*np.abs(beta).sum())
        diff=xx@(b-model.coef_)
        record={'selection':name,'relative_alpha':float(analysis.RATIOS[j]),'alpha':alpha,
            'gpu_primal':objective(b),'sklearn_primal':objective(model.coef_),
            'sklearn_dual_gap':float(model.dual_gap_),'sklearn_n_iter':int(model.n_iter_),
            'prediction_RMS_difference_response_SD':float(np.sqrt(np.mean(diff**2))),
            'prediction_RMS_difference_input_units':float(np.sqrt(np.mean(diff**2))*fit['y_scale']),
            'max_coefficient_difference':float(np.max(np.abs(b-model.coef_))),
            'GPU_active_gt_1e_6':int(np.sum(np.abs(b)>1e-6)),'sklearn_active_gt_1e_6':int(np.sum(np.abs(model.coef_)>1e-6)),
            'reference_seconds':time.monotonic()-t}
        assert record['gpu_primal']-record['sklearn_primal']<=analysis.CONFIG['solver']['tol']*1.0001
        assert model.dual_gap_<1.0001e-8
        records.append(record);print(json.dumps(record),flush=True)
    output={'target':'B cell','reference':'sklearn Lasso, cold start, tol=1e-8','passed':True,'seconds':time.monotonic()-start,'records':records}
    (ROOT/'logs/03_solver_crosscheck_B_cell.json').write_text(json.dumps(output,indent=2)+'\n')

if __name__=='__main__':
    with threadpool_limits(limits=2):main()
