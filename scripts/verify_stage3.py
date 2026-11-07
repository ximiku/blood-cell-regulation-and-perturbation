"""Independently audit splits, original-scale predictions, selection and artifacts."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json
import xml.etree.ElementTree as ET
import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from PIL import Image
from PyPDF2 import PdfReader
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
CHECKS=[]


def passed(name,detail):CHECKS.append({'check':name,'passed':True,'detail':detail})
def same(a,b,atol=1e-9):np.testing.assert_allclose(a,b,atol=atol,rtol=1e-9,equal_nan=True)
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()


def direct_dual_gap(x_standard,y_standard,betas,alphas):
    """Evaluate certificates from original training rows, without a Gram cache."""
    residual=y_standard-x_standard@betas
    rr=np.mean(residual**2,axis=0)
    corr=x_standard.T@residual/len(residual)
    factor=np.maximum(1,np.max(np.abs(corr),axis=0)/alphas)
    primal=.5*rr+alphas*np.sum(np.abs(betas),axis=0)
    dual=np.mean(y_standard*residual,axis=0)/factor-.5*rr/factor**2
    return primal-dual


def main():
    config=json.loads((ROOT/'config/analysis_config.json').read_text())
    old=json.loads((ROOT/'logs/stage1_artifact_manifest.json').read_text())
    h5=ROOT/'data/processed/analysis_matrices.h5'
    h5_hash=digest(h5)
    assert h5_hash==next(x['sha256'] for x in old['files'] if x['path']=='data/processed/analysis_matrices.h5')
    with h5py.File(h5) as f:
        X=f['model/X'][:].astype(float);Y=f['model/Y'][:]
        genes=f['model/genes'].asstr()[:];targets=f['model/target_cells'].asstr()[:];pert=f['model/perturbations'].asstr()[:]
    split=json.loads((ROOT/'data/splits/gene_cv_splits.json').read_text())
    assert split['gene_order_sha256']==hashlib.sha256('\n'.join(genes).encode()).hexdigest()
    assert len(targets)==21 and len(genes)==11922 and X.shape[1]==2717
    allgenes=set(range(len(genes)));count=np.zeros(len(genes),int)
    for entry in split['outer']:
        tr=set(entry['train']);te=set(entry['test']);assert tr.isdisjoint(te) and tr|te==allgenes
        count[entry['test']]+=1
        inner_count={g:0 for g in tr}
        for inner in entry['inner']:
            itr=set(inner['train']);iva=set(inner['validation'])
            assert itr.isdisjoint(iva) and itr|iva==tr and (itr|iva).isdisjoint(te)
            for g in iva:inner_count[g]+=1
        assert set(inner_count.values())=={1}
    assert np.all(count==1)
    full_count=np.zeros(len(genes),int)
    for entry in split['full_inner']:
        tr=set(entry['train']);va=set(entry['validation']);assert tr.isdisjoint(va) and tr|va==allgenes
        full_count[entry['validation']]+=1
    assert np.all(full_count==1)
    passed('frozen_input_and_split_isolation','Frozen HDF5 unchanged; five disjoint outer tests and all inner train/validation sets form the intended nested partitions.')

    metric_table=pd.read_csv(ROOT/'tables/03_model_metrics.csv')
    selection=pd.read_csv(ROOT/'tables/03_selected_models.csv')
    all_coef=pd.read_csv(ROOT/'tables/03_all_coefficients_annotated.csv')
    model_files=[];pred_files=[]
    for j,t in enumerate(targets):
        directory=ROOT/f'results/lasso/target_{j:02d}'
        meta=json.loads((directory/'completed.json').read_text());assert meta['input_sha256']==h5_hash
        assert meta['target']==t and meta['all_converged']
        z=dict(np.load(directory/'model_arrays.npz'));model_files.append(z)
        pred=pd.read_csv(directory/'oof_predictions.csv');pred_files.append(pred)
        assert pred.gene_symbol.tolist()==genes.tolist() and z['perturbations'].tolist()==pert.tolist()
        same(pred.observed,Y[:,j])
        for col,model in enumerate(['Lasso_min','Lasso_1se','Mean_baseline']):
            yp=pred[model].to_numpy();yy=Y[:,j]
            expected=[1-np.sum((yy-yp)**2)/np.sum((yy-yy.mean())**2),np.sqrt(np.mean((yy-yp)**2)),np.corrcoef(yy,yp)[0,1],spearmanr(yy,yp).statistic]
            saved=metric_table[(metric_table.target==t)&(metric_table.model==model)].iloc[0]
            same(saved[['R2','RMSE','Pearson','Spearman']].to_numpy(dtype=float),expected)
        curves=pd.read_csv(directory/'tuning_curves.csv')
        assert len(curves)==6*50
        for scope,c in curves.groupby('scope',sort=False):
            c=c.reset_index(drop=True)
            err=c[[f'inner_{i}_MSE' for i in range(1,4)]].to_numpy()
            mean=err.mean(axis=1);se=err.std(axis=1,ddof=1)/np.sqrt(3)
            same(mean,c.mean_MSE);same(se,c.SE_MSE)
            imin=int(np.argmin(mean));i1=int(np.flatnonzero(mean<=mean[imin]+se[imin])[0])
            assert c.index[c.selected_min].tolist()==[imin] and c.index[c.selected_1se].tolist()==[i1]
            selected=selection[(selection.target==t)&(selection.scope==scope)].set_index('model')
            same(selected.loc['Lasso_min','relative_alpha'],c.relative_alpha.iloc[imin])
            same(selected.loc['Lasso_1se','relative_alpha'],c.relative_alpha.iloc[i1])
            assert c.relative_alpha.iloc[i1]>=c.relative_alpha.iloc[imin]
        for col,model in enumerate(['Lasso_min','Lasso_1se']):
            c=all_coef[(all_coef.target==t)&(all_coef.model==model)].set_index('perturbation').loc[pert]
            same(c.coefficient_standardized,z['full_coef_standardized'][:,col])
            same(c.coefficient_input,z['full_coef_input'][:,col])
            assert np.array_equal(c.active_gt_1e_6.to_numpy(),np.abs(z['full_coef_standardized'][:,col])>1e-6)
    passed('pooled_metrics_and_tuning','All 63 pooled metric rows and 126 scope-wise min/1-SE choices match independent recomputation from saved predictions and individual inner-fold errors.')
    passed('candidate_mapping','Every saved coefficient retains the frozen 2,717-perturbation order and exact active-threshold classification.')

    # Reconstruct every inner validation loss directly from the original rows and
    # the saved GPU coefficients; do not reuse cached means or validation losses.
    contexts=[(f"outer_{entry['fold']:02d}",inner) for entry in split['outer'] for inner in entry['inner']]
    contexts += [('full',inner) for inner in split['full_inner']]
    all_curves=[pd.read_csv(ROOT/f'results/lasso/target_{j:02d}/tuning_curves.csv') for j in range(len(targets))]
    max_inner_error=0.
    for scope,entry in contexts:
        tr=np.array(entry['train']);va=np.array(entry['validation'])
        xm=X[tr].mean(axis=0);xs=X[tr].std(axis=0);xs[xs==0]=1
        ym=Y[tr].mean(axis=0);ys=Y[tr].std(axis=0)
        name=f"{scope}_inner_{entry['fold']:02d}"
        with np.load(ROOT/'results/lasso_gpu_paths'/f'{name}.npz') as z:coef=z['coef']
        pred=(((X[va]-xm)/xs)@coef.reshape(X.shape[1],-1)).reshape(len(va),len(targets),50)
        pred=pred*ys[None,:,None]+ym[None,:,None]
        mse=np.mean((Y[va,:,None]-pred)**2,axis=0)
        for j,c in enumerate(all_curves):
            saved=c.loc[c.scope==scope,f"inner_{entry['fold']}_MSE"].to_numpy()
            same(mse[j],saved);max_inner_error=max(max_inner_error,float(np.max(np.abs(mse[j]-saved))))
    passed('inner_losses_from_original_validation_rows',f'All 18 inner folds × 21 targets × 50 penalties reproduce from independently recomputed training moments and original validation rows; maximum MSE discrepancy {max_inner_error:.3g}.')

    max_prediction_error=0.;certificate_maxima=[]
    for k,entry in enumerate(split['outer']):
        tr=np.array(entry['train']);te=np.array(entry['test'])
        xt=X[tr];xm=xt.mean(axis=0);xs=xt.std(axis=0);xs[xs==0]=1
        ym=Y[tr].mean(axis=0);ys=Y[tr].std(axis=0)
        xz=(xt-xm)/xs
        alpha_max=np.max(np.abs(xz.T@((Y[tr]-ym)/ys)),axis=0)/len(tr)
        for j,t in enumerate(targets):
            z=model_files[j];p=pred_files[j];b=z['outer_coef_input'][k];intercept=z['outer_intercept_input'][k]
            same(b,z['outer_coef_standardized'][k]*ys[j]/xs[:,None])
            same(intercept,ym[j]-xm@b)
            actual=X[te]@b+intercept
            saved=p.loc[te,['Lasso_min','Lasso_1se']].to_numpy()
            same(actual,saved);max_prediction_error=max(max_prediction_error,float(np.max(np.abs(actual-saved))))
            same(p.loc[te,'Mean_baseline'],np.full(len(te),ym[j]))
            assert (p.loc[te,'outer_fold']==k+1).all()
            rows=selection[(selection.target==t)&(selection.scope==f'outer_{k+1:02d}')]
            same(rows.alpha,rows.relative_alpha*alpha_max[j])
        betas=np.concatenate([z['outer_coef_standardized'][k] for z in model_files],axis=1)
        penalties=np.array([selection[(selection.target==t)&(selection.scope==f'outer_{k+1:02d}')].set_index('model').loc[model,'alpha']
                    for t in targets for model in ['Lasso_min','Lasso_1se']])
        gaps=direct_dual_gap(xz,np.repeat((Y[tr]-ym)/ys,2,axis=1),betas,penalties)
        assert np.all(gaps>=-1e-8) and np.all(gaps<=config['solver']['tol']*1.0001)
        certificate_maxima.append(float(gaps.max()))
    passed('outer_models_and_predictions',f'All 105 outer folds reproduce both OOF Lasso predictions, training-fold means/scales, intercepts, baseline predictions and alpha_max; maximum saved-prediction error {max_prediction_error:.3g}.')

    xm=X.mean(axis=0);xs=X.std(axis=0);xs[xs==0]=1;ym=Y.mean(axis=0);ys=Y.std(axis=0)
    alpha_max=np.max(np.abs(((X-xm)/xs).T@((Y-ym)/ys)),axis=0)/len(X)
    for j,z in enumerate(model_files):
        same(z['full_x_mean'],xm);same(z['full_x_scale'],xs);same(z['full_y_mean'],ym[j]);same(z['full_y_scale'],ys[j]);same(z['full_alpha_max'],alpha_max[j])
        same(z['full_coef_input'],z['full_coef_standardized']*ys[j]/xs[:,None])
        same(z['full_intercept_input'],ym[j]-xm@z['full_coef_input'])
        same(z['full_path_standardized'][:,z['selected_indices']],z['full_coef_standardized'])
    betas=np.concatenate([z['full_coef_standardized'] for z in model_files],axis=1)
    penalties=np.array([selection[(selection.target==t)&(selection.scope=='full')].set_index('model').loc[model,'alpha']
                for t in targets for model in ['Lasso_min','Lasso_1se']])
    gaps=direct_dual_gap((X-xm)/xs,np.repeat((Y-ym)/ys,2,axis=1),betas,penalties)
    assert np.all(gaps>=-1e-8) and np.all(gaps<=config['solver']['tol']*1.0001)
    certificate_maxima.append(float(gaps.max()))
    passed('full_data_models','All 21 full-data preprocessing records, selected paths, raw coefficients and intercepts agree with the frozen original input.')
    passed('independent_original_data_dual_certificates',f'All 252 selected outer/full models pass a CPU primal-dual calculation directly from training X/y, independently of cached GPU Gram matrices; maximum gap {max(certificate_maxima):.9g}.')

    diagnostics=pd.read_csv(ROOT/'tables/03_solver_diagnostics.csv')
    assert diagnostics.converged.all()
    assert np.isfinite(diagnostics.dual_gap).all()
    assert (diagnostics.dual_gap<=diagnostics.dual_gap_tolerance*1.00002).all()
    assert (diagnostics.dual_gap>=-1e-8).all()
    passed('convergence',f'{len(diagnostics)} recorded fits satisfy the saved dual-gap tolerance; {int(diagnostics.retried.sum())} retried fits.')
    corr=pd.read_csv(ROOT/'tables/03_single_perturbation_correlations.csv')
    assert len(corr)==21*2717
    for j in [0,4,5,13,20]:
        c=corr[corr.target==targets[j]].set_index('perturbation')
        for p in [0,77,500,1500,2716]:
            same(c.loc[pert[p],'Pearson'],np.corrcoef(X[:,p],Y[:,j])[0,1])
            same(c.loc[pert[p],'Spearman'],spearmanr(X[:,p],Y[:,j]).statistic)
    passed('single_perturbation_ranking','All 57,057 target/perturbation pairs are saved; 25 pairs spanning targets and feature indices agree with scipy/NumPy correlation calculations.')

    figures=json.loads((ROOT/'figures/03_figure_manifest.json').read_text())
    for rec in figures:
        for ext in ['pdf','svg','png']:
            p=ROOT/'figures'/f"{rec['stem']}.{ext}"
            if ext=='pdf':
                pdf=PdfReader(str(p));assert len(pdf.pages)==1 and len(pdf.pages[0].extract_text())>50
            elif ext=='svg':ET.parse(p)
            else:
                with Image.open(p) as im:assert abs(im.info['dpi'][0]-450)<.1;im.verify()
    passed('figure_files',f'{len(figures)} figure groups parse successfully as PDF, SVG and 450-dpi PNG.')
    result={'passed':True,'generated_utc':datetime.now(timezone.utc).isoformat(),'checks':CHECKS}
    (ROOT/'logs/stage3_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
