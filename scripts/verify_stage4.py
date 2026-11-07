"""Independent read-back checks against original rows, splits and saved predictions."""
from pathlib import Path
import argparse,hashlib,json,time
import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1];CHECKS=[]
def record(name,detail):
    CHECKS.append({'check':name,'passed':True,'detail':detail})
    print(json.dumps(CHECKS[-1],ensure_ascii=False),flush=True)
def close(a,b,atol=2e-9):np.testing.assert_allclose(a,b,atol=atol,rtol=2e-8,equal_nan=True)

def numeric():
    cfg=json.loads((ROOT/'config/stage4_config.json').read_text());focus=cfg['focus_targets'];thr=cfg['active_threshold']
    h5=ROOT/'data/processed/analysis_matrices.h5'
    assert hashlib.sha256(h5.read_bytes()).hexdigest()==cfg['input_sha256']
    with h5py.File(h5) as h:
        X=h['model/X'][:].astype(float);t=h['model/target_cells'].asstr()[:].tolist()
        Y=h['model/Y'][:][:,[t.index(x) for x in focus]]
        genes=h['model/genes'].asstr()[:];pert=h['model/perturbations'].asstr()[:]
        lm=set(np.flatnonzero(h['model/is_landmark'][:]).tolist())
    with np.load(ROOT/'data/processed/stage4_input_metadata.npz') as z:scale=z['gene_response_scale']
    split=json.loads((ROOT/'data/splits/gene_cv_splits.json').read_text());contexts={}
    for e in split['outer']:
        contexts[f'outer_{e["fold"]:02d}_fit']=(e['train'],[])
        for s in e['inner']:contexts[f'outer_{e["fold"]:02d}_inner_{s["fold"]:02d}']=(s['train'],s['validation'])
    contexts['full_fit']=(list(range(len(X))),[])
    for s in split['full_inner']:contexts[f'full_inner_{s["fold"]:02d}']=(s['train'],s['validation'])
    jobs=json.loads((ROOT/'data/splits/stage4_jobs.json').read_text())['jobs']
    assert len(jobs)==172 and len({j['id'] for j in jobs})==172
    for dataset in ['raw','landmark','gene_scaled']:
        assert {j['name'] for j in jobs if j['kind']=='context' and j['dataset']==dataset}==set(contexts)
    assert sum(j['kind']=='stability' for j in jobs)==100
    rng=np.random.default_rng(cfg['stability']['seed']);samples=np.stack([np.sort(rng.choice(len(X),int(.8*len(X)),replace=False)) for _ in range(100)])
    with np.load(ROOT/'data/splits/stage4_stability_indices.npz') as z:assert np.array_equal(samples,z['train_indices'])
    with np.load(ROOT/'data/processed/stage4_stability_reference.npz') as z:fixed=z['relative_alphas'];full=z['full_coefficients']
    for job in jobs:
        if job['kind']=='context':
            expected=contexts[job['name']]
            if job['dataset']=='landmark':expected=tuple([g for g in v if g in lm] for v in expected)
            assert job['train']==expected[0] and job['validation']==expected[1]
        else:assert np.array_equal(job['train'],samples[int(job['name'][-3:])-1])
        assert not set(job['train'])&set(job['validation'])
    record('frozen_input_shared_splits_and_100_draws',{'jobs':len(jobs),'landmark':len(lm),'sample_size':samples.shape[1]})
    own=json.loads((ROOT/'data/splits/stage4_self_readout_jobs.json').read_text())['jobs']
    error_mse=0.;error_gap=0.;n_direct=0;n_certificates=0;max_iter=0
    for job in jobs+own:
        path=ROOT/'results/stage4'/f'{job["id"]}.npz'
        assert path.is_file() and path.with_suffix('.json').is_file()
        with np.load(path) as z:a={k:z[k] for k in z.files}
        tr=np.array(job['train']);va=np.array(job['validation'],dtype=int)
        assert np.array_equal(a['train'],tr) and np.array_equal(a['validation'],va)
        yy=Y/scale[:,None] if job['dataset']=='gene_scaled' else Y
        xm=X[tr].mean(axis=0);xs=X[tr].std(axis=0);xs[xs==0]=1
        ym=yy[tr].mean(axis=0);ys=yy[tr].std(axis=0)
        close(a['x_mean'],xm);close(a['x_scale'],xs);close(a['y_mean'],ym);close(a['y_scale'],ys)
        xx=(X[tr]-xm)/xs;yz=(yy[tr]-ym)/ys;cross=xx.T@yz/len(tr)
        meta=json.loads(path.with_suffix('.json').read_text())
        if job['kind'] in ['stability','self_readout']:assert meta['initialization']=='zero'
        for model,pars in meta['parameters'].items():
            b=a[model+'_coef'];assert np.isfinite(b).all()
            if model=='NonnegativeLasso':assert b.min()>=0
            if len(va):
                pred=(((X[va]-xm)/xs)@b.reshape(X.shape[1],-1)).reshape(len(va),4,-1)*ys[None,:,None]+ym[None,:,None]
                mse=np.mean((yy[va,:,None]-pred)**2,axis=0)
                error_mse=max(error_mse,float(np.max(np.abs(mse-a[model+'_validation_MSE']))));close(mse,a[model+'_validation_MSE'])
        for key in [k for k in a if k.endswith('diagnostics')]:
            d=a[key];gap=d[:,1,:];tol=job.get('tolerance',cfg['tol'])
            assert np.isfinite(d).all() and np.all(gap>=-1e-8) and np.all(gap<=tol*1.0001)
            close(gap,d[:,2,:]-d[:,3,:],atol=1e-12);n_certificates+=gap.size;max_iter=max(max_iter,int(d[:,4,:].max()))
            model=key.removesuffix('_diagnostics')
            if key=='diagnostics':
                b=a['coef'][:,:,None];indices=[0];mixes=[1.];positive=False
                close(a['relative_alphas'],fixed)
                close(d[0,0],np.max(np.abs(cross),axis=0)*fixed)
            else:
                b=a[model+'_coef'];indices=([0,25,49,50,75,99,100,125,149] if model=='ElasticNet' else [0,25,49])
                mixes=[meta['parameters'][model][i]['l1_ratio'] for i in indices];positive=model=='NonnegativeLasso'
                for j,par in enumerate(meta['parameters'][model]):
                    alpha_max=np.max(cross if positive else np.abs(cross),axis=0)/par['l1_ratio']
                    close(d[j,0],alpha_max*par['relative_alpha'])
            chosen=b[:,:,indices];pred=(xx@chosen.reshape(X.shape[1],-1)).reshape(len(tr),4,-1)
            for j,(index,mix) in enumerate(zip(indices,mixes)):
                beta=chosen[:,:,j];alpha=d[index,0];l1=alpha*mix;l2=alpha*(1-mix)
                residual=yz-pred[:,:,j];aug=np.mean(residual**2,axis=0)+l2*np.sum(beta**2,axis=0)
                corr=xx.T@residual/len(tr)-l2*beta
                dualscale=np.maximum(1,np.max(corr if positive else np.abs(corr),axis=0)/l1)
                primal=.5*aug+l1*np.abs(beta).sum(axis=0)
                dual=np.mean(yz*residual,axis=0)/dualscale-.5*aug/dualscale**2
                direct=primal-dual;error_gap=max(error_gap,float(np.max(np.abs(direct-d[index,1]))))
                close(direct,d[index,1]);assert np.all(direct<=tol*1.0001)
                n_direct+=4
    record('fold_preprocessing_all_validation_MSE_and_certificates',{'n_saved_certificates':n_certificates,
        'n_independent_direct_certificates':n_direct,'max_MSE_error':error_mse,'max_gap_error':error_gap,'max_iterations':max_iter})
    # Recompute every choice independently from the three stored fold errors.
    curves=pd.read_csv(ROOT/'tables/04_tuning_curves.csv');selections=pd.read_csv(ROOT/'tables/04_selected_models.csv')
    count=0
    for keys,f in curves.groupby(['dataset','target','model','scope'],sort=False):
        f=f.sort_values('parameter_index');mse=f[[f'fold_{i}_MSE' for i in range(1,4)]].to_numpy()
        means=mse.mean(axis=1);se=mse.std(axis=1,ddof=1)/np.sqrt(3);best=int(np.argmin(means));limit=means[best]+se[best]
        eligible=np.flatnonzero(means<=limit)
        if f.relative_alpha.notna().any():one=sorted(eligible,key=lambda i:(-f.iloc[i].relative_alpha,-f.iloc[i].l1_ratio))[0]
        elif f['lambda'].notna().any():one=eligible[np.argmax(f.iloc[eligible]['lambda'])]
        else:one=eligible[np.argmin(f.iloc[eligible].components)]
        close(f.mean_MSE,means);close(f.SE_MSE,se)
        for which,i in [('min',best),('1se',one)]:
            row=selections[(selections.dataset==keys[0])&(selections.target==keys[1])&(selections.model==keys[2])&(selections.scope==keys[3])&(selections.selection==which)]
            assert len(row)==1 and int(row.iloc[0].parameter_index)==int(f.iloc[i].parameter_index)
        count+=1
    record('all_minimum_and_one_SE_choices',{'pairs':count})
    metrics_rows=0
    for dataset in ['raw','landmark','gene_scaled']:
        frame=pd.read_csv(ROOT/'tables'/f'04_{dataset}_oof_predictions.csv',float_precision='round_trip')
        m=pd.read_csv(ROOT/'tables'/('04_model_metrics.csv' if dataset=='raw' else '04_sensitivity_metrics.csv'));m=m[m.dataset==dataset]
        for _,row in m.iterrows():
            f=frame[frame.target==row.target];obs=f.observed.to_numpy()
            col='Mean_baseline' if row.model=='Mean_baseline' else row.model+'_'+row.selection;p=f[col].to_numpy()
            calc=[1-np.sum((p-obs)**2)/np.sum((obs-obs.mean())**2),np.sqrt(np.mean((p-obs)**2)),np.corrcoef(obs,p)[0,1] if np.ptp(p)>0 else np.nan,spearmanr(obs,p).statistic if np.ptp(p)>0 else np.nan]
            close(row[['R2','RMSE','Pearson','Spearman']].to_numpy(float),calc);metrics_rows+=1
        for path in (ROOT/'results/stage4/final').glob(dataset+'_*.npz'):
            model=path.stem.removeprefix(dataset+'_')
            with np.load(path) as z:
                ids=z['gene_indices'];pred=z['oof_predictions'];fold=z['outer_fold']
                for k in range(5):
                    mask=fold==k+1;pp=(X[ids[mask]]@z['outer_coef_input'][k].reshape(X.shape[1],-1)).reshape(mask.sum(),4,2)+z['outer_intercept_input'][k]
                    close(pp,pred[mask])
    record('all_pooled_metrics_and_original_row_outer_predictions',{'metric_rows':metrics_rows})
    fold_table=pd.read_csv(ROOT/'tables/04_fold_metrics.csv')
    cached={d:pd.read_csv(ROOT/'tables'/f'04_{d}_oof_predictions.csv',float_precision='round_trip') for d in ['raw','landmark','gene_scaled']}
    for _,row in fold_table.iterrows():
        f=cached[row.dataset];f=f[(f.target==row.target)&(f.outer_fold==row.fold)]
        obs=f.observed.to_numpy();p=f[row.model+'_'+row.selection].to_numpy()
        calc=[1-np.sum((p-obs)**2)/np.sum((obs-obs.mean())**2),np.sqrt(np.mean((p-obs)**2)),np.corrcoef(obs,p)[0,1] if np.ptp(p)>0 else np.nan,spearmanr(obs,p).statistic if np.ptp(p)>0 else np.nan]
        close(row[['R2','RMSE','Pearson','Spearman']].to_numpy(float),calc)
    record('every_outer_fold_metric',{'rows':len(fold_table)})
    with np.load(ROOT/'results/stage4/final/stability_coefficients.npz') as z:b=z['coefficients']
    active=np.abs(b)>thr;b0=np.where(active,b,0);freq=active.mean(axis=0);n=active.sum(axis=0)
    consistency=np.maximum((b0>0).sum(axis=0),(b0<0).sum(axis=0))/np.maximum(n,1)
    median=np.median(np.abs(b0),axis=0);tab=pd.read_csv(ROOT/'tables/04_stability_all_perturbations.csv')
    for t,target in enumerate(focus):
        f=tab[tab.target==target].set_index('perturbation').loc[pert]
        close(f.selection_frequency,freq[:,t]);close(f.sign_consistency,consistency[:,t]);close(f.median_abs_coefficient,median[:,t])
        close(f.median_coefficient,np.median(b0[:,:,t],axis=0))
        expected=(freq[:,t]>=.8)&(consistency[:,t]>=.9)&(np.abs(full[:,t])>thr)
        assert np.array_equal(f.stable_priority,expected)
        order=f.reset_index().sort_values(['selection_frequency','median_abs_coefficient','full_abs_coefficient'],ascending=False,kind='stable')
        assert order.rank_stability.tolist()==list(range(1,len(pert)+1))
    record('all_stability_frequencies_signs_ranking_and_eligibility',{'coefficient_count':int(b.size),'candidate_rows':len(tab)})
    annotation=pd.read_csv(ROOT/'tables/01_perturbation_annotations.csv').set_index('perturbation_label')
    top=pd.read_csv(ROOT/'tables/04_top10_TF_candidates.csv')
    for target in focus:
        ranked=tab[tab.target==target].sort_values('rank_stability')
        expected=ranked[ranked.primary_TF_eligible & ranked.full_active].head(10).perturbation.tolist()
        assert top[top.target==target].perturbation.tolist()==expected
        for _,row in ranked.iterrows():assert bool(row.primary_TF_eligible)==bool(annotation.loc[row.perturbation,'primary_TF_eligible'])
    record('TF_identity_and_top10_rule',{'displayed_TF_entries':len(top)})
    own_results=pd.read_csv(ROOT/'tables/04_self_readout_results.csv')
    gene_ann=pd.read_csv(ROOT/'tables/01_shared_gene_annotations.csv')
    approved=dict(zip(gene_ann.gene_symbol,gene_ann.approved_symbol_annotation_only))
    with np.load(ROOT/'results/stage4/self_readout/baseline.npz') as z:base=z['coef']
    for _,row in own_results.iterrows():
        drop=np.array([i for i,g in enumerate(genes) if g==row.approved_symbol or approved.get(g)==row.approved_symbol],int)
        assert len(drop)==row.n_removed
        if not len(drop):assert row.status=='own_readout_absent';continue
        with np.load(ROOT/'results/stage4'/f'{row.job_id}.npz') as z:
            assert np.array_equal(z['train'],np.setdiff1d(np.arange(len(X)),drop))
            j=list(pert).index(row.perturbation);t=focus.index(row.target);value=z['coef'][j,t]
        close(row.deleted_coefficient,value);close(row.coefficient_change,value-base[j,t])
    record('own_readout_identity_matching_and_coefficient_changes',{'candidates':len(own_results)})
    assert json.loads((ROOT/'logs/stage4_solver_crosscheck.json').read_text())['passed']
    record('independent_CPU_model_references',{'record':'logs/stage4_solver_crosscheck.json'})

def figures():
    import xml.etree.ElementTree as ET
    from PIL import Image
    from PyPDF2 import PdfReader
    manifest=json.loads((ROOT/'figures/04_figure_manifest.json').read_text())
    count=0
    for item in manifest:
        base=ROOT/'figures'/item['stem']
        assert len(PdfReader(base.with_suffix('.pdf')).pages)==1
        ET.parse(base.with_suffix('.svg'))
        with Image.open(base.with_suffix('.png')) as im:
            assert im.info['dpi'][0]>=449;im.verify()
        count+=3
    record('publication_figure_files',{'files':count})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--figures-only',action='store_true');args=p.parse_args();start=time.monotonic()
    with threadpool_limits(limits=4):
        if args.figures_only:figures()
        else:numeric()
    out={'passed':True,'seconds':time.monotonic()-start,'checks':CHECKS}
    name='stage4_figure_validation.json' if args.figures_only else 'stage4_validation.json'
    (ROOT/'logs'/name).write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
