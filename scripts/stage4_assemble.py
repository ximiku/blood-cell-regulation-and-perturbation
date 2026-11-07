"""Nested selection, held-out predictions, stability ranking and self-readout jobs."""
from pathlib import Path
import argparse,json
import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'config/stage4_config.json').read_text())
FOCUS=CFG['focus_targets'];THR=CFG['active_threshold']
TABLES=ROOT/'tables';OUT=ROOT/'results/stage4/final'

def frame_save(name,frame):frame.to_csv(TABLES/name,index=False,encoding='utf-8-sig',float_format='%.17g')
def load_context(dataset,name):
    p=ROOT/'results/stage4'/dataset/name
    with np.load(p.with_suffix('.npz')) as z:a={k:z[k] for k in z.files}
    return a,json.loads(p.with_suffix('.json').read_text())
def metrics(y,p):
    sst=np.sum((y-y.mean())**2)
    return {'R2':float(1-np.sum((y-p)**2)/sst),'RMSE':float(np.sqrt(np.mean((y-p)**2))),
        'Pearson':float(np.corrcoef(y,p)[0,1]) if np.ptp(p)>0 else np.nan,
        'Spearman':float(spearmanr(y,p).statistic) if np.ptp(p)>0 else np.nan,'n_genes':len(y)}
def select(mse,parameters,model):
    means=mse.mean(axis=0);se=mse.std(axis=0,ddof=1)/np.sqrt(len(mse))
    best=int(np.argmin(means));allowed=np.flatnonzero(means<=means[best]+se[best])
    def key(j):
        p=parameters[j]
        if 'relative_alpha' in p:return (-p['relative_alpha'],-p.get('l1_ratio',1))
        if 'lambda' in p:return (-p['lambda'],)
        return (p['components'],)
    one=min(allowed,key=key)
    return [best,int(one)],means,se

def assemble_models(include_datasets=('raw','landmark','gene_scaled')):
    OUT.mkdir(parents=True,exist_ok=True)
    with h5py.File(ROOT/'data/processed/analysis_matrices.h5') as h:
        X=h['model/X'][:].astype(float);allY=h['model/Y'][:]
        alltargets=h['model/target_cells'].asstr()[:].tolist();focus=[alltargets.index(t) for t in FOCUS]
        genes=h['model/genes'].asstr()[:].astype(str);pert=h['model/perturbations'].asstr()[:].astype(str)
    with np.load(ROOT/'data/processed/stage4_input_metadata.npz') as z:
        scale=z['gene_response_scale'];landmark=z['landmark_indices']
    Yraw=allY[:,focus];split=json.loads((ROOT/'data/splits/gene_cv_splits.json').read_text())
    annot=pd.read_csv(TABLES/'01_perturbation_annotations.csv').rename(columns={'perturbation_label':'perturbation'})
    all_metrics=[];folds=[];choices=[];curves=[];coefs=[];datasets={}
    for dataset in include_datasets:
        Y=Yraw/scale[:,None] if dataset=='gene_scaled' else Yraw
        ids=landmark if dataset=='landmark' else np.arange(len(X));positions={int(g):j for j,g in enumerate(ids)}
        models=CFG['raw_models'] if dataset=='raw' else ['Lasso']
        prediction=pd.DataFrame({'gene_index':np.repeat(ids,4),'gene_symbol':np.repeat(genes[ids],4),
            'target':np.tile(FOCUS,len(ids)),'observed':Y[ids].reshape(-1)})
        outer_fold=np.zeros(len(ids),int);baseline=np.zeros((len(ids),4))
        for e in split['outer']:
            tr=np.array([i for i in e['train'] if i in positions]);te=np.array([i for i in e['test'] if i in positions])
            loc=[positions[int(g)] for g in te];outer_fold[loc]=e['fold'];baseline[loc]=Y[tr].mean(axis=0)
        prediction['outer_fold']=np.repeat(outer_fold,4);prediction['Mean_baseline']=baseline.reshape(-1)
        for model in models:
            pred=np.full((len(ids),4,2),np.nan);outer_b=[];outer_i=[];chosen_all=[]
            for scope in [f'outer_{k:02d}' for k in range(1,6)]+['full']:
                names=[f'{scope}_inner_{j:02d}' for j in range(1,4)]
                inners=[load_context(dataset,name) for name in names]
                fit,meta=load_context(dataset,scope+'_fit');pars=meta['parameters'][model]
                coef=fit[model+'_coef'];chosen=np.zeros((4,2),int)
                b=np.empty((X.shape[1],4,2));intercept=np.empty((4,2));std_b=np.empty_like(b)
                for t,target in enumerate(FOCUS):
                    mse=np.stack([a[model+'_validation_MSE'][t] for a,_ in inners])
                    selected,means,se=select(mse,pars,model);chosen[t]=selected
                    for j,param in enumerate(pars):
                        curves.append({'dataset':dataset,'target':target,'model':model,'scope':scope,'parameter_index':j,
                            **param,'mean_MSE':means[j],'SE_MSE':se[j],
                            **{f'fold_{k+1}_MSE':mse[k,j] for k in range(3)},
                            'selected_min':j==selected[0],'selected_1se':j==selected[1]})
                    for choice,j in zip(['min','1se'],selected):
                        ix=0 if choice=='min' else 1
                        std_b[:,t,ix]=coef[:,t,j]
                        b[:,t,ix]=coef[:,t,j]*fit['y_scale'][t]/fit['x_scale']
                        intercept[t,ix]=fit['y_mean'][t]-fit['x_mean']@b[:,t,ix]
                        choices.append({'dataset':dataset,'target':target,'model':model,'selection':choice,'scope':scope,
                            'parameter_index':j,**pars[j],'n_active':int(np.sum(np.abs(coef[:,t,j])>THR)),
                            'n_train':len(fit['train']),'inner_mean_MSE':means[j]})
                chosen_all.append(chosen)
                if scope!='full':
                    e=split['outer'][int(scope[-2:])-1];test=[i for i in e['test'] if i in positions];loc=[positions[i] for i in test]
                    pred[loc]=np.einsum('np,ptk->ntk',X[test],b)+intercept[None,:,:]
                    outer_b.append(b);outer_i.append(intercept)
                else:
                    np.savez_compressed(OUT/f'{dataset}_{model}.npz',gene_indices=ids,genes=genes[ids],targets=FOCUS,
                        perturbations=pert,oof_predictions=pred,outer_fold=outer_fold,full_coef_standardized=std_b,
                        full_coef_input=b,full_intercept_input=intercept,outer_coef_input=np.stack(outer_b),
                        outer_intercept_input=np.stack(outer_i),chosen_indices=np.stack(chosen_all))
                    for t,target in enumerate(FOCUS):
                        for j,choice in enumerate(['min','1se']):
                            f=pd.DataFrame({'dataset':dataset,'target':target,'model':model,'selection':choice,'perturbation':pert,
                                'coefficient_standardized':std_b[:,t,j],'coefficient_input':b[:,t,j],
                                'active':np.abs(std_b[:,t,j])>THR})
                            coefs.append(f.merge(annot,on='perturbation',validate='one_to_one'))
            assert np.isfinite(pred).all()
            for j,choice in enumerate(['min','1se']):
                prediction[model+'_'+choice]=pred[:,:,j].reshape(-1)
                for t,target in enumerate(FOCUS):
                    all_metrics.append({'dataset':dataset,'target':target,'model':model,'selection':choice,**metrics(Y[ids,t],pred[:,t,j])})
                    for k in range(1,6):
                        m=outer_fold==k
                        folds.append({'dataset':dataset,'target':target,'model':model,'selection':choice,'fold':k,**metrics(Y[ids[m],t],pred[m,t,j])})
        for t,target in enumerate(FOCUS):
            all_metrics.append({'dataset':dataset,'target':target,'model':'Mean_baseline','selection':'baseline',**metrics(Y[ids,t],baseline[:,t])})
        if dataset=='raw':
            for j,choice in enumerate(['min','1se']):
                values=[]
                for t,idx in enumerate(focus):
                    df=pd.read_csv(ROOT/f'results/lasso/target_{idx:02d}/oof_predictions.csv')
                    assert np.array_equal(df.gene_symbol.to_numpy(),genes)
                    yhat=df['Lasso_'+choice].to_numpy();values.append(yhat)
                    all_metrics.append({'dataset':dataset,'target':FOCUS[t],'model':'Lasso','selection':choice,**metrics(Y[:,t],yhat)})
                    for k in range(1,6):
                        m=outer_fold==k;folds.append({'dataset':dataset,'target':FOCUS[t],'model':'Lasso','selection':choice,'fold':k,**metrics(Y[m,t],yhat[m])})
                prediction['Lasso_'+choice]=np.stack(values,axis=1).reshape(-1)
        frame_save(f'04_{dataset}_oof_predictions.csv',prediction);datasets[dataset]=prediction
    reference=pd.read_csv(TABLES/'03_all_coefficients_annotated.csv')
    reference=reference[reference.target.isin(FOCUS)].copy();reference['dataset']='raw'
    reference['selection']=reference.model.str.replace('Lasso_','');reference['model']='Lasso';reference['active']=reference.active_gt_1e_6
    allcoef=pd.concat(coefs+[reference],ignore_index=True)
    frame_save('04_all_coefficients_annotated.csv',allcoef)
    metric=pd.DataFrame(all_metrics);frame_save('04_model_metrics.csv',metric[metric.dataset=='raw'])
    frame_save('04_sensitivity_metrics.csv',metric[metric.dataset!='raw'])
    frame_save('04_fold_metrics.csv',pd.DataFrame(folds));frame_save('04_selected_models.csv',pd.DataFrame(choices))
    frame_save('04_tuning_curves.csv',pd.DataFrame(curves))
    # Matched landmark test genes from the original, all-gene trained Lasso.
    records=[]
    if 'raw' in datasets:
        for target in FOCUS:
            full=datasets['raw'];sub=full[(full.target==target)&full.gene_index.isin(landmark)]
            for choice in ['min','1se']:
                records.append({'dataset':'raw_model_on_landmark_test_genes','target':target,'model':'Lasso','selection':choice,
                    **metrics(sub.observed.to_numpy(),sub['Lasso_'+choice].to_numpy())})
        frame_save('04_matched_landmark_reference_metrics.csv',pd.DataFrame(records))
    print(json.dumps({'event':'models_assembled','raw_metric_rows':int((metric.dataset=='raw').sum()),'coefficient_rows':len(allcoef)}),flush=True)

def assemble_stability():
    annot=pd.read_csv(TABLES/'01_perturbation_annotations.csv').rename(columns={'perturbation_label':'perturbation'})
    with np.load(ROOT/'data/processed/stage4_stability_reference.npz') as z:full=z['full_coefficients']
    with np.load(ROOT/'data/processed/stage4_input_metadata.npz') as z:pert=z['perturbations'];genes=z['genes']
    reps=[];diagnostics=[]
    for i in range(1,101):
        with np.load(ROOT/f'results/stage4/stability/repeat_{i:03d}.npz') as z:
            reps.append(z['coef']);d=z['diagnostics'][0]
            for t,target in enumerate(FOCUS):
                diagnostics.append({'repeat':i,'target':target,'alpha':d[0,t],'dual_gap':d[1,t],
                    'primal':d[2,t],'dual':d[3,t],'n_iter':int(d[4,t])})
    b=np.stack(reps);active=np.abs(b)>THR;bactive=np.where(active,b,0)
    count=active.sum(axis=0);freq=count/100
    sign=np.maximum((bactive>0).sum(axis=0),(bactive<0).sum(axis=0))/np.maximum(count,1)
    median=np.median(np.abs(bactive),axis=0);records=[]
    for t,target in enumerate(FOCUS):
        f=pd.DataFrame({'target':target,'perturbation':pert,'selection_frequency':freq[:,t],
            'sign_consistency':sign[:,t],'median_abs_coefficient':median[:,t],
            'median_coefficient':np.median(bactive[:,:,t],axis=0),
            'full_sign_agreement_among_selected':np.sum(active[:,:,t] & (np.sign(bactive[:,:,t])==np.sign(full[:,t])[None,:]),axis=0)/np.maximum(count[:,t],1),
            'full_coefficient':full[:,t],'full_abs_coefficient':np.abs(full[:,t]),
            'full_active':np.abs(full[:,t])>THR,'positive_frequency':(bactive[:,:,t]>0).mean(axis=0),
            'negative_frequency':(bactive[:,:,t]<0).mean(axis=0),
            'coefficient_q025':np.quantile(bactive[:,:,t],.025,axis=0),
            'coefficient_q975':np.quantile(bactive[:,:,t],.975,axis=0)})
        f['stable_priority']=(f.selection_frequency>=.8)&(f.sign_consistency>=.9)&f.full_active
        f=f.merge(annot,on='perturbation',validate='one_to_one').sort_values(
            ['selection_frequency','median_abs_coefficient','full_abs_coefficient'],ascending=False,kind='stable')
        f['rank_stability']=np.arange(1,len(f)+1);records.append(f)
    allf=pd.concat(records,ignore_index=True);tf=allf[allf.primary_TF_eligible].copy()
    tf['TF_stability_rank']=tf.groupby('target').cumcount()+1
    top=tf[tf.full_active].groupby('target',sort=False).head(10)
    frame_save('04_stability_all_perturbations.csv',allf);frame_save('04_stability_all_TFs.csv',tf)
    frame_save('04_stable_priority_TFs.csv',tf[tf.stable_priority]);frame_save('04_top10_TF_candidates.csv',top)
    frame_save('04_stability_diagnostics.csv',pd.DataFrame(diagnostics))
    np.savez_compressed(OUT/'stability_coefficients.npz',coefficients=b,targets=FOCUS,perturbations=pert)
    trial=ROOT/'results/stage4_exploratory/warmstart_stability'
    if trial.exists():
        old=[]
        for i in range(1,101):
            with np.load(trial/f'repeat_{i:03d}.npz') as z:old.append(z['coef'])
        old=np.stack(old);oldfreq=(np.abs(old)>THR).mean(axis=0)
        init=[]
        for t,target in enumerate(FOCUS):
            init.append(pd.DataFrame({'target':target,'perturbation':pert,
                'zero_start_frequency':freq[:,t],'full_coefficient_start_frequency':oldfreq[:,t],
                'frequency_change':freq[:,t]-oldfreq[:,t]}))
        frame_save('04_initialization_sensitivity.csv',pd.concat(init,ignore_index=True))
    # Match current approved identities, retain every matching readout row.
    gene_ann=pd.read_csv(TABLES/'01_shared_gene_annotations.csv')
    approved=dict(zip(gene_ann.gene_symbol,gene_ann.approved_symbol_annotation_only))
    own=[];jobs={};allids=np.arange(len(genes));jobs['baseline']={'id':'self_readout/baseline','name':'baseline',
        'kind':'self_readout','dataset':'raw','train':allids.tolist(),'validation':[],'tolerance':1e-6}
    for _,row in top.iterrows():
        symbol=row.approved_symbol
        drop=[i for i,g in enumerate(genes) if g==symbol or approved.get(g)==symbol]
        key='exclude_'+str(symbol)
        own.append({'target':row.target,'perturbation':row.perturbation,'approved_symbol':symbol,
            'readout_genes':';'.join(genes[drop]),'n_removed':len(drop),
            'job_id':'self_readout/'+key if drop else '', 'status':'planned' if drop else 'own_readout_absent'})
        if drop:
            jobs[key]={'id':'self_readout/'+key,'name':key,'kind':'self_readout','dataset':'raw',
                'train':np.setdiff1d(allids,drop).tolist(),'validation':[],'tolerance':1e-6}
    (ROOT/'data/splits/stage4_self_readout_jobs.json').write_text(json.dumps({'jobs':list(jobs.values())},indent=2)+'\n')
    frame_save('04_self_readout_match_plan.csv',pd.DataFrame(own))
    print(json.dumps({'event':'stability_assembled','stable_TFs':tf.groupby('target').stable_priority.sum().to_dict(),
                      'self_readout_jobs':len(jobs)},ensure_ascii=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stability-only',action='store_true');p.add_argument('--models-only',action='store_true')
    p.add_argument('--datasets',nargs='+',default=['raw','landmark','gene_scaled']);args=p.parse_args()
    if not args.stability_only:assemble_models(args.datasets)
    if not args.models_only:assemble_stability()
