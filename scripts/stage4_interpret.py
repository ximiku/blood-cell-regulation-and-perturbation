"""Candidate overlap, correlated predictors, scale sensitivity and own-readout refits."""
from pathlib import Path
import json
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits
from stage4_assemble import metrics,frame_save,FOCUS,THR

ROOT=Path(__file__).resolve().parents[1];TAB=ROOT/'tables'

def main():
    coef=pd.read_csv(TAB/'04_all_coefficients_annotated.csv',low_memory=False)
    stable=pd.read_csv(TAB/'04_stability_all_TFs.csv');top=pd.read_csv(TAB/'04_top10_TF_candidates.csv')
    with h5py.File(ROOT/'data/processed/analysis_matrices.h5') as f:
        X=f['model/X'][:].astype(float);allY=f['model/Y'][:]
        targets=f['model/target_cells'].asstr()[:].tolist();genes=f['model/genes'].asstr()[:]
        pert=f['model/perturbations'].asstr()[:]
    Y=allY[:,[targets.index(t) for t in FOCUS]];pindex={p:i for i,p in enumerate(pert)}
    with np.load(ROOT/'data/processed/stage4_input_metadata.npz') as z:scale=z['gene_response_scale']
    # Exploratory back-transform uses supplied scales containing target-cell expression.
    # Keep its provenance explicit; formal scale sensitivity uses the native response.
    rec=[]
    for dataset in ['raw','gene_scaled','landmark']:
        pred=pd.read_csv(TAB/f'04_{dataset}_oof_predictions.csv')
        for t,target in enumerate(FOCUS):
            sub=pred[pred.target==target];ix=sub.gene_index.to_numpy()
            for choice in ['min','1se']:
                yhat=sub['Lasso_'+choice].to_numpy()
                if dataset=='gene_scaled':yhat=yhat*scale[ix]
                rec.append({'dataset':dataset,'target':target,'selection':choice,
                    'evaluation_scale':'original_input_units',
                    'role':'descriptive_backtransform_with_observed_target_dependent_scales' if dataset=='gene_scaled' else 'original_response_evaluation',
                    'used_for_primary_model_ranking':dataset=='raw',**metrics(Y[ix,t],yhat)})
    frame_save('04_sensitivity_original_units.csv',pd.DataFrame(rec))
    agreement=[];detail=[]
    for target in FOCUS:
        raw=coef[(coef.dataset=='raw')&(coef.model=='Lasso')&(coef.selection=='1se')&(coef.target==target)].set_index('perturbation').loc[pert]
        for dataset,model in [('raw','ElasticNet'),('landmark','Lasso'),('gene_scaled','Lasso')]:
            alt=coef[(coef.dataset==dataset)&(coef.model==model)&(coef.selection=='1se')&(coef.target==target)].set_index('perturbation').loc[pert]
            rb=raw.coefficient_standardized.to_numpy();ab=alt.coefficient_standardized.to_numpy()
            r=np.abs(rb)>THR;a=np.abs(ab)>THR;both=r&a
            agreement.append({'target':target,'comparison':dataset+'_'+model,'n_active_reference':int(r.sum()),
                'n_active_alternative':int(a.sum()),'n_intersection':int(both.sum()),
                'active_Jaccard':float(both.sum()/np.sum(r|a)) if np.any(r|a) else np.nan,
                'sign_agreement_in_intersection':float(np.mean(np.sign(rb[both])==np.sign(ab[both]))) if both.any() else np.nan,
                'coefficient_Spearman':float(spearmanr(rb,ab).statistic) if np.ptp(rb)>0 and np.ptp(ab)>0 else np.nan,
                'comparison_status':'no_active_alternative' if not a.any() else 'active_alternative'})
            for _,row in top[top.target==target].iterrows():
                i=pindex[row.perturbation]
                detail.append({'target':target,'perturbation':row.perturbation,'comparison':dataset+'_'+model,
                    'reference_coefficient':rb[i],'alternative_coefficient':ab[i],'alternative_active':bool(a[i]),
                    'same_active_sign':bool(a[i] and np.sign(rb[i])==np.sign(ab[i]))})
    frame_save('04_candidate_set_agreement.csv',pd.DataFrame(agreement));frame_save('04_top10_candidate_sensitivity.csv',pd.DataFrame(detail))
    overlap=[]
    for left in FOCUS:
        ls=set(stable[(stable.target==left)&stable.stable_priority].approved_symbol.dropna())
        for right in FOCUS:
            rs=set(stable[(stable.target==right)&stable.stable_priority].approved_symbol.dropna())
            overlap.append({'left_target':left,'right_target':right,'left_n_TF_identities':len(ls),
                'right_n_TF_identities':len(rs),'intersection':len(ls&rs),'union':len(ls|rs),
                'Jaccard':len(ls&rs)/len(ls|rs) if ls|rs else 0,'shared_TFs':';'.join(sorted(ls&rs))})
    frame_save('04_stable_TF_overlap.csv',pd.DataFrame(overlap))
    # Predictor groups are descriptive, fitted after predictive evaluation is complete.
    Z=(X-X.mean(axis=0))/X.std(axis=0);R=Z.T@Z/len(X)
    groups=[];edges=[];group_metric=[]
    for threshold in [.7,.8,.9]:
        mask=np.abs(R)>=threshold;np.fill_diagonal(mask,False)
        n,labels=connected_components(csr_matrix(mask),directed=False)
        sizes=np.bincount(labels)
        for i,p in enumerate(pert):groups.append({'threshold':threshold,'perturbation':p,'group_id':int(labels[i]),'group_size':int(sizes[labels[i]])})
        if threshold==.8:
            ii,jj=np.where(np.triu(mask,1))
            edges=[{'perturbation_1':pert[i],'perturbation_2':pert[j],'Pearson':R[i,j]} for i,j in zip(ii,jj)]
            for target in FOCUS:
                for model in ['Lasso','ElasticNet']:
                    b=coef[(coef.dataset=='raw')&(coef.target==target)&(coef.model==model)&(coef.selection=='1se')].set_index('perturbation').loc[pert].coefficient_standardized.to_numpy()
                    for g in range(n):
                        ids=np.flatnonzero(labels==g);active=np.abs(b[ids])>THR
                        group_metric.append({'target':target,'model':model,'group_id':g,'group_size':len(ids),
                            'n_active':int(active.sum()),'sum_absolute_coefficients':float(np.abs(b[ids]).sum()),
                            'group_prediction_SD':float(np.sqrt(max(b[ids]@R[np.ix_(ids,ids)]@b[ids],0))),
                            'active_members':';'.join(pert[ids][active])})
    frame_save('04_correlation_groups.csv',pd.DataFrame(groups));frame_save('04_correlation_edges.csv',pd.DataFrame(edges,columns=['perturbation_1','perturbation_2','Pearson']))
    frame_save('04_group_model_comparison.csv',pd.DataFrame(group_metric))
    partners=[]
    for _,row in top.iterrows():
        i=pindex[row.perturbation];order=np.argsort(-np.abs(R[i]));order=[j for j in order if j!=i][:5]
        for rank,j in enumerate(order,1):partners.append({'target':row.target,'candidate':row.perturbation,
            'partner':pert[j],'partner_rank':rank,'Pearson':R[i,j],'above_primary_threshold':abs(R[i,j])>=.8})
    frame_save('04_candidate_correlated_partners.csv',pd.DataFrame(partners))
    # Own-readout comparison uses matched, tighter-tolerance baseline and deletion fits.
    match=pd.read_csv(TAB/'04_self_readout_match_plan.csv');self_records=[]
    with np.load(ROOT/'results/stage4/self_readout/baseline.npz') as z:baseline=z['coef']
    for _,row in match.iterrows():
        record=row.to_dict();t=FOCUS.index(row.target);i=pindex[row.perturbation]
        record['baseline_tight_coefficient']=baseline[i,t]
        if row.n_removed:
            with np.load(ROOT/'results/stage4'/f'{row.job_id}.npz') as z:b=z['coef'];tr=z['train'];diag=z['diagnostics'][0]
            before=baseline[i,t];after=b[i,t]
            record.update(status='completed',deleted_coefficient=after,coefficient_change=after-before,
                relative_coefficient_change=(after-before)/abs(before) if abs(before)>THR else np.nan,
                same_active_sign=bool(abs(after)>THR and np.sign(after)==np.sign(before)),
                full_marginal_Pearson=float(np.corrcoef(X[:,i],Y[:,t])[0,1]),
                deleted_marginal_Pearson=float(np.corrcoef(X[tr,i],Y[tr,t])[0,1]),
                dual_gap=diag[1,t],tol=1e-6)
        self_records.append(record)
    own_final=pd.DataFrame(self_records);frame_save('04_self_readout_results.csv',own_final)
    evidence=top[['target','perturbation','approved_symbol','selection_frequency','sign_consistency',
        'full_coefficient','stable_priority']].copy()
    for comparison in ['raw_ElasticNet','landmark_Lasso','gene_scaled_Lasso']:
        f=pd.DataFrame(detail);f=f[f.comparison==comparison][['target','perturbation','same_active_sign']]
        evidence=evidence.merge(f.rename(columns={'same_active_sign':comparison+'_same_active_sign'}),on=['target','perturbation'],validate='one_to_one')
    evidence=evidence.merge(own_final[['target','perturbation','status','same_active_sign','relative_coefficient_change']].rename(
        columns={'status':'self_readout_status','same_active_sign':'self_readout_same_active_sign'}),on=['target','perturbation'],validate='one_to_one')
    frame_save('04_candidate_evidence_matrix.csv',evidence)
    summary={'stable_TF_entries':stable.groupby('target').stable_priority.sum().astype(int).to_dict(),
        'top10_candidates':{t:top[top.target==t].perturbation.tolist() for t in FOCUS},
        'self_readout_completed':int(sum(r['status']=='completed' for r in self_records)),
        'self_readout_absent':int(sum(r['status']=='own_readout_absent' for r in self_records)),
        'correlation_group_counts':{str(th):int(pd.DataFrame(groups).loc[lambda f:f.threshold==th,'group_id'].nunique()) for th in [.7,.8,.9]}}
    (ROOT/'logs/stage4_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=4):main()
