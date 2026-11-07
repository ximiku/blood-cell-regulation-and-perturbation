"""Evaluate shared condition predictions, equal weights, fixed genes and direct-target exclusion."""
from pathlib import Path
import json,argparse
import numpy as np,pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1]
MODELS=['zero_change','training_condition_mean','Ridge_additive','direct_single_additive','GEARS','scGPT_blood','scFoundation_GEARS']
def correlation(x,y):
 if np.ptp(x)<1e-12 or np.ptp(y)<1e-12:return np.nan,np.nan
 return float(np.corrcoef(x,y)[0,1]),float(spearmanr(x,y).statistic)
def main(partial=False):
 P=R/'data/processed/stage6';split=pd.read_csv(R/'data/splits/stage6/condition_splits.csv');ix=np.load(P/'evaluation_indices.npz');rows=[];availability=[]
 for seed in [1,2,3]:
  expected=sorted(split[(split.seed==seed)&(split.split=='test')].condition)
  for model in MODELS:
   out=R/f'results/stage7/{model}/seed{seed}';p=out/'test_predictions.npz'
   if not p.exists():
    if not partial:raise FileNotFoundError(p)
    continue
   with np.load(p) as z:
    names=z['conditions'].astype(str);genes=z['gene_symbols'].astype(str);indices=z['gene_indices'];pred=z['predicted'].astype(np.float64);truth=z['observed'].astype(np.float64);ctrl=z['control_mean'].astype(np.float64)
   assert names.tolist()==expected and np.array_equal(indices,ix['common'])
   pdlt=pred-ctrl;tdlt=truth-ctrl;elig=np.isfinite(pred).all(1)
   assert model=='direct_single_additive' or elig.all()
   availability.append({'model':model,'seed':seed,'eligible':int(elig.sum()),'total':len(names)})
   for i,c in enumerate(names):
    if not elig[i]:continue
    info=split[(split.seed==seed)&(split.condition==c)].iloc[0];pert=set(c.split('+'))-{'ctrl'}
    for gene_set in ['shared','HVG1000']:
     base=np.ones(len(genes),dtype=bool) if gene_set=='shared' else np.isin(indices,ix[f'hvg1000_seed{seed}'])
     for exclude in [False,True]:
      mask=base&(~np.isin(genes,list(pert)) if exclude else True);a=pdlt[i,mask];b=tdlt[i,mask];pear,spear=correlation(a,b)
      rows.append({'model':model,'seed':seed,'condition':c,'group':info.evaluation_group,'official_group':info.official_subgroup,'gene_set':gene_set,'exclude_direct_targets':exclude,'n_genes':int(mask.sum()),'MSE':float(np.mean((a-b)**2)),'Pearson_delta':pear,'Spearman_delta':spear,'observed_delta_RMS':float(np.sqrt(np.mean(b*b)))})
 d=pd.DataFrame(rows);keys=['seed','condition','gene_set','exclude_direct_targets']
 for baseline,label in [('zero_change','zero'),('training_condition_mean','mean'),('Ridge_additive','ridge')]:
  ref=d[d.model==baseline][keys+['MSE']].rename(columns={'MSE':f'MSE_{label}'})
  d=d.merge(ref,on=keys,how='left',validate='many_to_one');d[f'improvement_vs_{label}']=1-d.MSE/d[f'MSE_{label}']
 # Aggregate MSE first; RMSE is sqrt of the condition-equal aggregate MSE.
 summaries=[]
 for group in ['all','seen0','seen1','seen2','unseen_single','direct_additive_eligible']:
  subset=d if group=='all' else d[d.group==group]
  if group=='direct_additive_eligible':
   eligible=d[d.model=='direct_single_additive'][keys].drop_duplicates();subset=d.merge(eligible,on=keys,how='inner',validate='many_to_one')
  for key,s in subset.groupby(['model','seed','gene_set','exclude_direct_targets'],sort=False):
   mse=float(s.MSE.mean());row=dict(zip(['model','seed','gene_set','exclude_direct_targets'],key));row.update(group=group,n_conditions=len(s),MSE=mse,RMSE=float(np.sqrt(mse)),Pearson_delta=float(s.Pearson_delta.mean()),Spearman_delta=float(s.Spearman_delta.mean()),n_defined_correlations=int(s.Pearson_delta.notna().sum()))
   for label in ['zero','mean','ridge']:row[f'improvement_vs_{label}']=1-mse/float(s[f'MSE_{label}'].mean())
   summaries.append(row)
 agg=pd.DataFrame(summaries);root=R/('tables/stage7_partial' if partial else 'tables');root.mkdir(exist_ok=True,parents=True)
 d.to_csv(root/'07_condition_metrics.csv',index=False);agg.to_csv(root/'07_seed_metrics.csv',index=False);pd.DataFrame(availability).to_csv(root/'07_coverage.csv',index=False)
 summary=agg.groupby(['model','gene_set','exclude_direct_targets','group']).agg(seeds=('seed','size'),MSE_mean=('MSE','mean'),MSE_sd=('MSE','std'),RMSE_mean=('RMSE','mean'),Pearson_delta_mean=('Pearson_delta','mean'),Spearman_delta_mean=('Spearman_delta','mean'),improvement_vs_zero_mean=('improvement_vs_zero','mean'),improvement_vs_mean_mean=('improvement_vs_mean','mean'),improvement_vs_ridge_mean=('improvement_vs_ridge','mean')).reset_index();summary.to_csv(root/'07_model_summary.csv',index=False)
 print(summary[(summary.gene_set=='shared')&(~summary.exclude_direct_targets)&(summary.group=='all')].to_string(index=False))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--allow-partial',action='store_true');main(p.parse_args().allow_partial)
