"""Validate and summarize the 20 fresh-tuning core resamples."""
from pathlib import Path
import json
import numpy as np,pandas as pd,h5py
R=Path(__file__).resolve().parents[1];cfg=json.loads((R/'config/stage4_config.json').read_text());targets=cfg['focus_targets'];saved_samples=np.load(R/'data/splits/stage4_stability_indices.npz')['train_indices'][:20];bs=[];fixed=[];selections=[];checks=[]
for rep in range(20):
 z=np.load(R/f'results/stage7/core_retuning/repeat_{rep+1:02d}.npz');meta=json.loads((R/f'results/stage7/core_retuning/repeat_{rep+1:02d}.json').read_text());assert meta['all_converged'] and np.array_equal(z['train_indices'],saved_samples[rep]);tr=set(z['train_indices'].tolist());vals=[]
 for k in [1,2,3]:
  inner=set(z[f'inner{k}_train'].tolist());val=set(z[f'inner{k}_validation'].tolist());assert not inner&val and inner|val==tr;vals.extend(val);diag=z[f'inner{k}_diagnostics'];assert np.isfinite(diag).all() and (diag[:,1,:]>=-1e-8).all() and (diag[:,1,:]<=cfg['tol']*1.00001).all()
 assert len(vals)==len(tr) and len(set(vals))==len(tr)
 err=z['cv_MSE'];mean=err.mean(0);se=err.std(0,ddof=1)/np.sqrt(3);mini=mean.argmin(1);one=np.array([np.flatnonzero(mean[j]<=mean[j,mini[j]]+se[j,mini[j]])[0] for j in range(4)]);assert np.array_equal(z['selected_indices'],np.stack([mini,one],1))
 for mode in ['min','1se']:
  diag=z['diagnostics_'+mode];assert np.isfinite(diag).all() and (diag[:,1,:]>=-1e-8).all() and (diag[:,1,:]<=cfg['tol']*1.00001).all()
 bs.append(z['coef_1se']);f=np.load(R/f'results/stage4/stability/repeat_{rep+1:03d}.npz');assert np.array_equal(f['train'],z['train_indices']);fixed.append(f['coef'])
 for j,t in enumerate(targets):selections.append({'repeat':rep+1,'target':t,'relative_alpha_min':z['relative_path'][mini[j]],'relative_alpha_1se':z['relative_path'][one[j]],'n_nonzero_min':int((np.abs(z['coef_min'][:,j])>cfg['active_threshold']).sum()),'n_nonzero_1se':int((np.abs(z['coef_1se'][:,j])>cfg['active_threshold']).sum()),'seconds_joint_four_targets':meta['seconds']})
 checks.append({'repeat':rep+1,'split_selection_and_all_certificates_passed':True})
B=np.stack(bs);F=np.stack(fixed);threshold=cfg['active_threshold'];annotations=pd.read_csv(R/'tables/01_perturbation_annotations.csv');rows=[]
with h5py.File(R/'data/processed/analysis_matrices.h5') as h:pert=h['model/perturbations'].asstr()[:];is_tf=h['model/is_TF'][:]
for j,t in enumerate(targets):
 for p in np.flatnonzero(is_tf):
  b=B[:,p,j];f=F[:,p,j];active=np.abs(b)>threshold;signs=np.sign(b[active]);consistency=max(np.sum(signs>0),np.sum(signs<0))/len(signs) if len(signs) else np.nan
  rows.append({'target':t,'perturbation':pert[p],'retuned_selection_frequency':float(active.mean()),'fixed_same20_selection_frequency':float((np.abs(f)>threshold).mean()),'retuned_sign_consistency':float(consistency),'retuned_median_coefficient':float(np.median(b)),'retuned_median_absolute_coefficient':float(np.median(np.abs(b)))})
d=pd.DataFrame(rows);d.to_csv(R/'tables/07_core_retuned_TF_stability.csv',index=False);pd.DataFrame(selections).to_csv(R/'tables/07_core_retuned_parameter_choices.csv',index=False);top=pd.read_csv(R/'tables/04_top10_TF_candidates.csv').merge(d,on=['target','perturbation'],validate='one_to_one');top.to_csv(R/'tables/07_core_top10_retuning_check.csv',index=False)
summary={'passed':True,'resamples':20,'targets':4,'inner_folds':3,'checks':checks,'display_top40_min_retuned_frequency':float(top.retuned_selection_frequency.min()),'display_top40_same_median_sign':int((np.sign(top.retuned_median_coefficient)==np.sign(top.full_coefficient)).sum())};(R/'logs/stage7_core_retuning_validation.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n');print(json.dumps(summary,ensure_ascii=False))
