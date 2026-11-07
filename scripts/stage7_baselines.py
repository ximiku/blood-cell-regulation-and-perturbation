"""Condition-equal baselines; select Ridge penalty using validation conditions only."""
from pathlib import Path
import json,os
import numpy as np,pandas as pd
R=Path(os.environ.get('STAGE6_ROOT','/root/HW2_20261107/solution'));P=R/'data/processed/stage6';Q=R/'data/processed/stage7'
conditions=json.loads((Q/'conditions.json').read_text());means=np.load(Q/'condition_means.npy').astype(np.float64);common=np.load(P/'evaluation_indices.npz')['common'];genes=pd.read_csv(P/'genes.csv').gene_symbol.to_numpy()[common];ci={c:i for i,c in enumerate(conditions)}
control=np.load(P/'control_mean.npy').astype(np.float64)[common];delta=means[:,common]-control;pg=json.loads((P/'perturbation_genes.json').read_text());pi={p:i for i,p in enumerate(pg)}
design=np.zeros((len(conditions),len(pg)))
for i,c in enumerate(conditions):
 for g in c.split('+'):
  if g!='ctrl':design[i,pi[g]]=1
for seed in [1,2,3]:
 split=json.loads((R/f'data/splits/stage6/seed{seed}.json').read_text())['conditions'];tr=[ci[c] for c in split['train'] if c!='ctrl'];va=[ci[c] for c in split['val'] if c!='ctrl'];te=[ci[c] for c in split['test'] if c!='ctrl'];rows=[]
 U,s,VT=np.linalg.svd(design[tr],full_matrices=False);projection=U.T@delta[tr];best=None
 for alpha in [.0001,.001,.01,.1,1,10,100,1000,10000]:
  coef=VT.T@((s/(s*s+alpha))[:,None]*projection);error=float(np.mean((design[va]@coef-delta[va])**2));rows.append({'alpha':alpha,'validation_MSE':error})
  if best is None or error<best[0]:best=(error,alpha,coef)
 preds={'zero_change':np.zeros_like(delta[te]),'training_condition_mean':np.repeat(delta[tr].mean(0,keepdims=True),len(te),axis=0),'Ridge_additive':design[te]@best[2]}
 direct=np.full_like(delta[te],np.nan)
 for j,i in enumerate(te):
  singles=[g+'+ctrl' for g in conditions[i].split('+') if g!='ctrl']
  if all(c in split['train'] for c in singles):direct[j]=sum(delta[ci[c]] for c in singles)
 preds['direct_single_additive']=direct
 for name,pred in preds.items():
  out=R/f'results/stage7/{name}/seed{seed}';out.mkdir(parents=True,exist_ok=True)
  np.savez_compressed(out/'test_predictions.npz',predicted=(pred+control).astype(np.float32),observed=means[te][:,common].astype(np.float32),conditions=np.array([conditions[i] for i in te]),gene_indices=common,gene_symbols=genes.astype(str),control_mean=control.astype(np.float32),eligible=np.isfinite(pred).all(1))
  result={'model':name,'seed':seed,'status':'completed','test_conditions':len(te),'eligible_conditions':int(np.isfinite(pred).all(1).sum()),'training_conditions_equal_weight':len(tr),'controls_excluded_from_training_mean':True}
  if name=='Ridge_additive':
   result.update({'alpha':best[1],'validation_MSE':best[0],'fit_intercept':False});pd.DataFrame(rows).to_csv(out/'penalty_validation.csv',index=False);np.savez_compressed(out/'ridge_coefficients.npz',coef=best[2],perturbation_genes=np.array(pg),gene_indices=common)
  (out/'completed.json').write_text(json.dumps(result,indent=2)+'\n')
 print('baselines complete',seed,flush=True)
