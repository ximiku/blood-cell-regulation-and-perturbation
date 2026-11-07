"""Check full-run provenance, selected checkpoints, train/test isolation and reported metrics."""
from pathlib import Path
import json,hashlib,zipfile,time
import numpy as np,pandas as pd,scipy.sparse as sp
R=Path(__file__).resolve().parents[1];P=R/'data/processed/stage6';Q=R/'data/processed/stage7';checks=[]
def check(name,value,**details):
 checks.append({'check':name,'passed':bool(value),**details});assert value,name
def main():
 t=time.time();cells=pd.read_csv(P/'cells.csv.gz');genes=pd.read_csv(P/'genes.csv');ix=np.load(P/'evaluation_indices.npz');common=ix['common'];cond=json.loads((Q/'conditions.json').read_text());means=np.load(Q/'condition_means.npy');ci={c:i for i,c in enumerate(cond)};control=np.load(P/'control_mean.npy');refs=np.load(Q/'reference_control_rows.npy');X=sp.load_npz(P/'expression_log1p.npz')
 expected_refs=np.random.default_rng(20261107).choice(np.flatnonzero(cells.condition.to_numpy()=='ctrl'),64,replace=False)
 check('shared_reference_pool_recomputed',np.array_equal(refs,expected_refs) and np.array_equal(X[refs].toarray(),np.load(Q/'reference_expression.npy')) and np.array_equal(cells.total_count.to_numpy()[refs],np.load(Q/'reference_total_counts.npy')))
 check('reference_controls_train_all_splits',all((cells.iloc[refs][f'split_seed{s}']=='train').all() for s in [1,2,3]))
 # Recompute all condition means from the frozen sparse matrix.
 codes=np.load(Q/'cell_condition_codes.npy');recomputed=np.stack([np.asarray(X[codes==i].mean(0)).ravel() for i in range(len(cond))]);check('condition_means_recomputed',np.array_equal(recomputed,means))
 del X,recomputed
 metrics=pd.read_csv(R/'tables/07_condition_metrics.csv');aggregate=pd.read_csv(R/'tables/07_seed_metrics.csv')
 gpu=['GEARS','scGPT_blood','scFoundation_GEARS'];baseline=['zero_change','training_condition_mean','Ridge_additive','direct_single_additive']
 for seed in [1,2,3]:
  split=json.loads((R/f'data/splits/stage6/seed{seed}.json').read_text())['conditions'];test=sorted(split['test']);val=sorted(split['val']);train=split['train'];ntrain=int((cells[f'split_seed{seed}']=='train').sum())
  check(f'seed{seed}_partition',not(set(train)&set(test) or set(train)&set(val) or set(val)&set(test)) and set.union(set(train),set(val),set(test))==set(cond))
  for model in baseline+gpu:
   out=R/f'results/stage7/{model}/seed{seed}';z=np.load(out/'test_predictions.npz');result=json.loads((out/'completed.json').read_text());prediction=z['predicted'];truth=z['observed'];names=z['conditions'].astype(str);eligible=np.isfinite(prediction).all(1)
   check(f'{model}_seed{seed}_prediction_alignment',names.tolist()==test and np.array_equal(z['gene_indices'],common) and np.array_equal(z['gene_symbols'],genes.gene_symbol.to_numpy()[common]) and np.array_equal(truth,means[[ci[c] for c in test]][:,common]) and np.array_equal(z['control_mean'],control[common]))
   expected_mse=np.mean((prediction[eligible].astype(np.float64)-truth[eligible])**2,axis=1);saved=metrics[(metrics.model==model)&(metrics.seed==seed)&(metrics.gene_set=='shared')&(~metrics.exclude_direct_targets)].set_index('condition').loc[names[eligible]].MSE.to_numpy()
   check(f'{model}_seed{seed}_MSE_recomputed',np.allclose(expected_mse,saved,atol=1e-12,rtol=1e-10))
   if model in gpu:
    h=pd.read_csv(out/'training_history.csv');bestidx=h.validation_condition_equal_MSE.idxmin();bestepoch=int(h.loc[bestidx,'epoch']);v=np.load(out/'validation_predictions.npz');vmse=float(np.mean((v['predicted'].astype(np.float64)-v['observed'])**2))
    check(f'{model}_seed{seed}_best_validation_selection',bestepoch==result['best_epoch'] and np.isclose(vmse,result['best_validation_MSE'],atol=1e-12) and v['conditions'].astype(str).tolist()==val and np.array_equal(v['observed'],means[[ci[c] for c in val]][:,common]))
    check(f'{model}_seed{seed}_training_and_early_stop',h.epoch.tolist()==list(range(1,len(h)+1)) and len(h)<=20 and bestepoch<=8 and (len(h)>=8 or len(h)-bestepoch>=5) and np.isfinite(h.mean_training_batch_loss).all() and np.isfinite(h.max_gradient_norm).all() and (h.training_rows_seen>=ntrain).all() and (h.training_rows_seen<=ntrain+1).all() and eligible.all() and result['status']=='completed')
    if model!='scGPT_blood':check(f'{model}_seed{seed}_loss_equivalence',json.loads((out/'loss_equivalence.json').read_text())['passed'])
    if model=='scFoundation_GEARS':check(f'{model}_seed{seed}_frozen_cache',json.loads((out/'run_config.json').read_text())['cache_check_max_abs_difference']<=.0002)
   if model=='zero_change':check(f'zero_seed{seed}_definition',np.array_equal(prediction,np.broadcast_to(control[common],prediction.shape)))
   if model=='training_condition_mean':
    tr=[ci[c] for c in train if c!='ctrl'];expected=(means[tr][:,common].astype(np.float64)-control[common]).mean(0)+control[common];check(f'mean_seed{seed}_condition_weighting',np.allclose(prediction,expected,atol=1e-6,rtol=1e-6))
   if model=='Ridge_additive':
    path=pd.read_csv(out/'penalty_validation.csv');check(f'ridge_seed{seed}_validation_selection',float(path.loc[path.validation_MSE.idxmin(),'alpha'])==result['alpha'])
   if model=='direct_single_additive':
    expected=np.array([all(g+'+ctrl' in train for g in c.split('+') if g!='ctrl') for c in test]);check(f'additive_seed{seed}_eligibility',np.array_equal(expected,eligible))
 for key,s in aggregate.groupby(['model','seed','gene_set','exclude_direct_targets','group']):
  m,seed,gs,ex,group=key;d=metrics[(metrics.model==m)&(metrics.seed==seed)&(metrics.gene_set==gs)&(metrics.exclude_direct_targets==ex)]
  if group=='direct_additive_eligible':
   allowed=metrics[(metrics.model=='direct_single_additive')&(metrics.seed==seed)].condition.unique();d=d[d.condition.isin(allowed)]
  elif group!='all':d=d[d.group==group]
  assert len(d)==s.iloc[0].n_conditions and np.isclose(d.MSE.mean(),s.iloc[0].MSE,atol=1e-12)
 check('all_condition_equal_aggregate_MSE',True,rows=len(aggregate))
 result={'passed':True,'checks':checks,'check_count':len(checks),'elapsed_seconds':time.time()-t,'review_mode':'self_checks_only_user_requested'};(R/'logs/stage7_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps({'passed':True,'checks':len(checks),'elapsed_seconds':time.time()-t}))
if __name__=='__main__':main()
