"""Prepare shared condition means, training-only graphs and fixed reference inputs."""
from pathlib import Path
import json,os,time
import numpy as np,pandas as pd,scipy.sparse as sp,torch
R=Path(os.environ.get('STAGE6_ROOT','/root/HW2_20261107/solution'));P=R/'data/processed/stage6';Q=R/'data/processed/stage7';C=Path('/root/HW2_20261107/cache/stage7')
Q.mkdir(parents=True,exist_ok=True);C.mkdir(parents=True,exist_ok=True)
t=time.time();torch.set_num_threads(8)
X=sp.load_npz(P/'expression_log1p.npz');cells=pd.read_csv(P/'cells.csv.gz');conditions=sorted(cells.condition.unique());code=np.array([conditions.index(c) for c in cells.condition],dtype=np.int64)
dense=np.lib.format.open_memmap(C/'expression_dense.npy',mode='w+',dtype='float32',shape=X.shape)
for i in range(0,X.shape[0],2048):dense[i:i+2048]=X[i:i+2048].toarray()
dense.flush();del dense
means=np.stack([np.asarray(X[code==i].mean(0)).ravel() for i in range(len(conditions))]);np.save(Q/'condition_means.npy',means);np.save(Q/'cell_condition_codes.npy',code)
(Q/'conditions.json').write_text(json.dumps(conditions)+'\n')
ctrl=np.flatnonzero(cells.condition.to_numpy()=='ctrl');refs=np.random.default_rng(20261107).choice(ctrl,64,replace=False);np.save(Q/'reference_control_rows.npy',refs)
cells.iloc[refs].to_csv(Q/'reference_controls.csv',index=False)
np.save(Q/'reference_expression.npy',X[refs].toarray());np.save(Q/'reference_total_counts.npy',np.load(P/'total_counts.npy')[refs])
nonzero=np.stack([np.asarray(X[code==i].sum(0)).ravel()>0 for i in range(len(conditions))]);np.save(Q/'condition_nonzero.npy',nonzero)
base=np.load(P/'gears_graphs_seed1.npz');np.savez_compressed(Q/'gears_graphs_seed1.npz',**dict(base))
(Q/'common_prepared.json').write_text(json.dumps({'cells':len(cells),'conditions':len(conditions),'reference_cells':64,'all_training_targets_retained':True,'control_pool_used_for_training_and_inference':True})+'\n')
for seed in [2,3]:
 mask=(cells[f'split_seed{seed}']=='train')&cells.condition.str.contains('ctrl');device='cuda:3'
 A=torch.tensor(X[mask.to_numpy()].toarray(),dtype=torch.float64,device=device);A-=A.mean(0,keepdim=True);norm=torch.sqrt((A*A).sum(0));cov=A.T@A;corr=torch.nan_to_num(cov/(norm[:,None]*norm[None,:]),nan=0.0).clamp(-1,1).abs();del A,cov
 vals,src=torch.topk(corr,21,dim=1);target=torch.arange(X.shape[1],device=device)[:,None].expand_as(src);valid=vals>.4
 edge=torch.stack([src[valid],target[valid]]).cpu().numpy();weight=vals[valid].float().cpu().numpy()
 np.savez_compressed(Q/f'gears_graphs_seed{seed}.npz',coexpression_edge_index=edge,coexpression_edge_weight=weight,go_edge_index=base['go_edge_index'],go_edge_weight=base['go_edge_weight'])
 del corr,vals,src,target;torch.cuda.empty_cache();print('graph',seed,int(mask.sum()),edge.shape[1],flush=True)
report={'passed':True,'elapsed_seconds':time.time()-t,'dense_cache':str(C/'expression_dense.npy'),'dense_cache_role':'temporary exact dense copy of archived sparse matrix; removable after jobs','reference_pool_seed':20261107,'reference_pool_size':64}
(Q/'prepared.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
