"""Freeze shared Norman inputs, official GEARS simulation splits, and a minimal pilot."""
from pathlib import Path
import os,json,hashlib,zipfile,shutil,ast,pickle,time
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,torch
R=Path(os.environ.get('STAGE6_ROOT','/root/HW2_20261107/solution'));E=R/'data/external/stage6';P=R/'data/processed/stage6';L=R/'logs/stage6';S=R/'data/splits/stage6'
for p in [P,L,S]:p.mkdir(parents=True,exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(2**22),b''):h.update(b)
 return h.hexdigest()
def savejson(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
start=time.time();raw=E/'norman_source.h5ad'
if not raw.exists():
 with zipfile.ZipFile(E/'norman_with_totalcounts.h5ad.zip') as z,z.open(z.namelist()[0]) as src,raw.open('wb') as dst:shutil.copyfileobj(src,dst,2**22)
a=ad.read_h5ad(raw);X=a.X.tocsr().astype(np.float32);genes=a.var.gene_name.astype(str).to_numpy();obs=a.obs.copy();obs['source_condition']=obs.condition.astype(str)
def canonical(s):
 p=sorted(x for x in str(s).split('+') if x!='ctrl');return '+'.join(p) if len(p)>1 else (p[0]+'+ctrl' if p else 'ctrl')
obs['condition']=obs.source_condition.map(canonical)
assert len(set(genes))==len(genes) and obs.index.is_unique and np.isfinite(X.data).all() and (X.data>=0).all()
assert np.isfinite(obs.total_count).all() and (obs.total_count>0).all()
with (E/'gene2go_all.pkl').open('rb') as f:go=pickle.load(f)
pertgenes=sorted({g for c in obs.condition.unique() for g in c.split('+') if g!='ctrl'})
missing_go=[g for g in pertgenes if g not in go or not go[g]]
keep=~obs.condition.map(lambda c:any(g in missing_go for g in c.split('+'))).to_numpy();excluded=obs.loc[~keep,['condition','source_condition']];excluded.to_csv(P/'excluded_cells.csv.gz',compression='gzip')
obs=obs.loc[keep].copy();X=X[keep];pertgenes=[g for g in pertgenes if g not in missing_go]
vocab=json.loads((R/'models/stage6/scGPT_blood/vocab.json').read_text());common=np.array([i for i,g in enumerate(genes) if g in vocab],dtype=np.int64)
foundation_genes=pd.read_csv(R/'vendor/stage6/scFoundation/OS_scRNA_gene_index.19264.tsv',sep='\t').gene_name.astype(str).to_numpy();assert np.array_equal(genes,foundation_genes)
assert all(g in genes for g in pertgenes)
g=pd.DataFrame({'gene_index':np.arange(len(genes)),'gene_symbol':genes,'scgpt_token_id':[vocab.get(v,-1) for v in genes],'shared_output':[v in vocab for v in genes],'scfoundation_index':np.arange(len(genes))});g.to_csv(P/'genes.csv',index=False)
# Execute the unmodified upstream class/function definitions without optional import side effects.
namespace={'np':np};uf=R/'vendor/stage6/GEARS/gears/utils.py';df=R/'vendor/stage6/GEARS/gears/data_utils.py'
for path,names in [(uf,{'parse_single_pert','parse_combo_pert','parse_any_pert'}),(df,{'DataSplitter'})]:
 tree=ast.parse(path.read_text());nodes=[n for n in tree.body if isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name in names];exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),namespace)
conditions=[];indices={'common':common};splitstats=[]
for seed in [1,2,3]:
 small=ad.AnnData(X=sp.csr_matrix((len(obs),1)),obs=obs[['condition']].copy());splitter=namespace['DataSplitter'](small,split_type='simulation');out,sub=splitter.split_data(seed=seed,train_gene_set_size=.75,combo_seen2_train_frac=.75)
 splits=out.obs['split'].astype(str);assert set(splits)=={'train','val','test'};obs[f'split_seed{seed}']=splits.to_numpy()
 groups={s:sorted(obs.loc[splits.to_numpy()==s,'condition'].unique().tolist()) for s in ['train','val','test']};assert not(set(groups['train'])&set(groups['test']) or set(groups['val'])&set(groups['test']) or set(groups['train'])&set(groups['val']))
 train_genes={x for c in groups['train'] for x in c.split('+') if x!='ctrl'}
 official={c:k for ss in sub.values() for k,cs in ss.items() for c in cs}
 for s,cs in groups.items():
  for c in cs:
   pg=[x for x in c.split('+') if x!='ctrl'];seen=sum(x in train_genes for x in pg);conditions.append({'seed':seed,'split':s,'condition':c,'n_cells':int((obs.condition==c).sum()),'n_perturbed_genes':len(pg),'seen_in_final_training':seen,'evaluation_group':'ctrl' if not pg else (f'seen{seen}' if len(pg)==2 else ('seen_single' if seen else 'unseen_single')),'official_subgroup':official.get(c,'train')})
 train=X[splits.to_numpy()=='train'][:,common];mean=np.asarray(train.mean(axis=0)).ravel();var=np.asarray(train.power(2).mean(axis=0)).ravel()-mean**2;ix=np.argsort(-var,kind='stable')[:1000];indices[f'hvg1000_seed{seed}']=common[ix]
 savejson(S/f'seed{seed}.json',{'seed':seed,'split':'GEARS simulation','conditions':groups,'official_subgroups':sub,'train_gene_set_size':.75,'combo_seen2_train_frac':.75})
 splitstats.append({'seed':seed,**{k:len(v) for k,v in groups.items()}})
np.savez_compressed(P/'evaluation_indices.npz',**indices);obs.index.name='cell_id';obs.to_csv(P/'cells.csv.gz',compression='gzip');ct=pd.DataFrame(conditions);ct.to_csv(S/'condition_splits.csv',index=False)
sp.save_npz(P/'expression_log1p.npz',X,compressed=True);np.save(P/'total_counts.npy',obs.total_count.to_numpy(dtype=np.float64));savejson(P/'perturbation_genes.json',pertgenes)
ctrl=np.flatnonzero(obs.condition.to_numpy()=='ctrl');assert np.all(obs.iloc[ctrl].split_seed1=='train');control_mean=np.asarray(X[ctrl].mean(axis=0)).ravel();np.save(P/'control_mean.npy',control_mean)
# Paired control inputs use only training controls. The same pairs feed every model.
rng=np.random.default_rng(20261107);pilot=[]
for split in ['train','val','test']:
 rows=ct[(ct.seed==1)&(ct.split==split)&(ct.condition!='ctrl')]
 selected=[]
 if split=='train':
  for n in [1,2]:selected+=rows[rows.n_perturbed_genes==n].sample(n=min(2,sum(rows.n_perturbed_genes==n)),random_state=20261107).condition.tolist()
 else:
  for group in ['seen0','seen1','seen2','unseen_single']:
   eligible=rows[rows.evaluation_group==group]
   if len(eligible):selected.append(eligible.iloc[0].condition)
 for c in selected:
  ids=rng.choice(np.flatnonzero(obs.condition.to_numpy()==c),2,replace=False)
  for i in ids:pilot.append({'split':split,'condition':c,'target_row':int(i),'target_cell':str(obs.index[i]),'control_row':int(rng.choice(ctrl))})
pd.DataFrame(pilot).to_csv(S/'pilot_pairs.csv',index=False)
# Native GEARS co-expression definition: training controls and training single perturbations.
torch.set_num_threads(8);device='cuda:0';mask=(obs.split_seed1=='train') & obs.condition.str.contains('ctrl');A=torch.tensor(X[mask.to_numpy()].toarray(),dtype=torch.float64,device=device);A-=A.mean(0,keepdim=True);norm=torch.sqrt((A*A).sum(0));cov=A.T@A;corr=torch.nan_to_num(cov/(norm[:,None]*norm[None,:]),nan=0.0).clamp(-1,1).abs();del A,cov
vals,src=torch.topk(corr,21,dim=1);target=torch.arange(len(genes),device=device)[:,None].expand_as(src);valid=vals>.4;coedge=torch.stack([src[valid],target[valid]]).cpu().numpy();cow=vals[valid].float().cpu().numpy();del corr,vals,src,target;torch.cuda.empty_cache()
# GO Jaccard edges >0.1, top21 incoming per target, restricted to the declared perturbation universe.
goedge=[];gow=[]
for j,t in enumerate(pertgenes):
 scores=[]
 for i,s in enumerate(pertgenes):
  score=len(go[t]&go[s])/len(go[t]|go[s])
  if score>.1:scores.append((i,score))
 for i,v in sorted(scores,key=lambda z:(-z[1],z[0]))[:21]:goedge.append((i,j));gow.append(v)
np.savez_compressed(P/'gears_graphs_seed1.npz',coexpression_edge_index=coedge,coexpression_edge_weight=cow,go_edge_index=np.asarray(goedge,dtype=np.int64).T,go_edge_weight=np.asarray(gow,dtype=np.float32))
report={'source_shape':[int(a.n_obs),int(a.n_vars)],'retained_shape':list(X.shape),'conditions':int(obs.condition.nunique()),'controls':len(ctrl),'perturbation_genes':len(pertgenes),'missing_GO_genes':missing_go,'excluded_cells':len(excluded),'common_output_genes':len(common),'foundation_exact_gene_order':True,'expression_scale':'publisher log1p normalized expression; no second normalization','total_counts':'publisher obs.total_count, preserved separately; appended only for scFoundation adapter','training_HVG_per_seed':1000,'splits':splitstats,'pilot_pairs':len(pilot),'graph_cells':int(mask.sum()),'coexpression_edges':coedge.shape[1],'go_edges':len(goedge),'coexpression_precision':'float64 GPU; top21 absolute Pearson >0.4; train controls+single only','GO_scope':'custom graph on all declared target genes; shared across GEARS and scFoundation','official_split_code_sha256':sha(df),'upstream_all_data_DE_annotations_used':False,'elapsed_seconds':time.time()-start,'passed':True}
savejson(L/'data_preparation.json',report);print(json.dumps(report),flush=True)
