"""Full cell-level training and condition-level validation using stage6 upstream models."""
from pathlib import Path
import os,sys,json,time,argparse,traceback,types
import numpy as np,pandas as pd,torch
from torch import nn
import stage6_pilot as pilot
R=pilot.R;P=R/'data/processed/stage6';Q=R/'data/processed/stage7';V=pilot.V
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n')
def atomic_torch(p,x):
 tmp=p.with_suffix('.tmp');torch.save(x,tmp);tmp.replace(p)
class CachedEncoder(nn.Module):
 def __init__(self,cache):super().__init__();self.cache=cache;self.indices=None
 def forward(self,x):
  assert self.indices is not None and len(self.indices)==len(x)
  return self.cache[self.indices]

def build(name,seed,device):
 genes=pd.read_csv(P/'genes.csv');n=len(genes);pg=json.loads((P/'perturbation_genes.json').read_text());meta={}
 if name=='scGPT_blood':return pilot.build_scgpt(json.loads((R/'models/stage6/scGPT_blood/vocab.json').read_text()))
 graphs=np.load(Q/f'gears_graphs_seed{seed}.npz');cfg={'num_genes':n,'num_perts':len(pg),'hidden_size':64 if name=='GEARS' else 512,'uncertainty':False,'num_go_gnn_layers':1,'num_gene_gnn_layers':1,'decoder_hidden_size':16,'no_perturb':False,'device':device,'G_coexpress':torch.tensor(graphs['coexpression_edge_index']),'G_coexpress_weight':torch.tensor(graphs['coexpression_edge_weight']),'G_go':torch.tensor(graphs['go_edge_index']),'G_go_weight':torch.tensor(graphs['go_edge_weight'])}
 path=V/'GEARS/gears/model.py'
 if name=='scFoundation_GEARS':
  pilot.foundation_encoder();path=V/'scFoundation/GEARS/gears/model.py';cfg.update({'cell_fitness_pred':False,'model_type':'maeautobin','mode':'v1','bin_set':'autobin_resolution_append','highres':0,'load_path':str(R/'models/stage6/scFoundation/models.ckpt'),'finetune_method':'frozen'})
 model=pilot.module(path,'formal_'+name).GEARS_Model(cfg)
 if name=='scFoundation_GEARS':
  for p in model.singlecell_model.parameters():p.requires_grad=False
  model.singlecell_model.eval();meta['frozen_parameters']=sum(p.numel() for p in model.singlecell_model.parameters())
 return model,meta

def main(args):
 name,seed,device=args.model,args.seed,args.device;out=R/f'results/stage7/{name}/seed{seed}';out.mkdir(parents=True,exist_ok=True)
 if (out/'completed.json').exists():print('Already completed',out,flush=True);return
 torch.set_num_threads(6);torch.cuda.set_device(device);torch.manual_seed(20261107+seed);np.random.seed(20261107+seed);torch.set_float32_matmul_precision('high');torch.backends.cuda.matmul.allow_tf32=True;torch.backends.cudnn.allow_tf32=True
 start=time.time();max_epochs=int(json.loads((R/'config/stage7_config.json').read_text())['max_epochs']);torch.cuda.reset_peak_memory_stats();rng=np.random.default_rng(20261107+seed)
 X=np.load('/root/HW2_20261107/cache/stage7/expression_dense.npy',mmap_mode='r');cells=pd.read_csv(P/'cells.csv.gz');genes=pd.read_csv(P/'genes.csv');n=len(genes);common=np.load(P/'evaluation_indices.npz')['common'];conds=json.loads((Q/'conditions.json').read_text());codes=np.load(Q/'cell_condition_codes.npy');means=np.load(Q/'condition_means.npy');nonzero=np.load(Q/'condition_nonzero.npy');split=json.loads((R/f'data/splits/stage6/seed{seed}.json').read_text())['conditions'];ci={c:i for i,c in enumerate(conds)}
 train_rows=np.flatnonzero(cells[f'split_seed{seed}'].to_numpy()=='train');refs=torch.tensor(np.load(Q/'reference_expression.npy'),device=device);totals=torch.tensor(np.load(Q/'reference_total_counts.npy'),device=device,dtype=torch.float32);control=torch.tensor(np.load(P/'control_mean.npy'),device=device)
 pg=json.loads((P/'perturbation_genes.json').read_text());pmap={p:i for i,p in enumerate(pg)};gmap={g:i for i,g in enumerate(genes.gene_symbol)}
 condition_perts={c:[pmap[g] for g in c.split('+') if g!='ctrl'] or [-1] for c in conds};condition_gene={c:[gmap[g] for g in c.split('+') if g!='ctrl'] for c in conds}
 filters={c:np.flatnonzero(nonzero[ci[c]]).tolist() for c in split['train']}
 model,meta=build(name,seed,device);model.to(device)
 if name=='scFoundation_GEARS':
  cachefile=Q/'scFoundation_reference_embeddings.npy';lock=Q/'scFoundation_cache.lock'
  if not cachefile.exists():
   try:fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd);owner=True
   except FileExistsError:owner=False
   if owner:
    tmp=Q/'scFoundation_reference_embeddings.partial.npy';cache=np.lib.format.open_memmap(tmp,mode='w+',dtype='float32',shape=(64,n,512));model.singlecell_model.eval()
    with torch.no_grad():
     for k in range(0,64,2):cache[k:k+2]=model.singlecell_model(torch.cat([refs[k:k+2],totals[k:k+2,None]],1)).cpu().numpy()
    cache.flush();del cache;tmp.replace(cachefile);save(Q/'scFoundation_cache.json',{'shape':[64,n,512],'dtype':'float32','frozen_eval':True,'input_rows':'reference_control_rows.npy','TF32':True});lock.unlink()
   else:
    while not cachefile.exists():time.sleep(3)
  cache=torch.tensor(np.load(cachefile),device=device)
  with torch.no_grad():
   actual=model.singlecell_model(torch.cat([refs[:2],totals[:2,None]],1));err=float((actual-cache[:2]).abs().max());assert torch.allclose(actual,cache[:2],atol=2e-4,rtol=2e-4),err
  meta['cache_check_max_abs_difference']=err;model.singlecell_model=CachedEncoder(cache)
  torch.cuda.empty_cache()
 tokens=torch.tensor(genes.scgpt_token_id.to_numpy()[common],device=device) if name=='scGPT_blood' else None
 common_position={int(g):i for i,g in enumerate(common)}
 ns={'torch':torch,'np':np};pilot.definitions(V/'GEARS/gears/utils.py',{'loss_fct'},ns)
 loss_fct=ns['loss_fct'];lr=1e-4 if name=='scGPT_blood' else 1e-3;wd=0 if name=='scGPT_blood' else 5e-4
 optim=torch.optim.Adam([p for p in model.parameters() if p.requires_grad],lr=lr,weight_decay=wd)
 scheduler=torch.optim.lr_scheduler.StepLR(optim,step_size=1,gamma=.9 if name=='scGPT_blood' else .5)
 batch={'GEARS':64,'scGPT_blood':64,'scFoundation_GEARS':24}[name];evalbatch={'GEARS':32,'scGPT_blood':16,'scFoundation_GEARS':12}[name]
 nonzero[ci['ctrl']]=True
 loss_masks=torch.tensor(nonzero,device=device,dtype=torch.float32)
 loss_counts=loss_masks.sum(1)
 vectorized_loss_verified=False
 def fast_gears_loss(pr,yy,cs):
  nonlocal vectorized_loss_verified
  ids=np.array([ci[c] for c in cs]);_,inv,counts=np.unique(ids,return_inverse=True,return_counts=True)
  mask=loss_masks[ids];weights=1/(torch.tensor(counts[inv],device=device)*loss_counts[ids]*len(counts))
  power=((pr-yy)**4*mask).sum(1)
  direction=((torch.sign(yy-control)-torch.sign(pr-control))**2*mask).sum(1)
  value=((power+.001*direction)*weights).sum()
  if not vectorized_loss_verified:
   reference=loss_fct(pr,yy,cs,ctrl=control,dict_filter=filters)
   assert torch.allclose(value,reference,atol=2e-6,rtol=2e-5),(float(value),float(reference))
   a=torch.autograd.grad(value,pr,retain_graph=True)[0];b=torch.autograd.grad(reference,pr,retain_graph=True)[0]
   assert torch.allclose(a,b,atol=1e-7,rtol=2e-5)
   save(out/'loss_equivalence.json',{'passed':True,'value_abs_difference':float((value-reference).abs().detach()),'gradient_max_abs_difference':float((a-b).abs().max()),'batch_cells':len(cs)})
   vectorized_loss_verified=True
  return value
 meta.update({'train_batch_size':batch,'evaluation_batch_size':evalbatch,'train_cells':len(train_rows),'learning_rate':lr,'weight_decay':wd,'scheduler_gamma':.9 if name=='scGPT_blood' else .5,'gradient_clip':'norm1' if name=='scGPT_blood' else 'value1','reference_pool_size':64,'TF32':True,'initialization':'fresh public pretraining or random, no pilot checkpoint','train_graph_seed':seed,'max_epochs':max_epochs,'resume_prior_seconds':args.prior_seconds})
 save(out/'run_config.json',meta)
 def state():return {k:v.detach().cpu() for k,v in model.state_dict().items() if not k.startswith('singlecell_model.')}
 def forward(cs,refidx,training=False):
  xx=refs[refidx];b=len(cs)
  if name=='scGPT_blood':
   idx=np.arange(len(common))
   if training:
    req=sorted({common_position[g] for c in cs for g in condition_gene[c]});rest=np.setdiff1d(idx,req);idx=np.sort(np.r_[req,rng.choice(rest,1536-len(req),replace=False)]).astype(np.int64)
   sub=common[idx];flags=torch.zeros((b,len(idx)),dtype=torch.long,device=device);position={int(g):i for i,g in enumerate(sub)}
   for i,c in enumerate(cs):
    for g in condition_gene[c]:flags[i,position[g]]=1
   with torch.autocast('cuda',dtype=torch.bfloat16):pr=model(tokens[idx].expand(b,-1),xx[:,sub],flags,src_key_padding_mask=torch.zeros((b,len(idx)),dtype=torch.bool,device=device))['mlm_output']
   return pr.float(),sub
  if name=='scFoundation_GEARS':model.singlecell_model.indices=refidx;xx=torch.cat([xx,totals[refidx,None]],1)
  data=types.SimpleNamespace(x=xx.reshape(-1,1),pert_idx=[condition_perts[c] for c in cs],batch=torch.arange(b,device=device).repeat_interleave(xx.shape[1]))
  return model(data),np.arange(n)
 def predict(names):
  model.eval();pred=[];t=time.time()
  with torch.no_grad():
   for c in names:
    total=torch.zeros(len(common),device=device,dtype=torch.float64)
    for k in range(0,64,evalbatch):
     ids=torch.arange(k,min(k+evalbatch,64),device=device);pr,idx=forward([c]*len(ids),ids,False)
     if name!='scGPT_blood':pr=pr[:,common]
     assert torch.isfinite(pr).all();total+=pr.double().sum(0)
    pred.append((total/64).float().cpu().numpy())
  return np.stack(pred),time.time()-t
 history=[];best=float('inf');bestepoch=0;begin=1;latest=out/'latest.pt'
 if latest.exists():
  ck=torch.load(latest,map_location=device,weights_only=False);model.load_state_dict(ck['model'],strict=False);optim.load_state_dict(ck['optimizer']);scheduler.load_state_dict(ck['scheduler']);history=ck['history'];best=ck['best'];bestepoch=ck['bestepoch'];begin=ck['epoch']+1;rng.bit_generator.state=ck['rng'];torch.set_rng_state(ck['torch_rng'].cpu());torch.cuda.set_rng_state(ck['cuda_rng'].cpu(),device)
 valnames=[c for c in split['val'] if c!='ctrl'];truth=means[[ci[c] for c in valnames]][:,common]
 for epoch in range(begin,max_epochs+1):
  t=time.time();epoch_lr=optim.param_groups[0]['lr'];model.train();order=rng.permutation(train_rows);losses=[];maxgrad=0;count=0
  for startrow in range(0,len(order),batch):
   row=order[startrow:startrow+batch]
   if len(row)==1: # BatchNorm needs two rows. One repeated training row is logged.
    row=np.r_[row,order[0]]
   cs=[conds[c] for c in codes[row]];ri=torch.tensor(rng.integers(0,64,len(row)),device=device);optim.zero_grad(set_to_none=True);pr,idx=forward(cs,ri,True);yy=torch.tensor(np.asarray(X[row][:,idx]),device=device)
   loss=((pr-yy)**2).mean() if name=='scGPT_blood' else fast_gears_loss(pr,yy,cs)
   if not torch.isfinite(loss):raise RuntimeError(f'nonfinite loss epoch {epoch} batch {startrow}')
   loss.backward();gn=torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad],1. if name=='scGPT_blood' else float('inf'))
   if name!='scGPT_blood':torch.nn.utils.clip_grad_value_([p for p in model.parameters() if p.grad is not None],1.)
   if not torch.isfinite(gn):raise RuntimeError('nonfinite gradient')
   optim.step();losses.append(float(loss.detach()));maxgrad=max(maxgrad,float(gn));count+=len(row)
   if startrow==0 or (startrow//batch+1)%250==0:
    save(out/'progress.json',{'model':name,'seed':seed,'epoch':epoch,'phase':'training','rows_processed':count,'training_rows':len(order),'epoch_elapsed_seconds':time.time()-t,'time_unix':time.time()})
  trainsec=time.time()-t;prediction,valsec=predict(valnames);val=float(np.mean((prediction.astype(np.float64)-truth)**2));improved=val<best
  if improved:
   best=val;bestepoch=epoch;atomic_torch(out/'best.pt',state());np.savez_compressed(out/'validation_predictions.npz',predicted=prediction,observed=truth,conditions=np.array(valnames),gene_indices=common)
  item={'epoch':epoch,'learning_rate':epoch_lr,'mean_training_batch_loss':float(np.mean(losses)),'validation_condition_equal_MSE':val,'training_seconds':trainsec,'validation_seconds':valsec,'best_epoch':bestepoch,'best_validation_MSE':best,'training_rows_seen':count,'max_gradient_norm':maxgrad}
  history.append(item);pd.DataFrame(history).to_csv(out/'training_history.csv',index=False);print(name,seed,item,flush=True)
  scheduler.step()
  atomic_torch(latest,{'model':state(),'optimizer':optim.state_dict(),'scheduler':scheduler.state_dict(),'history':history,'best':best,'bestepoch':bestepoch,'epoch':epoch,'rng':rng.bit_generator.state,'torch_rng':torch.get_rng_state(),'cuda_rng':torch.cuda.get_rng_state(device)})
  save(out/'progress.json',{'model':name,'seed':seed,'epoch':epoch,'phase':'epoch_completed',**item,'time_unix':time.time()})
  if epoch-bestepoch>=5:break
 model.load_state_dict(torch.load(out/'best.pt',map_location=device,weights_only=True),strict=False)
 testnames=[c for c in split['test'] if c!='ctrl'];prediction,testsec=predict(testnames);truth=means[[ci[c] for c in testnames]][:,common]
 np.savez_compressed(out/'test_predictions.npz',predicted=prediction,observed=truth,conditions=np.array(testnames),gene_indices=common,gene_symbols=genes.gene_symbol.to_numpy()[common].astype(str),control_mean=control[common].cpu().numpy())
 summary={'status':'completed','model':name,'seed':seed,'epochs':history[-1]['epoch'],'best_epoch':bestepoch,'best_validation_MSE':best,'test_condition_equal_MSE':float(np.mean((prediction.astype(np.float64)-truth)**2)),'test_conditions':len(testnames),'test_seconds':testsec,'total_seconds':args.prior_seconds+time.time()-start,'max_epochs_requested':max_epochs,'selection_max_epochs':max_epochs,'prechange_extra_epochs':max(0,len(history)-max_epochs),'peak_allocated_GiB':torch.cuda.max_memory_allocated()/2**30,'peak_reserved_GiB':torch.cuda.max_memory_reserved()/2**30,'all_predictions_finite':bool(np.isfinite(prediction).all()),'train_cells':len(train_rows),'training_reference_pool_size':64,'frozen_cache':name=='scFoundation_GEARS'}
 save(out/'completed.json',summary);print(summary,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--model',required=True);p.add_argument('--seed',type=int,required=True);p.add_argument('--device',required=True);p.add_argument('--prior-seconds',type=float,default=0);a=p.parse_args()
 try:main(a)
 except Exception as e:
  out=R/f'results/stage7/{a.model}/seed{a.seed}';out.mkdir(parents=True,exist_ok=True);save(out/'error.json',{'error':str(e),'traceback':traceback.format_exc(),'time_unix':time.time()});raise
