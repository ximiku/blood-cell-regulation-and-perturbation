"""Minimal real GPU training using pinned upstream architectures and shared pairs.
Runs 4 optimizer steps, then held-out predictions; not a model benchmark.
"""
from pathlib import Path
import os,sys,json,time,ast,types,importlib.util,importlib,logging,re,typing,argparse,traceback,hashlib
import numpy as np,pandas as pd,scipy.sparse as sp,torch
from torch import nn
from types import SimpleNamespace
R=Path(os.environ.get('STAGE6_ROOT','/root/HW2_20261107/solution'));P=R/'data/processed/stage6';V=R/'vendor/stage6';D=R/'results/stage6'
def js(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n')
def module(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def package(name,path):
 m=types.ModuleType(name);m.__path__=[str(path)];sys.modules[name]=m;return m

def definitions(path,names,ns,torch_load_cpu=False):
 tree=ast.parse(path.read_text());nodes=[x for x in tree.body if isinstance(x,(ast.ClassDef,ast.FunctionDef)) and x.name in names]
 if torch_load_cpu:
  for n in nodes:
   for c in ast.walk(n):
    if isinstance(c,ast.Call) and isinstance(c.func,ast.Attribute) and isinstance(c.func.value,ast.Name) and c.func.value.id=='torch' and c.func.attr=='load':c.keywords.extend([ast.keyword(arg='map_location',value=ast.Constant('cpu')),ast.keyword(arg='weights_only',value=ast.Constant(False))])
 exec(compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),str(path),'exec'),ns)

def build_scgpt(vocab):
 pkg=package('scgpt',V/'scGPT/scgpt');pkg.logger=logging.getLogger('scGPT');package('scgpt.model',V/'scGPT/scgpt/model');utils=package('scgpt.utils',V/'scGPT/scgpt/utils')
 utils.__dict__.update({**vars(typing),'np':np,'torch':torch,'re':re,'logger':pkg.logger})
 definitions(V/'scGPT/scgpt/utils/util.py',{'map_raw_id_to_vocab_id','load_pretrained','_rename_key'},utils.__dict__)
 M=importlib.import_module('scgpt.model.generation_model').TransformerGenerator
 cfg=json.loads((R/'models/stage6/scGPT_blood/args.json').read_text());m=M(ntoken=len(vocab),d_model=cfg['embsize'],nhead=cfg['nheads'],d_hid=cfg['d_hid'],nlayers=cfg['nlayers'],nlayers_cls=cfg['n_layers_cls'],n_cls=1,vocab=vocab,dropout=0,pad_token='<pad>',pad_value=0,pert_pad_id=0,n_input_bins=0,use_fast_transformer=False)
 weights=torch.load(R/'models/stage6/scGPT_blood/best_model.pt',map_location='cpu',weights_only=True);prefix=['encoder.','value_encoder.','transformer_encoder.'];utils.load_pretrained(m,weights,prefix=prefix,verbose=False)
 renamed={utils._rename_key(k,{r'self_attn\._impl\.Wqkv\.':'self_attn.in_proj_',r'self_attn\.Wqkv\.':'self_attn.in_proj_',r'self_attn\._impl\.out_proj\.':'self_attn.out_proj.'}):v for k,v in weights.items()}
 expected={k:v for k,v in m.state_dict().items() if any(k.startswith(p) for p in prefix)};assert all(k in renamed and torch.equal(v,renamed[k]) for k,v in expected.items()),'Incomplete pretrained encoder loading'
 return m,{'encoder_tensors_loaded':len(expected),'encoder_parameters_loaded':sum(v.numel() for v in expected.values()),'backend':'native PyTorch SDPA; upstream key translation','train_sequence_length':1536,'prediction_sequence':'all shared output genes','pretraining_config':cfg}

def foundation_encoder():
 package('modules',V/'scFoundation/GEARS/modules');pf=importlib.import_module('modules.performer_module');tr=importlib.import_module('modules.transformer');ma=importlib.import_module('modules.mae_autobin')
 ns={'torch':torch,'nn':nn,'np':np,'math':__import__('math'),'PerformerModule':pf.PerformerModule,'pytorchTransformerModule':tr.pytorchTransformerModule,'MaeAutobin':ma.MaeAutobin}
 definitions(V/'scFoundation/model/load.py',{'convertconfig'},ns)
 definitions(V/'scFoundation/GEARS/modules/encoders.py',{'next_16x','gatherData','getEncoerDecoderData','select_module','choose_binset','MAEAutobinencoder'},ns,torch_load_cpu=True)
 mod=types.ModuleType('modules.encoders');mod.MAEAutobinencoder=ns['MAEAutobinencoder'];sys.modules['modules.encoders']=mod

def main(name,device):
 start=time.time();out=D/name;out.mkdir(parents=True,exist_ok=True);torch.set_num_threads(6);torch.manual_seed(20261107);np.random.seed(20261107);torch.cuda.set_device(device);torch.cuda.reset_peak_memory_stats(device)
 X=sp.load_npz(P/'expression_log1p.npz');genes=pd.read_csv(P/'genes.csv');cells=pd.read_csv(P/'cells.csv.gz');pairs=pd.read_csv(R/'data/splits/stage6/pilot_pairs.csv');common=np.load(P/'evaluation_indices.npz')['common'];totals=np.load(P/'total_counts.npy');pg=json.loads((P/'perturbation_genes.json').read_text());pmap={p:i for i,p in enumerate(pg)};gmap={g:i for i,g in enumerate(genes.gene_symbol)};control=torch.tensor(np.load(P/'control_mean.npy'),device=device)
 n=X.shape[1];graphs=np.load(P/'gears_graphs_seed1.npz');meta={}
 if name=='scGPT_blood':
  vocab=json.loads((R/'models/stage6/scGPT_blood/vocab.json').read_text());model,meta=build_scgpt(vocab);gene_tokens=torch.tensor([vocab[genes.iloc[i].gene_symbol] for i in common],device=device);output_idx=common
 else:
  cfg={'num_genes':n,'num_perts':len(pg),'hidden_size':64 if name=='GEARS' else 512,'uncertainty':False,'num_go_gnn_layers':1,'num_gene_gnn_layers':1,'decoder_hidden_size':16,'no_perturb':False,'device':device,'G_coexpress':torch.tensor(graphs['coexpression_edge_index']),'G_coexpress_weight':torch.tensor(graphs['coexpression_edge_weight']),'G_go':torch.tensor(graphs['go_edge_index']),'G_go_weight':torch.tensor(graphs['go_edge_weight'])}
  path=V/'GEARS/gears/model.py'
  if name=='scFoundation_GEARS':
   foundation_encoder();path=V/'scFoundation/GEARS/gears/model.py';cfg.update({'cell_fitness_pred':False,'model_type':'maeautobin','mode':'v1','bin_set':'autobin_resolution_append','highres':0,'load_path':str(R/'models/stage6/scFoundation/models.ckpt'),'finetune_method':'frozen'})
  model=module(path,'hw2_'+name).GEARS_Model(cfg);output_idx=np.arange(n)
  if name=='scFoundation_GEARS':
   for p in model.singlecell_model.parameters():p.requires_grad=False
   model.singlecell_model.eval();meta['frozen_parameters']=sum(p.numel() for p in model.singlecell_model.parameters());meta['loading_adaptation']='Official encoder class, CPU checkpoint load with explicit weights_only=False; unrelated optional imports omitted'
  ns={'torch':torch,'np':np};definitions(V/'GEARS/gears/utils.py',{'loss_fct'},ns);loss_fct=ns['loss_fct']
 model.to(device);optim=torch.optim.Adam([p for p in model.parameters() if p.requires_grad],lr=1e-4 if name=='scGPT_blood' else 1e-3)
 witness=next(p for p in model.parameters() if p.requires_grad);before=witness.detach().clone();loaded_seconds=time.time()-start;rng=np.random.default_rng(20261107);train_rows=pairs[pairs.split=='train'].reset_index(drop=True);history=[]
 filters={c:np.flatnonzero(np.asarray(X[cells.condition.to_numpy()==c].sum(axis=0)).ravel()>0).tolist() for c in train_rows.condition.unique()}
 def forward(rows,training=False):
  cr=rows.control_row.to_numpy();tr=rows.target_row.to_numpy();xx=torch.tensor(X[cr].toarray(),device=device);yy=torch.tensor(X[tr].toarray(),device=device);conditions=rows.condition.tolist();batch=len(rows)
  if name=='scGPT_blood':
   idx=np.arange(len(common))
   if training and len(idx)>1536:
    required={int(np.flatnonzero(common==gmap[g])[0]) for c in conditions for g in c.split('+') if g!='ctrl'};rest=np.setdiff1d(idx,list(required));idx=np.sort(np.r_[list(required),rng.choice(rest,1536-len(required),replace=False)]).astype(int)
   sub=common[idx];flags=torch.zeros((batch,len(idx)),dtype=torch.long,device=device)
   for i,c in enumerate(conditions):
    for g in c.split('+'):
     if g!='ctrl':flags[i,torch.tensor(sub==gmap[g],device=device)]=1
   vals=xx[:,sub];tokens=gene_tokens[idx].expand(batch,-1);mask=torch.zeros_like(tokens,dtype=torch.bool)
   with torch.autocast('cuda',dtype=torch.bfloat16):pred=model(tokens,vals,flags,src_key_padding_mask=mask)['mlm_output']
   return pred.float(),yy[:,sub],sub
  inds=[[pmap[g] for g in c.split('+') if g!='ctrl'] or [-1] for c in conditions]
  if name=='scFoundation_GEARS':xx=torch.cat([xx,torch.tensor(totals[cr,None],device=device,dtype=xx.dtype)],dim=1)
  data=SimpleNamespace(x=xx.reshape(-1,1),pert_idx=inds,batch=torch.arange(batch,device=device).repeat_interleave(xx.shape[1]))
  pred=model(data);return pred,yy,np.arange(n)
 model.train()
 if name=='scFoundation_GEARS':model.singlecell_model.eval()
 for step in range(4):
  rows=train_rows.iloc[step*2:(step+1)*2];assert len(rows)==2
  torch.cuda.synchronize();t=time.time();optim.zero_grad(set_to_none=True);pred,yy,idx=forward(rows,True)
  loss=((pred-yy)**2).mean() if name=='scGPT_blood' else loss_fct(pred,yy,rows.condition.tolist(),ctrl=control,dict_filter=filters)
  assert torch.isfinite(loss);loss.backward();grad=float(torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad],1.0));assert np.isfinite(grad) and grad>0;optim.step();torch.cuda.synchronize()
  history.append({'step':step+1,'loss':float(loss.detach()),'gradient_norm_before_clip':grad,'seconds':time.time()-t,'batch_cells':2,'training_genes':len(idx)});print(name,history[-1],flush=True)
 assert not torch.equal(before,witness),'No trainable parameter update'
 frozen_ok=True
 if name=='scFoundation_GEARS':frozen_ok=all(p.grad is None and not p.requires_grad for p in model.singlecell_model.parameters());assert frozen_ok
 model.eval();predictions=[];target=[];details=[];t=time.time()
 with torch.no_grad():
  for split in ['val','test']:
   rows_all=pairs[pairs.split==split]
   for k in range(0,len(rows_all),2):
    rows=rows_all.iloc[k:k+2];pr,yy,idx=forward(rows,False);assert pr.shape==(len(rows),len(idx)) and torch.isfinite(pr).all()
    if name!='scGPT_blood':pr=pr[:,common];yy=yy[:,common]
    predictions.append(pr.cpu().numpy());target.append(yy.cpu().numpy());details.extend(rows.to_dict('records'))
 torch.cuda.synchronize();prediction_seconds=time.time()-t
 prediction=np.concatenate(predictions);truth=np.concatenate(target);np.savez_compressed(out/'pilot_predictions.npz',predicted=prediction,observed=truth,gene_indices=common,gene_symbols=genes.gene_symbol.to_numpy()[common].astype(str));pd.DataFrame(details).to_csv(out/'prediction_rows.csv',index=False);pd.DataFrame(history).to_csv(out/'training_steps.csv',index=False)
 # Persist only newly trained parameters/buffers; frozen public weights stay in their source checkpoint.
 if name=='scFoundation_GEARS':state={k:v.cpu() for k,v in model.state_dict().items() if not k.startswith('singlecell_model.')}
 else:state={k:v.cpu() for k,v in model.state_dict().items()}
 torch.save(state,out/'pilot_final.pt')
 result={'model':name,'status':'passed','device':str(device),'gpu_name':torch.cuda.get_device_name(device),'optimizer_steps':4,'training_cells':len(train_rows),'prediction_cells':len(details),'shared_output_genes':len(common),'data_loading_and_model_init_seconds':loaded_seconds,'training_seconds':sum(x['seconds'] for x in history),'prediction_seconds':prediction_seconds,'total_seconds':time.time()-start,'peak_allocated_GiB':torch.cuda.max_memory_allocated(device)/2**30,'peak_reserved_GiB':torch.cuda.max_memory_reserved(device)/2**30,'parameters_total':sum(p.numel() for p in model.parameters()),'parameters_trainable':sum(p.numel() for p in model.parameters() if p.requires_grad),'parameter_update_verified':True,'frozen_gradient_check':frozen_ok,'predictions_finite':True,'checkpoint':'pilot_final.pt','checkpoint_role':'four-step pilot, not formal trained/best model','metadata':meta};js(out/'pilot_summary.json',result);print(json.dumps(result,default=str),flush=True)
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--model',required=True,choices=['GEARS','scGPT_blood','scFoundation_GEARS']);a.add_argument('--device',required=True);args=a.parse_args()
 try:main(args.model,args.device)
 except Exception as e:
  p=D/args.model;p.mkdir(parents=True,exist_ok=True);js(p/'pilot_error.json',{'model':args.model,'error':str(e),'traceback':traceback.format_exc()});raise
