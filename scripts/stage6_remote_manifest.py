from pathlib import Path
import hashlib,json,subprocess,sys,importlib.metadata as im
R=Path('/root/HW2_20261107/solution');L=R/'logs/stage6';E=R/'environment';E.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**22),b''):h.update(b)
 return h.hexdigest()
records=[]
for rel in ['data/external/stage6','data/processed/stage6','data/splits/stage6','results/stage6','vendor/stage6']:
 for p in sorted((R/rel).rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:records.append({'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)})
for p in sorted((R/'scripts').glob('stage6*.py')):records.append({'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)})
(L/'remote_scientific_manifest.json').write_text(json.dumps({'files':records,'bytes':sum(x['bytes'] for x in records)},indent=2)+'\n')
weights=[]
for model in ['scGPT_blood','scFoundation']:
 for p in sorted((R/'models/stage6'/model).glob('*')):
  if p.is_file():weights.append({'model':model,'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)})
expect={'best_model.pt':'56e68896b1b4e8c708b695f1a60bb7560e08952ea76842f695109737beaeeb41','models.ckpt':'9f40bf324d3d0084c4b288d06f5af4fddd12206e2a3f022551d12e89e33a0ea9'}
for x in weights:
 if Path(x['path']).name in expect:assert x['sha256']==expect[Path(x['path']).name]
(L/'remote_public_weights_verified.json').write_text(json.dumps({'passed':True,'files':weights},indent=2)+'\n')
packages=['torch','numpy','scipy','pandas','scanpy','anndata','h5py','torch-geometric','einops','local-attention'];env={'python':sys.version,'executable':sys.executable,'packages':{p:im.version(p) for p in packages},'CUDA_VISIBLE_DEVICES':'all 4; trials on cuda:0,1,2','inherited_torch':'/opt/conda with system-site-packages; isolated additional dependencies'}
(E/'stage6_environment.json').write_text(json.dumps(env,indent=2)+'\n');(E/'requirements-stage6.txt').write_text(subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True));print(json.dumps({'scientific_files':len(records),'scientific_bytes':sum(x['bytes'] for x in records),'public_weights_verified':len(weights)}))
