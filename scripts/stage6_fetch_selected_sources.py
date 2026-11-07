from pathlib import Path
import urllib.request,json,hashlib,concurrent.futures,time
R=Path('/root/HW2_20261107/solution');V=R/'vendor/stage6';L=R/'logs/stage6';V.mkdir(parents=True,exist_ok=True);L.mkdir(parents=True,exist_ok=True)
def get(u):
 for k in range(3):
  try:return urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'HW2-reproducible-analysis'}),timeout=40).read()
  except Exception:
   if k==2:raise
   time.sleep(1)
def repo(name,owner,branch):
 meta=json.loads(get(f'https://api.github.com/repos/{owner}/commits/{branch}'));sha=meta['sha'];tree=json.loads(get(f'https://api.github.com/repos/{owner}/git/trees/{sha}?recursive=1'))
 prefixes={'GEARS':['gears/','demo/'],'scGPT':['scgpt/','tutorials/'],'scFoundation':['GEARS/gears/','GEARS/run_sh/','model/']}[name]
 keep=[x['path'] for x in tree['tree'] if x['type']=='blob' and (x['path'] in ['README.md','requirements.txt','setup.py','pyproject.toml','LICENSE','MODEL_LICENSE','OS_scRNA_gene_index.19264.tsv'] or (any(x['path'].startswith(p) for p in prefixes) and x['path'].endswith(('.py','.md','.sh','.json','.yaml','.yml'))))]
 if name=='scGPT':keep+=['tutorials/Tutorial_Perturbation.ipynb']
 if name=='scFoundation':keep+=['GEARS/README.md','GEARS/train.py']
 records=[]
 def fetch(path):
  url=f'https://raw.githubusercontent.com/{owner}/{sha}/{path}';b=get(url);p=V/name/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
  return {'path':path,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
 with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
  for rec in ex.map(fetch,dict.fromkeys(keep)):records.append(rec)
 return {'name':name,'repository':f'https://github.com/{owner}','commit':sha,'selection':'runtime source only; large repository data omitted','files':records}
records=[]
for item in [('GEARS','snap-stanford/GEARS','master'),('scGPT','bowang-lab/scGPT','main'),('scFoundation','biomap-research/scFoundation','main')]:
 try:
  r=repo(*item);records.append(r);print(item[0],len(r['files']),flush=True)
 except Exception as e:records.append({'name':item[0],'error':str(e)});print(item[0],str(e),flush=True)
 (L/'source_versions.json').write_text(json.dumps(records,indent=2)+'\n')
