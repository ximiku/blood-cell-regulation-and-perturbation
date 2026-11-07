from pathlib import Path
import urllib.request,json,tarfile,hashlib,concurrent.futures,io,time
R=Path('/root/HW2_20261107/solution');V=R/'vendor/stage6';L=R/'logs/stage6';V.mkdir(parents=True,exist_ok=True);L.mkdir(parents=True,exist_ok=True)
def get(u):
 req=urllib.request.Request(u,headers={'User-Agent':'HW2-reproducible-analysis'})
 return urllib.request.urlopen(req,timeout=60).read()
def repo(item):
 name,owner,branch=item;meta=json.loads(get(f'https://api.github.com/repos/{owner}/commits/{branch}'));sha=meta['sha'];u=f'https://codeload.github.com/{owner}/tar.gz/{sha}';b=get(u);dest=V/name
 dest.mkdir(exist_ok=True)
 with tarfile.open(fileobj=io.BytesIO(b),mode='r:gz') as tf:
  for member in tf.getmembers():
   parts=Path(member.name).parts[1:]
   if not parts:continue
   rel=Path(*parts)
   if rel.is_absolute() or '..' in rel.parts or member.issym() or member.islnk():raise RuntimeError('Unsupported archive member')
   out=dest/rel
   if member.isdir():out.mkdir(parents=True,exist_ok=True)
   elif member.isfile():out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(tf.extractfile(member).read())
 return {'name':name,'repository':f'https://github.com/{owner}','commit':sha,'archive_url':u,'archive_sha256':hashlib.sha256(b).hexdigest(),'archive_bytes':len(b)}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 items=list(ex.map(repo,[('GEARS','snap-stanford/GEARS','master'),('scGPT','bowang-lab/scGPT','main'),('scFoundation','biomap-research/scFoundation','main')]))
(L/'source_versions.json').write_text(json.dumps(items,indent=2)+'\n');print(json.dumps(items),flush=True)
f=json.loads(get('https://api.figshare.com/v2/articles/24049200'));(L/'figshare_article.json').write_text(json.dumps(f,indent=2)+'\n');print(json.dumps({'files':[{k:x.get(k) for k in ['id','name','size','download_url','computed_md5']} for x in f['files']]}),flush=True)
