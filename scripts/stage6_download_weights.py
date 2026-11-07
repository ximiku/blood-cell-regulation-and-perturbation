from pathlib import Path
import urllib.request,json,hashlib,concurrent.futures,time,os
R=Path(os.environ.get('STAGE6_ROOT','/root/HW2_20261107/solution'));P=Path(os.environ.get('STAGE6_WEIGHTS',str(R/'models/stage6/scGPT_blood')));P.mkdir(parents=True,exist_ok=True)
files={'args.json':'1y4UJVflGl-b2qm-fvpxIoQ3XcC2umjj0','vocab.json':'127FdcUyY1EM7rQfAS0YI4ms6LwjmnT9J','best_model.pt':'1MJaavaG0ZZkC_yPO4giGRnuCe3F1zt30'}
def fetch(item):
 name,ident=item;url=f'https://drive.usercontent.google.com/download?id={ident}&export=download&confirm=t';p=P/name;t=time.time();h=hashlib.sha256()
 with urllib.request.urlopen(url,timeout=60) as r,p.open('wb') as f:
  first=r.read(8192)
  if b'<html' in first.lower() or b'<!doctype html' in first.lower():raise RuntimeError(f'{name}: HTML instead of file')
  f.write(first);h.update(first)
  while True:
   b=r.read(2**20)
   if not b:break
   f.write(b);h.update(b)
 out={'file':name,'source_url':url,'official_folder':'https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU','bytes':p.stat().st_size,'sha256':h.hexdigest(),'download_seconds':time.time()-t};print(json.dumps(out),flush=True);return out
out=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 for fut in concurrent.futures.as_completed([ex.submit(fetch,x) for x in files.items()]):
  try:out.append(fut.result())
  except Exception as e:out.append({'error':str(e)});print(str(e),flush=True)
  (R/'logs/stage6/scgpt_weight_manifest.json').write_text(json.dumps(out,indent=2)+'\n')
