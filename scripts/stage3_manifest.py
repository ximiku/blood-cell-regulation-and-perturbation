"""Create a portable snapshot of scientific files in the remote task directory."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'logs/stage3_remote_manifest.json'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

if __name__=='__main__':
    rows=[]
    for p in sorted(ROOT.rglob('*')):
        rel=p.relative_to(ROOT)
        if not p.is_file() or p==OUTPUT or any(s.startswith('.') or s=='__pycache__' for s in rel.parts):continue
        rows.append({'path':str(rel),'bytes':p.stat().st_size,'sha256':sha(p)})
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),'remote_root':str(ROOT),'excluded':'Derived caches and this manifest','files':rows}
    OUTPUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'files':len(rows),'bytes':sum(x['bytes'] for x in rows)}))
