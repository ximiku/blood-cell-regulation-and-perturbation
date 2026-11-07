"""Snapshot the local core-report delivery without overwriting earlier stages."""
from pathlib import Path
from datetime import datetime,timezone
import json
from stage4_manifest import sha
root=Path(__file__).resolve().parents[1];out=root/'logs/stage5_artifact_manifest.json';rows=[]
for p in sorted(root.rglob('*')):
    rel=p.relative_to(root)
    if not p.is_file() or p==out or p.suffix=='.tmp' or 'Autosaved' in p.name:continue
    if any(x.startswith('.') or x=='__pycache__' for x in rel.parts):continue
    rows.append({'path':str(rel),'bytes':p.stat().st_size,'sha256':sha(p)})
out.write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),'root':str(root),'scope':'all visible files excluding this manifest, caches, temporary and autosaved files','files':rows},ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'files':len(rows),'bytes':sum(r['bytes'] for r in rows),'manifest':str(out)}))
