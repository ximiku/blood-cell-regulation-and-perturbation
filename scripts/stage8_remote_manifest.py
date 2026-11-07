"""Manifest new remote scientific artifacts and model files left on the instance."""
from pathlib import Path
import json,hashlib,time
R=Path('/root/HW2_20261107/solution');L=R/'logs';records=[];remote_models=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**22),b''):h.update(b)
 return h.hexdigest()
def entry(p):return {'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)}
q=json.loads((L/'stage7/queue_finished.json').read_text());ret=json.loads((L/'stage7/retuning_finished.json').read_text())
assert not q['active'] and not q['queued'] and not q['failed'] and len(q['completed'])==9
assert not ret['active'] and not ret['pending'] and not ret['failed'] and len(ret['completed'])==20
for rel in ['data/processed/stage7','results/stage7','results/stage7_exploratory','logs/stage7']:
 for p in sorted((R/rel).rglob('*')):
  if not p.is_file() or '__pycache__' in p.parts:continue
  if p.suffix=='.pt' and (p.name!='best.pt' or 'stage7_exploratory' in p.parts):remote_models.append(entry(p));continue
  records.append(entry(p))
for p in sorted((R/'scripts').glob('stage7*.py')):records.append(entry(p))
for p in [R/'scripts/stage8_remote_manifest.py',R/'config/stage7_config.json']:records.append(entry(p))
for p in sorted((R/'models/stage6').rglob('*')):
 if p.is_file():remote_models.append(entry(p))
result={'created_unix':time.time(),'all_nine_GPU_runs_completed':True,'all_20_core_retuning_resamples_completed':True,'files':records,'file_count':len(records),'bytes':sum(x['bytes'] for x in records),'scope':'stage7 new data, formal best weights, predictions, core retuning, exploratory non-model records, logs, executed sources and configuration'}
(L/'stage8_remote_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
(L/'stage8_remote_model_inventory.json').write_text(json.dumps({'files':remote_models,'file_count':len(remote_models),'bytes':sum(x['bytes'] for x in remote_models),'note':'Public pretraining and non-best/resume/aborted checkpoints remain remote; formal best checkpoints are archived locally.'},indent=2)+'\n')
print(json.dumps({'files':len(records),'bytes':result['bytes'],'remote_only_model_files':len(remote_models)}))
