from pathlib import Path
import json,csv,time,subprocess
R=Path('/root/HW2_20261107/solution');L=R/'logs/stage7'
def load(p):
 try:return json.loads(p.read_text())
 except (FileNotFoundError,json.JSONDecodeError):return {}
q=load(L/'queue_status.json');runs=[]
for p in sorted((R/'results/stage7').glob('*/seed*/run_config.json')):
 o=p.parent;h=[]
 if (o/'training_history.csv').exists():
  with (o/'training_history.csv').open() as f:h=list(csv.DictReader(f))
 progress=load(o/'progress.json');completed=load(o/'completed.json');error=load(o/'error.json')
 runs.append({'model':o.parent.name,'seed':o.name,'epochs_done':len(h),'last':h[-1] if h else None,'progress':progress,'completed':bool(completed),'error':error or None})
ret=load(L/'retuning_queue_status.json')
print(json.dumps({'time_unix':time.time(),'queued':q.get('queued'),'active':q.get('active'),'finished':q.get('completed'),'failed':q.get('failed'),'runs':runs,'retuning':{'pending':len(ret.get('pending',[])),'active':ret.get('active'),'completed':len(ret.get('completed',[])),'failed':ret.get('failed')},'gpu':subprocess.check_output(['nvidia-smi','--query-gpu=index,memory.used,utilization.gpu','--format=csv,noheader'],text=True)},ensure_ascii=False))
