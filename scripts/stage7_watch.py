"""Compact foreground status stream for the current authorized batch."""
from pathlib import Path
import json,csv,time
R=Path('/root/HW2_20261107/solution');L=R/'logs/stage7'
def load(p):
 try:return json.loads(p.read_text())
 except (FileNotFoundError,json.JSONDecodeError):return {}
while True:
 q=load(L/'queue_status.json');r=load(L/'retuning_queue_status.json');runs=[]
 for x in q.get('active',[]):
  out=R/f"results/stage7/{x['model']}/seed{x['seed']}";h=[]
  if (out/'training_history.csv').exists():
   with (out/'training_history.csv').open() as f:h=list(csv.DictReader(f))
  p=load(out/'progress.json');last=h[-1] if h else {};runs.append({'model':x['model'],'seed':x['seed'],'epochs':len(h),'best_epoch':last.get('best_epoch'),'best_val_MSE':last.get('best_validation_MSE'),'last_epoch_minutes':round((float(last.get('training_seconds',0))+float(last.get('validation_seconds',0)))/60,2),'current_rows':p.get('rows_processed'),'training_rows':p.get('training_rows')})
 row={'utc':time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime()),'GPU_completed':len(q.get('completed',[])),'GPU_queued':len(q.get('queued',[])),'active':runs,'GPU_failed':q.get('failed',[]),'retuning_completed':len(r.get('completed',[])),'retuning_failed':r.get('failed',[])}
 print(json.dumps(row),flush=True)
 if q.get('failed') or r.get('failed'):break
 if len(q.get('completed',[]))==9 and len(r.get('completed',[]))==20:break
 time.sleep(55)
