"""Use GPUs released by the main queue for the planned core resampling supplement."""
from pathlib import Path
import json,time,subprocess,os
R=Path('/root/HW2_20261107/solution');L=R/'logs/stage7';py='/root/HW2_20261107/envs/gpu_lasso/bin/python'
pending=list(range(20));active={};done=[];failed=[]
while pending or active:
 try:s=json.loads((L/'queue_status.json').read_text())
 except (FileNotFoundError,json.JSONDecodeError):time.sleep(5);continue
 busy={x['gpu'] for x in s['active']}
 if not s['queued'] and not s['failed']:
  for gpu in range(4):
   if gpu in busy or gpu in active or not pending:continue
   rep=pending.pop(0);log=(L/f'core_retuning_{rep+1:02d}.log').open('a');env=os.environ.copy();env.update({'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2'})
   p=subprocess.Popen([py,'-u',str(R/'scripts/stage7_core_retuning.py'),'--device',str(gpu),'--shard',str(rep),'--shards','20'],stdout=log,stderr=subprocess.STDOUT,env=env);active[gpu]={'process':p,'repeat':rep+1,'log':log,'start':time.time()}
 for gpu,item in list(active.items()):
  code=item['process'].poll()
  if code is not None:
   row={'repeat':item['repeat'],'gpu':gpu,'exit_code':code,'seconds':time.time()-item['start']};(done if code==0 else failed).append(row);item['log'].close();del active[gpu]
 report={'pending':pending,'active':[{'gpu':g,'repeat':x['repeat'],'pid':x['process'].pid} for g,x in active.items()],'completed':done,'failed':failed,'time_unix':time.time()};(L/'retuning_queue_status.json').write_text(json.dumps(report,indent=2)+'\n');time.sleep(5)
(L/'retuning_finished.json').write_text(json.dumps(report,indent=2)+'\n')
