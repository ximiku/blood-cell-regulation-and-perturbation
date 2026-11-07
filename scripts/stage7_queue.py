"""Four-GPU queue for nine independent, pre-specified runs; no extra agent or automation."""
from pathlib import Path
import subprocess,time,json,os
R=Path('/root/HW2_20261107/solution');L=R/'logs/stage7';L.mkdir(exist_ok=True);py='/root/HW2_20261107/envs/perturbation_sys/bin/python'
jobs=[('scFoundation_GEARS',1),('scGPT_blood',1),('GEARS',1),('GEARS',2),('scFoundation_GEARS',2),('scFoundation_GEARS',3),('scGPT_blood',2),('scGPT_blood',3),('GEARS',3)]
while not (R/'data/processed/stage7/prepared.json').exists():time.sleep(3)
base=subprocess.run([py,str(R/'scripts/stage7_baselines.py')],stdout=(L/'baselines.log').open('w'),stderr=subprocess.STDOUT)
if base.returncode:raise RuntimeError('baseline failure')
active={};done=[];failed=[];start=time.time()
while jobs or active:
 for gpu in range(4):
  if gpu not in active and jobs:
   name,seed=jobs.pop(0);log=(L/f'{name}_seed{seed}.log').open('a');env=os.environ.copy();env.update({'PYTHONHASHSEED':'0','OMP_NUM_THREADS':'6','MKL_NUM_THREADS':'6'})
   p=subprocess.Popen([py,'-u',str(R/'scripts/stage7_train.py'),'--model',name,'--seed',str(seed),'--device',f'cuda:{gpu}'],stdout=log,stderr=subprocess.STDOUT,env=env)
   active[gpu]={'p':p,'model':name,'seed':seed,'log':log,'started':time.time()}
 for gpu,item in list(active.items()):
  rc=item['p'].poll()
  if rc is not None:
   row={'model':item['model'],'seed':item['seed'],'gpu':gpu,'exit_code':rc,'seconds':time.time()-item['started']};(done if rc==0 else failed).append(row);item['log'].close();del active[gpu];print(row,flush=True)
 status={'queued':[{'model':m,'seed':s} for m,s in jobs],'active':[{'gpu':g,'model':x['model'],'seed':x['seed'],'pid':x['p'].pid,'elapsed_seconds':time.time()-x['started']} for g,x in active.items()],'completed':done,'failed':failed,'elapsed_seconds':time.time()-start,'time_unix':time.time()}
 (L/'queue_status.json').write_text(json.dumps(status,indent=2)+'\n')
 time.sleep(5)
(L/'queue_finished.json').write_text(json.dumps(status,indent=2)+'\n')
