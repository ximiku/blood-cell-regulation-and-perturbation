"""Adopt running jobs and finish at the user-approved budget without repeating epochs."""
from pathlib import Path
import os,signal,subprocess,time,json,csv
R=Path('/root/HW2_20261107/solution');L=R/'logs/stage7';py='/root/HW2_20261107/envs/perturbation_sys/bin/python'
original=json.loads((L/'queue_before_budget_change.json').read_text());cap=json.loads((R/'config/stage7_config.json').read_text())['max_epochs']
jobs=original['queued'];done=original['completed'];failed=original['failed'];start=original['time_unix']-original['elapsed_seconds'];active={}
def save(path,value):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)
def running(pid):
 try:return Path(f'/proc/{pid}/stat').read_text().split(') ')[1][0]!='Z'
 except FileNotFoundError:return False
def launch(gpu,name,seed,prior=0):
 log=(L/f'{name}_seed{seed}.log').open('a');env=os.environ.copy();env.update(PYTHONHASHSEED='0',OMP_NUM_THREADS='6',MKL_NUM_THREADS='6')
 p=subprocess.Popen([py,'-u',str(R/'scripts/stage7_train.py'),'--model',name,'--seed',str(seed),'--device',f'cuda:{gpu}','--prior-seconds',str(prior)],stdout=log,stderr=subprocess.STDOUT,env=env)
 return dict(model=name,seed=seed,pid=p.pid,p=p,log=log,started=time.time()-prior,adopted=False)
for x in original['active']:
 active[x['gpu']]=dict(model=x['model'],seed=x['seed'],pid=x['pid'],started=original['time_unix']-x['elapsed_seconds'],adopted=True)
while jobs or active:
 for gpu,item in list(active.items()):
  out=R/f"results/stage7/{item['model']}/seed{item['seed']}";history=out/'training_history.csv';latest=out/'latest.pt'
  if item['adopted'] and not (out/'completed.json').exists() and history.exists() and latest.exists():
   h=list(csv.DictReader(history.open()))
   if len(h)>=cap and latest.stat().st_mtime>=history.stat().st_mtime:
    # Checkpoint is atomically saved after CSV. A completed checkpoint is safe to resume.
    os.kill(item['pid'],signal.SIGTERM);os.kill(item['pid'],signal.SIGCONT)
    while running(item['pid']):time.sleep(.2)
    prior=time.time()-item['started'];active[gpu]=launch(gpu,item['model'],item['seed'],prior)
    print('budget_finalize',item['model'],item['seed'],len(h),flush=True);continue
  alive=running(item['pid']) if item['adopted'] else item['p'].poll() is None
  if not alive:
   ok=(out/'completed.json').exists();row=dict(model=item['model'],seed=item['seed'],gpu=gpu,exit_code=0 if ok else 1,seconds=time.time()-item['started']);(done if ok else failed).append(row)
   if 'log' in item:item['log'].close()
   del active[gpu];print(row,flush=True)
 for gpu in range(4):
  if gpu not in active and jobs:
   x=jobs.pop(0);active[gpu]=launch(gpu,x['model'],x['seed'])
 status=dict(queued=jobs,active=[dict(gpu=g,model=x['model'],seed=x['seed'],pid=x['pid'],elapsed_seconds=time.time()-x['started']) for g,x in active.items()],completed=done,failed=failed,elapsed_seconds=time.time()-start,time_unix=time.time(),max_epochs=cap)
 save(L/'queue_status.json',status);time.sleep(2)
save(L/'queue_finished.json',status)
