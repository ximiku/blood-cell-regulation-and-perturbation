"""Copy each finished experiment during remaining training; stop when both queues finish."""
from pathlib import Path
import subprocess,json,time
R=Path(__file__).resolve().parents[1];remote='/root/HW2_20261107/solution'
host='root@gpu.example.org';ssh=['ssh','-i','/path/to/ssh_private_key','-p','22','-o','BatchMode=yes','-o','ControlPath=/tmp/hw2-ssh-%C']
copied=set();last=None
def sync(rel,exclude=False):
 dest=R/rel;dest.mkdir(parents=True,exist_ok=True)
 cmd=['rsync','-a','--partial']+(['--exclude=latest.pt'] if exclude else [])+['-e',' '.join(ssh),f'{host}:{remote}/{rel}/',str(dest)+'/']
 subprocess.run(cmd,check=True)
while True:
 # Read independent JSON files in one SSH call, using explicit lengths via remote Python.
 cmd=f"/root/HW2_20261107/envs/perturbation_sys/bin/python -c \"import json,pathlib; p=pathlib.Path('{remote}/logs/stage7'); print(json.dumps([json.loads((p/'queue_status.json').read_text()),json.loads((p/'retuning_queue_status.json').read_text())]))\""
 raw=subprocess.run(ssh+[host,cmd],capture_output=True,text=True)
 if raw.returncode:print('temporary SSH retry',flush=True);time.sleep(30);continue
 try:q,ret=json.loads(raw.stdout)
 except json.JSONDecodeError:time.sleep(5);continue
 if q['failed'] or ret['failed']:raise RuntimeError({'main':q['failed'],'retuning':ret['failed']})
 for x in q['completed']:
  key=(x['model'],x['seed'])
  if key in copied:continue
  sync(f'results/stage7/{key[0]}/seed{key[1]}',True);copied.add(key);print('archived',key,flush=True)
 state=(len(q['completed']),len(ret['completed']))
 if state!=last:print('completed',state,'active',[(x['model'],x['seed']) for x in q['active']],flush=True);last=state
 if not q['active'] and not q['queued'] and len(q['completed'])==9 and not ret['active'] and not ret['pending'] and len(ret['completed'])==20:
  sync('results/stage7/core_retuning');sync('logs/stage7');print('all formal results copied',flush=True);break
 time.sleep(45)
