"""Read-only remote inventory. No remote installation or training is performed."""
from datetime import datetime, timezone
import json
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REMOTE_CODE = r'''
import importlib.metadata as md, json, os, pathlib, platform, shutil, subprocess
result = {'host_platform': platform.platform(), 'python': platform.python_version(),
          'cpu_count': os.cpu_count(), 'cwd': os.getcwd(), 'paths': {}}
for p in ['/root','/tmp','/app','/data','/workspace','/mnt','/data1','/data2','/public','/datasets']:
    if pathlib.Path(p).exists():
        d=shutil.disk_usage(p)
        result['paths'][p]={'total_bytes':d.total,'free_bytes':d.free,
                            'top_level_entries': sorted(os.listdir(p))[:30] if pathlib.Path(p).is_dir() else []}
result['memory_kb']={l.split(':')[0]:l.split(':')[1].strip() for l in pathlib.Path('/proc/meminfo').read_text().splitlines() if l.startswith(('MemTotal:','MemAvailable:'))}
result['cgroup_limits']={}
for key in ['memory.max','cpu.max','cpuset.cpus.effective']:
    p=pathlib.Path('/sys/fs/cgroup')/key
    if p.exists(): result['cgroup_limits'][key]=p.read_text().strip()
result['packages']={}
for name in ['torch','numpy','pandas','scipy','scikit-learn','h5py','anndata','scanpy','scgpt','cell-gears','pip']:
    try: result['packages'][name]=md.version(name)
    except md.PackageNotFoundError: result['packages'][name]=None
for name, command in {
    'gpus':['nvidia-smi','--query-gpu=name,memory.total,memory.free,driver_version','--format=csv,noheader'],
    'filesystems':['df','-h'],
    'tools':['sh','-c','command -v python3; command -v pip3; command -v conda; command -v rsync; command -v git']
}.items():
    p=subprocess.run(command,capture_output=True,text=True,timeout=20)
    result[name]={'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
print(json.dumps(result))
'''
command = ['ssh','-i',str(Path.home()/'.ssh/id_ed25519'),'-p','22',
           '-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=15',
           'root@gpu.example.org','python3 -c '+shlex.quote(REMOTE_CODE)]
p = subprocess.run(command,capture_output=True,text=True,timeout=60)
output = ROOT/'environment'
output.mkdir(parents=True,exist_ok=True)
if p.returncode:
    (output/'remote_inspection_error.txt').write_text(p.stderr)
    raise SystemExit(p.returncode)
data = json.loads(p.stdout)
data['inspected_utc'] = datetime.now(timezone.utc).isoformat()
data['connection'] = {'host':'gpu.example.org','port':22,'user':'root'}
(output/'remote_inventory.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(data,ensure_ascii=False,indent=2))
