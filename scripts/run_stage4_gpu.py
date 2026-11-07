"""Run the four-GPU stage-4 job list with paired-file recovery."""
from pathlib import Path
import os,json,subprocess,time,argparse
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--jobs',default='data/splits/stage4_jobs.json');args=p.parse_args()
    start=time.monotonic();jobs=json.loads((ROOT/args.jobs).read_text())['jobs']
    handles=[];env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',CUPY_CACHE_DIR=str(ROOT/'.cache/cupy'))
    tag=next((t for t in ['self_readout','stability_cold','size_control'] if t in args.jobs),'models_stability')
    (ROOT/'logs').mkdir(exist_ok=True)
    for i in range(4):
        log=(ROOT/'logs'/f'04_{tag}_gpu_{i}.log').open('w')
        cmd=[str(ROOT.parent/'envs/gpu_lasso/bin/python'),'-u',str(ROOT/'scripts/stage4_gpu.py'),
            '--device',str(i),'--shard',str(i),'--n-shards','4','--jobs',args.jobs]
        handles.append((i,subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env),log))
    failed=[]
    for i,proc,log in handles:
        status=proc.wait();log.close()
        if status:failed.append(i)
        print(json.dumps({'device':i,'exit_code':status}),flush=True)
    missing=[j['id'] for j in jobs if not all((ROOT/'results/stage4'/f'{j["id"]}.{ext}').is_file() for ext in ['npz','json'])]
    result={'seconds':time.monotonic()-start,'n_jobs':len(jobs),'failed_devices':failed,'missing':missing}
    (ROOT/'logs'/f'stage4_{tag}_execution.json').write_text(json.dumps(result,indent=2)+'\n')
    if failed or missing:raise RuntimeError(result)
