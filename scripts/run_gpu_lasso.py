"""Distribute independent training contexts across four GPUs, then assemble CV."""
from pathlib import Path
from datetime import datetime, timezone
import json,os,subprocess,time

ROOT=Path(__file__).resolve().parents[1]
TASK=ROOT.parent
GPU_PY=TASK/'envs/gpu_lasso/bin/python'
CPU_PY=TASK/'envs/statistics/bin/python'

if __name__=='__main__':
    start=time.monotonic()
    contexts=sorted(p.name.removesuffix('_moments.npz') for p in (ROOT/'.cache/lasso').glob('*_moments.npz'))
    cfg=json.loads((ROOT/'config/analysis_config.json').read_text())
    expected=['full_fit']+[f'full_inner_{j:02d}' for j in range(1,cfg['inner_folds']+1)]
    for k in range(1,cfg['outer_folds']+1):
        expected += [f'outer_{k:02d}_fit']+[f'outer_{k:02d}_inner_{j:02d}' for j in range(1,cfg['inner_folds']+1)]
    if set(contexts)!=set(expected):
        raise RuntimeError(f'Unexpected training contexts: missing={set(expected)-set(contexts)}, extra={set(contexts)-set(expected)}')
    def completed(c):
        return all((ROOT/'results/lasso_gpu_paths'/f'{c}.{ext}').exists() for ext in ['json','npz'])
    pending=[c for c in contexts if not completed(c)]
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',CUPY_CACHE_DIR=str(ROOT/'.cache/cupy'))
    jobs=[]
    for device in range(4):
        group=pending[device::4]
        if not group:continue
        log=(ROOT/'logs'/f'03_gpu_device_{device}.log').open('w')
        command=[str(GPU_PY),'-u',str(ROOT/'scripts/lasso_gpu_paths.py'),'--device',str(device),'--contexts',*group]
        proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env)
        jobs.append((device,group,proc,log))
    errors=[]
    for device,group,proc,log in jobs:
        code=proc.wait();log.close()
        print(json.dumps({'event':'device_done','device':device,'exit_code':code,'contexts':group}),flush=True)
        if code:errors.append(device)
    if errors:raise RuntimeError(f'GPU jobs failed: {errors}')
    if not all(completed(c) for c in contexts):raise RuntimeError('A context is missing its path or completion record.')
    report={'completed_utc':datetime.now(timezone.utc).isoformat(),'parallel_path_wall_seconds':time.monotonic()-start,
            'contexts':contexts,'n_contexts':len(contexts),'devices':[0,1,2,3],
            'first_context_separate_log':'03_gpu_first_context.log'}
    (ROOT/'logs/stage3_gpu_execution.json').write_text(json.dumps(report,indent=2)+'\n')
    command=[str(CPU_PY),'-u',str(ROOT/'scripts/lasso_analysis.py'),'--correlations','--threads','2',
             '--workers','4','--targets',*[str(i) for i in range(21)]]
    with (ROOT/'logs/03_target_assembly.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
    print('All paths and target outputs completed.',flush=True)
