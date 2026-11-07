"""Verify remote/local byte identity and open final deliverables before instance shutdown."""
from pathlib import Path
import json,hashlib,re,time,zipfile
import numpy as np,pandas as pd,h5py
R=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**22),b''):h.update(b)
 return h.hexdigest()
start=time.time();checks=[]
for rel in ['logs/stage3_remote_manifest.json','logs/stage4_remote_manifest.json','logs/stage6/remote/remote_scientific_manifest.json','logs/stage8_remote_manifest.json']:
 manifest=json.loads((R/rel).read_text()); records=manifest['files']
 if rel in ['logs/stage3_remote_manifest.json','logs/stage4_remote_manifest.json']:records=[x for x in records if x['path'].startswith(('data/','results/'))]
 for x in records:
  p=R/x['path'];assert p.exists() and p.stat().st_size==x['bytes'],x['path'];assert sha(p)==x['sha256'],x['path']
 checks.append({'manifest':rel,'files':len(records),'bytes':sum(x['bytes'] for x in records),'passed':True})
for x in json.loads((R/'data/input_manifest.json').read_text()):
 p=Path(x['absolute_path']);assert p.stat().st_size==x['bytes'] and sha(p)==x['sha256'],x['role']
for report in ['08_enhanced_report.md','07_gpu_benchmark_report.md','08_archive_handoff.md','08_results_brief.md']:
 p=R/'reports'/report;text=p.read_text();assert len(text)>300
 for link in re.findall(r'\]\(([^)]+)\)',text):
  if link.startswith(('http://','https://','#')):continue
  assert (p.parent/link.split('#')[0]).exists(),link
for name in ['07_model_summary.csv','07_seed_metrics.csv','07_condition_metrics.csv','07_core_top10_retuning_check.csv']:
 assert len(pd.read_csv(R/'tables'/name))>0
with h5py.File(R/'data/processed/analysis_matrices.h5','r') as f:
 assert f['model/X'].shape==(11922,2717) and np.isfinite(f['model/X'][0,:20]).all() and np.isfinite(f['model/Y'][0,:20]).all()
with h5py.File(R/'data/external/stage6/norman_source.h5ad','r') as f:
 assert 'X' in f and 'obs' in f
 values=f['X/data'][:100] if isinstance(f['X'],h5py.Group) else f['X'][0,:100]
 assert np.isfinite(values).all()
z=np.load(R/'data/processed/stage7/scFoundation_reference_embeddings.npy',mmap_mode='r');assert z.shape==(64,19264,512) and np.isfinite(z[0]).all()
for m in ['GEARS','scGPT_blood','scFoundation_GEARS']:
 for seed in [1,2,3]:
  z=np.load(R/f'results/stage7/{m}/seed{seed}/test_predictions.npz');assert np.isfinite(z['predicted']).all()
  assert zipfile.is_zipfile(R/f'results/stage7/{m}/seed{seed}/best.pt')
output={'passed':True,'remote_local_manifests':checks,'reports_opened':4,'key_tables_opened':4,'core_and_external_data_opened':True,'all_model_predictions_opened':True,'elapsed_seconds':time.time()-start}
(R/'logs/stage8_archive_validation.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n');print(json.dumps(output,ensure_ascii=False))
