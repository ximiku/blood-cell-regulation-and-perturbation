"""Verify report links, facts and preservation of approved scientific results."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
body=(R/'reports/05_core_report.md').read_text(); links=[]
for name in ['05_core_report.md','05_core_supplement.md','05_literature_notes.md','05_independent_review.md']:
    p=R/'reports'/name
    for ref in re.findall(r'\]\(([^)]+)\)',p.read_text()):
        if '://' not in ref:assert (p.parent/ref).exists(),(name,ref);links.append(ref)
old={x['path']:x for x in json.loads((R/'logs/stage4_artifact_manifest.json').read_text())['files']}
images=re.findall(r'!\[[^\]]*\]\(([^)]+)\)',body);assert len(images)==4
for ref in images:
    p=(R/'reports'/ref).resolve();assert sha(p)==old[str(p.relative_to(R.resolve()))]['sha256']
refs=pd.read_csv(R/'tables/05_references.csv').fillna('');assert len(refs)==16 and refs.id.nunique()==16
assert set(re.findall(r'\[(R\d\d)：',body))==set(refs.id)
source=pd.read_csv(R/'tables/04_stability_all_TFs.csv');ev=pd.read_csv(R/'tables/05_TF_literature_evidence.csv');assert len(ev)==7
for _,r in ev.iterrows():
    s=source[(source.target==r.target)&(source.perturbation==r.TF)].iloc[0]
    assert r.TF_stability_rank==s.TF_stability_rank and r.stable_priority==s.stable_priority
    np.testing.assert_allclose([r.selection_frequency,r.coefficient_standardized],[s.selection_frequency,s.full_coefficient],atol=1e-14)
b=ev[ev.TF=='BCL6'].iloc[0];assert b.TF_stability_rank==15 and b.self_readout_check_status=='not_in_prespecified_top10_check'
m=pd.read_csv(R/'tables/03_model_metrics.csv').query('model=="Lasso_min"').set_index('target')
for t in ['B cell','plasma cell','monocyte','erythrocyte']:assert f'{m.loc[t,"R2"]:.3f}' in body
assert '41、43、52、42' in body and '7 个稳定规范身份' in body
n=0
for rel,r in old.items():
    if rel.startswith(('results/','tables/01_','tables/02_','tables/03_','tables/04_')) or rel=='data/processed/analysis_matrices.h5':
        p=R/rel;assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],rel;n+=1
result={'generated_utc':datetime.now(timezone.utc).isoformat(),'passed':True,'new_model_fits':0,'local_links':len(links),'unchanged_main_figures':4,'references_cited':16,'candidate_cases_checked':7,'BCL6_scope_correct':True,'headline_metrics_correct':True,'earlier_scientific_files_unchanged':n}
(R/'logs/stage5_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(result,ensure_ascii=False))
