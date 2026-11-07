"""Check the stage-1 handoff and write a local artifact manifest."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b''): h.update(chunk)
    return h.hexdigest()


def main():
    summary = json.loads((ROOT/'logs/stage1_summary.json').read_text())
    checks = []
    for row in json.loads((ROOT/'data/input_manifest.json').read_text()):
        p = Path(row['absolute_path'])
        assert p.stat().st_size == row['bytes'] and digest(p) == row['sha256']
    checks.append('original_input_files_match_saved_hashes')
    for row in json.loads((ROOT/'data/references/source_manifest.json').read_text()):
        p=ROOT/row['file']
        assert p.stat().st_size == row['bytes'] and digest(p) == row['sha256']
    checks.append('official_reference_files_match_saved_hashes')
    data_path=ROOT/summary['model_matrix_file']
    assert digest(data_path)==summary['model_matrix_sha256']
    checks.append('processed_matrix_matches_saved_hash')
    ann=pd.read_csv(ROOT/'tables/01_perturbation_annotations.csv')
    ga=pd.read_csv(ROOT/'tables/01_shared_gene_annotations.csv')
    with h5py.File(data_path,'r') as f:
        genes=f['model/genes'].asstr()[:].tolist()
        perts=f['model/perturbations'].asstr()[:].tolist()
        targets=f['model/target_cells'].asstr()[:].tolist()
        blood=pd.DataFrame(f['blood/expression'][:],index=f['blood/genes'].asstr()[:],columns=f['blood/cell_types'].asstr()[:])
        X=f['model/X'][:]
        Y=f['model/Y'][:]
        assert list(X.shape)==summary['model_X_shape'] and list(Y.shape)==summary['model_Y_shape']
        assert len(genes)==len(set(genes)) and len(perts)==len(set(perts))
        assert np.isfinite(X).all() and np.isfinite(Y).all()
        assert (X.std(axis=0)>0).all()
        expected=blood.loc[genes,targets].subtract(blood.loc[genes,'hematopoietic stem cell'],axis=0).to_numpy()
        assert np.array_equal(Y,expected)
        checks.extend(['matrix_dimensions_unique_labels_finite_values','response_independently_reconstructed_from_saved_blood','no_constant_predictors'])
        assert perts==ann.perturbation_label.tolist()
        assert np.array_equal(f['model/is_TF'][:],ann.is_TF.to_numpy())
        ga=ga.loc[ga.in_model].set_index('gene_symbol').loc[genes]
        assert np.array_equal(f['model/is_landmark'][:],ga.is_landmark.eq(1).to_numpy())
        assert int(f['model/is_main_target'][:].sum())==20
        assert len(targets)==21 and 'hematopoietic stem cell' not in targets
        checks.extend(['TF_and_landmark_masks_match_feature_and_gene_order','twenty_main_and_one_supplementary_target'])
    required=['PLAN.md','WORK_LOG.md','README.md','reports/01_data_qc.md','environment/local_environment.json',
              'environment/remote_inventory.json','tables/01_TF_catalog_coverage.csv','tables/01_duplicate_gene_audit.csv']
    for name in required:
        p=ROOT/name
        assert p.exists() and p.stat().st_size>0
        if p.suffix=='.md': assert len(p.read_text())>100
    state=json.loads((ROOT/'stage_status.json').read_text())
    assert state['stage']==1 and state['status']=='awaiting_user_review' and state['next_stage_started'] is False
    checks.extend(['handoff_documents_present_and_readable','stage_two_has_not_started'])
    validation={'completed_utc':datetime.now(timezone.utc).isoformat(),'status':'passed','checks':checks,
                'training_started':False,'local_data_archive_verified':True}
    (ROOT/'logs/stage1_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    records=[]
    manifest=ROOT/'logs/stage1_artifact_manifest.json'
    for p in sorted(ROOT.rglob('*')):
        rel=p.relative_to(ROOT)
        if not p.is_file() or p==manifest or '.cache' in rel.parts or '__pycache__' in rel.parts:
            continue
        records.append({'path':str(rel),'bytes':p.stat().st_size,'sha256':digest(p)})
    manifest.write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),'files':records},ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'passed','checks':len(checks),'manifest_files':len(records),
                      'artifact_bytes':sum(r['bytes'] for r in records)},indent=2))


if __name__=='__main__': main()
