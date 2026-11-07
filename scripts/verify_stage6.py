"""Local scientific and archival checks for the completed minimal GPU trials."""
from pathlib import Path
import hashlib, json, time, zipfile
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py

R = Path(__file__).resolve().parents[1]
P = R / 'data/processed/stage6'
L = R / 'logs/stage6'
checks = []
def check(name, condition, **details):
    checks.append({'check': name, 'passed': bool(condition), **details})
    assert condition, name
def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(2**22), b''):
            h.update(block)
    return h.hexdigest()

start = time.time()
manifest = json.loads((L / 'remote/remote_scientific_manifest.json').read_text())
for entry in manifest['files']:
    p = R / entry['path']
    assert p.exists() and p.stat().st_size == entry['bytes'], entry['path']
    assert sha(p) == entry['sha256'], entry['path']
check('remote_to_local_path_bytes_sha256', True,
      files=len(manifest['files']), bytes=manifest['bytes'])

X = sp.load_npz(P / 'expression_log1p.npz')
cells = pd.read_csv(P / 'cells.csv.gz')
genes = pd.read_csv(P / 'genes.csv')
indices = np.load(P / 'evaluation_indices.npz')
common = indices['common']
check('matrix_dimensions_order_finite', X.shape == (84143, 19264)
      and len(cells) == X.shape[0] and len(genes) == X.shape[1]
      and cells.cell_id.is_unique and genes.gene_symbol.is_unique
      and np.array_equal(genes.gene_index, np.arange(X.shape[1]))
      and np.isfinite(X.data).all() and (X.data >= 0).all())
check('shared_gene_indices', len(common) == 18977
      and np.array_equal(common, np.flatnonzero(genes.shared_output)))
check('total_counts_preserved', np.array_equal(np.load(P/'total_counts.npy'), cells.total_count))
check('target_gene_coverage', set(json.loads((P/'perturbation_genes.json').read_text()))
      <= set(genes.loc[genes.shared_output, 'gene_symbol']))
split_table = pd.read_csv(R/'data/splits/stage6/condition_splits.csv')
for seed in [1, 2, 3]:
    labels = cells[f'split_seed{seed}']
    split = json.loads((R/f'data/splits/stage6/seed{seed}.json').read_text())['conditions']
    sets = {k: set(v) for k,v in split.items()}
    check(f'seed{seed}_condition_disjoint',
          not (sets['train'] & sets['val'] or sets['train'] & sets['test'] or sets['val'] & sets['test'])
          and set.union(*sets.values()) == set(cells.condition)
          and all(set(cells.loc[labels == k, 'condition']) == v for k,v in sets.items())
          and set(labels[cells.condition == 'ctrl']) == {'train'})
    train = X[labels.to_numpy() == 'train'][:, common]
    variance = np.asarray(train.power(2).mean(0)).ravel() - np.asarray(train.mean(0)).ravel()**2
    expected_hvg = common[np.argsort(-variance, kind='stable')[:1000]]
    check(f'seed{seed}_training_only_HVG_recomputed', np.array_equal(expected_hvg, indices[f'hvg1000_seed{seed}']))
    trained = {g for c in sets['train'] for g in c.split('+') if g != 'ctrl'}
    rows = split_table[split_table.seed == seed]
    check(f'seed{seed}_actual_seen_recomputed', all(
        int(row.seen_in_final_training) == sum(g in trained for g in row.condition.split('+') if g != 'ctrl')
        for row in rows.itertuples()))
    del train

ctrl = np.flatnonzero(cells.condition.to_numpy() == 'ctrl')
check('control_mean_recomputed', np.allclose(np.asarray(X[ctrl].mean(0)).ravel(), np.load(P/'control_mean.npy'), atol=1e-7, rtol=0))
pairs = pd.read_csv(R/'data/splits/stage6/pilot_pairs.csv')
check('pilot_pair_identity_and_training_controls', all(
    cells.iloc[r.target_row].cell_id == r.target_cell
    and cells.iloc[r.target_row].condition == r.condition
    and cells.iloc[r.target_row].split_seed1 == r.split
    and cells.iloc[r.control_row].condition == 'ctrl'
    and cells.iloc[r.control_row].split_seed1 == 'train' for r in pairs.itertuples()))

# Reopen the publisher H5AD and compare real pilot input/target rows after filtering.
with h5py.File(R/'data/external/stage6/norman_source.h5ad', 'r') as f:
    index_key = f['obs'].attrs['_index']
    if isinstance(index_key, bytes): index_key = index_key.decode()
    source_ids = f['obs'][index_key].asstr()[:]
    source_row = {c:i for i,c in enumerate(source_ids)}
    indptr = f['X/indptr'][:]
    for row in sorted(set(pairs.target_row) | set(pairs.control_row)):
        j = source_row[cells.iloc[row].cell_id]
        a,b = indptr[j:j+2]
        original = np.zeros(X.shape[1], dtype=np.float32)
        original[f['X/indices'][a:b]] = f['X/data'][a:b]
        assert np.array_equal(original, X[row].toarray()[0]), row
check('publisher_matrix_pilot_rows_exact', True)

expected_rows = pd.concat([pairs[pairs.split == 'val'], pairs[pairs.split == 'test']], ignore_index=True)
for model in ['GEARS', 'scGPT_blood', 'scFoundation_GEARS']:
    out = R/'results/stage6'/model
    result = json.loads((out/'pilot_summary.json').read_text())
    pred = np.load(out/'pilot_predictions.npz')
    rows = pd.read_csv(out/'prediction_rows.csv')
    steps = pd.read_csv(out/'training_steps.csv')
    check(model+'_shared_rows_genes_truth', rows.equals(expected_rows)
          and np.array_equal(pred['gene_indices'], common)
          and np.array_equal(pred['gene_symbols'], genes.gene_symbol.to_numpy()[common])
          and pred['predicted'].shape == (14,18977)
          and np.array_equal(pred['observed'], X[rows.target_row.to_numpy()][:,common].toarray()))
    check(model+'_finite_predictions_and_real_updates',
          np.isfinite(pred['predicted']).all() and len(steps) == 4
          and np.isfinite(steps.loss).all()
          and np.isfinite(steps.gradient_norm_before_clip).all()
          and (steps.gradient_norm_before_clip > 0).all()
          and result['parameter_update_verified'] and result['frozen_gradient_check']
          and result['status'] == 'passed')
    check(model+'_checkpoint_readable_zip', zipfile.is_zipfile(out/'pilot_final.pt'))
    with zipfile.ZipFile(out/'pilot_final.pt') as z:
        assert any(n.endswith('data.pkl') for n in z.namelist())

audit = {'passed': True, 'checks': checks, 'check_count': len(checks),
         'review_mode':'self_checks_only_user_requested', 'completed_utc': datetime.now(timezone.utc).isoformat(),
         'elapsed_seconds': time.time()-start}
(L/'local_validation.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({'passed':True, 'checks':len(checks), 'files':len(manifest['files']),
                  'bytes':manifest['bytes'], 'elapsed_seconds':audit['elapsed_seconds']}))
