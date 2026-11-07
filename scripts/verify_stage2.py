"""Read-only checks of stage-2 artifacts; writes only its validation record."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import xml.etree.ElementTree as ET
import h5py
import numpy as np
import pandas as pd
from PIL import Image
from PyPDF2 import PdfReader
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / 'tables'
checks = []


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(2**20), b''):
            h.update(chunk)
    return h.hexdigest()


def passed(name, detail):
    checks.append({'check': name, 'passed': True, 'detail': detail})


def same(a, b, atol=1e-9, rtol=1e-9):
    np.testing.assert_allclose(a, b, atol=atol, rtol=rtol)


def read(name, index=None):
    return pd.read_csv(TAB / name, index_col=index)


def zscore(a):
    return (a - a.mean(axis=1, keepdims=True)) / a.std(axis=1, ddof=1, keepdims=True)


def main():
    cfg = json.loads((ROOT / 'config/analysis_config.json').read_text())
    summary = json.loads((ROOT / 'logs/stage2_eda_summary.json').read_text())
    old = json.loads((ROOT / 'logs/stage1_artifact_manifest.json').read_text())
    h5 = ROOT / 'data/processed/analysis_matrices.h5'
    expected = next(f for f in old['files'] if f['path'] == str(h5.relative_to(ROOT)))
    assert sha(h5) == expected['sha256']
    for item in json.loads((ROOT / 'data/input_manifest.json').read_text()):
        assert sha(Path(item['absolute_path'])) == item['sha256']
    passed('frozen_inputs_unchanged', 'HDF5 and four original input files match the stage-1 SHA-256 records.')

    with h5py.File(h5, 'r') as f:
        a = f['blood/expression'][:]
        genes = f['blood/genes'].asstr()[:].tolist()
        cells = f['blood/cell_types'].asstr()[:].tolist()
        model_genes = f['model/genes'].asstr()[:].tolist()
    main_cells = [c for c in cells if c != cfg['supplementary_target']]
    a = a[:, [cells.index(c) for c in main_cells]]
    assert a.shape == (19994, 21)
    variance = np.var(a, axis=1, ddof=1)
    ranking = np.argsort(-variance, kind='stable')
    for n in [100, 500, 1000]:
        table = read(f'02_HVG_{n}.csv')
        assert table.gene_symbol.tolist() == [genes[i] for i in ranking[:n]]
        same(table.variance_21_cells.values, variance[ranking[:n]])
    z = zscore(a[ranking[:1000]])
    saved_z = read('02_HVG1000_zscores.csv', 0)
    assert saved_z.columns.tolist() == main_cells
    same(saved_z.values, z)
    same(read('02_heatmap_HVG100_gene_zscores.csv', 0).values, zscore(a[ranking[:100]]))
    passed('HVG_selection_and_scaling', 'All 100/500/1000 gene lists, variances and displayed z-scores match direct recomputation.')

    matrix = z.T - z.T.mean(axis=0)
    u, s, vt = np.linalg.svd(matrix, full_matrices=False)
    ratios = s**2 / np.sum(s**2)
    same(read('02_PCA_explained_variance.csv').explained_variance_ratio, ratios)
    coords = read('02_PCA_coordinates.csv', 0)
    assert coords.index.tolist() == main_cells
    # Singular vectors may have either sign. Compare their absolute inner products
    # and the reconstructed sample covariance, preserving all 10 saved components.
    expected_scores = u[:, :10] * s[:10]
    same(coords.values @ coords.values.T, expected_scores @ expected_scores.T, atol=1e-7)
    same(summary['pca_first2_cumulative'], ratios[:2].sum())
    same(summary['pca_first3_cumulative'], ratios[:3].sum())
    passed('PCA', 'Independent NumPy SVD agrees with saved explained variance and 10-dimensional sample geometry.')

    hsc = main_cells.index(cfg['reference_cell'])
    delta = a - a[:, [hsc]]
    rms = np.sqrt(np.mean(delta**2, axis=0))
    summary_table = read('02_HSC_difference_summary.csv', 0).loc[main_cells]
    same(summary_table.rms_all_19994_genes.values, rms)
    lookup = {g: i for i, g in enumerate(genes)}
    model_i = [lookup[g] for g in model_genes]
    same(summary_table.rms_model_11922_genes, np.sqrt(np.mean(delta[model_i]**2, axis=0)))
    focus_i = [main_cells.index(c) for c in cfg['focus_targets']]
    saved_delta = read('02_focus_HSC_differences.csv', 0)
    assert saved_delta.index.tolist() == genes
    same(saved_delta.values, delta[:, focus_i])
    markers = read('02_marker_panel_manifest.csv').gene_symbol.tolist()
    marker_z = zscore(a[[lookup[g] for g in markers]])
    same(read('02_marker_gene_zscore_HSC_difference.csv', 0).values, marker_z[:, focus_i] - marker_z[:, [hsc]])
    same(read('02_marker_input_expression.csv', 0).values, a[[lookup[g] for g in markers]])
    passed('HSC_differences', 'All-gene/model-gene RMS, full focus differences and all 29 marker values agree with frozen inputs.')

    dist = squareform(pdist(z.T, metric='euclidean'))
    corr = np.corrcoef(z.T)
    same(read('02_HVG1000_euclidean_distance.csv', 0).values, dist)
    same(read('02_HVG1000_gene_zscore_Pearson.csv', 0).values, corr)
    other = [i for i in range(len(main_cells)) if i != hsc]
    nearest = main_cells[other[np.argmin(dist[hsc, other])]]
    assert nearest == summary['hsc_nearest']['euclidean_nearest']
    assert main_cells[other[np.argmax(corr[hsc, other])]] == summary['hsc_nearest']['correlation_nearest']
    passed('distances_and_nearest_neighbors', 'Full Euclidean/correlation matrices and HSC nearest neighbor agree.')

    trees = {f'HVG{n}_z_Ward': linkage(zscore(a[ranking[:n]]).T, method='ward') for n in [100, 500, 1000]}
    trees['HVG1000_z_correlation_average'] = linkage(pdist(z.T, metric='correlation'), method='average')
    trees['HVG1000_input_Ward'] = linkage(a[ranking[:1000]].T, method='ward')
    sensitivity = read('02_clustering_sensitivity.csv')
    assignments = read('02_cluster_assignments.csv')
    for row in sensitivity.itertuples():
        reference = cut_tree(trees['HVG1000_z_Ward'], n_clusters=row.n_clusters).ravel()
        labels = cut_tree(trees[row.method], n_clusters=row.n_clusters).ravel()
        same(adjusted_rand_score(reference, labels), row.ARI_vs_HVG1000_z_Ward)
        saved = assignments[(assignments.method == row.method) & (assignments.n_clusters == row.n_clusters)].set_index('cell_type').loc[main_cells, 'cluster']
        assert adjusted_rand_score(labels, saved) == 1
    passed('clustering_sensitivity', 'All 30 method/k combinations and 630 cell assignments agree with independent tree reconstruction.')

    prism_xml = ET.parse(ROOT / 'prism/HSC_difference_data.pzfx').getroot()
    vals = [float(e.find('Subcolumn/d').text) for e in prism_xml.findall('Table/YColumn')]
    same(vals, rms[focus_i])
    same(pd.read_csv(ROOT / 'prism/HSC_difference_data.tsv', sep='\t').values[0], rms[focus_i])
    passed('Prism_data_only', 'XML/TSV values agree with RMS; this checks data, not native Prism rendering or schema compliance.')

    records = json.loads((ROOT / 'figures/figure_manifest.json').read_text())
    assert len(records) == 7
    for item in records:
        for ext in ['pdf', 'svg', 'png']:
            p = ROOT / 'figures' / f"{item['stem']}.{ext}"
            assert p.stat().st_size > 1000
            if ext == 'pdf':
                pdf = PdfReader(str(p))
                assert len(pdf.pages) == 1
                assert len(pdf.pages[0].extract_text()) > 40
            elif ext == 'svg':
                assert ET.parse(p).getroot().tag.endswith('svg')
            else:
                with Image.open(p) as img:
                    assert min(img.size) > 1000
                    assert abs(img.info['dpi'][0] - 450) < .1
                    img.verify()
    passed('figure_outputs', 'All 21 PDF/SVG/PNG files exist and parse; every PNG has 450 dpi metadata and >1000 px on each side.')

    assert summary['training_started'] is False
    passed('stage_scope', 'Stage-2 summary records no training; stage 3 requires user approval.')
    record = {'passed': True, 'checks': checks, 'generated_utc': datetime.now(timezone.utc).isoformat()}
    (ROOT / 'logs/stage2_validation.json').write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(record, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
