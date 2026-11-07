"""Stage 1: reproducible QC, reference annotation and aligned model inputs.

This script prepares data only. No EDA, model fitting or data-driven feature
selection is performed. Original expression labels remain the alignment key.
"""
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
HW2 = ROOT.parent
for name in ["tables", "reports", "environment", "logs", "data/processed", ".cache/matplotlib", ".cache/fontconfig"]:
    (ROOT / name).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / ".cache"))
import h5py
import numpy as np
import pandas as pd

CONFIG = json.loads((ROOT / "config/analysis_config.json").read_text())
REF = ROOT / "data/references"
TABLES = ROOT / "tables"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def save_table(name, frame):
    frame.to_csv(TABLES / name, index=False, encoding="utf-8-sig", float_format="%.12g")


def matrix_stats(name, frame):
    a = frame.to_numpy()
    q = np.quantile(a, [0, .01, .25, .5, .75, .99, 1])
    return {"dataset": name, "rows": frame.shape[0], "columns": frame.shape[1],
            "duplicate_rows_to_remove": int(frame.index.duplicated().sum()),
            "duplicate_columns": int(frame.columns.duplicated().sum()),
            "missing": int(np.isnan(a).sum()), "infinite": int(np.isinf(a).sum()),
            "zero_fraction": float((a == 0).mean()),
            **dict(zip(["minimum", "p01", "p25", "median", "p75", "p99", "maximum"], map(float, q)))}


def duplicates(name, frame):
    positions = defaultdict(list)
    for i, symbol in enumerate(frame.index):
        positions[symbol].append(i)
    records = []
    for symbol, ids in positions.items():
        if len(ids) < 2:
            continue
        ref = frame.iloc[ids[0]].to_numpy()
        for i in ids:
            records.append({"dataset": name, "gene_symbol": symbol, "data_row_1based": i + 1,
                            "action": "keep_first" if i == ids[0] else "drop_duplicate",
                            "max_abs_difference_from_first": float(np.max(np.abs(frame.iloc[i].to_numpy() - ref)))})
    return records


def markdown_table(frame):
    def cell(x):
        return str(x).replace("|", "/").replace("\n", " ")
    head = "| " + " | ".join(map(cell, frame.columns)) + " |\n"
    sep = "| " + " | ".join(["---"] * len(frame.columns)) + " |\n"
    return head + sep + "\n".join("| " + " | ".join(map(cell, row)) + " |" for row in frame.itertuples(index=False, name=None))


def main():
    blood_path = HW2 / "07-高维回归/数据/blood mean expr.csv"
    l1000_path = HW2 / "L1000 MCF7_OE gene expression.h5"
    inputs = []
    for role, path in [("blood", blood_path), ("L1000", l1000_path),
                       ("assignment", HW2 / "第二次作业2026.pdf"),
                       ("course_demo", HW2 / "07-高维回归/数据/Tabula blood demo3.ipynb")]:
        inputs.append({"role": role, "absolute_path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})
    save_json(ROOT / "data/input_manifest.json", inputs)
    blood_raw = pd.read_csv(blood_path, index_col=0)
    l1000_raw = pd.read_hdf(l1000_path, key="expression")
    for frame in [blood_raw, l1000_raw]:
        assert frame.index.notna().all() and frame.columns.notna().all()
        assert all(isinstance(x, str) and x and x == x.strip() for x in frame.index)
        assert frame.columns.is_unique
        assert np.isfinite(frame.to_numpy()).all()
    assert (blood_raw.to_numpy() >= 0).all()
    save_table("01_matrix_qc.csv", pd.DataFrame([matrix_stats("blood_raw", blood_raw), matrix_stats("L1000_raw", l1000_raw)]))
    save_table("01_duplicate_gene_audit.csv", pd.DataFrame(duplicates("blood", blood_raw) + duplicates("L1000", l1000_raw)))
    blood = blood_raw.loc[~blood_raw.index.duplicated(keep="first")].copy()
    l1000 = l1000_raw.loc[~l1000_raw.index.duplicated(keep="first")].copy()
    shared = blood.index[blood.index.isin(l1000.index)]
    all_zero = (blood.loc[shared] == 0).all(axis=1)
    genes = shared[~all_zero]
    X = l1000.loc[genes].copy()
    valid_predictors = X.nunique(dropna=False) > 1
    save_table("01_predictor_qc.csv", pd.DataFrame({
        "perturbation_label": X.columns,
        "minimum": X.min().values, "maximum": X.max().values,
        "mean": X.mean().values, "sd": X.std().values,
        "zero_fraction": (X == 0).mean().values,
        "included": valid_predictors.values,
        "reason": np.where(valid_predictors, "valid", "constant_across_model_genes")}))
    X = X.loc[:, valid_predictors]
    hsc = CONFIG["reference_cell"]
    targets = [c for c in blood.columns if c != hsc]
    Y = blood.loc[genes, targets].subtract(blood.loc[genes, hsc], axis=0)
    assert X.index.equals(Y.index)
    assert all(c in targets for c in CONFIG["focus_targets"])

    hgnc = pd.read_csv(REF / "hgnc_complete_set.txt", sep="\t", dtype=str, keep_default_na=False)
    hgnc = hgnc.loc[hgnc.status.eq("Approved")].set_index("symbol", drop=False)
    aliases = defaultdict(set)
    previous = defaultdict(set)
    for row in hgnc.itertuples(index=False):
        for old in row.prev_symbol.split("|"):
            if old: previous[old].add(row.symbol)
        for old in row.alias_symbol.split("|"):
            if old: aliases[old].add(row.symbol)

    def resolve(symbol):
        if symbol in hgnc.index: return symbol, "approved_exact"
        for method, lookup in [("previous_symbol_unique", previous), ("alias_unique", aliases)]:
            hits = lookup.get(symbol, set())
            if len(hits) == 1: return next(iter(hits)), method
            if len(hits) > 1: return "", "ambiguous_" + method.split("_unique")[0]
        return "", "unresolved"

    tf_db = pd.read_csv(REF / "HumanTFs_v1.01.csv", keep_default_na=False)
    tf_records = defaultdict(list)
    for _, row in tf_db.iterrows():
        approved, method = resolve(row["HGNC symbol"])
        if approved: tf_records[approved].append(row.to_dict())

    perturbations = pd.read_csv(REF / "GSE92742_Broad_LINCS_pert_info.txt.gz", sep="\t", dtype=str, keep_default_na=False)
    oe_names = set(perturbations.loc[perturbations.pert_type.eq("trt_oe"), "pert_iname"])
    annotations = []
    for label in X.columns:
        name = label
        parenthesized = name.startswith("(") and name.endswith(")")
        if parenthesized: name = name[1:-1]
        mutated = name.endswith("_MUT")
        if mutated: name = name[:-4]
        approved, method = resolve(name)
        records = tf_records.get(approved, [])
        tf_yes = [r for r in records if r["Is TF?"] == "Yes"]
        # Exact curated symbols remain informative when HGNC has no unambiguous mapping.
        if not approved:
            exact = tf_db.loc[tf_db["HGNC symbol"].eq(name)]
            tf_yes = exact.loc[exact["Is TF?"].eq("Yes")].to_dict("records")
        group = hgnc.loc[approved, "locus_group"] if approved else ""
        locus_type = hgnc.loc[approved, "locus_type"] if approved else ""
        if tf_yes: category = "TF"
        elif group == "protein-coding gene": category = "other_protein_coding"
        elif group == "non-coding RNA": category = "noncoding_RNA"
        elif group == "pseudogene": category = "pseudogene"
        elif group: category = "other_locus"
        else: category = "unresolved"
        annotations.append({"perturbation_label": label, "lookup_label": name,
                            "approved_symbol": approved, "mapping_method": method,
                            "parenthesized_label": parenthesized, "mutation_tag": mutated,
                            "construct_status": "MUT_label" if mutated else "no_MUT_label",
                            "locus_group": group, "locus_type": locus_type,
                            "annotation_class": category, "is_TF": bool(tf_yes),
                            "primary_TF_eligible": bool(tf_yes) and not mutated,
                            "TF_catalog_symbols": "|".join(sorted({r["HGNC symbol"] for r in tf_yes})),
                            "TF_assessment": "|".join(sorted({r["TF assessment"] for r in tf_yes})),
                            "DBD": "|".join(sorted({r["DBD"] for r in tf_yes})),
                            "GSE92742_OE_name_match": name in oe_names or approved in oe_names,
                            "name_match_scope": "annotation_only_not_signature_provenance"})
    ann = pd.DataFrame(annotations)
    save_table("01_perturbation_annotations.csv", ann)
    save_table("01_TF_perturbations.csv", ann.loc[ann.is_TF])
    duplicate_tf_symbols = ann.loc[ann.is_TF, "approved_symbol"].duplicated(keep=False)
    save_table("01_TF_shared_identity_labels.csv", ann.loc[ann.is_TF].loc[duplicate_tf_symbols])
    save_table("01_unresolved_perturbations.csv", ann.loc[ann.annotation_class.eq("unresolved")])
    category_counts = ann.groupby(["annotation_class", "construct_status"]).size().rename("n").reset_index()
    save_table("01_perturbation_categories.csv", category_counts)

    gi = pd.read_csv(REF / "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="\t", dtype={"pr_gene_symbol":str})
    gi_groups = {key: value for key, value in gi.groupby("pr_gene_symbol", sort=False)}
    gi_by_entrez = {str(key): value for key, value in gi.groupby("pr_gene_id", sort=False)}
    gi_canonical = defaultdict(list)
    for i, name in gi.pr_gene_symbol.items():
        approved, _ = resolve(name)
        if approved: gi_canonical[approved].append(i)
    gene_records = []
    for name in shared:
        approved, method = resolve(name)
        g = gi_groups.get(name)
        reference_method = "exact_original_symbol" if g is not None else "unresolved"
        if g is None and approved:
            entrez = hgnc.loc[approved, 'entrez_id']
            if entrez in gi_by_entrez:
                g = gi_by_entrez[entrez]
                reference_method = "HGNC_Entrez_ID"
            elif approved in gi_canonical:
                g = gi.loc[gi_canonical[approved]]
                reference_method = "HGNC_symbol_history"
        available = g is not None
        lm_values = sorted(g.pr_is_lm.unique().tolist()) if available else []
        bing_values = sorted(g.pr_is_bing.unique().tolist()) if available else []
        # Duplicate source symbols with conflicting measurement class stay unresolved.
        lm = int(lm_values[0]) if len(lm_values) == 1 else -1
        bing = int(bing_values[0]) if len(bing_values) == 1 else -1
        gene_records.append({"gene_symbol": name, "approved_symbol_annotation_only": approved,
                             "hgnc_mapping_method": method,
                             "in_model": bool(not all_zero.loc[name]),
                             "exclusion_reason": "all_blood_types_zero" if all_zero.loc[name] else "",
                             "GSE92742_gene_info_match": available,
                             "reference_mapping_method": reference_method,
                             "source_gene_symbols": "|".join(g.pr_gene_symbol.tolist()) if available else "",
                             "source_gene_ids": "|".join(map(str, g.pr_gene_id.tolist())) if available else "",
                             "is_landmark": lm, "is_BING": bing,
                             "reference_symbol_ambiguous": len(lm_values) > 1 or len(bing_values) > 1})
    gene_ann = pd.DataFrame(gene_records)
    save_table("01_shared_gene_annotations.csv", gene_ann)
    save_table("01_excluded_shared_genes.csv", gene_ann.loc[~gene_ann.in_model])
    model_gene_ann = gene_ann.loc[gene_ann.in_model].set_index("gene_symbol").loc[genes]
    assert (model_gene_ann.is_landmark >= 0).all(), "Resolve landmark mapping before freezing inputs"
    save_table("01_landmark_model_genes.csv", model_gene_ann.loc[model_gene_ann.is_landmark.eq(1)].reset_index())

    readout_symbols = {resolve(n)[0] or n for n in genes}
    all_blood_symbols = {resolve(n)[0] or n for n in blood.index}
    tf_coverage = []
    for _, row in tf_db.loc[tf_db["Is TF?"].eq("Yes")].iterrows():
        name = row["HGNC symbol"]
        approved, method = resolve(name)
        canonical = approved or name
        found = ann.loc[ann.is_TF & (ann.approved_symbol.eq(canonical) | ann.lookup_label.eq(name))]
        tf_coverage.append({"catalog_symbol": name, "approved_symbol": approved,
                            "mapping_method": method, "DBD": row["DBD"], "TF_assessment": row["TF assessment"],
                            "in_blood_readout": canonical in all_blood_symbols,
                            "in_model_readout": canonical in readout_symbols,
                            "available_perturbation_count": len(found),
                            "eligible_perturbation_count": int(found.primary_TF_eligible.sum()),
                            "perturbation_labels": "|".join(found.perturbation_label)})
    coverage = pd.DataFrame(tf_coverage)
    save_table("01_TF_catalog_coverage.csv", coverage)
    marker_names = ['EBF1','PRDM1','BCL6','GATA2','GATA3','FOXP3','MYB','FLI1','STAT5A',
                    'GATA1','KLF1','TAL1','SPI1','CEBPA','CEBPB','PAX5','XBP1','IRF4','IRF8','RUNX1','IKZF1']
    marker_panel = coverage.loc[coverage.catalog_symbol.isin(marker_names)].copy()
    marker_panel['panel_order'] = marker_panel.catalog_symbol.map({n:i for i,n in enumerate(marker_names)})
    marker_panel = marker_panel.sort_values('panel_order').drop(columns='panel_order')
    save_table("01_lineage_TF_coverage.csv", marker_panel)

    cell_rows = []
    for name in blood.columns:
        v = blood[name]
        role = "reference" if name == hsc else "supplementary_target" if name == CONFIG['supplementary_target'] else "main_target"
        cell_rows.append({"cell_type": name, "role": role, "focus": name in CONFIG['focus_targets'],
                          "n_genes": len(v), "minimum": float(v.min()), "maximum": float(v.max()),
                          "mean": float(v.mean()), "median": float(v.median()), "sd": float(v.std()),
                          "zero_count": int((v == 0).sum()), "zero_fraction": float((v == 0).mean()),
                          "n_nonzero_target_minus_HSC_model_genes": int((Y[name] != 0).sum()) if name in Y else 0})
    cells = pd.DataFrame(cell_rows)
    save_table("01_cell_quality.csv", cells)
    save_table("01_target_manifest.csv", cells.loc[cells.role.ne('reference'), ['cell_type','role','focus']])

    path = ROOT / "data/processed/analysis_matrices.h5"
    string_type = h5py.string_dtype('utf-8')
    with h5py.File(path, 'w') as f:
        f.attrs.update(schema_version=1, response='input_scale_target_minus_HSC',
                       observations='genes', alignment='exact_original_symbol_in_blood_order',
                       generated_utc=datetime.now(timezone.utc).isoformat(),
                       blood_unit='source_unit_unrecorded_preserved_input_scale',
                       l1000_processing_level='not_identified_from_local_HDF_metadata')
        def strings(key, values): f.create_dataset(key, data=list(values), dtype=string_type)
        f.create_dataset('blood/expression', data=blood.to_numpy(), compression='gzip', compression_opts=1, fletcher32=True)
        strings('blood/genes', blood.index)
        strings('blood/cell_types', blood.columns)
        f.create_dataset('model/X', data=X.to_numpy(), compression='gzip', compression_opts=1, fletcher32=True)
        f.create_dataset('model/Y', data=Y.to_numpy(), compression='gzip', compression_opts=1, fletcher32=True)
        strings('model/genes', genes)
        strings('model/perturbations', X.columns)
        strings('model/target_cells', Y.columns)
        f.create_dataset('model/is_landmark', data=model_gene_ann.is_landmark.eq(1).to_numpy())
        f.create_dataset('model/is_TF', data=ann.is_TF.to_numpy())
        f.create_dataset('model/primary_TF_eligible', data=ann.primary_TF_eligible.to_numpy())
        f.create_dataset('model/is_main_target', data=np.array([c != CONFIG['supplementary_target'] for c in Y.columns]))
    # Round-trip comparison checks that labels, masks and numerical arrays survived storage.
    with h5py.File(path, 'r') as f:
        assert list(f['model/genes'].asstr()[:]) == list(genes)
        assert list(f['model/perturbations'].asstr()[:]) == list(X.columns)
        assert list(f['model/target_cells'].asstr()[:]) == targets
        assert np.array_equal(f['model/X'][:], X.to_numpy())
        assert np.array_equal(f['model/Y'][:], Y.to_numpy())
        assert np.array_equal(f['blood/expression'][:], blood.to_numpy())

    packages = ['numpy','pandas','scipy','scikit-learn','h5py','tables','matplotlib','seaborn']
    versions = {name: metadata.version(name) for name in packages}
    environment = {'python_executable':sys.executable, 'python_version':platform.python_version(),
                   'platform':platform.platform(), 'packages':versions,
                   'strategy':'reuse_verified_runtime_pin_dependencies_no_system_environment_mutation'}
    save_json(ROOT / 'environment/local_environment.json', environment)
    (ROOT / 'environment/requirements-stage1.txt').write_text('\n'.join(f'{name}=={version}' for name,version in versions.items())+'\n')
    distribution_versions = sorted({(d.metadata['Name'],d.version) for d in metadata.distributions() if d.metadata['Name']})
    (ROOT / 'environment/package_inventory.txt').write_text('\n'.join(f'{n}=={v}' for n,v in distribution_versions)+'\n')

    summary = {
        'blood_raw_shape':list(blood_raw.shape), 'blood_unique_shape':list(blood.shape),
        'l1000_raw_shape':list(l1000_raw.shape), 'l1000_unique_shape':list(l1000.shape),
        'shared_gene_count':len(shared), 'excluded_all_blood_zero_count':int(all_zero.sum()),
        'model_X_shape':list(X.shape), 'model_Y_shape':list(Y.shape),
        'landmark_model_gene_count':int(model_gene_ann.is_landmark.eq(1).sum()),
        'gene_reference_mapping_counts':{str(k):int(v) for k,v in gene_ann.reference_mapping_method.value_counts().items()},
        'reference_landmark_count':int(gi.pr_is_lm.eq(1).sum()),
        'BING_model_gene_count':int(model_gene_ann.is_BING.eq(1).sum()),
        'TF_perturbation_count':int(ann.is_TF.sum()),
        'primary_TF_eligible_count':int(ann.primary_TF_eligible.sum()),
        'unique_TF_catalog_coverage':int(coverage.eligible_perturbation_count.gt(0).sum()),
        'mutation_label_count':int(ann.mutation_tag.sum()),
        'unresolved_perturbation_count':int(ann.annotation_class.eq('unresolved').sum()),
        'supplementary_target_zero_fraction':float(cells.set_index('cell_type').loc[CONFIG['supplementary_target'],'zero_fraction']),
        'model_matrix_file':str(path.relative_to(ROOT)), 'model_matrix_bytes':path.stat().st_size,
        'model_matrix_sha256':sha256(path),
        'model_training_started':False, 'next_stage':'02_EDA_requires_user_review',
        'checks': {'finite_values':True, 'unique_gene_and_feature_labels':True, 'deterministic_alignment':True,
                   'response_equals_target_minus_HSC':True, 'TF_mask_matches_predictor_order':True,
                   'landmark_mapping_resolved':True, 'HDF5_roundtrip_values_and_labels':True}
    }
    save_json(ROOT / 'logs/stage1_summary.json', summary)
    save_json(ROOT / 'stage_status.json', {'stage':1,'status':'awaiting_user_review','completed_utc':datetime.now(timezone.utc).isoformat(), 'next_stage':2,'next_stage_started':False})

    available_markers = marker_panel.loc[marker_panel.eligible_perturbation_count.gt(0),'catalog_symbol'].tolist()
    missing_markers = marker_panel.loc[marker_panel.eligible_perturbation_count.eq(0),'catalog_symbol'].tolist()
    category_md = markdown_table(category_counts)
    focus_md = markdown_table(cells.loc[cells.cell_type.isin(CONFIG['focus_targets']+[hsc,CONFIG['supplementary_target']]), ['cell_type','role','zero_fraction','median','maximum']].round(4))
    report = f'''# 第 1 阶段：环境与数据定稿

本阶段已完成数据质控、确定性对齐、TF/基因类别注释及 landmark 注释，形成后续 EDA 和建模的统一输入。模型训练尚未启动。

## 主要结果

| 项目 | 结果 |
|---|---|
| 血液原始矩阵 | {blood_raw.shape[0]:,} 基因行 × {blood_raw.shape[1]} 细胞类型 |
| 血液去重矩阵 | {blood.shape[0]:,} 基因 × {blood.shape[1]} 细胞类型 |
| L1000 原始矩阵 | {l1000_raw.shape[0]:,} 基因行 × {l1000_raw.shape[1]:,} 扰动 |
| 共同独立基因 | {len(shared):,} |
| 全部血液类型为零而移除的共同基因 | {int(all_zero.sum())} |
| 正式 X | {X.shape[0]:,} 基因 × {X.shape[1]:,} 扰动 |
| 正式 Y | {Y.shape[0]:,} 基因 × {Y.shape[1]} 目标 |
| 直接测量 landmark 共同基因 | {summary['landmark_model_gene_count']} / 官方 {summary['reference_landmark_count']} |
| TF 扰动 / 无 MUT 标签的 TF 扰动 | {summary['TF_perturbation_count']} / {summary['primary_TF_eligible_count']} |
| 原始数值缺失与无穷值 | 两份矩阵均为 0 |

## 数据处理与表达尺度

重复基因保留原始首条，并保存每个重复行与首条的最大数值差。血液表删除 {int(blood_raw.index.duplicated().sum())} 个重复行，L1000 删除 {int(l1000_raw.index.duplicated().sum())} 个重复行。原始文件保持完整；输入文件摘要见 `data/input_manifest.json`。

共同基因按血液数据的原始顺序排列，X 和 Y 使用完全相同的基因顺序。HGNC 别名用于注释，矩阵仍按原始基因符号精确对齐。

血液输入范围为 {blood_raw.to_numpy().min():.4g}–{blood_raw.to_numpy().max():.6g}，保留当前尺度，响应定义为目标细胞表达值减去 HSC 表达值。上游单位与变换未包含在课程文件中，因此报告称其为“输入尺度表达差值”。L1000 范围为 {l1000_raw.to_numpy().min():.4g}–{l1000_raw.to_numpy().max():.6g}，保留正负分数；HDF 仅保存表达矩阵，具体处理层级、重复、时间和剂量需要源数据记录才能确定。

GEO 官方说明区分直接测量与推断基因，并提供 gene_info 注释；本次用该注释构建 landmark 子集。在 {len(shared):,} 个共同基因中，{summary['gene_reference_mapping_counts'].get('exact_original_symbol',0):,} 个直接按原符号对应，另 {summary['gene_reference_mapping_counts'].get('HGNC_Entrez_ID',0)} 个通过 HGNC 的 Entrez ID 对应历史名称。例如 MCUB 对应旧名 CCDC109B、ELOA 对应 TCEB3。每个映射均保存源符号与源 ID，矩阵标签和数值保持一致。课程 HDF 的具体实验来源仍记录为未定。[GEO GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)

## 细胞类型与分析职责

{focus_md}

HSC 为参考；B 细胞、浆细胞、单核细胞和红细胞为重点。髓系树突状细胞零值率为 {summary['supplementary_target_zero_fraction']:.2%}，其模型结果将单列。其余 20 个目标用于主要跨类型候选比较。

## 扰动身份与 TF 覆盖

{category_md}

TF 判定来自 Human TFs v1.01 的 `Is TF?=Yes`，其中包含已知及可能的人类 TF。HGNC approved、唯一历史符号及唯一别名按顺序匹配；歧义保留待核实。括号只作为原始命名标志；`_MUT` 单独记录。无 MUT 标签表示命名中没有该标记，其构建体细节由后续来源核验补充。[Human TFs](https://humantfs.ccbr.utoronto.ca/download.php)

预先指定的谱系参照中，可用过表达扰动包括：{', '.join(available_markers)}。缺少对应扰动的参照包括：{', '.join(missing_markers)}。完整目录同时记录表达读出是否覆盖这些基因；表达读出与可选扰动分别列示。

291 个 TF 扰动条目对应 290 个目录 TF：`DDIT3` 与 `(DDIT3)` 是独立的两列，保留各自表达模式，并在后续解释中归到同一 TF 身份。EBF1 具有可用扰动列，其自身基因不在本次共同读出集合中；该扰动可由其他基因的响应参与建模，后续自身读出排除检查对它标为不适用。

类别为 unresolved 的条目保留在完整扰动字典中，进入全部扰动模型；它们在 TF 主表中的资格保持未确认。GSE92742 的 OE 名称命中仅作注释交叉核验，未据此推定课程 HDF 中每一列的实验来源。

## 运行环境与复现

本地使用 Python {platform.python_version()} 与已验证的运行环境；核心包精确版本见 `environment/requirements-stage1.txt`。数据准备脚本：`scripts/prepare_data.py`；官方参考文件与摘要：`data/references/source_manifest.json`。

已通过数值有限性、唯一标签、确定性对齐、响应定义、TF/landmark 掩码顺序与 HDF5 写后重读检查。正式输入保存为 `data/processed/analysis_matrices.h5`；该文件含全部去重血液表达、X/Y、基因/扰动/目标顺序及筛选掩码。

## 审核要点

本阶段沿用已批准的输入尺度、首条去重、全扰动建模及 mDC 附加分析规则。请审核最终输入范围、注释覆盖和数据处理说明。审核通过后进入第 2 阶段，制作热图、PCA、聚类对照和图形样板。
'''
    (ROOT / 'reports/01_data_qc.md').write_text(report)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
