# 核心报告补充材料

## S1．全目标 Lasso 筛查

| 目标 | 定位 | Lasso R² | RMSE | Pearson | Spearman |
|---|---|---|---|---|---|
| common myeloid progenitor | main | 0.1268 | 0.1601 | 0.3569 | 0.2862 |
| erythrocyte | main | 0.6093 | 0.2413 | 0.7809 | 0.6953 |
| platelet | main | 0.4370 | 0.3251 | 0.6620 | 0.5464 |
| macrophage | main | 0.5377 | 0.1979 | 0.7337 | 0.6308 |
| B cell | main | 0.4271 | 0.1595 | 0.6540 | 0.6107 |
| monocyte | main | 0.4835 | 0.1927 | 0.6960 | 0.6285 |
| natural killer cell | main | 0.4702 | 0.1757 | 0.6865 | 0.5801 |
| CD4-positive, alpha-beta T cell | main | 0.4638 | 0.1402 | 0.6813 | 0.6274 |
| CD8-positive, alpha-beta T cell | main | 0.4861 | 0.1630 | 0.6979 | 0.6062 |
| basophil | main | 0.4298 | 0.2416 | 0.6561 | 0.5123 |
| neutrophil | main | 0.5466 | 0.2954 | 0.7397 | 0.6028 |
| myeloid dendritic cell | supplementary | 0.1713 | 0.4654 | 0.4139 | 0.4784 |
| plasmacytoid dendritic cell | main | 0.3930 | 0.1475 | 0.6282 | 0.5470 |
| plasma cell | main | 0.4289 | 0.1620 | 0.6555 | 0.5829 |
| mature NK T cell | main | 0.4090 | 0.1657 | 0.6403 | 0.5569 |
| regulatory T cell | main | 0.3999 | 0.1879 | 0.6327 | 0.5396 |
| classical monocyte | main | 0.4944 | 0.2019 | 0.7037 | 0.6336 |
| non-classical monocyte | main | 0.4725 | 0.1946 | 0.6880 | 0.6047 |
| naive thymus-derived CD4-positive, alpha-beta T cell | main | 0.4010 | 0.1539 | 0.6336 | 0.5921 |
| intermediate monocyte | main | 0.5038 | 0.2078 | 0.7102 | 0.6354 |
| hematopoietic precursor cell | main | 0.3124 | 0.0898 | 0.5599 | 0.5288 |

以上是最低误差模型的留出预测。最低误差与 1-SE 的完整四项指标见[阶段 3 指标](../tables/03_model_metrics.csv)，逐折结果见[分折指标](../tables/03_fold_metrics.csv)，参数及系数见[模型选择](../tables/03_selected_models.csv)与[全扰动系数](../tables/03_all_coefficients_annotated.csv)。

## S2．四目标前 10 个 TF

| 目标 | 名次 | TF | 全数据标准化系数 | 选择频率 |
|---|---|---|---|---|
| B 细胞 | 1 | ZNF22 | 0.06950 | 100% |
| B 细胞 | 2 | DLX3 | 0.03359 | 100% |
| B 细胞 | 3 | RORC | 0.04461 | 100% |
| B 细胞 | 4 | ZNF174 | -0.03991 | 100% |
| B 细胞 | 5 | STAT4 | 0.03659 | 100% |
| B 细胞 | 6 | ZNF296 | 0.04026 | 100% |
| B 细胞 | 7 | KLF3 | 0.03591 | 100% |
| B 细胞 | 8 | HOXB5 | -0.03385 | 100% |
| B 细胞 | 9 | NFYB | -0.03517 | 100% |
| B 细胞 | 10 | ZNF114 | 0.03492 | 100% |
| 浆细胞 | 1 | OVOL2 | 0.07187 | 100% |
| 浆细胞 | 2 | NCOA3 | -0.06393 | 100% |
| 浆细胞 | 3 | ZNF653 | 0.06699 | 100% |
| 浆细胞 | 4 | (DDIT3) | 0.06726 | 100% |
| 浆细胞 | 5 | LHX8 | -0.05911 | 100% |
| 浆细胞 | 6 | ERF | 0.05265 | 100% |
| 浆细胞 | 7 | CLOCK | -0.04856 | 100% |
| 浆细胞 | 8 | ZNF558 | 0.03912 | 100% |
| 浆细胞 | 9 | TFAP2A | 0.04664 | 100% |
| 浆细胞 | 10 | TP53 | -0.03985 | 100% |
| 单核细胞 | 1 | PAX8 | 0.08841 | 100% |
| 单核细胞 | 2 | KLF6 | 0.08430 | 100% |
| 单核细胞 | 3 | TCF7L2 | 0.06733 | 100% |
| 单核细胞 | 4 | TFAP2A | 0.06727 | 100% |
| 单核细胞 | 5 | OVOL2 | 0.06151 | 100% |
| 单核细胞 | 6 | IRF3 | 0.05526 | 100% |
| 单核细胞 | 7 | HIF1A | 0.05968 | 100% |
| 单核细胞 | 8 | PRDM1 | 0.04847 | 100% |
| 单核细胞 | 9 | SIM2 | -0.04987 | 100% |
| 单核细胞 | 10 | SATB2 | 0.04007 | 100% |
| 红细胞 | 1 | MYB | 0.07820 | 100% |
| 红细胞 | 2 | IRF3 | 0.04999 | 100% |
| 红细胞 | 3 | NFYB | -0.05312 | 100% |
| 红细胞 | 4 | PREB | 0.04659 | 100% |
| 红细胞 | 5 | ELK3 | 0.04142 | 100% |
| 红细胞 | 6 | NCOA3 | -0.03989 | 100% |
| 红细胞 | 7 | CREB3 | 0.04000 | 100% |
| 红细胞 | 8 | RARA | -0.03710 | 100% |
| 红细胞 | 9 | OVOL2 | 0.03193 | 100% |
| 红细胞 | 10 | KLF6 | 0.03592 | 100% |

全部 TF 稳定性见[完整表](../tables/04_stability_all_TFs.csv)。稳定优先条目数为 41／43／52／42；规范身份层面两两重叠见[重叠表](../tables/04_stable_TF_overlap.csv)。浆细胞 BCL6 为第 15 名，因直接机制证据进入正文，其自身读出删除不在预定前 10 项检查中。

## S3．完整敏感性与解释口径

- [六类模型及双参数结果](../tables/04_model_metrics.csv)、[逐折结果](../tables/04_fold_metrics.csv)、[全部参数选择](../tables/04_selected_models.csv)。
- [逐项候选证据矩阵](../tables/04_candidate_evidence_matrix.csv)、[相关组](../tables/04_correlation_groups.csv)、[候选相关伙伴](../tables/04_candidate_correlated_partners.csv)。
- [landmark／标准化响应指标](../tables/04_sensitivity_metrics.csv)、[同一 landmark 测试集的全基因模型结果](../tables/04_matched_landmark_reference_metrics.csv)、[自身读出删除](../tables/04_self_readout_results.csv)。
- [三组等规模对照](../tables/04_size_control_metrics.csv)：每组 965 个非 landmark 基因，精确匹配外层数量，内层名单求交；训练与测试身份均改变，三组结果描述本轮抽样范围。

逐基因标准化使用包括目标类型的跨 21 类细胞标准差，形成一个新的响应定义；其回乘原单位的分数在源表中标为描述性复算，未用于主预测排名。条件稳定性从零初始化；先前初始化试算保存在独立探索目录。全部解释与七组图见[阶段 4 报告](04_multimodel_robustness_report.md)。

## S4．方法、文献与复现

[Lasso 方法](03_lasso_methods.md)、[多模型与稳健性方法](04_methods.md)、[环境与重建命令](../environment/README.md)包含完整参数、版本与计算入口。本阶段复用已有结果，未新增训练、重采样或 GPU 拓展。核心报告可由 `scripts/stage5_report.py` 结合既有表格与文献题录重建。

[文献证据表](../tables/05_TF_literature_evidence.csv)、[参考文献题录](../tables/05_references.csv)、[检索与证据范围](05_literature_notes.md)区分实验背景、统计发现及研究假设。完整题录核验记录位于 `data/references/stage5/`。全文访问受限的文献以已核对摘要与原文检索片段支持相应层级的概述，未声称完成所有论文的全文复核。

## S5．图形成果

核心报告采用已审核的 Python 图形，提供 PDF、SVG 与高分辨率 PNG。用户已取消后续 Prism 制图与迁移。已完成的 [Prism 模型比较 PDF](../figures/05_Prism_model_comparison.pdf)、[可编辑项目](../prism/stage5/05_Prism_model_comparison.pzfx)及迁移输入保留为附加历史成果，后续分析和报告采用现有绘图流程。

## S6．独立审核

[阶段 5 独立审核](05_independent_review.md)覆盖主线、数值、引用与作业要求；[阶段 4 审核](04_independent_review.md)保存模型结果的独立数值核验。用户已批准本阶段并要求进入第 6 阶段。
