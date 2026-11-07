# 第 1 阶段：环境与数据定稿

本阶段已完成数据质控、确定性对齐、TF/基因类别注释及 landmark 注释，形成后续 EDA 和建模的统一输入。模型训练尚未启动。

## 主要结果

| 项目 | 结果 |
|---|---|
| 血液原始矩阵 | 20,001 基因行 × 22 细胞类型 |
| 血液去重矩阵 | 19,994 基因 × 22 细胞类型 |
| L1000 原始矩阵 | 12,328 基因行 × 2,717 扰动 |
| 共同独立基因 | 11,940 |
| 全部血液类型为零而移除的共同基因 | 18 |
| 正式 X | 11,922 基因 × 2,717 扰动 |
| 正式 Y | 11,922 基因 × 21 目标 |
| 直接测量 landmark 共同基因 | 965 / 官方 978 |
| TF 扰动 / 无 MUT 标签的 TF 扰动 | 291 / 291 |
| 原始数值缺失与无穷值 | 两份矩阵均为 0 |

## 数据处理与表达尺度

重复基因保留原始首条，并保存每个重复行与首条的最大数值差。血液表删除 7 个重复行，L1000 删除 1 个重复行。原始文件保持完整；输入文件摘要见 `data/input_manifest.json`。

共同基因按血液数据的原始顺序排列，X 和 Y 使用完全相同的基因顺序。HGNC 别名用于注释，矩阵仍按原始基因符号精确对齐。

血液输入范围为 0–8.77847，保留当前尺度，响应定义为目标细胞表达值减去 HSC 表达值。上游单位与变换未包含在课程文件中，因此报告称其为“输入尺度表达差值”。L1000 范围为 -10–10，保留正负分数；HDF 仅保存表达矩阵，具体处理层级、重复、时间和剂量需要源数据记录才能确定。

GEO 官方说明区分直接测量与推断基因，并提供 gene_info 注释；本次用该注释构建 landmark 子集。在 11,940 个共同基因中，11,693 个直接按原符号对应，另 247 个通过 HGNC 的 Entrez ID 对应历史名称。例如 MCUB 对应旧名 CCDC109B、ELOA 对应 TCEB3。每个映射均保存源符号与源 ID，矩阵标签和数值保持一致。课程 HDF 的具体实验来源仍记录为未定。[GEO GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)

## 细胞类型与分析职责

| cell_type | role | zero_fraction | median | maximum |
| --- | --- | --- | --- | --- |
| hematopoietic stem cell | reference | 0.2921 | 0.041 | 4.5651 |
| erythrocyte | main_target | 0.2015 | 0.0006 | 8.7785 |
| B cell | main_target | 0.1081 | 0.0215 | 4.1891 |
| monocyte | main_target | 0.0265 | 0.017 | 5.0442 |
| myeloid dendritic cell | supplementary_target | 0.94 | 0.0 | 6.3168 |
| plasma cell | main_target | 0.1596 | 0.0226 | 4.4132 |

HSC 为参考；B 细胞、浆细胞、单核细胞和红细胞为重点。髓系树突状细胞零值率为 94.00%，其模型结果将单列。其余 20 个目标用于主要跨类型候选比较。

## 扰动身份与 TF 覆盖

| annotation_class | construct_status | n |
| --- | --- | --- |
| TF | no_MUT_label | 291 |
| noncoding_RNA | MUT_label | 7 |
| noncoding_RNA | no_MUT_label | 179 |
| other_locus | no_MUT_label | 7 |
| other_protein_coding | MUT_label | 2 |
| other_protein_coding | no_MUT_label | 1908 |
| pseudogene | no_MUT_label | 3 |
| unresolved | MUT_label | 3 |
| unresolved | no_MUT_label | 317 |

TF 判定来自 Human TFs v1.01 的 `Is TF?=Yes`，其中包含已知及可能的人类 TF。HGNC approved、唯一历史符号及唯一别名按顺序匹配；歧义保留待核实。括号只作为原始命名标志；`_MUT` 单独记录。无 MUT 标签表示命名中没有该标记，其构建体细节由后续来源核验补充。[Human TFs](https://humantfs.ccbr.utoronto.ca/download.php)

预先指定的谱系参照中，可用过表达扰动包括：EBF1, PRDM1, BCL6, GATA2, GATA3, FOXP3, MYB, FLI1, STAT5A。缺少对应扰动的参照包括：GATA1, KLF1, TAL1, SPI1, CEBPA, CEBPB, PAX5, XBP1, IRF4, IRF8, RUNX1, IKZF1。完整目录同时记录表达读出是否覆盖这些基因；表达读出与可选扰动分别列示。

291 个 TF 扰动条目对应 290 个目录 TF：`DDIT3` 与 `(DDIT3)` 是独立的两列，保留各自表达模式，并在后续解释中归到同一 TF 身份。EBF1 具有可用扰动列，其自身基因不在本次共同读出集合中；该扰动可由其他基因的响应参与建模，后续自身读出排除检查对它标为不适用。

类别为 unresolved 的条目保留在完整扰动字典中，进入全部扰动模型；它们在 TF 主表中的资格保持未确认。GSE92742 的 OE 名称命中仅作注释交叉核验，未据此推定课程 HDF 中每一列的实验来源。

## 运行环境与复现

本地使用 Python 3.12.7 与已验证的运行环境；核心包精确版本见 `environment/requirements-stage1.txt`。数据准备脚本：`scripts/prepare_data.py`；官方参考文件与摘要：`data/references/source_manifest.json`。

已通过数值有限性、唯一标签、确定性对齐、响应定义、TF/landmark 掩码顺序与 HDF5 写后重读检查。正式输入保存为 `data/processed/analysis_matrices.h5`；该文件含全部去重血液表达、X/Y、基因/扰动/目标顺序及筛选掩码。

## 审核要点

本阶段沿用已批准的输入尺度、首条去重、全扰动建模及 mDC 附加分析规则。请审核最终输入范围、注释覆盖和数据处理说明。审核通过后进入第 2 阶段，制作热图、PCA、聚类对照和图形样板。
