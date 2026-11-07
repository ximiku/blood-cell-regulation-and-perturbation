# 血液细胞调控候选与扰动响应预测

**Blood-cell regulation and perturbation prediction** · 基因组学数据统计分析课程项目

以 HSC 为参照，通过 L1000 扰动表达字典重建血液细胞差异，筛选稳定的转录因子候选；进一步在 Norman 2019 K562 CRISPRa 数据上比较线性基线、GEARS、scGPT 与 scFoundation＋GEARS。

## 阅读入口

- **[独立分析报告（PDF）](submission/分析报告-刘家希.pdf)**：论文式正文、六组综合图、方法、参考文献与完整附录。
- **[工作记录](submission/工作记录.md)**：任务过程、主要调整与核查记录。
- **[可编辑 LaTeX 工程](manuscript/)**：正文、附录、文献和图形；使用 XeLaTeX 编译。
- **[完整 solution 归档](https://github.com/ximiku/blood-cell-regulation-and-perturbation/releases/tag/solution-archive-2026-11-08)**：包含数据、预测、模型检查点、缓存、环境包及全部研究成果。

## 研究主线

全类型筛查中，红细胞表达差异的 Lasso 留出基因重建表现最突出（R² = 0.609）。MYB、负向 BCL6、KLF6／PRDM1 将稳定候选联系到红系、浆细胞与髓系表达程序。外部响应预测中，共同基因的线性基线表现较强，高变基因上的 GEARS 平均误差最低；完整比较与细胞背景讨论见报告。

## 目录与复现

| 目录 | 内容 |
|---|---|
| `scripts/`、`config/` | 分析与绘图代码、配置 |
| `tables/`、`results/` | 指标、候选、诊断与实验结果；大数组和检查点见 Release |
| `figures/`、`manuscript/` | 图形与论文工程 |
| `reports/`、`logs/` | 各阶段交付与运行记录 |
| `environment/`、`vendor/` | 环境记录与第三方代码 |
| `submission/` | 最终提交稿与工作记录 |

统计分析及 GPU 实验按阶段配置运行；环境版本见 `environment/`，完整流程见下方原始项目说明与 `WORK_LOG.md`。历史脚本中的个人目录和 GPU 接入地址已替换为占位符，需要按自己的环境设置。报告编译不需要 GPU。

普通 Git 仓库提供可浏览成果；`.gitignore` 列出的数据、检查点与缓存均保存在 Release 完整归档中。完整归档按 64 MiB 分卷（初次 1 GiB 上传遇到 HTTP 408 后调整），下载同一 Release 的所有 `solution.tar.gz.part*`、`SHA256SUMS` 和 `file_manifest.json` 后，在新目录执行：

```bash
shasum -a 256 -c SHA256SUMS
cat solution.tar.gz.part* | tar -xzf -
```

归档还原为 `solution/`。`file_manifest.json` 给出逐文件字节数和 SHA-256，`verification.json` 记录重新读取压缩包后的逐文件验证结果。

## 公开副本说明

原始本地成果保持不变。公开副本替换本机用户路径和历史 GPU 接入地址，省略 Finder 的 `.DS_Store` 元数据；变更见 `PUBLICATION_NOTES.json`。归档清单对应公开副本；历史阶段清单记录的是当时的计算成果。保留第三方许可证和数据来源，本仓库未对第三方数据、模型及依赖授予额外许可。

---

## 原项目详细说明

## 最新提交稿（2026-11-08）

- `submission/分析报告-刘家希.pdf`：35 页独立报告，正文 13 页、六组综合图及完整附录。
- `submission/工作记录.md`：执行记录与本人待补部分。
- `manuscript/Homework2_LaTeX_source.zip`：可独立编译的编辑工程。

# Homework 2 分析成果

**当前状态：第1–8阶段分析均已完成，增强报告与本地归档供用户审核。九项GPU正式实验和每目标20次重新调参补充全部完成，最终核验通过。Prism制图及后续独立subagent审核已按用户要求取消。**

## 阅读入口

- [一页结果说明](reports/08_results_brief.md)：先看主要结论和科研价值。
- [增强分析报告](reports/08_enhanced_report.md)：完整作业与外部模型拓展。
- [GPU实验完整结果](reports/07_gpu_benchmark_report.md)：全部比较、分层与案例。
- [最终归档与使用说明](reports/08_archive_handoff.md)：核验、复现及远端权重清单。

- [第 6 阶段 GPU 试点报告](reports/06_gpu_pilot_report.md)：三模型实测、统一数据和第 7 步运行设置。
- [模型来源与预训练记录](reports/06_pretraining_audit.md)：固定源码、权重、输入适配和复现入口。

- [核心分析报告](reports/05_core_report.md)：保留核心版本；围绕 Lasso 表达重建、代表 TF、生物学解释与研究连接组织。
- [核心报告补充材料](reports/05_core_supplement.md)：全目标结果、40 个前 10 候选及完整分析入口。
- [阶段 5 独立审核](reports/05_independent_review.md)：关键数值、作业覆盖、文献解释与主线审核。
- [TF 文献证据表](tables/05_TF_literature_evidence.csv)及[检索说明](reports/05_literature_notes.md)：统计证据、实验背景、支持范围与验证方向。
- [阶段 4 多模型与稳健性报告](reports/04_multimodel_robustness_report.md)：模型比较、稳定 TF、相关组、读出与尺度敏感性，以及补充等规模对照。
- [阶段 4 独立审核](reports/04_independent_review.md)：原始数据独立复算、图文审核与供用户参考的建议。
- [阶段 4 方法说明](reports/04_methods.md)：新增模型、共享验证、条件稳定性与敏感性设计。
- [阶段 3 Lasso 报告](reports/03_lasso_report.md)：21 个目标的折外表现、四类重点目标、系数路径与 TF 初筛。
- [阶段 3 独立审核](reports/03_independent_review.md)：独立重算、复现检查及供用户参考的审核建议。
- [阶段 3 方法说明](reports/03_lasso_methods.md)：嵌套交叉验证、1-SE 选择与 GPU 双精度求解。
- [阶段 1 质控报告](reports/01_data_qc.md)：最终矩阵、表达尺度、细胞类型、扰动注释和审核要点。
- [阶段 2 探索分析报告](reports/02_eda_report.md)：主要发现、首批图组、聚类敏感性与图形样式。
- [阶段 2 独立审核](reports/02_independent_review.md)：独立重算、图形检查与供用户参考的审核建议。
- [执行计划](PLAN.md)：全部阶段、统计设计、GPU 拓展及归档规则。
- [工作日志](WORK_LOG.md)：关键提示词、授权、执行与排错记录。
- [环境说明](environment/README.md)：本地运行入口、固定依赖与远端资源。

## 第7–8阶段成果

GEARS、scGPT blood和scFoundation＋GEARS各完成三组划分，并与四种基线统一比较。共同基因上Ridge加和的平均MSE为0.001970，与训练均值接近；训练高变基因上GEARS的平均MSE最低，为0.01974。按最终训练覆盖的分层及预先确定的PRDM1案例见增强报告。

每目标20次重新调参后，展示40个TF候选的最低选择频率为90%，40项系数中位数方向全部一致。88项正式实验自检、补充收敛检查及最终远端/本地SHA-256核验通过。九个最佳检查点、原始和处理数据、预测、代码、五组新增图及工作日志均在本地。

全部远端计算已经结束，GPU实例可由用户关闭。公共预训练及恢复/早期试算权重在远端保留的明细见最终归档说明。最终PDF转换和作业提交由用户完成。

## 第 6 阶段成果

GEARS、scGPT blood、scFoundation＋GEARS 均完成 4 次真实参数更新，输出同一组 14 个留出细胞、18,977 个基因的预测。统一输入含 84,143 个细胞、225 个条件和三个官方 simulation 划分。26 项自检通过，126 个远端科学文件（3,847,328,837 字节）已回存并通过 SHA-256 核验。公共预训练权重留在远端，本次试点检查点均在本地。正式效果比较已在第7步完成。

## 第 5 阶段成果

核心报告复用四张已验收主图，以 Lasso 为主线，将完整多模型比较放入附录。16 篇文献围绕红系 MYB、浆细胞 BCL6、髓系 KLF6／PRDM1、探索性候选及基础模型／编辑／递送研究连接组织。正文明确现有实测结果与未来模型实验的范围。

本阶段没有新增模型拟合，原 1,031 个科学结果文件的摘要保持不变。文档链接、引用覆盖、7 个代表文献案例的统计值及四张图均通过验收，记录见 `logs/stage5_validation.json`。报告可通过 `sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage5_report.py` 重建；最终 PDF 转换及提交由用户完成。

[Prism 模型比较样板](prism/stage5/05_Prism_model_comparison.pzfx)已完成标题、配色、原生 PDF／300 dpi TIFF 导出。用户现已取消后续 Prism 制图；已保存样板和六组输入保留为历史附加成果，核心报告继续使用已审核的原图。详见[样板记录](prism/stage5/README.md)及[迁移输入](prism/migration/README.md)。

## 第 4 阶段成果

四目标完成六类模型的共享 5×3 嵌套验证。最低误差 Ridge 的折外 R² 分别为 B 细胞 0.4350、浆细胞 0.4429、单核细胞 0.4958、红细胞 0.6148。100 次 80% 子抽样得到 41、43、52、42 个稳定优先 TF 扰动条目；四目标共同身份为 NFYB 和 OVOL2。自身读出删除的 34 个可匹配条目均保留方向，另 6 项明确记录读出缺失。

7 组图以 PDF、SVG、450 dpi PNG 保存于 `figures/04_*`，完整源表见 `tables/04_*`；40 个展示候选的证据集中在 [候选证据矩阵](tables/04_candidate_evidence_matrix.csv)。landmark、逐基因标准化和三组 965 基因探索对照各自保留调参路径、预测与指标。正式稳定性均从零初始化，早期初始化试算保留在独立探索目录。

17 项 CPU 求解器对照、9 项主体数值验收、补充对照验收及 21 个图形文件检查均通过。824 个远端科学文件按路径、字节数和 SHA-256 与本地核对一致，记录见 `logs/stage4_remote_archive_validation.json`。全本地交付快照另见 `logs/stage4_artifact_manifest.json`。

## 第 3 阶段成果

11,922 个共同基因、2,717 个扰动、21 个目标采用共享 5×3 嵌套验证。主要 20 个目标的最低误差 Lasso 折外 R² 中位数为 0.4504；B、浆、单核、红细胞分别为 0.4271、0.4289、0.4835、0.6093。最低误差与 1-SE 两套预测、系数和参数均已保存。

5 组图位于 `figures/03_*`，提供 PDF、SVG、450 dpi PNG；源表位于 `tables/03_*`，中文图注见 [图注文件](reports/03_figure_captions.md)。每目标的折外预测、调参曲线和系数在 `results/lasso/target_*`；24 个训练上下文的完整 50 点路径与对偶证书在 `results/lasso_gpu_paths/`。初筛候选按系数排序，重采样稳定性属于第 4 阶段。

10 项实质性验收全部通过，另完成 B 细胞两套参数的严格 scikit-learn 交叉核验。结果见 `logs/stage3_validation.json` 和 `logs/03_solver_crosscheck_B_cell.json`。远端同步与本地交付快照分别记录在 `logs/stage3_remote_archive_validation.json` 和 `logs/stage3_artifact_manifest.json`。

## 第 2 阶段成果

7 组图形提供 PDF、SVG、450 dpi PNG：表达分布、HVG100 热图、PCA、HSC 差异、聚类对照、mDC 补充与聚类敏感性。图形位于 `figures/02_*`，源数据位于 `tables/02_*`；中文图注见 [图形规范](reports/02_figure_captions.md)。

PCA 前两轴解释 51.54%、前三轴解释 63.67%；四类重点目标具有不同的参照基因模式。完整解释见报告。验收脚本独立核对输入摘要、HVG、PCA、差值、距离、30 组聚类设置、Prism 数据和图形文件，结果见 `logs/stage2_validation.json`。

Prism 已载入数据，创建原生图形时界面超时；数据表保留于 `prism/`，用户已同意后续补充原生样板。审核原文保留在报告中，处理记录见工作日志。全阶段本地快照及校验值见 `logs/stage2_artifact_manifest.json`。

## 第 1 阶段成果

| 内容 | 文件 |
|---|---|
| 去重血液表达、X/Y、基因/扰动/目标顺序及掩码 | `data/processed/analysis_matrices.h5` |
| 原始文件路径及 SHA-256 | `data/input_manifest.json` |
| 官方注释原文件、来源、时间及摘要 | `data/references/` |
| 全扰动身份与 TF 资格 | `tables/01_perturbation_annotations.csv` |
| 291 个 TF 扰动条目 | `tables/01_TF_perturbations.csv` |
| 完整 1,639 个 TF 目录覆盖情况 | `tables/01_TF_catalog_coverage.csv` |
| 谱系参照的表达读出与扰动覆盖 | `tables/01_lineage_TF_coverage.csv` |
| 共同基因、历史名称映射与 landmark | `tables/01_shared_gene_annotations.csv` |
| 重复基因处理和细胞/扰动 QC | `tables/01_*qc.csv`、`tables/01_cell_quality.csv`、`tables/01_duplicate_gene_audit.csv` |
| 数值摘要与验收结果 | `logs/stage1_summary.json`、`logs/stage1_validation.json` |
| 本阶段本地文件清单及 SHA-256 | `logs/stage1_artifact_manifest.json` |

## 数据接口

HDF5 的 `blood/expression` 为 19,994×22 去重血液矩阵，分别以 `blood/genes`、`blood/cell_types` 标识。

`model/X` 为 11,922×2,717 的 L1000 预测变量；`model/Y` 为 11,922×21 的目标减 HSC 表达差值。两者共用 `model/genes`；列名为 `model/perturbations` 和 `model/target_cells`。所有数组仍为输入尺度，折内标准化将在建模阶段执行。

`model/is_landmark` 对应基因维度，`model/is_TF` 和 `model/primary_TF_eligible` 对应扰动维度，`model/is_main_target` 对应目标维度。`DDIT3` 与 `(DDIT3)` 保留独立特征列，注释中映射到同一身份。

## 第 4 阶段复现

已归档结果可直接重建图表和报告：

```bash
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage4_figures.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/verify_stage4.py --figures-only
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage4_report.py
```

完整计算及数值验收顺序见[环境说明](environment/README.md)。每目标20次重新调参的重采样已在第7步完成。

## 第 3 阶段复现

使用已归档结果重建图表和执行数值验收：

```bash
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/lasso_figures.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/verify_stage3.py
```

从冻结矩阵重新拟合的环境安装与 `prepare → GPU 路径 → 结果汇总 → 图表 → 验收` 顺序见[环境说明](environment/README.md)。现有结果支持同输入、同配置的断点续跑；新实验在独立任务目录运行。

## 第 2 阶段复现

在课程根目录执行：

```bash
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/eda_figures.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/verify_stage2.py
```

第 1 阶段从原始文件重建的命令保留在工作日志中；日常使用已冻结的 HDF5。阶段清单是各次交付时的快照，后续报告或日志更新会产生新摘要，原始数据与冻结矩阵的摘要单独核验。

原始文件与本阶段分析成果均保存在本地。远端 `/root/HW2_20261107` 已安装独立统计和 GPU 求解环境，完成第 3、4 阶段计算；第7–8阶段也已完成，全部远端任务结束且本地归档已验收。
