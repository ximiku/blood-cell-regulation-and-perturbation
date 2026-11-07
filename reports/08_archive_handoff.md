# 第8步：本地归档与使用说明

远端与本地的路径、字节数和SHA-256核验全部通过，报告、结果表及处理数据已实际打开。

## 阅读入口

1. [一页结果说明](08_results_brief.md)：主要结论及科研价值。
2. [增强分析报告](08_enhanced_report.md)：核心作业、GPU拓展、主图和文献。
3. [GPU实验完整报告](07_gpu_benchmark_report.md)：模型比较、分层、案例及资源。
4. [核心补充材料](05_core_supplement.md)：完整统计模型和候选表。

## 已保存的成果

本地根目录：`/path/to/user/Desktop/CourseData/26_fall/Statistical_analysis_of_genomics_data/HW2/solution`。

- `data/`：课程输入摘要、冻结矩阵、Norman原始与处理数据、细胞/基因/扰动清单、三组划分和参考对照。
- `results/stage7/`：四种基线、三种GPU模型各三个划分的条件预测；九个正式最佳检查点；20次核心重新调参的全部路径、系数和诊断。
- `tables/07_*`：逐条件、逐划分、模型汇总、覆盖范围及候选重采样结果；`tables/08_*`：核心候选与Norman的衔接案例。
- `figures/07_*`：五组图的PDF、SVG和450dpi PNG及图形清单；[中文图注](07_figure_captions.md)。
- `scripts/`、`config/`、`vendor/stage6/`和`environment/`：执行代码、固定配置、原模型源码及依赖记录。
- [工作日志](../WORK_LOG.md)：用户授权、训练预算调整、排错及结果路径。

正式GPU实验自检共88项通过。20次重新调参也通过分折、参数选择和全部收敛证书检查。新增远端归档清单含232个文件、4,328,108,490字节；第3、4阶段冻结数据/结果及第6阶段清单在最终验收时重新核对。各清单存在共享文件，统计时按清单分别保留记录。

核验入口：[正式实验验收](../logs/stage7_validation.json)、[重新调参验收](../logs/stage7_core_retuning_validation.json)、[新增远端文件清单](../logs/stage8_remote_manifest.json)。最终归档结果保存于 `logs/stage8_archive_validation.json`，完整本地快照及SHA-256保存于 `logs/stage8_local_artifact_manifest.json`。

## 复现与训练预算

从已归档结果重建图表、报告及指标的命令见[环境说明](../environment/README.md)。用户在运行中将训练上限由20轮收紧至8轮；全部最终检查点来自前8轮。调整前已完成的GEARS种子1/2和scFoundation种子1额外轮次保留在CSV历史中。正式配置见[配置文件](../config/stage7_config.json)。

scFoundation正式检查点保存训练后的预测部分，冻结编码器的参考对照表示已在本地保存；处理新的输入可按来源记录重新获取公共权重。临时稠密表达缓存是归档稀疏矩阵的等值展开，可由准备脚本重建。

## 远端仍保留的模型文件

远端任务根目录为 `/root/HW2_20261107/solution`。公共预训练权重、恢复检查点和早期试算检查点保留在远端；九个正式最佳检查点均已在本地。完整路径、字节数、SHA-256及伴随配置文件见[远端模型清单](../logs/stage8_remote_model_inventory.json)，权重来源和版本见[预训练审计记录](06_pretraining_audit.md)。

| 远端保留文件（相对任务根目录） | 字节数 |
|---|---:|
| `results/stage7/GEARS/seed1/latest.pt` | 75,384,656 |
| `results/stage7/GEARS/seed2/latest.pt` | 75,384,400 |
| `results/stage7/GEARS/seed3/latest.pt` | 75,384,144 |
| `results/stage7/scFoundation_GEARS/seed1/latest.pt` | 630,163,930 |
| `results/stage7/scFoundation_GEARS/seed2/latest.pt` | 630,163,802 |
| `results/stage7/scFoundation_GEARS/seed3/latest.pt` | 630,163,546 |
| `results/stage7/scGPT_blood/seed1/latest.pt` | 470,074,893 |
| `results/stage7/scGPT_blood/seed2/latest.pt` | 470,074,893 |
| `results/stage7/scGPT_blood/seed3/latest.pt` | 470,074,893 |
| `results/stage7_exploratory/constant_lr/GEARS/seed1/best.pt` | 25,131,161 |
| `results/stage7_exploratory/constant_lr/GEARS/seed1/latest.pt` | 75,383,312 |
| `results/stage7_exploratory/constant_lr/GEARS/seed2/best.pt` | 25,131,161 |
| `results/stage7_exploratory/constant_lr/GEARS/seed2/latest.pt` | 75,383,312 |
| `results/stage7_exploratory/constant_lr/scFoundation_GEARS/seed1/best.pt` | 290,420,192 |
| `results/stage7_exploratory/constant_lr/scFoundation_GEARS/seed1/latest.pt` | 630,162,586 |
| `results/stage7_exploratory/constant_lr/scGPT_blood/seed1/best.pt` | 158,078,717 |
| `results/stage7_exploratory/constant_lr/scGPT_blood/seed1/latest.pt` | 470,073,677 |
| `models/stage6/scFoundation/models.ckpt` | 1,432,587,886 |
| `models/stage6/scGPT_blood/best_model.pt` | 155,981,354 |


最终PDF转换、AI智能体配置说明与作业提交由用户完成。实例由用户关闭；关闭提示在全部计算及最终归档核验结束后给出。
