# Stage 5 Prism 原生样板（已恢复并导出）

2026-11-07；GraphPad Prism 10.5.0 (673)。

## 成果

- `05_Prism_model_comparison.pzfx`：Prism 本机实际创建、保存的原生图形项目（含 Template 图形数据）。
- `../../figures/05_Prism_model_comparison.pdf`：Prism 原生矢量 PDF 导出。
- `../../figures/05_Prism_model_comparison.tiff`：Prism 原生 300 dpi TIFF（实测 1388 × 1009 像素） 导出（RGB、白底、压缩）。
- `05_Prism_model_comparison_preview.png`：由上述原生 PDF 光栅化的审核预览，并非重新绘图。
- `model_comparison_source.csv`：8 行原始指标；`model_comparison.tsv`：宽表输入；`model_comparison_data.pzfx`：脚本生成的纯数据载入文件。

## 数据和图形

源表为 `../../tables/04_model_metrics.csv`，选择 `dataset=raw`、`selection=min`、模型 Lasso/Ridge，目标 B cell、plasma cell、monocyte、erythrocyte。Prism 中通过 Create New Graph → Grouped → Individual values → Interleaved scatter → No line or error bar 创建。

四组各两个点，纵轴 0–0.7、间隔 0.1；无连接线、无误差条。英文类别与图例。主标题 Out-of-fold R²；Prism Colorblind Safe 配色：Lasso 黑色方块，Ridge 洋红三角形。两种不同符号同时提供非颜色区分。为避开自绘文字编辑控件的 AX 错误，指标通过主标题说明，隐藏默认的 YTitle；纵轴保留数字刻度。

每个点为一个目标、一个模型的合并折外 R²（11,922 个基因的折外预测汇总）；四类细胞不是重复样本，不添加 SD/SEM、置信区间或显著性标记。

| Target | Lasso | Ridge |
|---|---:|---:|
| B cell | 0.4270963427 | 0.4350207887 |
| Plasma cell | 0.4289249270 | 0.4428658098 |
| Monocyte | 0.4834976082 | 0.4958167708 |
| Erythrocyte | 0.6092764273 | 0.6148057794 |

## 操作记录

最初尝试中发生 `noWindowsAvailable`、自绘文字点击 `AXError.notImplemented` 及多次 `timeoutReached`，初次交付只有保存的草图。用户随后扩大 Prism 迁移范围，重新绑定应用后响应恢复，前次待完成的标题/隐藏设置已生效。已正常保存并通过 File → Export 导出 PDF、TIFF。已人工查看 PDF 预览，8 点、标签、范围完整，无默认占位标题。项目仍为兼容模式 `.pzfx`，可在 Prism 中继续编辑。
