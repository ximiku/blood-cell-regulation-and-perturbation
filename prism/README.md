# Prism 样板状态

2026-11-07，GraphPad Prism 10.5。

- 数据：四类重点目标相对 HSC 的跨基因 RMS，每个类型一个摘要值；完整分母为 19,994 个去重基因。
- `HSC_difference_data.tsv` 为通用表格输入；`HSC_difference_data.pzfx` 为脚本生成的 Prism 数据 XML。
- 本机界面已成功载入四列数据，列名及数值与源表一致。
- 打开 XML 时出现格式不完全兼容提示，继续打开后数据可见。创建 Column 图时 UI 工具超时，之后 Prism 界面读取继续超时。
- 本轮未获得可核验的 Prism 原生图形或 PDF/TIFF 导出。现有 `.pzfx` 只包含数据，尚无图形样式。
- 官方随附 XSD 在 lxml 构建 schema 时报告非确定性 content model，未取得 schema 验收；该错误与界面超时分别记录。

已完成的专业图形位于 `../figures/`，包括 7 组 PDF/SVG/450 dpi PNG。原生 Prism 样板待界面可用后补充。可直接将 TSV 四列导入 Column 表，使用各列单点展示、隐藏误差条，纵轴标记 `RMS difference vs HSC (input units)`，并沿用主图配色。

本文件描述本轮实际操作结果；无 Prism 图形成功导出的记录。
