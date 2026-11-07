# Homework 2 可编辑 LaTeX 工程

最终提交文件位于 `../submission/分析报告.pdf` 和 `../submission/工作记录.md`。报告共 35 页：正文 13 页、参考文献 1 页、附录 A–D 21 页。六组正文综合图、完整方法、文献和补充结果都在同一 PDF 中。

## 编辑与编译

- 主文件：`main.tex`；附录：`appendices.tex`；实际排版文献：`references.tex`。
- `references.bib` 保留结构化文献数据，当前采用显式数字文献条目，无需运行 BibTeX。
- 图形位于 `figures/`；正文和附录的 21 个 PDF 图文件是完整编译依赖。
- 采用 TeX Live 2025、XeLaTeX、ctexart 和 Fandol 中文字体。执行 `xelatex main.tex` 两次，或使用 `latexmk -xelatex main.tex`。
- 修改姓名学号：`main.tex` 首页；个人使用记录：提交目录中的工作记录末尾。
- 在线工程：https://latex.pku.edu.cn/project/6aef41a7e2551280bf30ff24 ，项目名 genomic-stat-HW2，已同步正文、附录、文献和全部编译图形。

## 文件与复现

`figure_coverage.csv` 记录原 24 组正式图的覆盖去向。`reproduction/` 保存本轮重组脚本；重新绘图需要原 solution 归档中的结果表和分析环境，脚本原始位置为 solution/scripts。正常编译与文字修改只需要本工程内的 TeX 和 PDF 图，不需要模型、GPU 或原始数据。新综合图同时保留 SVG 便于进一步编辑。

原在线模板已备份在本地 manuscript/template_backup/original_template.zip。旧增强报告和计算结果保持原有归档。本轮重新组织图版、叙事和版面，未重新拟合统计或 GPU 模型。

## 验收

最终 PDF 于 2026-11-08 使用在线 XeLaTeX 成功编译，35 页 A4。已检查全部页面的图表和文字布局，最终修订页再次渲染检查；引用标签完整。submission 目录仅包含分析报告 PDF 与工作记录 MD。完整性记录见 validation.json。
