# Prism 图形迁移输入

全部数值来自已审核的保存结果。脚本 `scripts/prepare_prism_migration.py` 生成可编辑数据表，尚未生成原生图形。正式图形由 Prism 建立、保存和导出。

1. PCA：两张 XY 图，PC1 对 PC2／PC3。
2. HVG100：100×21 热图，保留 Ward 行列顺序，色阶 ±2.5。
3. OOF：四目标各 11,922 个点，2×2 布局。
4. TF：各 10 个候选，中心为中位数，绝对上下限为 97.5%／2.5% 分位数，线性坐标。
5. 多模型：完整六模型及基线，最低误差／1-SE 两种选择。
6. 跨目标候选：带方向的选择频率热图。

字体 Arial，图内英文，中文图注见报告。详情与配色见 figure_specs.json。

## 当前进度

六组输入均通过独立核验，见 `../../reports/05_Prism_input_review.md`。`prepared_native/` 是替换官方模板 XML 数据表的待渲染输入，压缩 Template 未改动；这些文件需要 Prism 实际打开、调整坐标和导出后才作为成品使用。

官方自动化资料：
- [AppleScript 打开 Prism 脚本](https://www.graphpad.com/guides/prism/latest/user-guide/launching_prism_from_applescript.htm)
- [Prism 原生脚本命令](https://www.graphpad.com/support/faq/list-of-all-script-commands-for-graphpad-prism/)
- [原生图形导出参数](https://www.graphpad.com/guides/prism/latest/user-guide/exporting_a_graph_or_layout.htm)

用户已明确授权 AppleScript 激活 Prism 并打开本任务 `.pzc`。激活成功；脚本未生成预期输出，尚不能宣称批量运行已完成。原生界面的阻碍和恢复尝试记录于工作日志。
