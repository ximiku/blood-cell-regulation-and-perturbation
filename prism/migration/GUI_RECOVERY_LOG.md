# Prism 迁移 GUI 恢复记录

日期：2026-11-07；GraphPad Prism 10.5.0 (673)。

## 已确认成果

`../stage5/05_Prism_model_comparison.pzfx` 已实际在 Prism 创建并保存。`../../figures/05_Prism_model_comparison.pdf` 与 `.tiff` 已通过 Prism File → Export 原生导出。TIFF 实测 1388×1009 像素，约300 dpi。PDF审核预览位于 `../stage5/05_Prism_model_comparison_preview.png`。图为4目标、2模型、合并折外R²、无误差条、0–0.7范围、主标题Out-of-fold R²。

## 扩大迁移尝试

用户明确授权 AppleScript 仅激活Prism、打开本任务原生脚本及保存导出后，执行官方入口：activate + open POSIX file。首次沙箱内启动报 -10822；提升沙箱权限后成功返回脚本对象，无自动审批拒绝。

- `activate_and_open_script.applescript` 打开 `open_pca.pzc`：Prism前台恢复，窗口标题栏与窗口按钮恢复彩色，但PCA输入未加载。
- `activate_test.applescript` 打开 `activation_test.pzc`：诊断脚本采用ASCII/CRLF，包含CreateLog、写done标记、Open单表PCA；系统返回成功打开脚本对象，但 `activation_test_done.txt` 与 `.log` 未产生。因此不能声称原生脚本已经运行。
- 当前 File菜单的New、Open、Close、Save、Export、Run Script全部disabled。CUA只暴露旧的Untitled主窗口（RMS expression difference vs HSC），未暴露可点击Run/Cancel的脚本模态窗口。
- 按恢复要求尝试Escape一次、Command-period一次；均无变化。cua.getState未列出单独Prism模态表面或浏览器表面。
- 已关闭菜单，保留应用与所有原有文档；未强制退出、未丢弃旧Untitled内容。

## 完成边界

尚未创建新的XY/TF/热图原生模板；`*_data.pzfx` 是输入数据，prepared_native官方Template克隆是待Prism渲染和检查的项目，不能作为已导出成品。当前可审核的Prism成品仅模型比较图PDF/TIFF及其原生项目。

恢复建议：由用户查看Prism是否存在未被自动化接口暴露、被遮挡或在屏幕外的脚本对话框，取消/关闭该模态状态（保留未保存文档）后继续。下一步优先单表PCA创建一个XY散点模板，再复用数据；TF官方森林图模板需先把log横轴改linear并核验真实上下限；热图输入默认按钮367。
