# 运行环境与复现

## 本地分析与绘图

任务入口 `scripts/run_python.sh` 使用 `/opt/anaconda3/envs/base_clone1/bin/python`（Python 3.12.7），设置独立缓存、`PYTHONNOUSERSITE=1` 与任务 `.deps` 路径。主要版本为 NumPy 1.26.4、SciPy 1.13.1、pandas 2.2.2、scikit-learn 1.5.1、h5py 3.11.0、Matplotlib 3.9.2、seaborn 0.13.2。阶段 1、2 直接依赖见各自 requirements 文件；完整本地环境清单见 `package_inventory.txt`。

冻结输入为 `data/processed/analysis_matrices.h5`。日常建模直接使用该文件；`prepare_data.py` 属于阶段 1 原始输入重建，会重新生成交付状态与文件摘要。原始参考下载及 URL、日期、SHA-256 保存在 `data/references/`。

本地 Celer 0.7.4 安装在任务 `.deps` 中，用于曾经的 CPU 路径试算。正式阶段 3 结果使用下述 GPU 求解器。本机 Prism 10.5.0（build 100673）的原生样板经用户同意后续补充。

## 远端实际环境

用户授权实例为 `root@gpu.example.org:22`，使用本机 `/path/to/ssh_private_key`；私钥保持在本机。SSH 使用 ControlMaster 复用连接，避开网关并发握手限制。

硬件为 4×NVIDIA A100-SXM4-80GB，驱动 580.95.05，CUDA 13.0；容器 CPU 配额 52 核、内存配额 480 GiB。`remote_inventory.json` 保存首次资源清点。正式启动时确认 `/app` 是只读挂载，实际工作目录改为 **`/root/HW2_20261107`**，使用约 30 GiB 的容器系统盘。后续大数据阶段需再次检查可写空间。

| 环境 | 解释器 | 主要依赖与作用 |
|---|---|---|
| 统计 | `/root/HW2_20261107/envs/statistics/bin/python` | Python 3.10.12；NumPy 1.26.4、SciPy 1.13.1、scikit-learn 1.5.1、pandas 2.2.2、h5py 3.11.0；生成划分、训练折统计量、预测与评价 |
| GPU 路径 | `/root/HW2_20261107/envs/gpu_lasso/bin/python` | Python 3.10.12；CuPy 14.1.1（CUDA 13）、NumPy 2.2.6、cuda-pathfinder 1.8.3；float64 批量 FISTA |

`requirements-stage3-remote.txt` 与 `requirements-stage3-gpu.txt` 固定正式计算依赖。`remote_wheels/` 与 `gpu_solver_wheels/` 保存 Linux x86_64 / CPython 3.10 安装包，`stage3_wheel_manifest.json` 保存校验值。Celer 另用于早期 CPU 试算，正式 GPU 路线无需该包。阶段 3 未使用预训练权重。

系统 Python 缺少 ensurepip，实际采用 `python3 -m venv --without-pip` 创建上述环境，再将已保存的 pip wheel 加入临时 PYTHONPATH 引导安装。可在相同 Linux/Python/CUDA 条件下按如下方式重建；`HW2_TASK` 指向新的独立任务根目录，目录布局为 `solution/` 与 `envs/`：

```bash
export HW2_TASK=/root/HW2_20261107_reproduce
python3 -m venv --without-pip "$HW2_TASK/envs/statistics"
python3 -m venv --without-pip "$HW2_TASK/envs/gpu_lasso"
for HW2_PIP_WHEEL in "$HW2_TASK"/solution/environment/remote_wheels/pip-*.whl; do
  PYTHONPATH="$HW2_PIP_WHEEL" "$HW2_TASK/envs/statistics/bin/python" -m pip install --no-index "$HW2_PIP_WHEEL"
  PYTHONPATH="$HW2_PIP_WHEEL" "$HW2_TASK/envs/gpu_lasso/bin/python" -m pip install --no-index "$HW2_PIP_WHEEL"
done
"$HW2_TASK/envs/statistics/bin/python" -m pip install --no-index --find-links "$HW2_TASK/solution/environment/remote_wheels" -r "$HW2_TASK/solution/environment/requirements-stage3-remote.txt"
"$HW2_TASK/envs/gpu_lasso/bin/python" -m pip install --no-index --find-links "$HW2_TASK/solution/environment/gpu_solver_wheels" -r "$HW2_TASK/solution/environment/requirements-stage3-gpu.txt"
```

先将 `solution/config`、`solution/scripts`、冻结 HDF5、阶段 1 注释表和上述依赖包复制到新目录。全新重算使用空的 `results/lasso`、`results/lasso_gpu_paths`、`.cache/lasso` 目录；保留原交付目录便于比较。

## 第 3 阶段运行顺序

在准备好上述目录后运行：

```bash
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/lasso_analysis.py" --prepare --threads 4
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/run_gpu_lasso.py"
```

第一条命令创建固定的 5 外层 / 3 内层划分及 24 个训练上下文的均值、尺度、Gram 矩阵。第二条命令按 4 张 GPU 分配完整路径，检查每个上下文的 NPZ 与 JSON 成对完成，再自动执行单扰动相关和 21 个目标的预测、调参选择及指标汇总。正式输入、参数与种子见 `config/analysis_config.json`。

断点续跑按已有完成文件跳过路径与目标，适用于**同一冻结输入和配置**。变更输入、划分或算法参数时使用新任务目录，防止复用旧缓存。派生 `.cache/lasso` 可按第一条命令重建；完整 GPU 路径、系数和预测已归档本地。

将结果增量同步回本地后，在课程目录执行：

```bash
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/lasso_figures.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/verify_stage3.py
```

本地已有归档结果，可直接执行这两条命令重建图表和复核数值。`check_lasso_solver.py` 另提供 B 细胞最终参数的严格 scikit-learn 交叉核验，已完成的记录为 `logs/03_solver_crosscheck_B_cell.json`。这项核验耗时较长，常规图表重建无需重复。

阶段 3 所需远端计算已结束，全部科学结果已回存；实例继续保留给后续多模型与外部数据实验。项目全部远端工作结束、最终归档核验通过后，再提示用户关闭实例。

## 第 4 阶段运行顺序

第 4 阶段沿用上述统计与 CuPy 环境，无新增依赖。它以冻结 HDF5、第 3 阶段划分、四目标的完整模型文件和阶段 1/3 注释表为输入。正式设置另存 `config/stage4_config.json`；重复实验使用独立任务目录。

在该目录布局下，用统计环境运行准备脚本，再启动 GPU 路径与 100 次子抽样：

```bash
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/stage4_prepare.py"
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/run_stage4_gpu.py"
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/check_stage4_solvers.py"
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/stage4_assemble.py"
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/run_stage4_gpu.py" --jobs data/splits/stage4_self_readout_jobs.json
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/stage4_interpret.py"
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/verify_stage4.py"
```

自身读出任务由稳定性前 10 项生成，规范身份相同的删除任务合并。所有新任务保存在 `results/stage4/`，未改写第 3 阶段科学结果。`--models-only` 与 `--stability-only` 可分别汇总已完成的对应任务；正式交付包含两部分。数值验收直接读取原始训练行，所需 `.cache/stage4` 可由准备脚本重建。

当前代码直接从零系数执行正式稳定性拟合。早期全数据系数初始化的试算及其对照属于既有探索归档，新目录复现无需执行该试算。三组等规模非 landmark 补充对照依次运行：

```bash
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/stage4_size_control.py" --prepare
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/run_stage4_gpu.py" --jobs data/splits/stage4_size_control_jobs.json
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/stage4_size_control.py"
"$HW2_TASK/envs/statistics/bin/python" "$HW2_TASK/solution/scripts/verify_stage4_size_control.py"
```

将结果、源表、划分及日志同步回本地后，生成图表、检验图形文件并撰写报告：

```bash
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage4_figures.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/verify_stage4.py --figures-only
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage4_report.py
```

独立 subagent 审核在数值结果、图组与报告完成后进行，审核意见单独保存。最终入口状态、报告中的审核结论和工作日志随审核结果更新；生成阶段文件清单后提交用户审核。

## 第 6 阶段 GPU 模型环境

已完成三个模型最小试点。远端使用 `/root/HW2_20261107/envs/perturbation_sys/bin/python`；依赖与解释器见 [stage6_environment.json](stage6_environment.json) 和 [requirements-stage6.txt](requirements-stage6.txt)。数据与输入适配、预训练来源及命令见 [第6步来源说明](../reports/06_pretraining_audit.md)。正式第7步尚未开始。公共权重保留远端，本次试点检查点已回存。


## 第7–8阶段复现与归档

正式GPU实验在 `/root/HW2_20261107/solution` 内执行。公共模型和输入来源见 `reports/06_pretraining_audit.md`，固定源码位于 `vendor/stage6`。GEARS/scGPT/scFoundation运行环境为 `envs/perturbation_sys/bin/python`；核心Lasso重新调参使用 `envs/gpu_lasso/bin/python`。三个种子的条件清单在 `data/splits/stage6/seed{seed}.json`；统一设置见 `config/stage7_config.json`。

用户在运行中批准将训练上限由20轮收紧至8轮；全部最终选用检查点位于前8轮。调整前已完成的额外轮次保留在CSV训练历史中。`stage7_budget_queue.py` 记录本次从运行队列交接到新预算的过程；全新复现直接使用当前8轮设置的 `stage7_train.py` 即可。

新任务依次执行 `stage7_prepare.py`、`stage7_baselines.py`，再对三模型、三种子调用 `stage7_train.py --model MODEL --seed SEED --device cuda:GPU`。原始预处理和图构建说明沿用第6阶段记录。冻结编码器表示缓存与64个参考对照一一对应，缓存文件已归档；临时稠密表达矩阵为已保存稀疏输入的等值展开，可由准备脚本重新生成。

本地结果齐备后，在课程根目录运行：

```bash
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage7_retuning_summary.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage7_evaluate.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/verify_stage7.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage7_figures.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage8_report.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/verify_stage8_archive.py
```

图形输出为PDF、SVG和450dpi PNG，中文图注位于 `reports/07_figure_captions.md`。正式最佳检查点、预测和训练记录在 `results/stage7`；临时恢复检查点与公共预训练权重由最终远端模型清单单独列明。
