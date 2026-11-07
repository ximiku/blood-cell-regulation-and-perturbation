"""Build the stage-6 handoff from saved pilot results; launches no training."""
from pathlib import Path
import json
from datetime import datetime, timezone
import pandas as pd

R = Path(__file__).resolve().parents[1]
L = R/'logs/stage6'
models = ['GEARS','scGPT_blood','scFoundation_GEARS']
results = [json.loads((R/'results/stage6'/m/'pilot_summary.json').read_text()) for m in models]
audit = json.loads((L/'local_validation.json').read_text())
assert audit['passed']
data = json.loads((L/'remote/data_preparation.json').read_text())
splits = pd.read_csv(R/'data/splits/stage6/condition_splits.csv')
resources = pd.DataFrame([{k:r[k] for k in ['model','status','optimizer_steps','training_cells','prediction_cells',
    'shared_output_genes','training_seconds','prediction_seconds','total_seconds','peak_allocated_GiB',
    'peak_reserved_GiB','parameters_total','parameters_trainable']} for r in results])
resources.to_csv(R/'tables/06_pilot_resources.csv',index=False)
split_summary = splits.groupby(['seed','split']).agg(conditions=('condition','size'),cells=('n_cells','sum')).reset_index()
split_summary.to_csv(R/'tables/06_split_summary.csv',index=False)
seen = splits[splits.split=='test'].groupby(['seed','evaluation_group']).size().unstack(fill_value=0)
seen.to_csv(R/'tables/06_test_seen_groups.csv')
def save(path, obj):
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')

config = {
 'status':'executed_and_self_checked','random_seed':20261107,'split_seed':1,
 'source_dataset':'Norman 2019 GSE133344, scFoundation Figshare processed release with total_count',
 'retained_cells':84143,'conditions_including_control':225,'input_genes':19264,'common_output_genes':18977,
 'expression_scale':'publisher log1p normalized, retained without re-normalization',
 'pilot_pairs':'data/splits/stage6/pilot_pairs.csv','batch_size':2,'optimizer_steps':4,
 'training_cells':8,'heldout_prediction_cells':14,'gradient_clip_norm':1.0,
 'models':{
  'GEARS':{'hidden_size':64,'decoder_hidden_size':16,'learning_rate':.001,'loss':'upstream power4 + direction', 'precision':'float32'},
  'scGPT_blood':{'learning_rate':.0001,'loss':'MSE','train_sequence_length':1536,'prediction_genes':18977,
     'always_include_perturbed_genes':True,'precision':'bfloat16 autocast','pretrained_encoder_verified':True},
  'scFoundation_GEARS':{'hidden_size':512,'decoder_hidden_size':16,'learning_rate':.001,
     'loss':'upstream power4 + direction','mode':'v1','highres':0,'precision':'float32','encoder_frozen':True,
     'encoder_eval_mode':True,'append_original_total_count':True}},
 'graphs':{'GO':'100 declared targets; Jaccard >0.1; top21 incoming including self',
    'coexpression':'seed1 training controls and single perturbations only; absolute Pearson >0.4; top21 incoming'},
 'review_mode':'self_checks_only_user_requested','Prism':'cancelled_by_user'}
save(R/'config/stage6_config.json', config)
formal = {
 'status':'proposed_waiting_for_user_stage7_approval','execution_started':False,
 'task':'held-out perturbation-condition expression-change prediction',
 'seeds':[1,2,3],'split_manifests':'data/splits/stage6/seed{seed}.json',
 'models':['zero_change','training_condition_mean','Ridge_additive','direct_single_additive',*models],
 'max_epochs':20,'early_stopping_patience':5,
 'early_stopping_metric':'validation condition-equal mean MSE on expression changes across all shared output genes',
 'best_checkpoint':'minimum validation metric; fresh training, pilot weights not reused',
 'training_conditions':'all conditions in each training manifest; train controls only as reference inputs',
 'gene_sets':['18977 shared outputs','1000 high-variance shared genes selected using each training split'],
 'evaluation_groups':['seen0','seen1','seen2','unseen_single'],
 'seen_definition':'gene occurs in final training conditions; original GEARS labels preserved separately',
 'condition_aggregation':'equal weight per condition; true mean across all cells in that condition',
 'control_reference_for_evaluation':{'n_cells':64,'seed':20261107,'selection':'fixed training-control sample shared across all models and seeds',
    'prediction':'average predicted response across these controls','control_baseline':'mean over all training controls'},
 'metrics':['expression-change MSE','expression-change RMSE','Pearson on changes','Spearman on changes','relative MSE improvement over each baseline'],
 'supplement':'repeat metrics excluding directly perturbed genes; report undefined correlations explicitly',
 'direct_additive_eligibility':'all corresponding single-perturbation conditions present in training; report coverage',
 'Ridge_alpha_grid':[.0001,.001,.01,.1,1,10,100,1000,10000],
 'Ridge_selection':'validation condition-equal MSE',
 'initial_batch_sizes':{'GEARS':32,'scGPT_train':32,'scGPT_full_gene_eval':4,'scFoundation_GEARS':6},
 'learning_rates':{'GEARS':.001,'scGPT_blood':.0001,'scFoundation_GEARS':.001},
 'scGPT_training_sequence_length':1536,'scGPT_always_include_perturbed_genes':True,
 'scFoundation_frozen_encoder':True,
 'resource_adjustment':'batch sizes may be reduced for memory; log final values; no extra pilot rounds',
 'speed_strategy':'parallelize model/seed jobs across 4 GPUs; cached data; early stopping; vectorized condition summaries',
 'coexpression_graphs':'seed-specific training-only graphs; seed2/3 built at stage7 start',
 'runtime_estimate':'hour-scale; refine after first full epoch, minimal-trial timing is insufficient for a reliable total',
 'self_checks_only':True}
save(R/'config/stage7_proposed_config.json',formal)

resource_rows='\n'.join(f"| {r['model']} | 通过 | {r['total_seconds']:.1f} | {r['peak_allocated_GiB']:.2f} / {r['peak_reserved_GiB']:.2f} |" for r in results)
split_rows='\n'.join(f"| {seed} | {int((splits[(splits.seed==seed)&(splits.split=='train')]).shape[0])} | {int((splits[(splits.seed==seed)&(splits.split=='val')]).shape[0])} | {int((splits[(splits.seed==seed)&(splits.split=='test')]).shape[0])} | {seen.loc[seed,'seen0']} / {seen.loc[seed,'seen1']} / {seen.loc[seed,'seen2']} |" for seed in [1,2,3])
report=f'''# 第 6 步：外部扰动数据与 GPU 试点

**三个模型已全部跑通，可以进入正式实验。** 本步完成统一数据、三组划分、真实梯度更新和留出预测；当前成果回答的是模型流程能否运行，预测效果将在第 7 步用完整训练比较。

## 这一步与作业主线的联系

核心报告通过 Lasso 从表达差异中筛选候选调控因子。拓展实验进一步研究一个清楚的问题：**已知细胞受到哪些基因扰动，模型能否预测它随后发生的表达变化？** GEARS 提供专门的扰动预测模型，scGPT 和 scFoundation 用预训练表示加入比较。这使分析从血细胞状态差异，延伸到状态—干预—响应预测。

## 已完成的结果

三模型共享 8 个训练细胞、4 次参数更新，以及 14 个留出细胞的预测任务；每次更新的 batch size 为 2。预测统一为 **18,977 个基因**，基因顺序和真实表达逐项一致。梯度、参数更新和输出数值检查通过，scFoundation 的 119,251,772 个冻结参数保持无梯度。

| 模型 | 实际训练与预测 | 单次流程耗时（秒） | GPU 峰值已分配 / 已预留（GiB） |
|---|---|---:|---:|
{resource_rows}

表中耗时包括试点的数据载入、模型初始化、4 次更新、预测与检查点保存，设备为 A100 80 GB；首次环境、数据及权重下载耗时单独记录在工作日志。此表用于安排资源，正式性能比较采用完整训练后的测试集结果。完整数值见[资源表](../tables/06_pilot_resources.csv)。

## 统一数据和比较规则

使用 Norman 2019 K562 CRISPRa 的公开处理数据，来自 [scFoundation 数据发布页](https://doi.org/10.6084/m9.figshare.24049200)。原文件含 84,986 个细胞、19,264 个基因及原始总计数 `total_count`。表达矩阵沿用发布者的 log1p 归一化尺度，`total_count` 独立保留供 scFoundation 输入。

统一排除缺少 GO 注释的 IER5L 相关 2 个条件、843 个细胞后，得到 **84,143 个细胞、225 个条件（含对照）、100 个扰动基因**；其中对照 7,275 个。scFoundation 的 19,264 个基因顺序与发布矩阵完全匹配，scGPT 词表覆盖 18,977 个共同输出基因和全部 100 个扰动基因。所有方法共享细胞与条件，输入通过各自模型适配器组织。

| 划分种子 | 训练条件 | 验证条件 | 测试条件 | 测试双扰动 seen0 / seen1 / seen2 |
|---|---:|---:|---:|---:|
{split_rows}

各划分另含 25 个测试单扰动条件。seen0、seen1、seen2 分别表示组合的两个扰动基因中，有 0、1、2 个出现在最终训练条件中。同一条件的所有细胞保持在同一个集合。划分使用固定版本 GEARS 官方 `simulation` 实现；训练基因比例和 seen2 组合训练比例均为 0.75。

划出验证集后，三个种子分别有 12、11、7 个测试组合的实际 seen 分类与官方初始标签不同。表中按最终训练集计算，原始标签同时保留在[条件划分表](../data/splits/stage6/condition_splits.csv)。每个种子的前 1,000 个高变基因均由该训练集单独确定。

GEARS 与 scFoundation 使用相同图结构：GO 图覆盖这 100 个声明的扰动基因；共表达图来自训练对照与训练单扰动细胞。种子 1 图含 17,803 条共表达边和 1,297 条 GO 边。当前模型属于该明示 GO 范围下的实现；公开文件附带的全数据差异基因注释未参与建模。

## 第 7 步的具体运行设置（待审核）

研究主线是“基础模型能否改善未见扰动条件的响应预测”。比较保留零变化、训练条件均值、Ridge 加和、可用条件下的单扰动直接加和，以及三个 GPU 模型。每模型使用三个相同划分，最多 20 epochs，验证指标连续 5 次无改善时早停，选验证误差最小的检查点。

评价以条件等权汇总，比较表达变化的误差、相关性和相对基线改善；分别报告全部共同基因、训练集前 1,000 个高变基因，以及去除直接扰动基因的结果。预测时固定使用同一组 64 个训练对照细胞作为输入并求平均，实际响应使用该条件的全部细胞均值；基线对照均值使用全部训练对照。64 个参考细胞使各模型共享推断输入，并控制推断开销。训练保留各自训练清单中的全部细胞。

四张 A100 按模型和种子并行安排。复用已完成的输入与环境，直接开始正式训练，在首个完整 epoch 后更新预计总耗时。初始 batch size 为 GEARS 32、scGPT 训练 32／全基因预测 4、scFoundation 6；显存需要时调整并记录。scFoundation 冻结编码器，scGPT 使用 blood 通用预训练编码器，正式训练从通用权重或随机初始化开始。

完整可审阅设置见[第 7 步配置](../config/stage7_proposed_config.json)，当前尚未运行第 7 步。权重来源、适配细节与预训练数据记录见[来源说明](06_pretraining_audit.md)。

## 本地交付与自检

{audit['check_count']} 项自检全部通过：输入矩阵可读且数值有限、条件隔离、训练集高变基因复算、实际 seen 分类复算、源 H5AD 试点行逐值匹配、预测基因与细胞顺序、真实表达、梯度和参数更新，以及检查点文件可读性。按路径、字节数和 SHA-256 核对 **126 个远端科学文件，共 3,847,328,837 字节**，本地与远端完全一致。

- [自检记录](../logs/stage6/local_validation.json)与[远端文件清单](../logs/stage6/remote/remote_scientific_manifest.json)。
- 数据：`data/external/stage6/`、`data/processed/stage6/`、`data/splits/stage6/`。
- 三模型预测、训练记录与试点检查点：`results/stage6/`；四次更新的检查点标记为试点用途。
- [实际试点配置](../config/stage6_config.json)、[环境](../environment/stage6_environment.json)、[依赖列表](../environment/requirements-stage6.txt)。
- 公共预训练权重留在远端 `models/stage6/`，本地保留来源、版本与校验值；本次新训练的试点检查点均已回存。

本步按用户要求采用必要自检，取消后续 Prism 制图与独立 subagent 审核。实例当前没有本任务的训练进程，保留用于待审核的第 7 步。
'''
(R/'reports/06_gpu_pilot_report.md').write_text(report)

versions=json.loads((L/'remote/source_versions.json').read_text())
version_lines='\n'.join(f"| {v['name']} | [`{v['commit'][:12]}`]({v['repository']}/tree/{v['commit']}) |" for v in versions)
weights=json.loads((L/'remote/remote_public_weights_verified.json').read_text())['files']
weight_lines='\n'.join(f"| {w['path']} | {w['bytes']} | `{w['sha256']}` |" for w in weights)
provenance=f'''# 第 6 步：模型、权重与预训练数据来源

## 固定源码

| 项目 | 本次实际版本 |
|---|---|
{version_lines}

本地 `vendor/stage6/` 保存执行所需源码及许可证，文件摘要纳入远端归档清单。GEARS、scGPT、scFoundation 均使用原项目模型类。运行环境为 Python 3.12.13、PyTorch 2.12.0+cu130；依赖完整版本见[环境记录](../environment/stage6_environment.json)。

## 通用预训练权重

**scGPT blood：** 从[官方仓库](https://github.com/bowang-lab/scGPT)所链接的 [blood 权重目录](https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU)获取 `best_model.pt`、`args.json`、`vocab.json`。官方描述为 10.3M 血液与骨髓细胞预训练。下载使用通用 blood 权重，扰动模型中新建任务相关模块；153 个编码器张量、37,927,424 个参数逐项核对加载值一致。

**scFoundation：** 模型架构及加载逻辑来自[官方仓库](https://github.com/biomap-research/scFoundation)，官方描述为 50M 人类细胞预训练。本次权重从 [genbio-ai 公共镜像固定版本](https://huggingface.co/genbio-ai/scFoundation/tree/cb434153a1acfacd215eefc956ea445f7cc39cc3)下载，SHA-256 与该版本 LFS 标识一致。该一致性核验对应镜像发布物；尚未与原发布方文件逐字节比对。[镜像元数据](../logs/stage6/scfoundation_mirror_metadata.json)与[下载记录](../logs/stage6/scfoundation_weight_manifest.json)已保存。模型许可证随源码与远端权重保存。

已读取的发布材料未提供可逐细胞核对的预训练成员清单，因此两套预训练权重与 Norman 数据的样本重叠状态记录为“未核实”。报告将准确区分通用预训练来源与本任务划分，保留这一适用范围说明。

| 远端路径 | 字节数 | SHA-256 |
|---|---:|---|
{weight_lines}

## 数据来源

Norman GSE133344 使用 scFoundation 发布的含总计数版本，来源 [Figshare](https://doi.org/10.6084/m9.figshare.24049200)，[下载地址](https://ndownloader.figshare.com/files/44477939)。压缩文件为 614,079,631 字节，MD5 `f9c76065fafa43303dae3f80fe4587a8` 与发布者记录一致。解压 H5AD 为 2,184,186,178 字节，归档清单另记 SHA-256。

表达值为发布者 log1p 归一化矩阵，原始 UMI 总计数在 `obs.total_count` 中；本次保留两者，按模型输入格式适配。GO 映射来自 GEARS 所用 [Harvard Dataverse 发布文件](https://dataverse.harvard.edu/api/access/datafile/6153417)。

## 兼容适配与任务训练

scGPT 通过官方键名转换加载预训练编码器，并采用 PyTorch 原生注意力后端。训练每次抽取 1,536 个基因，包含该批全部直接扰动基因；预测覆盖 18,977 个共同基因。训练使用 bfloat16 自动混合精度，任务损失为表达 MSE。血液预训练输入配置与扰动任务配置分别保留，任务直接使用公开扰动教程的连续表达输入路径。

GEARS 保留官方模型和 power4＋方向损失，隐藏维数 64。scFoundation＋GEARS 保留官方 `v1`、`autobin_resolution_append` 模式及隐藏维数 512；将同一个对照细胞的总计数附在 19,264 个基因之后。编码器处于冻结及 eval 模式，训练 72,583,107 个预测部分参数。两者的损失过滤基因由对应训练条件确定。

模型导入采用最小模块命名空间或抽取原源码定义，避开无关训练器和可选依赖的初始化；scFoundation 公共检查点的读取显式指定 CPU 及 `weights_only=False`。所有适配代码保存在[试点脚本](../scripts/stage6_pilot.py)，完整训练前的科学比较配置保存在[第 7 步设置](../config/stage7_proposed_config.json)。

## 复现入口

本地验收已有成果：

```bash
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/verify_stage6.py
sh HW2/solution/scripts/run_python.sh HW2/solution/scripts/stage6_report.py
```

远端可用环境：`/root/HW2_20261107/envs/perturbation_sys/bin/python`；任务目录：`/root/HW2_20261107/solution`。在已有目录中直接运行试点脚本会覆盖同名试点输出；完整复现实验应复制到独立目录，并通过 `STAGE6_ROOT` 指定。数据准备入口为 `scripts/stage6_prepare.py`，三个模型分别由 `scripts/stage6_pilot.py --model GEARS|scGPT_blood|scFoundation_GEARS --device cuda:0` 调用。正式第 7 步的完整训练驱动在阶段获准后编写。
'''
(R/'reports/06_pretraining_audit.md').write_text(provenance)
print('Wrote stage 6 reports, tables, actual configuration, and proposed stage 7 configuration.')
