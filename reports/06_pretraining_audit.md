# 第 6 步：模型、权重与预训练数据来源

## 固定源码

| 项目 | 本次实际版本 |
|---|---|
| GEARS | [`f374e43e197b`](https://github.com/snap-stanford/GEARS/tree/f374e43e197b295016d80395d7a54ddb81cc6769) |
| scGPT | [`cebd6fae655b`](https://github.com/bowang-lab/scGPT/tree/cebd6fae655b9c585a4807daa3ac31bb764f06b4) |
| scFoundation | [`397631c495ed`](https://github.com/biomap-research/scFoundation/tree/397631c495eddf9ad6644fc00c6ea8139e651245) |

本地 `vendor/stage6/` 保存执行所需源码及许可证，文件摘要纳入远端归档清单。GEARS、scGPT、scFoundation 均使用原项目模型类。运行环境为 Python 3.12.13、PyTorch 2.12.0+cu130；依赖完整版本见[环境记录](../environment/stage6_environment.json)。

## 通用预训练权重

**scGPT blood：** 从[官方仓库](https://github.com/bowang-lab/scGPT)所链接的 [blood 权重目录](https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU)获取 `best_model.pt`、`args.json`、`vocab.json`。官方描述为 10.3M 血液与骨髓细胞预训练。下载使用通用 blood 权重，扰动模型中新建任务相关模块；153 个编码器张量、37,927,424 个参数逐项核对加载值一致。

**scFoundation：** 模型架构及加载逻辑来自[官方仓库](https://github.com/biomap-research/scFoundation)，官方描述为 50M 人类细胞预训练。本次权重从 [genbio-ai 公共镜像固定版本](https://huggingface.co/genbio-ai/scFoundation/tree/cb434153a1acfacd215eefc956ea445f7cc39cc3)下载，SHA-256 与该版本 LFS 标识一致。该一致性核验对应镜像发布物；尚未与原发布方文件逐字节比对。[镜像元数据](../logs/stage6/scfoundation_mirror_metadata.json)与[下载记录](../logs/stage6/scfoundation_weight_manifest.json)已保存。模型许可证随源码与远端权重保存。

已读取的发布材料未提供可逐细胞核对的预训练成员清单，因此两套预训练权重与 Norman 数据的样本重叠状态记录为“未核实”。报告将准确区分通用预训练来源与本任务划分，保留这一适用范围说明。

| 远端路径 | 字节数 | SHA-256 |
|---|---:|---|
| models/stage6/scGPT_blood/args.json | 902 | `8db3c600fec669fcba2619dad82810267d1f1490f64b1d2375ca0c814586a936` |
| models/stage6/scGPT_blood/best_model.pt | 155981354 | `56e68896b1b4e8c708b695f1a60bb7560e08952ea76842f695109737beaeeb41` |
| models/stage6/scGPT_blood/vocab.json | 761427 | `9b8df99caf1da1d58f0b7c7f871c64b75f2c5fba256d1f18fe14b29c4529e64c` |
| models/stage6/scFoundation/MODEL_LICENSE | 2310 | `bdcd01b3213efb813f4461ef60c2a30d56cf07b63446d234283d8c0547ee3cbb` |
| models/stage6/scFoundation/models.ckpt | 1432587886 | `9f40bf324d3d0084c4b288d06f5af4fddd12206e2a3f022551d12e89e33a0ea9` |

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
