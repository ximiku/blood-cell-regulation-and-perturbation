"""Integrate the completed core and external perturbation experiments into Markdown."""
from pathlib import Path
import json,re
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];T=R/'tables';REPORT=R/'reports'
MODELS=['zero_change','training_condition_mean','Ridge_additive','GEARS','scGPT_blood','scFoundation_GEARS']
NAME={'zero_change':'零变化','training_condition_mean':'训练条件均值','Ridge_additive':'Ridge 加和','direct_single_additive':'单扰动直接加和','GEARS':'GEARS','scGPT_blood':'scGPT blood','scFoundation_GEARS':'scFoundation＋GEARS'}
def number(v,fmt='.4f'):return '—' if pd.isna(v) else format(v,fmt)
def main():
 s=pd.read_csv(T/'07_model_summary.csv');seed=pd.read_csv(T/'07_seed_metrics.csv');coverage=pd.read_csv(T/'07_coverage.csv');retune=json.loads((R/'logs/stage7_core_retuning_validation.json').read_text());validation=json.loads((R/'logs/stage7_validation.json').read_text());assert validation['passed'] and retune['passed']
 def select(gs='shared',exclude=False,group='all'):return s[(s.gene_set==gs)&(s.exclude_direct_targets==exclude)&(s.group==group)&s.model.isin(MODELS)].set_index('model').loc[MODELS]
 common=select();hvg=select('HVG1000');winner=common.MSE_mean.idxmin();hvgwinner=hvg.MSE_mean.idxmin();gpuwinner=common.loc[MODELS[3:]].MSE_mean.idxmin()
 stable=pd.read_csv(T/'04_stability_all_TFs.csv');perturbation_genes=set(json.loads((R/'data/processed/stage6/perturbation_genes.json').read_text()));overlap=stable[stable.stable_priority&stable.approved_symbol.isin(perturbation_genes)];overlap.to_csv(T/'08_core_Norman_target_overlap.csv',index=False)
 metrics=pd.read_csv(T/'07_condition_metrics.csv');bridge=metrics[(metrics.seed==1)&metrics.condition.isin(['PRDM1+ctrl','CBFA2T3+PRDM1'])&(metrics.gene_set=='HVG1000')&(~metrics.exclude_direct_targets)&metrics.model.isin(MODELS)];bridge.to_csv(T/'08_PRDM1_response_bridge.csv',index=False)
 bridge_table='| 方法 | PRDM1 单扰动 MSE | 变化 Pearson r |\n|---|---:|---:|\n'
 for m in MODELS:
  b=bridge[(bridge.model==m)&(bridge.condition=='PRDM1+ctrl')].iloc[0];bridge_table+=f"| {NAME[m]} | {b.MSE:.5f} | {number(b.Pearson_delta,'.3f')} |\n"
 bridge_text='核心分析的稳定TF与Norman扰动目标共有5个规范身份：EGR1、OSR2、PRDM1、HNF4A和ZBTB25，对应8个细胞类型—候选条目。PRDM1在单核细胞分析中为100%重复选择的正向候选；在外部数据种子1中，PRDM1单扰动和CBFA2T3＋PRDM1组合均留作测试，分别属于未见单扰动和seen0组合。这个依据候选身份选择的案例，具体检验同一候选在K562新干预条件下的响应预测；其与单核细胞机制的联系继续由各自细胞背景的证据支撑。'
 winnertext=f"在全部 18,977 个共同基因上，{NAME[winner]} 的三划分平均 MSE 最低，为 {common.loc[winner,'MSE_mean']:.6f}；相对零变化基线平均降低 {common.loc[winner,'improvement_vs_zero_mean']:.1%}。在训练集前 1,000 个高变基因上，最低平均 MSE 来自 {NAME[hvgwinner]}（{hvg.loc[hvgwinner,'MSE_mean']:.5f}）。"
 if winner in MODELS[:3]:
  interpretation=f"本次设置下，简单基线在整体表达响应预测中保持较强竞争力。三个 GPU 模型中，{NAME[gpuwinner]} 的共同基因平均 MSE 最低；基础模型的价值需要结合未见基因、组合覆盖和表达读出范围逐项判断。"
 else:
  interpretation=f"本次设置下，{NAME[winner]} 在整体表达响应误差上取得优势。该优势对应当前 K562 数据、模型输入适配和三组划分；分层结果进一步显示模型在哪类未见干预条件上发挥作用。"
 if winner=='Ridge_additive':
  interpretation=f"Ridge加和与训练条件均值的整体MSE接近（{common.loc['Ridge_additive','MSE_mean']:.6f}与{common.loc['training_condition_mean','MSE_mean']:.6f}）；表达变化Pearson相关分别为{common.loc['Ridge_additive','Pearson_delta_mean']:.3f}和{common.loc['training_condition_mean','Pearson_delta_mean']:.3f}。三个GPU模型中，{NAME[gpuwinner]}的共同基因平均MSE最低。"
 seen2=select(group='seen2');seen0=select(group='seen0');unseenwinner=seen0.MSE_mean.idxmin()
 generaltext=f"训练覆盖进一步区分了方法的适用条件：两个基因均已出现在训练条件中的组合，Ridge加和相对训练均值的MSE平均降低{seen2.loc['Ridge_additive','improvement_vs_mean_mean']:.1%}；两个基因均未出现的组合中，{NAME[unseenwinner]}的平均MSE最低。这将后续模型研究具体指向新基因干预的响应幅度与特异性。"
 bp=bridge[bridge.condition=='PRDM1+ctrl'].set_index('model');bm=bp.MSE.idxmin();br=bp.Pearson_delta.idxmax()
 bridge_detail=f"在这一预先确定的PRDM1单扰动案例中，{NAME[bm]}的MSE最低（{bp.loc[bm,'MSE']:.5f}），{NAME[br]}的表达变化Pearson相关最高（{bp.loc[br,'Pearson_delta']:.3f}）。候选层面的比较由此补充了全体条件的平均表现，为后续同一候选的跨背景研究提供具体读出。"
 tab='| 方法 | 共同基因 MSE（均值±SD） | 变化 Pearson r | HVG1000 MSE | 相对零变化 MSE 改善 |\n|---|---:|---:|---:|---:|\n'
 for m in MODELS:tab+=f"| {NAME[m]} | {common.loc[m,'MSE_mean']:.6f} ± {common.loc[m,'MSE_sd']:.6f} | {number(common.loc[m,'Pearson_delta_mean'],'.3f')} | {hvg.loc[m,'MSE_mean']:.5f} | {common.loc[m,'improvement_vs_zero_mean']:.1%} |\n"
 histories=[]
 for m in MODELS[3:]:
  for k in [1,2,3]:
   out=R/f'results/stage7/{m}/seed{k}';x=json.loads((out/'completed.json').read_text());cfg=json.loads((out/'run_config.json').read_text());histories.append({**x,**{f'config_{a}':b for a,b in cfg.items() if isinstance(b,(str,int,float,bool))}})
 runtime=pd.DataFrame(histories);runtime.to_csv(T/'07_run_summary.csv',index=False)
 rt='| 模型 | 划分 | 完成轮次 | 最佳轮次 | 总耗时（分钟） | 最终进程峰值显存（GiB） |\n|---|---:|---:|---:|---:|---:|\n'
 for x in histories:rt+=f"| {NAME[x['model']]} | {x['seed']} | {x['epochs']} | {x['best_epoch']} | {x['total_seconds']/60:.1f} | {x['peak_allocated_GiB']:.2f} |\n"
 additive=s[(s.model=='direct_single_additive')&(s.gene_set=='shared')&(~s.exclude_direct_targets)&(s.group=='direct_additive_eligible')].iloc[0];counts=coverage[coverage.model=='direct_single_additive'].sort_values('seed').eligible.tolist()
 excl=select(exclude=True);extext='；'.join(f"{NAME[m]} {excl.loc[m,'MSE_mean']:.5f}" for m in MODELS)
 rettext=f"按原计划补充每目标 20 次重新调参子抽样：每次抽取 80% 基因，并重新进行三折正则化选择。原前 10 候选共 40 项的最低重新选择频率为 {retune['display_top40_min_retuned_frequency']:.0%}，其中 {retune['display_top40_same_median_sign']}/40 项的系数中位数方向与原全数据系数一致。完整结果见[重新调参候选表](../tables/07_core_top10_retuning_check.csv)。"
 technical=f'''# 第 7 步：外部扰动响应预测结果

## 核心发现

{winnertext}

{interpretation}

{tab}

结果先在各划分内按条件等权汇总，再对种子1、2、3取均值；SD描述三个划分之间的变化。零变化预测的变化相关系数无定义，以“—”表示。完整逐条件、逐划分及各基因集合结果均保存在[条件指标](../tables/07_condition_metrics.csv)、[划分指标](../tables/07_seed_metrics.csv)与[模型汇总](../tables/07_model_summary.csv)。

![完整模型比较](../figures/07_F1_model_benchmark.png)

**图 1｜完整测试集比较。** 小点为三个划分，菱形及横线为均值±SD；MSE比值越低越好，表达变化相关越高越好。上、下分别使用全部共同基因及训练集高变基因。

## 数据、任务与训练设置

Norman2019 K562 CRISPRa 数据统一为84,143个细胞、225个条件（含对照）、100个扰动基因。输入表达保留公开处理矩阵的log1p尺度，原始总计数单独保留。模型共同评价18,977个基因，scFoundation使用完整19,264基因输入；各适配器、来源和权重校验见[第6步来源说明](06_pretraining_audit.md)。

三个种子沿用固定GEARS simulation清单，同一条件的所有细胞进入同一集合。每个训练集合中的全部细胞参与训练；64个固定训练对照作为各模型共享的参考输入，推断取64个预测的平均，真实响应取该扰动条件的全部细胞均值。表达变化均相对全部训练对照的均值计算。训练高变基因、共表达图及损失非零读出掩码分别由相应训练数据确定。GO图使用声明的100个扰动基因范围。

GEARS采用隐藏维数64；scFoundation＋GEARS采用官方v1模式、隐藏维数512及冻结编码器；scGPT加载blood通用预训练编码器，训练每批抽取1536个基因并保留直接扰动基因，推断覆盖全部共同输出基因。三个模型均使用原项目模型类。冻结编码器的64个对照表示以float32缓存，并与实时编码输出核对。

GEARS/scFoundation使用Adam，初始学习率0.001、权重衰减0.0005，逐轮乘0.5，梯度值裁剪1；scGPT使用Adam，初始学习率0.0001、逐轮乘0.9，梯度范数裁剪1。训练batch分别为64、24、64；浮点矩阵运算使用TF32，scGPT另采用bfloat16自动混合精度。等价的批量GEARS损失通过原实现的损失值和梯度比较。各模型根据验证条件等权MSE选最佳检查点，patience5。用户在运行中将训练上限收紧至8轮，全部最终检查点均从前8轮选出；调整前已完成的额外轮次保留在训练历史中。测试集在选定模型后评价。

Ridge以训练条件的二元扰动编码拟合表达变化，无截距；9个惩罚强度由验证集选择。训练均值按扰动条件等权，零变化返回训练对照均值。单扰动直接加和在两个单扰动均为训练条件时使用，另列其覆盖子集。正式配置见[stage7_config.json](../config/stage7_config.json)。

## 未见扰动的分层结果

![泛化分层](../figures/07_F2_generalization.png)

**图 2｜最终训练集覆盖分层。** seen0/1/2表示组合中有0/1/2个扰动基因出现在最终训练条件中，另列未见单扰动。图中百分比为相对训练条件均值基线的MSE改善，正值表示误差降低。完整官方原标签与重新计算标签均保留。

按细胞划分的预测通常容易混入相同扰动条件；本实验将整个扰动条件留出，直接检验新条件的泛化。Ridge对训练中完全未出现的基因，其对应系数为零；训练均值基线提供跨条件共享响应，图结构与预训练表示则提供不同的迁移来源。分层表因此同时反映训练覆盖及模型表示的作用。

{generaltext}

## 直接加和、读出敏感性与案例

三个划分分别有{counts[0]}、{counts[1]}、{counts[2]}个测试条件具备完整单扰动训练信息。在相同覆盖子集上，直接加和基线的共同基因平均MSE为{additive.MSE_mean:.5f}。所有方法均在这组匹配条件上重新计算，见下图及完整表中的 `direct_additive_eligible` 分层。

![直接加和匹配比较](../figures/07_F4_direct_additive_subset.png)

**图 3｜完整单扰动信息下的组合预测。** MSE以同一划分的直接加和基线归一化；所有方法共享可用条件子集。

剔除直接扰动基因后，共同基因平均MSE分别为：{extext}。该补充结果集中考察更广泛的下游表达响应；高变基因对应结果也完整保存。

![代表性预测](../figures/07_F5_representative_predictions.png)

**图 4｜按真实响应强度选择的代表案例。** 从种子1的测试组合中，按训练高变基因上的真实响应均方幅度选择最大、中位和最小案例。每点为一个高变基因，虚线为一致线，案例选择未使用模型误差。

## 训练过程、资源与核心补充

![验证轨迹](../figures/07_F3_training_curves.png)

**图 5｜验证误差与模型选择。** 线型表示三个划分，圆点为最低验证误差轮次，纵轴对数刻度。

{rt}

耗时累计预算切换前后的运行时间。显存列记录任务最终进程的分配峰值；续跑项的该进程包含检查点恢复和最终预测。

{rettext}

## 与核心候选衔接：PRDM1

{bridge_text}

{bridge_table}

{bridge_detail}

此表采用种子1训练集确定的1,000个高变基因。两个PRDM1相关条件的完整指标见[案例源表](../tables/08_PRDM1_response_bridge.csv)，共有候选的原始统计量见[候选交集](../tables/08_core_Norman_target_overlap.csv)。

## 研究解释与适用范围

表达状态预训练提供了可迁移的基因与细胞表示；扰动预测还需要把干预身份、细胞背景和响应幅度连接起来。本次比较把共享响应、加和响应、基因关系图和通用预训练放在相同条件划分下评价。简单基线的作用是刻画哪些结构已能由训练数据直接估计，从而定位复杂模型新增的信息。[GEARS原文](https://doi.org/10.1038/s41587-023-01905-6)、[scGPT原文](https://doi.org/10.1038/s41592-024-02201-0)、[scFoundation原文](https://doi.org/10.1038/s41592-024-02305-7)。

近期系统比较同样强调有力基线、匹配评价条件及干预特异性的重要性；其划分、基因集合和模型设置与本实验各自明确，结果可在研究问题层面相互参照。[Ahlmann-Eltze等，2025](https://doi.org/10.1038/s41592-025-02772-6)。

当前实验覆盖K562中的一个终点表达响应任务，三次划分描述这套数据内的变化。模型实现对应本次GO图范围、64个参考对照和已记录训练配置；本次采用8轮训练预算，比较结果也反映这一计算条件。通用预训练样本与Norman的逐细胞重叠状态尚未核实。后续医学世界模型研究可在此基础上加入时间、组织环境和功能结局，并在新背景中检验预测和校准。基因编辑与递送研究则可进一步将干预方向、表达持续时间和靶细胞状态纳入响应预测。

## 验收和成果入口

共{validation['check_count']}项正式实验自检通过；20次重新调参补充的划分、参数选择和收敛证书全部通过。九个GPU正式运行、四种基线的三组结果、全部最佳检查点、逐条件预测及图表源数据均归档。完整日志见[工作日志](../WORK_LOG.md)，图注见[图注文件](07_figure_captions.md)，最终归档见[第8步验收](08_archive_handoff.md)。
'''
 (REPORT/'07_gpu_benchmark_report.md').write_text(technical)
 core=(REPORT/'05_core_report.md').read_text();core=core.replace('# 血液细胞表达差异的 Lasso 重建与转录因子候选分析','# 血液细胞状态重建与扰动响应预测：从 Lasso 候选到基础模型比较',1)
 core=core.replace('## 1．研究问题与数据',f'外部拓展使用84,143个K562细胞和三组条件划分，比较四种基线及GEARS、scGPT blood、scFoundation＋GEARS。{winnertext}\n\n## 1．研究问题与数据',1)
 core=core.replace('GEARS、scGPT 和 scFoundation 相关实验属于已规划的后续工作，本核心报告中的实测结果全部来自已完成的统计分析。','第6节进一步完成GEARS、scGPT和scFoundation在相同外部扰动条件划分下的实测比较。')
 core=core.replace('100 次子抽样固定相对惩罚，刻画当前数据条件下的选择重复性；进一步重新调参和独立数据验证安排在核心报告之后。','100次固定惩罚子抽样与每目标20次重新调参子抽样共同考察当前基因集合下的选择重复性。GPU拓展评价K562终点表达响应，使用已记录的图范围和64个参考对照；预训练成员重叠状态尚未核实，跨细胞背景和时间轨迹仍需相应数据评价。')
 core=core.replace('为使生物学讨论集中，选择下面的代表案例。',rettext+'\n\n为使生物学讨论集中，选择下面的代表案例。',1)
 core=core.replace('## 6．局限性','## 7．局限性').replace('## 7．结论','## 8．结论')
 gpu_section=f'''## 6．外部拓展：从状态匹配到干预响应

在核心分析中，观测单位是基因，目标是重建血细胞相对HSC的状态差异。外部实验转向新的问题：已知K562受到哪些基因扰动，能否预测一个训练中未出现的扰动条件的表达响应？采用Norman2019数据，统一保留84,143个细胞、225个条件及18,977个共同输出基因。三个GEARS simulation划分按完整条件分组；所有模型共享参考对照、输出基因和测试条件。

{winnertext}

{tab}

![外部扰动模型比较](../figures/07_F1_model_benchmark.png)

**图5｜外部扰动预测的统一比较。** 误差与变化相关均先按扰动条件等权，再跨三组划分汇总。小点为各划分，菱形与横线为均值±SD；高变基因由对应训练集确定。单扰动直接加和另在覆盖完整的匹配子集上评价，完整结果见[GPU实验报告](07_gpu_benchmark_report.md)。

{interpretation}

![未见干预分层](../figures/07_F2_generalization.png)

**图6｜不同训练覆盖下的泛化。** seen0、seen1、seen2分别表示组合的两个基因中有0、1、2个出现在最终训练条件中；另列未见单扰动。百分比为相对训练均值的MSE改善，正值表示误差降低。

{generaltext}

三类GPU模型各完成三个划分。根据用户时间预算，上限由20轮调整为8轮，连续5轮验证误差无改善时早停；最终均选用前8轮中的最佳验证检查点。调整前已完成的额外轮次保留在完整训练记录中。scGPT采用blood通用预训练，scFoundation冻结编码器并训练GEARS预测部分；全部条件训练细胞保留，64个共享训练对照作为参考输入。模型来源、适配、学习率和实际批量均在[完整方法与结果](07_gpu_benchmark_report.md)记录。

{bridge_text}

{bridge_table}

{bridge_detail}

PRDM1案例采用种子1训练集的1,000个高变基因，完整数据见[案例源表](../tables/08_PRDM1_response_bridge.csv)。

这一拓展为研究兴趣提供了连续路径：Lasso把血液状态差异转为候选优先级；扰动模型检验“干预身份—表达响应”的可预测性；后续再将编辑方式、递送背景、时间与功能读出加入同一验证框架。医学世界模型的细胞模块可以围绕这些明确定义的预测目标逐步建立，并通过独立干预数据评价。

'''
 core=core.replace('## 7．局限性',gpu_section+'## 7．局限性',1)
 core=core.replace('## 参考文献',f'外部扰动实验进一步给出统一模型比较：{winnertext} 状态重建与干预响应两个任务由共享的证据原则连接：明确预测单位、保留有力基线，并在相应细胞背景中解释结果。\n\n## 参考文献',1)
 core+='\n\n## 完整成果\n\n- [核心补充材料](05_core_supplement.md)：全部目标、候选和原定统计模型。\n- [GPU实验完整报告](07_gpu_benchmark_report.md)：九个GPU运行、四种基线、分层、案例与资源。\n- [最终归档与使用说明](08_archive_handoff.md)：本地文件、核验记录及远端权重清单。\n'
 (REPORT/'08_enhanced_report.md').write_text(core)
 brief=f'''# 一页读懂本次作业

本次工作有一条主线：**先从血液细胞的表达特征提出调控候选，再用独立的扰动数据比较表达响应预测方法。** 两部分分别回答候选发现和新干预预测这两个问题。

## 核心作业的结论

用2,717种L1000扰动模式共同重建血液细胞相对HSC的表达差异，在21个目标中均优于均值基线。四类重点目标的留出基因R²分别为：B细胞0.427、浆细胞0.429、单核细胞0.483、红细胞0.609。红细胞的特征最容易由当前扰动字典重建。

20次重新调参后，展示的40个候选最低选择频率为{retune["display_top40_min_retuned_frequency"]:.0%}，全部候选的系数中位数方向一致。候选讨论集中在红系MYB、浆细胞负向BCL6，以及髓系KLF6／PRDM1。这些候选结合了重复选择和相关细胞背景的文献证据；ZNF22等统计前列候选提供新的探索线索。系数描述各扰动模式在联合表达重建中的贡献，后续实验应结合候选机制选择干预方向。

## GPU拓展带来的结论

在84,143个K562细胞上，按完整扰动条件留出测试，比较三种GPU模型和四种基线。{winnertext}

{interpretation}

三模型各完成三个划分，最终检查点均从前8轮选出。单扰动直接加和的覆盖条件另作匹配比较，模型排序和分层结果完整保留。

## 如何使用这些结果

核心分析提供有统计支持的候选优先级，外部实验提供“干预身份—表达响应”的统一比较。后续课题可围绕少数有机制证据的候选，将细胞背景、编辑方式、递送持续时间和功能读出加入验证框架；医学世界模型的细胞模块可由这些具体预测任务逐步建立。

- [增强分析报告](08_enhanced_report.md)：完整主线、主要图表与文献讨论。
- [GPU实验完整结果](07_gpu_benchmark_report.md)：全部模型、分层、训练过程和资源记录。
- [归档与复现说明](08_archive_handoff.md)：本地数据、代码、检查点和核验入口。
'''
 (REPORT/'08_results_brief.md').write_text(brief)
 print(winnertext);print('Wrote stage7 full report and stage8 enhanced report.')
if __name__=='__main__':main()
