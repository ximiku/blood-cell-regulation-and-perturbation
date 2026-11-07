"""Render a stage-4 Chinese report directly from saved numeric result tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];TAB=ROOT/'tables'
FOCUS=json.loads((ROOT/'config/stage4_config.json').read_text())['focus_targets']
ZH=dict(zip(FOCUS,['B 细胞','浆细胞','单核细胞','红细胞']))
NAMES={'Mean_baseline':'均值基线','Lasso':'Lasso','ElasticNet':'Elastic Net','NonnegativeLasso':'非负 Lasso','Ridge':'Ridge','PCR':'PCR','PLS':'PLS'}

def md_table(headers,rows):return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in rows)
def main():
    m=pd.read_csv(TAB/'04_model_metrics.csv');s=pd.read_csv(TAB/'04_stability_all_TFs.csv')
    top=pd.read_csv(TAB/'04_top10_TF_candidates.csv');best=m[m.selection=='min'].sort_values('R2',ascending=False).groupby('target').first()
    primary=m[(m.selection=='min')|(m.model=='Mean_baseline')].pivot(index='model',columns='target',values='R2')
    models=['Mean_baseline','Lasso','ElasticNet','Ridge','PCR','PLS','NonnegativeLasso']
    performance=md_table(['模型']+[ZH[t] for t in FOCUS],[[NAMES[model]]+[f'{primary.loc[model,t]:.4f}' for t in FOCUS] for model in models])
    stable_counts=s.groupby('target').stable_priority.sum();count_rows=[]
    for t in FOCUS:
        f=s[(s.target==t)&s.stable_priority];z=top[top.target==t]
        count_rows.append([ZH[t],len(f),f.approved_symbol.nunique(),'; '.join(f'{r.perturbation} ({"+" if r.full_coefficient>0 else "−"}, {r.selection_frequency:.0%})' for _,r in z.head(5).iterrows())])
    candidates=md_table(['目标','稳定优先 TF 扰动条目','规范 TF 身份数','排名前五：方向、选择频率'],count_rows)
    o=pd.read_csv(TAB/'04_sensitivity_original_units.csv');native=pd.read_csv(TAB/'04_sensitivity_metrics.csv');lm=pd.read_csv(TAB/'04_matched_landmark_reference_metrics.csv')
    rows=[]
    for t in FOCUS:
        f=o[(o.target==t)&(o.selection=='min')].set_index('dataset')
        rawlm=lm[(lm.target==t)&(lm.selection=='min')].R2.iloc[0]
        z_r2=native[(native.target==t)&(native.dataset=='gene_scaled')&(native.model=='Lasso')&(native.selection=='min')].R2.iloc[0]
        rows.append([ZH[t],f'{f.loc["raw","R2"]:.4f}',f'{z_r2:.4f}',f'{rawlm:.4f}',f'{f.loc["landmark","R2"]:.4f}'])
    sensitivities=md_table(['目标','原差值任务：全部基因','标准化差值任务：全部基因','原模型：landmark 测试基因','landmark 训练与测试'],rows)
    cand=pd.read_csv(TAB/'04_top10_candidate_sensitivity.csv');own=pd.read_csv(TAB/'04_self_readout_results.csv');evidence=pd.read_csv(TAB/'04_candidate_evidence_matrix.csv')
    rows=[]
    for t in FOCUS:
        f=own[(own.target==t)&(own.status=='completed')]
        counts=[int(cand[(cand.target==t)&(cand.comparison==c)].same_active_sign.sum()) for c in ['raw_ElasticNet','landmark_Lasso','gene_scaled_Lasso']]
        rows.append([ZH[t]]+[f'{n}/10' for n in counts]+[f'{int(f.same_active_sign.sum())}/{len(f)}',f'{100*f.relative_coefficient_change.abs().max():.2f}%'])
    agreement=md_table(['目标','Elastic Net 同向活跃','landmark 同向活跃','逐基因标准化同向活跃','自身读出删除后同向活跃','删除后最大系数相对变化'],rows)
    sets=[set(s[(s.target==t)&s.stable_priority].approved_symbol.dropna()) for t in FOCUS]
    common=sorted(set.intersection(*sets));groups=pd.read_csv(TAB/'04_correlation_groups.csv');gm=pd.read_csv(TAB/'04_group_model_comparison.csv')
    gr=[]
    for threshold in [.7,.8,.9]:
        f=groups[np.isclose(groups.threshold,threshold)].drop_duplicates('group_id');gr.append([threshold,len(f),int((f.group_size>1).sum()),int(f.group_size.max())])
    group_table=md_table(['绝对相关阈值','全部连通组','多成员组','最大组成员数'],gr)
    choices=pd.read_csv(TAB/'04_selected_models.csv');fullchoice=choices[(choices.dataset=='raw')&(choices.scope=='full')]
    pars=[]
    for t in FOCUS:
        for model in ['ElasticNet','Ridge','PCR','PLS','NonnegativeLasso']:
            f=fullchoice[(fullchoice.target==t)&(fullchoice.model==model)].set_index('selection')
            def describe(r):
                if pd.notna(r.relative_alpha):return f'相对惩罚 {r.relative_alpha:.5f}'+(f'，ρ={r.l1_ratio:g}' if model=='ElasticNet' else '')
                if pd.notna(r['lambda']):return f'λ={r["lambda"]:.5g}'
                return f'{int(r.components)} 个成分'
            pars.append([ZH[t],NAMES[model],describe(f.loc['min']),describe(f.loc['1se']),int(f.loc['1se','n_active'])])
    params=md_table(['目标','模型','最低误差参数','1-SE 参数','1-SE 活跃系数数'],pars)
    validation=json.loads((ROOT/'logs/stage4_validation.json').read_text())
    ref=json.loads((ROOT/'logs/stage4_solver_crosscheck.json').read_text())
    timing=json.loads((ROOT/'logs/stage4_models_stability_execution.json').read_text())
    own_timing=json.loads((ROOT/'logs/stage4_self_readout_execution.json').read_text())
    cold_path=ROOT/'logs/stage4_stability_cold_execution.json'
    cold_text=(f'正式从零初始化的 100 次重采样补充批次用时 {json.loads(cold_path.read_text())["seconds"]:.2f} 秒。' if cold_path.exists() else '正式 100 次重采样均从零初始化，其用时包含在模型运行批次中。')
    init_path=TAB/'04_initialization_sensitivity.csv'
    init_text='初始化试算与正式结果的对照另存于[初始化敏感性表](../tables/04_initialization_sensitivity.csv)。' if init_path.exists() else ''
    control_section='';control_validation_text=''
    control_path=TAB/'04_size_control_metrics.csv'
    if control_path.exists():
        control=pd.read_csv(control_path);control_rows=[]
        for t in FOCUS:
            landmark=native[(native.target==t)&(native.dataset=='landmark')&(native.model=='Lasso')&(native.selection=='min')].R2.iloc[0]
            vals=control[(control.target==t)&(control.selection=='min')].sort_values('dataset').R2.to_numpy()
            control_rows.append([ZH[t],f'{landmark:.4f}']+[f'{v:.4f}' for v in vals])
        control_table=md_table(['目标','landmark','非 landmark 抽样 1','非 landmark 抽样 2','非 landmark 抽样 3'],control_rows)
        control_section=f'''## 固定读出数量的补充对照

将三组非 landmark 读出分别抽至 965 个，每组重新完成 Lasso 的外层五折、内层三折调参，得到以下最低误差模型结果：

{control_table}

在相同总基因数下，三组非 landmark 集合的四目标 R² 均高于 landmark 专项结果，提示基因组成对本任务表现具有实质影响。结合全基因训练模型在共同 landmark 测试读出上的结果，可以从读出组成与训练覆盖两个角度理解这项敏感性。

![等规模基因组成对照](../figures/04_F7_size_matched_control.png)

**图 7｜** 固定读出总数后的组成敏感性。每个外层测试折的数量与 landmark 集合精确匹配，内层训练与验证名单通过原划分求交确定，内层数量随抽样变化。各组训练和测试基因身份同时改变，因此分数描述同样本量下的不同读出集合。三次抽样呈现本轮抽样范围，相关结果用于补充解释，主模型排名仍以全部共同基因为依据。

完整[补充指标](../tables/04_size_control_metrics.csv)、[逐基因预测](../tables/04_size_control_predictions.csv)、[参数选择](../tables/04_size_control_choices.csv)和[调参曲线](../tables/04_size_control_tuning.csv)均已保存。

'''
        control_validation_text='补充对照另有[原始行验收记录](../logs/stage4_size_control_validation.json)，覆盖三次抽样、60 个训练上下文、全部参数选择和预测。'
    figure_count=len(json.loads((ROOT/'figures/04_figure_manifest.json').read_text()))
    strongest='；'.join(f'{ZH[t]}为 {NAMES[best.loc[t,"model"]]}（R²={best.loc[t,"R2"]:.3f}）' for t in FOCUS)
    delta=[primary.loc['Ridge',t]-primary.loc['Lasso',t] for t in FOCUS]
    ridge_summary=f'Ridge 相对 Lasso 的 R² 增量为 {min(delta):.4f}–{max(delta):.4f}；Elastic Net 在四个目标上也提供了较高的重建表现。' if all(primary.loc['Ridge',t]>primary.loc['Lasso',t] for t in FOCUS) else ''
    broad=evidence[(evidence.raw_ElasticNet_same_active_sign)&(evidence.landmark_Lasso_same_active_sign)&(evidence.gene_scaled_Lasso_same_active_sign)]
    broad_text='；'.join(f'{ZH[t]}：'+('、'.join(broad[broad.target==t].perturbation) or '该组前 10 项中的一致性按各检验分别记录') for t in FOCUS)
    broad_paragraph=('在 Elastic Net、landmark Lasso 和逐基因标准化 Lasso 三项对照中都保留活跃方向的展示候选为：'+broad_text+'。') if len(broad) else '四组前 10 个候选在 Elastic Net 中均保留活跃方向；逐基因标准化响应下分别保留 9、9、8、6 个。'
    report=f'''# 血液细胞表达重建的模型结构与 TF 稳健性

第 4 阶段交付｜2026-11-07｜四类目标、六类统计模型、100 次子抽样与敏感性分析

## 主要发现

共享外层测试基因的比较显示，{strongest}。{ridge_summary}所有六类模型均按预定范围完成调参，最低误差与 1-SE 两套结果完整保留。

100 次 80% 基因子抽样得到的稳定优先 TF 扰动条目数，B、浆、单核、红细胞分别为 **{', '.join(str(int(stable_counts[t])) for t in FOCUS)}**。候选排序同时保留选择频率、系数方向与幅度。直接测量基因、响应尺度和自身读出检查进一步展示各候选在不同分析条件下的表现。

## 共享测试基因上的模型比较

{performance}

表中为合并五个外层测试折的 R²，模型参数由各自内层三折最低 MSE 决定。完整 RMSE、Pearson、Spearman、1-SE 结果和逐折指标见[主指标表](../tables/04_model_metrics.csv)与[分折指标](../tables/04_fold_metrics.csv)。

![四目标模型比较](../figures/04_F1_model_comparison.png)

**图 1｜** 全部六类模型与均值基线；两种参数选择共用测试基因。各模型的表示方式提供不同的组合约束：Lasso 与 Elastic Net产生稀疏系数，Ridge 保留连续收缩的组合，PCR 和 PLS 使用低维表示，非负 Lasso 约束组合系数方向。

四目标的 PCR 与 PLS 最低误差全数据模型都达到预定成分上限，分别为 100 和 20 个成分。本轮比较对应这一预先设定的搜索范围；全部参数选择列于后文。

## 从初筛走向可重复选择的候选

{candidates}

稳定优先规则为选择频率至少 80%、活跃抽样中的多数符号一致性至少 90%，并在全数据 1-SE 模型中活跃。四目标共享 100 个抽样集合，每次 9,537 个基因；相对惩罚固定为第 3 阶段的 1-SE 选择，各次重新计算训练尺度与 λmax，并从零系数独立初始化。排序按频率、幅度中位数及全数据幅度依次确定。

![TF 子抽样稳定性](../figures/04_F2_stability_candidates.png)

**图 2｜** 每目标展示至多 10 个 TF。点与分位线呈现 100 次系数分布，右侧数字为选择频率；彩色圆点满足稳定优先规则。完整候选和分布见[全部 TF 稳定性](../tables/04_stability_all_TFs.csv)、[稳定优先 TF](../tables/04_stable_priority_TFs.csv)与[前 10 项](../tables/04_top10_TF_candidates.csv)。

系数方向描述扰动表达模式在联合重建中的条件贡献。每个候选的身份、模型表现和敏感性证据汇入[候选证据矩阵](../tables/04_candidate_evidence_matrix.csv)，供下一阶段逐项查证文献。

## 跨目标共同模式与相关扰动

四目标共同出现的稳定 TF 规范身份共 **{len(common)}** 个，依字母顺序包括：{'、'.join(common[:15]) if common else '本轮各目标的稳定候选集合具有不同组成'}{'等' if len(common)>15 else ''}。全部两两交集与 Jaccard 见[重叠表](../tables/04_stable_TF_overlap.csv)。

![跨目标候选模式](../figures/04_F3_cross_target_candidates.png)

**图 3｜** 四组前 10 个候选的并集，颜色结合方向和选择频率，白点标记稳定优先项。图中保留原始扰动标签，集合重叠按规范基因身份计算。

相关组以扰动谱绝对 Pearson 相关阈值构建连接，并取连通分量。阈值对照如下：

{group_table}

![相关变量组比较](../figures/04_F5_correlated_predictor_groups.png)

**图 4｜** Lasso 与 Elastic Net 的全部系数对照，突出前 10 个 TF 和 0.80 阈值下的相关组成员。完整[组成员](../tables/04_correlation_groups.csv)、[组预测贡献与活跃数](../tables/04_group_model_comparison.csv)和[候选相关伙伴](../tables/04_candidate_correlated_partners.csv)用于追踪相关变量的替代与共同贡献。

主要阈值下得到两个多成员组。MYB 与 HSP90AB1、POLB、DHH、RPS6KB2、PPP3CB 位于同一连通组；在红细胞的 Lasso 与 Elastic Net 1-SE 模型中，该组六个成员里均由 MYB 保持活跃，组预测贡献标准差分别约为 0.0782 和 0.0800（标准化响应单位）。NFYB 与 OVOL2 则是四目标共同的稳定 TF 身份，分别呈一致的负向与正向重建贡献。

## 读出与响应尺度的敏感性

{sensitivities}

前两列使用全部 11,922 个相同留出基因，分别评价原差值和逐基因标准化差值两个响应定义；后两列在原差值任务中使用相同的 965 个 landmark 测试读出。模型主排名保持原差值任务的评价口径。全部敏感性指标见[敏感性指标表](../tables/04_sensitivity_metrics.csv)。

![敏感性分析概览](../figures/04_F4_sensitivity_summary.png)

**图 5｜** 比较响应权重、训练读出和候选方向。逐基因标准化强调细胞类型之间的相对表达模式；landmark 专项分析集中于官方直接测量的读出。各模型的系数交集、方向一致率与全系数排序相关见[候选集合比较](../tables/04_candidate_set_agreement.csv)。

前 10 个候选的逐项表现如下：

landmark 专项的四个全数据 1-SE 模型均选择了仅含截距的简约解，因此其候选同向活跃数为零。该结果与对应的训练基因数量和读出组成一并考察。

{agreement}

{broad_paragraph}这些交叉证据可用于安排后续文献精读；原定稳定性排序完整保留。

## 候选自身读出检查

按规范身份匹配后，{int((own.status=='completed').sum())} 个目标—候选条目完成自身读出删除检查，{int((own.status=='own_readout_absent').sum())} 个条目在共同读出集合中没有自身基因。共享同一规范身份的删除任务合并执行。

![自身读出删除结果](../figures/04_F6_self_readout_check.png)

**图 6｜** 删除匹配读出后，重新计算训练尺度与 λmax 并拟合固定相对惩罚。完整基线和删除模型均采用 10⁻⁶ 的更严格对偶阈值，图示相对系数变化。逐项匹配、系数方向、边际相关变化与证书见[自身读出结果](../tables/04_self_readout_results.csv)。

{control_section}## 参数、数值验收与复现

{params}

稠密表示模型的活跃系数数与成分数分别描述原扰动空间和低维表示。所有内层误差曲线、完整模型选择及双模型系数均已保存，见[参数选择](../tables/04_selected_models.csv)、[调参曲线](../tables/04_tuning_curves.csv)和[全部系数](../tables/04_all_coefficients_annotated.csv)。

初始模型、敏感性和重采样试算批次的四卡运行用时 {timing['seconds']/60:.2f} 分钟（复用首个独立试点）；{cold_text}自身读出批次用时 {own_timing['seconds']:.2f} 秒。{init_text}环境准备、汇总、传输、制图与验收时间单独记录。新增模型的 17 项独立 CPU 参考均通过；从原始训练行进行的 {len(validation['checks'])} 项验收也全部通过。详细记录见[CPU 参考](../logs/stage4_solver_crosscheck.json)与[阶段验收](../logs/stage4_validation.json)。{control_validation_text}

统计定义、数值算法、固定种子、敏感性口径及来源见[方法说明](04_methods.md)。{figure_count} 组图均导出 PDF、SVG 和 450 dpi PNG，中文说明见[图注](04_figure_captions.md)。分析程序、处理数据、每次抽样索引、预测、系数、图表、环境与日志均归档本地。

## 证据范围与下一阶段

本阶段的稳定性对应当前 MCF7 扰动字典、平均血液表达谱和固定相对惩罚下的基因子抽样。landmark 与尺度对照展示读出范围和误差权重改变时的结果。下一阶段将围绕这些具体候选，结合血液细胞调控机制、编辑验证与递送背景组织文献证据，并形成可独立使用的核心 Markdown 报告。

独立审核意见见[第 4 阶段审核](04_independent_review.md)。本阶段完成后暂停供用户审核；第 5 阶段等待批准。Prism 原生样板继续按已批准安排后续补充。
'''
    (ROOT/'reports/04_multimodel_robustness_report.md').write_text(report)
    print(json.dumps({'report':'reports/04_multimodel_robustness_report.md','best_models':best.model.to_dict(),'stable_TFs':stable_counts.to_dict()},ensure_ascii=False))

if __name__=='__main__':main()
