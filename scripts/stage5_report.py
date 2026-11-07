"""Build the core report from approved results and a curated literature evidence file.

This script performs no model fitting and does not modify stage 1-4 results.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT/'tables'
FOCUS = ['B cell','plasma cell','monocyte','erythrocyte']
ZH = dict(zip(FOCUS,['B 细胞','浆细胞','单核细胞','红细胞']))


def table(headers, rows):
    return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)


def references():
    items=json.loads((ROOT/'data/references/stage5/pubmed_exact_records.json').read_text())['result']['summary']['PubmedArticleSet']['PubmedArticle']
    by_id={r['MedlineCitation']['PMID']:r for r in items}
    order=[('R01','29195078'),('R03','20686118'),('R04','16484593'),('R05','15240705'),('R06','24385430'),('R07','27539853'),('R08','7542362'),('R09','33283989'),('R10','37499029'),('R12','37592036'),('R13','38409223'),('R14','38844628')]
    rows=[]
    for rid,pmid in order:
        r=by_id[pmid];a=r['MedlineCitation']['Article'];ids=r['PubmedData']['ArticleIdList']['ArticleId'];ids=ids if isinstance(ids,list) else [ids]
        authors=a['AuthorList']['Author'];authors=authors if isinstance(authors,list) else [authors]
        rows.append({'id':rid,'pmid':pmid,'title':a['ArticleTitle'].rstrip('.'),'authors':'; '.join(x.get('LastName','')+' '+x.get('Initials','') for x in authors),
          'year':a['Journal']['JournalIssue']['PubDate']['Year'],'journal':a['Journal']['Title'],'doi':next(x for x in ids if x.startswith('10.')),
          'url':f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/','evidence_access':'PubMed exact-ID metadata; abstract and primary-source indexed passages; some full-text pages inaccessible'})
    rows += [
      {'id':'R02','pmid':'29425488','title':'The Human Transcription Factors','authors':'Lambert SA; Jolma A; Campitelli LF; et al.','year':'2018','journal':'Cell','doi':'10.1016/j.cell.2018.01.029','url':'https://pubmed.ncbi.nlm.nih.gov/29425488/','evidence_access':'PubMed indexed record and official Human TFs download page'},
      {'id':'R11','pmid':'31395745','title':'Exploring genetic interaction manifolds constructed from rich single-cell phenotypes','authors':'Norman TM; Horlbeck MA; Replogle JM; et al.','year':'2019','journal':'Science','doi':'10.1126/science.aax4438','url':'https://pubmed.ncbi.nlm.nih.gov/31395745/','evidence_access':'Primary article indexed abstract and PMC passages'},
      {'id':'R15','pmid':'','title':'Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines','authors':'Ahlmann-Eltze C; Huber W; Anders S','year':'2025','journal':'Nature Methods','doi':'10.1038/s41592-025-02772-6','url':'https://www.nature.com/articles/s41592-025-02772-6','evidence_access':'Publisher full text: abstract and main results inspected'},
      {'id':'R16','pmid':'','title':'TxPert: using multiple knowledge graphs for prediction of transcriptomic perturbation effects','authors':'Wenkel F; Tu W; Masschelein C; et al.','year':'2026','journal':'Nature Biotechnology','doi':'10.1038/s41587-026-03113-4','url':'https://www.nature.com/articles/s41587-026-03113-4','evidence_access':'Publisher full text: abstract, introduction, metric design, architecture inspected; published 2026-05-01'}]
    rows=sorted(rows,key=lambda r:r['id'])
    (ROOT/'data/references/stage5/references.json').write_text(json.dumps({'checked_date':'2026-11-07','references':rows},ensure_ascii=False,indent=2)+'\n')
    pd.DataFrame(rows).to_csv(TAB/'05_references.csv',index=False,encoding='utf-8-sig')
    return {r['id']:r for r in rows}


def main():
    refs=references()
    def cite(rid):return f'[{rid}：{refs[rid]["authors"].split(";")[0]}，{refs[rid]["year"]}]({refs[rid]["url"]})'
    metrics=pd.read_csv(TAB/'03_model_metrics.csv');stable=pd.read_csv(TAB/'04_stability_all_TFs.csv')
    top=pd.read_csv(TAB/'04_top10_TF_candidates.csv');models=pd.read_csv(TAB/'04_model_metrics.csv')
    evidence=pd.read_csv(TAB/'04_candidate_evidence_matrix.csv');overlap=pd.read_csv(TAB/'04_stable_TF_overlap.csv')
    m=metrics[metrics.model=='Lasso_min'].set_index('target');base=metrics[metrics.model.str.contains('baseline',case=False)].set_index('target')
    counts=stable.groupby('target').stable_priority.sum()
    performance=table(['目标','Lasso R²','Pearson r','RMSE','相对均值基线 RMSE 降低','稳定优先 TF 条目'],[[ZH[t],f'{m.loc[t,"R2"]:.3f}',f'{m.loc[t,"Pearson"]:.3f}',f'{m.loc[t,"RMSE"]:.4f}',f'{100*(1-m.loc[t,"RMSE"]/base.loc[t,"RMSE"]):.1f}%',int(counts[t])] for t in FOCUS])
    chosen=[('erythrocyte','MYB','R03;R04','红系机制支持','人 CD34+ 细胞中的 MYB 下调影响红系命运，并涉及 KLF1、LMO2；小鼠研究显示其作用随红系阶段变化。','将 MYB 相关匹配优先联系红系发育程序，检验阶段与表达持续时间。'),
      ('plasma cell','BCL6','R05','浆细胞机制支持','小鼠 Bcl6 直接抑制 Prdm1；缺失后浆细胞发育加快。','负向匹配与解除 BCL6 抑制的方向相呼应；在 B 到浆细胞背景核对转变。'),
      ('monocyte','KLF6','R06;R07','邻近髓系背景支持','人/鼠巨噬细胞研究支持 KLF6 的炎症调控；后续工作涉及 PRDM1–BCL6 轴。','优先考察炎症功能状态，并区分其与单核细胞谱系形成的关系。'),
      ('monocyte','PRDM1','R07','邻近髓系背景支持','巨噬细胞 KLF6–PRDM1–BCL6 轴提供与 KLF6 共同出现的研究线索。','在对应背景检验两者共同响应与先后关系；本次未直接估计调控边。'),
      ('B cell','ZNF22','','本次定向检索尚未定位直接证据','预定基因—B/淋巴/造血检索返回空；保留统计排名。','作为探索性候选，在独立 B 细胞数据中优先复核。'),
      ('plasma cell','OVOL2','','共同模式候选','本次定向检索命中一篇镉暴露/糖代谢研究，未作为血液分化机制证据。','比较其跨细胞共同表达程序与目标特有程序。'),
      ('erythrocyte','NFYB','','共同模式候选','本阶段按四目标稳定交集讨论；未新增 NFYB 特异机制结论。','比较跨目标贡献，查明是否主要追踪共同状态特征。')]
    literature=[]
    for target,gene,rids,level,support,direction in chosen:
        row=stable[(stable.target==target)&(stable.perturbation==gene)].iloc[0]
        match=evidence[(evidence.target==target)&(evidence.perturbation==gene)]
        literature.append({'target':target,'TF':gene,'TF_stability_rank':int(row.TF_stability_rank),'selection_frequency':row.selection_frequency,'coefficient_standardized':row.full_coefficient,
          'stable_priority':bool(row.stable_priority),'source_ids':rids,'evidence_scope':level,'literature_support_zh':support,'interpretation_and_validation_zh':direction,
          'self_readout_check_status':match.iloc[0].self_readout_status if len(match) else 'not_in_prespecified_top10_check',
          'source_urls':';'.join(refs[i]['url'] for i in rids.split(';') if i)})
    pd.DataFrame(literature).to_csv(TAB/'05_TF_literature_evidence.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
    candidate_table=table(['目标','代表 TF','本目标 TF 稳定性名次','系数方向','抽样选择频率','讨论定位'],[[ZH[r['target']],r['TF'],r['TF_stability_rank'],'正' if r['coefficient_standardized']>0 else '负',f'{r["selection_frequency"]:.0%}',r['evidence_scope']] for r in literature])
    model_table=table(['模型']+[ZH[t] for t in FOCUS],[[name]+[f'{models[(models.target==t)&(models.model==key)&(models.selection=="min")].R2.iloc[0]:.4f}' for t in FOCUS] for name,key in [('Lasso','Lasso'),('Ridge','Ridge'),('Elastic Net','ElasticNet'),('PCR','PCR'),('PLS','PLS'),('非负 Lasso','NonnegativeLasso')]])
    pair=overlap[(overlap.left_target=='B cell')&(overlap.right_target=='plasma cell')].iloc[0]
    def short_authors(value):
        names=value.split('; ')
        return ', '.join(names[:3])+(', et al' if len(names)>3 else '')
    bib='\n\n'.join(f'{rid}. {short_authors(r["authors"])}. **{r["title"]}.** *{r["journal"]}*. {r["year"]}. [DOI](https://doi.org/{r["doi"]})'+('；GEARS 于 2023 年在线发表，卷期年份为 2024。' if rid=='R12' else '') for rid,r in refs.items())
    report=f'''# 血液细胞表达差异的 Lasso 重建与转录因子候选分析

## 摘要

不同血液细胞的功能差异伴随着不同的基因表达程序。本研究以造血干细胞（HSC）为参考，将目标细胞的表达差值与 L1000 过表达扰动谱进行匹配，筛选能够共同重建目标模式的转录因子（TF）候选。对 11,922 个共同基因和 2,717 个扰动建立 Lasso 模型，在 21 个非 HSC 类型中进行筛查，重点讨论 B 细胞、浆细胞、单核细胞和红细胞。四目标留出基因上的 R² 分别为 **0.427、0.429、0.483、0.609**，红细胞表现最佳。100 次基因子抽样得到 41、43、52、42 个稳定优先 TF 条目。文献支持将红系 MYB、浆细胞负向 BCL6、单核细胞 KLF6／PRDM1 作为重点解释对象；ZNF22 等统计前列候选提供进一步研究线索。研究形成了从表达模式匹配、重复选择到机制解释的候选优先级，可为细胞背景明确的干预验证和后续扰动预测提供起点。

## 1．研究问题与数据

研究问题是：**哪些 L1000 扰动的表达变化模式，与目标血液细胞相对于 HSC 的表达差异相匹配？** 血液平均表达谱给出目标模式，L1000 提供已观测的扰动模式，回归模型估计各模式在联合重建中的权重。候选 TF 的价值来自统计匹配与相关生物学证据的结合。

作业提供 `bloodmeanexpr.csv` 和 `L1000MCF7_OEgeneexpression.h5`。后者为 MCF7 背景的过表达扰动字典。L1000 平台直接测量一组 landmark 基因，并推断更广的转录组读出；本研究使用课程给定矩阵，并通过官方注释识别其中的直接测量基因。{cite('R01')}

| 数据检查 | 结果与处理 |
|---|---|
| 血液矩阵 | 原 20,001 行、22 类平均谱；重复基因保留首条后为 19,994 个基因 |
| L1000 矩阵 | 原 12,328 行、2,717 个扰动；去重后为 12,327 个基因 |
| 基因对齐 | 按原始符号精确匹配，共 11,940 个；去除血液全部类型中均为零的 18 个，最终 11,922 个 |
| 缺失与数值 | 冻结输入均为有限值；血液表达保留输入尺度，L1000 保留正负扰动分数 |
| TF 与直接测量读出 | 291 个 TF 扰动条目、290 个规范身份；共同基因中 965 个为 landmark |
| 类型范围 | 所有 21 个非 HSC 类型均建模；零值率 94.00% 的髓系树突状细胞（mDC）单列附加结果 |

TF 资格依据 Human TFs 目录，并以 HGNC 核对身份。不同构建体保留原标签，例如 `DDIT3` 和 `(DDIT3)` 为独立模型列；跨目标集合比较采用规范基因身份。{cite('R02')}

## 2．表达结构与模型方法

探索分析使用 HSC 和 20 个主要目标的平均表达谱。按跨细胞类型方差选择高变基因，逐基因标准化；前 100 个用于热图，前 1,000 个用于 PCA 和聚类。主聚类采用欧氏距离和 Ward 连接，并以相关距离作对照。PCA 前两个主成分解释 **51.54%** 的变异，HSC 与造血前体在两种距离下均最接近。

![血液平均表达谱的 PCA](../figures/02_F3_PCA.png)

**图 1｜血液细胞表达结构。** 每个点代表一类平均表达谱，重点目标使用统一颜色。A、B 展示前两个及第一／第三主成分，C 为解释方差。图中基因筛选与标准化用于探索展示，回归使用完整共同基因并在训练折内预处理。

![高变基因热图](../figures/02_F2_HVG100_heatmap.png)

**图 2｜前 100 个高变基因的表达模式。** 行为基因，列为细胞类型，颜色为逐基因 z-score。树由图示矩阵的欧氏距离与 Ward 方法生成。该集合中 56 个基因带有核糖体相关分组注释，提示较强的跨类型变异同时包含共同细胞功能程序。

目标相对 HSC 的参照表达变化具有不同特征：B 细胞突出 MS4A1、CD79A、CD79B；浆细胞突出 JCHAIN、MZB1、SDC1；单核细胞突出 LYZ、S100A8、S100A9；红细胞突出 HBB、HBA1、HBA2、ALAS2。这些来自当前矩阵的观察构成回归解释的背景。

对每个目标定义：

\\[
y_g=E_{{g,\\mathrm{{target}}}}-E_{{g,\\mathrm{{HSC}}}},\\qquad
y_g\\approx\\beta_0+\\sum_{{p=1}}^{{2717}}X_{{gp}}\\beta_p.
\\]

这里每一行是一个基因，每一列是一个扰动的表达模式。Lasso 通过 L1 惩罚筛选贡献模式，目标函数为平方误差的一半加系数绝对值惩罚。正、负系数分别表示该模式在联合组合中的同向和反向贡献。

评价采用固定种子 20261107、外层五折与内层三折。训练折内拟合 X／y 的均值与尺度、确定正则化路径；预测回到原响应单位后评价。内层比较 50 个相对惩罚强度，由 λmax 至 0.001λmax。最低误差参数用于预测表现；全数据独立调参后采用 1-SE 规则选择更简约的候选模型。每个基因均有一次外层留出预测，主要报告 R²、RMSE、Pearson 和 Spearman 相关。

四目标进一步各做 100 次无放回 80% 基因子抽样，每次 9,537 个基因，固定全数据 1-SE 相对惩罚并重新计算训练尺度、λmax，从零系数拟合。选择频率至少 80%、活跃抽样中符号一致性至少 90%、全数据系数非零的 TF 记为稳定优先候选；排序依次依据频率、系数幅度中位数和全数据幅度。

## 3．主要结果：可重建的表达差异与目标特异候选

### 3.1 红细胞表达差异的重建表现最好

全部 21 个目标的 Lasso 均优于训练折均值预测。20 个主要目标的 R² 中位数为 0.450，范围为 0.127–0.609；红细胞为最高，其后为中性粒细胞和巨噬细胞。重点结果如下。

{performance}

R² 评价被留出基因上的表达差异重建。红细胞的 R²=0.609，说明当前字典能捕获其较大部分表达差异结构。四目标相对均值基线的 RMSE 降低 24.3%–37.5%，支持进一步利用模型开展候选分析。

![四目标 Lasso 留出预测](../figures/03_F2_focus_OOF_predictions.png)

**图 3｜四目标的留出基因预测。** 横、纵轴分别为观测与预测表达差值，六边形颜色表示基因数，虚线为完全一致。多数基因接近低变化区，大幅变化的预测存在向中心收缩。每个点的预测均来自未使用该基因训练的外层模型。

### 3.2 候选具有重复选择性，也体现细胞类型差异

四目标分别得到 41、43、52、42 个稳定优先 TF 条目。B 细胞与浆细胞共有 {int(pair.intersection)} 个稳定规范身份，Jaccard={pair.Jaccard:.3f}；四类目标共同出现的稳定身份为 **NFYB 和 OVOL2**。这表明同一扰动字典可以同时识别共享模式与目标相关模式。

![TF 候选的子抽样稳定性](../figures/04_F2_stability_candidates.png)

**图 4｜每目标前 10 个 TF 的条件稳定性。** 点与分位线呈现 100 次系数分布，右侧标注选择频率；方向与全数据候选保持明确对应。图示前 10 项均在 100 次抽样中被选中，完整排序和各项敏感性保存在源表。分位线描述本次子抽样分布。

为使生物学讨论集中，选择下面的代表案例。统计排序保持原值；BCL6 位于浆细胞第 15 名，因其直接机制证据而列入讨论，图 4 的前 10 项展示范围保持不变。

{candidate_table}

### 3.3 模型与读出对照支持怎样的解释

Ridge 在四目标的 R² 为 0.435、0.443、0.496、0.615，较 Lasso 增加 0.0055–0.0139；Elastic Net 的结果也接近 Lasso。本报告以 Lasso 的表达重建和候选解释为主线，全部六类模型的原定比较见附录。

前 10 个 TF 在 Elastic Net 中均保留活跃方向。删除候选自身表达读出后，34 个可匹配条目仍保留方向，最大系数相对变化为 7.05%；其余 6 项缺少自身读出。该检查支持这些候选的联合贡献由更广的表达模式承载。

读出选择对结果有实质影响：仅用 965 个 landmark 基因重新训练，四目标 R² 为 0.005、0.007、0.034、−0.003；使用全部基因训练，再在相同 landmark 留出读出上评价，则为 0.316、0.324、0.328、0.224。三组等数量非 landmark 对照的表现均高于 landmark 专项结果，提示训练覆盖与读出组成需要共同考察。逐基因标准化也会改变候选集合，相关结果完整保留在附录中。

## 4．生物学解释：从匹配模式到可检验的机制

### 4.1 红系 MYB：与发育程序相呼应的统计首位候选

MYB 为红细胞的第 1 位 TF，系数为正，100 次抽样均被选中。Bianchi 等在人 CD34+ 造血细胞中研究了 MYB 下调，提出其通过 KLF1 和 LMO2 参与红系与巨核系命运选择；这一证据为当前红系匹配提供了较直接的机制联系。{cite('R03')}

进一步的红系研究显示，c-Myb 参与谱系承诺、向红细胞母细胞阶段推进及终末细胞分裂，其作用具有发育阶段特点。因此，本结果最适合被解释为 **MYB 相关表达程序与红系差异相匹配**，后续验证应关注目标细胞所处阶段与表达持续时间。{cite('R04')}

### 4.2 B／浆细胞：BCL6 的方向线索与 ZNF22 的探索价值

浆细胞模型中 BCL6 为负向稳定候选，选择频率 100%。Tunyaplin 等发现 Bcl-6 可以直接结合并抑制 Prdm1，小鼠 Bcl6 缺失后浆细胞发育加快。当前负向匹配与解除 BCL6 对浆细胞程序的抑制这一方向相呼应，是连接统计系数与分化机制的清晰案例。{cite('R05')}

B 细胞统计首位为 ZNF22。此次围绕 ZNF22 与 B 细胞、淋巴细胞和造血的定向检索尚未定位直接机制研究，因此将它保留为探索候选。作为参照，EBF1 有小鼠 B 细胞发育的遗传证据，但在当前 B 细胞全数据模型中系数为零。两者提示候选匹配同时受可用扰动背景、目标状态和联合变量组合影响；下一步可将 ZNF22 与 EBF1 这一已知参照共同放到 B 细胞背景中检验。{cite('R08')}

### 4.3 单核细胞：KLF6 与 PRDM1 指向髓系功能状态

单核细胞中 KLF6 和 PRDM1 分列第 2 和第 8 位，均为正向、100% 选择频率。Date 等在人和小鼠巨噬细胞体系中支持 KLF6 参与炎症极化和相关转录调控，为它与髓系表达模式的联系提供了实验基础。{cite('R06')}

后续研究提出巨噬细胞中的 KLF6–PRDM1–BCL6 调控轴。这使两个候选共同出现具有可检验的机制线索：优先考察它们是否反映髓系炎症功能状态，再检验这种状态与当前单核细胞平均谱之间的联系。该文献证据来自巨噬细胞背景，本研究没有直接估计三者之间的调控边。{cite('R07')}

### 4.4 共享候选用于区分共同状态与谱系特征

NFYB 在四目标中均呈负向稳定贡献，OVOL2 均呈正向稳定贡献。本数据的高变基因又包含大量核糖体相关基因，因此可以提出一个后续假设：跨目标共享候选可能部分反映共同的生长、代谢或状态变化程序。该假设需要在独立数据中通过表达模块和细胞状态读出检验。当前保留它们的原始统计排序，并与 MYB、KLF6 等具有具体背景证据的候选分别讨论。

## 5．后续验证与科研连接

### 5.1 先检验候选，再扩展模型

后续验证采用三个连续问题：候选能否在独立血液数据中重现；在合适细胞背景中改变候选后，是否出现预测方向的表达程序；这种变化是否伴随预期的细胞状态或功能变化。分析层面优先使用含供体信息的数据，按供体评价重复性，并对阶段、增殖和应激状态作分层。

实验设计可首先聚焦红系 MYB、B／浆细胞 BCL6 和探索候选 ZNF22、髓系 KLF6／PRDM1。分别在相关发育或功能背景中比较方向明确的调节，设置非靶向和已知谱系参照，结合独立生物重复、单细胞转录组及谱系／功能读出。基因编辑与短期表达调节可帮助区分持续调控和暂时响应；具体干预方向由候选自身的机制证据和背景验证共同确定。

### 5.2 与基因编辑和递送研究的联系

血液系统已有通过调控元件编辑改变红系输出的实例：Frangoul 等报道在人自体造血干／祖细胞中编辑 BCL11A 的红系增强子，以提高胎儿血红蛋白。该研究说明，转录调控研究可以从候选 TF 进一步落到细胞特异调控元件与功能输出。{cite('R09')}

Breda 等研究了靶向 HSC 的 mRNA 脂质纳米颗粒递送，为在体调节造血细胞提供了实验研究方向。结合本作业，递送方向最值得关注的变量是细胞靶向性、表达持续时间和干预后状态分布；对阶段依赖的调控因子，这些因素应与候选身份共同评价。本作业的候选表可用于提出研究问题，具体应用需在对应模型中验证。{cite('R10')}

### 5.3 与基础模型、扰动预测和医学世界模型的联系

本研究已有一套可解释的表达模式匹配基线。面向基础模型，下一层问题是：模型能否利用细胞状态和干预信息预测新的表达响应？Norman 等的单细胞组合扰动研究提供了这一类数据；GEARS 将基因关系图与扰动响应预测结合，scGPT 和 scFoundation 则从大规模细胞数据学习表示。{cite('R11')}；{cite('R12')}；{cite('R13')}；{cite('R14')}

近期研究强调评价设计的重要性。2025 年 Ahlmann-Eltze 等在所测试任务中发现简单基线具有很强竞争力；2026 年 TxPert 工作进一步结合多种关系图，并强调匹配批次对照与扰动身份检索评价。两项工作共同提示，研究重点应放在清楚定义新干预／新背景的泛化问题、保持有力基线，以及检验模型是否抓住干预特异信号。{cite('R15')}；{cite('R16')}

对医学世界模型，本报告提出一个具体连接：将细胞层面的“当前状态—干预—响应”作为可验证模块，再逐步加入时间、组织环境和功能结局。当前血液平均谱与 MCF7 字典可用于候选提出；外部单细胞扰动数据承担响应预测检验。GEARS、scGPT 和 scFoundation 相关实验属于已规划的后续工作，本核心报告中的实测结果全部来自已完成的统计分析。

## 6．局限性

本研究使用每类细胞的一列平均表达，缺少供体与单细胞层面的生物重复，因此表达比较采用描述性口径。扰动字典来自 MCF7，匹配候选向血液细胞迁移需要相应背景证据。评价单位为留出基因，基因间关联及 L1000 推断读出影响评价空间；landmark 结果展示了这一敏感性。线性组合描述当前表达重建，候选系数对应联合贡献。100 次子抽样固定相对惩罚，刻画当前数据条件下的选择重复性；进一步重新调参和独立数据验证安排在核心报告之后。

## 7．结论

L1000 扰动表达模式能够重建血液目标细胞相对 HSC 的一部分表达差异，四类重点目标中以红细胞表现最佳。稳定性筛查给出了具有细胞差异和共享成分的 TF 候选，MYB、BCL6、KLF6／PRDM1 分别提供红系、B／浆细胞和髓系的机制讨论切入点。当前成果的核心价值是形成一条可追溯的研究路径：**用表达重建提出候选，用重复选择确定优先级，再用细胞背景和文献组织可检验的生物学问题。**

## 参考文献

{bib}

## 附录：完整预定模型比较

下表保留四目标所有预定模型的最低误差折外 R²。PCR／PLS 的全数据最低误差模型分别达到本轮 100／20 成分上限，比较对应已设定的搜索范围。

{model_table}

全部 21 个目标的四项指标、两套参数模型、40 个展示候选、逐项敏感性与复现入口见[核心报告补充材料](05_core_supplement.md)。[TF 文献证据表](../tables/05_TF_literature_evidence.csv)保存统计证据、文献背景和验证方向；[工作记录](../WORK_LOG.md)保存关键指令、操作与审核记录。
'''
    (ROOT/'reports/05_core_report.md').write_text(report)
    allmetrics=table(['目标','定位','Lasso R²','RMSE','Pearson','Spearman'],[[r.target,r.analysis_role,f'{r.R2:.4f}',f'{r.RMSE:.4f}',f'{r.Pearson:.4f}',f'{r.Spearman:.4f}'] for _,r in metrics[metrics.model=='Lasso_min'].iterrows()])
    top_table=table(['目标','名次','TF','全数据标准化系数','选择频率'],[[ZH[r.target],int(r.TF_stability_rank),r.perturbation,f'{r.full_coefficient:.5f}',f'{r.selection_frequency:.0%}'] for _,r in top.iterrows()])
    supplement=f'''# 核心报告补充材料

## S1．全目标 Lasso 筛查

{allmetrics}

以上是最低误差模型的留出预测。最低误差与 1-SE 的完整四项指标见[阶段 3 指标](../tables/03_model_metrics.csv)，逐折结果见[分折指标](../tables/03_fold_metrics.csv)，参数及系数见[模型选择](../tables/03_selected_models.csv)与[全扰动系数](../tables/03_all_coefficients_annotated.csv)。

## S2．四目标前 10 个 TF

{top_table}

全部 TF 稳定性见[完整表](../tables/04_stability_all_TFs.csv)。稳定优先条目数为 41／43／52／42；规范身份层面两两重叠见[重叠表](../tables/04_stable_TF_overlap.csv)。浆细胞 BCL6 为第 15 名，因直接机制证据进入正文，其自身读出删除不在预定前 10 项检查中。

## S3．完整敏感性与解释口径

- [六类模型及双参数结果](../tables/04_model_metrics.csv)、[逐折结果](../tables/04_fold_metrics.csv)、[全部参数选择](../tables/04_selected_models.csv)。
- [逐项候选证据矩阵](../tables/04_candidate_evidence_matrix.csv)、[相关组](../tables/04_correlation_groups.csv)、[候选相关伙伴](../tables/04_candidate_correlated_partners.csv)。
- [landmark／标准化响应指标](../tables/04_sensitivity_metrics.csv)、[同一 landmark 测试集的全基因模型结果](../tables/04_matched_landmark_reference_metrics.csv)、[自身读出删除](../tables/04_self_readout_results.csv)。
- [三组等规模对照](../tables/04_size_control_metrics.csv)：每组 965 个非 landmark 基因，精确匹配外层数量，内层名单求交；训练与测试身份均改变，三组结果描述本轮抽样范围。

逐基因标准化使用包括目标类型的跨 21 类细胞标准差，形成一个新的响应定义；其回乘原单位的分数在源表中标为描述性复算，未用于主预测排名。条件稳定性从零初始化；先前初始化试算保存在独立探索目录。全部解释与七组图见[阶段 4 报告](04_multimodel_robustness_report.md)。

## S4．方法、文献与复现

[Lasso 方法](03_lasso_methods.md)、[多模型与稳健性方法](04_methods.md)、[环境与重建命令](../environment/README.md)包含完整参数、版本与计算入口。本阶段复用已有结果，未新增训练、重采样或 GPU 拓展。核心报告可由 `scripts/stage5_report.py` 结合既有表格与文献题录重建。

[文献证据表](../tables/05_TF_literature_evidence.csv)、[参考文献题录](../tables/05_references.csv)、[检索与证据范围](05_literature_notes.md)区分实验背景、统计发现及研究假设。完整题录核验记录位于 `data/references/stage5/`。全文访问受限的文献以已核对摘要与原文检索片段支持相应层级的概述，未声称完成所有论文的全文复核。

## S5．图形成果

核心报告采用已审核的 Python 图形，提供 PDF、SVG 与高分辨率 PNG。用户已取消后续 Prism 制图与迁移。已完成的 [Prism 模型比较 PDF](../figures/05_Prism_model_comparison.pdf)、[可编辑项目](../prism/stage5/05_Prism_model_comparison.pzfx)及迁移输入保留为附加历史成果，后续分析和报告采用现有绘图流程。

## S6．独立审核

[阶段 5 独立审核](05_independent_review.md)覆盖主线、数值、引用与作业要求；[阶段 4 审核](04_independent_review.md)保存模型结果的独立数值核验。用户已批准本阶段并要求进入第 6 阶段。
'''
    (ROOT/'reports/05_core_supplement.md').write_text(supplement)
    print(json.dumps({'report':'reports/05_core_report.md','references':len(refs),'evidence_rows':len(literature),'new_model_fits':0},ensure_ascii=False))


if __name__=='__main__':main()
