# 第 5 阶段文献检索与证据范围

检索日期：2026-11-07。围绕已冻结候选、基础模型、编辑和递送兴趣定向检索，正文收录 16 篇文献。本轮直接复用科学结果，未增加模型计算或扩展为全领域系统综述。

## 检索与筛选

使用 Life Sciences Literature 的 ncbi-entrez-skill 脚本查询 PubMed，并通过网页检索核对原始论文、PubMed、PMC、期刊官网。题录优先以准确 PMID 批量核对；近似题名自动扩展结果只作发现线索。

| 范围 | 代表查询 | 结果 |
|---|---|---|
| MYB 红系 | c-myb supports erythropoiesis；c-Myb erythropoiesis threshold levels | 人 CD34+ 研究 R03、小鼠红系阶段研究 R04 |
| B／浆细胞 | BCL6 represses BLIMP1 plasma cell differentiation；B-cell differentiation lacking EBF | Bcl6–Prdm1 研究 R05、EBF 缺失研究 R08 |
| 髓系 | KLF6 macrophage polarization；KLF6 BCL6 PRDM1 | R06/R07，明确为巨噬细胞及相关体系 |
| ZNF22 | `ZNF22 AND (B-cell OR lymphocyte OR hematopoiesis)` | 本次查询为空；保留 checked_sources，结论限定为此次定向检索 |
| OVOL2 | `OVOL2 AND (B-cell OR plasma-cell OR hematopoiesis)` | PMID 42748818 为镉暴露与糖代谢研究，排除为血液分化机制证据 |
| 模型与方法 | L1000；Human Transcription Factors；GEARS；scGPT；scFoundation；linear perturbation baselines | R01/R02/R12–R15；另核验 2026-05-01 发表的 TxPert（R16） |
| 编辑／递送 | CRISPR-Cas9 sickle cell beta-thalassemia；in vivo hematopoietic stem cell mRNA delivery | R09/R10，保持研究背景层面 |
| 组合扰动 | Norman genetic interaction manifolds 2019 | R11，用于说明后续响应预测的数据类型 |

## 题录与实质性证据

`data/references/stage5/pubmed_exact_records.json` 保存准确作者、题名、年份、DOI、PMID 及来源信息。**插件将长字符串截断到约 240 字符，AbstractText 不是完整摘要。** 该文件用于题录核验；正文机制概述另依据本次检索返回的原始论文索引片段、可访问摘要和期刊正文。

部分 PubMed／PMC 打开时返回验证页，部分 Blood／ScienceDirect 正文返回 403。访问失败不作为支持证据，也未声称完成全部论文的全文精读。

| 文献 | 实际核对材料及位置 | 支持范围 |
|---|---|---|
| R03 | [PubMed 原文索引](https://pubmed.ncbi.nlm.nih.gov/20686118/)及 Blood 检索摘要的 CD34、KLF1、LMO2、lineage decision 片段 | 人造血细胞 MYB 下调与红系／巨核系命运的联系 |
| R04 | [原始摘要索引](https://pubmed.ncbi.nlm.nih.gov/16484593/)中的 commitment、erythroblast、terminal cell divisions 结果段 | MYB 的阶段特点 |
| R05 | [PubMed 完整摘要](https://pubmed.ncbi.nlm.nih.gov/15240705/)两段结果与准确 DOI | Bcl6 结合／抑制 Prdm1；缺失小鼠和脾 B 细胞的浆细胞发育变化 |
| R06 | [PMC 原文索引](https://pmc.ncbi.nlm.nih.gov/articles/PMC4036156/)的 Discussion central findings 段；[PubMed 题录](https://pubmed.ncbi.nlm.nih.gov/24385430/) | 人／鼠巨噬细胞 KLF6 炎症调控；未扩写为 HSC 向单核分化的直接证明 |
| R07 | [PubMed 结果索引](https://pubmed.ncbi.nlm.nih.gov/27539853/)及[PMC 摘要索引](https://pmc.ncbi.nlm.nih.gov/articles/PMC5076533/)明确涉及 KLF6 提高 PRDM1、降低 BCL6 | KLF6／PRDM1 共同入选的邻近髓系线索；本模型未估计网络边 |
| R08 | [Nature 摘要片段](https://www.nature.com/articles/376263a0)和[PubMed](https://pubmed.ncbi.nlm.nih.gov/7542362/) | 小鼠 EBF 缺失与 B 细胞发育阻断 |
| R09 | [NEJM 原文检索段](https://www.nejm.org/doi/full/10.1056/NEJMoa2031054)，BCL11A erythroid-specific enhancer；PMID 33283989 | 调控元件编辑与红系功能输出的联系 |
| R10 | [PMC 原始论文索引](https://pmc.ncbi.nlm.nih.gov/articles/PMC10567133/)及[PubMed](https://pubmed.ncbi.nlm.nih.gov/37499029/) | HSC 靶向 mRNA 递送研究方向，无配方和操作参数 |
| R11–R14 | 原始论文摘要／期刊索引与题录；DOI 见参考表 | 组合扰动、GEARS、scGPT／scFoundation 的任务定位；本项目相应实验尚未执行 |
| R15 | [Nature Methods 全文](https://www.nature.com/articles/s41592-025-02772-6) Abstract、Main 的双扰动／未见单扰动比较、基线段 | 所测试任务中简单基线具有竞争力 |
| R16 | [Nature Biotechnology 全文](https://www.nature.com/articles/s41587-026-03113-4) Abstract、Results 的批次对照／检索评价和模型结构；日期 2026-05-01 | 多关系图与评价设计的新方向；未自动追加模型实验 |

R01 平台论文和 R02 官方 TF 目录用于方法来源，前期注释原件及校验保留于阶段 1。医学世界模型部分是本报告提出的研究连接，以细胞状态、干预、响应和未来时间维度组织可验证问题。

## 候选范围与审核意见落实

40 个前 10 候选在主图与补充表完整保留。BCL6 为浆细胞第 15 名，因机制证据进入正文，其自身读出删除不在预定检查范围。ZNF22 和共享候选保持探索定位；NFYB 本轮未做单独机制专题检索，共享状态解释为待检验假设。

独立审核指出题录摘要截断问题后，本说明明确了文件用途、实际核验位置与访问范围。统计证据、文献背景和验证方向见[候选证据表](../tables/05_TF_literature_evidence.csv)。
