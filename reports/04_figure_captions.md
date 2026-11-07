# 第 4 阶段图注

## 04_F1_model_comparison

四目标共享五折基因划分的模型比较。圆点为最低内层 MSE 选择，菱形为 1-SE 选择；灰色短线为训练均值基线。指标由全部折外预测合并计算，连线连接同一模型的两种选择。六种统计模型全部保留。

[PDF](../figures/04_F1_model_comparison.pdf) · [SVG](../figures/04_F1_model_comparison.svg)

## 04_F2_stability_candidates

四目标按选择频率、系数绝对值中位数和全数据系数幅度排序的前 10 个 TF。每次无放回抽取 80% 基因，共 100 次；固定全数据 1-SE 相对惩罚，重新计算训练尺度和 λmax。点为系数中位数，线段为子抽样系数的 2.5%–97.5% 分位数，右侧为选择频率。彩色圆点满足稳定优先规则；灰色叉号为其余展示候选。分位线描述重采样分布。

[PDF](../figures/04_F2_stability_candidates.pdf) · [SVG](../figures/04_F2_stability_candidates.svg)

## 04_F3_cross_target_candidates

四目标前 10 个 TF 的并集。红蓝表示全数据系数的正负方向，颜色强度表示该目标中 100 次子抽样的选择频率；白点标记满足稳定优先规则。每个扰动保留原始构建体标签。完整稳定 TF 身份层面的交集与 Jaccard 表另存。

[PDF](../figures/04_F3_cross_target_candidates.pdf) · [SVG](../figures/04_F3_cross_target_candidates.svg)

## 04_F4_sensitivity_summary

A：分别评价原差值与逐基因标准化差值两个响应定义，报告各自响应单位下的折外 R²；模型主排名使用原差值任务。B：在相同的 965 个 landmark 留出读出上，比较全基因训练与 landmark 专项训练。C：前 10 个 TF 在 Elastic Net、landmark Lasso、逐基因标准化 Lasso 的全数据 1-SE 模型中保持活跃且方向相同的比例。三种候选对照以不同形状区分。

[PDF](../figures/04_F4_sensitivity_summary.pdf) · [SVG](../figures/04_F4_sensitivity_summary.svg)

## 04_F5_correlated_predictor_groups

Lasso 与 Elastic Net 全数据 1-SE 标准化系数逐项对照。灰点为全部 2,717 个扰动，彩点为稳定性排序的前 10 个 TF；黑色空心圈标记绝对 Pearson 相关至少 0.80 的多成员连通组成员，虚线为两模型系数相同。各目标标注展示 TF 中最正和最负的两项。完整相关组成员、预测贡献、活跃数和 0.70/0.90 阈值对照保存在源表。

[PDF](../figures/04_F5_correlated_predictor_groups.pdf) · [SVG](../figures/04_F5_correlated_predictor_groups.svg)

## 04_F6_self_readout_check

逐候选删除与其 HGNC 规范身份匹配的自身表达读出行，在固定全数据 Lasso 1-SE 相对惩罚下重新拟合。横轴为删除前后系数差除以基线绝对系数；基线和删除模型均采用 10⁻⁶ 的更严格对偶间隙阈值。圆点表示保留活跃方向，叉号表示方向或活跃状态改变；无自身读出的候选单列。完整 Pearson 变化、匹配行与证书另存。

[PDF](../figures/04_F6_self_readout_check.pdf) · [SVG](../figures/04_F6_self_readout_check.svg)

## 04_F7_size_matched_control

各 965 个读出的 landmark 集合与三次非 landmark 随机抽样对照。各集合独立完成共享划分下的 5×3 Lasso 调参，展示最低误差模型。非 landmark 抽样精确匹配每个外层测试折的 landmark 数量，内层通过与原名单求交保留划分关系；内层数量随抽样变化。各集合的训练和测试基因身份均不同。三次抽样用于描述固定总数下基因组成的敏感性，点间范围为本次抽样范围。

[PDF](../figures/04_F7_size_matched_control.pdf) · [SVG](../figures/04_F7_size_matched_control.svg)
