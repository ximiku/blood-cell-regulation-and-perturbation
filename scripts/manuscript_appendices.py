"""Generate numbered appendix summaries directly from frozen result tables."""
from pathlib import Path
import pandas as pd,numpy as np,json
R=Path(__file__).resolve().parents[1];T=R/'tables';O=R/'manuscript'
def rd(n):return pd.read_csv(T/(n+'.csv'))
def esc(s):
 s=str(s).replace('β', r'$\beta$')
 for a,b in [('&',r'\&'),('%',r'\%'),('_',r'\_'),('#',r'\#')]:s=s.replace(a,b)
 return s
S=rd('02_cell_dictionary').set_index('cell_type').short_label.to_dict();focus=['B cell','plasma cell','monocyte','erythrocyte'];AB=dict(zip(focus,['B','Plasma','Mono','Ery']))
mods=['zero_change','training_condition_mean','Ridge_additive','GEARS','scGPT_blood','scFoundation_GEARS','direct_single_additive'];names=dict(zip(mods,['Zero','Train mean','Ridge add.','GEARS','scGPT','scF+GEARS','Direct add.']))
def val(x,n=3):return '---' if pd.isna(x) else f'{x:.{n}f}'
parts=[]
def add(s):parts.append(s)
def page(title):add(r'\clearpage\subsection{'+title+'}')
def section(letter,title):add(r'\clearpage\section{'+title+r'}\label{app:'+letter+r'}\setcounter{figure}{0}\setcounter{table}{0}\renewcommand{\thefigure}{'+letter+r'\arabic{figure}}\renewcommand{\thetable}{'+letter+r'\arabic{table}}')
def table(caption,headers,rows,align=None,size='small'):
 if align is None:align='l'+'r'*(len(headers)-1)
 add('{\\'+size+r'\linespread{1}\selectfont\setlength{\tabcolsep}{4pt}\renewcommand{\arraystretch}{1.08}\begin{longtable}{'+align+'}\caption{'+caption+r'}\\\toprule '+' & '.join(headers)+r'\\\midrule\endfirsthead\toprule '+' & '.join(headers)+r'\\\midrule\endhead')
 for r in rows:add(' & '.join(str(v) for v in r)+r'\\')
 add(r'\bottomrule\end{longtable}}')
def fig(stem,label,cap,height):add(r'\plate{'+stem+'}{'+label+'}{'+cap+'}{'+str(height)+'}')
# Bibliography uses full source metadata; author truncation is conventional.
refs=rd('05_references');b=[r'\begin{thebibliography}{99}\small']
for _,r in refs.iterrows():
 au=r.authors.split(';');au='; '.join(au[:3])+('; et al.' if len(au)>3 else '')
 b.append(r'\bibitem{'+r.id+'} '+esc(au)+'. '+esc(r.title)+'. '+esc(r.journal)+', '+str(r.year)+'. DOI: '+r'\nolinkurl{'+r.doi+'}.')
b.append(r'\end{thebibliography}');(O/'references.tex').write_text('\n'.join(b))
# Also provide an editable bibliographic database.
bib=[]
for _,r in refs.iterrows():bib.append('@article{'+r.id+',\n title={'+r.title+'},\n author={'+r.authors.replace(';',' and')+'},\n journal={'+r.journal+'},\n year={'+str(r.year)+'},\n doi={'+r.doi+'}\n}')
(O/'references.bib').write_text('\n\n'.join(bib))
section('A','数据与探索分析')
add('本附录给出原始矩阵质控、细胞覆盖及完整探索图。除 mDC 的敏感性分析外，主要 EDA 覆盖 HSC 和 20 个目标。')
q=rd('01_matrix_qc');table('原始输入矩阵的质量与尺度。缺失和无穷值均为零。',['输入','行数','列数','删除重复行','最小值','最大值'],[[esc(r.dataset),int(r.rows),int(r.columns),int(r.duplicate_rows_to_remove),val(r.minimum),val(r.maximum)] for _,r in q.iterrows()])
add('血液去重后 19,994 个基因，L1000 去重后 12,327 个。模型交集去除全血液零值行后为 11,922 个基因；全部 2,717 个扰动列保留，291 个为 TF 条目。重复处理逐条保留审计，矩阵按原符号对齐，HGNC 映射用于注释。')
fig('02_F1_expression_distributions','a:dist','表达覆盖与分布。A：全部 22 类型的零值比例。B：HSC 与重点目标的跨基因累积分布；symlog 为显示坐标，输入表达尺度保持。',.56)
page('完整高变基因表达图')
fig('Supp_HVG100','a:hvg','完整 HVG100 逐基因 z-score。A、B 按原 Ward 行聚类顺序分别显示前、后 50 个基因；两块共享同一细胞列顺序和色标。列顺序来自图示 HVG100 的 Ward 聚类；完整列树见顶部。图中 56 个基因带有核糖体相关分组注释。',.83)
page('主成分与谱系差异的补充视图')
fig('02_F3_PCA','a:pca','完整 PCA：A 为 PC1/PC2，B 为 PC1/PC3，C 为解释方差。前两轴累计 51.54\%，前三轴 63.67\%。所有点为细胞类型平均谱。',.33)
fig('02_F4_HSC_differences','a:markers','相对 HSC 的差异。A 为四目标差值 RMS；B 为预设 29 个谱系标记在逐基因标准化后与 HSC 的差值，保留全部可用标记。',.44)
page('聚类距离与 mDC 的影响')
fig('02_F5_clustering_comparison','a:cluster','HVG1000 下欧氏距离/Ward 与 Pearson 相关距离/平均连接的对照。HSC 在两者中均最接近造血前体。',.37)
fig('02_S1_mDC_sensitivity','a:mdc','mDC 敏感性：A 将 mDC 投影到固定主要 PCA 空间；B 将全部 22 类型重新拟合。mDC 的高零值率及其对全局坐标的影响单独显示。',.33)
page('全部细胞的覆盖与聚类稳定性')
q=rd('01_cell_quality');table('细胞覆盖。零值率以去重后的 19,994 个基因为分母。',['细胞类型','零值率（\%）','均值','最大值'],[[esc(S[r.cell_type]),val(100*r.zero_fraction,1),val(r['mean']),val(r.maximum)] for _,r in q.iterrows()],align='lrrr')
fig('02_S2_clustering_sensitivity','a:clustsens','聚类敏感性。改变 HVG 数目、表达尺度和距离定义，在聚类数 3–8 下以 ARI 对照主聚类；ARI 衡量分组一致性。',.28)
# B
section('B','统计模型与拟合诊断')
add('本附录的 r 和 $\\rho$ 分别为 Pearson 与 Spearman 相关。最低误差模型用于预测比较，1-SE 用于简约候选解释；模型规模均来自全数据最终拟合。')
m=rd('03_model_metrics');sel=rd('03_selected_models');se=sel[sel.scope=='full'].set_index(['target','model']);order=m[m.model=='Lasso_min'].sort_values('R2',ascending=False).target.tolist();rows=[]
for t in order:
 a=m[(m.target==t)&(m.model=='Lasso_min')].iloc[0];b=m[(m.target==t)&(m.model=='Lasso_1se')].iloc[0];c=m[(m.target==t)&(m.model=='Mean_baseline')].iloc[0]
 rows.append([esc(S[t]),val(a.R2),val(b.R2),val(a.RMSE,4),val(a.Pearson),val(a.Spearman),str(int(se.loc[(t,'Lasso_min'),'n_active_gt_1e_6']))+'/'+str(int(se.loc[(t,'Lasso_1se'),'n_active_gt_1e_6']))])
table('全部目标 Lasso 的折外性能与模型规模。min/1SE 列保留两套选择。',['目标','$R^2_{min}$','$R^2_{1SE}$','RMSE$_{min}$','$r_{min}$',r'$\rho_{min}$','活跃数 min/1SE'],rows)
# remaining complete baseline/1se metrics compact table
rows=[]
for t in order:
 a=m[(m.target==t)&(m.model=='Lasso_1se')].iloc[0];c=m[(m.target==t)&(m.model=='Mean_baseline')].iloc[0];rows.append([esc(S[t]),val(a.RMSE,4),val(a.Pearson),val(a.Spearman),val(c.R2,4),val(c.RMSE,4)])
table('1-SE 与均值基线的补充指标。基线由各训练折均值构成。',['目标','RMSE$_{1SE}$','$r_{1SE}$',r'$\rho_{1SE}$','$R^2_{mean}$','RMSE$_{mean}$'],rows,size='footnotesize')
page('六类统计模型的完整比较')
mm=rd('04_model_metrics');rows=[]
for t in focus:
 for model in ['Lasso','ElasticNet','Ridge','NonnegativeLasso','PCR','PLS','Mean_baseline']:
  z=mm[(mm.target==t)&(mm.model==model)]
  for _,r in z.sort_values('selection',ascending=False).iterrows():rows.append([AB[t],{'NonnegativeLasso':'Nonneg. Lasso','ElasticNet':'Elastic Net','Mean_baseline':'Mean'}.get(model,model),esc(r.selection),val(r.R2),val(r.RMSE,4),val(r.Pearson),val(r.Spearman)])
table('四目标所有模型与选择规则。各行均采用相同外层基因划分。',['目标','模型','选择','$R^2$','RMSE','$r$',r'$\rho$'],rows,size='footnotesize')
page('全部重点目标的留出预测')
fig('03_F2_focus_OOF_predictions','b:oof','四目标最低误差 Lasso 折外预测。每个基因仅使用其所属测试折的预测，全部 11,922 个基因均进入图中。六边形颜色为基因数的对数色阶，各面板独立；虚线为一致线。',.78)
add('预测在低幅度区域集中，大幅变化呈现收缩。评价以整套折外预测合并计算，因而图示代表所有共同基因的重建情况。模型预测单位与目标相对 HSC 的输入尺度差值相同。')
page('调参曲线与单次拟合初筛')
fig('03_F3_tuning_curves','b:tuning','全数据三折调参曲线：均值 MSE 及一个标准误，蓝虚线与朱红点线分别为最低误差与 1-SE。该图解释最终候选模型选择，独立外层评价见表 B1。',.43)
i=rd('03_focus_TF_top10_1se');rows=[]
for k in range(10):
 row=[k+1]
 for t in focus:
  rr=i[i.target==t].iloc[k];row.append(esc(rr.perturbation)+f' ({rr.coefficient_standardized:+.3f})')
 rows.append(row)
table('单次全数据 1-SE 拟合按系数幅度的前十 TF。括号为标准化系数，排序与稳定性排序用途不同。',['名次','B','Plasma','Mono','Ery'],rows,align='rllll',size='footnotesize')
page('完整系数路径的代表变量')
fig('03_F4_coefficient_paths','b:paths','每目标全数据 1-SE 模型中幅度最大的十二个活跃扰动，显示完整正则化路径。横轴由强到弱惩罚，蓝线为最低误差，朱红线为 1-SE。正负号贯穿路径保留，变量标签为原始扰动标签。',.78)
add('全部 2,717 个扰动参与拟合。单扰动相关性按与目标差值的绝对 Pearson 相关排序，用于直观参照；联合回归系数还取决于其他扰动模式，故两种排序分别保留。求解过程记录了收敛及边界参数，冻结结果通过数值验收。')
# C
section('C','候选与稳健性证据')
add('B、Plasma、Mono、Ery 对应四目标。频率为选中比例；固定抽样为 100 次，重新调参为 20 次。所有展示条目的固定选择频率及活跃符号一致性均为 100\%。分位数描述子抽样系数分布。')
top=rd('07_core_top10_retuning_check');selfr=rd('04_self_readout_results').set_index(['target','perturbation']);rows=[]
for t in focus:
 for _,r in top[top.target==t].iterrows():
  sr=selfr.loc[(t,r.perturbation)];v=val(100*sr.relative_coefficient_change,2) if sr.status=='completed' else '---'
  rows.append([AB[t],int(r.TF_stability_rank),esc(r.perturbation),val(r.full_coefficient),val(r.median_coefficient),f'[{r.coefficient_q025:.3f}, {r.coefficient_q975:.3f}]',f'{100*r.retuned_selection_frequency:.0f}',v])
table('40 个展示候选的完整统计及自身读出删除结果。最后一列为相对系数变化百分比，--- 为没有匹配自身读出。',['目标','排名','TF','全数据','中位数','2.5\%–97.5\%','重选\%','删除变化\%'],rows,align='lrlrrlrr',size='footnotesize')
page('全部稳定身份与交集')
s=rd('04_stable_priority_TFs')
for t in focus:
 ids=sorted(s[s.target==t].approved_symbol.dropna().unique());add(r'\noindent\textbf{'+esc(S[t])+f'（{len(ids)} 个规范身份）：'+r'} '+esc('、'.join(ids))+r'。\par\medskip')
o=rd('04_stable_TF_overlap');rows=[]
for j,t in enumerate(focus):
 for u in focus[j+1:]:
  r=o[(o.left_target==t)&(o.right_target==u)].iloc[0];rows.append([AB[t]+' / '+AB[u],int(r.intersection),val(r.Jaccard),esc(r.shared_TFs.replace(';',', '))])
table('成对交集。构建体按规范身份合并后计算。',['目标对','交集数','Jaccard','共享身份'],rows,align='lrrp{10cm}',size='footnotesize')
page('统计候选的文献证据与研究定位')
e=rd('05_TF_literature_evidence');rows=[]
for _,r in e.iterrows():
 refs='、'.join(r.source_ids.split(';')) if isinstance(r.source_ids,str) else '定向检索'
 rows.append([esc(S.get(r.target,r.target))+' / '+esc(r.TF),str(r.TF_stability_rank),esc(r.evidence_scope),esc(r.literature_support_zh)+r' '+esc(r.interpretation_and_validation_zh)+' '+esc(refs)])
table('统计证据—文献背景—验证方向。R 编号对应参考文献序号。',['目标 / TF','排名','证据背景','文献依据与后续问题'],rows,align='p{2.6cm}rp{2.6cm}p{9cm}',size='footnotesize')
add('文献核对以原始研究的摘要、可访问全文及官方索引为依据。部分全文访问受限，表中结论限定在已核对内容。BCL6 因机制证据列入讨论，未被追加到预设的前十自身读出删除实验；EBF1 作为已知 B 细胞参照保留其零系数结果。')
ed=rd('04_correlation_edges');table('相关扰动边（主阈值 $|r|\geq0.8$）。用于观察可替代模式。',['扰动 1','扰动 2','Pearson'],[[esc(r.perturbation_1),esc(r.perturbation_2),val(r.Pearson)] for _,r in ed.iterrows()],align='llr',size='footnotesize')
page('Lasso 与 Elastic Net 的相关变量结构')
fig('04_F5_correlated_predictor_groups','c:correlation','Lasso 与 Elastic Net 的全扰动系数对照及相关组。各点为扰动，重点前十 TF 高亮；连通组按相关阈值定义，保留组内不同成员。两种模型的活跃方向一致性需与候选组成和预测性能共同阅读。',.73)
a=rd('04_candidate_set_agreement');a=a[a.comparison=='raw_ElasticNet'];table('全活跃扰动集合的一致性。',['目标','Lasso 数','EN 数','交集','Jaccard','系数 Spearman'],[[AB[r.target],r.n_active_reference,r.n_active_alternative,r.n_intersection,val(r.active_Jaccard),val(r.coefficient_Spearman)] for _,r in a.iterrows()],size='footnotesize')
page('完整尺度与读出敏感性汇总')
add('原尺度差值、逐基因标准化差值和 landmark 子集具有各自响应与读出定义。表中保留各任务的原生评价尺度；同一 landmark 测试集上的训练覆盖比较在正文图 4D。等数量非 landmark 对照检验读出组成的影响。')
se=rd('04_sensitivity_metrics');base=rd('04_model_metrics');ref=rd('04_matched_landmark_reference_metrics');sz=rd('04_size_control_metrics');rows=[]
for label,d in [('Raw',base[(base.model=='Lasso')]),('Gene-scaled',se[se.dataset=='gene_scaled']),('Landmark',se[se.dataset=='landmark']),('All → landmark',ref),('Size match 1',sz[sz.dataset=='size_matched_01']),('Size match 2',sz[sz.dataset=='size_matched_02']),('Size match 3',sz[sz.dataset=='size_matched_03'])]:
 for t in focus:
  v=d[(d.target==t)&(d.model=='Lasso')];v1=v[v.selection=='min'];v2=v[v.selection=='1se']
  if v1.empty:continue
  r=v1.iloc[0];rr=v2.iloc[0];rows.append([label,AB[t],val(r.R2),val(rr.R2),val(r.RMSE,4),val(r.Pearson),val(r.Spearman)])
table('全部尺度与训练读出比较。MSE 选择相同，列出两套 $R^2$ 与最低误差模型其他指标。',['任务','目标','$R^2_{min}$','$R^2_{1SE}$','RMSE$_{min}$','$r_{min}$',r'$\rho_{min}$'],rows,align='llrrrrr',size='footnotesize')
ev=rd('04_candidate_evidence_matrix');rows=[]
for t in focus:
 d=ev[ev.target==t];rows.append([AB[t],f'{d.raw_ElasticNet_same_active_sign.sum()}/10',f'{d.gene_scaled_Lasso_same_active_sign.sum()}/10',f'{d.landmark_Lasso_same_active_sign.sum()}/10',f'{(d.self_readout_status=="completed").sum()}/10'])
table('前十候选的活跃方向保留与自身读出覆盖。',['目标','Elastic Net','Gene-scaled','Landmark','自身读出可匹配'],rows,align='lllll',size='footnotesize')
# D
section('D','外部扰动实验')
add('本附录覆盖划分、训练、完整汇总与代表条件。统一输出 18,977 个共同基因，HVG1000 每次由训练集确定；单个条件等权，三个种子再等权。模型缩写 scF+GEARS 指冻结 scFoundation 编码器的任务模型。')
table('条件划分与最终训练覆盖。各测试集另含 25 个未见单扰动。',['种子','训练','验证','测试','Seen0','Seen1','Seen2','训练细胞'],[[1,98,18,109,18,56,10,44158],[2,104,28,93,15,37,16,45661],[3,104,24,97,10,48,14,45696]])
run=rd('07_run_summary');table('九项正式运行。轮数为实际记录；所有最佳检查点在前八轮。耗时含训练与预测，显存为最终进程的分配峰值。',['模型','种子','实际轮数','最佳轮','验证 MSE','分钟','显存 GiB'],[[names[r.model],r.seed,r.epochs,r.best_epoch,val(r.best_validation_MSE,5),val(r.total_seconds/60,1),val(r.peak_allocated_GiB,1)] for _,r in run.iterrows()],align='lrrrrrr',size='footnotesize')
fig('07_F3_training_curves','d:training','所有模型与种子的验证曲线。点标最佳检查点，灰区标示预算变更前已完成的第八轮之后历史；这些额外轮次没有进入最终检查点选择。',.43)
page('总体指标与直接靶基因排除')
sm=rd('07_model_summary');rows=[]
for ex in [False,True]:
 for gs in ['shared','HVG1000']:
  for mod in mods[:-1]:
   r=sm[(sm.model==mod)&(sm.gene_set==gs)&(sm.exclude_direct_targets==ex)&(sm.group=='all')].iloc[0];rows.append(['排除' if ex else '包含','Shared' if gs=='shared' else 'HVG',names[mod],val(r.MSE_mean,5),val(r.MSE_sd,5),val(r.RMSE_mean,4),val(r.Pearson_delta_mean),val(r.Spearman_delta_mean)])
table('全部方法的总体指标。SD 为三个划分间的标准差；相关无定义记为 ---。',['直接靶基因','读出','模型','MSE','SD','RMSE','$r$',r'$\rho$'],rows,align='lllrrrrr',size='footnotesize')
add('三次划分不是独立生物重复，SD 描述划分差异。相关系数按有定义的条件汇总。对零变化或全未见组成下的恒定预测，相关系数无定义，但 MSE／RMSE 仍可计算。')
page('泛化分层的完整误差与模式相关性')
rows=[]
for gs in ['shared','HVG1000']:
 for gr in ['unseen_single','seen0','seen1','seen2']:
  for mod in mods[:-1]:
   d=sm[(sm.model==mod)&(sm.gene_set==gs)&(~sm.exclude_direct_targets)&(sm.group==gr)]
   if d.empty:continue
   r=d.iloc[0];rows.append(['Shared' if gs=='shared' else 'HVG',gr.replace('unseen_single','Single'),names[mod],val(r.MSE_mean,5),val(r.RMSE_mean,4),val(r.Pearson_delta_mean),val(r.Spearman_delta_mean)])
table('分层评价保留全部方法。Seen0/1/2 按最终训练条件重新核对。',['读出','条件组','模型','MSE','RMSE','$r$',r'$\rho$'],rows,align='lllrrrr',size='footnotesize')
page('直接加和的匹配比较与权重来源')
fig('07_F4_direct_additive_subset','d:direct','直接单扰动加和的可用子集。三个种子分别有 10、16、14 个可匹配组合；全部七种方法在相同条件子集内评价。比值相对直接加和计算，点与误差线保留划分间波动。',.43)
rows=[]
for mod in mods:
 for gs in ['shared','HVG1000']:
  d=sm[(sm.model==mod)&(sm.gene_set==gs)&(~sm.exclude_direct_targets)&(sm.group=='direct_additive_eligible')]
  if d.empty:continue
  r=d.iloc[0];rows.append([names[mod],'Shared' if gs=='shared' else 'HVG',val(r.MSE_mean,5),val(r.RMSE_mean,4),val(r.Pearson_delta_mean),val(r.Spearman_delta_mean)])
table('匹配子集完整指标。',['模型','读出','MSE','RMSE','$r$',r'$\rho$'],rows,align='llrrrr',size='footnotesize')
add(r'{\footnotesize 源码版本：GEARS \texttt{f374e43e197b}；scGPT \texttt{cebd6fae655b}；scFoundation \texttt{397631c495ed}。scGPT 使用官方 blood 通用权重，scFoundation 使用公共镜像固定版本 \texttt{cb434153a1ac}，文件 SHA-256 与该镜像发布标识一致。运行环境为 Python 3.12.13、PyTorch 2.12.0+cu130。预训练数据与 Norman 的样本重叠状态未核实。}\par')
page('不同响应幅度的完整代表预测')
fig('07_F5_representative_predictions','d:cases','种子 1 测试组合中，按真实 HVG 响应均方根选取最大、中位和最小响应条件。每行同一条件，保留 Ridge、GEARS、scGPT 与 scFoundation＋GEARS 全部预测；同一行采用共同坐标范围。该选择覆盖不同信号强度，不依据模型误差筛选。',.84)
page('PRDM1 的完整响应及指标')
fig('Supp_PRDM1','d:prdm','PRDM1 衔接案例的完整方法预测。上排单扰动，下排 CBFA2T3＋PRDM1 组合；每种方法均显示同一套训练 HVG1000。各排共享坐标范围，虚线为一致线。图中使用连续 A–L 标识。',.42)
br=rd('08_PRDM1_response_bridge');rows=[]
for cond in ['PRDM1+ctrl','CBFA2T3+PRDM1']:
 for mod in mods[:-1]:
  r=br[(br.condition==cond)&(br.model==mod)].iloc[0];rows.append(['单扰动' if cond=='PRDM1+ctrl' else '组合',names[mod],val(r.MSE,5),val(r.Pearson_delta),val(r.Spearman_delta),f'{r.improvement_vs_zero*100:.1f}'])
table('PRDM1 及其组合的所有方法。种子 1；改进为相对零变化 MSE 降低百分比。',['条件','模型','MSE','$r$',r'$\rho$','改进\%'],rows,align='llrrrr',size='footnotesize')
add('核心稳定 TF 与 Norman 靶标交集为 EGR1、OSR2、PRDM1、HNF4A、ZBTB25，共八个目标细胞—候选条目。PRDM1 单扰动在种子 1 为未见单扰动，组合为 seen0；在种子 2、3 中，PRDM1 单扰动处于训练集。直接单扰动加和不覆盖这两个种子 1 条件。')
# Coverage inline compact end; preserves self-contained provenance.
page('原正式图的内容覆盖')
coverage=[('02_F1','表达分布与零值','A1'),('02_F2','完整 HVG100','A2'),('02_F3','PCA 全部视图','图1B；A3'),('02_F4','HSC 差异及29标记','图1C；A4'),('02_F5','距离与聚类对照','图1D；A5'),('02_S1','mDC 投影和重拟合','A6'),('02_S2','聚类敏感性','A7'),('03_F1','21目标表现及模型规模','图2A、B；表B1、B2'),('03_F2','四目标完整预测','图2D；B1'),('03_F3','全数据调参曲线','B2'),('03_F4','系数路径','B3'),('03_F5','单次拟合TF初筛','表B4'),('04_F1','六模型两种选择','图2C；表B3'),('04_F2','40候选抽样分布','图3A；表C1'),('04_F3','跨目标候选热图','图3B、C；表C2'),('04_F4','尺度、landmark和方向','图4C、D；表C6、C7'),('04_F5','相关扰动组和EN比较','C1；表C4、C5'),('04_F6','自身读出删除','图4B；表C1'),('04_F7','等量非landmark对照','表C6'),('07_F1','共同/HVG全方法比较','图5A、B；表D3'),('07_F2','泛化分层','图5C、D；表D4'),('07_F3','九项训练曲线','D1；表D2'),('07_F4','直接加和匹配比较','D2；表D5'),('07_F5','响应强中弱代表条件','D3')]
table('24 组原正式图的独立内容去向。A1 等在此为附录图号，分图内部均使用连续字母。',['原图组','内容','新报告位置'],[[esc(a),b,c] for a,b,c in coverage],align='lp{6cm}p{7cm}',size='footnotesize')
add('相同信息通过共享坐标、对齐视图或数值表合并。正文图 6 及附录 D4 另外补全了 PRDM1 的全部方法响应。完整图表、方法、文献与执行记录可分别从本 PDF 和配套 MD 工作记录独立阅读。')
(O/'appendices.tex').write_text('\n'.join(parts).split(r'\clearpage\subsection{原正式图的内容覆盖}')[0])
pd.DataFrame(coverage,columns=['original_group','content','new_location']).to_csv(O/'figure_coverage.csv',index=False)
print('Appendices generated',len(parts))
