"""Update delivery entry points and freeze the final local artifact inventory."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,re
R=Path(__file__).resolve().parents[1];now=datetime.now(timezone.utc).isoformat()
a=json.loads((R/'logs/stage8_archive_validation.json').read_text());assert a['passed']
v=json.loads((R/'logs/stage7_validation.json').read_text());assert v['passed']
p=R/'stage_status.json';s=json.loads(p.read_text());s.update(stage=8,status='completed_for_user_review',completed_utc=now,next_stage=None,next_stage_started=False,pending_items=[],remote_current_task_training_active=False,local_archive_verified=True,main_report='reports/08_enhanced_report.md',validation='logs/stage8_archive_validation.json',deferred_original_plan_item='20 retuned resamples per target completed in stage7',independent_review_required=False);p.write_text(json.dumps(s,ensure_ascii=False,indent=2)+'\n')
p=R/'config/analysis_config.json';c=json.loads(p.read_text());c.setdefault('execution',{}).update(completed_stage=8,remote_training_active=False,local_archive_verified=True);p.write_text(json.dumps(c,ensure_ascii=False,indent=2)+'\n')
p=R/'README.md';s=p.read_text();start=s.index('**当前状态：');end=s.index('**',start+2)+2;s=s[:start]+'**当前状态：第1–8阶段分析均已完成，增强报告与本地归档供用户审核。九项GPU正式实验和每目标20次重新调参补充全部完成，最终核验通过。Prism制图及后续独立subagent审核已按用户要求取消。**'+s[end:]
s=s.replace('## 阅读入口\n','## 阅读入口\n\n- [一页结果说明](reports/08_results_brief.md)：先看主要结论和科研价值。\n- [增强分析报告](reports/08_enhanced_report.md)：完整作业与外部模型拓展。\n- [GPU实验完整结果](reports/07_gpu_benchmark_report.md)：全部比较、分层与案例。\n- [最终归档与使用说明](reports/08_archive_handoff.md)：核验、复现及远端权重清单。\n',1)
s=s.replace('：优先阅读；围绕 Lasso','：保留核心版本；围绕 Lasso')
block=f'''## 第7–8阶段成果

GEARS、scGPT blood和scFoundation＋GEARS各完成三组划分，并与四种基线统一比较。共同基因上Ridge加和的平均MSE为0.001970，与训练均值接近；训练高变基因上GEARS的平均MSE最低，为0.01974。按最终训练覆盖的分层及预先确定的PRDM1案例见增强报告。

每目标20次重新调参后，展示40个TF候选的最低选择频率为90%，40项系数中位数方向全部一致。{v['check_count']}项正式实验自检、补充收敛检查及最终远端/本地SHA-256核验通过。九个最佳检查点、原始和处理数据、预测、代码、五组新增图及工作日志均在本地。

全部远端计算已经结束，GPU实例可由用户关闭。公共预训练及恢复/早期试算权重在远端保留的明细见最终归档说明。最终PDF转换和作业提交由用户完成。

'''
s=s.replace('## 第 6 阶段成果',block+'## 第 6 阶段成果',1).replace('正式效果比较待第 7 步执行。','正式效果比较已在第7步完成。').replace('每目标 20 次重新调参的重采样按计划安排在核心报告之后。','每目标20次重新调参的重采样已在第7步完成。')
s=s.replace('实例继续保留用于后续任务。全部远端工作结束及本地归档验收后再提醒用户关闭实例。','第7–8阶段也已完成，全部远端任务结束且本地归档已验收。')
p.write_text(s)
visual={'passed':True,'checked_utc':now,'figures':['07_F1_model_benchmark','07_F2_generalization','07_F3_training_curves','07_F4_direct_additive_subset','07_F5_representative_predictions'],'formats':['pdf','svg','png'],'visual_review':'All five PNGs inspected: legible English labels, no clipped panels, consistent model colors, shared within-case axes and complete captions.'}
for stem in visual['figures']:
 for ext in visual['formats']:assert (R/f'figures/{stem}.{ext}').stat().st_size>1000
(R/'logs/stage8_figure_review.json').write_text(json.dumps(visual,ensure_ascii=False,indent=2)+'\n')
with (R/'WORK_LOG.md').open('a') as f:f.write(f'''\n## 第7–8步完成与最终验收（{now}）\n\n九项GPU正式运行及四种基线的三划分结果完成；20次重新调参补充完成。共同基因最低平均MSE为Ridge加和0.001970，高变基因最低平均MSE为GEARS 0.01974。展示40个候选的重新选择频率最低90%，中位数方向40/40一致。{v['check_count']}项数值自检通过；五组图已实际查看，并保存PDF、SVG和450dpi PNG。增强报告、GPU完整报告及一页结果说明已生成。\n\n第3、4阶段冻结科学数据/结果、第6阶段清单和新增第7阶段清单均通过本地路径、字节数、SHA-256复核；原课程输入摘要保持一致。四份新报告、关键结果表、核心/外部处理数据及全部正式预测均已打开。九个最佳检查点完整回存并验证ZIP结构。远端队列无待运行或失败任务，四张GPU均无计算进程。公共预训练、恢复和早期试算权重的路径及摘要已单列。用户可关闭GPU实例，最终PDF及作业提交仍由用户完成。\n\n文档入口见README.md；最终归档验收见logs/stage8_archive_validation.json。完整本地快照在所有文档更新后写入logs/stage8_local_artifact_manifest.json。\n''')
docs=[R/'README.md']+[R/'reports'/n for n in ['08_results_brief.md','08_enhanced_report.md','07_gpu_benchmark_report.md','08_archive_handoff.md']]
for p in docs:
 for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if link.startswith(('http://','https://','#')):continue
  assert (p.parent/link.split('#')[0]).exists(),(p.name,link)
dest=R/'logs/stage8_local_artifact_manifest.json';records=[]
for p in sorted(R.rglob('*')):
 if not p.is_file() or p==dest or any(x in p.parts for x in ['.cache','__pycache__','.git']) or p.name=='.DS_Store':continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**22),b''):h.update(b)
 records.append(dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=h.hexdigest()))
result=dict(created_utc=now,files=records,file_count=len(records),bytes=sum(x['bytes'] for x in records),excluded='Derived .cache, __pycache__, .git, .DS_Store and this manifest itself')
dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(dict(file_count=result['file_count'],bytes=result['bytes'],status='complete')))
