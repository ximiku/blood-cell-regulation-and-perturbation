"""Publication-style plots of the complete perturbation benchmark."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
R=Path(__file__).resolve().parents[1];T=R/'tables';F=R/'figures';F.mkdir(exist_ok=True)
ORDER=['zero_change','training_condition_mean','Ridge_additive','GEARS','scGPT_blood','scFoundation_GEARS']
LABEL={'zero_change':'No change','training_condition_mean':'Training mean','Ridge_additive':'Additive Ridge','direct_single_additive':'Direct additive','GEARS':'GEARS','scGPT_blood':'scGPT blood','scFoundation_GEARS':'scFoundation + GEARS'}
COLOR={'zero_change':'#9AA5B1','training_condition_mean':'#525D68','Ridge_additive':'#0072B2','direct_single_additive':'#56B4E9','GEARS':'#009E73','scGPT_blood':'#D55E00','scFoundation_GEARS':'#8E61A9'}
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.65,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
manifest=[]
def save(fig,stem,caption):
 for ext in ['png','pdf','svg']:fig.savefig(F/f'{stem}.{ext}',dpi=450,bbox_inches='tight',pad_inches=.08)
 plt.close(fig);manifest.append({'stem':stem,'caption_zh':caption,'formats':['png','pdf','svg']})
def panel(ax,letter):ax.text(-.12,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=12)
def dots(ax,data,metric,models,reference=None):
 for i,m in enumerate(models):
  s=data[data.model==m].sort_values('seed');v=s[metric].to_numpy();finite=v[np.isfinite(v)]
  if not len(finite):continue
  ax.scatter(v,i+np.linspace(-.12,.12,len(v)),s=17,color=COLOR[m],alpha=.65,zorder=3,edgecolors='white',linewidths=.3)
  ax.errorbar(np.nanmean(v),i,xerr=np.nanstd(v,ddof=1),fmt='D',ms=4.5,color=COLOR[m],elinewidth=1.3,capsize=2,zorder=4)
 ax.set_yticks(range(len(models)),[LABEL[m] for m in models]);ax.invert_yaxis();ax.grid(axis='x',lw=.5,color='#E7EBEF');ax.set_axisbelow(True)
 if reference is not None:ax.axvline(reference,color='#9BA4AD',ls='--',lw=.75)
def main():
 s=pd.read_csv(T/'07_seed_metrics.csv');d=pd.read_csv(T/'07_condition_metrics.csv');summary=pd.read_csv(T/'07_model_summary.csv')
 main=s[(s.group=='all')&(~s.exclude_direct_targets)&s.model.isin(ORDER)].copy();main['MSE_ratio_to_zero']=1-main.improvement_vs_zero
 fig,axes=plt.subplots(2,2,figsize=(9,6.7))
 for row,gs in enumerate(['shared','HVG1000']):
  sub=main[main.gene_set==gs];dots(axes[row,0],sub,'MSE_ratio_to_zero',ORDER,1);dots(axes[row,1],sub,'Pearson_delta',ORDER)
  axes[row,0].set_xlabel('Condition-equal MSE / no-change MSE');axes[row,1].set_xlabel('Pearson correlation of expression changes')
  for col in [0,1]:axes[row,col].set_title('18,977 shared genes' if gs=='shared' else '1,000 training-variable genes',loc='left',pad=12);panel(axes[row,col],'ABCD'[2*row+col])
  axes[row,1].set_xlim(min(-.05,float(sub.Pearson_delta.min())-.05),1)
 fig.subplots_adjust(left=.20,right=.98,bottom=.09,top=.94,wspace=.95,hspace=.55)
 save(fig,'07_F1_model_benchmark','完整测试集上的模型比较。小圆点为三个划分，菱形及横线为三个划分的均值±标准差。A、C为条件等权MSE相对零变化基线的比值，越小越好；B、D为相对对照的表达变化Pearson相关，越大越好。零变化预测的相关系数无定义，保留为空。每个划分按条件等权，再跨划分汇总；训练高变基因由相应训练集确定。')

 groups=['seen0','seen1','seen2','unseen_single'];models=ORDER[2:];fig,axes=plt.subplots(1,2,figsize=(9,3.7));arrays=[]
 for gs in ['shared','HVG1000']:
  arr=np.array([[summary[(summary.model==m)&(summary.gene_set==gs)&(~summary.exclude_direct_targets)&(summary.group==g)].improvement_vs_mean_mean.iloc[0]*100 for g in groups] for m in models]);arrays.append(arr)
 vmax=max(10,float(np.nanmax(np.abs(arrays))));norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
 for ax,arr,title,letter in zip(axes,arrays,['18,977 shared genes','1,000 training-variable genes'],'AB'):
  im=ax.imshow(arr,cmap='RdBu',norm=norm,aspect='auto');ax.set_xticks(range(4),['Seen 0','Seen 1','Seen 2','Unseen\nsingle']);ax.set_yticks(range(4),[LABEL[m] for m in models]);ax.set_title(title,loc='left',pad=13);panel(ax,letter)
  for i in range(4):
   for j in range(4):ax.text(j,i,f'{arr[i,j]:+.1f}%',ha='center',va='center',color='white' if abs(arr[i,j])>.62*vmax else '#17212B',fontsize=8)
  ax.set_xticks(np.arange(-.5,4,1),minor=True);ax.set_yticks(np.arange(-.5,4,1),minor=True);ax.grid(which='minor',color='white',lw=1);ax.tick_params(which='minor',bottom=False,left=False)
 fig.subplots_adjust(left=.20,right=.87,bottom=.18,top=.87,wspace=.85);cax=fig.add_axes([.91,.22,.018,.59]);fig.colorbar(im,cax=cax,label='MSE improvement over training mean (%)')
 save(fig,'07_F2_generalization','按最终训练条件中的基因覆盖分层。数值为各划分条件等权MSE相对训练条件均值基线的改善百分比，再对三个划分求均值；正值表示误差降低。Seen0/1/2表示组合中有0/1/2个基因出现在最终训练条件中，另列未见单扰动。')

 fig,axes=plt.subplots(1,3,figsize=(10,3.3));histrows=[]
 for ax,m,letter in zip(axes,ORDER[3:],'ABC'):
  for seed,style in zip([1,2,3],['-','--',':']):
   h=pd.read_csv(R/f'results/stage7/{m}/seed{seed}/training_history.csv');h['model']=m;h['seed']=seed;histrows.append(h);ax.plot(h.epoch,h.validation_condition_equal_MSE,color=COLOR[m],ls=style,lw=1.4,label=f'Seed {seed}');j=h.validation_condition_equal_MSE.idxmin();ax.scatter(h.loc[j,'epoch'],h.loc[j,'validation_condition_equal_MSE'],color=COLOR[m],s=24,edgecolors='white',lw=.4,zorder=3)
  ax.axvspan(8.5,max(8.5,max(pd.read_csv(R/f'results/stage7/{m}/seed{k}/training_history.csv').epoch.max() for k in [1,2,3])+.2),color='#EEF1F4',zorder=0);ax.set_title(LABEL[m],loc='left',pad=13);panel(ax,letter);ax.set_xlabel('Epoch');ax.set_ylabel('Validation MSE');ax.set_yscale('log');ax.grid(color='#E7EBEF',lw=.5);ax.legend(frameon=False)
 fig.subplots_adjust(left=.08,right=.99,bottom=.19,top=.84,wspace=.45);pd.concat(histrows).to_csv(T/'07_training_curves.csv',index=False)
 save(fig,'07_F3_training_curves','三类GPU模型的完整验证轨迹，按扰动条件等权计算MSE；每种线型表示一个划分，圆点为最低验证误差轮次。纵轴为对数刻度。连续5轮无改善时早停；用户将预算调整为8轮，灰区保留调整前已经完成的额外轮次。所有最终检查点均来自前8轮。')

 sub=s[(s.group=='direct_additive_eligible')&(s.gene_set=='shared')&(~s.exclude_direct_targets)].copy();ref=sub[sub.model=='direct_single_additive'].set_index('seed').MSE;sub['MSE_ratio_to_additive']=[v/ref.loc[k] for v,k in zip(sub.MSE,sub.seed)]
 fig,ax=plt.subplots(figsize=(6.4,3.8));dots(ax,sub,'MSE_ratio_to_additive',ORDER[:3]+['direct_single_additive']+ORDER[3:],1);ax.set_xlabel('Condition-equal MSE / direct-additive MSE');ax.set_title('Matched conditions with both training singles available',loc='left',pad=13);fig.subplots_adjust(left=.30,right=.98,bottom=.18,top=.87)
 save(fig,'07_F4_direct_additive_subset','直接加和基线可用条件上的匹配比较。每个划分均要求组合的两个单扰动均在训练集中，所有模型在该相同子集上重新评价。小点为三个划分，菱形与横线为均值±标准差。适用条件数保存在覆盖表中，与全测试集排序分别呈现。')

 # Choose representative conditions by observed response amplitude before inspecting model errors.
 seed=1;base=np.load(R/'results/stage7/Ridge_additive/seed1/test_predictions.npz');ix=np.load(R/'data/processed/stage6/evaluation_indices.npz');mask=np.isin(base['gene_indices'],ix['hvg1000_seed1']);names=base['conditions'].astype(str);delta=base['observed']-base['control_mean'];combos=np.flatnonzero(np.array([len(c.split('+'))==2 and 'ctrl' not in c for c in names]));rank=combos[np.argsort(np.mean(delta[combos][:,mask]**2,axis=1))];chosen=[rank[-1],rank[len(rank)//2],rank[0]];case_models=['Ridge_additive',*ORDER[3:]];fig,axes=plt.subplots(3,4,figsize=(10,7.4));caserows=[]
 case_predictions={m:np.load(R/f'results/stage7/{m}/seed1/test_predictions.npz')['predicted']-base['control_mean'] for m in case_models}
 row_limits={j:max(float(np.max(np.abs(delta[j,mask]))),max(float(np.max(np.abs(case_predictions[m][j,mask]))) for m in case_models))*1.07 for j in chosen}
 for col,m in enumerate(case_models):
  z=np.load(R/f'results/stage7/{m}/seed1/test_predictions.npz');pr=z['predicted']-z['control_mean'];assert np.array_equal(z['conditions'],names)
  for row,j in enumerate(chosen):
   ax=axes[row,col];x=delta[j,mask];y=pr[j,mask];limit=row_limits[j];ax.scatter(x,y,s=5,color=COLOR[m],alpha=.40,edgecolors='none',rasterized=True);ax.plot([-limit,limit],[-limit,limit],ls='--',color='#87929C',lw=.65);ax.axhline(0,lw=.4,color='#CCD2D8');ax.axvline(0,lw=.4,color='#CCD2D8');ax.set_xlim(-limit,limit);ax.set_ylim(-limit,limit);pear=float(np.corrcoef(x,y)[0,1]);mse=float(np.mean((x-y)**2));ax.text(.04,.95,f'r = {pear:.2f}\nMSE = {mse:.3g}',transform=ax.transAxes,va='top',fontsize=7)
   if row==0:ax.set_title(LABEL[m],fontsize=8,pad=10)
   if col==0:ax.set_ylabel(names[j]+'\nPredicted change',fontsize=7)
   if row==2:ax.set_xlabel('Observed change')
   caserows.append({'condition':names[j],'selection':['largest','median','smallest'][row]+'_observed_HVG_response','model':m,'seed':1,'Pearson_delta':pear,'MSE':mse})
 fig.subplots_adjust(left=.12,right=.99,bottom=.09,top=.94,wspace=.35,hspace=.35);pd.DataFrame(caserows).to_csv(T/'07_representative_cases.csv',index=False)
 save(fig,'07_F5_representative_predictions','种子1测试组合中，按真实高变基因响应均方幅度选取最大、中位、最小三个案例，选择过程不使用模型误差。每点为一个训练高变基因；横轴为真实表达变化，纵轴为预测变化，虚线为一致线。模型与案例的完整数值保存于源表。')

 (R/'figures/07_figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
 (R/'reports/07_figure_captions.md').write_text('# GPU 拓展图注\n\n'+'\n\n'.join(f"## {x['stem']}\n\n{x['caption_zh']}" for x in manifest)+'\n')
 print('Exported',len(manifest),'figure groups in PDF, SVG and 450 dpi PNG.')
if __name__=='__main__':main()
