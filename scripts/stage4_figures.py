"""Stage-4 comparisons and candidate robustness in a consistent publication style."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm

ROOT=Path(__file__).resolve().parents[1];TAB=ROOT/'tables';FIG=ROOT/'figures'
FOCUS=json.loads((ROOT/'config/stage4_config.json').read_text())['focus_targets']
D=pd.read_csv(TAB/'02_cell_dictionary.csv').set_index('cell_type');SHORT=D.short_label.to_dict();COLORS=D.color.to_dict()
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.65,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,
 'savefig.facecolor':'white','xtick.major.width':.6,'ytick.major.width':.6})
MANIFEST=[]
def panel(ax,letter):ax.annotate(letter,(0,1),xycoords='axes fraction',xytext=(-25,11),textcoords='offset points',fontsize=12,fontweight='bold')
def save(fig,name,caption):
    for ext in ['pdf','svg','png']:fig.savefig(FIG/f'{name}.{ext}',dpi=450,bbox_inches='tight',pad_inches=.08)
    plt.close(fig);MANIFEST.append({'stem':name,'caption_zh':caption,'formats':['pdf','svg','png']})

def main():
    metric=pd.read_csv(TAB/'04_model_metrics.csv');top=pd.read_csv(TAB/'04_top10_TF_candidates.csv')
    order=['Mean_baseline','PCR','PLS','Ridge','NonnegativeLasso','Lasso','ElasticNet']
    names={'Mean_baseline':'Mean baseline','NonnegativeLasso':'Nonnegative Lasso','ElasticNet':'Elastic Net'}
    fig,axes=plt.subplots(2,2,figsize=(7.6,6.1))
    for ax,t,letter in zip(axes.ravel(),FOCUS,'ABCD'):
        f=metric[metric.target==t]
        for k,model in enumerate(order):
            m=f[f.model==model]
            if model=='Mean_baseline':ax.scatter(m.R2,k,marker='|',s=55,color='#67727D')
            else:
                vals=m.set_index('selection').R2
                ax.plot([vals['min'],vals['1se']],[k,k],color='#B7BEC5',lw=1)
                ax.scatter(vals['min'],k-.10,s=23,color='#0072B2',zorder=3)
                ax.scatter(vals['1se'],k+.10,s=21,color='#D55E00',marker='D',zorder=3)
        ax.set_yticks(range(len(order)),[names.get(m,m) for m in order]);ax.invert_yaxis();ax.axvline(0,color='#C1C8CE',lw=.7)
        ax.set_xlabel('Out-of-fold R²');ax.set_title(SHORT[t],color=COLORS[t],loc='left',pad=14);panel(ax,letter)
        ax.grid(axis='x',color='#ECEFF2',lw=.6);ax.set_axisbelow(True)
    handles=[Line2D([],[],marker='o',ls='',color='#0072B2',label='Minimum validation MSE'),Line2D([],[],marker='D',ls='',color='#D55E00',label='1-SE selection')]
    fig.legend(handles=handles,loc='lower center',ncol=2,frameon=False)
    fig.subplots_adjust(left=.17,right=.98,bottom=.14,top=.92,wspace=.70,hspace=.48)
    save(fig,'04_F1_model_comparison','四目标共享五折基因划分的模型比较。圆点为最低内层 MSE 选择，菱形为 1-SE 选择；灰色短线为训练均值基线。指标由全部折外预测合并计算，连线连接同一模型的两种选择。六种统计模型全部保留。')

    fig,axes=plt.subplots(2,2,figsize=(7.8,7.0))
    for ax,t,letter in zip(axes.ravel(),FOCUS,'ABCD'):
        f=top[top.target==t].copy();yy=np.arange(len(f));median=f.median_coefficient.to_numpy()
        for i,(_,row) in enumerate(f.iterrows()):
            color=COLORS[t] if row.stable_priority else '#9AA4AD'
            ax.plot([row.coefficient_q025,row.coefficient_q975],[i,i],color=color,lw=1.6,alpha=.65)
            ax.scatter(median[i],i,color=color,s=23,zorder=3,marker='o' if row.stable_priority else 'x')
            ax.text(1.025,i,f'{row.selection_frequency:.0%}',transform=ax.get_yaxis_transform(),va='center',fontsize=6.7)
        ax.axvline(0,color='#A5AFB8',lw=.8,ls='--');ax.set_yticks(yy,f.perturbation);ax.invert_yaxis()
        ax.set_xlabel('Standardized coefficient across subsamples');ax.set_title(SHORT[t],color=COLORS[t],loc='left',pad=14)
        ax.text(1.025,1.015,'Selected',transform=ax.transAxes,fontsize=6.7);panel(ax,letter)
    fig.subplots_adjust(left=.12,right=.91,bottom=.09,top=.93,wspace=.74,hspace=.47)
    save(fig,'04_F2_stability_candidates','四目标按选择频率、系数绝对值中位数和全数据系数幅度排序的前 10 个 TF。每次无放回抽取 80% 基因，共 100 次；固定全数据 1-SE 相对惩罚，重新计算训练尺度和 λmax。点为系数中位数，线段为子抽样系数的 2.5%–97.5% 分位数，右侧为选择频率。彩色圆点满足稳定优先规则；灰色叉号为其余展示候选。分位线描述重采样分布。')

    stable=pd.read_csv(TAB/'04_stability_all_TFs.csv');labels=[]
    for t in FOCUS:
        for p in top[top.target==t].perturbation:
            if p not in labels:labels.append(p)
    matrix=np.zeros((len(labels),4));qualify=np.zeros_like(matrix,bool)
    for j,t in enumerate(FOCUS):
        f=stable[stable.target==t].set_index('perturbation').loc[labels]
        matrix[:,j]=f.selection_frequency.to_numpy()*np.sign(f.full_coefficient.to_numpy());qualify[:,j]=f.stable_priority
    fig,ax=plt.subplots(figsize=(5.9,max(6.6,.20*len(labels)+1.5)))
    im=ax.imshow(matrix,aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1)
    ax.set_xticks(range(4),[SHORT[t] for t in FOCUS]);ax.set_yticks(range(len(labels)),labels)
    for i,j in zip(*np.where(qualify)):ax.scatter(j,i,s=9,color='white',edgecolors='#333333',linewidths=.25)
    ax.set_xticks(np.arange(-.5,4,1),minor=True);ax.set_yticks(np.arange(-.5,len(labels),1),minor=True)
    ax.grid(which='minor',color='white',lw=.5);ax.tick_params(which='minor',left=False,bottom=False)
    cb=fig.colorbar(im,ax=ax,fraction=.05,pad=.06);cb.set_label('Selection frequency × full-model sign')
    ax.set_title('Shared and target-specific candidate patterns',loc='left',pad=14)
    fig.subplots_adjust(left=.20,right=.90,bottom=.08,top=.94)
    save(fig,'04_F3_cross_target_candidates','四目标前 10 个 TF 的并集。红蓝表示全数据系数的正负方向，颜色强度表示该目标中 100 次子抽样的选择频率；白点标记满足稳定优先规则。每个扰动保留原始构建体标签。完整稳定 TF 身份层面的交集与 Jaccard 表另存。')

    native=pd.read_csv(TAB/'04_sensitivity_metrics.csv');original=pd.read_csv(TAB/'04_sensitivity_original_units.csv')
    matched=pd.read_csv(TAB/'04_matched_landmark_reference_metrics.csv');sensitivity=pd.read_csv(TAB/'04_top10_candidate_sensitivity.csv')
    fig,axes=plt.subplots(1,3,figsize=(10.1,3.4))
    for k,t in enumerate(FOCUS):
        f=original[(original.target==t)&(original.selection=='min')].set_index('dataset')
        a=f.loc['raw','R2'];b=native[(native.target==t)&(native.dataset=='gene_scaled')&(native.model=='Lasso')&(native.selection=='min')].R2.iloc[0];axes[0].plot([a,b],[k,k],color=COLORS[t],lw=1.1)
        axes[0].scatter(a,k,color=COLORS[t],s=25);axes[0].scatter(b,k,color=COLORS[t],s=25,marker='s',facecolors='none')
        a=matched[(matched.target==t)&(matched.selection=='min')].R2.iloc[0]
        b=native[(native.target==t)&(native.dataset=='landmark')&(native.model=='Lasso')&(native.selection=='min')].R2.iloc[0]
        axes[1].plot([a,b],[k,k],color=COLORS[t],lw=1.1);axes[1].scatter(a,k,color=COLORS[t],s=25);axes[1].scatter(b,k,color=COLORS[t],s=25,marker='s',facecolors='none')
    for ax in axes[:2]:ax.set_yticks(range(4),[SHORT[t] for t in FOCUS]);ax.invert_yaxis();ax.axvline(0,color='#C3CAD0',lw=.6);ax.set_xlabel('Out-of-fold R² (input units)')
    axes[0].set_xlabel('Out-of-fold R² (each response definition)')
    comparisons=['raw_ElasticNet','landmark_Lasso','gene_scaled_Lasso'];markers=['o','s','^']
    for i,c in enumerate(comparisons):
        for j,t in enumerate(FOCUS):
            val=sensitivity[(sensitivity.comparison==c)&(sensitivity.target==t)].same_active_sign.mean()
            axes[2].scatter(val,j+(i-1)*.15,color=COLORS[t],marker=markers[i],s=25)
    axes[2].set_yticks(range(4),[SHORT[t] for t in FOCUS]);axes[2].invert_yaxis();axes[2].set_xlim(-.04,1.07);axes[2].set_xlabel('Top-10 candidates retaining active sign')
    for ax,title,letter in zip(axes,['Response definitions','Landmark test genes','Candidate agreement'],'ABC'):ax.set_title(title,loc='left',pad=13);panel(ax,letter)
    axes[0].legend(handles=[Line2D([],[],marker='o',ls='',color='#666666',label='Raw-response fit'),Line2D([],[],marker='s',ls='',color='#666666',mfc='none',label='Gene-scaled fit')],loc='upper left',bbox_to_anchor=(-.1,-.22),frameon=False,fontsize=6.5)
    axes[1].legend(handles=[Line2D([],[],marker='o',ls='',color='#666666',label='All-gene training'),Line2D([],[],marker='s',ls='',color='#666666',mfc='none',label='Landmark training')],loc='upper left',bbox_to_anchor=(-.1,-.22),frameon=False,fontsize=6.5)
    axes[2].legend(handles=[Line2D([],[],marker=m,ls='',color='#666666',label=l) for m,l in zip(markers,['Elastic Net','Landmark Lasso','Gene-scaled Lasso'])],loc='upper left',bbox_to_anchor=(-.1,-.22),frameon=False,fontsize=6.5)
    fig.subplots_adjust(left=.08,right=.98,top=.84,bottom=.30,wspace=.72)
    save(fig,'04_F4_sensitivity_summary','A：分别评价原差值与逐基因标准化差值两个响应定义，报告各自响应单位下的折外 R²；模型主排名使用原差值任务。B：在相同的 965 个 landmark 留出读出上，比较全基因训练与 landmark 专项训练。C：前 10 个 TF 在 Elastic Net、landmark Lasso、逐基因标准化 Lasso 的全数据 1-SE 模型中保持活跃且方向相同的比例。三种候选对照以不同形状区分。')

    group=pd.read_csv(TAB/'04_correlation_groups.csv');coefficient=pd.read_csv(TAB/'04_all_coefficients_annotated.csv',low_memory=False)
    grouped=set(group[np.isclose(group.threshold,.8)&(group.group_size>1)].perturbation)
    fig,axes=plt.subplots(2,2,figsize=(7.5,6.5))
    for ax,t,letter in zip(axes.ravel(),FOCUS,'ABCD'):
        f=coefficient[(coefficient.dataset=='raw')&(coefficient.target==t)&(coefficient.selection=='1se')&coefficient.model.isin(['Lasso','ElasticNet'])].pivot(index='perturbation',columns='model',values='coefficient_standardized')
        ax.scatter(f.Lasso,f.ElasticNet,s=5,color='#ACB6BE',alpha=.45,linewidths=0,rasterized=True)
        selected=f.loc[top[top.target==t].perturbation]
        ax.scatter(selected.Lasso,selected.ElasticNet,s=20,color=COLORS[t],zorder=3,edgecolors='white',linewidths=.3)
        related=f.loc[list(grouped)]
        ax.scatter(related.Lasso,related.ElasticNet,s=34,facecolors='none',edgecolors='#424D56',linewidths=.6,zorder=4)
        low=min(f.min().min(),-.02)*1.12;high=max(f.max().max(),.02)*1.12
        ax.plot([low,high],[low,high],ls='--',lw=.7,color='#8A959F');ax.set_xlim(low,high);ax.set_ylim(low,high)
        ax.axhline(0,color='#D2D7DC',lw=.5);ax.axvline(0,color='#D2D7DC',lw=.5)
        ax.set_xlabel('Lasso standardized coefficient');ax.set_ylabel('Elastic Net standardized coefficient')
        ax.set_title(SHORT[t],color=COLORS[t],loc='left',pad=14);panel(ax,letter)
        for label in [selected.Lasso.idxmin(),selected.Lasso.idxmax()]:
            ax.annotate(label,(selected.loc[label,'Lasso'],selected.loc[label,'ElasticNet']),xytext=(5,-10 if selected.loc[label,'Lasso']<0 else 6),textcoords='offset points',fontsize=6.5)
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color='#ACB6BE',label='All perturbations'),Line2D([],[],marker='o',ls='',color='#0072B2',label='Top-10 TFs'),Line2D([],[],marker='o',ls='',color='#424D56',mfc='none',label='Correlated-group members (|r| ≥ 0.80)')],loc='lower center',ncol=3,frameon=False,fontsize=6.5)
    fig.subplots_adjust(left=.10,right=.97,top=.92,bottom=.14,wspace=.43,hspace=.49)
    save(fig,'04_F5_correlated_predictor_groups','Lasso 与 Elastic Net 全数据 1-SE 标准化系数逐项对照。灰点为全部 2,717 个扰动，彩点为稳定性排序的前 10 个 TF；黑色空心圈标记绝对 Pearson 相关至少 0.80 的多成员连通组成员，虚线为两模型系数相同。各目标标注展示 TF 中最正和最负的两项。完整相关组成员、预测贡献、活跃数和 0.70/0.90 阈值对照保存在源表。')

    own=pd.read_csv(TAB/'04_self_readout_results.csv');fig,axes=plt.subplots(2,2,figsize=(7.7,6.9))
    for ax,t,letter in zip(axes.ravel(),FOCUS,'ABCD'):
        labels=top[top.target==t].perturbation.tolist();f=own[own.target==t].set_index('perturbation').loc[labels]
        for k,(label,row) in enumerate(f.iterrows()):
            if row.status=='completed':ax.scatter(100*row.relative_coefficient_change,k,s=24,color=COLORS[t],marker='o' if row.same_active_sign else 'x')
            else:ax.text(.02,k,'readout absent',transform=ax.get_yaxis_transform(),fontsize=6.5,va='center',color='#7D8790')
        ax.axvline(0,color='#A7B0B8',ls='--',lw=.7);ax.set_yticks(range(len(labels)),labels);ax.set_ylim(len(labels)-.5,-.5)
        ax.set_xlabel('Coefficient change / |baseline| (%)');ax.set_title(SHORT[t],color=COLORS[t],loc='left',pad=14);panel(ax,letter)
        lo,hi=ax.get_xlim();extent=max(abs(lo),abs(hi),1);ax.set_xlim(-extent*1.1,extent*1.1)
    fig.subplots_adjust(left=.12,right=.98,top=.93,bottom=.10,wspace=.55,hspace=.46)
    save(fig,'04_F6_self_readout_check','逐候选删除与其 HGNC 规范身份匹配的自身表达读出行，在固定全数据 Lasso 1-SE 相对惩罚下重新拟合。横轴为删除前后系数差除以基线绝对系数；基线和删除模型均采用 10⁻⁶ 的更严格对偶间隙阈值。圆点表示保留活跃方向，叉号表示方向或活跃状态改变；无自身读出的候选单列。完整 Pearson 变化、匹配行与证书另存。')
    control_path=TAB/'04_size_control_metrics.csv'
    if control_path.exists():
        control=pd.read_csv(control_path)
        fig,ax=plt.subplots(figsize=(6.4,3.8))
        for k,t in enumerate(FOCUS):
            landmark=native[(native.target==t)&(native.dataset=='landmark')&(native.model=='Lasso')&(native.selection=='min')].R2.iloc[0]
            ax.scatter(landmark,k,s=49,marker='s',facecolors='none',edgecolors=COLORS[t],linewidths=1.3,zorder=3)
            vals=control[(control.target==t)&(control.selection=='min')].sort_values('dataset').R2.to_numpy()
            for rep,value in enumerate(vals):
                ax.scatter(value,k+(rep-1)*.14,s=27,color=COLORS[t],marker=['o','^','D'][rep],zorder=3)
        ax.axvline(0,color='#C3CAD0',lw=.8);ax.set_yticks(range(4),[SHORT[t] for t in FOCUS]);ax.set_ylim(3.45,-.45)
        ax.set_xlabel('Out-of-fold R² (each sampled gene set)');ax.set_title('Fixed gene count: 965 readouts per set',loc='left',pad=16)
        ax.legend(handles=[Line2D([],[],marker='s',ls='',color='#555555',mfc='none',label='Landmark genes')]+[Line2D([],[],marker=m,ls='',color='#555555',label=f'Non-landmark draw {i+1}') for i,m in enumerate(['o','^','D'])],loc='upper left',bbox_to_anchor=(0,-.20),ncol=2,frameon=False,fontsize=7)
        fig.subplots_adjust(left=.19,right=.96,top=.85,bottom=.28)
        save(fig,'04_F7_size_matched_control','各 965 个读出的 landmark 集合与三次非 landmark 随机抽样对照。各集合独立完成共享划分下的 5×3 Lasso 调参，展示最低误差模型。非 landmark 抽样精确匹配每个外层测试折的 landmark 数量，内层通过与原名单求交保留划分关系；内层数量随抽样变化。各集合的训练和测试基因身份均不同。三次抽样用于描述固定总数下基因组成的敏感性，点间范围为本次抽样范围。')
    (FIG/'04_figure_manifest.json').write_text(json.dumps(MANIFEST,ensure_ascii=False,indent=2)+'\n')
    text='# 第 4 阶段图注\n\n'+'\n\n'.join(f'## {item["stem"]}\n\n{item["caption_zh"]}\n\n[PDF](../figures/{item["stem"]}.pdf) · [SVG](../figures/{item["stem"]}.svg)' for item in MANIFEST)+'\n'
    (ROOT/'reports/04_figure_captions.md').write_text(text)
    print(json.dumps({'figures':len(MANIFEST),'files':3*len(MANIFEST)}))

if __name__=='__main__':main()
