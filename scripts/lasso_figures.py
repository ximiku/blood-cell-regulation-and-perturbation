"""Assemble stage-3 numerical tables and publication-oriented figure files."""
from pathlib import Path
from datetime import datetime, timezone
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap

ROOT=Path(__file__).resolve().parents[1]
TAB=ROOT/'tables';FIG=ROOT/'figures'
CFG=json.loads((ROOT/'config/analysis_config.json').read_text())
FOCUS=CFG['focus_targets'];MDC=CFG['supplementary_target'];THR=CFG['stability']['nonzero_threshold_standardized']
DICT=pd.read_csv(TAB/'02_cell_dictionary.csv').set_index('cell_type')
SHORT=DICT.short_label.to_dict();COLORS=DICT.color.to_dict()
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.65,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,
 'savefig.facecolor':'white','xtick.major.width':.6,'ytick.major.width':.6})
FIGURES=[]
DENSITY_CMAP=ListedColormap(plt.get_cmap('Blues')(np.linspace(.23,1,256)))


def csv(name,frame):
    frame.to_csv(TAB/name,index=False,encoding='utf-8-sig',float_format='%.12g')


def panel(ax,letter):
    ax.annotate(letter,(0,1),xycoords='axes fraction',xytext=(-22,10),textcoords='offset points',fontsize=12,fontweight='bold',va='bottom')


def save(fig,name,caption):
    for ext in ['pdf','svg','png']:
        fig.savefig(FIG/f'{name}.{ext}',dpi=450,bbox_inches='tight',pad_inches=.07)
    plt.close(fig);FIGURES.append({'stem':name,'caption_zh':caption,'formats':['pdf','svg','png']})


def main():
    dirs=sorted((ROOT/'results/lasso').glob('target_*'))
    assert len(dirs)==21 and all((d/'completed.json').exists() for d in dirs)
    annotation=pd.read_csv(TAB/'01_perturbation_annotations.csv')
    annotation=annotation.rename(columns={'perturbation_label':'perturbation'})
    metrics=[];folds=[];selections=[];coefs=[];runtimes=[];diagnostics=[];by_target={}
    for directory in dirs:
        meta=json.loads((directory/'completed.json').read_text());t=meta['target']
        by_target[t]=directory;runtimes.append({'target':t,'seconds':meta['seconds'],'n_fits':meta['n_fits'],'n_retried':meta['n_retried'],'all_converged':meta['all_converged']})
        metric=pd.read_csv(directory/'pooled_metrics.csv');metric['analysis_role']='supplementary' if t==MDC else 'main';metrics.append(metric)
        for file,holder in [('fold_metrics.csv',folds),('selected_models.csv',selections),('solver_diagnostics.csv',diagnostics)]:
            f=pd.read_csv(directory/file);f.insert(0,'target',t);holder.append(f)
        with np.load(directory/'model_arrays.npz') as z:
            for j,model in enumerate(['Lasso_min','Lasso_1se']):
                f=pd.DataFrame({'target':t,'model':model,'perturbation':z['perturbations'],
                     'coefficient_standardized':z['full_coef_standardized'][:,j],
                     'coefficient_input':z['full_coef_input'][:,j],
                     'effect_per_X_SD_in_response_units':z['full_coef_standardized'][:,j]*z['full_y_scale'],
                     'active_gt_1e_6':np.abs(z['full_coef_standardized'][:,j])>THR,
                     'outer_selection_fraction':np.mean(np.abs(z['outer_coef_standardized'][:,:,j])>THR,axis=0)})
                f['abs_standardized']=f.coefficient_standardized.abs()
                f['sign']=np.sign(f.coefficient_standardized).astype(int)
                f=f.merge(annotation,on='perturbation',validate='one_to_one',sort=False)
                f['rank_abs_standardized']=f.abs_standardized.rank(method='first',ascending=False).astype(int)
                coefs.append(f)
    metrics=pd.concat(metrics,ignore_index=True);folds=pd.concat(folds,ignore_index=True)
    selection=pd.concat(selections,ignore_index=True);coef=pd.concat(coefs,ignore_index=True)
    diagnostic=pd.concat(diagnostics,ignore_index=True)
    csv('03_model_metrics.csv',metrics);csv('03_fold_metrics.csv',folds);csv('03_selected_models.csv',selection)
    csv('03_all_coefficients_annotated.csv',coef);csv('03_runtime.csv',pd.DataFrame(runtimes))
    csv('03_solver_diagnostics.csv',diagnostic)
    tf=coef[(coef.primary_TF_eligible)&(coef.active_gt_1e_6)].copy()
    tf=tf.sort_values(['target','model','abs_standardized'],ascending=[True,True,False],kind='stable')
    tf['TF_rank']=tf.groupby(['target','model']).cumcount()+1
    csv('03_TF_candidates_all.csv',tf)
    csv('03_TF_top10_1se.csv',tf[(tf.model=='Lasso_1se')&(tf.TF_rank<=10)])
    csv('03_focus_TF_top10_1se.csv',tf[(tf.model=='Lasso_1se')&(tf.TF_rank<=10)&(tf.target.isin(FOCUS))])
    correlations=pd.read_csv(TAB/'03_single_perturbation_correlations.csv')
    correlations=correlations.merge(annotation,on='perturbation',validate='many_to_one',sort=False)
    csv('03_single_perturbation_correlations_annotated.csv',correlations)
    corr_top=correlations.sort_values(['target','rank_abs_Pearson']).groupby('target').head(20)
    csv('03_single_perturbation_top20.csv',corr_top)

    # F1: complete target screen and the sparsity tradeoff; mDC remains visible as supplemental.
    pivot=metrics.pivot(index='target',columns='model',values='R2')
    order=pivot.drop(index=MDC).sort_values('Lasso_min',ascending=False).index.tolist()+[MDC]
    count=selection[selection.scope=='full'].pivot(index='target',columns='model',values='n_active_gt_1e_6')
    fig,axes=plt.subplots(1,2,figsize=(7.5,6.4),gridspec_kw={'width_ratios':[1.15,.85]})
    yy=np.arange(len(order));symbols={'Lasso_min':'o','Lasso_1se':'D','Mean_baseline':'|'}
    for ax in axes:ax.axhspan(19.5,20.5,color='#F1F3F5',zorder=0)
    for model,shift,color in [('Lasso_min',-.13,'#0072B2'),('Lasso_1se',.13,'#D55E00'),('Mean_baseline',0,'#66717B')]:
        legend_label={'Lasso_min':'Lasso (minimum MSE)','Lasso_1se':'Lasso (1-SE)','Mean_baseline':'Mean baseline'}[model]
        axes[0].scatter(pivot.loc[order,model],yy+shift,s=18,marker=symbols[model],color=color,label=legend_label,zorder=3)
    axes[0].axvline(0,color='#AAB4BD',lw=.7);axes[0].set_yticks(yy,[SHORT[t]+(' [supp.]' if t==MDC else '') for t in order]);axes[0].invert_yaxis()
    axes[0].set_xlabel('Out-of-fold R²');axes[0].set_title('Reconstruction across target cell types',loc='left',pad=13)
    axes[0].legend(frameon=False,loc='lower left',bbox_to_anchor=(-.05,-.20),ncol=1)
    for model,shift,color in [('Lasso_min',-.13,'#0072B2'),('Lasso_1se',.13,'#D55E00')]:
        axes[1].scatter(count.loc[order,model],yy+shift,s=18,marker=symbols[model],color=color)
    axes[1].set_yticks(yy,[]);axes[1].invert_yaxis();axes[1].set_xlim(left=0)
    axes[1].set_xlabel('Active perturbations (full-data fit)');axes[1].set_title('Model size',loc='left',pad=13)
    for ax,letter in zip(axes,'AB'):panel(ax,letter)
    fig.subplots_adjust(left=.25,right=.98,top=.92,bottom=.21,wspace=.25)
    save(fig,'03_F1_target_screen','全部 21 个非 HSC 目标的外层五折折外 R²，以及全数据重新调参后两套模型的活跃扰动数。蓝色为最低内层 MSE 模型，朱红色为 1-SE 模型，灰色为训练折均值基线。mDC 位于底部阴影区，作为附加目标。活跃阈值为标准化系数绝对值大于 10⁻⁶；全数据模型大小用于描述，不参与折外性能计算。')

    # F2: every held-out gene is represented, using density rather than selective points.
    fig,axes=plt.subplots(2,2,figsize=(7.5,7.0))
    for ax,t,letter in zip(axes.ravel(),FOCUS,'ABCD'):
        pred=pd.read_csv(by_target[t]/'oof_predictions.csv');y=pred.observed.to_numpy();p=pred.Lasso_min.to_numpy()
        hb=ax.hexbin(y,p,gridsize=42,mincnt=1,bins='log',cmap=DENSITY_CMAP,linewidths=0,rasterized=True)
        lim=(min(y.min(),p.min()),max(y.max(),p.max()));pad=(lim[1]-lim[0])*.06
        lo,hi=lim[0]-pad,lim[1]+pad;ax.plot([lo,hi],[lo,hi],ls='--',lw=.8,color='#687681');ax.set_xlim(lo,hi);ax.set_ylim(lo,hi)
        m=metrics[(metrics.target==t)&(metrics.model=='Lasso_min')].iloc[0]
        ax.text(.04,.96,f"R² = {m.R2:.3f}\nRMSE = {m.RMSE:.3f}\nr = {m.Pearson:.3f}",transform=ax.transAxes,va='top',fontsize=7)
        ax.set_title(SHORT[t],color=COLORS[t],loc='left',pad=13);ax.set_xlabel('Observed target − HSC');ax.set_ylabel('Out-of-fold prediction')
        cb=fig.colorbar(hb,ax=ax,fraction=.045,pad=.03);cb.set_label('Genes per hexagon',fontsize=7);cb.ax.tick_params(labelsize=6)
        panel(ax,letter)
    fig.subplots_adjust(left=.09,right=.95,bottom=.09,top=.93,wspace=.52,hspace=.43)
    save(fig,'03_F2_focus_OOF_predictions','四个重点目标的最低误差 Lasso 折外预测。每个基因使用其所属外层测试折的预测，全部 11,922 个基因进入密度图。横纵轴为输入尺度差值，虚线为 y=x；颜色表示每个六边形的基因数（对数色阶，各面板独立）。R²、RMSE 和 Pearson r 由全部折外预测合并计算。')

    # F3: explanatory full-data tuning curves, explicitly separate from outer-test evaluation.
    fig,axes=plt.subplots(2,2,figsize=(7.5,5.7))
    for ax,t,letter in zip(axes.ravel(),FOCUS,'ABCD'):
        c=pd.read_csv(by_target[t]/'tuning_curves.csv');c=c[c.scope=='full'];x=c.relative_alpha.to_numpy()
        ax.plot(x,c.mean_MSE,color=COLORS[t],lw=1.3);ax.fill_between(x,c.mean_MSE-c.SE_MSE,c.mean_MSE+c.SE_MSE,color=COLORS[t],alpha=.16,lw=0)
        amin=c.loc[c.selected_min,'relative_alpha'].iloc[0];a1=c.loc[c.selected_1se,'relative_alpha'].iloc[0]
        ax.axvline(amin,color='#0072B2',ls='--',lw=.9);ax.axvline(a1,color='#D55E00',ls=':',lw=1.2)
        ax.set_xscale('log');ax.set_xlabel('Relative penalty (λ / λmax)');ax.set_ylabel('Validation MSE (input units²)')
        ax.set_title(SHORT[t],color=COLORS[t],loc='left',pad=13);panel(ax,letter)
    handles=[Line2D([],[],color='#0072B2',ls='--',label='Minimum MSE'),Line2D([],[],color='#D55E00',ls=':',label='1-SE selection')]
    fig.legend(handles=handles,loc='lower center',ncol=2,frameon=False)
    fig.subplots_adjust(left=.11,right=.98,top=.91,bottom=.16,wspace=.39,hspace=.54)
    save(fig,'03_F3_tuning_curves','全数据最终模型的三折调参曲线，展示 50 个相对惩罚强度下原尺度 MSE 的折均值及 ±1 个标准误。每个内层训练折单独计算 λmax；蓝虚线为最低误差选择，朱红点线为最强惩罚且落在最低误差加一个标准误范围内的选择。本图用于解释全数据候选模型的选择；外层测试结果独立见图 03_F1/03_F2。')

    # F4: twelve largest 1-SE coefficients, with full paths retained in the model archive.
    fig,axes=plt.subplots(2,2,figsize=(7.5,6.5))
    for ax,t,letter in zip(axes.ravel(),FOCUS,'ABCD'):
        with np.load(by_target[t]/'model_arrays.npz') as z:
            paths=z['full_path_standardized'];ratios=z['relative_alphas'];labels=z['perturbations'];chosen=z['selected_indices']
            rank=np.argsort(-np.abs(paths[:,chosen[1]]),kind='stable')[:12]
            rank=[i for i in rank if abs(paths[i,chosen[1]])>THR]
            if not rank:rank=np.argsort(-np.abs(paths[:,chosen[0]]),kind='stable')[:12]
            palette=plt.get_cmap('tab20').colors
            for k,i in enumerate(rank):ax.plot(ratios,paths[i],lw=1.0,color=palette[k],label=labels[i])
            ax.axvline(ratios[chosen[0]],color='#0072B2',ls='--',lw=.8)
            ax.axvline(ratios[chosen[1]],color='#D55E00',ls=':',lw=1)
        ax.axhline(0,color='#BFC7CE',lw=.6);ax.set_xscale('log');ax.invert_xaxis()
        ax.set_xlabel('Relative penalty (strong → weak)');ax.set_ylabel('Standardized coefficient')
        ax.set_title(SHORT[t],color=COLORS[t],loc='left',pad=13);panel(ax,letter)
        ax.legend(frameon=False,ncol=3,fontsize=5.8,loc='upper center',bbox_to_anchor=(.5,-.23),handlelength=1.3,columnspacing=.7)
    fig.subplots_adjust(left=.1,right=.98,top=.93,bottom=.17,wspace=.38,hspace=.92)
    save(fig,'03_F4_coefficient_paths','各目标全数据 1-SE 模型中按标准化系数绝对值排序的前 12 个活跃扰动的完整 Lasso 路径。若 1-SE 模型为空，则显示最低误差模型前 12 个系数。横轴由强到弱惩罚；蓝虚线和朱红点线分别标记最低误差与 1-SE 选择。图中保留正负方向；全部 2,717 个扰动的路径存于各目标 model_arrays.npz。')

    # F5: initial TF candidates, no stability claim before stage 4.
    fig,axes=plt.subplots(2,2,figsize=(7.5,6.1))
    for ax,t,letter in zip(axes.ravel(),FOCUS,'ABCD'):
        c=tf[(tf.target==t)&(tf.model=='Lasso_1se')].head(10)
        yy=np.arange(len(c));ax.barh(yy,c.coefficient_standardized,color=COLORS[t],height=.64)
        ax.set_yticks(yy,c.perturbation);ax.invert_yaxis();ax.axvline(0,color='#687681',lw=.7)
        ax.set_xlabel('Standardized coefficient (1-SE fit)');ax.set_title(SHORT[t],color=COLORS[t],loc='left',pad=13);panel(ax,letter)
        if c.empty:ax.text(.5,.5,'No active TF perturbation',transform=ax.transAxes,ha='center')
    fig.subplots_adjust(left=.13,right=.98,bottom=.10,top=.92,wspace=.62,hspace=.51)
    save(fig,'03_F5_initial_TF_candidates','四个重点目标全数据 1-SE Lasso 的 TF 初筛结果，每目标至多 10 个，按标准化系数绝对值排序，保留符号与原扰动标签。TF 资格来自第 1 阶段冻结注释。该图展示条件拟合系数，后续第 4 阶段再以重复抽样评价选择频率和符号一致性。')

    (FIG/'03_figure_manifest.json').write_text(json.dumps(FIGURES,ensure_ascii=False,indent=2)+'\n')
    (ROOT/'reports/03_figure_captions.md').write_text('# 第 3 阶段图注\n\n'+'\n\n'.join(f"## {f['stem']}\n\n{f['caption_zh']}" for f in FIGURES)+'\n')
    summary={'generated_utc':datetime.now(timezone.utc).isoformat(),'targets':len(dirs),'focus_targets':FOCUS,
             'main_targets':20,'supplementary_target':MDC,'metrics':metrics.to_dict('records'),
             'runtime':runtimes,'n_solver_fits':len(diagnostic),'n_nonconverged':int((~diagnostic.converged).sum()),
             'n_retried':int(diagnostic.retried.sum()),'n_low_alpha_boundary_selected':int(selection.selection_at_lowest_alpha.sum()),
             'full_data_model_sizes':selection[selection.scope=='full'].to_dict('records'),
             'figure_count':len(FIGURES),'stability_analysis_started':False}
    (ROOT/'logs/stage3_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['metrics','runtime','full_data_model_sizes']},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
