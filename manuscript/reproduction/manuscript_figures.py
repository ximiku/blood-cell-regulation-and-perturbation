"""Recompose frozen HW2 outputs for a standalone paper; no model fitting."""
from pathlib import Path
import json, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import dendrogram,leaves_list
R=Path(__file__).resolve().parents[1]; T=R/'tables'; O=R/'manuscript'; F=O/'figures';F.mkdir(exist_ok=True)
def read(n):return pd.read_csv(T/(n+'.csv'))
D=read('02_cell_dictionary').set_index('cell_type'); C=D.color.to_dict(); S=D.short_label.to_dict()
focus=['B cell','plasma cell','monocyte','erythrocyte']; ab=['B','Plasma','Mono','Ery']
mods=['zero_change','training_condition_mean','Ridge_additive','GEARS','scGPT_blood','scFoundation_GEARS']
MN=dict(zip(mods,['Zero change','Train mean','Ridge additive','GEARS','scGPT blood','scFoundation + GEARS']))
MC=dict(zip(mods,['#9AA5B1','#525D68','#0072B2','#009E73','#D55E00','#8E61A9']))
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'axes.titlesize':10,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
def panel(ax,l,title=''):
 ax.text(-.11,1.08,l,transform=ax.transAxes,fontweight='bold',fontsize=13,va='bottom');ax.set_title(title,loc='left',pad=12)
def save(f,n):
 for e in ['pdf','svg','png']:f.savefig(F/(n+'.'+e),dpi=320,bbox_inches='tight',pad_inches=.10)
 plt.close(f);print(n,flush=True)
def grid(ax):ax.grid(axis='x',color='#E8ECF0',lw=.6);ax.set_axisbelow(True)
def legend_cells(fig,y=0):
 fig.legend([Line2D([],[],marker='o',ls='',color=C[t]) for t in focus],[S[t] for t in focus],loc='lower center',bbox_to_anchor=(.5,y),ncol=4,frameon=False)

# 1: frame, PCA, aligned marker heatmap/tree and target amplitudes.
f=plt.figure(figsize=(7.3,7.7));g=f.add_gridspec(3,2,height_ratios=[.6,2.2,2.4],hspace=.65,wspace=.5)
a=f.add_subplot(g[0,:]);a.axis('off');panel(a,'A')
for x,title,txt in [(0,'Cell states','Blood target − HSC'),(.37,'Candidate discovery','L1000 dictionary + sparse regression'),(.76,'Response prediction','K562 perturbations')]:
 a.text(x,.72,title,weight='bold',color='#153A52');a.text(x,.15,txt,fontsize=8)
a.annotate('',(.35,.52),(.25,.52),arrowprops={'arrowstyle':'->','color':'#A4AFB8'});a.annotate('',(.74,.52),(.67,.52),arrowprops={'arrowstyle':'->','color':'#A4AFB8'})
a=f.add_subplot(g[1,0]);p=read('02_PCA_coordinates').set_index('cell_type');pv=read('02_PCA_explained_variance')
for t,r in p.iterrows():a.scatter(r.PC1,r.PC2,c=C[t],s=42 if t in focus else 22,edgecolors='white',lw=.5)
for t,off in zip(['hematopoietic stem cell']+focus,[(-25,-18),(9,13),(-52,-12),(9,10),(-10,-18)]):
 a.annotate(S[t],p.loc[t,['PC1','PC2']],xytext=off,textcoords='offset points',fontsize=8,color=C[t],arrowprops={'arrowstyle':'-','lw':.5,'color':'#A4AFB8'})
a.set(xlabel=f'PC1 ({pv.iloc[0].explained_variance_ratio:.1%})',ylabel=f'PC2 ({pv.iloc[1].explained_variance_ratio:.1%})');a.margins(.25);panel(a,'B','Blood-cell expression space')
a=f.add_subplot(g[1,1]);diff=read('02_HSC_difference_summary').set_index('cell_type').loc[focus]
a.barh(range(4),diff.rms_all_19994_genes,color=[C[t] for t in focus],height=.55);a.set_yticks(range(4),[S[t] for t in focus]);a.invert_yaxis();a.set_xlabel('RMS difference from HSC (input units)');panel(a,'C','Magnitude of state differences')
sub=g[2,:].subgridspec(2,1,height_ratios=[.22,1],hspace=.02);tree=f.add_subplot(sub[0]);heat=f.add_subplot(sub[1]);link=read('02_linkage_ward').values;den=dendrogram(link,ax=tree,no_labels=True,color_threshold=0,above_threshold_color='#697985');tree.axis('off');panel(tree,'D','Lineage markers aligned to cell clustering')
cells=read('02_HVG1000_zscores').columns[1:].tolist();order=[cells[i] for i in den['leaves']]
b=read('02_marker_input_expression').set_index('gene_symbol');genes=['CD34','MS4A1','CD79A','JCHAIN','MZB1','LYZ','S100A8','HBB','HBA1','ALAS2'];v=b.loc[genes,order];z=v.sub(v.mean(axis=1),axis=0).div(v.std(axis=1),axis=0)
im=heat.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-2.5,vmax=2.5);heat.set_yticks(range(len(genes)),genes);heat.set_xticks(range(len(order)),[S[t] for t in order],rotation=60,ha='right');heat.tick_params(length=0);cbax=heat.inset_axes([1.025,0,.018,1]);f.colorbar(im,cax=cbax,label='Gene z-score')
f.subplots_adjust(left=.12,right=.96,top=.96,bottom=.18);save(f,'Fig1')

# 2: all types and sizes share rows; common-axis model comparisons; marginal scatter.
f=plt.figure(figsize=(7.3,7.0));g=f.add_gridspec(2,2,height_ratios=[1.4,1],hspace=.42,wspace=.6)
sub=g[0,:].subgridspec(1,2,width_ratios=[1.3,1],wspace=.12);a=f.add_subplot(sub[0]);b=f.add_subplot(sub[1],sharey=a)
m=read('03_model_metrics').pivot(index='target',columns='model',values='R2');mdc='myeloid dendritic cell';order=m.sort_values('Lasso_min',ascending=False).index.tolist();sel=read('03_selected_models');sizes=sel[sel.scope=='full'].pivot(index='target',columns='model',values='n_active_gt_1e_6')
for model,dy,mark,col in [('Lasso_min',-.15,'o','#0072B2'),('Lasso_1se',.15,'D','#D55E00')]:
 a.scatter(m.loc[order,model],np.arange(21)+dy,s=18,marker=mark,color=col,label='Minimum MSE' if model=='Lasso_min' else '1-SE');b.scatter(sizes.loc[order,model],np.arange(21)+dy,s=18,marker=mark,color=col)
a.scatter(m.loc[order,'Mean_baseline'],range(21),s=20,marker='|',color='#888888',label='Mean baseline');a.set_yticks(range(21),[S[t]+(' *' if t==mdc else '') for t in order]);a.invert_yaxis();b.tick_params(labelleft=False);a.set_xlabel('Out-of-fold R²');b.set_xlabel('Active perturbations');a.legend(frameon=False,ncol=3,fontsize=7,loc='upper left',bbox_to_anchor=(-.5,-.14));panel(a,'A','All targets');panel(b,'B','Model size');grid(a);grid(b)
a=f.add_subplot(g[1,0]);met=read('04_model_metrics');names=['PCR','PLS','NonnegativeLasso','Lasso','ElasticNet','Ridge']
for j,t in enumerate(focus):
 for ss,mark,sz in [('min','o',20),('1se','D',14)]:
  z=met[(met.target==t)&(met.selection==ss)].set_index('model').loc[names];a.scatter(z.R2,np.arange(6)+(j-1.5)*.16+(0 if ss=='min' else .045),marker=mark,s=sz,color=C[t],facecolors=C[t] if ss=='min' else 'none')
a.set_yticks(range(6),['PCR','PLS','Nonneg. Lasso','Lasso','Elastic Net','Ridge']);a.invert_yaxis();a.set_xlabel('Out-of-fold R²');panel(a,'C','Six models on shared splits');grid(a)
sub=g[1,1].subgridspec(2,2,width_ratios=[1,.18],height_ratios=[.18,1],hspace=.03,wspace=.03);a=f.add_subplot(sub[1,0]);top=f.add_subplot(sub[0,0],sharex=a);right=f.add_subplot(sub[1,1],sharey=a)
pred=read('04_raw_oof_predictions');pred=pred[pred.target=='erythrocyte'];x=pred.observed;y=pred.Lasso_min
hb=a.hexbin(x,y,gridsize=36,bins='log',mincnt=1,cmap='Blues',linewidths=0,rasterized=True);lim=(min(x.min(),y.min())-.1,max(x.max(),y.max())+.1);a.plot(lim,lim,'--',c='#788691',lw=.8);a.set(xlim=lim,ylim=lim,xlabel='Observed target − HSC',ylabel='Predicted difference');top.hist(x,bins=60,color=C['erythrocyte']);right.hist(y,bins=60,orientation='horizontal',color=C['erythrocyte']);top.axis('off');right.axis('off');panel(top,'D','Erythrocyte');a.text(.04,.94,'R² = 0.609\nr = 0.781',transform=a.transAxes,va='top',fontsize=8)
legend_cells(f,-.005);f.subplots_adjust(left=.25,right=.96,top=.94,bottom=.12);save(f,'Fig2')

# 3: representative genes aligned coefficients/frequency/evidence, cross-cell structure and overlap.
f=plt.figure(figsize=(7.3,7.3));g=f.add_gridspec(2,2,height_ratios=[.85,1.4],wspace=.6,hspace=.6);sub=g[0,:].subgridspec(1,3,width_ratios=[1.5,.6,1],wspace=.1)
a=f.add_subplot(sub[0]);b=f.add_subplot(sub[1],sharey=a);e=f.add_subplot(sub[2],sharey=a);st=read('04_stability_all_TFs')
reps=[('erythrocyte','MYB','Erythroid'),('plasma cell','BCL6','Plasma differentiation'),('monocyte','KLF6','Macrophage state'),('monocyte','PRDM1','Macrophage state'),('B cell','ZNF22','Exploratory'),('erythrocyte','NFYB','Shared signature'),('plasma cell','OVOL2','Shared signature')]
for i,(t,p,ev) in enumerate(reps):
 r=st[(st.target==t)&(st.perturbation==p)].iloc[0];a.plot([r.coefficient_q025,r.coefficient_q975],[i,i],c=C[t]);a.scatter(r.median_coefficient,i,color=C[t],s=22);b.scatter(100*r.selection_frequency,i,color=C[t],s=22);e.text(0,i,ev,va='center',fontsize=8)
a.set_yticks(range(7),[p+' · '+S[t] for t,p,_ in reps]);a.invert_yaxis();a.axvline(0,c='#AAAAAA',lw=.7);a.set_xlabel('Subsample coefficient');b.set_xlim(75,105);b.set_xlabel('Selected (%)');b.tick_params(labelleft=False);e.axis('off');panel(a,'A','Statistical candidates and evidence');e.set_title('Evidence context',loc='left');grid(a)
top=read('04_top10_TF_candidates');genes=list(dict.fromkeys(top.perturbation));a=f.add_subplot(g[1,0]);m=st[st.perturbation.isin(genes)].pivot(index='perturbation',columns='target',values='selection_frequency').loc[genes,focus];sg=st.pivot(index='perturbation',columns='target',values='full_coefficient').loc[genes,focus];im=a.imshow(m*np.sign(sg),aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1);a.set_yticks(range(len(genes)),genes,fontsize=7);a.set_xticks(range(4),ab);panel(a,'B','Candidate programs');f.colorbar(im,ax=a,pad=.04,fraction=.05,label='Selection frequency × sign')
a=f.add_subplot(g[1,1]);over=read('04_stable_TF_overlap');v=over.pivot(index='left_target',columns='right_target',values='Jaccard').loc[focus,focus];im=a.imshow(v,cmap='Blues',vmin=0,vmax=1);a.set_xticks(range(4),ab);a.set_yticks(range(4),ab)
for i in range(4):
 for j in range(4):a.text(j,i,f'{v.iloc[i,j]:.2f}',ha='center',va='center',color='white' if i==j else '#183040',fontsize=9)
panel(a,'C','Stable identity overlap');a.text(.5,-.32,'Shared across all four:\nNFYB (−), OVOL2 (+)',ha='center',transform=a.transAxes,fontsize=10);f.colorbar(im,ax=a,fraction=.05,pad=.05,label='Jaccard')
f.subplots_adjust(left=.25,right=.94,bottom=.08,top=.93);save(f,'Fig3')

# 4: meaningful shared axes, different response definitions separated.
f=plt.figure(figsize=(7.3,6.2));g=f.add_gridspec(2,2,hspace=.65,wspace=.6);a=f.add_subplot(g[0,0]);r=read('07_core_top10_retuning_check')
for j,t in enumerate(focus):
 d=r[r.target==t];a.scatter(d.TF_stability_rank+(j-1.5)*.15,d.retuned_selection_frequency*100,color=C[t],s=22,alpha=.7)
a.set(xlabel='Within-target candidate rank',ylabel='Selected after retuning (%)',ylim=(87,102),xticks=[1,5,10]);panel(a,'A','Retuned stability (20 repeats)')
a=f.add_subplot(g[0,1]);rr=read('04_self_readout_results');d=rr[rr.status=='completed']
for j,t in enumerate(focus):
 v=d[d.target==t].relative_coefficient_change.to_numpy()*100;a.scatter(v,j+np.linspace(-.12,.12,len(v)),color=C[t],s=22)
a.axvline(0,c='#999999',lw=.7);a.set_yticks(range(4),[S[t] for t in focus]);a.invert_yaxis();a.set_xlabel('Coefficient change after deletion (%)');panel(a,'B','Candidate self-readout removal');grid(a)
a=f.add_subplot(g[1,0]);ev=read('04_candidate_evidence_matrix');cols=['raw_ElasticNet_same_active_sign','gene_scaled_Lasso_same_active_sign','landmark_Lasso_same_active_sign'];labels=['Elastic Net','Gene-scaled response','Landmark readout']
for j,t in enumerate(focus):a.plot([ev[ev.target==t][col].sum()/10*100 for col in cols],np.arange(3)+(j-1.5)*.13,'o',color=C[t],ms=4)
a.set_yticks(range(3),labels);a.invert_yaxis();a.set(xlabel='Top-10 directions retained (%)',xlim=(-5,105));panel(a,'C','Candidate agreement');grid(a)
a=f.add_subplot(g[1,1]);met=read('04_sensitivity_metrics');ref=read('04_matched_landmark_reference_metrics')
for j,t in enumerate(focus):
 x=met[(met.target==t)&(met.dataset=='landmark')&(met.selection=='min')].R2.iloc[0];y=ref[(ref.target==t)&(ref.selection=='min')].R2.iloc[0];a.plot([x,y],[j,j],c=C[t],lw=1);a.scatter(x,j,c=C[t],marker='s',s=26);a.scatter(y,j,facecolors='none',edgecolors=C[t],s=34)
a.set_yticks(range(4),ab);a.invert_yaxis();a.set_xlabel('R² on the same 965 landmark test genes');panel(a,'D','Readout coverage');grid(a);a.legend([Line2D([],[],marker='s',ls='',color='#555555'),Line2D([],[],marker='o',ls='',color='#555555',mfc='white')],['Landmark training','All-gene training'],loc='upper left',bbox_to_anchor=(0,-.25),ncol=1,fontsize=7,frameon=False)
a.get_legend().remove();legend_cells(f,-.01);f.subplots_adjust(left=.2,right=.96,bottom=.12,top=.92);save(f,'Fig4')

# 5: common/HVG metrics on shared axes, aligned generalization heatmaps.
f=plt.figure(figsize=(7.3,6.5));g=f.add_gridspec(2,2,height_ratios=[1.2,1],hspace=.68,wspace=.48);sm=read('07_model_summary');full=sm[(sm.group=='all')&(~sm.exclude_direct_targets)];a=f.add_subplot(g[0,0]);b=f.add_subplot(g[0,1])
for k,m in enumerate(mods):
 for gs,dy,mk in [('shared',-.12,'o'),('HVG1000',.12,'D')]:
  d=full[(full.model==m)&(full.gene_set==gs)]
  if d.empty:continue
  r=d.iloc[0];a.scatter(1-r.improvement_vs_zero_mean,k+dy,color=MC[m],marker=mk,s=30);b.scatter(r.Pearson_delta_mean,k+dy,color=MC[m],marker=mk,s=30)
for ax in [a,b]:ax.set_yticks(range(6),[MN[m] for m in mods]);ax.invert_yaxis();grid(ax)
b.tick_params(labelleft=False);a.axvline(1,c='#AAAAAA',ls='--',lw=.8);a.set_xlabel('MSE / zero-change MSE');b.set_xlabel('Pearson r (expression change)');panel(a,'A','Prediction error');panel(b,'B','Response pattern');a.legend([Line2D([],[],marker='o',ls='',color='#555'),Line2D([],[],marker='D',ls='',color='#555')],['Shared genes','Training HVG1000'],frameon=False,fontsize=8,loc='lower left',bbox_to_anchor=(0,-.3))
groups=['unseen_single','seen0','seen1','seen2'];mm=mods[2:]
for j,gs in enumerate(['shared','HVG1000']):
 a=f.add_subplot(g[1,j]);d=sm[(sm.gene_set==gs)&(~sm.exclude_direct_targets)];v=d[d.model.isin(mm)&d.group.isin(groups)].pivot(index='model',columns='group',values='improvement_vs_mean_mean').loc[mm,groups]*100
 im=a.imshow(v,aspect='auto',cmap='RdBu',norm=TwoSlopeNorm(vmin=-110,vcenter=0,vmax=70));a.set_xticks(range(4),['Single','Seen 0','Seen 1','Seen 2'],rotation=35,ha='right');a.set_yticks(range(4),[MN[m] for m in mm] if j==0 else [])
 for ii in range(4):
  for jj in range(4):a.text(jj,ii,f'{v.iloc[ii,jj]:.0f}',ha='center',va='center',fontsize=8,color='white' if abs(v.iloc[ii,jj]+20)>55 else '#172F42')
 panel(a,chr(67+j),'Shared genes' if j==0 else 'Training HVG1000')
f.subplots_adjust(left=.25,right=.96,bottom=.17,top=.93);cb=f.colorbar(im,ax=f.axes[-2:],orientation='horizontal',fraction=.07,pad=.28);cb.set_label('MSE improvement over training-condition mean (%)');save(f,'Fig5')

# 6: prospective gene choice with distinct prediction objectives.
f=plt.figure(figsize=(7.3,6.5));g=f.add_gridspec(3,2,height_ratios=[.55,1,1.3],hspace=.8,wspace=.58);a=f.add_subplot(g[0,:]);a.axis('off');panel(a,'A')
for x,txt in [(0,'Blood monocyte\nPRDM1: rank 8, positive'),(.39,'K562 CRISPRa\nPRDM1 held out in seed 1'),(.79,'Gene-level responses\nSingle + combination')]:a.text(x,.5,txt,fontsize=9,va='center',color='#173E53')
a.annotate('',(.37,.5),(.3,.5),arrowprops={'arrowstyle':'->'});a.annotate('',(.77,.5),(.7,.5),arrowprops={'arrowstyle':'->'})
br=read('08_PRDM1_response_bridge');a=f.add_subplot(g[1,0]);b=f.add_subplot(g[1,1]);conditions=['PRDM1+ctrl','CBFA2T3+PRDM1']
for k,m in enumerate(mods):
 for cond,dy,mk in zip(conditions,[-.12,.12],['o','D']):
  r=br[(br.model==m)&(br.condition==cond)].iloc[0];a.scatter(r.MSE,k+dy,color=MC[m],marker=mk,s=25);b.scatter(r.Pearson_delta,k+dy,color=MC[m],marker=mk,s=25)
a.set_yticks(range(6),[MN[m] for m in mods]);b.set_yticks(range(6),[])
for ax in [a,b]:ax.invert_yaxis();grid(ax)
a.set_xlabel('HVG1000 MSE');b.set_xlabel('Pearson correlation of change');panel(a,'B','Single and combination errors');panel(b,'C','Response pattern agreement');b.legend([Line2D([],[],color='#555',marker='o',ls=''),Line2D([],[],color='#555',marker='D',ls='')],['PRDM1','CBFA2T3 + PRDM1'],frameon=False,fontsize=7,loc='lower center',bbox_to_anchor=(-.35,1.15),ncol=2)
idx=np.load(R/'data/processed/stage6/evaluation_indices.npz');hg=idx['hvg1000_seed1']
for j,m in enumerate(['GEARS','scFoundation_GEARS']):
 a=f.add_subplot(g[2,j]);z=np.load(R/f'results/stage7/{m}/seed1/test_predictions.npz',allow_pickle=True);i=list(z['conditions']).index('PRDM1+ctrl');mask=np.isin(z['gene_indices'],hg);x=(z['observed'][i]-z['control_mean'])[mask];y=(z['predicted'][i]-z['control_mean'])[mask];a.scatter(x,y,s=5,c=MC[m],alpha=.35,rasterized=True);lo=min(x.min(),y.min());hi=max(x.max(),y.max());a.plot([lo,hi],[lo,hi],ls='--',c='#88949E',lw=.7);a.set(xlabel='Observed expression change',ylabel='Predicted expression change');panel(a,chr(68+j),MN[m]);r=br[(br.model==m)&(br.condition=='PRDM1+ctrl')].iloc[0];a.text(.03,.95,f'MSE = {r.MSE:.3f}\nr = {r.Pearson_delta:.3f}',transform=a.transAxes,va='top',fontsize=8)
f.subplots_adjust(left=.25,right=.96,bottom=.09,top=.94);save(f,'Fig6')

# Copy complete legacy diagnostic panels, preserving vector graphics and original scales.
for p in (R/'figures').glob('*.pdf'):
 if p.name.startswith(('02_','03_','04_','07_')):shutil.copy2(p,F/p.name)
