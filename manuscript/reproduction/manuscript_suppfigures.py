from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage,leaves_list,dendrogram
R=Path(__file__).resolve().parents[1];F=R/'manuscript/figures';T=R/'tables'
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
def save(f,n):
 for ext in ['pdf','svg','png']:f.savefig(F/(n+'.'+ext),bbox_inches='tight',dpi=320,pad_inches=.1)
 plt.close(f)
S=pd.read_csv(T/'02_cell_dictionary.csv').set_index('cell_type').short_label.to_dict();z=pd.read_csv(T/'02_heatmap_HVG100_gene_zscores.csv').set_index('gene_symbol');row=leaves_list(linkage(z,method='ward',optimal_ordering=True));col=linkage(z.T,method='ward',optimal_ordering=True);cols=leaves_list(col)
f=plt.figure(figsize=(7.4,8.8));g=f.add_gridspec(2,2,height_ratios=[.07,1],hspace=.015,wspace=.38)
for j in range(2):
 a=f.add_subplot(g[0,j]);dendrogram(col,ax=a,no_labels=True,color_threshold=0,above_threshold_color='#687985');a.axis('off');a.text(-.1,1.3,'AB'[j],transform=a.transAxes,fontweight='bold',fontsize=13)
 a=f.add_subplot(g[1,j]);v=z.iloc[row[j*50:(j+1)*50],cols];im=a.imshow(v,aspect='auto',cmap='RdBu_r',vmin=-2.5,vmax=2.5);a.set_yticks(range(50),v.index,fontsize=8);a.set_xticks(range(21),[S[t] for t in v.columns],rotation=90,fontsize=7);a.tick_params(length=0,pad=2)
f.subplots_adjust(left=.10,right=.98,top=.95,bottom=.2);cax=f.add_axes([.36,.04,.3,.015]);f.colorbar(im,cax=cax,orientation='horizontal',label='Gene z-score');save(f,'Supp_HVG100')
mods=['zero_change','training_condition_mean','Ridge_additive','GEARS','scGPT_blood','scFoundation_GEARS'];names=['Zero','Train mean','Ridge additive','GEARS','scGPT','scF + GEARS'];colors=['#9AA5B1','#525D68','#0072B2','#009E73','#D55E00','#8E61A9'];conds=['PRDM1+ctrl','CBFA2T3+PRDM1'];idx=np.load(R/'data/processed/stage6/evaluation_indices.npz');preds={m:np.load(R/f'results/stage7/{m}/seed1/test_predictions.npz') for m in mods};f,ax=plt.subplots(2,6,figsize=(10,4.1))
for row,cond in enumerate(conds):
 limit=0
 for m in mods:
  z=preds[m];i=list(z['conditions']).index(cond);mask=np.isin(z['gene_indices'],idx['hvg1000_seed1']);x=(z['observed'][i]-z['control_mean'])[mask];y=(z['predicted'][i]-z['control_mean'])[mask];limit=max(limit,np.abs(x).max(),np.abs(y).max())
 limit*=1.06
 for col,m in enumerate(mods):
  a=ax[row,col];z=preds[m];i=list(z['conditions']).index(cond);x=(z['observed'][i]-z['control_mean'])[mask];y=(z['predicted'][i]-z['control_mean'])[mask];a.scatter(x,y,s=3,c=colors[col],alpha=.45,rasterized=True);a.plot([-limit,limit],[-limit,limit],'--',c='#9CA6AF',lw=.65);a.set(xlim=(-limit,limit),ylim=(-limit,limit));a.text(-.15,1.13,chr(65+row*6+col),transform=a.transAxes,fontweight='bold',fontsize=10)
  if row==0:a.set_title(names[col],fontsize=8)
  if col==0:a.set_ylabel(('PRDM1' if row==0 else 'CBFA2T3 + PRDM1')+'\nPredicted change',fontsize=8)
  else:a.tick_params(labelleft=False)
  if row==1:a.set_xlabel('Observed',fontsize=8)
f.subplots_adjust(left=.07,right=.99,top=.88,bottom=.13,wspace=.2,hspace=.5);save(f,'Supp_PRDM1')
