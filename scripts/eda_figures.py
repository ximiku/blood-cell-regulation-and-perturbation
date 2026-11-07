"""Stage 2: blood-expression EDA, reproducible figures and Prism source data."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME',str(ROOT/'.cache'))
import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, cut_tree
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
import seaborn as sns

CONFIG=json.loads((ROOT/'config/analysis_config.json').read_text())
FIG=ROOT/'figures'; TAB=ROOT/'tables'; PRISM=ROOT/'prism'; REPORT=ROOT/'reports'
for p in [FIG,TAB,PRISM,REPORT,ROOT/'logs']: p.mkdir(parents=True,exist_ok=True)
HSC=CONFIG['reference_cell']; MDC=CONFIG['supplementary_target']; FOCUS=CONFIG['focus_targets']
SHORT={HSC:'HSC','common myeloid progenitor':'CMP','erythrocyte':'Erythrocyte','platelet':'Platelet',
       'macrophage':'Macrophage','B cell':'B cell','monocyte':'Monocyte','natural killer cell':'NK cell',
       'CD4-positive, alpha-beta T cell':'CD4 T','CD8-positive, alpha-beta T cell':'CD8 T',
       'basophil':'Basophil','neutrophil':'Neutrophil',MDC:'mDC','plasmacytoid dendritic cell':'pDC',
       'plasma cell':'Plasma cell','mature NK T cell':'NK T','regulatory T cell':'Treg',
       'classical monocyte':'Classical mono.','non-classical monocyte':'Non-classical mono.',
       'naive thymus-derived CD4-positive, alpha-beta T cell':'Naive CD4 T',
       'intermediate monocyte':'Intermediate mono.','hematopoietic precursor cell':'Precursor'}
COLORS={HSC:'#40464D','B cell':'#0072B2','plasma cell':'#AA4499','monocyte':'#E69F00','erythrocyte':'#D55E00',MDC:'#222222'}
OTHER='#BFC7CE'
CMAP=LinearSegmentedColormap.from_list('expression_diverging',['#2864A1','#FAFAFA','#B73745'])
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
                     'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.65,
                     'lines.linewidth':1.3,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
                     'savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,
                     'xtick.major.width':.6,'ytick.major.width':.6})
np.random.seed(CONFIG['seed'])
figure_records=[]


def save_csv(name, frame, index=False):
    frame.to_csv(TAB/name,index=index,encoding='utf-8-sig',float_format='%.12g')


def panel(ax, letter):
    ax.annotate(letter,xy=(0,1),xycoords='axes fraction',xytext=(-22,10),
                textcoords='offset points',fontsize=12,fontweight='bold',va='bottom')


def save_figure(fig, stem, caption):
    for ext in ['pdf','svg','png']:
        fig.savefig(FIG/f'{stem}.{ext}',dpi=450 if ext=='png' else 300,bbox_inches='tight',pad_inches=.06)
    figure_records.append({'stem':stem,'caption_zh':caption,'formats':['pdf','svg','png']})
    plt.close(fig)


def row_z(frame):
    scale=frame.std(axis=1,ddof=1)
    assert (scale>0).all()
    return frame.sub(frame.mean(axis=1),axis=0).div(scale,axis=0)


def write_prism(rms):
    root=ET.Element('GraphPadPrismFile',PrismXMLVersion='5.00')
    created=ET.SubElement(root,'Created')
    ET.SubElement(created,'OriginalVersion',CreatedByProgram='HW2 data export',CreatedByVersion='1.0',DateTime=datetime.now(timezone.utc).isoformat())
    seq=ET.SubElement(root,'TableSequence',Selected='1');ET.SubElement(seq,'Ref',ID='Table0',Selected='1')
    table=ET.SubElement(root,'Table',ID='Table0',XFormat='none',YFormat='replicates',Replicates='1',TableType='OneWay',EVFormat='AsteriskAfterNumber')
    ET.SubElement(table,'Title').text='RMS expression difference vs HSC'
    for c in FOCUS:
        column=ET.SubElement(table,'YColumn',Width='100',Decimals='5',Subcolumns='1')
        ET.SubElement(column,'Title').text=SHORT[c]
        sub=ET.SubElement(column,'Subcolumn');ET.SubElement(sub,'d').text=f'{rms[c]:.12g}'
    ET.indent(root)
    ET.ElementTree(root).write(PRISM/'HSC_difference_data.pzfx',encoding='utf-8',xml_declaration=True)
    pd.DataFrame({SHORT[c]:[rms[c]] for c in FOCUS}).to_csv(PRISM/'HSC_difference_data.tsv',sep='\t',index=False,float_format='%.12g')


def pca_scatter(ax, coords, variance, x=0, y=1, extra=None, title=None):
    for c in coords.index:
        focus=c in FOCUS or c==HSC
        ax.scatter(coords.loc[c,f'PC{x+1}'],coords.loc[c,f'PC{y+1}'],
                   s=42 if focus else 22,color=COLORS.get(c,OTHER),
                   marker='*' if c==HSC else 'o',linewidth=.55,edgecolor='white',zorder=4 if focus else 2)
    # Fixed labels emphasize the predeclared targets; complete cell labels are in the heatmap and tables.
    if y==1:
        offsets={HSC:(-36,-20),'B cell':(18,20),'plasma cell':(-58,-15),'monocyte':(15,15),'erythrocyte':(-7,-20),'platelet':(12,-1),'neutrophil':(10,0)}
    else:
        offsets={HSC:(-42,22),'B cell':(13,16),'plasma cell':(-61,-1),'monocyte':(16,8),'erythrocyte':(-5,-21),'platelet':(-6,14),'neutrophil':(10,-7)}
    for c,offset in offsets.items():
        if c not in coords.index: continue
        ax.annotate(SHORT[c],(coords.loc[c,f'PC{x+1}'],coords.loc[c,f'PC{y+1}']),xytext=offset,textcoords='offset points',
                    fontsize=7,color=COLORS.get(c,'#66717B'),ha='left',va='center',
                    arrowprops={'arrowstyle':'-','lw':.5,'color':'#9EA7AF'},zorder=5)
    if extra is not None:
        ax.scatter(extra[x],extra[y],marker='D',s=43,facecolors='none',edgecolors='#222222',linewidth=1.0,zorder=6)
        ax.annotate('mDC',extra[[x,y]],xytext=(10,-8),textcoords='offset points',fontsize=7)
    ax.axhline(0,color='#DDE2E7',lw=.55,zorder=0);ax.axvline(0,color='#DDE2E7',lw=.55,zorder=0)
    ax.set_xlabel(f'PC{x+1} ({variance[x]*100:.1f}%)');ax.set_ylabel(f'PC{y+1} ({variance[y]*100:.1f}%)')
    ax.margins(x=.28,y=.23)
    if title: ax.set_title(title,loc='left',pad=13)


def main():
    input_path=ROOT/'data/processed/analysis_matrices.h5'
    with h5py.File(input_path,'r') as f:
        blood=pd.DataFrame(f['blood/expression'][:],index=f['blood/genes'].asstr()[:],columns=f['blood/cell_types'].asstr()[:])
        model_genes=f['model/genes'].asstr()[:].tolist()
    main_cells=[c for c in blood.columns if c!=MDC]
    main=blood[main_cells]
    variance=main.var(axis=1,ddof=1).sort_values(ascending=False,kind='mergesort')
    assert (variance.iloc[:1000]>0).all()
    ranked=pd.DataFrame({'gene_symbol':variance.index,'variance_21_cells':variance.values,
                         'rank':np.arange(1,len(variance)+1),
                         'ribosomal_prefix':variance.index.str.match(r'^RP[SL]\d'),
                         'globin_gene':variance.index.isin(['HBA1','HBA2','HBB','HBD','HBG1','HBG2'])})
    hgnc=pd.read_csv(ROOT/'data/references/hgnc_complete_set.txt',sep='\t',low_memory=False).set_index('symbol')
    ranked['HGNC_gene_group']=ranked.gene_symbol.map(hgnc.gene_group).fillna('')
    ranked['HGNC_ribosomal_annotation']=ranked.HGNC_gene_group.str.contains('ribosomal',case=False)
    save_csv('02_gene_variance.csv',ranked)
    for n in [100,500,1000]:save_csv(f'02_HVG_{n}.csv',ranked.head(n))
    Z=row_z(main.loc[variance.index[:1000]])
    ward=linkage(Z.T,method='ward',metric='euclidean',optimal_ordering=True)
    corr_link=linkage(pdist(Z.T,metric='correlation'),method='average',optimal_ordering=True)
    order=[main_cells[i] for i in leaves_list(ward)]
    pca=PCA(svd_solver='full')
    scores=pca.fit_transform(Z.T)
    coords=pd.DataFrame(scores[:,:10],index=main_cells,columns=[f'PC{i}' for i in range(1,11)])
    save_csv('02_PCA_coordinates.csv',coords.rename_axis('cell_type').reset_index())
    pvar=pd.DataFrame({'PC':np.arange(1,len(pca.explained_variance_ratio_)+1),
                       'explained_variance_ratio':pca.explained_variance_ratio_,
                       'cumulative':np.cumsum(pca.explained_variance_ratio_)})
    save_csv('02_PCA_explained_variance.csv',pvar)
    save_csv('02_PCA_loadings.csv',pd.DataFrame(pca.components_[:5].T,index=Z.index,columns=[f'PC{i}' for i in range(1,6)]).rename_axis('gene_symbol').reset_index())
    save_csv('02_cell_dictionary.csv',pd.DataFrame({'cell_type':blood.columns,'short_label':[SHORT[c] for c in blood.columns],
             'color':[COLORS.get(c,OTHER) for c in blood.columns], 'main_EDA':[c!=MDC for c in blood.columns], 'focus':[c in FOCUS for c in blood.columns]}))
    save_csv('02_HVG1000_zscores.csv',Z.rename_axis('gene_symbol').reset_index())
    corr=Z.corr();dist=pd.DataFrame(squareform(pdist(Z.T)),index=main_cells,columns=main_cells)
    save_csv('02_HVG1000_gene_zscore_Pearson.csv',corr.rename_axis('cell_type').reset_index())
    save_csv('02_HVG1000_euclidean_distance.csv',dist.rename_axis('cell_type').reset_index())
    save_csv('02_linkage_ward.csv',pd.DataFrame(ward,columns=['left_node','right_node','distance','count']))
    save_csv('02_linkage_correlation_average.csv',pd.DataFrame(corr_link,columns=['left_node','right_node','distance','count']))
    nearest=[]
    for c in main_cells:
        other=[a for a in main_cells if a!=c]
        nearest.append({'cell_type':c,'euclidean_nearest':dist.loc[c,other].idxmin(),'euclidean_distance':dist.loc[c,other].min(),
                        'correlation_nearest':corr.loc[c,other].idxmax(),'correlation':corr.loc[c,other].max()})
    save_csv('02_nearest_cell_profiles.csv',pd.DataFrame(nearest))

    # Figure 1: zero fraction and full-distribution ECDF, without extra expression transformation.
    zeros=(blood==0).mean().sort_values(ascending=False,kind='mergesort')
    fig,axes=plt.subplots(1,2,figsize=(7.5,5.3),gridspec_kw={'width_ratios':[1.05,1]})
    ax=axes[0];yy=np.arange(len(zeros))
    ax.barh(yy,100*zeros.values,color=[COLORS.get(c,OTHER) for c in zeros.index],height=.65)
    ax.set_yticks(yy,[SHORT[c] for c in zeros.index]);ax.invert_yaxis();ax.set_xlim(0,105)
    ax.set_xlabel('Zero-valued genes (%)');ax.set_title('Expression coverage',loc='left',pad=13)
    for i,(c,val) in enumerate(zeros.items()):
        if c in FOCUS+[HSC,MDC]:ax.text(val*100+1.2,i,f'{val*100:.1f}',va='center',fontsize=6.5)
    panel(ax,'A')
    ax=axes[1]
    ecdf=[]
    for c in [HSC]+FOCUS:
        values=np.sort(blood[c].to_numpy());fraction=np.arange(1,len(values)+1)/len(values)
        ax.plot(values,fraction,color=COLORS[c],label=SHORT[c],lw=1.5)
        ecdf.append(pd.DataFrame({'cell_type':c,'expression':values,'cumulative_fraction':fraction}))
    ax.set_xscale('symlog',linthresh=.001);ax.set_xlim(0,10);ax.set_ylim(0,1.01)
    ax.set_xticks([0,.001,.01,.1,1,10],['0','0.001','0.01','0.1','1','10'])
    ax.set_xlabel('Expression (input scale; symlog axis)');ax.set_ylabel('Cumulative fraction of genes')
    ax.set_title('Expression distributions',loc='left',pad=13);ax.legend(frameon=False,loc='lower right')
    panel(ax,'B');fig.subplots_adjust(left=.19,right=.98,bottom=.13,top=.91,wspace=.49)
    save_csv('02_focus_expression_ECDF.csv',pd.concat(ecdf,ignore_index=True))
    save_figure(fig,'02_F1_expression_distributions','A：全部 22 类细胞的零值比例，分母为 19,994 个去重基因。B：HSC 与四个重点目标的基因表达累积分布；symlog 仅为横轴显示方式，原始表达数值保持输入尺度。曲线描述跨基因分布。')

    # Figure 2: HVG100 heatmap, row/column trees both based on the displayed matrix.
    Z100=row_z(main.loc[variance.index[:100]])
    hrow=linkage(Z100,method='ward',metric='euclidean',optimal_ordering=True)
    hcol=linkage(Z100.T,method='ward',metric='euclidean',optimal_ordering=True)
    display=Z100.rename(columns=SHORT)
    strip=pd.Series([COLORS.get(c,OTHER) for c in main_cells],index=display.columns,name='Targets')
    g=sns.clustermap(display,row_linkage=hrow,col_linkage=hcol,col_colors=strip,cmap=CMAP,
                    vmin=-2.5,vmax=2.5,figsize=(7.7,10.9),dendrogram_ratio=(.08,.07),
                    colors_ratio=.018,cbar_pos=(.04,.91,.025,.06),xticklabels=True,yticklabels=True,
                    cbar_kws={'ticks':[-2.5,0,2.5]})
    g.ax_heatmap.tick_params(axis='y',labelsize=5.8,length=0,pad=2)
    g.ax_heatmap.tick_params(axis='x',labelsize=7,length=0,pad=3)
    plt.setp(g.ax_heatmap.get_xticklabels(),rotation=60,ha='right')
    g.ax_heatmap.set_xlabel('');g.ax_heatmap.set_ylabel('')
    g.ax_cbar.set_title('Gene\nz-score',fontsize=7,pad=4)
    g.fig.suptitle('High-variance expression programs',x=.51,y=1.005,fontsize=11,fontweight='bold')
    g.fig.legend(handles=[Patch(facecolor=COLORS[c],label=SHORT[c]) for c in [HSC]+FOCUS]+[Patch(facecolor=OTHER,label='Other types')],
                 loc='lower center',bbox_to_anchor=(.5,-.055),ncol=3,frameon=False,fontsize=7)
    save_csv('02_heatmap_HVG100_gene_zscores.csv',Z100.rename_axis('gene_symbol').reset_index())
    save_csv('02_heatmap_HVG100_order.csv',pd.DataFrame({'cell_type':[main_cells[i] for i in leaves_list(hcol)]}))
    save_figure(g.fig,'02_F2_HVG100_heatmap','以除 mDC 外的 21 类细胞计算逐基因方差，选前 100 个基因；每行按这 21 类细胞标准化为 z-score。行与列均按显示矩阵的欧氏距离/Ward 法聚类。色阶固定为 ±2.5，超出部分在显示中截断；原数值完整保存。上方色条标出 HSC 和预先指定的四个重点目标。')

    # Figure 3: two PCA views plus variance explained; labels focus on predeclared targets.
    fig=plt.figure(figsize=(7.5,5.7));gs=fig.add_gridspec(2,2,height_ratios=[3,1.05],hspace=.58,wspace=.36)
    ax=fig.add_subplot(gs[0,0]);pca_scatter(ax,coords,pca.explained_variance_ratio_,0,1,title='Principal components 1 and 2');panel(ax,'A')
    ax=fig.add_subplot(gs[0,1]);pca_scatter(ax,coords,pca.explained_variance_ratio_,0,2,title='Principal components 1 and 3');panel(ax,'B')
    ax=fig.add_subplot(gs[1,:]);r=pca.explained_variance_ratio_[:10]*100
    ax.bar(np.arange(1,11),r,color=['#40464D']*3+['#BFC7CE']*7,width=.65)
    ax.set_xticks(np.arange(1,11));ax.set_xlabel('Principal component');ax.set_ylabel('Variance (%)');ax.set_ylim(0,max(r)*1.22)
    for i,v in enumerate(r[:3],1):ax.text(i,v+.5,f'{v:.1f}',ha='center',fontsize=7)
    panel(ax,'C');fig.subplots_adjust(left=.1,right=.98,bottom=.10,top=.93)
    save_figure(fig,'02_F3_PCA','前 1,000 个高变基因逐行标准化后，对 21 类细胞平均表达谱做 PCA。每个点是一类细胞的平均谱；A、B 分别显示 PC1–PC2 与 PC1–PC3，C 显示各轴解释方差。彩色点为预先指定的四个重点目标，HSC 为深灰星形，其余类型为浅灰。完整细胞坐标表提供所有点的身份。')

    # Figure 4: global difference magnitude and a prespecified lineage marker panel.
    delta=main.subtract(main[HSC],axis=0)
    rms=np.sqrt(delta.pow(2).mean())
    delta_rows=[]
    for c in main_cells:
        delta_rows.append({'cell_type':c,'rms_all_19994_genes':float(rms[c]),
                           'mean_abs_difference_all_genes':float(delta[c].abs().mean()),
                           'rms_model_11922_genes':float(np.sqrt(delta.loc[model_genes,c].pow(2).mean())),
                           'mean_abs_difference_model_genes':float(delta.loc[model_genes,c].abs().mean())})
    save_csv('02_HSC_difference_summary.csv',pd.DataFrame(delta_rows))
    save_csv('02_focus_HSC_differences.csv',delta[FOCUS].rename_axis('gene_symbol').reset_index())
    marker_groups={'HSC / progenitor':['CD34','KIT','HLF','GATA2'],
                   'B cell':['MS4A1','CD79A','CD79B','CD74','HLA-DRA'],
                   'Plasma cell':['MZB1','SDC1','JCHAIN','PRDM1','XBP1','DERL3'],
                   'Monocyte':['LYZ','LST1','FCN1','S100A8','S100A9','CTSS'],
                   'Erythroid':['GYPA','HBA1','HBA2','HBB','ALAS2','AHSP','KLF1','GATA1']}
    marker_records=[];present=[]
    for group,names in marker_groups.items():
        for gene in names:
            ok=gene in main.index and main.loc[gene].std(ddof=1)>0
            marker_records.append({'group':group,'gene_symbol':gene,'available_variable':ok,'selection':'prespecified_reference_panel'})
            if ok:present.append(gene)
    marker_Z=row_z(main.loc[present]);marker_delta=marker_Z[FOCUS].subtract(marker_Z[HSC],axis=0)
    save_csv('02_marker_panel_manifest.csv',pd.DataFrame(marker_records))
    save_csv('02_marker_input_expression.csv',main.loc[present].rename_axis('gene_symbol').reset_index())
    save_csv('02_marker_gene_zscore_HSC_difference.csv',marker_delta.rename_axis('gene_symbol').reset_index())
    fig=plt.figure(figsize=(7.5,6.8));gs=fig.add_gridspec(1,2,width_ratios=[1,1.15],wspace=.65)
    ax=fig.add_subplot(gs[0,0]);sorted_rms=rms.drop(HSC).sort_values()
    yy=np.arange(len(sorted_rms));ax.hlines(yy,0,sorted_rms.values,color='#D8DFE5',lw=1)
    ax.scatter(sorted_rms.values,yy,c=[COLORS.get(c,OTHER) for c in sorted_rms.index],s=31,edgecolors='white',linewidth=.4)
    ax.set_yticks(yy,[SHORT[c] for c in sorted_rms.index]);ax.set_xlim(0,.42)
    ax.set_xlabel('RMS difference vs HSC\n(input units; all 19,994 genes)');ax.set_title('Magnitude of expression differences',loc='left',pad=13)
    panel(ax,'A')
    ax=fig.add_subplot(gs[0,1]);sns.heatmap(marker_delta,cmap=CMAP,vmin=-3,vmax=3,ax=ax,cbar=True,
                         cbar_kws={'label':'Change in gene z-score vs HSC','shrink':.6,'ticks':[-3,0,3]},
                         xticklabels=[SHORT[c] for c in FOCUS],yticklabels=True,linewidths=.25,linecolor='white')
    ax.set_xticklabels([SHORT[c] for c in FOCUS],rotation=45,ha='right');ax.tick_params(axis='y',labelsize=7,length=0)
    ax.set_xlabel('');ax.set_ylabel('');ax.set_title('Prespecified lineage-marker panel',loc='left',pad=13)
    cursor=0
    for names in marker_groups.values():
        cursor+=sum(gene in present for gene in names)
        if cursor<len(present):ax.axhline(cursor,color='#6B7780',lw=.8)
    panel(ax,'B');fig.subplots_adjust(left=.19,right=.97,bottom=.16,top=.92)
    save_figure(fig,'02_F4_HSC_differences','A：相对 HSC 的表达差值均方根（RMS），使用全部 19,994 个去重基因，衡量整体差异幅度。B：预先指定的谱系参照基因；先用 21 类细胞计算逐基因 z-score，再展示各目标减 HSC 的变化。分隔线对应 HSC/祖细胞、B、浆细胞、单核和红系面板。色阶 ±3；细节数值保存于表格。这些摘要描述平均表达谱，图中不设置生物重复误差条。')
    write_prism(rms)

    # Figure 5: clustering comparison using the same HVG1000 and gene z-scores.
    fig,axes=plt.subplots(1,2,figsize=(7.5,5.4))
    for ax,L,title,letter in zip(axes,[ward,corr_link],['Euclidean / Ward','Correlation / average'],['A','B']):
        dendrogram(L,orientation='left',labels=[SHORT[c] for c in main_cells],ax=ax,
                   color_threshold=0,above_threshold_color='#56626D',leaf_font_size=7)
        ax.set_title(title,loc='left',pad=13);ax.set_xlabel('Linkage distance')
        for label in ax.get_yticklabels():
            reverse={v:k for k,v in SHORT.items()};c=reverse[label.get_text()]
            label.set_color(COLORS.get(c,'#4E5962'));label.set_fontweight('bold' if c in FOCUS+[HSC] else 'normal')
        ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0);panel(ax,letter)
    fig.subplots_adjust(left=.06,right=.82,bottom=.12,top=.91,wspace=1.1)
    save_figure(fig,'02_F5_clustering_comparison','两种聚类使用相同的前 1,000 个高变基因及相同的逐基因标准化矩阵。A 为欧氏距离/Ward，B 为 Pearson 相关距离/平均连接。叶节点顺序由各自算法确定。颜色强调 HSC 与四个重点目标；分支长度分别采用各方法自己的距离尺度。')

    # Sensitivity: fixed-vs-refitted mDC inclusion, and prespecified HVG counts / preprocessing.
    mdc_z=(blood.loc[Z.index,MDC]-main.loc[Z.index].mean(axis=1))/main.loc[Z.index].std(axis=1,ddof=1)
    mdc_projection=pca.transform(pd.DataFrame([mdc_z.to_numpy()],columns=Z.index,index=[MDC]))[0]
    all_var=blood.var(axis=1,ddof=1).sort_values(ascending=False,kind='mergesort')
    Zall=row_z(blood.loc[all_var.index[:1000]])
    pca_all=PCA(svd_solver='full');all_scores=pca_all.fit_transform(Zall.T)
    all_coords=pd.DataFrame(all_scores[:,:3],index=blood.columns,columns=['PC1','PC2','PC3'])
    save_csv('02_PCA_all22_refit_coordinates.csv',all_coords.rename_axis('cell_type').reset_index())
    save_csv('02_mDC_fixed_projection.csv',pd.DataFrame([mdc_projection[:3]],columns=['PC1','PC2','PC3']))
    fig,axes=plt.subplots(1,2,figsize=(7.5,3.6))
    pca_scatter(axes[0],coords,pca.explained_variance_ratio_,0,1,extra=mdc_projection,title='mDC projected into fixed 21-type PCA');panel(axes[0],'A')
    pca_scatter(axes[1],all_coords.drop(index=MDC),pca_all.explained_variance_ratio_,0,1,extra=all_scores[list(blood.columns).index(MDC)],title='PCA refitted with all 22 types');panel(axes[1],'B')
    fig.subplots_adjust(left=.1,right=.97,bottom=.2,top=.86,wspace=.4)
    save_figure(fig,'02_S1_mDC_sensitivity','A：使用主分析已经固定的 HVG、逐基因均值/标准差和 PCA 轴，将 mDC 投影到同一空间，空心菱形为 mDC。B：纳入全部 22 类细胞后重新选择 HVG、标准化并拟合 PCA，主成分轴由该数据重新确定。两面板用于显示高零值细胞对整体结构的影响。')
    sensitivity=[];cluster_assignment=[]
    alternatives={}
    for n in [100,500,1000]:
        zn=row_z(main.loc[variance.index[:n]])
        alternatives[f'HVG{n}_z_Ward']=linkage(zn.T,method='ward',optimal_ordering=True)
    alternatives['HVG1000_z_correlation_average']=corr_link
    alternatives['HVG1000_input_Ward']=linkage(main.loc[variance.index[:1000]].T,method='ward',optimal_ordering=True)
    for k in range(3,9):
        reference=cut_tree(ward,n_clusters=k).ravel()
        for method,L in alternatives.items():
            labels=cut_tree(L,n_clusters=k).ravel()
            sensitivity.append({'method':method,'n_clusters':k,'ARI_vs_HVG1000_z_Ward':adjusted_rand_score(reference,labels)})
            for c,label in zip(main_cells,labels):cluster_assignment.append({'method':method,'n_clusters':k,'cell_type':c,'cluster':int(label)})
    sensitivity=pd.DataFrame(sensitivity);save_csv('02_clustering_sensitivity.csv',sensitivity)
    save_csv('02_cluster_assignments.csv',pd.DataFrame(cluster_assignment))
    fig,ax=plt.subplots(figsize=(6.7,2.7))
    method_labels={'HVG1000_z_Ward':'1,000 genes | z-score | Ward (reference)',
                   'HVG500_z_Ward':'500 genes | z-score | Ward',
                   'HVG100_z_Ward':'100 genes | z-score | Ward',
                   'HVG1000_z_correlation_average':'1,000 genes | z-score | correlation / average',
                   'HVG1000_input_Ward':'1,000 genes | input scale | Ward'}
    sensitivity_grid=sensitivity.pivot(index='method',columns='n_clusters',values='ARI_vs_HVG1000_z_Ward')
    sensitivity_grid=sensitivity_grid.loc[list(method_labels)].rename(index=method_labels)
    sns.heatmap(sensitivity_grid,
                cmap='Blues',vmin=0,vmax=1,annot=True,fmt='.2f',linewidths=.5,ax=ax,
                cbar_kws={'label':'Adjusted Rand index','ticks':[0,.5,1]})
    ax.set_xlabel('Number of clusters');ax.set_ylabel('');ax.set_yticklabels(ax.get_yticklabels(),rotation=0)
    ax.set_title('Clustering sensitivity to gene count, scaling and distance',loc='left',pad=12)
    save_figure(fig,'02_S2_clustering_sensitivity','预先固定 k=3–8，比较 HVG100/500/1000、标准化方式及距离/连接方法。ARI 以 HVG1000、逐基因 z-score、欧氏距离/Ward 为参考，衡量同一批细胞类型的分组一致性。完整结果并列展示，供评估图形结构对分析选择的敏感程度。')

    caption_text='# 第 2 阶段图注与图形规范\n\n'
    for item in figure_records:caption_text+=f"## {item['stem']}\n\n{item['caption_zh']}\n\n"
    caption_text+='统一字体 Arial；主图字体约 7–9 pt、HVG100 热图基因标签 5.8 pt；颜色固定；矢量 PDF/SVG 和 450 dpi PNG。色阶截断仅用于显示，源数据完整保存。Prism 数据源中的四个数值为跨基因 RMS 摘要，每种细胞类型一个值。\n'
    (REPORT/'02_figure_captions.md').write_text(caption_text)
    (FIG/'figure_manifest.json').write_text(json.dumps(figure_records,ensure_ascii=False,indent=2)+'\n')
    top100=ranked.head(100)
    summary={'main_cell_types':len(main_cells),'supplementary_cell':MDC,'blood_genes':len(blood),
             'pca_HVG':1000,'pca_variance_ratio_first3':pca.explained_variance_ratio_[:3].tolist(),
             'pca_first2_cumulative':float(pca.explained_variance_ratio_[:2].sum()),
             'pca_first3_cumulative':float(pca.explained_variance_ratio_[:3].sum()),
             'focus_rms_all_genes':{c:float(rms[c]) for c in FOCUS},
             'HVG100_ribosomal_prefix_count':int(top100.ribosomal_prefix.sum()),
             'HVG100_HGNC_ribosomal_count':int(top100.HGNC_ribosomal_annotation.sum()),
             'HVG100_globin_count':int(top100.globin_gene.sum()),
             'marker_panel_requested':len(marker_records),'marker_panel_available':len(present),
             'marker_missing':[r['gene_symbol'] for r in marker_records if not r['available_variable']],
             'hsc_nearest':next(r for r in nearest if r['cell_type']==HSC),
             'pca_all22_variance_first3':pca_all.explained_variance_ratio_[:3].tolist(),
             'clustering_sensitivity':sensitivity.to_dict('records'),
             'figure_count':len(figure_records),'training_started':False,
             'generated_utc':datetime.now(timezone.utc).isoformat()}
    (ROOT/'logs/stage2_eda_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='clustering_sensitivity'},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
