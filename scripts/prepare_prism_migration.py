"""Reuse saved scientific tables as lossless Prism inputs; no model fitting."""
from pathlib import Path
import xml.etree.ElementTree as E
from datetime import datetime, timezone
import json, hashlib
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
R=Path(__file__).resolve().parents[1];P=R/'prism/migration';P.mkdir(exist_ok=True)
T=R/'tables';D=pd.read_csv(T/'02_cell_dictionary.csv').set_index('cell_type');short=D.short_label.to_dict()
focus=['B cell','plasma cell','monocyte','erythrocyte'];spec=[]
def root():
    r=E.Element('GraphPadPrismFile',PrismXMLVersion='5.00'); c=E.SubElement(r,'Created');E.SubElement(c,'OriginalVersion',CreatedByProgram='HW2 lossless data export',CreatedByVersion='1.0',DateTime=datetime.now(timezone.utc).isoformat());E.SubElement(r,'TableSequence',Selected='1');return r

def values(parent,seq):
    s=E.SubElement(parent,'Subcolumn')
    for v in seq:E.SubElement(s,'d').text='' if pd.isna(v) else (format(float(v),'.17g') if isinstance(v,(float,int,np.number)) else str(v))

def table(r,title,frame,kind='TwoWay',x=None,interval=False):
    i=len(r.findall('Table'));tid=f'Table{i}';E.SubElement(r.find('TableSequence'),'Ref',ID=tid,Selected='1' if i==0 else '0')
    t=E.SubElement(r,'Table',ID=tid,XFormat='numbers' if x is not None else 'none',YFormat='upper-lower-limits' if interval else 'replicates',TableType=kind,EVFormat='AsteriskAfterNumber')
    if not interval:t.set('Replicates','1')
    E.SubElement(t,'Title').text=title
    rows=E.SubElement(t,'RowTitlesColumn',Width='130');values(rows,frame.index)
    if x is not None:
        c=E.SubElement(t,'XColumn',Width='100',Decimals='8',Subcolumns='1');E.SubElement(c,'Title').text=x.name;values(c,x)
    if interval:
        c=E.SubElement(t,'YColumn',Width='240',Decimals='8',Subcolumns='3');E.SubElement(c,'Title').text='Median and 2.5-97.5 percentiles'
        for name in ['median_coefficient','coefficient_q975','coefficient_q025']:values(c,frame[name])
    else:
        for name in frame.columns:
            c=E.SubElement(t,'YColumn',Width='100',Decimals='8',Subcolumns='1');E.SubElement(c,'Title').text=str(name);values(c,frame[name])
    return t

def save(r,stem):
    E.indent(r);E.ElementTree(r).write(P/f'{stem}.pzfx',encoding='utf-8',xml_declaration=True)

def csv(f,stem):f.to_csv(P/f'{stem}.tsv',sep='\t',float_format='%.17g')

pc=pd.read_csv(T/'02_PCA_coordinates.csv').set_index('cell_type');var=pd.read_csv(T/'02_PCA_explained_variance.csv')
r=root();categories=['hematopoietic stem cell']+focus+['Other cell types']
for y in ['PC2','PC3']:
    frame=pd.DataFrame(index=[short[c] for c in pc.index])
    for c in categories:
        sel=~pc.index.isin(categories[:-1]) if c=='Other cell types' else pc.index==c
        frame[short.get(c,c)]=np.where(sel,pc[y],np.nan)
    x=pd.Series(pc.PC1.values,index=frame.index,name='PC1 (28.9%)')
    table(r,f'PCA PC1 vs {y}',frame,kind='XY',x=x);csv(pd.concat([x,frame],axis=1),f'01_PCA_{y}')
save(r,'01_PCA_data');spec.append({'project':'01_PCA_data.pzfx','graph_type':'XY scatter, no connecting lines','table_count':2,'rows_per_table':21,'x_title':'PC1 (28.9%)','y_titles':['PC2 (22.6%)','PC3 (12.1%)'],'palette':D.color.to_dict()})

z=pd.read_csv(T/'02_heatmap_HVG100_gene_zscores.csv').set_index('gene_symbol');cols=pd.read_csv(T/'02_heatmap_HVG100_order.csv').cell_type.tolist()
# Restore original Ward row order from saved matrix; float precision checked independently.
order=leaves_list(linkage(z.to_numpy(),method='ward',metric='euclidean',optimal_ordering=True));z=z.iloc[order][cols];z.columns=[short[c] for c in cols]
r=root();E.SubElement(r,'DefGraphButton',ButtonID='367');table(r,'HVG100 expression patterns',z);save(r,'02_HVG100_heatmap_data');csv(z,'02_HVG100_heatmap')
spec.append({'project':'02_HVG100_heatmap_data.pzfx','graph_type':'Grouped heatmap','rows':100,'columns':21,'color_limits':[-2.5,2.5],'note':'Data retain full z-scores. Display clipping only. Original Ward order; no native dendrogram assumed.'})

pred=pd.read_csv(T/'04_raw_oof_predictions.csv');r=root()
for j,t in enumerate(focus):
    f=pred[pred.target==t].sort_values('gene_index');x=f.observed.reset_index(drop=True);x.name='Observed expression difference'
    frame=pd.DataFrame({'Lasso prediction':f.Lasso_min.to_numpy()},index=f.gene_symbol);x.index=frame.index
    table(r,short[t],frame,kind='XY',x=x);csv(pd.concat([x,frame],axis=1),f'03_OOF_{j+1}')
save(r,'03_OOF_data');spec.append({'project':'03_OOF_data.pzfx','graph_type':'XY scatter, no lines; small points','table_count':4,'rows_per_table':11922,'x_title':'Observed expression difference','y_title':'Predicted expression difference','layout':'2 x 2','note':'All genes retained. Add y=x if supported. Original hexbin changed to native scatter.'})

c=pd.read_csv(T/'04_top10_TF_candidates.csv');r=root()
for j,t in enumerate(focus):
    f=c[c.target==t].set_index('perturbation');table(r,short[t],f,interval=True);csv(f[['median_coefficient','coefficient_q975','coefficient_q025','selection_frequency']],f'04_TF_{j+1}')
save(r,'04_TF_stability_data');spec.append({'project':'04_TF_stability_data.pzfx','graph_type':'Grouped horizontal interval/forest, linear coefficient axis','table_count':4,'rows_per_table':10,'interval_order':['median','upper absolute limit','lower absolute limit'],'layout':'2 x 2','note':'All displayed TFs selected in 100/100 subsamples. Lines are empirical 2.5-97.5 percentiles, not confidence intervals.'})

m=pd.read_csv(T/'04_model_metrics.csv');r=root();models=['Mean_baseline','PCR','PLS','Ridge','NonnegativeLasso','Lasso','ElasticNet']
for t in focus:
    f=m[(m.target==t)&(m.dataset=='raw')].pivot(index='model',columns='selection',values='R2').reindex(models)
    f=f.rename(index={'Mean_baseline':'Mean baseline','NonnegativeLasso':'Nonnegative Lasso','ElasticNet':'Elastic Net'},columns={'min':'Minimum MSE','1se':'1-SE'})
    table(r,short[t],f);csv(f,'05_Models_'+str(focus.index(t)+1))
save(r,'05_All_models_data');spec.append({'project':'05_All_models_data.pzfx','graph_type':'Grouped points, no error bars','table_count':4,'layout':'2 x 2','y_title':'Out-of-fold R2','note':'Complete six models plus baseline; each scalar pooled OOF, no replicate error bars.'})

stable=pd.read_csv(T/'04_stability_all_TFs.csv');labels=list(dict.fromkeys(c.perturbation));f=pd.DataFrame(index=labels)
for t in focus:
    st=stable[stable.target==t].set_index('perturbation').loc[labels];f[short[t]]=st.selection_frequency*np.sign(st.full_coefficient)
r=root();E.SubElement(r,'DefGraphButton',ButtonID='367');table(r,'Shared candidate patterns',f);save(r,'06_Cross_target_heatmap_data');csv(f,'06_Cross_target_heatmap');spec.append({'project':'06_Cross_target_heatmap_data.pzfx','graph_type':'Grouped heatmap','rows':len(f),'columns':4,'color_limits':[-1,1],'note':'Signed selection frequency; original white stable-priority marker documented in separate source table.'})

(P/'figure_specs.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n')
(P/'README.md').write_text('# Prism 图形迁移输入\n\n全部数值来自已审核的保存结果。脚本 `scripts/prepare_prism_migration.py` 生成可编辑数据表，尚未生成原生图形。正式图形由 Prism 建立、保存和导出。\n\n1. PCA：两张 XY 图，PC1 对 PC2／PC3。\n2. HVG100：100×21 热图，保留 Ward 行列顺序，色阶 ±2.5。\n3. OOF：四目标各 11,922 个点，2×2 布局。\n4. TF：各 10 个候选，中心为中位数，绝对上下限为 97.5%／2.5% 分位数，线性坐标。\n5. 多模型：完整六模型及基线，最低误差／1-SE 两种选择。\n6. 跨目标候选：带方向的选择频率热图。\n\n字体 Arial，图内英文，中文图注见报告。详情与配色见 figure_specs.json。\n')
files=[]
for p in sorted(P.iterdir()):
    if p.is_file():files.append({'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(R/'logs/stage5_prism_input_manifest.json').write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),'new_model_fits':0,'files':files},ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'projects':len(spec),'tables':sum(len(E.parse(P/s['project']).getroot().findall('Table')) for s in spec),'output':str(P)}))
