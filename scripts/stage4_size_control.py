"""Exploratory 965-gene non-landmark controls with matched outer-fold counts."""
from pathlib import Path
import argparse,json
import h5py
import numpy as np
import pandas as pd
from stage4_assemble import select,metrics,frame_save,FOCUS,load_context
ROOT=Path(__file__).resolve().parents[1]

def prepare():
    with np.load(ROOT/'data/processed/stage4_input_metadata.npz') as z:
        lm=set(z['landmark_indices'].tolist());n=len(z['genes'])
    split=json.loads((ROOT/'data/splits/gene_cv_splits.json').read_text())
    base=json.loads((ROOT/'data/splits/stage4_jobs.json').read_text())['jobs']
    base=[j for j in base if j['dataset']=='raw' and j['kind']=='context' and j['name'].startswith('outer_')]
    alljobs=[];samples={}
    for rep in range(1,4):
        seed=20261107+4100+rep;rng=np.random.default_rng(seed);chosen=[]
        for e in split['outer']:
            n_lm=sum(g in lm for g in e['test']);pool=[g for g in e['test'] if g not in lm]
            chosen.extend(rng.choice(pool,n_lm,replace=False).tolist())
        chosen=sorted(chosen);keep=set(chosen);dataset=f'size_matched_{rep:02d}'
        samples[dataset]={'seed':seed,'gene_indices':chosen}
        for c in base:
            job=dict(c);job.update(id=dataset+'/'+c['name'],dataset=dataset,
                train=[g for g in c['train'] if g in keep],validation=[g for g in c['validation'] if g in keep])
            alljobs.append(job)
    for name,obj in [('stage4_size_control_jobs.json',{'jobs':alljobs}),('stage4_size_control_samples.json',samples)]:
        (ROOT/'data/splits'/name).write_text(json.dumps(obj,indent=2)+'\n')
    cfgpath=ROOT/'config/stage4_config.json';cfg=json.loads(cfgpath.read_text())
    cfg['exploratory_size_control']={'repetitions':3,'n_genes':965,'sampling':'non-landmark genes; exact landmark count in each outer test fold',
        'seeds':[20265208,20265209,20265210],'model':'Lasso, full 5x3 nested tuning','scope':'training-set size/composition comparison; not part of prespecified model ranking'}
    cfgpath.write_text(json.dumps(cfg,ensure_ascii=False,indent=2)+'\n')
    print({'jobs':len(alljobs),'sample_sizes':{k:len(v['gene_indices']) for k,v in samples.items()}})

def assemble():
    with h5py.File(ROOT/'data/processed/analysis_matrices.h5') as h:
        X=h['model/X'][:].astype(float);t=h['model/target_cells'].asstr()[:].tolist();Y=h['model/Y'][:][:,[t.index(x) for x in FOCUS]]
        genes=h['model/genes'].asstr()[:].astype(str)
    split=json.loads((ROOT/'data/splits/gene_cv_splits.json').read_text());samples=json.loads((ROOT/'data/splits/stage4_size_control_samples.json').read_text())
    metric=[];predictions=[];choices=[];curves=[]
    for dataset,sample in samples.items():
        ids=np.array(sample['gene_indices']);pos={g:i for i,g in enumerate(ids)};pred=np.zeros((len(ids),4,2));base=np.zeros((len(ids),4));fold=np.zeros(len(ids),int)
        for e in split['outer']:
            scope=f'outer_{e["fold"]:02d}';fit,meta=load_context(dataset,scope+'_fit');tr=fit['train']
            inners=[load_context(dataset,scope+f'_inner_{j:02d}')[0] for j in range(1,4)];pars=meta['parameters']['Lasso']
            te=np.array([g for g in e['test'] if g in pos]);loc=[pos[g] for g in te];fold[loc]=e['fold'];base[loc]=Y[tr].mean(axis=0)
            for t,target in enumerate(FOCUS):
                mse=np.stack([a['Lasso_validation_MSE'][t] for a in inners]);selected,means,se=select(mse,pars,'Lasso')
                for j,param in enumerate(pars):curves.append({'dataset':dataset,'target':target,'scope':scope,'parameter_index':j,
                    **param,'mean_MSE':means[j],'SE_MSE':se[j],**{f'fold_{k+1}_MSE':mse[k,j] for k in range(3)}})
                for k,(choice,j) in enumerate(zip(['min','1se'],selected)):
                    b=fit['Lasso_coef'][:,t,j];pred[loc,t,k]=((X[te]-fit['x_mean'])/fit['x_scale'])@b*fit['y_scale'][t]+fit['y_mean'][t]
                    choices.append({'dataset':dataset,'target':target,'scope':scope,'selection':choice,'parameter_index':j,**pars[j]})
        for t,target in enumerate(FOCUS):
            for k,choice in enumerate(['min','1se']):metric.append({'dataset':dataset,'target':target,'selection':choice,'model':'Lasso',**metrics(Y[ids,t],pred[:,t,k])})
            metric.append({'dataset':dataset,'target':target,'selection':'baseline','model':'Mean_baseline',**metrics(Y[ids,t],base[:,t])})
        f=pd.DataFrame({'dataset':dataset,'gene_index':np.repeat(ids,4),'gene_symbol':np.repeat(genes[ids],4),
            'target':np.tile(FOCUS,len(ids)),'outer_fold':np.repeat(fold,4),'observed':Y[ids].ravel(),
            'Lasso_min':pred[:,:,0].ravel(),'Lasso_1se':pred[:,:,1].ravel(),'Mean_baseline':base.ravel()})
        predictions.append(f)
    frame_save('04_size_control_metrics.csv',pd.DataFrame(metric));frame_save('04_size_control_predictions.csv',pd.concat(predictions,ignore_index=True))
    frame_save('04_size_control_choices.csv',pd.DataFrame(choices));frame_save('04_size_control_tuning.csv',pd.DataFrame(curves))
    print(pd.DataFrame(metric).query('selection=="min"').pivot(index='dataset',columns='target',values='R2').to_string())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');args=p.parse_args()
    prepare() if args.prepare else assemble()
