"""Freeze stage-4 inputs and jobs while preserving the approved stage-3 outputs."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')

def main():
    config=json.loads((ROOT/'config/analysis_config.json').read_text())
    with h5py.File(ROOT/'data/processed/analysis_matrices.h5') as f:
        X=f['model/X'][:].astype('float64'); allY=f['model/Y'][:]
        targets=f['model/target_cells'].asstr()[:].tolist()
        focus=[targets.index(t) for t in config['focus_targets']]
        genes=f['model/genes'].asstr()[:].astype(str);pert=f['model/perturbations'].asstr()[:].astype(str)
        landmark=np.flatnonzero(f['model/is_landmark'][:])
        blood=f['blood/expression'][:];bg=f['blood/genes'].asstr()[:].tolist()
        cells=f['blood/cell_types'].asstr()[:].tolist()
        included=[i for i,c in enumerate(cells) if c!=config['supplementary_target']]
        scale=blood[:,included].std(axis=1,ddof=1)[[bg.index(g) for g in genes]]
    assert np.all(scale>0)
    cache=ROOT/'.cache/stage4';cache.mkdir(parents=True,exist_ok=True)
    np.save(cache/'X.npy',X)
    np.save(cache/'Y_raw.npy',allY[:,focus]);np.save(cache/'Y_gene_scaled.npy',allY[:,focus]/scale[:,None])
    np.savez_compressed(ROOT/'data/processed/stage4_input_metadata.npz',genes=genes,perturbations=pert,
        focus_targets=np.array(config['focus_targets']),focus_indices=focus,landmark_indices=landmark,
        gene_response_scale=scale,scale_cell_types=np.array([cells[i] for i in included]))
    split=json.loads((ROOT/'data/splits/gene_cv_splits.json').read_text())
    contexts=[]
    for entry in split['outer']:
        k=entry['fold'];contexts.append({'name':f'outer_{k:02d}_fit','train':entry['train'],'validation':[]})
        for inner in entry['inner']:
            contexts.append({'name':f'outer_{k:02d}_inner_{inner["fold"]:02d}',
                'train':inner['train'],'validation':inner['validation']})
    contexts.append({'name':'full_fit','train':list(range(len(X))),'validation':[]})
    for inner in split['full_inner']:
        contexts.append({'name':f'full_inner_{inner["fold"]:02d}','train':inner['train'],'validation':inner['validation']})
    lm=set(landmark.tolist());jobs=[]
    for dataset in ['raw','landmark','gene_scaled']:
        for context in contexts:
            job={'id':dataset+'/'+context['name'],'kind':'context','dataset':dataset,**context}
            if dataset=='landmark':
                job['train']=[i for i in job['train'] if i in lm]
                job['validation']=[i for i in job['validation'] if i in lm]
            jobs.append(job)
    ratios=[];full=[]
    for i in focus:
        with np.load(ROOT/f'results/lasso/target_{i:02d}/model_arrays.npz') as z:
            ratios.append(float(z['relative_alphas'][z['selected_indices'][1]]))
            full.append(z['full_coef_standardized'][:,1])
    np.savez_compressed(ROOT/'data/processed/stage4_stability_reference.npz',relative_alphas=ratios,
                        full_coefficients=np.stack(full,axis=1))
    rng=np.random.default_rng(config['seed']+4000)
    samples=np.stack([np.sort(rng.choice(len(X),int(np.floor(.8*len(X))),replace=False)) for _ in range(100)])
    np.savez_compressed(ROOT/'data/splits/stage4_stability_indices.npz',train_indices=samples,
                        seed=config['seed']+4000)
    for i,tr in enumerate(samples):
        jobs.append({'id':f'stability/repeat_{i+1:03d}','name':f'repeat_{i+1:03d}',
            'kind':'stability','dataset':'raw','train':tr.tolist(),'validation':[]})
    design={'created_utc':datetime.now(timezone.utc).isoformat(),'focus_targets':config['focus_targets'],
      'input_sha256':hashlib.sha256((ROOT/'data/processed/analysis_matrices.h5').read_bytes()).hexdigest(),
      'n_genes':len(X),'n_perturbations':X.shape[1],'landmark_genes':len(landmark),
      'raw_models':['ElasticNet','NonnegativeLasso','Ridge','PCR','PLS'],
      'reference_models':['Lasso','Mean_baseline'],'relative_path':config['lasso_relative_alpha'],
      'elastic_net_l1_ratios':config['elastic_net_l1_ratio'],'ridge':config['ridge_mean_loss_lambda'],
      'pcr_components':config['pcr_components'],'pls_components':config['pls_components'],
      'tol':config['solver']['tol'],'max_iter':config['solver']['retry_max_iter'],
      'active_threshold':config['stability']['nonzero_threshold_standardized'],
      'one_se_order':{'ElasticNet':'largest relative penalty, then largest l1_ratio',
        'NonnegativeLasso':'largest relative penalty','Lasso':'largest relative penalty',
        'Ridge':'largest mean-loss lambda','PCR':'fewest components','PLS':'fewest components'},
      'stability':{'repetitions':100,'fraction':.8,'seed':config['seed']+4000,
        'initialization':'zero coefficients, independently for every subsample',
        'selection':'full-data Lasso 1-SE relative penalty held fixed; alpha_max and scaling recomputed per subsample',
        'frequency_threshold':.8,'sign_consistency_threshold':.9,
        'sign_consistency':'majority sign among active resamples',
        'median_absolute_coefficient':'all 100 resamples, including inactive values as zero',
        'ranking':'selection frequency, median absolute standardized coefficient, full coefficient magnitude'},
      'gene_scaled_response':'target-HSC divided by gene sample SD across 21 blood types excluding mDC',
      'landmark_splits':'intersection with the saved outer and inner gene lists',
      'self_readout':'remove matched own-gene rows, refit fixed relative Lasso 1-SE, inspect coefficient/sign and marginal correlation',
      'correlation_groups':{'primary_abs_Pearson':.8,'sensitivity':[.7,.9],'definition':'connected components of thresholded perturbation correlations'},
      'extended_retuning':'20 resamples per target after stage 5, as planned'}
    dump(ROOT/'config/stage4_config.json',design);dump(ROOT/'data/splits/stage4_jobs.json',{'jobs':jobs})
    config['execution']['current_authorized_stage']=4;dump(ROOT/'config/analysis_config.json',config)
    dump(ROOT/'stage_status.json',{'stage':4,'status':'in_progress','started_utc':design['created_utc'],
        'previous_stage_user_approved':True,'independent_review_required':True,'next_stage':5,
        'next_stage_started':False,'pending_items':[{'item':'Prism native sample','status':'user_approved_deferral'}]})
    print(json.dumps({'jobs':len(jobs),'shape':X.shape,'landmark':len(landmark),'focus':design['focus_targets']}))

if __name__=='__main__':main()
