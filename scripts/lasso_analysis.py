"""Stage 3: shared nested gene CV, fold-local Lasso paths, and recoverable outputs."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
from sklearn.linear_model import lasso_path
from sklearn.model_selection import KFold
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / '.cache/lasso'
CONFIG = json.loads((ROOT / 'config/analysis_config.json').read_text())
RATIOS = np.geomspace(CONFIG['lasso_relative_alpha']['maximum'], CONFIG['lasso_relative_alpha']['minimum'], CONFIG['lasso_relative_alpha']['count'])
THRESHOLD = CONFIG['stability']['nonzero_threshold_standardized']


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(2**20), b''):
            h.update(chunk)
    return h.hexdigest()


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')
    temp.replace(path)


def csv(path, frame):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, float_format='%.14g')


def load_data():
    with h5py.File(ROOT / 'data/processed/analysis_matrices.h5', 'r') as f:
        return (f['model/X'][:].astype('float64'), f['model/Y'][:],
                f['model/genes'].asstr()[:], f['model/perturbations'].asstr()[:],
                f['model/target_cells'].asstr()[:])


def prepare():
    """Only training rows determine cached moments, Gram matrices and X'y."""
    start = time.monotonic()
    CACHE.mkdir(parents=True, exist_ok=True)
    X, Y, genes, perturbations, targets = load_data()
    np.save(CACHE / 'X.npy', X)
    np.save(CACHE / 'Y.npy', Y)
    n = len(genes)
    indices = np.arange(n)
    split = {'seed': CONFIG['seed'], 'gene_order_sha256': hashlib.sha256('\n'.join(genes).encode()).hexdigest(), 'outer': [], 'full_inner': []}
    contexts = {}
    assignments = pd.DataFrame({'gene_index': indices, 'gene_symbol': genes, 'outer_test_fold': -1, 'full_inner_validation_fold': -1})
    outer_cv = KFold(n_splits=CONFIG['outer_folds'], shuffle=True, random_state=CONFIG['seed'])
    for k, (train, test) in enumerate(outer_cv.split(indices), 1):
        assignments.loc[test, 'outer_test_fold'] = k
        assignments[f'outer_{k:02d}_inner_validation_fold'] = -1
        entry = {'fold': k, 'train': train.tolist(), 'test': test.tolist(), 'inner': []}
        contexts[f'outer_{k:02d}_fit'] = train
        for j, (itr, iva) in enumerate(KFold(CONFIG['inner_folds'], shuffle=True, random_state=CONFIG['seed'] + k).split(train), 1):
            tr, va = train[itr], train[iva]
            entry['inner'].append({'fold': j, 'train': tr.tolist(), 'validation': va.tolist()})
            contexts[f'outer_{k:02d}_inner_{j:02d}'] = tr
            assignments.loc[va, f'outer_{k:02d}_inner_validation_fold'] = j
        split['outer'].append(entry)
    for j, (tr, va) in enumerate(KFold(CONFIG['inner_folds'], shuffle=True, random_state=CONFIG['seed']).split(indices), 1):
        split['full_inner'].append({'fold': j, 'train': tr.tolist(), 'validation': va.tolist()})
        contexts[f'full_inner_{j:02d}'] = tr
        assignments.loc[va, 'full_inner_validation_fold'] = j
    contexts['full_fit'] = indices
    dump(ROOT / 'data/splits/gene_cv_splits.json', split)
    csv(ROOT / 'data/splits/gene_cv_assignments.csv', assignments)
    context_records = []
    for name, tr in contexts.items():
        t0 = time.monotonic()
        xx = np.array(X[tr], order='F')
        xm = xx.mean(axis=0)
        xs = xx.std(axis=0, ddof=0)
        constant = xs == 0
        xs[constant] = 1
        xx -= xm
        xx /= xs
        yy = Y[tr]
        ym = yy.mean(axis=0)
        ys = yy.std(axis=0, ddof=0)
        assert np.all(ys > 0)
        yz = (yy - ym) / ys
        xy = xx.T @ yz
        if CONFIG['solver'].get('lasso_engine','sklearn_gram') in ['sklearn_gram','gpu_fista']:
            gram = np.ascontiguousarray(xx.T @ xx)
            np.save(CACHE / f'{name}_gram.npy', gram)
        np.savez(CACHE / f'{name}_moments.npz', train=tr, x_mean=xm, x_scale=xs, y_mean=ym, y_scale=ys, xy=xy, constant=constant)
        rec = {'context': name, 'n_train': len(tr), 'n_constant_features': int(constant.sum()), 'seconds': time.monotonic()-t0}
        context_records.append(rec)
        print(json.dumps({'event': 'prepared', **rec}), flush=True)
    dump(CACHE / 'metadata.json', {'input_sha256': sha(ROOT / 'data/processed/analysis_matrices.h5'), 'targets': targets.tolist(), 'genes': genes.tolist(), 'perturbations': perturbations.tolist(), 'n_samples': n, 'n_features': X.shape[1], 'seed': CONFIG['seed']})
    csv(ROOT / 'tables/03_preprocessing_contexts.csv', pd.DataFrame(context_records))
    dump(ROOT / 'logs/stage3_preparation.json', {'generated_utc': now(), 'seconds': time.monotonic()-start, 'n_contexts': len(contexts), 'all_preprocessing_fitted_on_context_train_only': True, 'gram_cache_role': 'derived_disposable_acceleration_cache'})


def fit_context(X, y, target_index, name, ratios=RATIOS):
    start = time.monotonic()
    with np.load(CACHE / f'{name}_moments.npz') as z:
        m = {key: z[key] for key in z.files}
    tr = m['train']
    engine = CONFIG['solver'].get('lasso_engine', 'sklearn_gram')
    if engine == 'gpu_fista':
        with np.load(ROOT / 'results/lasso_gpu_paths' / f'{name}.npz') as z:
            idx=[int(np.flatnonzero(np.isclose(z['ratios'],r,rtol=1e-12,atol=1e-15))[0]) for r in ratios]
            coef=z['coef'][:,target_index,:][:,idx]
            gap=z['gap'][target_index,idx];iters=z['n_iter'][target_index,idx]
            amax=float(z['alpha_max'][target_index]);threshold=float(z['dual_tolerance'][target_index])
        diagnostics=[{'context':name,'solver':engine,'iteration_unit':'proximal_gradient_steps',
            'relative_alpha':float(r),'alpha':float(r*amax),'alpha_max':amax,
            'dual_gap':float(gap[j]),'dual_gap_tolerance':threshold,'n_iter':int(iters[j]),
            'retried':bool(iters[j]>CONFIG['solver']['max_iter']),
            'converged':bool(-1e-8<=gap[j]<=threshold*1.00001),
            'n_nonzero':int(np.count_nonzero(coef[:,j])),
            'n_active_gt_1e_6':int(np.sum(np.abs(coef[:,j])>THRESHOLD))} for j,r in enumerate(ratios)]
        assert all(d['converged'] for d in diagnostics)
        return {'coef':coef,'x_mean':m['x_mean'],'x_scale':m['x_scale'],
                'y_mean':float(m['y_mean'][target_index]),'y_scale':float(m['y_scale'][target_index]),
                'alpha_max':amax,'diagnostics':diagnostics,'seconds':time.monotonic()-start,
                'warning_messages':[]}
    xx = np.asarray((X[tr] - m['x_mean']) / m['x_scale'], order='F')
    ym, ys = float(m['y_mean'][target_index]), float(m['y_scale'][target_index])
    yy = np.ascontiguousarray((y[tr]-ym)/ys)
    xy = np.ascontiguousarray(m['xy'][:, target_index])
    gram = np.load(CACHE / f'{name}_gram.npy', mmap_mode='r') if engine == 'sklearn_gram' else None
    alpha_max = float(np.max(np.abs(xy)) / len(tr))
    assert alpha_max > 0
    alphas = alpha_max * np.asarray(ratios)
    solver = CONFIG['solver']
    def solve(penalties, limit, initial=None):
        if engine == 'celer':
            from celer import celer_path
            return celer_path(xx, yy, pb='lasso', alphas=penalties, coef_init=initial,
                    tol=solver['tol'], max_iter=solver['celer_max_working_set_iterations'],
                    max_epochs=limit, return_n_iter=True)
        return lasso_path(xx, yy, alphas=penalties, precompute=gram, Xy=xy, coef_init=initial,
                    copy_X=False, check_input=False, tol=solver['tol'], max_iter=limit,
                    return_n_iter=True, selection='cyclic')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', ConvergenceWarning)
        _, coef, gap, iters = solve(alphas, solver['max_iter'])
    threshold = solver['tol'] * np.dot(yy, yy) / len(yy)
    retry = np.zeros(len(alphas), dtype=bool)
    iters = np.asarray(iters)
    for j in np.flatnonzero((gap > threshold * 1.00001) | ~np.isfinite(gap)):
        retry[j] = True
        with warnings.catch_warnings(record=True) as extra:
            warnings.simplefilter('always', ConvergenceWarning)
            _, c, g, ni = solve([alphas[j]], solver['retry_max_iter'], coef[:, j].copy())
        caught.extend(extra)
        coef[:, j], gap[j], iters[j] = c[:, 0], g[0], ni[0]
    converged = np.isfinite(gap) & (gap <= threshold * 1.00001) & (gap >= -1e-8)
    diagnostics = [{'context': name, 'solver': engine,
                    'iteration_unit': 'working_set_updates' if engine == 'celer' else 'coordinate_sweeps',
                    'relative_alpha': float(ratios[j]), 'alpha': float(alphas[j]),
                    'alpha_max': alpha_max, 'dual_gap': float(gap[j]), 'dual_gap_tolerance': float(threshold),
                    'n_iter': int(iters[j]), 'retried': bool(retry[j]), 'converged': bool(converged[j]),
                    'n_nonzero': int(np.count_nonzero(coef[:, j])),
                    'n_active_gt_1e_6': int(np.sum(np.abs(coef[:, j]) > THRESHOLD))}
                   for j in range(len(alphas))]
    if not np.all(converged):
        raise RuntimeError(f'{name}: nonconverged alphas {alphas[~converged]}')
    print(json.dumps({'event':'path_done','target_index':target_index,'context':name,
                      'solver':engine,'seconds':time.monotonic()-start,'retried':int(retry.sum()),
                      'warning_count':len(caught),'max_dual_gap':float(gap.max())}),flush=True)
    return {'coef': coef, 'x_mean': m['x_mean'], 'x_scale': m['x_scale'], 'y_mean': ym, 'y_scale': ys,
            'alpha_max': alpha_max, 'diagnostics': diagnostics, 'seconds': time.monotonic()-start,
            'warning_messages': [str(w.message) for w in caught]}


def predict(X, fit):
    # Exact inverse scaling also supplies portable raw-input model coefficients.
    beta = fit['coef'] * fit['y_scale'] / fit['x_scale'][:, None]
    intercept = fit['y_mean'] - fit['x_mean'] @ beta
    return X @ beta + intercept


def tune(X, y, target_index, inner, prefix):
    errors, diagnostics, times = [], [], []
    for item in inner:
        fit = fit_context(X, y, target_index, f"{prefix}_inner_{item['fold']:02d}")
        va = np.array(item['validation'])
        prediction = predict(X[va], fit)
        errors.append(np.mean((y[va, None] - prediction)**2, axis=0))
        diagnostics.extend(fit['diagnostics'])
        times.append(fit['seconds'])
    errors = np.array(errors)
    mean = errors.mean(axis=0)
    se = errors.std(axis=0, ddof=1) / np.sqrt(errors.shape[0])
    imin = int(np.argmin(mean))
    i1se = int(np.flatnonzero(mean <= mean[imin] + se[imin])[0])
    return {'errors': errors, 'mean': mean, 'se': se, 'min_index': imin, 'one_se_index': i1se,
            'diagnostics': diagnostics, 'seconds': float(sum(times))}


def metrics(y, pred):
    residual = y - pred
    r2 = 1 - np.sum(residual**2) / np.sum((y - y.mean())**2)
    variable = np.ptp(pred) > 1e-14
    return {'R2': float(r2), 'RMSE': float(np.sqrt(np.mean(residual**2))),
            'Pearson': float(np.corrcoef(y, pred)[0,1]) if variable else np.nan,
            'Spearman': float(spearmanr(y, pred).statistic) if variable else np.nan}


def run_target(target_index, threads=1):
    with threadpool_limits(limits=threads):
        return _run_target(target_index)


def _run_target(target_index):
    start = time.monotonic()
    meta = json.loads((CACHE / 'metadata.json').read_text())
    target = meta['targets'][target_index]
    slug = f'target_{target_index:02d}'
    out = ROOT / 'results/lasso' / slug
    out.mkdir(parents=True, exist_ok=True)
    X = np.load(CACHE / 'X.npy', mmap_mode='r')
    y = np.load(CACHE / 'Y.npy', mmap_mode='r')[:, target_index]
    split = json.loads((ROOT / 'data/splits/gene_cv_splits.json').read_text())
    oof = np.full((len(y), 3), np.nan)
    outer_fold = np.zeros(len(y), dtype=int)
    diagnostics, tuning_rows, selected, fold_metrics, raw_betas, intercepts, std_betas = [], [], [], [], [], [], []
    def add_tuning(cv, label):
        diagnostics.extend(cv['diagnostics'])
        for j, ratio in enumerate(RATIOS):
            tuning_rows.append({'scope': label, 'relative_alpha': ratio, 'mean_MSE': cv['mean'][j], 'SE_MSE': cv['se'][j],
                                **{f'inner_{k+1}_MSE': cv['errors'][k,j] for k in range(CONFIG['inner_folds'])},
                                'selected_min': j == cv['min_index'], 'selected_1se': j == cv['one_se_index']})
    for item in split['outer']:
        k = item['fold']; label = f'outer_{k:02d}'
        cv = tune(X, y, target_index, item['inner'], label)
        add_tuning(cv, label)
        chosen = sorted(set([cv['one_se_index'], cv['min_index']]))
        fit = fit_context(X, y, target_index, f'{label}_fit', RATIOS[chosen])
        diagnostics.extend(fit['diagnostics'])
        mapping = [chosen.index(cv['min_index']), chosen.index(cv['one_se_index'])]
        te = np.array(item['test']); tr = np.array(item['train'])
        oof[te,:2] = predict(X[te], fit)[:,mapping]
        oof[te,2] = y[tr].mean()
        outer_fold[te] = k
        raw = fit['coef'][:,mapping] * fit['y_scale'] / fit['x_scale'][:,None]
        raw_betas.append(raw);std_betas.append(fit['coef'][:,mapping])
        intercepts.append(fit['y_mean'] - fit['x_mean'] @ raw)
        for label_model, idx, col in [('Lasso_min',cv['min_index'],0),('Lasso_1se',cv['one_se_index'],1)]:
            selected.append({'scope':label,'model':label_model,'relative_alpha':RATIOS[idx],'alpha':RATIOS[idx]*fit['alpha_max'],
                             'n_active_gt_1e_6':int(np.sum(np.abs(fit['coef'][:,mapping[col]])>THRESHOLD)),
                             'selection_at_lowest_alpha':bool(idx == len(RATIOS)-1)})
        for col, model in enumerate(['Lasso_min','Lasso_1se','Mean_baseline']):
            fold_metrics.append({'fold':k,'model':model,'n_test':len(te), **metrics(y[te],oof[te,col])})
        print(json.dumps({'event':'outer_done','target':target,'fold':k,'seconds_elapsed':time.monotonic()-start}),flush=True)
    cv = tune(X,y,target_index,split['full_inner'],'full')
    add_tuning(cv,'full')
    full = fit_context(X,y,target_index,'full_fit')
    diagnostics.extend(full['diagnostics'])
    indices = [cv['min_index'],cv['one_se_index']]
    coef = full['coef'][:,indices]
    raw = coef * full['y_scale'] / full['x_scale'][:,None]
    intercept = full['y_mean'] - full['x_mean'] @ raw
    for j, model in enumerate(['Lasso_min','Lasso_1se']):
        selected.append({'scope':'full','model':model,'relative_alpha':RATIOS[indices[j]],'alpha':RATIOS[indices[j]]*full['alpha_max'],
                         'n_active_gt_1e_6':int(np.sum(np.abs(coef[:,j])>THRESHOLD)),
                         'selection_at_lowest_alpha':bool(indices[j] == len(RATIOS)-1)})
    np.savez_compressed(out/'model_arrays.npz',target=target,perturbations=np.asarray(meta['perturbations']),
         full_coef_standardized=coef,full_coef_input=raw,full_intercept_input=intercept,
         full_x_mean=full['x_mean'],full_x_scale=full['x_scale'],full_y_mean=full['y_mean'],full_y_scale=full['y_scale'],
         full_path_standardized=full['coef'],relative_alphas=RATIOS,full_alpha_max=full['alpha_max'],
         selected_indices=np.asarray(indices),outer_coef_input=np.stack(raw_betas),outer_intercept_input=np.stack(intercepts),
         outer_coef_standardized=np.stack(std_betas))
    pred_frame=pd.DataFrame({'gene_symbol':meta['genes'],'outer_fold':outer_fold,'observed':y,
                            'Lasso_min':oof[:,0],'Lasso_1se':oof[:,1],'Mean_baseline':oof[:,2]})
    csv(out/'oof_predictions.csv',pred_frame)
    csv(out/'tuning_curves.csv',pd.DataFrame(tuning_rows))
    csv(out/'solver_diagnostics.csv',pd.DataFrame(diagnostics))
    csv(out/'selected_models.csv',pd.DataFrame(selected))
    csv(out/'fold_metrics.csv',pd.DataFrame(fold_metrics))
    metric_rows=[{'target':target,'model':model,**metrics(y,oof[:,j])} for j,model in enumerate(['Lasso_min','Lasso_1se','Mean_baseline'])]
    csv(out/'pooled_metrics.csv',pd.DataFrame(metric_rows))
    result={'target_index':target_index,'target':target,'completed_utc':now(),'seconds':time.monotonic()-start,
            'timing_scope':'assembly_after_shared_GPU_paths' if CONFIG['solver'].get('lasso_engine')=='gpu_fista' else 'target_fit_and_assembly',
            'metrics':metric_rows,'n_fits':len(diagnostics),'n_retried':sum(d['retried'] for d in diagnostics),
            'all_converged':all(d['converged'] for d in diagnostics),'input_sha256':meta['input_sha256']}
    dump(out/'completed.json',result)
    print(json.dumps({'event':'target_done',**result},ensure_ascii=False),flush=True)
    return result


def single_correlations():
    X,Y,genes,perturbations,targets=load_data()
    def standard(a):
        z=a-a.mean(axis=0)
        return z/np.sqrt(np.sum(z*z,axis=0))
    pearson=standard(X).T @ standard(Y)
    spearman=standard(rankdata(X,axis=0)).T @ standard(rankdata(Y,axis=0))
    annot=pd.read_csv(ROOT/'tables/01_perturbation_annotations.csv')
    records=[]
    for k,t in enumerate(targets):
        f=pd.DataFrame({'target':t,'perturbation':perturbations,'Pearson':pearson[:,k],'Spearman':spearman[:,k]})
        f['rank_abs_Pearson']=f.Pearson.abs().rank(method='first',ascending=False).astype(int)
        f['rank_positive_Pearson']=f.Pearson.rank(method='first',ascending=False).astype(int)
        records.append(f)
    csv(ROOT/'tables/03_single_perturbation_correlations.csv',pd.concat(records,ignore_index=True))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--prepare',action='store_true');p.add_argument('--targets',nargs='*',type=int)
    p.add_argument('--workers',type=int,default=1);p.add_argument('--threads',type=int,default=1)
    p.add_argument('--correlations',action='store_true')
    args=p.parse_args()
    with threadpool_limits(limits=args.threads):
        if args.prepare:prepare()
        if args.correlations:single_correlations()
    if args.targets is not None:
        todo=[i for i in args.targets if not (ROOT/f'results/lasso/target_{i:02d}/completed.json').exists()]
        if args.workers==1:
            for i in todo:run_target(i,args.threads)
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                jobs={pool.submit(run_target,i,args.threads):i for i in todo}
                for job in as_completed(jobs):job.result()
