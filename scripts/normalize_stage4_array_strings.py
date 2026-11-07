"""Make archived HDF5-derived string labels readable without object/pickle arrays."""
from pathlib import Path
import h5py
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    with h5py.File(ROOT/'data/processed/analysis_matrices.h5') as f:
        genes=f['model/genes'].asstr()[:].astype(str);pert=f['model/perturbations'].asstr()[:].astype(str)
    paths=[ROOT/'data/processed/stage4_input_metadata.npz']+list((ROOT/'results/stage4/final').glob('*.npz'))
    for p in paths:
        with np.load(p) as z:
            a={k:z[k] for k in z.files if k not in ['genes','perturbations']}
            if 'genes' in z.files:a['genes']=genes[a['gene_indices']] if 'gene_indices' in a else genes
            if 'perturbations' in z.files:a['perturbations']=pert
        temp=p.with_suffix('.npz.tmp')
        with temp.open('wb') as f:np.savez_compressed(f,**a)
        temp.replace(p)
        with np.load(p) as z:assert all(z[k].dtype.kind!='O' for k in z.files)
    print({'normalized_archives':len(paths)})
