import json
import numpy as np
import anndata as ad
import scipy.sparse as sp
from pathlib import Path

DATA = Path('10_DFU_Dataset/dfu_outcome_dataset/GSE165816_sc.h5ad')
a = ad.read_h5ad(DATA)
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
obs = a.obs.copy()
gi = {g:i for i,g in enumerate(a.var_names.astype(str))}
tc = obs['tissue_class'].astype(str).values

def gate_mask(genes):
    cols = [gi[g] for g in genes if g in gi]
    sub = np.asarray(X[:, cols].todense(), dtype=np.float64)
    subz = (sub - sub.mean(0))/(sub.std(0)+1e-9)
    fs = subz.mean(1)
    return fs >= np.quantile(fs, .75)

g1 = gate_mask(['DCN','LUM','PDGFRA','COL1A1','COL1A2','FAP','THY1'])
g2 = gate_mask(['DCN','LUM','PDGFRA','VIM'])
sel1 = set(np.where((tc=='wound_bed') & g1)[0].tolist())
sel2 = set(np.where((tc=='wound_bed') & g2)[0].tolist())
inter = len(sel1 & sel2); union = len(sel1 | sel2)
out = {'primary_n': len(sel1), 'cf_n': len(sel2), 'intersection': inter, 'union': union,
       'jaccard': inter/union, 'pct_primary_in_cf': inter/len(sel1), 'pct_cf_in_primary': inter/len(sel2)}
json.dump(out, open('revision_work/reanalysis/S66_gate_overlap.json','w'), indent=2)
print(out)
