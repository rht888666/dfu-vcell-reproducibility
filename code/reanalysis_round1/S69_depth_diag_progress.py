# -*- coding: utf-8 -*-
"""S69: residual-depth diagnostics for Table 5 + GSE231643 progress score restoration."""
import sys, io, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np, pandas as pd
import scipy.sparse as sp
import anndata as ad
from pathlib import Path
from itertools import combinations

DATA = Path('10_DFU_Dataset/dfu_outcome_dataset/GSE165816_sc.h5ad')
OUT = Path('revision_work/reanalysis')

GATES = {
    'primary_S25_gate': ['DCN','LUM','PDGFRA','COL1A1','COL1A2','FAP','THY1'],
    'collagen_free_gate': ['DCN','LUM','PDGFRA','VIM'],
}
PROG = {
 'F1':['COL13A1','COL23A1','SFRP2','CRABP1','APCDD1'],
 'F5':['SCN7A','SOX10','PLP1','NGFR','S100B'],
 'F6':['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
 'F7':['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
 'F7e':['COL1A1','COL1A2','COL3A1','FN1','POSTN','MMP2'],
}

def rank_resid(s, rc):
    rs = pd.Series(s).rank().values
    rs = rs - rs.mean()
    beta = (rc*(rs)).sum()/((rc**2).sum()+1e-12)
    return rs - beta*rc

def lin_resid(s, logu):
    lc = logu - logu.mean()
    beta = (lc*(s-s.mean())).sum()/((lc**2).sum()+1e-12)
    return s - beta*lc

def gate_mask(X, gi, genes):
    cols = [gi[g] for g in genes if g in gi]
    sub = np.asarray(X[:, cols].todense(), dtype=np.float64)
    subz = (sub - sub.mean(0))/(sub.std(0)+1e-9)
    fs = subz.mean(1)
    return fs >= np.quantile(fs, .75)

def pearson(x, y):
    x = x - x.mean(); y = y - y.mean()
    return float((x*y).sum()/np.sqrt((x**2).sum()*(y**2).sum()+1e-30))

print('loading h5ad...')
a = ad.read_h5ad(DATA)
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
obs = a.obs.copy()
gi = {g:i for i,g in enumerate(a.var_names.astype(str))}
tc = obs['tissue_class'].astype(str).values

out = {'gates': {}}
for gname, ggenes in GATES.items():
    isfib = gate_mask(X, gi, ggenes)
    sel = np.where((tc=='wound_bed') & isfib)[0]
    Xsel = X[sel]
    umi = np.asarray(Xsel.sum(1)).ravel()
    logu = np.log1p(umi)
    allp = sorted({g for v in PROG.values() for g in v})
    cols = [gi[g] for g in allp]
    mat = np.asarray(Xsel[:, cols].todense(), dtype=np.float64)
    z = (mat - mat.mean(0))/(mat.std(0)+1e-9)
    gpos = {g:i for i,g in enumerate(allp)}
    lib = obs['sample_id'].astype(str).values[sel]
    gr = {}
    for k in ['F7','F7e']:
        idx = [gpos[g] for g in PROG[k] if g in gpos]
        s = z[:, idx].mean(1)
        schemes = {'raw': s, 'lin': lin_resid(s, logu), 'rank': rank_resid(s, logu - logu.mean())}
        d = {}
        for sn, sv in schemes.items():
            cell_r = pearson(sv, logu)
            # library-level
            lm, ld = {}, {}
            for u in np.unique(lib):
                m = lib == u
                lm[u] = sv[m].mean(); ld[u] = logu[m].mean()
            kk = sorted(lm)
            lib_r = pearson(np.array([lm[u] for u in kk]), np.array([ld[u] for u in kk]))
            d[sn] = {'cell_r': round(cell_r, 4), 'lib_r': round(lib_r, 4)}
        gr[k] = d
        print(gname, k, json.dumps(d))
    out['gates'][gname] = gr

# --- GSE231643 progress score restoration ---
dt = pd.read_csv('revision_work/github_update/dfu-vcell-reproducibility/processed_tables/S20_donor_trajectory.csv')
TP = ['D0','D1','D7','D30']
S = {}
for dnr in sorted(dt.donor.unique()):
    sub = dt[dt.donor==dnr].set_index('tp').loc[TP]
    S[dnr] = sub[['F6','F7']].values
mean_traj = np.mean([S[d] for d in S], axis=0)

def project_progress(points, poly):
    seglen = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    cum = np.concatenate([[0], np.cumsum(seglen)]); total = cum[-1]
    progs = []
    for p in points:
        best = (np.inf, None)
        for kk in range(len(poly)-1):
            A, B = poly[kk], poly[kk+1]; AB = B-A; L2 = AB@AB
            t = 0.0 if L2 < 1e-12 else np.clip((p-A)@AB/L2, 0, 1)
            pr = A + t*AB; dist = np.linalg.norm(p-pr)
            prog = (cum[kk] + t*seglen[kk])/total
            if dist < best[0]: best = (dist, prog)
        progs.append(min(max(best[1],0.0),1.0))
    return np.array(progs)

df = pd.read_csv('revision_work/github_update/dfu-vcell-reproducibility/processed_tables/S12_GSE231643_samples.csv')
pts = df[['F6','F7']].values
prog = project_progress(pts, mean_traj)
df['progress'] = prog
h = df[df.group=='Healer']['progress'].values
n = df[df.group=='NonHealer']['progress'].values
diff = float(h.mean() - n.mean())
# exact two-sided permutation C(8,3)=56
vals = df['progress'].values
groups = df['group'].values
cnt = 0; tot = 0
for idx in combinations(range(8), 3):
    mask = np.zeros(8, bool); mask[list(idx)] = True
    # enumerate nonhealer assignments of size 3
    d2 = vals[~mask].mean() - vals[mask].mean()
    if abs(d2) >= abs(diff) - 1e-12: cnt += 1
    tot += 1
p = cnt/tot
print('GSE231643 progress: diff', round(diff,3), 'p', round(p,4), f'({cnt}/{tot})')
print('  healer mean', round(float(h.mean()),3), 'nonhealer mean', round(float(n.mean()),3))
out['GSE231643_progress'] = {'diff': round(diff,4), 'p_two_sided_exact': round(p,4),
    'n_perm': tot, 'healer_mean': round(float(h.mean()),4), 'nonhealer_mean': round(float(n.mean()),4),
    'caveat': 'cross-dataset projection; z-space computed within GSE231643, not calibrated to the acute-wound reference'}
json.dump(out, open(OUT/'S69_depth_diag_progress.json','w'), indent=1)
print('saved S69 json')
