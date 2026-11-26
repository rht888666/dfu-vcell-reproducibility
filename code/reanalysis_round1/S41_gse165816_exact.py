# -*- coding: utf-8 -*-
"""
S41: GSE165816 exact two-sided reanalysis.
- primary gate + collagen-free gate (upper quartile, as S25/S37)
- wound_bed cells; z-scored raw counts within gate-selected cells
- raw / linear-OLS-residualized / rank-residualized (pooled across selected cells) scores
- library-level (sample_id) aggregation
- EXACT two-sided permutation tests by exhaustive enumeration of C(14,5)=2002 label assignments
- composite progress projection onto acute-wound reference polyline (heuristic, cross-dataset z units)
- depth diagnostics at library level
"""
import json
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.stats import rankdata, spearmanr

OUT = Path('revision_work/reanalysis'); OUT.mkdir(parents=True, exist_ok=True)
DATA = Path('10_DFU_Dataset/dfu_outcome_dataset/GSE165816_sc.h5ad')

PROG = {
 'F6':['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
 'F7':['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
 'F7e':['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
 'F1':['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
 'F5':['SCN7A','SOX10','PLP1','NGFR','S100B'],
}
GATES = {
 'primary_S25_gate': ['DCN','LUM','PDGFRA','COL1A1','COL1A2','FAP','THY1'],
 'collagen_free_gate': ['DCN','LUM','PDGFRA','VIM'],
}

def exact_two_sided(a, b):
    """Exact two-sided permutation test of mean difference by exhaustive enumeration.
    a=Healer values, b=NonHealer values. Enumerates all ways to choose len(b) of n units as group b."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    obs = a.mean() - b.mean()
    allv = np.concatenate([a, b]); n = len(allv); nb = len(b)
    total = 0; cnt = 0
    idx = np.arange(n)
    for comb in combinations(idx, nb):
        mask = np.zeros(n, bool); mask[list(comb)] = True
        d = allv[~mask].mean() - allv[mask].mean()
        total += 1
        if abs(d) >= abs(obs) - 1e-12: cnt += 1
    return float(obs), float(cnt/total), int(total)

def rank_resid(s, logu):
    rs = rankdata(s, method='average'); ru = rankdata(logu, method='average')
    rc = ru - ru.mean()
    beta = (rc*(rs-rs.mean())).sum()/((rc**2).sum()+1e-12)
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

# reference polyline from acute wound donor means (S20 donor trajectory, residualized z units)
dt = pd.read_csv('revision_work/github_update/dfu-vcell-reproducibility/processed_tables/S20_donor_trajectory.csv')
TP = ['D0','D1','D7','D30']
S = {}
for dnr in sorted(dt.donor.unique()):
    sub = dt[dt.donor==dnr].set_index('tp').loc[TP]
    S[dnr] = sub[['F6','F7']].values
mean_traj = np.mean([S[d] for d in S], axis=0)  # 4x2 in (F6,F7)

def project_progress(points, poly):
    """points: (n,2); polyline reference; returns clipped arc-length progress in [0,1] and distances."""
    seglen = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    cum = np.concatenate([[0], np.cumsum(seglen)]); total = cum[-1]
    progs, dists = [], []
    for p in points:
        best = (np.inf, None)
        for k in range(len(poly)-1):
            A, B = poly[k], poly[k+1]; AB = B-A; L2 = AB@AB
            t = 0.0 if L2 < 1e-12 else np.clip((p-A)@AB/L2, 0, 1)
            pr = A + t*AB; dist = np.linalg.norm(p-pr)
            prog = (cum[k] + t*seglen[k])/total
            if dist < best[0]: best = (dist, prog)
        dists.append(best[0]); progs.append(min(max(best[1],0.0),1.0))
    return np.array(progs), np.array(dists)

print('loading h5ad...')
a = ad.read_h5ad(DATA)
print('loaded', a.shape)
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
obs = a.obs.copy()
gi = {g:i for i,g in enumerate(a.var_names.astype(str))}
tc = obs['tissue_class'].astype(str).values

results = {'dataset':'GSE165816','method':'z-scored raw UMI counts within gate-selected wound-bed fibroblasts; raw / linear-OLS / rank residualization against log1p(UMI) computed POOLED across all gate-selected cells; library (sample_id) means; exact two-sided permutation by exhaustive enumeration','gates':{}}

for gname, ggenes in GATES.items():
    print('gate', gname)
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
    grp = obs['group'].astype(str).values[sel]
    lib = obs['sample_id'].astype(str).values[sel]
    subj = obs['subject_hint'].astype(str).values[sel]
    scores = {}
    for k, gs in PROG.items():
        idx = [gpos[g] for g in gs if g in gpos]
        s = z[:, idx].mean(1)
        scores[k] = {'raw': s, 'lin': lin_resid(s, logu), 'rank': rank_resid(s, logu)}
    # library table
    rows = []
    for unit in np.unique(lib):
        m = lib == unit
        if m.sum() < 30: continue
        row = {'library': unit, 'subject_hint': pd.Series(subj[m]).mode()[0],
               'group': pd.Series(grp[m]).mode()[0], 'n_cells': int(m.sum()), 'umi': float(umi[m].mean())}
        for k in PROG:
            for sc in ['raw','lin','rank']:
                row[f'{k}_{sc}'] = float(scores[k][sc][m].mean())
        rows.append(row)
    P = pd.DataFrame(rows)
    P.to_csv(OUT/f'S41_library_table_{gname}.csv', index=False)
    gr = {'n_selected_cells': int(len(sel)), 'n_libraries': int(len(P)),
          'group_counts': {g:int((P.group==g).sum()) for g in sorted(P.group.unique())},
          'tests': {}, 'progress': {}, 'depth': {}}
    hmask = P.group=='Healer'; nmask = P.group=='NonHealer'
    print('  libraries H/NH:', hmask.sum(), nmask.sum())
    for k in PROG:
        for sc in ['raw','lin','rank']:
            col = f'{k}_{sc}'
            o,p2,tot = exact_two_sided(P.loc[hmask,col].values, P.loc[nmask,col].values)
            gr['tests'][col] = {'diff': o, 'p_two_sided_exact': p2, 'n_permutations': tot,
                                'healer_mean': float(P.loc[hmask,col].mean()), 'nonhealer_mean': float(P.loc[nmask,col].mean())}
    # UMI difference test
    o,p2,tot = exact_two_sided(P.loc[hmask,'umi'].values, P.loc[nmask,'umi'].values)
    gr['tests']['umi'] = {'diff': o, 'p_two_sided_exact': p2, 'healer_mean': float(P.loc[hmask,'umi'].mean()), 'nonhealer_mean': float(P.loc[nmask,'umi'].mean())}
    # composite progress on rank-residualized (F6,F7) library means, clipped [0,1]
    pts = P[['F6_rank','F7_rank']].values
    prog, dist = project_progress(pts, mean_traj)
    P['progress'] = prog; P['progress_dist'] = dist
    o,p2,tot = exact_two_sided(P.loc[hmask,'progress'].values, P.loc[nmask,'progress'].values)
    gr['progress'] = {'diff': o, 'p_two_sided_exact': p2,
                      'healer_mean': float(P.loc[hmask,'progress'].mean()), 'nonhealer_mean': float(P.loc[nmask,'progress'].mean()),
                      'n_clipped_0': int((prog==0).sum()), 'n_clipped_1': int((prog==1).sum()),
                      'note':'heuristic cross-dataset projection in z units; clipped to [0,1]'}
    # depth correlations at library level
    for col in ['F7_raw','F7e_raw','F6_raw','progress']:
        v = P[col].values
        r_all = spearmanr(v, P['umi'].values)
        dful = P[P.group.isin(['Healer','NonHealer'])]
        r_dfu = spearmanr(dful[col].values, dful['umi'].values)
        gr['depth'][col] = {'rho_all_foot_libraries': float(r_all[0]), 'p_all': float(r_all[1]),
                            'rho_dfu_only': float(r_dfu[0]), 'p_dfu': float(r_dfu[1])}
    r = spearmanr(P['F7_raw'].values, np.log(P['umi'].values))
    gr['depth']['F7_raw_vs_logumi'] = {'rho_all': float(r[0]), 'p_all': float(r[1])}
    results['gates'][gname] = gr
    P.to_csv(OUT/f'S41_library_table_{gname}.csv', index=False)

with open(OUT/'S41_gse165816_exact.json','w') as f:
    json.dump(results, f, indent=2)
print('saved S41')
g = results['gates']['primary_S25_gate']
for col in ['F7_raw','F7_rank','F7e_raw','F7e_rank','F6_rank','F1_rank','F5_rank','umi']:
    t = g['tests'][col]
    print(col, 'diff', round(t['diff'],4), 'p2', round(t['p_two_sided_exact'],4))
print('progress', g['progress'])
print('depth', json.dumps(g['depth'], indent=1)[:800])
