# -*- coding: utf-8 -*-
"""
S42: CRISPRa axis audit and recomputation.
Recomputes Axis A / Axis-M7 / F7e-5gene pure structural vector rankings
using the exact LOO-cosine + axis-shuffle permutation method of the archived
04_insilico_screen script. Verifies whether manuscript Table 6 values reproduce.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path('revision_work/github_update/dfu-vcell-reproducibility')
OUT = Path('revision_work/reanalysis')
OUT.mkdir(parents=True, exist_ok=True)

d = np.load(ROOT/'data/perturb_V3.npz', allow_pickle=True)
M = d['M'].astype(np.float64)
tfs = np.array([str(x) for x in d['tfs']])
genes = np.array([str(x) for x in d['genes']])
gene_idx = {g: i for i, g in enumerate(genes)}
print('matrix', M.shape)

def build(pos, neg):
    v = np.zeros(len(genes))
    for g in pos:
        if g in gene_idx: v[gene_idx[g]] = 1.0
    for g in neg:
        if g in gene_idx: v[gene_idx[g]] = -1.0
    return v

AXES = {
    'AxisA_contractile_inclusive': build(['CTHRC1','ACTA2','TAGLN','COL1A1','COL1A2','COL3A1','FN1','LOX'],['MMP3','CCL2']),
    'AxisM7_measured_matrix': build(['COL1A1','COL1A2','COL3A1','FN1','LOX'],['MMP3','CCL2']),
    'F7e5_pure_structural': build(['COL1A1','COL1A2','COL3A1','FN1','LOX'],[]),
}

def screen(axis_vec, n_perm=2000, seed=20261010):
    axis_self = set(np.where(axis_vec != 0)[0].tolist())
    loo = np.zeros(M.shape[0])
    for i in range(M.shape[0]):
        a_i = axis_vec.copy()
        if tfs[i] in gene_idx:
            j = gene_idx[tfs[i]]
            if j in axis_self: a_i[j] = 0.0
        loo[i] = float(np.dot(M[i], a_i) / ((np.linalg.norm(M[i])+1e-12)*(np.linalg.norm(a_i)+1e-12)))
    rng = np.random.default_rng(seed)
    null_cos = np.empty((n_perm, M.shape[0]))
    norm_X = np.linalg.norm(M, axis=1) + 1e-12
    for k in range(n_perm):
        ap = axis_vec.copy()
        rng.shuffle(ap)
        null_cos[k] = (M @ ap) / (norm_X * (np.linalg.norm(ap)+1e-12))
    mu = null_cos.mean(axis=0); sd = null_cos.std(axis=0) + 1e-12
    z = (loo - mu) / sd
    return loo, z, mu, sd

out = {'n_tf_rows': int(M.shape[0]), 'n_genes': int(M.shape[1]), 'axes': {}}
allz = {}
for name, v in AXES.items():
    loo, z, mu, sd = screen(v)
    order = np.argsort(-loo)
    ranks = np.empty(len(loo), dtype=int); ranks[order] = np.arange(1, len(loo)+1)
    allz[name] = (loo, z, ranks)
    key_tfs = ['EGR3','EOMES','EGR2','EBF2','ZNF827','SATB1','CEBPB','KLF4','KLF5','PLAGL1','GATA6']
    kt = {}
    for t in key_tfs:
        if t in tfs:
            i = list(tfs).index(t)
            kt[t] = {'cos': float(loo[i]), 'z': float(z[i]), 'rank': int(ranks[i])}
    out['axes'][name] = {
        'nonzero_genes': sorted([genes[i] for i in np.where(v!=0)[0]]),
        'top10': [{'tf': tfs[i], 'cos': float(loo[i]), 'z': float(z[i])} for i in order[:10]],
        'key_tfs': kt,
    }
    print(f'\n== {name} ==')
    print('top5:', [(tfs[i], round(float(loo[i]),4), round(float(z[i]),2)) for i in order[:5]])
    print('EGR3:', kt.get('EGR3'))

# LOGO for AxisM7: drop each nonzero entry in turn
v = AXES['AxisM7_measured_matrix']
nz = np.where(v != 0)[0]
logo = {}
for j in nz:
    g = genes[j]
    v2 = v.copy(); v2[j] = 0.0
    loo, z, mu, sd = screen(v2)
    order = np.argsort(-loo)
    ranks = np.empty(len(loo), dtype=int); ranks[order] = np.arange(1, len(loo)+1)
    i = list(tfs).index('EGR3')
    logo[g] = {'EGR3_cos': float(loo[i]), 'EGR3_rank': int(ranks[i]), 'EGR3_z': float(z[i])}
    print(f'LOGO drop {g}: EGR3 rank {ranks[i]}, cos {loo[i]:.4f}')
out['axisM7_LOGO'] = logo

# EGR3 perturbation profile values for key genes
i = list(tfs).index('EGR3')
prof = {}
for g in ['COL1A1','COL1A2','COL3A1','FN1','LOX','ACTA2','TAGLN','CTHRC1','MMP3','CCL2','POSTN']:
    if g in gene_idx:
        prof[g] = float(M[i, gene_idx[g]])
    else:
        prof[g] = None
out['EGR3_profile'] = prof
print('\nEGR3 profile:', prof)
# also EBF2 and ZNF827 for contractile comparison
for t in ['EBF2','ZNF827','EOMES','KLF4']:
    ii = list(tfs).index(t)
    out[f'{t}_profile'] = {g: (float(M[ii, gene_idx[g]]) if g in gene_idx else None) for g in prof}

with open(OUT/'S42_crispra_axis_audit.json','w') as f:
    json.dump(out, f, indent=2)
# full ranking csvs
rows = {'TF': tfs}
for name,(loo,z,ranks) in allz.items():
    rows[f'cos_{name}'] = loo; rows[f'z_{name}'] = z; rows[f'rank_{name}'] = ranks
pd.DataFrame(rows).to_csv(OUT/'S42_full_rankings.csv', index=False)
print('\nsaved S42 outputs')
