# -*- coding: utf-8 -*-
"""Reviewer-requested fibroblast-gating sensitivity audit for GSE165816.
Recomputes the primary gate exactly as in S25/S25b and a collagen-free alternative gate.
Outputs library-level rank-residualized F7/F7e contrasts.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.stats import rankdata, spearmanr

ROOT=Path('.')
DATA=ROOT/'10_DFU_Dataset/dfu_outcome_dataset/GSE165816_sc.h5ad'
OUT=ROOT/'revision_work/gating_sensitivity_GSE165816.json'

PROG={
 'F6':['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
 'F7':['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
 'F7e':['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
 'F1':['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
 'F5':['SCN7A','SOX10','PLP1','NGFR','S100B'],
}
PRIMARY=['DCN','LUM','PDGFRA','COL1A1','COL1A2','FAP','THY1']
COLFREE=['DCN','LUM','PDGFRA','VIM']

def perm_test(a,b,nperm=10000,seed=11):
    rng=np.random.default_rng(seed)
    obs=float(a.mean()-b.mean())
    allv=np.concatenate([a,b]); na=len(a); cnt=0
    for _ in range(nperm):
        pm=rng.permutation(allv)
        if (pm[:na].mean()-pm[na:].mean()) >= obs: cnt+=1
    return obs,float((cnt+1)/(nperm+1))

def rank_resid(s,logu):
    rs=rankdata(s,method='average'); ru=rankdata(logu,method='average')
    rc=ru-ru.mean()
    beta=(rc*(rs-rs.mean())).sum()/((rc**2).sum()+1e-12)
    return rs-beta*rc

def upper_quartile_gate(X,obs,gi,genes):
    cols=[gi[g] for g in genes if g in gi]
    sub=np.asarray(X[:,cols].todense(),dtype=np.float64)
    subz=(sub-sub.mean(0))/(sub.std(0)+1e-9)
    fs=subz.mean(1)
    return fs>=np.quantile(fs,.75),fs

def compute_for_gate(X,obs,gi,gate_name,gate_genes):
    isfib,fs=upper_quartile_gate(X,obs,gi,gate_genes)
    tc=obs['tissue_class'].astype(str).values
    sel=np.where((tc=='wound_bed')&isfib)[0]
    Xsel=X[sel]
    umi=np.asarray(Xsel.sum(1)).ravel()
    logu=np.log1p(umi)
    # z-score target genes within selected wound-bed fibroblasts, matching S25 scaling logic
    allprog=sorted({g for v in PROG.values() for g in v})
    cols=[gi[g] for g in allprog]
    mat=np.asarray(Xsel[:,cols].todense(),dtype=np.float64)
    mean=mat.mean(0); sd=mat.std(0)+1e-9
    z=(mat-mean)/sd
    gpos={g:i for i,g in enumerate(allprog)}
    grp=obs['group'].astype(str).values[sel]
    lib=obs['sample_id'].astype(str).values[sel]
    subj=obs['subject_hint'].astype(str).values[sel]
    scores={}
    resscores={}
    for k,gs in PROG.items():
        idx=[gpos[g] for g in gs if g in gpos]
        s=z[:,idx].mean(1)
        scores[k]=s
        resscores[k]=rank_resid(s,logu)
    rows=[]
    for unit in np.unique(lib):
        m=lib==unit
        if m.sum()<30: continue
        row={'library':unit,'subject_hint':pd.Series(subj[m]).mode()[0],
             'group':pd.Series(grp[m]).mode()[0], 'n_cells':int(m.sum()),'umi':float(umi[m].mean())}
        for k in PROG:
            row[k]=float(resscores[k][m].mean())
        rows.append(row)
    P=pd.DataFrame(rows)
    out={'gate':gate_name,'gate_genes':gate_genes,'n_selected_cells':int(len(sel)),'n_libraries':int(len(P)),
         'groups':{g:int((P.group==g).sum()) for g in sorted(P.group.unique())},
         'tests':{},'library_table':P.to_dict(orient='records')}
    for k in ['F7','F7e','F6','F1','F5']:
        h=P[P.group=='Healer'][k].values; n=P[P.group=='NonHealer'][k].values
        if len(h) and len(n):
            d,p=perm_test(h,n)
            out['tests'][k]={'healer_mean':float(h.mean()),'nonhealer_mean':float(n.mean()),'diff':d,'p_one_sided':p,
                             'n_healer_libraries':int(len(h)),'n_nonhealer_libraries':int(len(n))}
    return out

print('loading h5ad')
a=ad.read_h5ad(DATA)
print('loaded',a.shape)
X=a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
obs=a.obs.copy()
gi={g:i for i,g in enumerate(a.var_names.astype(str))}
res={'dataset':'GSE165816','method':'cell-level program scores z-scored within gate-selected wound-bed fibroblasts; average-rank residualization against rank(log1p UMI); library means; one-sided permutation test (10,000 permutations, seed 11)','gates':[]}
for name,genes in [('primary_S25_gate',PRIMARY),('collagen_free_sensitivity_gate',COLFREE)]:
    print('computing',name)
    res['gates'].append(compute_for_gate(X,obs,gi,name,genes))
    print(res['gates'][-1]['tests'])
OUT.write_text(json.dumps(res,indent=2),encoding='utf-8')
print('wrote',OUT)
