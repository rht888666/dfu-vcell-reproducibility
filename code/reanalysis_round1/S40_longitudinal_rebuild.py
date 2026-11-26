# -*- coding: utf-8 -*-
"""
S40: longitudinal rebuild - single source of truth for Table 3 and Figure 1.
Same loading/scoring as S20 (z-scored raw counts, linear depth residualization),
plus ablation panels F7noCT (F7 minus ACTA2/TAGLN) and F7epCT (F7e plus ACTA2/TAGLN).
Outputs donor x timepoint x program table and per-panel angular dispersion.
"""
import json, tarfile, zipfile, gzip
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
import scipy.sparse as sp

OUT = Path('revision_work/reanalysis'); OUT.mkdir(parents=True, exist_ok=True)
G = Path('10_DFU_Dataset/gse241132')

M = {
 'F1': ['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
 'F5': ['SCN7A','SOX10','PLP1','NGFR','S100B'],
 'F6': ['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
 'F6i': ['IL6','CXCL1','CCL2','CXCL8','VEGFA'],
 'F7': ['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
 'F7noCT': ['LRRC15','SFRP4','ASPN','POSTN','CTHRC1'],
 'F7e': ['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
 'F7epCT': ['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX','ACTA2','TAGLN'],
}
ORDER = ['Skin','Wound1','Wound7','Wound30']
LABEL = {'Skin':'D0','Wound1':'D1','Wound7':'D7','Wound30':'D30'}

print('reading metadata')
meta = pd.read_csv(G/'cell_metadata.txt.gz', sep='\t')
print('meta', meta.shape)

# locate sample matrices (reuse S20 logic)
ex = G/'extracted'
files = [f for f in ex.rglob('*') if f.is_file()]
zips = sorted([f for f in files if f.suffix=='.zip'])
if zips:
    unz = ex/'unzipped'
    if not unz.exists():
        unz.mkdir(parents=True, exist_ok=True)
        for z in zips:
            with zipfile.ZipFile(z) as zf:
                for m in zf.namelist():
                    if '__MACOSX' in m or m.endswith('/'): continue
                    zf.extract(m, unz)
    files = [f for f in unz.rglob('*') if f.is_file()]
mtxs = sorted([f for f in files if 'matrix.mtx' in f.name])
samples = {}
for m in mtxs:
    d = m.parent
    bcs = list(d.glob('*barcodes.tsv*')); gns = list(d.glob('*features.tsv*'))+list(d.glob('*genes.tsv*'))
    if bcs and gns:
        samples[d.name] = {'mtx': m, 'genes': gns[0], 'barcodes': bcs[0]}
print('samples', len(samples))

ALLGENES = sorted({g for v in M.values() for g in v})

def read_10x_mtx(mtx, genes_f, barcodes_f, wanted):
    g = pd.read_csv(genes_f, sep='\t', header=None)
    gname = (g[1] if g.shape[1]>=2 else g[0]).astype(str).values
    wanted=set(wanted); keep_rows={}
    for i,nm in enumerate(gname):
        if nm in wanted and nm not in keep_rows: keep_rows[nm]=i+1
    kept=list(keep_rows.keys())
    row2col={r:k for k,r in enumerate(keep_rows.values())}
    bcs=pd.read_csv(barcodes_f,sep='\t',header=None)[0].astype(str).values
    rows,cols,vals=[],[],[]
    with gzip.open(mtx,'rt') as fh:
        for line in fh:
            if line.startswith('%'): continue
            parts=line.split()
            if len(parts)<3: continue
            r=int(parts[0])
            if r in row2col:
                rows.append(row2col[r]); cols.append(int(parts[1])-1); vals.append(float(parts[2]))
    Mm=sp.csr_matrix((vals,(cols,rows)),shape=(len(bcs),len(kept)))
    df=pd.DataFrame(Mm.toarray(),columns=kept)
    return df,bcs

frames=[]
for sid,spec in samples.items():
    df,bcs=read_10x_mtx(spec['mtx'],spec['genes'],spec['barcodes'],ALLGENES)
    df['barcode']=bcs; df['sample']=sid
    frames.append(df)
    print(' loaded',sid,df.shape)
expr=pd.concat(frames,ignore_index=True)
expr['bc_key']=expr.apply(lambda r: str(r['barcode']) if str(r['barcode']).startswith(str(r['sample'])+'_') else f"{r['sample']}_{r['barcode']}",axis=1)
meta['bc_key']=meta['barcode'].astype(str)
merged=expr.merge(meta[['bc_key','Patient','Condition','newMainCellTypes','nCount_RNA']],on='bc_key',how='inner')
merged=merged[merged['newMainCellTypes']=='Fibroblast'].copy()
print('fibroblasts',len(merged))

merged['logumi']=np.log1p(merged['nCount_RNA'].astype(float))
logu=merged['logumi'].values; logu_c=logu-logu.mean()

def score_program(df,genes):
    gs=[g for g in genes if g in df.columns]
    sub=df[gs].astype(float)
    z=(sub-sub.mean(0))/(sub.std(0)+1e-9)
    return z.mean(1).values

scored={}
for prog,genes in M.items():
    s=score_program(merged,genes)
    lc=logu_c
    beta=(lc*(s-s.mean())).sum()/((lc**2).sum()+1e-12)
    scored[prog]=s-beta*lc  # linear depth residualization, same as S20
    merged[f'rr_{prog}']=scored[prog]

# donor x timepoint table
rows=[]
for (d,t),sub in merged.groupby(['Patient','Condition']):
    if len(sub)<20: continue
    r={'donor':d,'timepoint':t,'tp':LABEL[t],'n':int(len(sub)),'logumi':float(sub['logumi'].mean())}
    for prog in M:
        r[prog]=float(np.nanmean(sub[f'rr_{prog}']))
    rows.append(r)
ddf=pd.DataFrame(rows)
ddf.to_csv(OUT/'S40_donor_timepoint_programs.csv',index=False)
print(ddf.round(3).to_string())

# angular dispersion per matrix panel: D0->D30 displacement angle in (F6, panel)
def circ_stats(angles_deg):
    th=np.deg2rad(angles_deg)
    R=np.abs(np.mean(np.exp(1j*th)))
    spread=np.degrees(np.sqrt(2*(1-R)))
    return float(R),float(spread)
ang={}
TP=['D0','D1','D7','D30']
for panel in ['F7','F7noCT','F7e','F7epCT']:
    per={}
    for d in sorted(ddf.donor.unique()):
        sub=ddf[ddf.donor==d].set_index('tp').loc[TP]
        dF6=sub['F6'].iloc[-1]-sub['F6'].iloc[0]
        dM=sub[panel].iloc[-1]-sub[panel].iloc[0]
        per[d]={'dF6':float(dF6),'dM':float(dM),'angle':float(np.degrees(np.arctan2(dM,dF6)))}
    angles=[per[d]['angle'] for d in per]
    R,spread=circ_stats(angles)
    ang[panel]={'per_donor':per,'angles':angles,'range':float(max(angles)-min(angles)),
                'R':R,'circular_spread':spread,'mean_angle':float(np.mean(angles))}
    print(panel,'angles',[round(a,2) for a in angles],'range',round(max(angles)-min(angles),2),'R',round(R,4))

# random 7-gene subsample reference for angular range (contextualize F7noCT vs random)
rng=np.random.default_rng(20261006)
allpanel=['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN']
# note: interpretation aid only; computed across 200 random 5-gene subsets of the 7 F7 genes
res=[]
for _ in range(200):
    gs=list(rng.choice(allpanel,size=5,replace=False))
    s=score_program(merged,gs)
    lc=logu_c
    beta=(lc*(s-s.mean())).sum()/((lc**2).sum()+1e-12)
    rr=s-beta*lc
    merged['tmp']=rr
    angles=[]
    for d in sorted(ddf.donor.unique()):
        sub=merged[merged['Patient']==d]
        v0=np.nanmean(sub.loc[sub['Condition']=='Skin','tmp'])
        v3=np.nanmean(sub.loc[sub['Condition']=='Wound30','tmp'])
        f60=np.nanmean(sub.loc[sub['Condition']=='Skin','rr_F6'])
        f63=np.nanmean(sub.loc[sub['Condition']=='Wound30','rr_F6'])
        angles.append(float(np.degrees(np.arctan2(v3-v0,f63-f60))))
    res.append(max(angles)-min(angles))
res=np.array(res)
ang['random_5of7_range']={'mean':float(res.mean()),'p95':float(np.percentile(res,95)),
                          'frac_le_F7e':float((res<=ang['F7e']['range']).mean()),
                          'frac_le_F7noCT':float((res<=ang['F7noCT']['range']).mean())}
print('random 5-of-7 range mean',round(res.mean(),2),'frac<=F7e range',round((res<=ang['F7e']['range']).mean(),3))

# cell-level spearman per donor for F6/F7e (for table)
from scipy.stats import spearmanr
tmap={t:i for i,t in enumerate(ORDER)}
merged['t_ord']=merged['Condition'].map(tmap)
corr={}
for prog in ['F6','F6i','F7','F7noCT','F7e','F7epCT','F1','F5']:
    per=[]
    for d in sorted(merged.Patient.unique()):
        sub=merged[merged.Patient==d]
        r,_=spearmanr(sub['t_ord'],sub[f'rr_{prog}'])
        per.append(float(r))
    corr[prog]=per
out={'angular':ang,'correlations':corr,'n_cells':int(len(merged)),
     'donor_cellcounts':{d:{LABEL[t]:int(len(merged[(merged.Patient==d)&(merged.Condition==t)])) for t in ORDER} for d in sorted(merged.Patient.unique())}}
with open(OUT/'S40_longitudinal_rebuild.json','w') as f:
    json.dump(out,f,indent=2)
print('saved S40')
print('cellcounts',out['donor_cellcounts'])
