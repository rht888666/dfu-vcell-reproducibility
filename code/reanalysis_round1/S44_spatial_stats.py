# -*- coding: utf-8 -*-
"""S44: spatial per-donor compartment statistics + ablation panels + 27-gene candidate list export."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path('revision_work/reanalysis')
SP = pd.read_csv('revision_work/github_update/dfu-vcell-reproducibility/processed_tables/S26_spatial_scores.csv')
print('spots', len(SP), 'cols', len(SP.columns))

PROG = {
 'F6':['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
 'F7':['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
 'F7e':['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
 'F7noCT':['LRRC15','SFRP4','ASPN','POSTN','CTHRC1'],
 'F7epCT':['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX','ACTA2','TAGLN'],
 'F1':['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
 'F5':['SCN7A','SOX10','PLP1','NGFR','S100B'],
}
# gene-level depth-residualized z columns (GZ approach from S26b): recompute here
zcols=[c for c in SP.columns if c.endswith('_z')]
GENES=[c[:-2] for c in zcols]
logn=np.log1p(SP['total_counts'].values.astype(float))
lc=logn-logn.mean(); den=(lc**2).sum()+1e-12
GZ={}
for g in GENES:
    v=SP[g+'_z'].values.astype(float)
    beta=(lc*(v-v.mean())).sum()/den
    GZ[g]=v-beta*lc
GZ=pd.DataFrame(GZ)

avail={}
sc={}
for p,gs in PROG.items():
    use=[g for g in gs if g in GZ.columns]
    avail[p]=use
    if use:
        sc[p]=GZ[use].mean(1).values
print('gene availability:',{k:len(v) for k,v in avail.items()})

# candidate panel gene list actually used
allg=sorted({g for v in PROG.values() for g in v if g in GZ.columns})
print('candidate panel union size (excl F7noCT dup logic):', len(allg), allg)

out={'n_spots':int(len(SP)),'gene_availability':avail,'panel_genes':allg,
     'compartment_by_donor':{},'smooth_muscle_scores':{},'fb1_contrasts':{}}

donors=sorted(SP['Donor'].dropna().unique())
ann=SP['AnnoType'].astype(str).values
for comp in ['Smooth_muscle','FB1']:
    per={}
    for d in donors:
        m=(SP['Donor']==d).values&(ann==comp)
        row={'n_spots':int(m.sum())}
        for p in sc:
            row[p]=float(np.nanmean(sc[p][m])) if m.sum()>0 else None
        per[d]=row
    out['compartment_by_donor'][comp]=per
    # pooled across donors too
    m=(ann==comp)
    out['smooth_muscle_scores' if comp=='Smooth_muscle' else 'fb1_contrasts'][comp+'_pooled']={
        p: float(np.nanmean(sc[p][m])) for p in sc}
    print(f'\n{comp}:')
    for d in donors:
        print(' ',d, {p:round(per[d][p],3) if per[d][p] is not None else None for p in ['F7','F7noCT','F7e','F7epCT']}, 'n=',per[d]['n_spots'])

# FB1 vs unweighted mean of other fibroblast categories (descriptive, matching S35 observed definition)
FIBRO=['FB1','FB2','Sweatgland_FB','Dermis1','Dermis2','Papillary_dermis']
cats=pd.Categorical(ann)
catmean={}
for c in FIBRO:
    m=(ann==c)
    if m.sum()>0:
        catmean[c]={p: float(np.nanmean(sc[p][m])) for p in sc}
out['category_means']=catmean
fb1desc={}
others=[c for c in FIBRO if c!='FB1' and c in catmean]
for p in sc:
    fb1desc[p]={'fb1':catmean['FB1'][p],
                'other_fib_unweighted_mean':float(np.mean([catmean[c][p] for c in others])),
                'delta': catmean['FB1'][p]-float(np.mean([catmean[c][p] for c in others]))}
out['fb1_descriptive_contrast']=fb1desc
print('\nFB1 descriptive deltas:', {p:round(fb1desc[p]['delta'],3) for p in ['F7','F7e','F7noCT','F1','F5']})

# temporal per-donor spearman for key panels incl ablation (cond ordinal)
from scipy.stats import spearmanr
ord_map={'Skin':0,'Wound1':1,'Wound7':2,'Wound30':3}
SP['cond_ord']=SP['Condition'].map(ord_map)
temp={}
for p in ['F7','F7noCT','F7e','F7epCT','F6','F1','F5']:
    rhos=[]
    for d in donors:
        m=(SP['Donor']==d).values
        if m.sum()<50: continue
        r,_=spearmanr(SP.loc[m,'cond_ord'], sc[p][m])
        rhos.append(float(r))
    temp[p]={'donor_rhos':rhos,'mean':float(np.mean(rhos))}
out['temporal']=temp
print('\ntemporal donor-mean rho:',{p:round(temp[p]['mean'],3) for p in temp})

with open(OUT/'S44_spatial_stats.json','w') as f:
    json.dump(out,f,indent=2)
print('saved S44')
