# -*- coding: utf-8 -*-
"""S43: exact two-sided permutation tests for GSE231643 and GSE223964 from archived library-level CSVs."""
import json
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

OUT = Path('revision_work/reanalysis')
TAB = Path('revision_work/github_update/dfu-vcell-reproducibility/processed_tables')

def exact_two(a, b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    obs=a.mean()-b.mean()
    allv=np.concatenate([a,b]); n=len(allv); nb=len(b)
    idx=np.arange(n); tot=0; cnt=0
    for comb in combinations(idx, nb):
        m=np.zeros(n,bool); m[list(comb)]=True
        d=allv[~m].mean()-allv[m].mean()
        tot+=1
        if abs(d)>=abs(obs)-1e-12: cnt+=1
    return float(obs), float(cnt/tot), tot

out={}

# ---- GSE231643: 5 Healer vs 3 NonHealer, program scores raw and depth-residualized (_res)
df=pd.read_csv(TAB/'S12_GSE231643_samples.csv')
res={'n_healer':int((df.group=='Healer').sum()),'n_nonhealer':int((df.group=='NonHealer').sum()),
     'unit':'deposited GEO libraries','tests':{},'umi':{}}
for prog in ['F6','F7','F7e','F1','F5','F6i']:
    for sfx in ['', '_res']:
        col=prog+sfx
        if col not in df.columns: continue
        h=df.loc[df.group=='Healer',col].values; nn=df.loc[df.group=='NonHealer',col].values
        o,p2,tot=exact_two(h,nn)
        res['tests'][col]={'diff':o,'p_two_sided_exact':p2,'n_perm':tot,
                           'healer_mean':float(h.mean()),'nonhealer_mean':float(nn.mean())}
o,p2,tot=exact_two(df.loc[df.group=='Healer','umi'].values, df.loc[df.group=='NonHealer','umi'].values)
res['umi']={'diff':o,'p_two_sided_exact':p2,'healer_mean':float(df.loc[df.group=='Healer','umi'].mean()),
            'nonhealer_mean':float(df.loc[df.group=='NonHealer','umi'].mean())}
out['GSE231643']=res
print('GSE231643:')
for col in ['F6','F6_res','F7e','F7e_res','F7','F7_res','F1_res','F5_res']:
    t=res['tests'].get(col)
    if t: print(f"  {col}: diff {t['diff']:+.4f} p2 {t['p_two_sided_exact']:.4f}")
print('  umi:',res['umi'])

# ---- GSE223964: 5 diabetic vs 3 non-diabetic ulcer (disease status, NOT healing outcome)
df2=pd.read_csv(TAB/'S21_GSE223964_samples.csv')
res2={'n_diabetic':int((df2.condition=='diabetic').sum()),'n_nondiabetic':int((df2.condition=='non diabetic').sum()),
      'unit':'libraries','contrast':'diabetic vs non-diabetic ulcer (disease status, not healing outcome)','tests':{}}
for prog in ['F6','F7','F7e','F1','F5','F6i']:
    if prog not in df2.columns: continue
    h=df2.loc[df2.condition=='diabetic',prog].values; nn=df2.loc[df2.condition=='non diabetic',prog].values
    o,p2,tot=exact_two(h,nn)
    res2['tests'][prog]={'diff':o,'p_two_sided_exact':p2,'n_perm':tot,
                         'diabetic_mean':float(h.mean()),'nondiabetic_mean':float(nn.mean())}
out['GSE223964']=res2
print('GSE223964 (diabetic vs non-diabetic):')
for col,t in res2['tests'].items():
    print(f"  {col}: diff {t['diff']:+.4f} p2 {t['p_two_sided_exact']:.4f}")

with open(OUT/'S43_cohort_exact_tests.json','w') as f:
    json.dump(out,f,indent=2)
print('saved S43')
