# -*- coding: utf-8 -*-
import pandas as pd, numpy as np, json
OUT='17_DFU_VirtualCell_10plus'

def spearmanr(x,y):
    xr=pd.Series(x).rank(method='average').values
    yr=pd.Series(y).rank(method='average').values
    r=float(np.corrcoef(xr,yr)[0,1])
    return r, np.nan
SP=pd.read_csv(OUT+'/S26_spatial_scores.csv')
ED=pd.read_csv(OUT+'/S29_edge_dist.csv')
SP['bc']=SP['barcode'].astype(str); ED['bc']=ED['barcode'].astype(str)
M=SP.merge(ED[['bc','h5_stem','edge_dist']],on=['bc','h5_stem'],how='inner')
zcols=[c for c in M.columns if c.endswith('_z')]
log_n=np.log1p(M['total_counts'].values.astype(float)); lc=log_n-log_n.mean(); den=(lc**2).sum()+1e-12
GZ={}
for c in zcols:
 v=M[c].values.astype(float); beta=(lc*(v-v.mean())).sum()/den; GZ[c[:-2]]=v-beta*lc
GZ=pd.DataFrame(GZ,index=M.index)
MARK={'F6':['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],'F7':['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],'F7e':['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],'F1':['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],'F5':['SCN7A','SOX10','PLP1','NGFR','S100B']}
res={}
for cond in ['Wound7']:
 res[cond]={}
 for prog,genes in MARK.items():
  gs=[g for g in genes if g in GZ.columns]
  res[cond][prog]={}
  for donor in sorted(M['Donor'].dropna().unique()):
   mm=(M['Condition'].values==cond)&(M['Donor'].values==donor)
   if mm.sum()<30: continue
   r,p=spearmanr(M.loc[mm,'edge_dist'].values,GZ.loc[mm,gs].mean(1).values)
   secs=sorted(M.loc[mm,'h5_stem'].unique())
   res[cond][prog][donor]={'rho':float(r),'p_cellspot_descriptive':float(p),'n_spots':int(mm.sum()),'n_sections':len(secs),'sections':secs}
# section-level sign consistency for F6 Wound7
sec=[]
for (donor,h5),g in M[M['Condition'].eq('Wound7')].groupby(['Donor','h5_stem']):
 if len(g)<30: continue
 idx=g.index
 r,p=spearmanr(g['edge_dist'].values,GZ.loc[idx,[x for x in MARK['F6'] if x in GZ.columns]].mean(1).values)
 sec.append({'donor':donor,'section':h5,'rho':float(r),'p_descriptive':float(p),'n_spots':int(len(g))})
out={'wound7_donor_level':res,'wound7_f6_section_level':sec,'n_sections_total':len(sec),'n_sections_f6_negative':sum(x['rho']<0 for x in sec)}
print(json.dumps(out,indent=2))
open('revision_work/spatial_donor_sensitivity.json','w').write(json.dumps(out,indent=2))
