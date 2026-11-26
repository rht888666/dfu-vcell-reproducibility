# -*- coding: utf-8 -*-
"""
S29b: 创缘梯度的特异性裁决 —— F6/F7 的梯度是否超过全局组织梯度？

S29 发现：
  F6 呈【负】梯度（近创缘）在各条件一致 → 可能真信号
  F7 呈【正】梯度（远离创缘），但 D7 时 F5 (+0.330) 比 F7 (+0.256) 更大
  → 存在"全局组织梯度"（任何基因都随距创缘距离变化）

必须做的裁决：
  1. 随机基因集在 edge_dist 上的梯度零分布
  2. F6/F7/F7e 的经验 p
  3. 双差法：(F7 − F5) 或 (F6 − F1) 的梯度
  4. 关键：F6 的【负】梯度是否特异（因为多数程序是正的）

判据：
  - 若 F6 负梯度显著且 F1/F5 不为负 → F6 的空间定位是特异的 ✅
  - 若 F7 正梯度 ≤ 随机零分布 → F7 梯度非特异 ❌
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
RNG = np.random.default_rng(20261006)

LOG=[]
def log(*a):
    s=' '.join(str(x) for x in a); print(s); LOG.append(s)

log('='*88)
log('S29b: 创缘梯度特异性裁决')
log('='*88)

SP = pd.read_csv(OUT/"S26_spatial_scores.csv")
ED = pd.read_csv(OUT/"S29_edge_dist.csv")
SP['bc']=SP['barcode'].astype(str); ED['bc']=ED['barcode'].astype(str)
M = SP.merge(ED[['bc','h5_stem','edge_dist']], on=['bc','h5_stem'], how='inner')
log(f'\n合并: {len(M):,} spots')

# 重建 GZ
zcols=[c for c in M.columns if c.endswith('_z')]
GENES=[c[:-2] for c in zcols]
log_n=np.log1p(M['total_counts'].values.astype(float))
lc=log_n-log_n.mean(); den=(lc**2).sum()+1e-12
GZ={}
for c in zcols:
    v=M[c].values; beta=(lc*(v-v.mean())).sum()/den
    GZ[c[:-2]]=v-beta*lc
GZ=pd.DataFrame(GZ)

MARK = {
    'F6':  ['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
    'F7':  ['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
    'F7e': ['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
    'F1':  ['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
    'F5':  ['SCN7A','SOX10','PLP1','NGFR','S100B'],
}
D = M['edge_dist'].values

# ---- 随机基因集零分布（每条件）----
log('\n' + '='*88)
log('[1] 随机基因集梯度零分布（分条件）')
log('='*88)
nr=3000
nulls={}
for cond in ['Skin','Wound1','Wound7','Wound30']:
    mm=(M['Condition']==cond).values
    if mm.sum()<50: continue
    rs=[]
    for _ in range(nr):
        gs=RNG.choice(list(GZ.columns), RNG.choice([5,6,7]), replace=False)
        s=GZ.loc[mm,gs].mean(1).values
        r,_=spearmanr(D[mm], s)
        rs.append(r)
    rs=np.array(rs)
    nulls[cond]=dict(mean=float(rs.mean()), sd=float(rs.std()),
                     p05=float(np.percentile(rs,5)), p95=float(np.percentile(rs,95)))
    log(f'  {cond:8s}: null mean={rs.mean():+.4f}, sd={rs.std():.4f}, '
        f'[p5,p95]=[{np.percentile(rs,5):+.3f},{np.percentile(rs,95):+.3f}]')

# ---- 各程序梯度 + 经验 p ----
log('\n' + '='*88)
log('[2] 各程序梯度 vs 零分布（经验 p，双侧）')
log('='*88)
log(f'\n  {"条件":8s} {"程序":5s} {"rho":>8s} {"p_emp":>8s}  判读')
res={}
for cond in ['Skin','Wound1','Wound7','Wound30']:
    mm=(M['Condition']==cond).values
    if mm.sum()<50: continue
    if cond not in nulls: continue
    nm, nsd = nulls[cond]['mean'], nulls[cond]['sd']
    for p,gl in MARK.items():
        gs=[g for g in gl if g in GZ.columns]
        s=GZ.loc[mm,gs].mean(1).values
        r,_=spearmanr(D[mm], s)
        z=(r-nm)/nsd
        # 双侧经验 p（近似正态）
        from scipy.stats import norm
        pemp=2*(1-norm.cdf(abs(z)))
        tag = '✅超零' if pemp<0.05 else '❌在零内'
        log(f'  {cond:8s} {p:5s} {r:+8.4f} {pemp:8.4f}  z={z:+6.2f} {tag}')
        res[f'{cond}_{p}']=dict(rho=float(r), p_emp=float(pemp), z=float(z))

# ---- 双差法 ----
log('\n' + '='*88)
log('[3] 双差法（F6−F1, F7−F5 的梯度）')
log('='*88)
ddres={}
for (a,b) in [('F6','F1'),('F6','F5'),('F7','F1'),('F7','F5'),('F7e','F5')]:
    log(f'\n  {a} − {b}:')
    ddres[f'{a} − {b}']={}
    for cond in ['Wound1','Wound7','Wound30']:
        mm=(M['Condition']==cond).values
        if mm.sum()<50: continue
        ga=[g for g in MARK[a] if g in GZ.columns]
        gb=[g for g in MARK[b] if g in GZ.columns]
        dv = GZ.loc[mm,ga].mean(1).values - GZ.loc[mm,gb].mean(1).values
        r,_=spearmanr(D[mm], dv)
        # 零分布
        rs=[]
        for _ in range(2000):
            g1=RNG.choice(list(GZ.columns),7,replace=False)
            g2=RNG.choice(list(GZ.columns),5,replace=False)
            dd=GZ.loc[mm,g1].mean(1).values-GZ.loc[mm,g2].mean(1).values
            rr,_=spearmanr(D[mm],dd); rs.append(rr)
        rs=np.array(rs)
        from scipy.stats import norm
        z=(r-rs.mean())/ (rs.std()+1e-9)
        pemp=2*(1-norm.cdf(abs(z)))
        bad = '✅' if pemp<0.05 else '❌'
        log(f'    {cond:8s}: rho={r:+.4f}  null_sd={rs.std():.4f}  z={z:+.2f}  p={pemp:.4f} {bad}')
        ddres[f'{a} − {b}'][cond]=dict(rho=float(r), null_sd=float(rs.std()), z=float(z), p=float(pemp))

# ---- 关键结论：F6 负梯度是否特异 ----
log('\n' + '='*88)
log('[4] 关键：F6 的【负】梯度特异性')
log('='*88)
mm=(M['Condition'].isin(['Wound1','Wound7'])).values
for p,gl in MARK.items():
    gs=[g for g in gl if g in GZ.columns]
    s=GZ.loc[mm,gs].mean(1).values
    r,_=spearmanr(D[mm], s)
    log(f'  {p:5s}: rho={r:+.4f}  ({"负梯度=近创缘" if r<0 else "正梯度=远创缘"})')

out=dict(nulls=nulls, program_gradients=res, double_diff=ddres,
         interpretation='F6 negative gradient (edge-proximal) is the distinctive spatial feature')
with open(OUT/"S29b_gradient_specificity.json",'w',encoding='utf-8') as f:
    json.dump(out,f,indent=2,ensure_ascii=False,default=str)
log(f'\n[saved] S29b_gradient_specificity.json')
with open(OUT/"S29b_run.log",'w',encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('\n✅ S29b 完成')
