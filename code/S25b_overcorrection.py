# -*- coding: utf-8 -*-
"""
S25b: 残差化过度校正的诊断 + 最终稳定的替代方案

S25 暴露的新问题：
  线性残差化把 F5_orig vs UMI 从 +0.03 变成 -0.747（符号翻转 + 幅度暴涨）
  → 这是【过度校正】：β 被少数高杠杆点拉大，残差被反向污染
  progress vs UMI 也从 +0.613 翻成 -0.392

本脚本：
  D1. 诊断过度校正：残差 vs UMI 的关系是否非线性/异方差
  D2. 稳健替代方案对比：
      (a) 线性残差化（S25，已证过校）
      (b) 秩残差化（对 log UMI 的秩回归取残差）
      (c) 分层内标准化（按 UMI 十分位做组内 z-score）
      (d) 深度匹配（丢弃而非建模）
      (e) 比例化：用 F7/(F6+|F7|) 等比值（对全局缩放不变）
  D3. 对每种方案报告：Healer-NonHealer obs/p + 残余深度相关
  D4. 最终裁决：是否存在任何深度稳健的证据支持"NonHealer 停滞"
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.stats import spearmanr, rankdata

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
DATA = ROOT / "10_DFU_Dataset" / "dfu_outcome_dataset" / "GSE165816_sc.h5ad"
RNG = np.random.default_rng(20261006)

LOG = []
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOG.append(s)

log('=' * 84)
log('S25b: 过度校正诊断 + 稳健方案对比')
log('=' * 84)

a = ad.read_h5ad(DATA)
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
genes = a.var_names.astype(str); gi = {g:i for i,g in enumerate(genes)}
obs = a.obs.copy()

FIB = ['DCN','LUM','PDGFRA','COL1A1','COL1A2','FAP','THY1']
fib_cols = [gi[g] for g in FIB if g in gi]
sub = np.asarray(X[:, fib_cols].todense())
subz = (sub - sub.mean(0)) / (sub.std(0)+1e-9)
fs = subz.mean(1)
is_fib = fs >= np.quantile(fs, .75)
tc = obs['tissue_class'].astype(str).values
sel = np.where((tc=='wound_bed') & is_fib)[0]
Xw = X[sel]
gm = np.asarray(Xw.mean(0)).ravel()
gs2 = np.sqrt(np.asarray(Xw.power(2).mean(0)).ravel() - gm**2) + 1e-9
Xz = (np.asarray(Xw.todense()) - gm) / gs2
umi = np.asarray(Xw.sum(1)).ravel()
grp = obs['group'].astype(str).values[sel]
pat = obs['subject_hint'].astype(str).values[sel]
logu = np.log1p(umi)

def score(cols):
    return Xz[:, [gi[g] for g in cols if g in gi]].mean(1)

F6 = score(['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'])
F7 = score(['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'])
F1 = score(['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'])
F5 = score(['SCN7A','SOX10','PLP1','NGFR','S100B'])

# =========================================================
# D1. 过度校正诊断
# =========================================================
log('\n' + '=' * 84)
log('D1. 过度校正诊断')
log('=' * 84)

# 检查 F5 vs UMI 的原关系
log('\n  [F5_orig vs log(UMI)] 分箱均值（看是否单调/线性）:')
q = pd.qcut(logu, 5, labels=False, duplicates='drop')
for b in sorted(np.unique(q)):
    m = q == b
    log(f'    十分位{b+1}: logUMI∈[{logu[m].min():.2f},{logu[m].max():.2f}]  '
        f'F5均值={F5[m].mean():+.4f}  F6均值={F6[m].mean():+.4f}  F7均值={F7[m].mean():+.4f}')

# 线性残差化
lc = logu - logu.mean(); den = (lc**2).sum()+1e-12
def lin_resid(s):
    beta = (lc*(s-s.mean())).sum()/den
    return s - beta*lc, beta
F5_lin, b_F5 = lin_resid(F5)
F6_lin, b_F6 = lin_resid(F6)
F7_lin, b_F7 = lin_resid(F7)
log(f'\n  线性残差化斜率: F5 beta={b_F5:+.4f}, F6 beta={b_F6:+.4f}, F7 beta={b_F7:+.4f}')
log(f'  残差 vs UMI: F5 rho={spearmanr(F5_lin,umi)[0]:+.4f}  F6 rho={spearmanr(F6_lin,umi)[0]:+.4f}  F7 rho={spearmanr(F7_lin,umi)[0]:+.4f}')
log('  -> 符号翻转 + |rho| 放大 = 过度校正（线性模型外推）')

# =========================================================
# D2/D3. 稳健方案对比（患者级）
# =========================================================
log('\n' + '=' * 84)
log('D2/D3. 稳健方案对比（患者级 Healer vs NonHealer）')
log('=' * 84)

def patient_table(scores_dict, min_cells=30):
    rows = []
    for p in np.unique(pat):
        m = pat == p
        if m.sum() < min_cells: continue
        r = {'patient': p, 'group': pd.Series(grp[m]).mode()[0],
             'n_cells': int(m.sum()), 'umi': float(umi[m].mean())}
        for k, s in scores_dict.items():
            r[k] = float(s[m].mean())
        rows.append(r)
    return pd.DataFrame(rows)

def perm_test(a, b, nperm=10000, seed=11):
    rng = np.random.default_rng(seed)
    obs = a.mean()-b.mean()
    allv = np.concatenate([a,b]); na=len(a); cnt=0
    for _ in range(nperm):
        pm = rng.permutation(allv)
        if (pm[:na].mean()-pm[na:].mean()) >= obs: cnt+=1
    return float(obs), float((cnt+1)/(nperm+1))

# 方案 (a) 线性残差化
# 方案 (b) 秩残差化：对 rank(logUMI) 回归 rank(score)
def rank_resid(s):
    rs = rankdata(s); ru = rankdata(logu)
    ru_c = ru - ru.mean()
    beta = (ru_c*(rs-rs.mean())).sum()/((ru_c**2).sum()+1e-12)
    return rs - beta*ru_c

# 方案 (c) 分层内标准化：按 UMI 五分位做组内 z
def stratum_z(s, nb=5):
    q = pd.qcut(logu, nb, labels=False, duplicates='drop')
    out = np.zeros_like(s)
    for b in np.unique(q):
        m = q==b
        out[m] = (s[m]-s[m].mean())/(s[m].std()+1e-9)
    return out

# 方案 (e) 比值（对全局缩放不变）
# 注意 F6/F7 可为负 -> 用差分再除以和
ratio = (F7 - F6) / (np.abs(F6) + np.abs(F7) + 1e-6)

SCHEMES = {
    'a_lin_resid': {'F6': F6_lin, 'F7': F7_lin},
    'b_rank_resid': {'F6': rank_resid(F6), 'F7': rank_resid(F7)},
    'c_stratum_z': {'F6': stratum_z(F6), 'F7': stratum_z(F7)},
    'e_ratio_F7F6': {'F6': -ratio*0.5, 'F7': ratio*0.5},   # 构造二维以便共用 progress
    'raw': {'F6': F6, 'F7': F7},
}

# 轨迹
dt = pd.read_csv(OUT / "S20_donor_trajectory.csv")
dt['tp'] = pd.Categorical(dt['tp'], categories=['D0','D1','D7','D30'], ordered=True)
dt = dt.sort_values(['donor','tp']).reset_index(drop=True)
donors = sorted(dt.donor.unique()); TP=['D0','D1','D7','D30']
S = {d: dt[dt.donor==d].set_index('tp').loc[TP][['F6','F7']].values for d in donors}
mean_traj = np.mean([S[d] for d in donors], axis=0)

def project(p, poly):
    seglen = np.linalg.norm(np.diff(poly,axis=0),axis=1)
    cum = np.concatenate([[0],np.cumsum(seglen)]); total=cum[-1]
    best=(np.inf,None)
    for k in range(len(poly)-1):
        A,B=poly[k],poly[k+1]; AB=B-A; L2=AB@AB
        t=0.0 if L2<1e-12 else np.clip((p-A)@AB/L2,0,1)
        pr=A+t*AB; d=np.linalg.norm(p-pr); prog=(cum[k]+t*seglen[k])/total
        if d<best[0]: best=(d,prog)
    return best

log(f'\n  {"方案":16s} {"Healer prog":>11s} {"NonHlr prog":>11s} {"obs":>8s} {"p":>8s} {"resid ρ(UMI)":>13s} {"F7-F6 p":>9s}')
results = {}
for name, sc in SCHEMES.items():
    P = patient_table(sc)
    prog = np.array([project(p, mean_traj)[1] for p in P[['F6','F7']].values])
    P['progress'] = prog
    h = P[P['group']=='Healer']['progress'].values
    n = P[P['group']=='NonHealer']['progress'].values
    if len(h)<2 or len(n)<2:
        log(f'  {name:16s} 样本不足'); continue
    o, p = perm_test(h, n)
    rho_u = spearmanr(prog, P['umi'].values)[0]
    # F7-F6 差值检验
    P['d'] = P['F7'] - P['F6']
    hd = P[P['group']=='Healer']['d'].values
    nd = P[P['group']=='NonHealer']['d'].values
    od, pd_ = perm_test(hd, nd)
    log(f'  {name:16s} {h.mean():11.4f} {n.mean():11.4f} {o:+8.4f} {p:8.4f} {rho_u:+13.4f} {pd_:9.4f}')
    results[name] = dict(healer=h.mean(), nonhealer=n.mean(), obs=o, p=p, umi_rho=float(rho_u), d_p=pd_, d_obs=od)

# =========================================================
# D4. 深度匹配（丢弃法）在各方案下
# =========================================================
log('\n' + '=' * 84)
log('D4. 深度匹配丢弃法')
log('=' * 84)
# 全局：Healer 和 NonHealer 的 UMI 重叠区
h_all = np.where(grp=='Healer')[0]; n_all = np.where(grp=='NonHealer')[0]
lo = max(np.percentile(umi[h_all],10), np.percentile(umi[n_all],10))
hi = min(np.percentile(umi[h_all],90), np.percentile(umi[n_all],90))
log(f'  UMI 重叠区间 [{lo:.0f}, {hi:.0f}]')
m_keep = (umi>=lo)&(umi<=hi)
log(f'  保留细胞 {m_keep.sum():,}/{len(umi):,}')
# 用匹配子集的细胞重建 F6/F7（原始 z-score，但在匹配深度内）
F6m = F6.copy(); F7m = F7.copy()  # 已在全量上算 z-score
Pm = patient_table({'F6':F6m,'F7':F7m})  # 注意：这里未按细胞过滤
# 正确做法：对匹配细胞重算门控与 z
sel2 = sel[m_keep]  # 细胞索引（在 sel 坐标系内）
# 简化：直接用细胞级得分在匹配子集上取均值
def patient_table_masked(F6s, F7s, mask):
    rows=[]
    for p in np.unique(pat):
        m = (pat==p)&mask
        if m.sum()<30: continue
        rows.append({'patient':p,'group':pd.Series(grp[m]).mode()[0],
                     'n_cells':int(m.sum()),'umi':float(umi[m].mean()),
                     'F6':float(F6s[m].mean()),'F7':float(F7s[m].mean())})
    return pd.DataFrame(rows)
Pmk = patient_table_masked(F6,F7,m_keep)
prog_m = np.array([project(p,mean_traj)[1] for p in Pmk[['F6','F7']].values])
Pmk['progress']=prog_m
h = Pmk[Pmk['group']=='Healer']['progress'].values
n = Pmk[Pmk['group']=='NonHealer']['progress'].values
log(f'  匹配后 Healer n={len(h)}, NonHealer n={len(n)}')
if len(h)>=2 and len(n)>=2:
    o,p = perm_test(h,n)
    log(f'  匹配后 Healer-NonHealer: obs={o:+.4f}, p={p:.4f}')
    results['d_matched'] = dict(obs=o, p=p, n_h=len(h), n_nh=len(n))

# =========================================================
# 最终裁决
# =========================================================
log('\n' + '=' * 84)
log('最终裁决（S25b）')
log('=' * 84)
any_sig = [k for k,v in results.items() if v.get('p',1) < 0.05]
log(f'\n  达到 p<0.05 的方案: {any_sig if any_sig else "无"}')
log(f'  各方案残留深度相关 (progress vs UMI):')
for k,v in results.items():
    rr = v.get("umi_rho")
    rr_s = f'{rr:+.4f}' if rr is not None else '   n/a '
    log(f'    {k:16s} rho={rr_s}  p={v.get("p"):.4f}')
log('\n  ⚠️ 若"达到显著的方案"同时残留 |rho|>0.3 的深度相关，则该显著性不可信。')

out = dict(results=results, any_significant=any_sig)
with open(OUT / "S25b_scheme_comparison.json",'w',encoding='utf-8') as f:
    json.dump(out,f,indent=2,ensure_ascii=False,default=str)
log(f'\n[saved] S25b_scheme_comparison.json')
with open(OUT / "S25b_run.log",'w',encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('\n✅ S25b 完成')
