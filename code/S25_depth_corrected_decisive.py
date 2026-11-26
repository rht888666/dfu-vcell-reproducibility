# -*- coding: utf-8 -*-
"""
S25: 深度校正后的决定性重分析 —— "NonHealer 停滞"假说的最终裁决

背景（诚实链条）：
  S24b: 用 S10_patient_scores.csv（未做深度校正）得到 Healer vs NonHealer progress p=0.0245
  S24c: 发现 progress 与 UMI rho=+0.613 -> 深度混淆
  S24d: 裁决 DEPTH_CONFOUNDED；根因 = S10 输出的得分是【原始 z-score 均值】，未残差化
  S25 : 从 h5ad 重建【细胞级深度残差化】的得分 -> 重新做全部检验

方法：
  1. 重读 GSE165816_sc.h5ad，wound_bed 成纤维门控（与 S10 完全一致）
  2. 对每个程序得分，在【细胞级】对 log(UMI) 做残差化（与 S14 一致）
  3. 患者级 pseudobulk（残差化后）
  4. 重新做：
     R1. 分组均值表（残差化）
     R2. Healer vs NonHealer 置换检验
     R3. 轨迹 progress 投影 + 置换
     R4. 深度匹配子集
     R5. 阴性对照 F1/F5
  5. 与未校正版本并排对比，诚实给出"哪些结论存活、哪些死亡"

判据：
  - 若残差化后 Healer>NonHealer 仍显著 -> 信号是真实的，S24b 是"对的结论、脆弱的路径"
  - 若消失 -> 该结论死亡，必须从论文中撤下（第 9 次自我推翻的确认）
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.stats import spearmanr, mannwhitneyu

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
log('S25: 深度校正后的决定性重分析')
log('=' * 84)

PROG = {
    'F6_orig': ['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
    'F7_orig': ['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
    'F7_ecm':  ['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
    'F1_orig': ['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
    'F5_orig': ['SCN7A','SOX10','PLP1','NGFR','S100B'],
}

log(f'\n[load] {DATA.name}')
a = ad.read_h5ad(DATA)
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
genes = a.var_names.astype(str)
gi = {g: i for i, g in enumerate(genes)}
obs = a.obs.copy()
log(f'  {X.shape[0]:,} 细胞 × {X.shape[1]:,} 基因')

# 成纤维门控（与 S10 一致）
FIB = ['DCN','LUM','PDGFRA','COL1A1','COL1A2','FAP','THY1']
fib_cols = [gi[g] for g in FIB if g in gi]
sub = np.asarray(X[:, fib_cols].todense())
subz = (sub - sub.mean(0)) / (sub.std(0) + 1e-9)
fs = subz.mean(1)
is_fib = fs >= np.quantile(fs, .75)

tc = obs['tissue_class'].astype(str).values
wb = tc == 'wound_bed'
sel = np.where(wb & is_fib)[0]
log(f'  wound_bed 成纤维: {len(sel):,}')

Xw = X[sel]
gm = np.asarray(Xw.mean(0)).ravel()
gs2 = np.sqrt(np.asarray(Xw.power(2).mean(0)).ravel() - gm**2) + 1e-9
Xz = (np.asarray(Xw.todense()) - gm) / gs2
cell_umi = np.asarray(Xw.sum(1)).ravel()
grp = obs['group'].astype(str).values[sel]
pat = obs['subject_hint'].astype(str).values[sel]
logu = np.log1p(cell_umi)

# ---- 计算原始得分 与 深度残差化得分 ----
SC_raw, SC_res = {}, {}
log('\n[score] 计算原始 vs 深度残差化得分')
log('  程序          rho(raw,UMI)   rho(resid,UMI)')
for k, vs in PROG.items():
    cols = [gi[g] for g in vs if g in gi]
    if not cols:
        continue
    s = Xz[:, cols].mean(1)
    SC_raw[k] = s
    # 细胞级对 log UMI 残差化
    lc = logu - logu.mean()
    den = (lc**2).sum() + 1e-12
    beta = (lc * (s - s.mean())).sum() / den
    SC_res[k] = s - beta * lc
    r_raw = spearmanr(s, cell_umi)[0]
    r_res = spearmanr(SC_res[k], cell_umi)[0]
    log(f'  {k:12s}  {r_raw:+.4f}       {r_res:+.4f}')

# ---- 患者级 pseudobulk ----
def build_patient_scores(SC, min_cells=30):
    rows = []
    for p in np.unique(pat):
        m = pat == p
        if m.sum() < min_cells:
            continue
        r = {'patient': p, 'group': pd.Series(grp[m]).mode()[0],
             'n_cells': int(m.sum()), 'umi': float(cell_umi[m].mean())}
        for k in SC:
            r[k] = float(SC[k][m].mean())
        rows.append(r)
    return pd.DataFrame(rows)

P_raw = build_patient_scores(SC_raw)
P_res = build_patient_scores(SC_res)
log(f'\n  患者数: raw={len(P_raw)}, resid={len(P_res)}')

# =========================================================
# R1. 分组均值对比
# =========================================================
log('\n' + '=' * 84)
log('R1. 分组均值对比（原始 vs 残差化）')
log('=' * 84)
ORDER = ['Healthy','DiabeticNoDFU','Healer','NonHealer']
for name, P in [('RAW(未校正)', P_raw), ('RESID(深度残差化)', P_res)]:
    log(f'\n  [{name}]')
    t = P.groupby('group')[['F6_orig','F7_orig','F7_ecm','F1_orig','F5_orig']].mean().reindex(ORDER)
    log(t.round(4).to_string())

# =========================================================
# R2. Healer vs NonHealer 检验（每个程序）
# =========================================================
log('\n' + '=' * 84)
log('R2. Healer vs NonHealer 置换检验（原始 vs 残差化）')
log('=' * 84)

def perm_test(a, b, nperm=10000, seed=11, alternative='greater'):
    rng = np.random.default_rng(seed)
    if alternative == 'greater':
        obs = a.mean() - b.mean()
    else:
        obs = abs(a.mean() - b.mean())
    allv = np.concatenate([a, b]); na = len(a); cnt = 0
    for _ in range(nperm):
        pm = rng.permutation(allv)
        d = pm[:na].mean() - pm[na:].mean()
        if alternative == 'greater':
            if d >= obs: cnt += 1
        else:
            if abs(d) >= obs: cnt += 1
    return float(obs), float((cnt + 1) / (nperm + 1))

log(f'\n  {"程序":12s} {"RAW obs":>10s} {"RAW p":>8s} | {"RESID obs":>10s} {"RESID p":>8s}')
res_R2 = {}
for k in ['F6_orig','F7_orig','F7_ecm','F1_orig','F5_orig']:
    if k not in P_raw.columns: continue
    h_r = P_raw[P_raw['group']=='Healer'][k].values
    n_r = P_raw[P_raw['group']=='NonHealer'][k].values
    h_s = P_res[P_res['group']=='Healer'][k].values
    n_s = P_res[P_res['group']=='NonHealer'][k].values
    or_, pr = perm_test(h_r, n_r)
    os_, ps = perm_test(h_s, n_s)
    log(f'  {k:12s} {or_:+10.4f} {pr:8.4f} | {os_:+10.4f} {ps:8.4f}')
    res_R2[k] = dict(raw_obs=or_, raw_p=pr, resid_obs=os_, resid_p=ps)

# =========================================================
# R3. 轨迹 progress 投影（残差化）
# =========================================================
log('\n' + '=' * 84)
log('R3. 轨迹 progress 投影（用残差化得分）')
log('=' * 84)

dt = pd.read_csv(OUT / "S20_donor_trajectory.csv")
dt['tp'] = pd.Categorical(dt['tp'], categories=['D0','D1','D7','D30'], ordered=True)
dt = dt.sort_values(['donor','tp']).reset_index(drop=True)
donors = sorted(dt.donor.unique()); TP = ['D0','D1','D7','D30']
S = {d: dt[dt.donor==d].set_index('tp').loc[TP][['F6','F7']].values for d in donors}
mean_traj = np.mean([S[d] for d in donors], axis=0)

def project(p, poly):
    seglen = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    cum = np.concatenate([[0], np.cumsum(seglen)]); total = cum[-1]
    best = (np.inf, None)
    for k in range(len(poly)-1):
        A, B = poly[k], poly[k+1]; AB = B-A; L2 = AB@AB
        t = 0.0 if L2 < 1e-12 else np.clip((p-A)@AB/L2, 0, 1)
        pr = A + t*AB; dist = np.linalg.norm(p-pr)
        prog = (cum[k] + t*seglen[k])/total
        if dist < best[0]: best = (dist, prog)
    return best

for name, P in [('RAW', P_raw), ('RESID', P_res)]:
    prog = np.array([project(p, mean_traj)[1] for p in P[['F6_orig','F7_orig']].values])
    P['progress'] = prog
    log(f'\n  [{name}] 分组 progress:')
    for g in ORDER:
        s = P[P['group']==g]
        if len(s)==0: continue
        log(f'    {g:14s} n={len(s):2d} progress={s["progress"].mean():.4f}')
    h = P[P['group']=='Healer']['progress'].values
    n = P[P['group']=='NonHealer']['progress'].values
    o, p = perm_test(h, n)
    rho_umi = spearmanr(prog, P['umi'].values)[0]
    log(f'    Healer-NonHealer: obs={o:+.4f}, p={p:.4f}   |  progress vs UMI rho={rho_umi:+.4f}')
    res_R2[f'progress_{name}'] = dict(obs=o, p=p, umi_rho=float(rho_umi))

# =========================================================
# R4. 深度匹配子集（残差化）
# =========================================================
log('\n' + '=' * 84)
log('R4. Healer vs NonHealer（残差化 + 深度细分层）')
log('=' * 84)
P_res_sorted = P_res.sort_values('umi')
med = P_res['umi'].median()
for st, mask in [('high UMI', P_res['umi']>=med), ('low UMI', P_res['umi']<med)]:
    sub = P_res[mask]
    h = sub[sub['group']=='Healer']['progress'].values
    n = sub[sub['group']=='NonHealer']['progress'].values
    if len(h)>=2 and len(n)>=2:
        o, p = perm_test(h, n)
        log(f'  {st:9s}: H n={len(h)}, NH n={len(n)}, obs={o:+.4f}, p={p:.4f}')

# =========================================================
# 最终裁决
# =========================================================
log('\n' + '=' * 84)
log('最终裁决')
log('=' * 84)
pr = res_R2.get('progress_RESID', {})
if pr.get('p', 1) < 0.05:
    verdict = 'SIGNAL_SURVIVES_DEPTH_CORRECTION'
else:
    verdict = 'SIGNAL_DIES_UNDER_DEPTH_CORRECTION'
log(f'\n  残差化后 Healer-NonHealer progress: obs={pr.get("obs"):+.4f}, p={pr.get("p"):.4f}')
log(f'  -> {verdict}')
log(f'\n  未校正 progress vs UMI rho = {res_R2.get("progress_RAW",{}).get("umi_rho"):+.4f}')
log(f'  残差化 progress vs UMI rho = {pr.get("umi_rho"):+.4f}')

out = dict(
    R1=dict(
        raw=P_raw.groupby('group')[['F6_orig','F7_orig','F7_ecm']].mean().reindex(ORDER).round(4).to_dict(),
        resid=P_res.groupby('group')[['F6_orig','F7_orig','F7_ecm']].mean().reindex(ORDER).round(4).to_dict(),
    ),
    R2=res_R2,
    verdict=verdict,
)
with open(OUT / "S25_depth_corrected.json", 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False, default=str)
P_res.to_csv(OUT / "S25_patient_scores_depthresid.csv", index=False)
log(f'\n[saved] S25_depth_corrected.json + S25_patient_scores_depthresid.csv')
with open(OUT / "S25_run.log", 'w', encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('\n✅ S25 完成')
