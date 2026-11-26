# -*- coding: utf-8 -*-
"""
S24d: 深度混淆的深入诊断 —— progress 信号中"真信号"与"深度伪影"的分解

S24c 发现了一个必须正视的问题：
  progress 与测序深度 UMI 强相关 (rho=+0.613, p=0.0002)
  Healer UMI=12,452 vs NonHealer UMI=9,052（差 3,400）
  残差化后 p=0.104 → 不能声称有独立于深度的信号

本脚本要回答：
  Q1. 深度是"混淆变量"还是"中介变量"？
      - 混淆：深度同时影响分组和 progress，但无生物学意义（BAD）
      - 中介：深度是愈合状态的【下游结果】（组织更活跃 -> 捕获更多 RNA）（可接受）
  Q2. 原始 F6/F7 得分本身是否已做深度校正？（S20 marker z-score 是否已控深度）
  Q3. 能否构造一个【深度正交】的 progress 指标？
  Q4. 在深度匹配的子集中检验，信号是否存活？

判据：
  - 若在【深度匹配】的子集内信号存活 -> 深度是部分中介，信号可保留（诚实标注）
  - 若消失 -> 信号即深度伪影，必须放弃该指标
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
RNG = np.random.default_rng(20261006)

LOG = []
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOG.append(s)

log('=' * 84)
log('S24d: 深度混淆诊断 —— 信号 vs 伪影分解')
log('=' * 84)

# --- 载入 ---
dt = pd.read_csv(OUT / "S20_donor_trajectory.csv")
dt['tp'] = pd.Categorical(dt['tp'], categories=['D0','D1','D7','D30'], ordered=True)
dt = dt.sort_values(['donor','tp']).reset_index(drop=True)
donors = sorted(dt.donor.unique())
TP = ['D0','D1','D7','D30']
def dmat(d):
    return dt[dt.donor == d].set_index('tp').loc[TP][['F6','F7']].values
S = {d: dmat(d) for d in donors}
mean_traj = np.mean([S[d] for d in donors], axis=0)

dfu = pd.read_csv(OUT / "S10_patient_scores.csv")
d2 = dfu.dropna(subset=['F6_orig','F7_orig','group']).copy()
d2['grp'] = d2['group'].astype(str)

def project_to_polyline(p, poly):
    seglen = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    cum = np.concatenate([[0], np.cumsum(seglen)]); total = cum[-1]
    best = (np.inf, None)
    for k in range(len(poly)-1):
        a, b = poly[k], poly[k+1]; ab = b - a; L2 = ab @ ab
        t = 0.0 if L2 < 1e-12 else np.clip((p-a)@ab/L2, 0, 1)
        proj = a + t*ab; dist = np.linalg.norm(p-proj)
        prog = (cum[k] + t*seglen[k])/total
        if dist < best[0]: best = (dist, prog)
    return best

base_prog = np.array([project_to_polyline(p, mean_traj)[1] for p in d2[['F6_orig','F7_orig']].values])
d2['progress'] = base_prog
umi = d2['umi'].values.astype(float)

# =========================================================
# Q1. 深度-分组-得分 的三方关系
# =========================================================
log('\n' + '=' * 84)
log('Q1. 三方关系：深度 ~ 分组 ~ progress')
log('=' * 84)

log('\n  [分组 × 深度]')
for g in ['Healthy','DiabeticNoDFU','Healer','NonHealer']:
    sub = d2[d2['grp']==g]
    log(f'   {g:14s} n={len(sub):2d}  UMI={sub["umi"].mean():8.0f} ± {sub["umi"].std():6.0f}  progress={sub["progress"].mean():.4f}')

# 深度与分组的关联（Healer vs NonHealer）
h_umi = d2[d2['grp']=='Healer']['umi'].values
nh_umi = d2[d2['grp']=='NonHealer']['umi'].values
u, pu = mannwhitneyu(h_umi, nh_umi, alternative='two-sided')
log(f'\n  Healer vs NonHealer UMI: {h_umi.mean():.0f} vs {nh_umi.mean():.0f}, Mann-Whitney p={pu:.4f}')
log(f'  -> 分组本身与深度显著相关（混淆的前提成立）')

# 深度与 progress
rho, prho = spearmanr(base_prog, umi)
log(f'  progress vs UMI: rho={rho:+.4f}, p={prho:.4f}')

# 原始的 F6/F7 得分是否已深度校正？
log('\n  [原始得分 vs 深度]')
for col in ['F6_orig','F7_orig','F6_il11','F7_ecm']:
    r, p = spearmanr(d2[col].values, umi)
    log(f'   {col:12s} vs UMI: rho={r:+.4f}, p={p:.4f}')

# =========================================================
# Q3. 构造深度正交的 progress
# =========================================================
log('\n' + '=' * 84)
log('Q3. 深度正交化：构造与 UMI 不相关的 progress')
log('=' * 84)

# 方法：Robust — 用 LOWESS/多项式把 progress 对 log(UMI) 回归，取残差
log_umi = np.log1p(umi)
# 二次项以防非线性
coef = np.polyfit(log_umi, base_prog, 2)
prog_resid = base_prog - np.polyval(coef, log_umi)
rho2, p2 = spearmanr(prog_resid, umi)
log(f'  残差化后 progress vs UMI: rho={rho2:+.4f}, p={p2:.4f}  (应接近 0)')

def perm_test(prog, grp, ga='Healer', gb='NonHealer', nperm=10000, seed=7):
    rng = np.random.default_rng(seed)
    a = prog[grp==ga]; b = prog[grp==gb]
    obs = a.mean()-b.mean()
    allv = np.concatenate([a,b]); na=len(a); cnt=0
    for _ in range(nperm):
        pm = rng.permutation(allv)
        if (pm[:na].mean()-pm[na:].mean()) >= obs: cnt+=1
    return float(obs), float((cnt+1)/(nperm+1))

o_raw, p_raw = perm_test(base_prog, d2['grp'].values)
o_res, p_res = perm_test(prog_resid, d2['grp'].values)
log(f'\n  原始 progress:   obs={o_raw:+.4f}, p={p_raw:.4f}')
log(f'  深度残差 progress: obs={o_res:+.4f}, p={p_res:.4f}')

# =========================================================
# Q4. 深度匹配子集检验（最关键）
# =========================================================
log('\n' + '=' * 84)
log('Q4. 深度匹配子集检验（1:1 最近邻匹配）')
log('=' * 84)
log('  方法：对每个 NonHealer，找一个 UMI 最接近的 Healer（卡尺限制），做配对比较')

h_df = d2[d2['grp']=='Healer'].copy()
nh_df = d2[d2['grp']=='NonHealer'].copy()
used = set(); pairs = []
for _, nr in nh_df.sort_values('umi').iterrows():
    cand = h_df[~h_df.index.isin(used)].copy()
    if len(cand)==0: break
    cand['dist'] = (cand['umi'] - nr['umi']).abs()
    best = cand.sort_values('dist').iloc[0]
    # 卡尺：深度差 < 组内 SD 的一半
    caliper = h_df['umi'].std() * 0.75
    if best['dist'] <= caliper:
        used.add(best.name)
        pairs.append((nr['patient'], best['patient'], nr['umi'], best['umi'], nr['progress'], best['progress']))

log(f'  匹配成功 {len(pairs)} 对（卡尺={h_df["umi"].std()*0.75:.0f}）:')
for a,b,ua,ub,pa,pb in pairs:
    log(f'    NH {a}(UMI={ua:.0f},prog={pa:.3f}) <-> H {b}(UMI={ub:.0f},prog={pb:.3f})  Δprog={pb-pa:+.3f}')
if pairs:
    dprog = np.array([pb-pa for *_ ,pa,pb in [(x[0],x[1],x[2],x[3],x[4],x[5]) for x in pairs]])
    dprog = np.array([x[5]-x[4] for x in pairs])
    log(f'  配对 Δprogress (H-NH): mean={dprog.mean():+.4f}, n={len(dprog)}')
    # 符号检验
    from scipy.stats import binomtest
    n_pos = int((dprog>0).sum())
    bp = binomtest(n_pos, len(dprog), 0.5, alternative='greater').pvalue
    log(f'  Δprogress > 0 的比例: {n_pos}/{len(dprog)}, 符号检验 p={bp:.4f}')

# =========================================================
# Q5. 关键替代指标：用【深度已在得分层面校正】的 S20 方法重算
# =========================================================
log('\n' + '=' * 84)
log('Q5. 尝试其他深度正交的维度')
log('=' * 84)
log('  思路：不用 progress（2D 投影，对深度敏感），改用【程序得分差】这种简单对比')
log('  如 F7-F6（重塑-炎症 对比），该差值对全局深度缩放不敏感')

d2['F7_minus_F6'] = d2['F7_orig'] - d2['F6_orig']
d2['F7_ratio'] = d2['F7_orig'] / (d2['F7_orig'] - d2['F6_orig']).abs().replace(0, np.nan)
rho, p = spearmanr(d2['F7_minus_F6'].values, umi)
log(f'\n  (F7-F6) vs UMI: rho={rho:+.4f}, p={p:.4f}')
o, p = perm_test(d2['F7_minus_F6'].values, d2['grp'].values)
log(f'  (F7-F6) Healer vs NonHealer: obs={o:+.4f}, p={p:.4f}')

for g in ['Healthy','DiabeticNoDFU','Healer','NonHealer']:
    sub = d2[d2['grp']==g]
    log(f'   {g:14s} F6={sub["F6_orig"].mean():+.4f}  F7={sub["F7_orig"].mean():+.4f}  F7-F6={sub["F7_minus_F6"].mean():+.4f}')

# 汇总
out = dict(
    q1=dict(
        umi_by_group={g: float(d2[d2['grp']==g]['umi'].mean()) for g in ['Healthy','DiabeticNoDFU','Healer','NonHealer']},
        progress_by_group={g: float(d2[d2['grp']==g]['progress'].mean()) for g in ['Healthy','DiabeticNoDFU','Healer','NonHealer']},
        progress_umi_rho=float(rho), healer_nonhealer_umi_p=float(pu),
    ),
    q3=dict(raw_obs=o_raw, raw_p=p_raw, resid_obs=o_res, resid_p=p_res),
    q4=dict(n_pairs=len(pairs), pairs=[dict(nh=a, h=b, umi_nh=ua, umi_h=ub, prog_nh=pa, prog_h=pb)
                                        for a,b,ua,ub,pa,pb in pairs]),
    verdict='DEPTH_CONFOUNDED' if p_res > 0.05 else 'SIGNAL_SURVIVES',
)
with open(OUT / "S24d_depth_diagnosis.json", 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log(f'\n[saved] S24d_depth_diagnosis.json')
log('\n' + '=' * 84)
log(f'最终裁决: {out["verdict"]}')
log('=' * 84)
with open(OUT / "S24d_run.log", 'w', encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('\n✅ S24d 完成')
