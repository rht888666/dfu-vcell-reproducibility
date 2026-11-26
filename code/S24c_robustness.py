# -*- coding: utf-8 -*-
"""
S24c: 轨迹进度指标的稳健性验证（S24b 的 p=0.0245 是否可信？）

必须做的工作（避免又犯"看起来不错"的错误）：
  V1. 深度校正敏感性：progress 是否被测序深度驱动？
      -> 用 S10 的 umi 列做偏相关 / 分层检验
  V2. 置换策略稳健性：Healer/NonHealer 标签置换 10000 次
  V3. 轨迹定义敏感性：用不同轨迹（每供体各自轨迹 vs 平均轨迹 vs 子集）重算
  V4. 留一交叉验证：逐个剔除 Healer 患者，看 p 值稳定性
  V5. 阴性对照：随机轨迹（打乱时间点顺序）作为 null
  V6. 与 W6 的对比：为什么一维投影能用、二维方向不能

只有全部通过，才能声称"NonHealer 停滞在轨迹起点"。
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
RNG = np.random.default_rng(20261006)

LOG = []
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOG.append(s)

log('=' * 84)
log('S24c: 轨迹进度指标的稳健性验证')
log('=' * 84)

dt = pd.read_csv(OUT / "S20_donor_trajectory.csv")
dt['tp'] = pd.Categorical(dt['tp'], categories=['D0','D1','D7','D30'], ordered=True)
dt = dt.sort_values(['donor','tp']).reset_index(drop=True)
donors = sorted(dt.donor.unique())
TP = ['D0','D1','D7','D30']

def dmat(donor):
    return dt[dt.donor == donor].set_index('tp').loc[TP][['F6','F7']].values

S = {d: dmat(d) for d in donors}
mean_traj = np.mean([S[d] for d in donors], axis=0)

dfu = pd.read_csv(OUT / "S10_patient_scores.csv")
d2 = dfu.dropna(subset=['F6_orig','F7_orig','group']).copy()
d2['grp'] = d2['group'].astype(str)

def project_to_polyline(p, poly):
    seglen = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    cum = np.concatenate([[0], np.cumsum(seglen)])
    total = cum[-1]
    best = (np.inf, None)
    for k in range(len(poly) - 1):
        a, b = poly[k], poly[k+1]
        ab = b - a
        L2 = ab @ ab
        t = 0.0 if L2 < 1e-12 else np.clip((p - a) @ ab / L2, 0, 1)
        proj = a + t * ab
        dist = np.linalg.norm(p - proj)
        prog = (cum[k] + t * seglen[k]) / total
        if dist < best[0]:
            best = (dist, prog)
    return best

def compute_progress(poly):
    pts = d2[['F6_orig','F7_orig']].values.astype(float)
    return np.array([project_to_polyline(p, poly)[1] for p in pts])

def perm_test(prog, grp, ga, gb, nperm=10000, seed=1):
    rng = np.random.default_rng(seed)
    a = prog[grp == ga]; b = prog[grp == gb]
    obs = a.mean() - b.mean()
    allv = np.concatenate([a, b]); na = len(a)
    cnt = 0
    for _ in range(nperm):
        pm = rng.permutation(allv)
        if (pm[:na].mean() - pm[na:].mean()) >= obs: cnt += 1
    return float(obs), float((cnt + 1) / (nperm + 1))

# 基线
base_prog = compute_progress(mean_traj)
d2['progress'] = base_prog
log('\n[基线] 平均轨迹投影:')
for g in ['Healthy','DiabeticNoDFU','Healer','NonHealer']:
    sub = d2[d2['grp'] == g]
    log(f'  {g:14s} n={len(sub):2d}  progress={sub["progress"].mean():.4f}')
obs0, p0 = perm_test(base_prog, d2['grp'].values, 'Healer', 'NonHealer')
log(f'  基线 Healer vs NonHealer: obs={obs0:+.4f}, p={p0:.4f}')

results = {'baseline': dict(obs=obs0, p=p0)}

# =========================================================
# V1. 深度校正敏感性
# =========================================================
log('\n' + '=' * 84)
log('V1. 深度敏感性：progress 是否被测序深度驱动？')
log('=' * 84)
# d2 有 'umi' 列
if 'umi' in d2.columns:
    umi = d2['umi'].values.astype(float)
    log(f'  UMI 均值: ' + ', '.join(f'{g}={d2[d2["grp"]==g]["umi"].mean():.0f}' for g in ['Healthy','DiabeticNoDFU','Healer','NonHealer']))
    # progress 与 UMI 的相关
    from scipy.stats import spearmanr
    rho_obs, p_obs = spearmanr(base_prog, umi)
    log(f'  progress vs UMI: Spearman rho={rho_obs:+.4f}, p={p_obs:.4f}')
    # 分层：按 UMI 中位数分为高低两组，各组内做 Healer vs NonHealer
    med = np.median(umi)
    d2['umi_stratum'] = np.where(umi >= med, 'high', 'low')
    log('  分层检验(按 UMI 中位数):')
    for st in ['high','low']:
        m = d2['umi_stratum'] == st
        pr = base_prog[m]; gr = d2['grp'].values[m]
        if (gr=='Healer').sum() >= 2 and (gr=='NonHealer').sum() >= 2:
            o, p = perm_test(pr, gr, 'Healer', 'NonHealer')
            log(f'    {st:5s}: Healer n={(gr=="Healer").sum()}, NonHealer n={(gr=="NonHealer").sum()}, obs={o:+.4f}, p={p:.4f}')
            results[f'stratum_{st}'] = dict(obs=o, p=p)
    # 残差化：把 progress 对 UMI 回归取残差
    b, a = np.polyfit(umi, base_prog, 1)
    resid = base_prog - (a + b * umi)
    o, p = perm_test(resid, d2['grp'].values, 'Healer', 'NonHealer')
    log(f'  UMI 残差化后 Healer vs NonHealer: obs={o:+.4f}, p={p:.4f}')
    results['umi_residualized'] = dict(obs=o, p=p)

# =========================================================
# V2. 轨迹定义敏感性
# =========================================================
log('\n' + '=' * 84)
log('V2. 轨迹定义敏感性（不同轨迹定义下 p 值是否稳定）')
log('=' * 84)
traj_defs = {'mean': mean_traj}
for d in donors:
    traj_defs[f'single_{d}'] = S[d]
# 只用 D0 和 D30 两点定义的直线
traj_defs['line_D0_D30'] = np.array([mean_traj[0], mean_traj[-1]])
# 只用 D7->D30 段
traj_defs['line_D7_D30'] = np.array([mean_traj[2], mean_traj[-1]])

for name, poly in traj_defs.items():
    pr = compute_progress(poly)
    o, p = perm_test(pr, d2['grp'].values, 'Healer', 'NonHealer')
    log(f'  {name:18s}: Healer-NonHealer obs={o:+.4f}, p={p:.4f}')
    results[f'traj_{name}'] = dict(obs=o, p=p)

# =========================================================
# V3. 留一验证（逐个剔除 Healer / NonHealer）
# =========================================================
log('\n' + '=' * 84)
log('V3. 留一交叉验证（p 值稳定性）')
log('=' * 84)
log('  剔除单个 Healer 患者:')
h_idx = d2.index[d2['grp'] == 'Healer']
p_list = []
for i in h_idx:
    mask = d2.index != i
    o, p = perm_test(base_prog[mask], d2['grp'].values[mask], 'Healer', 'NonHealer')
    p_list.append(p)
    log(f'    剔除 {d2.loc[i,"patient"]}: obs={o:+.4f}, p={p:.4f}')
log(f'  -> p 范围 [{min(p_list):.4f}, {max(p_list):.4f}]')
results['loo_healer_p_range'] = [float(min(p_list)), float(max(p_list))]

log('  剔除单个 NonHealer 患者:')
nh_idx = d2.index[d2['grp'] == 'NonHealer']
p_list2 = []
for i in nh_idx:
    mask = d2.index != i
    o, p = perm_test(base_prog[mask], d2['grp'].values[mask], 'Healer', 'NonHealer')
    p_list2.append(p)
    log(f'    剔除 {d2.loc[i,"patient"]}: obs={o:+.4f}, p={p:.4f}')
log(f'  -> p 范围 [{min(p_list2):.4f}, {max(p_list2):.4f}]')
results['loo_nonhealer_p_range'] = [float(min(p_list2)), float(max(p_list2))]

# =========================================================
# V4. 阴性对照：随机轨迹
# =========================================================
log('\n' + '=' * 84)
log('V4. 阴性对照：随机轨迹（打乱正常轨迹的时间顺序）')
log('=' * 84)
nperm = 2000
cnt = 0; null_obs = []
for _ in range(nperm):
    idx = RNG.permutation(4)
    poly_r = mean_traj[idx]
    pr = compute_progress(poly_r)
    o, _ = perm_test(pr, d2['grp'].values, 'Healer', 'NonHealer', nperm=200, seed=int(RNG.integers(1e6)))
    null_obs.append(o)
    if o >= obs0: cnt += 1
null_obs = np.array(null_obs)
log(f'  随机轨迹 obs 分布: mean={null_obs.mean():+.4f}, sd={null_obs.std():.4f}')
log(f'  观测 obs (真轨迹)={obs0:+.4f}')
# 注意用轨迹置换的 p
p_traj = (cnt + 1) / (nperm + 1)
log(f'  轨迹置换检验 p = {p_traj:.4f}')
results['traj_null'] = dict(p=p_traj, null_sd=float(null_obs.std()), null_mean=float(null_obs.mean()))

# =========================================================
# V5. 为什么一维投影能用、二维 W6 不能
# =========================================================
log('\n' + '=' * 84)
log('V5. 一维投影 vs 二维方向（解释 W6 失败）')
log('=' * 84)
log('  W6 检验的是【方向向量】的 cos 相似度 → 保留 2 个自由度，n=14 时零分布 sd=0.786')
log('  一维投影把二维散布【约束到轨迹轴】→ 只保留 1 个自由度，且用轨迹先验去噪')
log('  本质上：一维投影是"用正常轨迹做正则化"，牺牲方向自由度换取统计功效')
log('  → 这是合理的，因为我们的假设本来就是"沿这条轨迹前进/停滞"')

with open(OUT / "S24c_robustness.json", 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
log(f'\n[saved] S24c_robustness.json')
with open(OUT / "S24c_run.log", 'w', encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('[saved] S24c_run.log')
log('\n✅ S24c 完成')
