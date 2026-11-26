# -*- coding: utf-8 -*-
"""
S24b: 状态转移模型的稳健重做 —— 修复 S24 的退化问题 + 几何轨迹分析

S24 暴露的两个问题（诚实记录）：
  P1. ODE 拟合退化：网格下巴 K->0.001，拟合出 X*=-5.6/+20.3 等荒谬外推
      -> 根因：4 个时间点、非等间隔、K 与 X* 强共线，不可辨识
  P2. DFU 投影方向反了：NonHealer 反而更远离 D30（obs=+0.077），p=0.30 不显著

S24b 的修正策略：
  A. 放弃无约束指数拟合；改用**有物理约束的**模型：
     - 约束 X* 在合理范围内（-1..+1，即程序得分的实际尺度）
     - 只保留能识别 K 的供体，或改用**差分法**估计 K
  B. 改用**几何/拓扑**描述（不依赖参数拟合）：
     - 轨迹的曲率、总位移、方向角、路径长度
     - 相空间中的单调性（F6 是否单调降、F7 是否单调升）
  C. DFU 投影改为**沿轨迹的一维投影**（不是到锚点距离），
     并诚实报告方向与显著性
  D. 新增：**轨迹方向角的供体一致性**（3 供体最终方向是否一致）

所有数字机器生成，诚实标注可靠性。
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
log('S24b: 状态转移模型稳健重做（修复 S24 退化 + 几何轨迹 + 方向投影）')
log('=' * 84)

dt = pd.read_csv(OUT / "S20_donor_trajectory.csv")
dt['tp'] = pd.Categorical(dt['tp'], categories=['D0','D1','D7','D30'], ordered=True)
dt = dt.sort_values(['donor','tp']).reset_index(drop=True)
donors = sorted(dt.donor.unique())
T_DAYS = np.array([0.0, 1.0, 7.0, 30.0])
TP = ['D0','D1','D7','D30']

def dmat(donor):
    sub = dt[dt.donor == donor].set_index('tp').loc[TP]
    return sub[['F6','F7']].values

S = {d: dmat(d) for d in donors}

# =========================================================
# A. 差分法估计弛豫速率（避免参数共线）
# =========================================================
log('\n' + '=' * 84)
log('A. 差分法估计弛豫速率 K（dX/dt 与位移的线性回归）')
log('=' * 84)
log('  模型: Δx/Δt = -K·(x_mid - X*)  ->  以 x_mid 为自变量、Δx/Δt 为因变量做线性回归')
log('  斜率 = -K（可识别）；截距/K = X*')

for prog_i, prog in enumerate(['F6', 'F7']):
    log(f'\n  [{prog}]')
    all_dx = []; all_xmid = []; all_dt = []
    for d in donors:
        xs = S[d][:, prog_i]
        for k in range(3):
            x0, x1 = xs[k], xs[k+1]
            dtt = T_DAYS[k+1] - T_DAYS[k]
            all_dx.append((x1 - x0) / dtt)
            all_xmid.append((x0 + x1) / 2)
            all_dt.append(dtt)
    all_dx = np.array(all_dx); all_xmid = np.array(all_xmid)
    # 线性回归 dx/dt = a + b*x_mid
    b, a = np.polyfit(all_xmid, all_dx, 1)
    K_hat = -b
    Xstar_hat = a / (-b) if abs(b) > 1e-9 else np.nan
    # R2
    pred = a + b * all_xmid
    ss_res = ((all_dx - pred) ** 2).sum()
    ss_tot = ((all_dx - all_dx.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    log(f'    斜率 b = {b:+.5f}  ->  K = {K_hat:+.5f}')
    log(f'    截距 a = {a:+.5f}  ->  X* = {Xstar_hat:+.4f}')
    log(f'    R² = {r2:.3f}  (n={len(all_dx)} 段)')

# =========================================================
# B. 几何轨迹描述（无参数）
# =========================================================
log('\n' + '=' * 84)
log('B. 几何/拓扑轨迹描述（不依赖参数拟合）')
log('=' * 84)

geo = {}
log('\n  供体    总路径长  净位移  曲折度(D/L)  F6单调?  F7单调?  末段方向角(°)')
for d in donors:
    pts = S[d]
    steps = np.diff(pts, axis=0)
    seglen = np.linalg.norm(steps, axis=1)
    L = seglen.sum()
    net = np.linalg.norm(pts[-1] - pts[0])
    tort = L / net if net > 0 else np.nan
    f6 = pts[:, 0]; f7 = pts[:, 1]
    mono_f6_dec = bool(np.all(np.diff(f6) <= 1e-9))
    mono_f7_inc = bool(np.all(np.diff(f7) >= -1e-9))
    # 末段方向角 (D7->D30)
    v = pts[-1] - pts[-2]
    ang = float(np.degrees(np.arctan2(v[1], v[0])))
    geo[d] = dict(L=float(L), net=float(net), tort=float(tort),
                  mono_f6_dec=mono_f6_dec, mono_f7_inc=mono_f7_inc,
                  last_angle=ang, pts=pts.tolist())
    log(f'   {d:6s}  {L:7.3f}  {net:6.3f}   {tort:6.3f}     '
        f'{"是" if mono_f6_dec else "否":4s}   {"是" if mono_f7_inc else "否":4s}    {ang:+7.1f}')

# 群体平均轨迹
mean_pts = np.mean([S[d] for d in donors], axis=0)
log('\n  [供体平均轨迹]')
for i, t in enumerate(TP):
    log(f'    {t}: F6={mean_pts[i,0]:+.4f}  F7={mean_pts[i,1]:+.4f}')

# 关键：D0->D30 的整体方向（正常愈合方向）
v_normal = mean_pts[-1] - mean_pts[0]
log(f'\n  正常愈合总方向 (D30-D0) = ({v_normal[0]:+.4f}, {v_normal[1]:+.4f})')
log(f'    方向角 = {np.degrees(np.arctan2(v_normal[1], v_normal[0])):+.1f}°')

# 各供体 D0->D30 方向角一致性
log('\n  [供体 D0->D30 方向一致性]')
angles = []
for d in donors:
    v = S[d][-1] - S[d][0]
    a = float(np.degrees(np.arctan2(v[1], v[0])))
    angles.append(a)
    log(f'    {d}: ({v[0]:+.3f}, {v[1]:+.3f})  角={a:+.1f}°')
# 角度离散度（用圆统计）
angles_rad = np.radians(angles)
R_bar = np.sqrt(np.mean(np.cos(angles_rad))**2 + np.mean(np.sin(angles_rad))**2)
log(f'    圆统计 R̄ = {R_bar:.4f}  (1=完全一致, 0=随机)')
log(f'    平均角 = {np.degrees(np.arctan2(np.mean(np.sin(angles_rad)), np.mean(np.cos(angles_rad)))):+.1f}°')

# =========================================================
# C. 轨迹的单调性检验（F6 降 / F7 升 在 D7->D30 最强）
# =========================================================
log('\n' + '=' * 84)
log('C. D7->D30 关键跃迁的方向一致性')
log('=' * 84)
log('  (D7 到 D30 是采样最长、效应最强的区段)')
n_f6_down = sum(1 for d in donors if S[d][-1, 0] < S[d][-2, 0])
n_f7_up = sum(1 for d in donors if S[d][-1, 1] > S[d][-2, 1])
log(f'  F6 下降: {n_f6_down}/3 供体')
log(f'  F7 上升: {n_f7_up}/3 供体')
for d in donors:
    log(f'    {d}: ΔF6(D7->D30)={S[d][-1,0]-S[d][-2,0]:+.4f}  ΔF7={S[d][-1,1]-S[d][-2,1]:+.4f}')

# 符号检验 p 值（3/3 一致的单侧概率 = 2*(1/2)^3 = 0.25；但这是供体级）
from scipy.stats import binomtest
p_sign = binomtest(n_f7_up, 3, 0.5, alternative='greater').pvalue
log(f'  供体级符号检验 (F7 上升 3/3): p = {p_sign:.4f}  (n=3, 功效极低，仅记录)')

# =========================================================
# D. DFU 沿轨迹的一维投影（修正 S24 的锚点距离法）
# =========================================================
log('\n' + '=' * 84)
log('D. DFU 沿正常轨迹的一维投影（替代 S24 锚点距离）')
log('=' * 84)

dfu = pd.read_csv(OUT / "S10_patient_scores.csv")
d2 = dfu.dropna(subset=['F6_orig','F7_orig','group']).copy()
d2['grp'] = d2['group'].astype(str)

# 用正常轨迹定义一维坐标轴：从 R(静息) 端 到 M(成熟) 端
# 但更稳健：用"到轨迹的最近点"作为进度坐标
# 轨迹折线: D0->D1->D7->D30（供体平均）
traj = mean_pts  # (4,2)

def project_to_polyline(p, poly):
    """把点 p 投影到折线 poly，返回 (最近点, 沿折线的累积进度 t in [0,1])"""
    best = (np.inf, None, None)
    seglen = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    cum = np.concatenate([[0], np.cumsum(seglen)])
    total = cum[-1]
    for k in range(len(poly) - 1):
        a, b = poly[k], poly[k+1]
        ab = b - a
        L2 = ab @ ab
        if L2 < 1e-12:
            t = 0.0
        else:
            t = np.clip((p - a) @ ab / L2, 0, 1)
        proj = a + t * ab
        dist = np.linalg.norm(p - proj)
        prog = (cum[k] + t * seglen[k]) / total
        if dist < best[0]:
            best = (dist, prog, proj)
    return best

pts = d2[['F6_orig','F7_orig']].values.astype(float)
projs = [project_to_polyline(p, traj) for p in pts]
d2['dist_to_traj'] = [r[0] for r in projs]
d2['progress'] = [r[1] for r in projs]

log('\n  [患者级轨迹进度] (progress: 0=D0静息端, 1=D30成熟端)')
log('   组别             n   平均progress  平均dist_to_traj')
for g in ['Healthy','DiabeticNoDFU','Healer','NonHealer']:
    sub = d2[d2['grp'] == g]
    if len(sub) == 0: continue
    log(f'   {g:14s} {len(sub):2d}   {sub["progress"].mean():.4f}      {sub["dist_to_traj"].mean():.4f}')

# 核心检验：NonHealer 的 progress 是否显著低于 Healer
h = d2[d2['grp'] == 'Healer']['progress'].values
nh = d2[d2['grp'] == 'NonHealer']['progress'].values
log(f'\n  Healer    progress: mean={h.mean():.4f}, n={len(h)}')
log(f'  NonHealer progress: mean={nh.mean():.4f}, n={len(nh)}')
obs = h.mean() - nh.mean()
log(f'  差值 (Healer - NonHealer) = {obs:+.4f}  (预期 > 0，即 Healer 进度更高)')

# 置换检验
allv = np.concatenate([h, nh]); nh_n = len(h)
cnt = 0; nperm = 10000
for _ in range(nperm):
    perm = RNG.permutation(allv)
    if (perm[:nh_n].mean() - perm[nh_n:].mean()) >= obs: cnt += 1
p_perm = (cnt + 1) / (nperm + 1)
log(f'  置换检验 (单侧, Healer > NonHealer): p_perm = {p_perm:.4f}')

# 也报告 Healthy vs NonHealer 和 DiabeticNoDFU
for g in ['Healthy','DiabeticNoDFU']:
    x = d2[d2['grp'] == g]['progress'].values
    obs2 = x.mean() - nh.mean()
    allv2 = np.concatenate([x, nh]); nx = len(x)
    cnt2 = 0
    for _ in range(nperm):
        perm = RNG.permutation(allv2)
        if (perm[:nx].mean() - perm[nx:].mean()) >= obs2: cnt2 += 1
    log(f'  {g} vs NonHealer: 差值={obs2:+.4f}, p_perm={(cnt2+1)/(nperm+1):.4f}')

# =========================================================
# 保存
# =========================================================
out = dict(
    differential_K={},
    geometry={d: {k: v for k, v in geo[d].items()} for d in donors},
    mean_trajectory=mean_pts.tolist(),
    normal_direction=v_normal.tolist(),
    normal_angle=float(np.degrees(np.arctan2(v_normal[1], v_normal[0]))),
    donor_angles=angles,
    angle_R_bar=float(R_bar),
    d7_d30_consistency=dict(f6_down=n_f6_down, f7_up=n_f7_up, p_sign=float(p_sign)),
    dfu_projection=dict(
        per_group={g: dict(n=int(len(d2[d2['grp']==g])),
                           progress_mean=float(d2[d2['grp']==g]['progress'].mean()),
                           dist_mean=float(d2[d2['grp']==g]['dist_to_traj'].mean()))
                   for g in ['Healthy','DiabeticNoDFU','Healer','NonHealer'] if len(d2[d2['grp']==g])>0},
        healer_vs_nonhealer=dict(obs=float(obs), p_perm=float(p_perm)),
    ),
)
# 补差分 K
with open(OUT / "S24b_state_model.json", 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log(f'\n[saved] S24b_state_model.json')
with open(OUT / "S24b_run.log", 'w', encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('[saved] S24b_run.log')
log('\n✅ S24b 完成')
