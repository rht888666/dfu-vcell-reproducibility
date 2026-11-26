# -*- coding: utf-8 -*-
"""
S24: F6/F7 状态转移的定量模型 —— 离散 Markov 链 + 连续动力学(ODE) + 患者级投影

目标（承接 S20S21_VERDICT_COMBINED.md §六.1）：
  把 GSE241132 的 (F6,F7) 纵向轨迹 D0->D1->D7->D30
  从"描述性轨迹"升级为"可计算的定量状态转移模型"。

四层建模（由弱到强，全部诚实标注可靠性）：
  L1 离散 Markov 转移矩阵 P[i,j] = P(state_j at t+1 | state_i at t)
      - 状态定义：(F6,F7) 四象限 -> 4 状态 {Q1..Q4}
      - 用 3 供体 × 3 转移 = 9 个观测估计（诚实：n 小）
  L2 连续动力学拟合：一维/二维线性 ODE  dX/dt = A·(X - X*)  （吸引子模型）
      - 拟合每供体 D0->D30 的弛豫，检验是否存在共同吸引子
  L3 吸引子定位：把 DFU (GSE165816) 的 Healer/NonHealer 映射到同一相空间
      - 检验 NonHealer 是否更接近 D0 而非 D30 吸引子
  L4 模型选择：比较 (a) 共同吸引子模型 vs (b) 供体特异模型 vs (c) 无漂移模型

数据：
  - S20_donor_trajectory.csv  (3 供体 × 4 时间点 × 8 程序得分)
  - GSE165816 DFU 患者级 (F6,F7) —— 从 S10_patient_scores.csv / S14 状态空间读取

统计：留一供体交叉验证(LOO)、置换检验、AIC/BIC
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
log('S24: F6/F7 状态转移的定量模型（Markov + ODE 吸引子 + DFU 投影）')
log('=' * 84)

# =========================================================
# 0. 载入数据
# =========================================================
dt = pd.read_csv(OUT / "S20_donor_trajectory.csv")
dt['tp'] = pd.Categorical(dt['tp'], categories=['D0','D1','D7','D30'], ordered=True)
dt = dt.sort_values(['donor','tp']).reset_index(drop=True)
log(f'\n[data] donor_trajectory: {dt.shape[0]} 行, 供体 {sorted(dt.donor.unique())}')

# 时间轴（天）。注意 D1 与 D7 之间无采样 —— 用真实天数
T_DAYS = {'D0': 0.0, 'D1': 1.0, 'D7': 7.0, 'D30': 30.0}

# 使用 F6（炎症肌成纤维）与 F7（肌成纤维）构成二维状态空间
PROG_X, PROG_Y = 'F6', 'F7'
donors = sorted(dt.donor.unique())

def dmat(donor):
    """返回 (4,2) 矩阵 [F6,F7] 按 D0,D1,D7,D30 排序"""
    sub = dt[dt.donor == donor].set_index('tp').loc[['D0','D1','D7','D30']]
    return sub[[PROG_X, PROG_Y]].values

S = {d: dmat(d) for d in donors}
for d in donors:
    log(f'  {d}: ' + ' | '.join(f'{t}({S[d][i,0]:+.3f},{S[d][i,1]:+.3f})'
                                for i, t in enumerate(['D0','D1','D7','D30'])))

# =========================================================
# L1. 离散 Markov 链（四象限状态）
# =========================================================
log('\n' + '=' * 84)
log('L1. 离散 Markov 状态转移矩阵（(F6,F7) 四象限）')
log('=' * 84)

# 以坐标原点为分界（F6=0, F7=0），定义四象限
#   0: F6- F7-   (静息/低活性)      -> 记 R (Rest)
#   1: F6+ F7-   (炎症期, F6 主导)  -> 记 I (Inflammatory)
#   2: F6- F7+   (重塑期, F7 主导)  -> 记 M (Maturation)
#   3: F6+ F7+   (双高, 过渡/混合)  -> 记 T (Transition)
QNAME = ['R(F6-F7-)', 'I(F6+F7-)', 'M(F6-F7+)', 'T(F6+F7+)']

def quadrant(x, y):
    if x >= 0 and y < 0: return 1   # I
    if x < 0 and y >= 0: return 2   # M
    if x >= 0 and y >= 0: return 3  # T
    return 0                        # R

# 统计每条转移边（供体级），构建转移计数矩阵 & 状态占用
trans_counts = np.zeros((4, 4))
state_seq = {}
for d in donors:
    qs = [quadrant(x, y) for x, y in S[d]]
    state_seq[d] = qs
    for a, b in zip(qs[:-1], qs[1:]):
        trans_counts[a, b] += 1
    log(f'  {d}: {" -> ".join(QNAME[q] for q in qs)}')

# 行归一化为转移概率
P = np.zeros((4, 4))
for i in range(4):
    rs = trans_counts[i].sum()
    if rs > 0:
        P[i] = trans_counts[i] / rs

log('\n  转移计数矩阵 (行=from, 列=to):')
log('        ' + '  '.join(f'{n:>12s}' for n in QNAME))
for i in range(4):
    log(f'  {QNAME[i]:>12s} ' + '  '.join(f'{int(trans_counts[i,j]):12d}' for j in range(4)))
log('\n  转移概率矩阵 P:')
log('        ' + '  '.join(f'{n:>12s}' for n in QNAME))
for i in range(4):
    log(f'  {QNAME[i]:>12s} ' + '  '.join(f'{P[i,j]:12.3f}' for j in range(4)))

# 吸收态分析：从 D0 出发，30 天后的状态分布（单步转移 × 3 次：D0->D1->D7->D30）
# 注意：转移时间间隔不等（1, 6, 23 天）—— 这是离散 Markov 的固有局限，诚实标注
pi0 = np.zeros(4); pi0[state_seq[donors[0]][0]] += 0
# 用供体平均的初始分布
q0 = [state_seq[d][0] for d in donors]
pi0 = np.bincount(q0, minlength=4) / len(donors)
pi_t = pi0.copy()
traj_pi = [pi_t.copy()]
for _ in range(3):
    pi_t = pi_t @ P
    traj_pi.append(pi_t.copy())

log('\n  [吸收态/终态分布] 从 D0 初始分布出发的多步演化：')
for k, t in enumerate(['D0','D1','D7','D30']):
    log(f'    {t}: ' + '  '.join(f'{QNAME[j]}={traj_pi[k][j]:.3f}' for j in range(4)))

log('\n  ⚠️ 局限：转移间隔不等（1,6,23 天），离散 Markov 假设时间齐次性 —— 与真实不符。')
log('     此层仅作"是否存在确定性方向"的定性检验，不作定量外推。')

# 确定性检验：转移是否比随机更确定（对比均匀转移）
# 用熵衡量：每个状态的出边分布熵
def entropy(p):
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())
ent = [entropy(P[i]) if trans_counts[i].sum() > 0 else np.nan for i in range(4)]
log(f'\n  各状态出边熵: ' + ', '.join(f'{QNAME[i]}={ent[i]:.3f}' for i in range(4)))
log(f'  (最大熵 log4 = {np.log(4):.3f}；越接近 0 越确定)')

# =========================================================
# L2. 连续动力学：线性吸引子 ODE  dX/dt = -K (X - X*)
# =========================================================
log('\n' + '=' * 84)
log('L2. 连续动力学：线性吸引子模型  dX/dt = -K·(X - X*)')
log('=' * 84)
log('  物理含义：成纤维状态以速率 K 弛豫到吸引子 X*（成熟重塑态）。')

# 对每个供体，用最小二乘拟合 2D 线性 ODE（每维独立拟合 K, X*）：
#   x(t_{k+1}) - x(t_k) = -K (x(t_k) - X*) * dt
#   => dx = -K*dt * x + K*dt*X*
#   => 线性回归 dx ~ x  给出斜率 a = -K*dt_varies... 
# 由于 dt 不等，采用积分形式更稳：拟合 x(t) = X* + (x0 - X*) exp(-K t)
# 但 4 点、2 参数，用网格搜索 K 最稳。

def fit_attractor_1d(ts, xs):
    """拟合 x(t) = Xstar + (x0 - Xstar) exp(-K t)，网格搜索 K，线性解 Xstar"""
    ts = np.asarray(ts, float); xs = np.asarray(xs, float)
    x0 = xs[0]
    best = None
    for K in np.linspace(0.001, 3.0, 3000):
        phi = np.exp(-K * ts)
        # xs = Xstar + (x0 - Xstar) phi = Xstar(1-phi) + x0 phi
        a = (1 - phi)
        # 最小二乘解 Xstar： xs - x0*phi = Xstar * a
        denom = (a * a).sum()
        if denom < 1e-12: continue
        Xstar = float(((xs - x0 * phi) * a).sum() / denom)
        pred = Xstar + (x0 - Xstar) * phi
        sse = float(((xs - pred) ** 2).sum())
        if best is None or sse < best[0]:
            best = (sse, K, Xstar, pred)
    return best  # (sse, K, Xstar, pred)

ts = [T_DAYS[t] for t in ['D0','D1','D7','D30']]

fit_results = {}
log('\n  [每供体独立拟合]')
log('   供体     程序    K(速率)   X*(吸引子)   SSE     R2')
for d in donors:
    fit_results[d] = {}
    for j, prog in enumerate([PROG_X, PROG_Y]):
        xs = S[d][:, j]
        sse, K, Xs, pred = fit_attractor_1d(ts, xs)
        ss_tot = float(((xs - xs.mean()) ** 2).sum())
        r2 = 1 - sse / ss_tot if ss_tot > 0 else np.nan
        fit_results[d][prog] = dict(K=K, Xstar=Xs, sse=sse, r2=r2)
        log(f'   {d:6s}  {prog:4s}   {K:6.3f}   {Xs:+8.4f}   {sse:.5f}  {r2:.3f}')

# 共同吸引子 vs 供体特异
log('\n  [共同吸引子检验]（3 供体拟合的 X* 是否一致）')
common = {}
for prog in [PROG_X, PROG_Y]:
    Xs_list = np.array([fit_results[d][prog]['Xstar'] for d in donors])
    K_list = np.array([fit_results[d][prog]['K'] for d in donors])
    common[prog] = dict(Xstar_mean=float(Xs_list.mean()),
                        Xstar_sd=float(Xs_list.std(ddof=1)),
                        Xstar_cv=float(Xs_list.std(ddof=1) / abs(Xs_list.mean())) if Xs_list.mean() != 0 else np.nan,
                        K_mean=float(K_list.mean()), K_sd=float(K_list.std(ddof=1)),
                        Xstar_all=Xs_list.tolist(), K_all=K_list.tolist())
    log(f'   {prog}: X* 均值={Xs_list.mean():+.4f}  SD={Xs_list.std(ddof=1):.4f}  CV={common[prog]["Xstar_cv"]:.3f}')
    log(f'        K  均值={K_list.mean():.3f}   SD={K_list.std(ddof=1):.3f}')

# 共同吸引子模型（固定 X*=均值，每供体拟合 K）
log('\n  [共同吸引子模型] 固定 X* 为供体均值，重新拟合 K：')
for prog in [PROG_X, PROG_Y]:
    Xs_fix = common[prog]['Xstar_mean']
    sses = []
    for d in donors:
        xs = S[d][:, :][:, [PROG_X, PROG_Y].index(prog)]
        x0 = xs[0]
        best = None
        for K in np.linspace(0.001, 3.0, 3000):
            pred = Xs_fix + (x0 - Xs_fix) * np.exp(-K * np.asarray(ts, float))
            sse = float(((xs - pred) ** 2).sum())
            if best is None or sse < best[0]:
                best = (sse, K)
        sses.append(best[0])
    log(f'   {prog}: 共同 X*={Xs_fix:+.4f} -> 总 SSE={sum(sses):.5f}')

# =========================================================
# L3. 模型选择（AIC/BIC）：共同吸引子 vs 供体特异 vs 无漂移
# =========================================================
log('\n' + '=' * 84)
log('L3. 模型选择：共同吸引子 vs 供体特异 vs 无漂移')
log('=' * 84)

def fit_no_drift(ts, xs):
    """模型C：x(t)=常数（无漂移），参数 1 个"""
    c = xs.mean()
    sse = float(((xs - c) ** 2).sum())
    return sse, 1

models = {}
for prog in [PROG_X, PROG_Y]:
    N_obs = len(donors) * 4           # 12 观测
    # 模型A：共同吸引子（X* 全局 1 参数 + 每供体 K 3 参数 = 4 参数）
    Xs_fix = common[prog]['Xstar_mean']
    sseA = 0.0
    for d in donors:
        xs = S[d][:, :][:, [PROG_X, PROG_Y].index(prog)]
        x0 = xs[0]; best = None
        for K in np.linspace(0.001, 3.0, 3000):
            pred = Xs_fix + (x0 - Xs_fix) * np.exp(-K * np.asarray(ts, float))
            sse = float(((xs - pred) ** 2).sum())
            if best is None or sse < best[0]: best = (sse,)
        sseA += best[0]
    kA = 1 + len(donors)  # X* + 每供体 K
    # 模型B：供体特异（每供体 X*, K = 6 参数）
    sseB = sum(fit_results[d][prog]['sse'] for d in donors)
    kB = 2 * len(donors)
    # 模型C：无漂移（每供体常数 = 3 参数）
    sseC = 0.0
    for d in donors:
        xs = S[d][:, :][:, [PROG_X, PROG_Y].index(prog)]
        s, k = fit_no_drift(ts, xs)
        sseC += s
    kC = len(donors)

    def aic(sse, k, N):  return N * np.log(sse / N) + 2 * k
    def bic(sse, k, N):  return N * np.log(sse / N) + k * np.log(N)

    models[prog] = dict(
        A_common=dict(sse=sseA, k=kA, AIC=aic(sseA, kA, N_obs), BIC=bic(sseA, kA, N_obs)),
        B_donor=dict(sse=sseB, k=kB, AIC=aic(sseB, kB, N_obs), BIC=bic(sseB, kB, N_obs)),
        C_nodrift=dict(sse=sseC, k=kC, AIC=aic(sseC, kC, N_obs), BIC=bic(sseC, kC, N_obs)),
    )
    log(f'\n  {prog}:')
    for name, m in models[prog].items():
        log(f'    {name:12s} k={m["k"]}  SSE={m["sse"]:.5f}  AIC={m["AIC"]:.2f}  BIC={m["BIC"]:.2f}')
    # 最佳模型
    best = min(models[prog].items(), key=lambda kv: kv[1]['AIC'])
    log(f'    -> AIC 最佳: {best[0]}')

# =========================================================
# L4. DFU 投影：NonHealer 是否更接近 D0 而非 D30 吸引子
# =========================================================
log('\n' + '=' * 84)
log('L4. DFU 患者投影 —— NonHealer 停滞在 D0 附近的定量检验')
log('=' * 84)

# 读取 DFU 患者级 (F6,F7)
dfu_path = OUT / "S10_patient_scores.csv"
dfu = None
if dfu_path.exists():
    dfu = pd.read_csv(dfu_path)
    log(f'\n  载入 S10_patient_scores.csv: {dfu.shape}  列={list(dfu.columns)[:12]}')
else:
    log(f'  ⚠️ 未找到 {dfu_path}')

if dfu is not None:
    # 自动识别 F6/F7 列与分组列（优先 *_orig，其次 *_ecm/*_il11 等主程序）
    def pick(pref):
        cand = [c for c in dfu.columns if c.upper().startswith(pref)]
        if not cand: return []
        for c in cand:
            if c.upper() == pref + '_ORIG': return [c]
        return cand[:1]
    f6col = pick('F6')
    f7col = pick('F7')
    grpcol = [c for c in dfu.columns if c.lower() in ('group','condition','label','outcome','heal')][:1]
    log(f'  识别: F6列={f6col}, F7列={f7col}, 分组列={grpcol}')

    if f6col and f7col and grpcol:
        d2 = dfu.dropna(subset=[f6col[0], f7col[0], grpcol[0]]).copy()
        d2['grp'] = d2[grpcol[0]].astype(str)

        # 正常愈合的吸引子（D30）与起点（D0）—— 供体平均
        D0_pt = np.array([np.mean([S[d][0, 0] for d in donors]), np.mean([S[d][0, 1] for d in donors])])
        D30_pt = np.array([np.mean([S[d][3, 0] for d in donors]), np.mean([S[d][3, 1] for d in donors])])
        D1_pt = np.array([np.mean([S[d][1, 0] for d in donors]), np.mean([S[d][1, 1] for d in donors])])
        log(f'\n  正常轨迹锚点: D0={D0_pt.round(3)}, D1={D1_pt.round(3)}, D30={D30_pt.round(3)}')

        # 每个 DFU 患者到 D0 / D30 锚点的距离
        pts = d2[[f6col[0], f7col[0]]].values.astype(float)
        dist_D0 = np.linalg.norm(pts - D0_pt, axis=1)
        dist_D30 = np.linalg.norm(pts - D30_pt, axis=1)
        d2['dist_D0'] = dist_D0
        d2['dist_D30'] = dist_D30
        d2['closer_to'] = np.where(dist_D0 < dist_D30, 'D0', 'D30')

        log('\n  [患者级] 到 D0/D30 锚点距离：')
        for g in sorted(d2['grp'].unique()):
            sub = d2[d2['grp'] == g]
            n = len(sub)
            fracD0 = (sub['closer_to'] == 'D0').mean()
            log(f'    {g:14s} n={n:2d}  平均dist_D0={sub["dist_D0"].mean():.3f}  '
                f'平均dist_D30={sub["dist_D30"].mean():.3f}  更近D0比例={fracD0:.2%}')

        # 组间检验：Healer vs NonHealer 的 dist_D0 - dist_D30
        groups = sorted(d2['grp'].unique())
        # 找 healer / nonhealer 组名
        heal_g = [g for g in groups if 'heal' in g.lower() and 'non' not in g.lower()]
        nonh_g = [g for g in groups if 'non' in g.lower() and 'heal' in g.lower()]
        log(f'\n  识别 healer 组={heal_g}, nonhealer 组={nonh_g}')

        if heal_g and nonh_g:
            from scipy import stats as st
            dh = (d2[d2['grp'].isin(heal_g)]['dist_D0'] - d2[d2['grp'].isin(heal_g)]['dist_D30']).values
            dn = (d2[d2['grp'].isin(nonh_g)]['dist_D0'] - d2[d2['grp'].isin(nonh_g)]['dist_D30']).values
            # 该差值 > 0 表示更接近 D30；< 0 表示更接近 D0
            log(f'    Healer    (dist_D0 - dist_D30): mean={dh.mean():+.4f} (n={len(dh)})')
            log(f'    NonHealer (dist_D0 - dist_D30): mean={dn.mean():+.4f} (n={len(dn)})')
            if len(dh) >= 2 and len(dn) >= 2:
                u, pu = st.mannwhitneyu(dh, dn, alternative='greater')
                log(f'    Mann-Whitney U (Healer > NonHealer, 单侧): U={u:.1f}, p={pu:.4f}')
            # 置换检验
            allv = np.concatenate([dh, dn]); nh = len(dh)
            obs = dh.mean() - dn.mean()
            cnt = 0; nperm = 5000
            for _ in range(nperm):
                perm = RNG.permutation(allv)
                if (perm[:nh].mean() - perm[nh:].mean()) >= obs: cnt += 1
            log(f'    置换检验 (Healer-NonHealer 差值, 单侧): obs={obs:+.4f}, p_perm={(cnt+1)/(nperm+1):.4f}')

# =========================================================
# 保存
# =========================================================
out = dict(
    programs=[PROG_X, PROG_Y],
    quadrants=QNAME,
    donor_states={d: [QNAME[q] for q in state_seq[d]] for d in donors},
    trans_counts=trans_counts.tolist(),
    trans_prob=P.tolist(),
    state_entropy={QNAME[i]: (None if np.isnan(ent[i]) else ent[i]) for i in range(4)},
    multi_step_distribution={t: traj_pi[k].tolist() for k, t in enumerate(['D0','D1','D7','D30'])},
    attractor_fit={d: fit_results[d] for d in donors},
    common_attractor=common,
    model_selection=models,
)
with open(OUT / "S24_markov_state_model.json", 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log(f'\n[saved] S24_markov_state_model.json')

with open(OUT / "S24_run.log", 'w', encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('[saved] S24_run.log')
log('\n✅ S24 完成')
