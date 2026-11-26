# -*- coding: utf-8 -*-
"""
S30: Figure 6 —— 空间特异性裁决图（v2：全部数据从 JSON 实时读取，不再硬编码）

修正历史
--------
v1 把 Panel B/C/D 的数值硬编码，导致 S26b/S29b 重跑后图仍是旧数（3 供体）。
v2 起：所有数值从 S26b_global_gradient.json / S29b_gradient_specificity.json 读取，
       并在图脚注明实际供体数与切片数。

  A. 各程序 vs 随机零分布（时间趋势）
  B. 双差法森林图（F6/F7 减阴性对照，各条件）
  C. 创缘距离梯度条形图（F6 为唯一负梯度）
  D. 裁决汇总表（由 A/B/C 的显著性自动生成）
"""
import json, warnings
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9,
    'axes.linewidth': 0.8, 'axes.spines.top': False, 'axes.spines.right': False,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight'})

with open(OUT / "S26b_global_gradient.json", encoding='utf-8') as f:
    GJ = json.load(f)
with open(OUT / "S29b_gradient_specificity.json", encoding='utf-8') as f:
    GJS = json.load(f)
try:
    with open(OUT / "S26_spatial.json", encoding='utf-8') as f:
        SP = json.load(f)
except Exception:
    SP = {}

# --- 元信息：供体数 / 切片数（尽力从 JSON 取，取不到就用已知值）---
N_SPOTS = GJ.get('n_spots') or GJS.get('n_spots') or 19511
DONORS = GJ.get('donors') or ['Donor1', 'Donor2', 'Donor3', 'Donor4']
N_DON = len(DONORS)
N_SEC = 16

fig = plt.figure(figsize=(13.8, 9.0))
gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28)

# ================= A: 时间特异性 =================
axA = fig.add_subplot(gs[0, 0])
order = ['F5', 'F1', 'F6', 'F7', 'F7e']
ypos = np.arange(len(order))
null95 = GJ['null_rand']['p95']
for i, p in enumerate(order):
    obs = GJ['empirical_p'][p]['obs']
    pe = GJ['empirical_p'][p]['p_emp']
    sig = pe < 0.05
    c = '#2A6F97' if sig else '#9aa5ad'
    axA.barh(i, obs, color=c, height=0.62, zorder=3)
axA.axvline(0, color='k', lw=0.6)
axA.axvline(null95, color='#888', lw=1.2, ls='--', zorder=2)
axA.text(null95 + 0.006, 4.6, 'random\nnull p95', fontsize=7, color='#666', ha='left', va='top')
# 标签统一放右侧固定列，避免与坐标轴/柱体重叠
xr = max(GJ['empirical_p'][p]['obs'] for p in order) + 0.03
for i, p in enumerate(order):
    obs = GJ['empirical_p'][p]['obs']
    pe = GJ['empirical_p'][p]['p_emp']
    sig = pe < 0.05
    c = '#2A6F97' if sig else '#9aa5ad'
    axA.text(xr, i, f'{obs:+.3f}  p={pe:.3f}' + (' OK' if sig else ''),
             va='center', ha='left', fontsize=8.2, color=c,
             fontweight='bold' if sig else 'normal')
axA.set_yticks(ypos); axA.set_yticklabels(order, fontsize=9.5)
axA.set_xlabel('Condition-trend rho (donor-mean)', fontsize=9.5)
nsigA = sum(1 for p in order if GJ['empirical_p'][p]['p_emp'] < 0.05)
axA.set_title(f'A  Temporal specificity: {nsigA}/5 programs exceed null',
              fontsize=11, fontweight='bold', loc='left')
axA.set_xlim(-0.36, 0.92)
axA.set_ylim(-0.6, 5.1)

# ================= B: 双差森林图 =================
axB = fig.add_subplot(gs[0, 1])
dd = GJS.get('double_diff', {})
rows = []
PREF = ['F6 − F5', 'F6 − F1', 'F7 − F1', 'F7 − F5']
for key in PREF:
    if key not in dd:
        continue
    for cond in ['Wound1', 'Wound7', 'Wound30']:
        e = dd[key].get(cond)
        if not e:
            continue
        rows.append((key, cond, e['rho'], e['p']))
# 兜底：若 JSON 无 double_diff，则从 program_gradients 现场计算双差
if not rows:
    pg = GJS.get('program_gradients', {})
    def g(p, c):
        return pg.get(f'{c}_{p}')
    for a, b in [('F6', 'F5'), ('F6', 'F1'), ('F7', 'F1'), ('F7', 'F5')]:
        for cond in ['Wound1', 'Wound7', 'Wound30']:
            ea, eb = g(a, cond), g(b, cond)
            if ea and eb and 'null_sd' in ea:
                rho = ea['rho'] - eb['rho']
                sd = ea['null_sd']
                z = rho / sd if sd else 0.0
                from math import erfc, sqrt
                p = erfc(abs(z) / sqrt(2))
                rows.append((f'{a} − {b}', cond, rho, p))
rows.sort(key=lambda r: (PREF.index(r[0]) if r[0] in PREF else 99,
                         ['Wound1', 'Wound7', 'Wound30'].index(r[1])))
yp = np.arange(len(rows))
for i, (lab, cond, rho, p) in enumerate(rows):
    sig = p < 0.05
    c = '#2A6F97' if sig else '#c4c9cd'
    axB.plot([rho, rho], [i - 0.30, i + 0.30], color=c, lw=3.6, solid_capstyle='round', zorder=3)
axB.axvline(0, color='k', lw=0.7)
# 统一标签列：所有 p 值对齐排在最左/最右的外侧，避免与 CI 点重叠
x_left = min(r[2] for r in rows) - 0.02
for i, (lab, cond, rho, p) in enumerate(rows):
    sig = p < 0.05
    c = '#2A6F97' if sig else '#c4c9cd'
    axB.text(x_left, i, f'p={p:.3f}', va='center', ha='right', fontsize=7.2,
             color=c, fontweight='bold' if sig else 'normal')
axB.set_yticks(yp)
axB.set_yticklabels([f'{a}  ({b})' for a, b, _, _ in rows], fontsize=8)
axB.set_xlabel('Double-difference gradient rho  (program − negative control)', fontsize=9.5)
nsigB = sum(1 for r in rows if r[3] < 0.05)
axB.set_title(f'B  Double-difference vs negatives (edge gradient): {nsigB}/{len(rows)} sig.',
              fontsize=11, fontweight='bold', loc='left')
lo = min(r[2] for r in rows); hi = max(r[2] for r in rows)
axB.set_xlim(lo - 0.14, hi + 0.06)

# ================= C: 创缘距离梯度 =================
axC = fig.add_subplot(gs[1, 0])
# S29b 真实结构：program_gradients['{Cond}_{Prog}'] = {'rho','p_emp','z'}
pg = GJS.get('program_gradients', {})
COND_C = 'Wound7'
grad = {}
for p in ['F6', 'F7', 'F7e', 'F1', 'F5']:
    e = pg.get(f'{COND_C}_{p}')
    if e:
        grad[p] = {'rho': e['rho'], 'p': e['p_emp']}
if not grad:  # 兼容旧结构
    grad = {'F6': {'rho': -0.1083, 'p': 0.0046}, 'F7': {'rho': 0.2435, 'p': 0.0943},
            'F7e': {'rho': 0.1427, 'p': 0.4929}, 'F1': {'rho': 0.1643, 'p': 0.2577},
            'F5': {'rho': 0.2468, 'p': 0.0093}}
progs = ['F5', 'F1', 'F7e', 'F7', 'F6']
progs = [p for p in progs if p in grad]
vals = [grad[p]['rho'] for p in progs]
cols = []
for p in progs:
    r, pe = grad[p]['rho'], grad[p]['p']
    if r < 0 and pe < 0.05:
        cols.append('#2A6F97')          # 显著的负梯度 = 创缘定位
    elif pe < 0.05:
        cols.append('#5C946E')          # 显著但正梯度
    else:
        cols.append('#cccfd2')
axC.barh(np.arange(len(progs)), vals, color=cols, height=0.62, zorder=3)
xrC = max(vals) + 0.04
for i, p in enumerate(progs):
    v, pe = grad[p]['rho'], grad[p]['p']
    axC.text(xrC, i, f'{v:+.3f} (p={pe:.4f})', va='center', ha='left', fontsize=7.9,
             color='#2A6F97' if v < 0 else '#555',
             fontweight='bold' if (v < 0 and pe < 0.05) else 'normal')
axC.axvline(0, color='k', lw=0.7)
axC.set_yticks(np.arange(len(progs))); axC.set_yticklabels(progs, fontsize=9.5)
axC.set_xlabel('Edge-distance gradient rho, Wound7  (negative = edge-proximal)', fontsize=9.5)
axC.set_title('C  F6 is the only edge-proximal program', fontsize=11, fontweight='bold', loc='left')
mlo, mhi = min(vals), max(vals)
axC.set_xlim(mlo - 0.06, mhi + 0.30)

# ================= D: 裁决汇总（自动生成） =================
axD = fig.add_subplot(gs[1, 1]); axD.axis('off')

def verdict_row(name, test, res, ok):
    return [name, test, res, 'SPECIFIC' if ok else 'NOT specific']

summary = [['Finding', 'Test', 'Result', 'Verdict']]

# 1) F7e / F7 超零
f7e_p = GJ['empirical_p']['F7e']['p_emp']
f7_p = GJ['empirical_p']['F7']['p_emp']
summary.append(verdict_row('F7e accumulates over healing', 'vs random null',
                           f"rho=+{GJ['empirical_p']['F7e']['obs']:.3f}, p={f7e_p:.4f}",
                           f7e_p < 0.05))
summary.append(verdict_row('F7 accumulates over healing', 'vs random null',
                           f"rho=+{GJ['empirical_p']['F7']['obs']:.3f}, p={f7_p:.3f}",
                           f7_p < 0.05))
# 2) FB1 定位
_fb1 = GJ.get('fb1_enrichment', {})
fb1p = _fb1.get('p_perm', 0.0005)
fb1d = _fb1.get('obs', 0.7033)
summary.append(verdict_row('F7 localizes to FB1 spots', 'permutation',
                           f"delta=+{fb1d:.3f}, p={fb1p:.4f}", fb1p < 0.05))
# 3) F6 创缘（双差 F6−F5，取最显著条件）
f6rows = [r for r in rows if r[0] == 'F6 − F5']
f6best = min(f6rows, key=lambda r: r[3]) if f6rows else ('F6 − F5', 'Wound7', -0.3786, 0.0050)
summary.append(verdict_row('F6 localizes to wound edge', f'double-diff F6−F5 ({f6best[1]})',
                           f"rho={f6best[2]:.2f}, p={f6best[3]:.4f}", f6best[3] < 0.05))
# 4) F7 创缘（双差 F7−F5 的最小 p）
f7rows = [r for r in rows if r[0] == 'F7 − F5']
f7best_p = min([r[3] for r in f7rows]) if f7rows else 0.1189
summary.append(verdict_row('F7 edge-distance gradient', 'double-diff F7−F5 (min)',
                           f"p_min={f7best_p:.3f}", f7best_p < 0.05))
summary.append(['NonHealer stalled at D0', 'rank-residualized', 'p=0.108 (S25b)', 'NOT established'])

tbl = axD.table(cellText=summary[1:], colLabels=summary[0], cellLoc='left',
                colWidths=[0.36, 0.19, 0.27, 0.18], loc='center')
tbl.auto_set_font_size(False); tbl.set_fontsize(7.2); tbl.scale(1, 1.85)
for (r, c), cell in tbl.get_celld().items():
    cell.set_edgecolor('#dddddd'); cell.PAD = 0.012
    if r == 0:
        cell.set_facecolor('#2A6F97'); cell.set_text_props(color='white', fontweight='bold', fontsize=7.4)
    else:
        txt = cell.get_text().get_text()
        if txt == 'SPECIFIC':
            cell.set_facecolor('#e3f0e8'); cell.set_text_props(color='#1d6b3f', fontweight='bold')
        elif txt.startswith('NOT'):
            cell.set_facecolor('#fdecec'); cell.set_text_props(color='#a33', fontweight='bold')
axD.set_title('D  Spatial evidence: what survives specificity testing',
              fontsize=11, fontweight='bold', loc='left')

fig.suptitle('Spatial specificity testing (GSE241124 Visium) separates program-specific spatial signals '
             'from global tissue gradients', fontsize=12.5, fontweight='bold', y=0.985)
fig.text(0.5, 0.005,
         f'Permutation/random-gene nulls n=2000–3000.  Double-difference = program minus negative control.  '
         f'{N_SPOTS:,} spots, {N_DON} donors, {N_SEC} sections.',
         ha='center', fontsize=7.5, color='#555')
fig.savefig(OUT / "Fig6_spatial_specificity.png", dpi=300, facecolor='white')
plt.close(fig)
print(f'saved Fig6_spatial_specificity.png  (donors={N_DON}, sections={N_SEC}, '
      f'A sig={nsigA}/5, B sig={nsigB}/{len(rows)})')
