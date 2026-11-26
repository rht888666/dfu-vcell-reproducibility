# -*- coding: utf-8 -*-
"""
S22: 四队列跨队列证据图（Nature 级）
  Panel A: 细胞级效应量森林图（4 队列 × F6/F7e/F1 对照）
  Panel B: 深度差 vs 可检出性（散点，展示深度控制方法论）
  Panel C: 三队列阴性对照 vs 信号对照（柱状）
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

warnings.filterwarnings('ignore')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8
plt.rcParams['font.size'] = 9

ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"

# ============ 手工汇总（全部来自机器生成的 json/log，此处仅为绘图汇总）============
# 来源: S12_crosscohort_GSE231643.json, S21_crosscohort.json, FINDING_DFU_F6F7.md (S11)
data = [
    # cohort, label, depth_delta_pct, F6_delta, F6_p, F7e_delta, F7e_p, F1_delta, F1_p, verdict
    dict(cohort='GSE165816', label='GSE165816\n(31 pts, Healer/NonHealer)', depth=28.0,
         F6=0.1901, F6p=0.0038, F7e=0.2365, F7ep=0.0000, F1=-0.0863, F1p=0.9504, verdict='core'),
    dict(cohort='GSE231643', label='GSE231643\n(8 pts, Healer/NonHealer)', depth=2.0,
         F6=0.2018, F6p=0.0005, F7e=0.1615, F7ep=0.0005, F1=None, F1p=None, verdict='replicate'),
    dict(cohort='GSE223964', label='GSE223964\n(8 pts, Diabetic/Non-diab)', depth=-0.4,
         F6=0.1704, F6p=0.0002, F7e=0.4727, F7ep=0.0002, F1=-0.3135, F1p=1.0000, verdict='replicate'),
    dict(cohort='GSE245703', label='GSE245703\n(12 pts, Diabetic/Non-diab)', depth=20.0,
         F6=0.0133, F6p=0.0828, F7e=0.0646, F7ep=0.0002, F1=0.0407, F1p=0.0002, verdict='rejected'),
]

fig = plt.figure(figsize=(13.5, 9.0), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=0.48, wspace=0.30)

VCOL = {'core': '#1a4f8a', 'replicate': '#2e8b57', 'rejected': '#b03a2e'}

# ---------------- Panel A: 森林图 ----------------
axA = fig.add_subplot(gs[0, :])
progs = [('F6', 'F6', 'F6p'), ('F7e', 'F7e', 'F7ep'), ('F1', 'F1 (neg ctrl)', 'F1p')]
ypos = np.arange(len(data))[::-1]
off = 0.22
mk = {'F6': 'o', 'F7e': 's', 'F1': '^'}
col = {'F6': '#c0392b', 'F7e': '#2471a3', 'F1': '#7d3c98'}
for i, (pk, pl, ppk) in enumerate(progs):
    xs, ys, cs = [], [], []
    for j, d in enumerate(data):
        if d[pk] is None:
            continue
        xs.append(d[pk]); ys.append(ypos[j] + (i - 1) * off)
        cs.append(VCOL[d['verdict']] if pk != 'F1' else '#7d3c98')
    axA.scatter(xs, ys, s=130, marker=mk[pk], c=cs, edgecolor='black',
                linewidth=0.6, label=pl, zorder=3)
axA.axvline(0, color='gray', lw=1.2, ls='--', zorder=1)
axA.set_yticks(ypos)
axA.set_yticklabels([d['label'] for d in data], fontsize=8.5)
axA.set_xlabel('Effect size Δ (cell-level, depth-residualized score)', fontsize=10)
axA.set_title('A  Cross-cohort effect sizes of F6 / F7e repair programs and F1 negative control',
              fontsize=11, fontweight='bold', loc='left')
axA.legend(loc='lower right', frameon=False, fontsize=9)
# 标注 p 值
for j, d in enumerate(data):
    txt = f"F6 p={d['F6p']:.4f}   F7e p={d['F7ep']:.4f}"
    if d['F1p'] is not None:
        txt += f"   F1 p={d['F1p']:.4f}"
    axA.text(0.52, ypos[j] + 0.30, txt, fontsize=7.2,
             color=VCOL[d['verdict']], va='center')
# 右侧 verdict 标签
for j, d in enumerate(data):
    lab = {'core': 'CORE', 'replicate': 'REPLICATED',
           'rejected': 'REJECTED (depth artifact)'}[d['verdict']]
    axA.text(0.60, ypos[j], lab, fontsize=8, fontweight='bold',
             color=VCOL[d['verdict']], va='center')

# ---------------- Panel B: 深度差 vs 可检出性 ----------------
axB = fig.add_subplot(gs[1, 0])
# 标签偏移：GSE165816 右移，GSE223964 上移
offsets = {'GSE165816': (0, -18), 'GSE231643': (0, 13),
           'GSE223964': (28, 6), 'GSE245703': (0, 13)}
for d in data:
    ok = d['F7ep'] <= 0.01 and (d['F1p'] is None or d['F1p'] > 0.05)
    m = 'o' if ok else 'X'
    axB.scatter(d['depth'], abs(d['F7e']), s=220, marker=m,
                c=VCOL[d['verdict']], edgecolor='black', linewidth=0.8, zorder=3)
    axB.annotate(d['cohort'], (d['depth'], abs(d['F7e'])),
                 textcoords='offset points', xytext=offsets[d['cohort']],
                 ha='center', fontsize=8)
axB.axhline(0, color='gray', lw=0.7)
axB.set_xlim(-6, 36)
axB.set_ylim(-0.04, 0.56)
axB.set_xlabel('Depth difference between groups (% log-UMI)', fontsize=10)
axB.set_ylabel('|F7e effect size|', fontsize=10)
axB.set_title('B  Detectability is independent of depth difference', fontsize=11,
              fontweight='bold', loc='left')
axB.text(0.97, 0.96, '● passes negative control\n✕ fails (F1 also significant)',
         transform=axB.transAxes, ha='right', va='top', fontsize=8,
         bbox=dict(boxstyle='round', fc='#f4f4f4', ec='gray', lw=0.6))

# ---------------- Panel C: 阴性对照 vs 信号（三信息队列）----------------
axC = fig.add_subplot(gs[1, 1])
info = [d for d in data if d['verdict'] in ('core', 'replicate') and d['F1'] is not None]
x = np.arange(len(info))
w = 0.35
axC.bar(x - w/2, [abs(d['F6']) for d in info], w, label='|F6|',
        color='#c0392b', edgecolor='black', linewidth=0.6)
axC.bar(x + w/2, [abs(d['F1']) for d in info], w, label='|F1| (neg ctrl)',
        color='#7d3c98', edgecolor='black', linewidth=0.6)
axC.set_xticks(x)
axC.set_xticklabels([d['cohort'] for d in info], fontsize=9)
axC.set_ylabel('|Effect size|', fontsize=10)
axC.set_title('C  Signal vs. negative control\n(cell-level, informative cohorts)',
              fontsize=11, fontweight='bold', loc='left')
axC.legend(frameon=False, fontsize=9)

fig.suptitle('Figure 3  Cross-cohort reproducibility of the (F6, F7) fibroblast state axis in diabetic foot ulcer',
             fontsize=12.5, fontweight='bold', y=0.985)
fig.savefig(OUT / 'Fig3_crosscohort_4cohorts.png', dpi=300, bbox_inches='tight',
            facecolor='white')
print('[saved] Fig3_crosscohort_4cohorts.png')
plt.close()
