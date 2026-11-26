# -*- coding: utf-8 -*-
"""
S23: 纵向轨迹图（Figure 4）—— 这是"虚拟细胞状态转移"的中心视觉
  Panel A: F6 / F6i 时间轨迹（炎症期，D1 达峰后下降）
  Panel B: F7 / F7e 时间轨迹（重塑期，单调上升至 D30）
  Panel C: F6 vs F7 的相位图（D0→D1→D7→D30 的轨迹箭头）—— 核心创新图
  Panel D: 阴性对照 F1/F5 + 第三方 FB-I 亚型占比（独立验证时序）
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['font.size'] = 9

ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"

don = pd.read_csv(OUT / 'S20_donor_trajectory.csv')
LABEL = {'Skin': 'D0', 'Wound1': 'D1', 'Wound7': 'D7', 'Wound30': 'D30'}
ORD = ['Skin', 'Wound1', 'Wound7', 'Wound30']
don['tp'] = don['timepoint'].map(LABEL)
don['t_ord'] = don['timepoint'].map({t: i for i, t in enumerate(ORD)})
don = don.sort_values(['donor', 't_ord'])

fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.6), dpi=300)
axA, axB, axC, axD = axes.ravel()

DCOL = {'PWH26': '#1f77b4', 'PWH27': '#ff7f0e', 'PWH28': '#2ca02c'}
XL = ['D0', 'D1', 'D7', 'D30']

# ---------------- Panel A: F6 / F6i (炎症期) ----------------
for d, g in don.groupby('donor'):
    g = g.sort_values('t_ord')
    axA.plot([LABEL[t] for t in g['timepoint']], g['F6'], '-o', color=DCOL[d],
             lw=1.8, ms=6, alpha=0.85, label=d)
    axA.plot([LABEL[t] for t in g['timepoint']], g['F6i'], '--s', color=DCOL[d],
             lw=1.2, ms=4, alpha=0.5)
meanF6 = don.groupby('timepoint')['F6'].mean().reindex(ORD)
axA.plot(XL, meanF6.values, '-', color='black', lw=3, alpha=0.65, label='mean')
axA.axhline(0, color='gray', lw=0.7, ls=':')
axA.set_title('A  F6 inflammatory-myofibroblast program\npeaks at D1 then declines',
              fontsize=10.5, fontweight='bold', loc='left')
axA.set_ylabel('F6 score (depth-residualized)', fontsize=9.5)
axA.legend(frameon=False, fontsize=8, ncol=2)
axA.annotate('inflammation phase', xy=(0.5, 0.04), xycoords='axes fraction',
             ha='center', fontsize=8, color='#555', style='italic')

# ---------------- Panel B: F7 / F7e (重塑期) ----------------
for d, g in don.groupby('donor'):
    g = g.sort_values('t_ord')
    axB.plot([LABEL[t] for t in g['timepoint']], g['F7'], '-o', color=DCOL[d],
             lw=1.8, ms=6, alpha=0.85, label=d)
    axB.plot([LABEL[t] for t in g['timepoint']], g['F7e'], '--s', color=DCOL[d],
             lw=1.2, ms=4, alpha=0.5)
meanF7 = don.groupby('timepoint')['F7'].mean().reindex(ORD)
axB.plot(XL, meanF7.values, '-', color='black', lw=3, alpha=0.65, label='mean')
axB.axhline(0, color='gray', lw=0.7, ls=':')
axB.set_title('B  F7 myofibroblast program\nrises monotonically to D30',
              fontsize=10.5, fontweight='bold', loc='left')
axB.set_ylabel('F7 score (depth-residualized)', fontsize=9.5)
axB.legend(frameon=False, fontsize=8, ncol=2)
axB.annotate('remodeling phase', xy=(0.5, 0.04), xycoords='axes fraction',
             ha='center', fontsize=8, color='#555', style='italic')

# ---------------- Panel C: (F6,F7) 相位图 —— 核心 ----------------
mp = don.groupby('timepoint')[['F6', 'F7']].mean().reindex(ORD)
pts = mp.values
axC.plot(pts[:, 0], pts[:, 1], '-', color='#8e44ad', lw=2.6, zorder=2, alpha=0.85)
colors_tp = ['#95a5a6', '#e67e22', '#f1c40f', '#c0392b']
for i, (x, y) in enumerate(pts):
    axC.scatter(x, y, s=320, c=colors_tp[i], edgecolor='black', linewidth=1.3, zorder=4)
    axC.annotate(XL[i], (x, y), textcoords='offset points', xytext=(11, 7),
                 fontsize=12, fontweight='bold')
# 供体级小点
for d, g in don.groupby('donor'):
    g = g.sort_values('t_ord')
    axC.scatter(g['F6'], g['F7'], s=40, c=DCOL[d], alpha=0.45, edgecolor='none')
# 方向箭头
axC.annotate('', xy=pts[3], xytext=pts[2],
             arrowprops=dict(arrowstyle='-|>', color='#7d3c98', lw=2.4, alpha=0.9))
axC.axhline(0, color='gray', lw=0.8, ls='--')
axC.axvline(0, color='gray', lw=0.8, ls='--')
axC.set_xlabel('F6 (inflammatory myofibroblast)', fontsize=9.5)
axC.set_ylabel('F7 (myofibroblast)', fontsize=9.5)
axC.set_title('C  Repair-state trajectory in (F6, F7) space\nduring normal human wound healing',
              fontsize=10.5, fontweight='bold', loc='left')
axC.text(0.03, 0.96, 'large dots = cohort mean\nsmall dots = individual donors',
         transform=axC.transAxes, va='top', fontsize=7.8,
         bbox=dict(boxstyle='round', fc='#fbfbfb', ec='gray', lw=0.6))

# ---------------- Panel D: 阴性对照 + 第三方亚型 ----------------
subp = pd.read_csv(OUT / 'S20a_fib_subtype_time.csv', index_col=0)
xin = np.arange(4)
axD.plot(xin, don.groupby('timepoint')['F1'].mean().reindex(ORD).values,
         '-o', color='#7d3c98', lw=2, ms=6, label='F1 (negative ctrl)')
axD.plot(xin, don.groupby('timepoint')['F5'].mean().reindex(ORD).values,
         '-^', color='#16a085', lw=2, ms=6, label='F5 (negative ctrl)')
axD.set_xticks(xin); axD.set_xticklabels(XL)
axD.set_ylabel('Score (depth-residualized)', fontsize=9.5)
axD.axhline(0, color='gray', lw=0.7, ls=':')
axD2 = axD.twinx()
fb1 = subp.loc['FB-I'].reindex(ORD).values
axD2.plot(xin, fb1, '-D', color='#c0392b', lw=2.4, ms=7, label='FB-I subtype (%)')
axD2.set_ylabel('FB-I subtype fraction (%)', fontsize=9.5, color='#c0392b')
axD2.tick_params(axis='y', labelcolor='#c0392b')
axD.set_title('D  Negative controls show no trend;\nthird-party FB-I annotation confirms timing',
              fontsize=10.5, fontweight='bold', loc='left')
h1, l1 = axD.get_legend_handles_labels(); h2, l2 = axD2.get_legend_handles_labels()
axD.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc='lower left')

fig.suptitle('Figure 4  The F6\u2192F7 repair-state transition is the normal trajectory of human wound healing',
             fontsize=12.5, fontweight='bold', y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.965])
fig.savefig(OUT / 'Fig4_longitudinal_trajectory.png', dpi=300, bbox_inches='tight',
            facecolor='white')
print('[saved] Fig4_longitudinal_trajectory.png')
plt.close()
