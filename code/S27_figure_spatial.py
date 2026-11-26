# -*- coding: utf-8 -*-
"""
S27: 空间证据配图（Figure 5）—— Nature 级 4 面板
  A. 四条件下 F6/F7/F7e 得分（空间版纵向轨迹，含供体散点）
  B. 各程序的趋势 rho vs 随机基因集零分布（展示 F7/F7e 超零，F1/F5 在零内）
  C. F7 在 17 个空间注释类别的得分（展示 FB1 特异富集）
  D. F7 在 FB1 的富集置换零分布
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
RNG = np.random.default_rng(20261006)

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9,
    'axes.linewidth': 0.8, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

SP = pd.read_csv(OUT / "S26_spatial_scores.csv")
ord_map = {'Skin':0,'Wound1':1,'Wound7':2,'Wound30':3}
SP['cond_ord'] = SP['Condition'].map(ord_map)
COND = ['Skin','Wound1','Wound7','Wound30']
LAB = ['D0\nSkin','D1\nWound','D7\nWound','D30\nWound']
donors = sorted(SP['Donor'].dropna().unique())
DCOL = {'Donor1':'#3B6FB6','Donor2':'#D1495B','Donor3':'#5C946E','Donor4':'#C77DFF'}

with open(OUT/"S26b_global_gradient.json", encoding='utf-8') as f:
    GJ = json.load(f)
with open(OUT/"S26_spatial.json", encoding='utf-8') as f:
    SJ = json.load(f)

fig = plt.figure(figsize=(14, 9.2))
gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.24)

# --- A ---
axA = fig.add_subplot(gs[0,0])
progs = ['F6','F7','F7e']
pcol = {'F6':'#D1495B', 'F7':'#2A6F97', 'F7e':'#5C946E'}
x = np.arange(4)
for p in progs:
    y = [SP[SP['Condition']==c][p].mean() for c in COND]
    axA.plot(x, y, '-o', color=pcol[p], lw=2.2, ms=7, label=p, zorder=3)
# 供体散点
for d in donors:
    for p in progs:
        yv = [SP[(SP['Condition']==c)&(SP['Donor']==d)][p].mean() for c in COND]
        axA.plot(x, yv, 'o', color=DCOL.get(d,'gray'), ms=3, alpha=0.45, zorder=2)
axA.axhline(0, color='k', lw=0.6, ls=':', alpha=0.5)
axA.set_xticks(x); axA.set_xticklabels(LAB, fontsize=8.5)
axA.set_ylabel('Program score (spot-level z)', fontsize=9.5)
axA.set_title('A  Spatial longitudinal trajectory', fontsize=11, fontweight='bold', loc='left')
axA.legend(frameon=False, fontsize=9, loc='upper left', ncol=3)
# 标注（R15：数值从 JSON 实时读取，避免硬编码过期）
import json as _json
with open(OUT / "S26b_global_gradient.json", encoding='utf-8') as _f:
    _GJ = _json.load(_f)
_ep = _GJ['empirical_p']
axA.annotate(f"F6 declines with healing\n"
             f"(4/4 donors, rho={_ep['F6']['obs']:+.2f}, n.s.)",
             xy=(3, SP[SP['Condition'] == 'Wound30']['F6'].mean()),
             xytext=(1.35, -0.62), fontsize=7.2, color='#D1495B',
             arrowprops=dict(arrowstyle='->', color='#D1495B', lw=0.8))
axA.annotate(f"F7e accumulates (4/4 donors)\n"
             f"rho={_ep['F7e']['obs']:+.2f}, p={_ep['F7e']['p_emp']:.3f}",
             xy=(3, SP[SP['Condition'] == 'Wound30']['F7e'].mean()),
             xytext=(1.55, 0.90), fontsize=7.2, color='#2A6F97',
             arrowprops=dict(arrowstyle='->', color='#2A6F97', lw=0.8))

# --- B ---
axB = fig.add_subplot(gs[0,1])
null_m, null_s = GJ['null_rand']['mean'], GJ['null_rand']['sd']
xs = np.linspace(null_m-4*null_s, null_m+4*null_s, 300)
axB.fill_between(xs, 0, np.exp(-(xs-null_m)**2/(2*null_s**2)), color='#CCCCCC', alpha=0.65, zorder=1)
axB.axvline(GJ['null_rand']['p95'], color='gray', lw=1.0, ls='--', zorder=2)
order = ['F6','F5','F1','F7','F7e']
ypos = np.arange(len(order))
for i, p in enumerate(order):
    obs = GJ['empirical_p'][p]['obs']
    pe = GJ['empirical_p'][p]['p_emp']
    sig = pe < 0.05
    c = '#2A6F97' if sig else '#999999'
    axB.plot([obs, obs], [i-0.28, i+0.28], color=c, lw=3.2, solid_capstyle='round', zorder=4)
    axB.text(obs, i+0.36, f'p={pe:.3f}', ha='center', fontsize=7.3, color=c)
    axB.plot(null_m, i, '|', color='#777', ms=11, mew=1.4, zorder=5)
axB.set_yticks(ypos); axB.set_yticklabels(order, fontsize=9)
axB.set_xlabel('Condition-trend rho (donor-mean)', fontsize=9)
axB.set_title('B  F7/F7e exceed random-gene null', fontsize=11, fontweight='bold', loc='left')
axB.set_xlim(min(GJ['empirical_p'][p]['obs'] for p in order)-0.08, max(GJ['empirical_p'][p]['obs'] for p in order)+0.10)
axB.legend(handles=[Patch(facecolor='#2A6F97', label='exceeds null (p<0.05)'),
                    Patch(facecolor='#999999', label='within null (n.s.)'),
                    Patch(facecolor='#CCCCCC', label='random-gene null ±4sd')],
           frameon=False, fontsize=8, loc='lower right')

# --- C ---
axC = fig.add_subplot(gs[1,0])
AN = SP.groupby('AnnoType').agg(n=('F7','size'), F7=('F7','mean'), F6=('F6','mean')).sort_values('F7')
AN = AN[AN['n']>=100]
yp = np.arange(len(AN))
axC.barh(yp, AN['F7'].values, color=['#2A6F97' if v>0.3 else '#A9C0D4' for v in AN['F7'].values], height=0.72)
axC.set_yticks(yp); axC.set_yticklabels(AN.index, fontsize=8)
axC.axvline(0, color='k', lw=0.6)
axC.set_xlabel('Mean F7 score (spot-level)', fontsize=9)
axC.set_title('C  F7 localizes to FB1 spots', fontsize=11, fontweight='bold', loc='left')
for i,(v,n) in enumerate(zip(AN['F7'].values, AN['n'].values)):
    if v > 0.3:
        axC.text(v+0.015, i, f'{v:.2f} (n={n})', va='center', fontsize=7.3, color='#2A6F97', fontweight='bold')

# --- D ---
axD = fig.add_subplot(gs[1,1])
# FB1 富集零分布
is_fb1 = (SP['AnnoType']=='FB1').values
zcols = [c for c in SP.columns if c.endswith('_z')]
log_n = np.log1p(SP['total_counts'].values.astype(float))
lc = log_n-log_n.mean(); den=(lc**2).sum()+1e-12
GZ = {}
for c in zcols:
    v = SP[c].values; beta=(lc*(v-v.mean())).sum()/den
    GZ[c[:-2]] = v-beta*lc
GZ = pd.DataFrame(GZ)
F7g = ['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN']
f7v = GZ[[g for g in F7g if g in GZ.columns]].mean(1).values
obs = f7v[is_fb1].mean()-f7v[~is_fb1].mean()
nullfb = []
for _ in range(3000):
    gs_ = RNG.choice(list(GZ.columns), 7, replace=False)
    s = GZ[gs_].mean(1).values
    nullfb.append(s[is_fb1].mean()-s[~is_fb1].mean())
nullfb = np.array(nullfb)
axD.hist(nullfb, bins=45, color='#CCCCCC', edgecolor='none', density=True, label='random 7-gene sets')
axD.axvline(obs, color='#2A6F97', lw=2.6, label=f'F7 (obs={obs:+.3f}, p={GJ["fb1_enrichment"]["p_specific"]:.3f})')
axD.set_xlabel('FB1 enrichment (Δscore, FB1 − non-FB1)', fontsize=9)
axD.set_ylabel('Density', fontsize=9)
axD.set_title('D  FB1 enrichment is spatially specific', fontsize=11, fontweight='bold', loc='left')
axD.legend(frameon=False, fontsize=8.5)

fig.suptitle('Spatial transcriptomics (GSE241124, 4 donors × 4 timepoints, 16 sections, 22,915 spots) '
             'places the fibroblast state transition in tissue space',
             fontsize=12.0, fontweight='bold', y=0.985)
fig.text(0.5, 0.005,
         'Visium spots ≈ 20 cells; scores are tissue-region programs, not single-cell states.  '
         'Depth-residualized on log total-counts.  Permutation n=2000–3000.',
         ha='center', fontsize=7.5, color='#555')

fig.savefig(OUT/"Fig5_spatial_validation.png", dpi=300, facecolor='white')
plt.close(fig)
print('saved Fig5_spatial_validation.png')
