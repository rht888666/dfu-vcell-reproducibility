# -*- coding: utf-8 -*-
"""S50: regenerate Figures 1-4 from harmonized reanalysis outputs (single source of truth)."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 9, 'axes.linewidth': 0.8, 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 150})
OUT = Path('revision_work/reanalysis/figures_v2'); OUT.mkdir(parents=True, exist_ok=True)
TAB = Path('revision_work/github_update/dfu-vcell-reproducibility/processed_tables')
REA = Path('revision_work/reanalysis')
M = json.load(open(REA / 'MASTER_numbers.json'))
DC = {'PWH26': '#1f77b4', 'PWH27': '#d62728', 'PWH28': '#2ca02c'}
DAYS = ['D0', 'D1', 'D7', 'D30']
XC = {'D0': 0, 'D1': 1, 'D7': 7, 'D30': 30}

# ================= FIGURE 1: longitudinal =================
df = pd.read_csv(REA / 'S40_donor_timepoint_programs.csv')

def dmat(prog):
    return {d: [df[(df.donor == d) & (df.tp == t)][prog].values[0] for t in DAYS] for d in DC}

fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.3))
# A state space F6 vs F7e
ax = axes[0]
for d in DC:
    f6 = dmat('F6')[d]; f7 = dmat('F7e')[d]
    for i in range(3):
        ax.annotate('', xy=(f6[i + 1], f7[i + 1]), xytext=(f6[i], f7[i]),
                    arrowprops=dict(arrowstyle='->', color=DC[d], lw=1.4, alpha=0.85))
    ax.scatter(f6, f7, c=DC[d], s=28, zorder=5, label=d)
    for i, t in enumerate(DAYS):
        ax.annotate(t, (f6[i], f7[i]), textcoords='offset points', xytext=(4, 4), fontsize=7, color=DC[d])
ax.axhline(0, color='grey', lw=0.5, ls=':'); ax.axvline(0, color='grey', lw=0.5, ls=':')
ax.set_xlabel('F6 inflammatory score (z)'); ax.set_ylabel('F7e matrix score (z)')
ax.set_title('A  Donor state-space trajectories', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=8)
# B time courses
ax = axes[1]
for prog, col, lab in [('F6', '#e41a1c', 'F6 (inflammatory)'), ('F7e', '#08519c', 'F7e (matrix)')]:
    means = []
    for t in DAYS:
        vals = [df[(df.donor == d) & (df.tp == t)][prog].values[0] for d in DC]
        means.append(np.mean(vals))
        ax.scatter([XC[t]] * 3, vals, color=col, s=12, alpha=0.55, zorder=3)
    ax.plot([XC[t] for t in DAYS], means, color=col, lw=2.2, marker='o', ms=5, label=lab, zorder=4)
ax.axhline(0, color='grey', lw=0.5, ls=':')
ax.set_xlabel('Day after wounding'); ax.set_ylabel('program score (z)')
ax.set_title('B  Donor-level time courses', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=8)
ax.set_xticks([0, 1, 7, 30])
# C D0->D30 change per program
ax = axes[2]
progs = ['F1', 'F5', 'F6', 'F6i', 'F7', 'F7noCT', 'F7e']
ylab = {'F7noCT': 'F7\u2212CT'}
for yi, prog in enumerate(progs):
    deltas = []
    for d in DC:
        v = dmat(prog)[d]; deltas.append(v[3] - v[0])
        ax.scatter(deltas[-1], yi, color=DC[d], s=22, zorder=3)
    ax.scatter(np.mean(deltas), yi, marker='D', color='black', s=34, zorder=4)
ax.axvline(0, color='grey', lw=0.7, ls='--')
ax.set_yticks(range(len(progs))); ax.set_yticklabels([ylab.get(p, p) for p in progs])
ax.set_xlabel('\u0394 score D30 \u2212 D0 (z, per donor)')
ax.set_title('C  Net day-0\u219230 displacement', loc='left', fontsize=10)
# D dispersion dumbbell
ax = axes[3]
ang = M['longitudinal']['angular']
order = ['F7', 'F7noCT', 'F7e', 'F7epCT']
labels = ['F7 (+ACTA2/TAGLN)', 'F7 \u2212 ACTA2/TAGLN', 'F7e (core matrix)', 'F7e + ACTA2/TAGLN']
colors = ['#cb181d', '#fb6a4a', '#08519c', '#6baed6']
for yi, (p, lab, c) in enumerate(zip(order, labels, colors)):
    a = ang[p]
    ax.scatter(a['angles'], [yi] * 3, color=c, s=30, alpha=0.75, zorder=3)
    lo, hi = min(a['angles']), max(a['angles'])
    ax.plot([lo, hi], [yi, yi], color=c, lw=2.5, zorder=2)
    ax.text(hi + 1.2, yi, f"range {a['range']:.1f}\u00b0", va='center', fontsize=7.5)
ax.set_yticks(range(4)); ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel('donor resultant angle (\u00b0)')
ax.set_title('D  Angular dispersion & ablation', loc='left', fontsize=10)
ax.set_xlim(98, 145)
plt.tight_layout()
plt.savefig(OUT / 'Figure1_longitudinal.png', dpi=300); plt.savefig(OUT / 'Figure1_longitudinal.pdf')
plt.close()
print('fig1 done')

# ================= FIGURE 2: spatial =================
SP = pd.read_csv(TAB / 'S26_spatial_scores.csv')
PROG = {'F6': ['CXCL8', 'MMP3', 'CHI3L1', 'CXCL1', 'CCL2'],
        'F7': ['LRRC15', 'SFRP4', 'ASPN', 'POSTN', 'CTHRC1', 'ACTA2', 'TAGLN'],
        'F7e': ['COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'POSTN', 'LOX'],
        'F7noCT': ['LRRC15', 'SFRP4', 'ASPN', 'POSTN', 'CTHRC1'],
        'F7epCT': ['COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'POSTN', 'LOX', 'ACTA2', 'TAGLN']}
logn = np.log1p(SP['total_counts'].values.astype(float))
lc = logn - logn.mean(); den = (lc ** 2).sum() + 1e-12
GENES = sorted({g for v in PROG.values() for g in v})
GZ = {}
for g in GENES:
    v = SP[g + '_z'].values.astype(float)
    beta = (lc * (v - v.mean())).sum() / den
    GZ[g] = v - beta * lc
for p, gs in PROG.items():
    SP[p + '_gz'] = np.mean([GZ[g] for g in gs], axis=0)
DON = {'Donor1': '#1b9e77', 'Donor2': '#d95f02', 'Donor3': '#7570b3', 'Donor4': '#e7298a'}
CONDS = ['Skin', 'Wound1', 'Wound7', 'Wound30']
FIBCAT = ['FB1', 'FB2', 'Dermis1', 'Dermis2', 'Papillary_dermis', 'Wound_edges', 'Sweatgland_FB']

def fibsel(mask):
    return mask & (SP.AnnoType.isin(FIBCAT))

fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.3))
# A temporal donor-level
ax = axes[0]
isfib = SP.AnnoType.isin(FIBCAT)
for prog, col in [('F7', '#cb181d'), ('F7e', '#08519c')]:
    dmean = {cd: [] for cd in CONDS}
    for dn in DON:
        vals = []
        for cd in CONDS:
            sel = SP[(SP.Donor == dn) & (SP.Condition == cd) & isfib]
            vals.append(sel[prog + '_gz'].mean() if len(sel) else np.nan)
        ax.plot(range(4), vals, color=col, alpha=0.35, lw=1)
        for cd, v in zip(CONDS, vals):
            dmean[cd].append(v)
    mv = [np.nanmean(dmean[cd]) for cd in CONDS]
    ax.plot(range(4), mv, color=col, lw=2.4, marker='o', ms=5, label=prog)
ax.set_xticks(range(4)); ax.set_xticklabels(['D0', 'D1', 'D7', 'D30'], fontsize=8)
ax.axhline(0, color='grey', lw=0.5, ls=':')
ax.set_ylabel('depth-residualized z'); ax.set_title('A  Temporal dynamics (fibroblast compartments)', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=8)
# B compartment heatmap at Wound7
ax = axes[1]
cats = ['Wound_edges', 'FB1', 'FB2', 'Dermis1', 'Dermis2', 'Papillary_dermis', 'Smooth_muscle', 'Immune']
hm = np.full((len(cats), 3), np.nan)
w7 = SP[SP.Condition == 'Wound7']
for ri, cat in enumerate(cats):
    for ci, prog in enumerate(['F6', 'F7', 'F7e']):
        vals = []
        for dn in DON:
            sel = w7[(w7.Donor == dn) & (w7.AnnoType == cat)]
            if len(sel) >= 5: vals.append(sel[prog + '_gz'].mean())
        if vals: hm[ri, ci] = np.mean(vals)
im = ax.imshow(hm, cmap='RdBu_r', vmin=-0.6, vmax=0.6, aspect='auto')
ax.set_xticks(range(3)); ax.set_xticklabels(['F6', 'F7', 'F7e'])
ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats, fontsize=7.5)
for ri in range(len(cats)):
    for ci in range(3):
        if not np.isnan(hm[ri, ci]):
            ax.text(ci, ri, f'{hm[ri, ci]:+.2f}', ha='center', va='center', fontsize=7)
ax.set_title('B  Compartment scores at day 7', loc='left', fontsize=10)
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
# C distance gradient day 7
ax = axes[2]
try:
    ed = pd.read_csv(TAB / 'S29_edge_dist.csv')
    mrg = SP.merge(ed[['barcode', 'h5_stem', 'edge_dist']], on=['barcode', 'h5_stem'], how='left')
    w7d = mrg[(mrg.Condition == 'Wound7') & mrg.edge_dist.notna() & mrg.AnnoType.isin(FIBCAT)].copy()
    bins = np.quantile(w7d.edge_dist, [0, .25, .5, .75, 1.0])
    bins = np.unique(bins)
    w7d['db'] = pd.cut(w7d.edge_dist, bins, include_lowest=True)
    for prog, col in [('F7', '#cb181d'), ('F7e', '#08519c')]:
        for dn in DON:
            sub = w7d[w7d.Donor == dn]
            g = sub.groupby('db', observed=True)[prog + '_gz'].mean()
            if len(g) >= 3: ax.plot(range(len(g)), g.values, color=col, alpha=0.35, lw=1)
        g = w7d.groupby('db', observed=True)[prog + '_gz'].mean()
        ax.plot(range(len(g)), g.values, color=col, lw=2.4, marker='o', ms=5, label=prog)
    nb = len(bins) - 1
    ax.set_xticks(range(nb)); ax.set_xticklabels(['nearest', 'near', 'far', 'farthest'][:nb], fontsize=7.5)
except Exception as e:
    ax.text(0.5, 0.5, 'distance data unavailable', transform=ax.transAxes, ha='center')
    print('panel2C fallback:', e)
ax.axhline(0, color='grey', lw=0.5, ls=':')
ax.set_xlabel('distance to wound edge (quartile)'); ax.set_ylabel('depth-residualized z')
ax.set_title('C  Day-7 gradient from wound edge', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=8)
# D smooth muscle ablation per donor
ax = axes[3]
sm = M['spatial']['smooth_muscle_per_donor']
panels = ['F7', 'F7noCT', 'F7e', 'F7epCT']
x = np.arange(4); w = 0.2
for i, (dn, col) in enumerate(DON.items()):
    vals = [sm.get(dn, {}).get(p, np.nan) if isinstance(sm.get(dn), dict) else np.nan for p in panels]
    ax.bar(x + (i - 1.5) * w, vals, w, color=col, label=dn)
ax.axhline(0, color='black', lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(['F7', 'F7\u2212CT', 'F7e', 'F7e+CT'], fontsize=8)
ax.set_ylabel('smooth-muscle spot score (z)')
ax.set_title('D  Smooth-muscle enrichment & ablation', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=7, ncol=2)
plt.tight_layout()
plt.savefig(OUT / 'Figure2_spatial.png', dpi=300); plt.savefig(OUT / 'Figure2_spatial.pdf')
plt.close()
print('fig2 done')

# ================= FIGURE 3: clinical validation =================
s41 = pd.read_csv(REA / 'S41_library_table_primary_S25_gate.csv')
s12 = pd.read_csv(TAB / 'S12_GSE231643_samples.csv')
s21 = pd.read_csv(TAB / 'S21_GSE223964_samples.csv')
g16 = M['gse165816']['primary_S25_gate']
co23 = M['cohorts']['GSE231643']
co22 = M['cohorts']['GSE223964']
fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.3))
# A depth confound scatter
ax = axes[0]
for grp, col in [('Healer', '#2171b5'), ('NonHealer', '#cb181d')]:
    sub = s41[s41.group == grp]
    ax.scatter(np.log10(sub.umi), sub.F7_raw, c=col, s=34, alpha=0.85, label=grp)
ax.set_xlabel('log10 UMI per library'); ax.set_ylabel('F7 raw score (z)')
ax.set_title('A  GSE165816: score vs sequencing depth', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=8)
# B forest: exact two-sided effects
ax = axes[1]
tests = g16['tests']
rows = [('F7 raw', tests['F7_raw']), ('F7 linear-resid', tests['F7_lin']), ('F7 rank-resid', tests['F7_rank']),
        ('F7e raw', tests['F7e_raw'])]
for yi, (lab, t) in enumerate(rows):
    p = t['p_two_sided_exact']; d = t['diff']
    if 'rank' in lab:
        d = d / 36119.0
    c = '#08519c' if p < 0.05 else '#999999'
    ax.scatter(d, yi, color=c, s=60, zorder=3)
    ax.text(d + 0.02 if d >= 0 else d - 0.02, yi, f"p={p:.3f}", va='center',
            ha='left' if d >= 0 else 'right', fontsize=7.5)
ax.axvline(0, color='black', lw=0.8, ls='--')
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
ax.set_xlabel('Healer \u2212 NonHealer (z or scaled rank units)')
ax.set_title('B  GSE165816 exact two-sided tests', loc='left', fontsize=10)
# C GSE231643 programs raw/res
ax = axes[2]
rows = []
for prog in ['F6', 'F7e', 'F1', 'F5']:
    for sfx, sl in [('', 'raw'), ('_res', 'depth-resid')]:
        key = prog + sfx
        if key in co23['tests']:
            rows.append((f'{prog} {sl}', co23['tests'][key]))
for yi, (lab, t) in enumerate(rows):
    p = t['p_two_sided_exact']
    c = '#08519c' if p < 0.05 else '#999999'
    ax.scatter(t['diff'], yi, color=c, s=60, zorder=3)
    ax.text(t['diff'] + 0.02 if t['diff'] >= 0 else t['diff'] - 0.02, yi, f"p={p:.3f}", va='center',
            ha='left' if t['diff'] >= 0 else 'right', fontsize=7.5)
ax.axvline(0, color='black', lw=0.8, ls='--')
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
ax.set_xlabel('Healer \u2212 NonHealer (z)')
ax.set_title('C  GSE231643 exact two-sided tests', loc='left', fontsize=10)
# D GSE223964 + library strip for GSE231643 F6
ax = axes[3]
rows = [(prog, co22['tests'][prog]) for prog in ['F6', 'F7e', 'F1'] if prog in co22['tests']]
for yi, (lab, t) in enumerate(rows):
    p = t['p_two_sided_exact']
    c = '#08519c' if p < 0.05 else '#999999'
    ax.scatter(t['diff'], yi, color=c, s=60, zorder=3)
    ax.text(t['diff'] + 0.02 if t['diff'] >= 0 else t['diff'] - 0.02, yi, f"p={p:.3f}", va='center',
            ha='left' if t['diff'] >= 0 else 'right', fontsize=7.5)
ax.axvline(0, color='black', lw=0.8, ls='--')
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
ax.set_xlabel('diabetic \u2212 non-diabetic ulcer (z)')
ax.set_title('D  GSE223964 exact two-sided tests', loc='left', fontsize=10)
plt.tight_layout()
plt.savefig(OUT / 'Figure3_clinical.png', dpi=300); plt.savefig(OUT / 'Figure3_clinical.pdf')
plt.close()
print('fig3 done')

# ================= FIGURE 4: CRISPRa =================
rk = pd.read_csv(REA / 'S42_full_rankings.csv')
cr = M['crispra']
fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.3))
# A scatter cos AxisA vs AxisM7
ax = axes[0]
top = rk.nsmallest(60, 'rank_AxisA_contractile_inclusive')
ax.scatter(top.cos_AxisA_contractile_inclusive, top.cos_AxisM7_measured_matrix, s=14, color='#999999', alpha=0.6)
keys = ['EGR3', 'EBF2', 'ZNF827']
for tf in keys:
    r = rk[rk.TF == tf].iloc[0]
    ax.scatter(r.cos_AxisA_contractile_inclusive, r.cos_AxisM7_measured_matrix, s=55,
               color='#cb181d' if tf == 'EGR3' else '#08519c', zorder=5)
    ax.annotate(tf, (r.cos_AxisA_contractile_inclusive, r.cos_AxisM7_measured_matrix),
                textcoords='offset points', xytext=(5, 5), fontsize=8, weight='bold')
lim = max(top.cos_AxisA_contractile_inclusive.max(), top.cos_AxisM7_measured_matrix.max()) * 1.1
ax.plot([0, lim], [0, lim], color='grey', lw=0.6, ls=':')
ax.set_xlabel('cosine to Axis-A (contractile-inclusive)')
ax.set_ylabel('cosine to Axis-M7 (measured matrix)')
ax.set_title('A  CRISPRa alignment, two axis definitions', loc='left', fontsize=10)
# B dumbbell key TFs across 3 vectors
ax = axes[1]
tfs = ['EGR3', 'EBF2', 'ZNF827', 'ZNF232', 'KLF4', 'FOSL2']
vecs = ['cos_AxisA_contractile_inclusive', 'cos_AxisM7_measured_matrix', 'cos_F7e5_pure_structural']
vlabs = ['Axis-A', 'Axis-M7', 'F7e5']
for yi, tf in enumerate(tfs):
    r = rk[rk.TF == tf]
    if not len(r): continue
    r = r.iloc[0]
    vals = [r[v] for v in vecs]
    ax.plot(vals, [yi] * 3, color='#bbbbbb', lw=1.5, zorder=2)
    for xi, (v, vv) in enumerate(zip(vals, vecs)):
        ax.scatter(v, yi, s=45, color=['#cb181d', '#08519c', '#2ca02c'][xi], zorder=3)
ax.set_yticks(range(len(tfs))); ax.set_yticklabels(tfs, fontsize=8)
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, ms=7, label=l)
           for c, l in zip(['#cb181d', '#08519c', '#2ca02c'], vlabs)]
ax.legend(handles=handles, frameon=False, fontsize=7.5, loc='lower right')
ax.set_xlabel('cosine similarity to axis')
ax.set_title('B  Key TF alignment across axis definitions', loc='left', fontsize=10)
# C EGR3 perturbation profile
ax = axes[2]
ep = cr['egr3_profile']
unmeas=[k for k,v in ep.items() if v is None]
items = sorted([(k,v) for k,v in ep.items() if v is not None], key=lambda kv: kv[1])
genes = [k for k, v in items]
vals = [v for k, v in items]
cols = ['#cb181d' if g in ('ACTA2', 'TAGLN') else ('#08519c' if g in ('COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'POSTN', 'LOX') else '#bbbbbb') for g in genes]
ax.barh(range(len(genes)), vals, color=cols)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6.5)
ax.axvline(0, color='black', lw=0.8)
ax.set_xlabel('EGR3 perturbation effect (z)')
ax.set_title('C  EGR3 perturbation profile (measured panel genes)', loc='left', fontsize=10)
ax.text(0.02,0.02,'not measured: '+', '.join(unmeas),transform=ax.transAxes,fontsize=6.5,color='#666666')
# D LOGO + axis sensitivity for EGR3
ax = axes[3]
logo = cr['logo']
lv = [v['EGR3_cos'] for v in logo.values() if isinstance(v, dict) and 'EGR3_cos' in v]
ranks = [v.get('EGR3_rank') for v in logo.values() if isinstance(v, dict)]
ax.scatter(np.zeros(len(lv)) + np.linspace(-0.12, 0.12, len(lv)), lv, s=42, color='#08519c', alpha=0.8)
ax.axhline(0, color='grey', lw=0.6, ls=':')
ax.set_xlim(-0.5, 0.5); ax.set_xticks([])
ax.set_ylabel('cosine to Axis-M7 (leave-one-out)')
mn, mx = min(lv), max(lv)
ax.set_title(f"D  EGR3 leave-one-out stability (rank 1 in 7/7; cos {mn:.3f}\u2013{mx:.3f})", loc='left', fontsize=9)
plt.tight_layout()
plt.savefig(OUT / 'Figure4_crispra.png', dpi=300); plt.savefig(OUT / 'Figure4_crispra.pdf')
plt.close()
print('fig4 done')
print('ALL FIGURES DONE')
