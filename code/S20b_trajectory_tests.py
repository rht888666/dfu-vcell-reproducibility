# -*- coding: utf-8 -*-
"""
S20b: 纵向轨迹的关键补充检验
  1. 双侧置换 p（修正 S20 中单侧检验的解读歧义）
  2. F6 早期达峰 vs F7 晚期达峰 —— 时序正交性检验（核心！）
  3. 用 GSE241132 的 FB 亚型注释（第三方）交叉验证时序
  4. 正常愈合"修复方向向量" vs DFU 位移的置换零分布（W6 的显著性）
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
G = ROOT / "10_DFU_Dataset" / "gse241132"
RNG = np.random.default_rng(7)

ORDER = ['Skin', 'Wound1', 'Wound7', 'Wound30']
LABEL = {'Skin': 'D0', 'Wound1': 'D1', 'Wound7': 'D7', 'Wound30': 'D30'}
M = {
    'F1': ['COL13A1', 'COL23A1', 'SFRP2', 'CRABP1', 'WNT5A'],
    'F2': ['PI16', 'CD34', 'MFAP5', 'DPT', 'COL15A1', 'SERPINF1'],
    'F3': ['CCL19', 'CD74', 'HLA-DRA', 'CCL21', 'TNFSF13B', 'IL7R'],
    'F5': ['SCN7A', 'SOX10', 'PLP1', 'NGFR', 'S100B'],
    'F6': ['CXCL8', 'MMP3', 'CHI3L1', 'CXCL1', 'CCL2'],
    'F6i': ['IL6', 'CXCL1', 'CCL2', 'CXCL8', 'VEGFA'],
    'F7': ['LRRC15', 'SFRP4', 'ASPN', 'POSTN', 'CTHRC1', 'ACTA2', 'TAGLN'],
    'F7e': ['COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'POSTN', 'LOX'],
}
ALLG = sorted({g for v in M.values() for g in v} | {'DCN', 'LUM', 'PDGFRA', 'COL1A1', 'COL1A2', 'FAP', 'THY1'})

LOG = []
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); LOG.append(s)

log('=' * 80)
log('S20b: 纵向轨迹关键补充检验')
log('=' * 80)


# ---- 重新载入（流式，只取 marker）----
def read_10x_mtx(mtx, genes_f, barcodes_f, wanted):
    import gzip
    g = pd.read_csv(genes_f, sep='\t', header=None)
    gname = (g[1] if g.shape[1] >= 2 else g[0]).astype(str).values
    want = set(wanted); keep_rows = {}
    for i, nm in enumerate(gname):
        if nm in want and nm not in keep_rows:
            keep_rows[nm] = i + 1
    kept = list(keep_rows.keys()); row2col = {r: k for k, r in enumerate(keep_rows.values())}
    bcs = pd.read_csv(barcodes_f, sep='\t', header=None)[0].astype(str).values
    rows, cols, vals = [], [], []
    with gzip.open(mtx, 'rt') as fh:
        for line in fh:
            if line.startswith('%'):
                continue
            p = line.split()
            if len(p) < 3:
                continue
            r = int(p[0])
            if r in row2col:
                rows.append(row2col[r]); cols.append(int(p[1]) - 1); vals.append(float(p[2]))
    import scipy.sparse as sp
    Mm = sp.csr_matrix((vals, (cols, rows)), shape=(len(bcs), len(kept)))
    return pd.DataFrame(Mm.toarray(), columns=kept), bcs


unz = G / 'extracted' / 'unzipped'
samples = {}
for d in sorted(unz.iterdir()):
    if not d.is_dir():
        continue
    mt = list(d.glob('*matrix.mtx*'))
    gn = list(d.glob('*features.tsv*')) + list(d.glob('*genes.tsv*'))
    bc = list(d.glob('*barcodes.tsv*'))
    if mt and gn and bc:
        samples[d.name] = (mt[0], gn[0], bc[0])
log(f'样本: {sorted(samples.keys())}')

frames = []
for sid, (mt, gn, bc) in samples.items():
    df, bcs = read_10x_mtx(mt, gn, bc, ALLG)
    df['sample'] = sid
    df['barcode'] = bcs
    frames.append(df)
expr = pd.concat(frames, ignore_index=True)
expr['bc_key'] = expr['sample'].astype(str) + '_' + expr['barcode'].astype(str)

meta = pd.read_csv(G / 'cell_metadata.txt.gz', sep='\t')
meta['bc_key'] = meta['barcode'].astype(str)
mg = expr.merge(meta[['bc_key', 'Patient', 'Condition', 'newCellTypes', 'newMainCellTypes', 'nCount_RNA']],
                on='bc_key', how='inner')
log(f'合并: {len(mg):,}')
fib = mg[mg['newMainCellTypes'] == 'Fibroblast'].copy()
log(f'成纤维: {len(fib):,}')

fib['logumi'] = np.log1p(fib['nCount_RNA'].astype(float))
logu = fib['logumi'].values; logu_c = logu - logu.mean()
tmap = {t: i for i, t in enumerate(ORDER)}
fib['t_ord'] = fib['Condition'].map(tmap)


def zscore(df, genes):
    gs = [g for g in genes if g in df.columns]
    sub = df[gs].astype(float)
    return ((sub - sub.mean(0)) / (sub.std(0) + 1e-9)).mean(1).values


def resid(v):
    ok = ~np.isnan(v)
    b = np.polyfit(logu_c[ok], v[ok], 1)
    return v - (b[0] * logu_c + b[1])


for prog, genes in M.items():
    fib[f'sc_{prog}'] = resid(zscore(fib, genes))


# =====================================================
# 1. 双侧置换 p
# =====================================================
log('\n' + '=' * 80)
log('【1】双侧置换检验（供体内时间重排，nperm=2000）')
log('=' * 80)

def perm_two_sided(df, col, nperm=2000):
    v = df[col].values; ok = ~np.isnan(v)
    obs, _ = spearmanr(df.loc[ok, 't_ord'], v[ok])
    cnt = 0
    donors = df['Patient'].unique()
    for _ in range(nperm):
        d2 = df.copy()
        for dd in donors:
            idx = d2.index[d2['Patient'] == dd]
            d2.loc[idx, 't_ord'] = RNG.permutation(d2.loc[idx, 't_ord'].values)
        vv = d2[col].values; oo = ~np.isnan(vv)
        r, _ = spearmanr(d2.loc[oo, 't_ord'], vv[oo])
        if abs(r) >= abs(obs):
            cnt += 1
    return float(obs), float((cnt + 1) / (nperm + 1))

two = {}
for prog in M:
    rho, p = perm_two_sided(fib, f'sc_{prog}')
    two[prog] = {'rho': rho, 'p_two_sided': p}
    log(f'  {prog:6s} rho={rho:+.4f}  p_two_sided={p:.4f}')


# =====================================================
# 2. 时序正交性：F6 是否早于 F7 达峰？
# =====================================================
log('\n' + '=' * 80)
log('【2】F6/F7 时序正交性（核心）：各时间点均值 + 达峰位置')
log('=' * 80)
tp = fib.groupby('Condition')[['sc_F6', 'sc_F6i', 'sc_F7', 'sc_F7e', 'sc_F1', 'sc_F5']].mean()
tp = tp.reindex([t for t in ORDER if t in tp.index])
log(tp.round(4).to_string())

peakF6 = tp['sc_F6'].idxmax(); peakF7 = tp['sc_F7'].idxmax()
peakF6e = tp['sc_F6i'].idxmax(); peakF7e = tp['sc_F7e'].idxmax()
log(f'\n  F6 达峰: {LABEL[peakF6]}   F6i 达峰: {LABEL[peakF6e]}')
log(f'  F7 达峰: {LABEL[peakF7]}   F7e 达峰: {LABEL[peakF7e]}')
log(f'  🔴 F6 早于 F7 达峰? {tmap[peakF6] < tmap[peakF7]}')

# 供体级：每个供体的 F6/F7 达峰时间
log('\n  供体级达峰:')
dpeak = {}
for d, g in fib.groupby('Patient'):
    t = g.groupby('Condition')[['sc_F6', 'sc_F7']].mean()
    t = t.reindex([x for x in ORDER if x in t.index])
    dpeak[d] = {'F6_peak': LABEL[t['sc_F6'].idxmax()], 'F7_peak': LABEL[t['sc_F7'].idxmax()]}
    log(f'    {d}: F6@{dpeak[d]["F6_peak"]}  F7@{dpeak[d]["F7_peak"]}')
n_f6_early = sum(1 for d in dpeak if tmap[[k for k, v in LABEL.items() if v == dpeak[d]['F6_peak']][0]] <
                 tmap[[k for k, v in LABEL.items() if v == dpeak[d]['F7_peak']][0]])
log(f'  → {n_f6_early}/3 供体 F6 早于 F7 达峰')


# =====================================================
# 3. 第三方 FB 亚型注释的时序（独立于我们的 marker）
# =====================================================
log('\n' + '=' * 80)
log('【3】第三方 FB 亚型注释的时间轨迹（原作者注释，独立验证）')
log('=' * 80)
sub = fib.pivot_table(index='newCellTypes', columns='Condition', values='bc_key', aggfunc='count').fillna(0)
sub = sub.reindex(columns=[c for c in ORDER if c in sub.columns])
log('细胞数:')
log(sub.astype(int).to_string())
subp = sub / sub.sum(0) * 100
log('\n占比 %:')
log(subp.round(2).to_string())
# FB-I 的时间相关性
for st in subp.index:
    v = subp.loc[st].values
    rho, _ = spearmanr(range(len(v)), v)
    log(f'  {st:10s} 时间 rho = {rho:+.3f}')


# =====================================================
# 4. W6 置换零分布
# =====================================================
log('\n' + '=' * 80)
log('【4】W6 修复方向向量余弦相似度的置换零分布')
log('=' * 80)
v_normal = np.array([tp.loc['Wound30', 'sc_F6'] - tp.loc['Skin', 'sc_F6'],
                     tp.loc['Wound30', 'sc_F7'] - tp.loc['Skin', 'sc_F7']])
v_normal_u = v_normal / np.linalg.norm(v_normal)
log(f'  正常愈合方向 (F6,F7): {v_normal.round(4)}')

sd = pd.read_csv(OUT / 'S14_state_space_patients.csv')
grp = 'group'
h = sd[sd[grp].astype(str).str.contains('Heal') & ~sd[grp].astype(str).str.contains('Non')]
nh = sd[sd[grp].astype(str).str.contains('Non')]
v_dfu = np.array([h['F6'].mean() - nh['F6'].mean(), h['F7'].mean() - nh['F7'].mean()])
v_dfu_u = v_dfu / np.linalg.norm(v_dfu)
cos_obs = float(np.dot(v_dfu_u, v_normal_u))
log(f'  DFU 位移 (F6,F7): {v_dfu.round(4)}')
log(f'  cos = {cos_obs:+.4f}')

# 零分布：随机打乱 Healer/NonHealer 标签，重算 DFU 位移方向
subj_F = sd[['F6', 'F7']].values
y_true = np.array([1 if 'Heal' in str(x) and 'Non' not in str(x) else 0 for x in sd[grp]])
cnt = 0; cos_dist = []
for _ in range(5000):
    yp = RNG.permutation(y_true)
    vp = subj_F[yp == 1].mean(0) - subj_F[yp == 0].mean(0)
    nu = np.linalg.norm(vp)
    if nu < 1e-9:
        continue
    c = float(np.dot(vp / nu, v_normal_u))
    cos_dist.append(c)
    if c >= cos_obs:
        cnt += 1
cos_dist = np.array(cos_dist)
p_w6 = (cnt + 1) / (len(cos_dist) + 1)
log(f'  置换零分布: mean={cos_dist.mean():+.3f}, sd={cos_dist.std():.3f}, '
    f'p95={np.percentile(cos_dist,95):+.3f}')
log(f'  🔴 W6 单侧 p_perm = {p_w6:.4f}   (cos_obs={cos_obs:+.4f})')

# =====================================================
# 输出
# =====================================================
RES = {'two_sided_perm': two,
       'timecourse': tp.round(6).to_dict(),
       'peaks': {'F6': LABEL[peakF6], 'F6i': LABEL[peakF6e],
                 'F7': LABEL[peakF7], 'F7e': LABEL[peakF7e]},
       'donor_peaks': dpeak,
       'n_donors_F6_before_F7': int(n_f6_early),
       'fib_subtype_timecourse_pct': subp.round(4).to_dict(),
       'W6': {'v_normal': v_normal.tolist(), 'v_dfu': v_dfu.tolist(),
              'cos_obs': cos_obs, 'p_perm': float(p_w6),
              'null_mean': float(cos_dist.mean()), 'null_sd': float(cos_dist.std()),
              'n_perm': int(len(cos_dist))}}
with open(OUT / 'S20b_trajectory_tests.json', 'w', encoding='utf-8') as f:
    json.dump(RES, f, indent=2, ensure_ascii=False, default=str)
Path(OUT / 'S20b_run.log').write_text('\n'.join(LOG), encoding='utf-8')
log('\n[saved] S20b_trajectory_tests.json / S20b_run.log')
log('✅ S20b 完成')
