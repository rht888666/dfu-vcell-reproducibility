# -*- coding: utf-8 -*-
"""
S20: GSE241132 纵向创面愈合 —— F6/F7 状态转移轨迹 + 跨体系方向一致性

预注册假设（写脚本时固定，不得事后修改）：
  W1: F6（炎症肌成纤维）在 D1/D7 达峰（炎症期）
  W2: F7（肌成纤维）在 D7/D30 主导（重塑期）
  W3: 存在 F1 -> F6 -> F7 的中间态序列
  W4: 正常愈合中 F6/F7 程序显著上升（普遍规律）
  W5: 阴性对照 F1/F5 无同向趋势
  W6: 用正常愈合 (F6,F7) 时间轨迹定义"修复方向向量"，
      投影 DFU (GSE165816) Healer-NonHealer 位移，检验同向性

方法：S9-S15 相同 marker 定义 + 门控 + 深度残差校正
统计：细胞级置换（供体内重排）+ 供体级符号一致性 + 深度残差
"""
import json, warnings, tarfile
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
G = ROOT / "10_DFU_Dataset" / "gse241132"
RNG = np.random.default_rng(20261006)

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
ORDER = ['Skin', 'Wound1', 'Wound7', 'Wound30']
LABEL = {'Skin': 'D0', 'Wound1': 'D1', 'Wound7': 'D7', 'Wound30': 'D30'}

LOG = []
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOG.append(s)

log('=' * 80)
log('S20: GSE241132 纵向创面愈合 —— F6/F7 状态转移轨迹')
log('=' * 80)


# =========================================================
# 1. 读 metadata（细胞注释）
# =========================================================
meta = pd.read_csv(G / 'cell_metadata.txt.gz', sep='\t')
log(f'\n[metadata] {meta.shape[0]:,} 细胞 × {meta.shape[1]} 列')
meta['is_fib'] = meta['newMainCellTypes'] == 'Fibroblast'
fib_meta = meta[meta['is_fib']].copy()
log(f'[metadata] Fibroblast: {len(fib_meta):,}')
log(f'[metadata] 供体: {sorted(meta["Patient"].unique())}')
log(f'[metadata] 时间点: {sorted(meta["Condition"].unique())}')
log(f'[metadata] 成纤维亚型: {sorted(fib_meta["newCellTypes"].dropna().unique())}')


# =========================================================
# 2. 加载表达矩阵
# =========================================================
def find_sample_matrices():
    """返回 {sample_id: {'mtx':.., 'genes':.., 'barcodes':..}} 或 {sample: h5path}"""
    ex = G / "extracted"
    tar = G / "GSE241132_RAW.tar"
    files = [f for f in ex.rglob("*") if f.is_file()] if ex.exists() else []
    if not files and tar.exists():
        sz = tar.stat().st_size
        log(f'\n[extract] 解压 GSE241132_RAW.tar ({sz:,} bytes, 预期 775,905,280)')
        if sz < 775_000_000:
            log(f'  ⚠️ 文件不完整！差 {775905280 - sz:,} bytes')
        ex.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar) as t:
            t.extractall(ex)
        files = [f for f in ex.rglob("*") if f.is_file()]
    log(f'[extract] 文件数: {len(files)}  后缀: {dict(Counter(f.suffix for f in files))}')

    # --- 关键：GSE241132 的 tar 内是【每样本一个 zip】，zip 内是 10x mtx ---
    zips = sorted([f for f in files if f.suffix == '.zip'])
    if zips:
        import zipfile
        unz = ex / "unzipped"
        unz.mkdir(parents=True, exist_ok=True)
        for z in zips:
            try:
                with zipfile.ZipFile(z) as zf:
                    # 只解压非 __MACOSX 的成员
                    for m in zf.namelist():
                        if '__MACOSX' in m or m.endswith('/'):
                            continue
                        zf.extract(m, unz)
            except Exception as e:
                log(f'   [zip fail] {z.name}: {e}')
        files = [f for f in unz.rglob("*") if f.is_file()]
        log(f'[extract] 解 zip 后文件数: {len(files)}  后缀: {dict(Counter(f.suffix for f in files))}')

    samples = {}
    # 模式 A: 10x mtx (matrix.mtx.gz + features/genes + barcodes)
    mtxs = sorted([f for f in files if 'matrix.mtx' in f.name])
    for m in mtxs:
        d = m.parent
        bcs = [f for f in d.glob('*barcodes.tsv*')]
        gns = [f for f in d.glob('*features.tsv*')] + [f for f in d.glob('*genes.tsv*')]
        if bcs and gns:
            samples[d.name] = {'mtx': m, 'genes': gns[0], 'barcodes': bcs[0]}
    # 模式 B: 每个样本一个 h5
    h5s = sorted([f for f in files if f.suffix in ('.h5', '.h5ad')])
    for h in h5s:
        samples[h.stem] = {'h5': h}
    return samples, files


samples, all_files = find_sample_matrices()
log(f'\n[matrix] 识别到 {len(samples)} 个样本矩阵')
for k, v in list(samples.items())[:20]:
    log(f'   {k}: {list(v.keys())}')

if not samples:
    log('\n⚠️ 未识别到样本矩阵。列出目录树：')
    for f in all_files[:60]:
        log('   ', f.relative_to(G), f.stat().st_size)
    Path(OUT / 'S20_inventory.txt').write_text('\n'.join(LOG), encoding='utf-8')
    raise SystemExit('需要人工检查 tar 结构')

RESULT = {'meta': {
    'n_cells': int(meta.shape[0]),
    'n_fib': int(len(fib_meta)),
    'donors': sorted(meta['Patient'].unique().tolist()),
    'timepoints': sorted(meta['Condition'].unique().tolist()),
    'fib_subtypes': sorted(fib_meta['newCellTypes'].dropna().unique().tolist()),
    'n_sample_matrices': len(samples),
    'sample_keys': list(samples.keys()),
}}


# =========================================================
# 3. 加载每个样本表达（基因 × 细胞），只保留 marker 基因
# =========================================================
ALLGENES = sorted({g for v in M.values() for g in v} | {'DCN', 'LUM', 'PDGFRA', 'FAP', 'THY1'})
log(f'\n[genes] 目标基因数（marker 并集）: {len(ALLGENES)}')


def read_10x_mtx(mtx, genes_f, barcodes_f, wanted=None):
    """流式读取 10x matrix.mtx.gz，只保留 wanted 基因（内存友好）。
    返回 DataFrame(cells × kept_genes) 与 barcodes。"""
    import gzip
    g = pd.read_csv(genes_f, sep='\t', header=None)
    gname = (g[1] if g.shape[1] >= 2 else g[0]).astype(str).values
    # 目标基因 → 行号集合（1-based），重复名取首个
    wanted = set(wanted or [])
    keep_rows = {}
    for i, nm in enumerate(gname):
        if nm in wanted and nm not in keep_rows:
            keep_rows[nm] = i + 1  # mtx 是 1-based
    kept = list(keep_rows.keys())
    row2col = {r: k for k, r in enumerate(keep_rows.values())}

    bcs = pd.read_csv(barcodes_f, sep='\t', header=None)[0].astype(str).values

    # 逐行流式解析：格式 'row col val'
    rows, cols, vals = [], [], []
    with gzip.open(mtx, 'rt') as fh:
        for line in fh:
            if line.startswith('%'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            r = int(parts[0])
            if r in row2col:
                rows.append(row2col[r])
                cols.append(int(parts[1]) - 1)
                vals.append(float(parts[2]))
    n_cells = len(bcs)
    n_keep = len(kept)
    M = sp.csr_matrix((vals, (cols, rows)), shape=(n_cells, n_keep))
    df = pd.DataFrame(M.toarray(), columns=kept)
    return df, bcs


def read_h5(hp):
    import anndata as ad
    a = ad.read_h5ad(str(hp)) if str(hp).endswith('.h5ad') else ad.read_10x_h5(str(hp))
    return a


# 统一组装成 pandas DataFrame（细胞 × marker gene），流式读取只保留 marker
frames = []
for sid, spec in samples.items():
    try:
        if 'mtx' in spec:
            df, bcs = read_10x_mtx(spec['mtx'], spec['genes'], spec['barcodes'], wanted=ALLGENES)
            df['barcode'] = bcs
            df['sample'] = sid
        elif 'h5' in spec:
            a = read_h5(spec['h5'])
            keep = [g for g in ALLGENES if g in set(a.var_names)]
            sub = a[:, keep].to_df()
            sub['barcode'] = a.obs_names.values
            sub['sample'] = sid
            df = sub
        else:
            df = pd.read_csv(spec['csv'], sep=None, engine='python', index_col=0)
            df = df.T
            keep = [g for g in ALLGENES if g in df.columns]
            df = df[keep]
            df['barcode'] = df.index
            df['sample'] = sid
        frames.append(df)
        log(f'   [ok] {sid}: {df.shape}')
    except Exception as e:
        log(f'   [fail] {sid}: {type(e).__name__}: {e}')

if not frames:
    Path(OUT / 'S20_inventory.txt').write_text('\n'.join(LOG), encoding='utf-8')
    raise SystemExit('所有样本加载失败')

expr = pd.concat(frames, ignore_index=True)
log(f'\n[expr] 合并后: {expr.shape}')

# barcode 构造键：matrix 里是裸 barcode（AAAC...-1），metadata 是 <SAMPLE>_<barcode>
# 样本键 = 'PWH26D0' 等，正好是 metadata barcode 的前缀
def make_key(row):
    b = str(row['barcode'])
    s = str(row['sample'])
    # 若已含前缀则直接用；否则拼接
    if b.startswith(s + '_'):
        return b
    return f'{s}_{b}'

expr['bc_key'] = expr.apply(make_key, axis=1)
meta['bc_key'] = meta['barcode'].astype(str)

n_before = len(expr)
merged = expr.merge(meta[['bc_key', 'orig.ident', 'Patient', 'Condition',
                          'newMainCellTypes', 'nCount_RNA']],
                    on='bc_key', how='inner', suffixes=('', '_meta'))
log(f'[merge] {n_before:,} × metadata → {len(merged):,} 匹配')
RESULT['merge'] = {'before': int(n_before), 'after': int(len(merged))}

if len(merged) == 0:
    # 兜底：去掉 metadata 前缀后再匹配
    log('  [retry] 前缀拼接失败，尝试 barcode 后缀匹配...')
    expr['bc_key2'] = expr['barcode'].astype(str)
    meta['bc_key2'] = meta['barcode'].astype(str).str.split('_', n=1).str[-1]
    merged = expr.merge(meta[['bc_key2', 'Patient', 'Condition', 'newMainCellTypes', 'nCount_RNA']],
                        left_on='bc_key2', right_on='bc_key2', how='inner')
    log(f'  [merge-retry] {len(merged):,}')
    RESULT['merge']['after_retry'] = int(len(merged))

merged = merged[merged['newMainCellTypes'] == 'Fibroblast'].copy()
log(f'[expr] 成纤维细胞: {len(merged):,}')
RESULT['merge']['n_fib_expr'] = int(len(merged))

if len(merged) < 200:
    log('⚠️ 成纤维细胞过少，可能 barcode 未匹配')
    log('样本示例:', expr['barcode'].head().tolist())
    log('meta 示例:', meta['barcode'].head().tolist())
    Path(OUT / 'S20_inventory.txt').write_text('\n'.join(LOG), encoding='utf-8')
    raise SystemExit('barcode 匹配需人工检查')


# =========================================================
# 4. 打分（score = 平均 log1p CPM，去批次=按供体 z-score）
# =========================================================
merged['logumi'] = np.log1p(merged['nCount_RNA'].astype(float))
logu = merged['logumi'].values
logu_c = logu - logu.mean()

def score_program(df, genes):
    gs = [g for g in genes if g in df.columns]
    if not gs:
        return np.full(len(df), np.nan), []
    sub = df[gs].astype(float)
    # 每个基因在全体细胞内 z-score（消除基因表达量差异），再平均
    z = (sub - sub.mean(0)) / (sub.std(0) + 1e-9)
    return z.mean(1).values, gs

for prog, genes in M.items():
    s, used = score_program(merged, genes)
    merged[f'sc_{prog}'] = s
    if np.isnan(s).all():
        log(f'   ⚠️ {prog}: 0/{len(genes)} 基因命中')
    else:
        miss = [g for g in genes if g not in merged.columns]
        log(f'   {prog}: {len(used)}/{len(genes)} 命中' + (f'  缺 {miss}' if miss else ''))

# 深度残差校正（对每个程序的分数）
SC = [f'sc_{p}' for p in M]
for c in SC:
    v = merged[c].values.copy()
    ok = ~np.isnan(v)
    if ok.sum() < 10:
        merged[c + '_res'] = v
        continue
    b = np.polyfit(logu_c[ok], v[ok], 1)
    merged[c + '_res'] = v - (b[0] * logu_c + b[1])

RESULT['marker_hits'] = {
    p: {'used': [g for g in genes if g in merged.columns],
        'missing': [g for g in genes if g not in merged.columns]}
    for p, genes in M.items()
}


# =========================================================
# 5. 时间轨迹（细胞级 + 供体级）
# =========================================================
log('\n' + '=' * 80)
log('【轨迹 1】各时间点的程序得分（细胞级均值，深度残差校正后）')
log('=' * 80)

tp_rows = []
for prog in M:
    row = {'program': prog}
    for t in ORDER:
        v = merged.loc[merged['Condition'] == t, f'sc_{prog}_res']
        row[LABEL[t]] = float(np.nanmean(v))
    tp_rows.append(row)
tp_df = pd.DataFrame(tp_rows).set_index('program')
log(tp_df.round(4).to_string())
tp_df.to_csv(OUT / 'S20_trajectory_scores.csv', encoding='utf-8-sig')

log('\n【轨迹 2】供体级（每个供体在每时间点的均值，只保留有成纤维细胞的时间点）')
donor_rows = []
for (d, t), sub in merged.groupby(['Patient', 'Condition']):
    if len(sub) < 20:
        continue
    r = {'donor': d, 'timepoint': t, 'n': len(sub), 'logumi': float(sub['logumi'].mean())}
    for prog in M:
        r[prog] = float(np.nanmean(sub[f'sc_{prog}_res']))
    donor_rows.append(r)
donor_df = pd.DataFrame(donor_rows)
donor_df['tp'] = donor_df['timepoint'].map(LABEL)
donor_df.to_csv(OUT / 'S20_donor_trajectory.csv', index=False, encoding='utf-8-sig')
log(donor_df.round(4).to_string())


# =========================================================
# 6. 趋势检验：Spearman rho（时间 vs 分数）+ 供体级符号一致性
# =========================================================
from scipy.stats import spearmanr

log('\n' + '=' * 80)
log('【趋势检验】时间序数 (0,1,2,3) 与程序得分的 Spearman 相关')
log('=' * 80)
tmap = {t: i for i, t in enumerate(ORDER)}
merged['t_ord'] = merged['Condition'].map(tmap)

trend = {}
for prog in M:
    v = merged[f'sc_{prog}_res'].values
    ok = ~np.isnan(v)
    rho, p = spearmanr(merged.loc[ok, 't_ord'], v[ok])
    # 供体内符号一致性
    donor_rhos = []
    for d in merged['Patient'].unique():
        sub = merged[merged['Patient'] == d]
        vv = sub[f'sc_{prog}_res'].values
        oo = ~np.isnan(vv)
        if oo.sum() > 50 and sub.loc[oo, 't_ord'].nunique() >= 3:
            rr, _ = spearmanr(sub.loc[oo, 't_ord'], vv[oo])
            donor_rhos.append(rr)
    n_pos = sum(1 for r in donor_rhos if r > 0)
    n_neg = sum(1 for r in donor_rhos if r < 0)
    trend[prog] = {'rho_cell': float(rho), 'p_cell': float(p),
                   'donor_rhos': [float(x) for x in donor_rhos],
                   'n_pos': n_pos, 'n_neg': n_neg}
    log(f'  {prog:6s} rho={rho:+.4f}  p={p:.3g}   供体: +{n_pos} / -{n_neg}   {["%.3f" % x for x in donor_rhos]}')
RESULT['trend'] = trend


# =========================================================
# 7. 置换检验（供体内时间标签重排，控制个体效应）
# =========================================================
log('\n' + '=' * 80)
log('【置换检验】供体内时间标签重排（严格零模型：保持供体结构）')
log('=' * 80)

def within_donor_perm_trend(df, col, nperm=2000):
    obs_rho, _ = spearmanr(df.loc[~df[col].isna(), 't_ord'], df.loc[~df[col].isna(), col])
    cnt = 0
    donors = df['Patient'].unique()
    for _ in range(nperm):
        d2 = df.copy()
        for d in donors:
            idx = d2.index[d2['Patient'] == d]
            d2.loc[idx, 't_ord'] = RNG.permutation(d2.loc[idx, 't_ord'].values)
        v = d2[col].values; ok = ~np.isnan(v)
        r, _ = spearmanr(d2.loc[ok, 't_ord'], v[ok])
        if r >= obs_rho:
            cnt += 1
    return obs_rho, (cnt + 1) / (nperm + 1)

perm_res = {}
for prog in M:
    rho, p = within_donor_perm_trend(merged, f'sc_{prog}_res')
    perm_res[prog] = {'rho': float(rho), 'p_perm': float(p)}
    log(f'  {prog:6s} rho={rho:+.4f}  p_perm={p:.4f}')
RESULT['perm_trend'] = perm_res


# =========================================================
# 8. W6: 修复方向向量 —— 正常愈合 (F6,F7) 轨迹 → 投影 DFU 位移
# =========================================================
log('\n' + '=' * 80)
log('【W6】修复方向向量：正常愈合轨迹 vs DFU Healer-NonHealer 位移')
log('=' * 80)

# 正常愈合方向 = (D30 - D0) in (F6res, F7res) 标准化空间
d0 = merged[merged['Condition'] == 'Skin']
d30 = merged[merged['Condition'] == 'Wound30']
vec_D0 = np.array([np.nanmean(d0['sc_F6_res']), np.nanmean(d0['sc_F7_res'])])
vec_D30 = np.array([np.nanmean(d30['sc_F6_res']), np.nanmean(d30['sc_F7_res'])])
v_normal = vec_D30 - vec_D0
v_normal_unit = v_normal / (np.linalg.norm(v_normal) + 1e-12)
log(f'  正常愈合 F6/F7 位移 (D30-D0): {v_normal.round(4)}')
log(f'  单位方向: {v_normal_unit.round(4)}')

# 读入 DFU (GSE165816) 的患者级分数（S14 产物）
dfu_csv = OUT / 'S14_state_space_patients.csv'
w6 = {'v_normal': v_normal.tolist(), 'v_normal_unit': v_normal_unit.tolist()}
if dfu_csv.exists():
    sd = pd.read_csv(dfu_csv)
    log(f'  DFU 患者级文件列: {list(sd.columns)}')
    # 兼容不同列名
    f6c = [c for c in sd.columns if 'F6' in c][0]
    f7c = [c for c in sd.columns if 'F7' in c][0]
    grp = [c for c in sd.columns if 'group' in c.lower() or 'Group' in c][0]
    h = sd[sd[grp].astype(str).str.contains('Heal', case=False) & ~sd[grp].astype(str).str.contains('Non', case=False)]
    nh = sd[sd[grp].astype(str).str.contains('Non', case=False)]
    if len(h) and len(nh):
        v_dfu = np.array([h[f6c].mean() - nh[f6c].mean(), h[f7c].mean() - nh[f7c].mean()])
        v_dfu_unit = v_dfu / (np.linalg.norm(v_dfu) + 1e-12)
        cos = float(np.dot(v_dfu_unit, v_normal_unit))
        log(f'  DFU Healer-NonHealer 位移: {v_dfu.round(4)}  (n={len(h)} vs {len(nh)})')
        log(f'  🔴 余弦相似度 cos(正常愈合方向, DFU 愈合位移) = {cos:+.4f}')
        w6.update({'v_dfu': v_dfu.tolist(), 'v_dfu_unit': v_dfu_unit.tolist(),
                   'cos': cos, 'n_healer': int(len(h)), 'n_nonhealer': int(len(nh))})
else:
    log(f'  ⚠️ 未找到 {dfu_csv}')
RESULT['w6'] = w6


# =========================================================
# 9. 输出
# =========================================================
with open(OUT / 'S20_trajectory.json', 'w', encoding='utf-8') as f:
    json.dump(RESULT, f, indent=2, ensure_ascii=False)
Path(OUT / 'S20_run.log').write_text('\n'.join(LOG), encoding='utf-8')
log('\n[saved] S20_trajectory.json / S20_run.log / S20_*.csv')
log('✅ S20 完成')
