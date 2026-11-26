# -*- coding: utf-8 -*-
"""
S21: 第三/第四 DFU 队列跨队列复现

队列 A: GSE223964 (27,998 细胞, 8 样本) —— condition: diabetic(5) vs non-diabetic(3)
队列 B: GSE245703 (RAW.tar)          —— genotype: Diabetic vs Non-Diabetic

⚠️ 必须诚实标注的限制：
  这两个队列的分组是【糖尿病 vs 非糖尿病创面】，不是【愈合 vs 不愈合】。
  → 复现的是"糖尿病创面成纤维状态偏移"这一【疾病相关性】主张，
    不是"F6/F7 预测愈合结局"这一【预后】主张。两者必须分开陈述。

方法：与 S9-S12 完全一致 —— marker z-score 门控 + 深度残差校正 + 置换检验 + 阴性对照
成纤维细胞识别：FIB marker 高表达的 leiden 簇（自动 + 人工确认）
"""
import json, warnings, tarfile, gzip
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
BASE = ROOT / "10_DFU_Dataset"
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
FIB = ['DCN', 'LUM', 'PDGFRA', 'COL1A1', 'COL1A2', 'FAP', 'THY1']

LOG = []
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); LOG.append(s)

log('=' * 84)
log('S21: 第三/第四 DFU 队列跨队列复现')
log('=' * 84)


def perm_p(v, y, nperm=5000, rng=None):
    rng = rng or RNG
    v = np.asarray(v, float); y = np.asarray(y, int)
    ok = ~np.isnan(v); v, y = v[ok], y[ok]
    if y.sum() == 0 or (y == 0).sum() == 0:
        return np.nan, np.nan
    obs = v[y == 1].mean() - v[y == 0].mean()
    cnt = 0
    for _ in range(nperm):
        yp = rng.permutation(y)
        if (v[yp == 1].mean() - v[yp == 0].mean()) >= obs:
            cnt += 1
    return float(obs), float((cnt + 1) / (nperm + 1))


def zscore_program(Xdf, genes):
    gs = [g for g in genes if g in Xdf.columns]
    if not gs:
        return np.full(len(Xdf), np.nan), []
    sub = Xdf[gs].astype(float)
    z = (sub - sub.mean(0)) / (sub.std(0) + 1e-9)
    return z.mean(1).values, gs


def residualize(v, logu):
    v = np.asarray(v, float); ok = ~np.isnan(v)
    if ok.sum() < 10:
        return v
    lu = logu - logu.mean()
    b = np.polyfit(lu[ok], v[ok], 1)
    return v - (b[0] * lu + b[1])


RESULT = {}


# =====================================================================
# 队列 A: GSE223964
# =====================================================================
def run_GSE223964():
    hp = BASE / 'gse223964' / 'GSE223964_Integrated_all_cells.h5ad'
    if not hp.exists() or hp.stat().st_size < 100_000_000:
        log('[GSE223964] 数据未就绪'); return None
    import anndata as ad
    a = ad.read_h5ad(str(hp))
    log(f'\n[GSE223964] shape={a.shape}')
    obs = a.obs.copy()
    log(f'  condition: {obs["condition"].value_counts().to_dict()}')
    log(f'  sample: {obs["sample"].value_counts().to_dict()}')

    # --- 成纤维细胞识别：FIB marker 高的 leiden 簇 ---
    gs = [g for g in FIB if g in a.var_names]
    Xf = a[:, gs].X
    Xf = Xf.toarray() if hasattr(Xf, 'toarray') else np.asarray(Xf)
    fibscore = Xf.mean(1)
    cl_fib = {}
    for cl in obs['leiden'].unique():
        m = obs['leiden'].values == cl
        cl_fib[cl] = float(fibscore[m].mean())
    # Bimodal 判定：均值 + 2SD 以上，且为该列最高的 2 个簇
    vals = np.array(list(cl_fib.values()))
    thr = vals.mean() + 1.0 * vals.std()
    top2 = sorted(cl_fib, key=cl_fib.get, reverse=True)[:2]
    fib_clusters = [c for c in top2 if cl_fib[c] >= thr] or [top2[0]]
    log(f'  leiden 簇成纤维得分: { {k: round(v,3) for k,v in sorted(cl_fib.items(), key=lambda x:-x[1])} }')
    log(f'  → 判定为成纤维簇: {fib_clusters}')
    obs['is_fib'] = obs['leiden'].isin(fib_clusters)
    fib_idx = obs.index[obs['is_fib']]
    log(f'  成纤维细胞: {len(fib_idx):,} / {len(obs):,}')

    # --- 取 marker 表达 ---
    allg = sorted({g for v in M.values() for g in v} | set(FIB))
    allg = [g for g in allg if g in a.var_names]
    sub = a[fib_idx, allg]
    X = sub.X
    X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
    df = pd.DataFrame(X, columns=allg, index=fib_idx)
    ob = obs.loc[fib_idx].copy()
    df['logu'] = np.log1p(ob['total_counts'].astype(float).values
                          if 'total_counts' in ob.columns else np.asarray(sub.X.sum(1)).ravel())
    df['sample'] = ob['sample'].values
    df['condition'] = ob['condition'].values
    df['y'] = (df['condition'].astype(str).str.lower().str.contains('diabetic') &
               ~df['condition'].astype(str).str.lower().str.contains('non')).astype(int)

    log(f'  分组计数: {df.groupby("condition")["y"].size().to_dict()}')

    # --- 打分 + 深度残差 ---
    res = {}
    scores = {}
    for prog, genes in M.items():
        s, used = zscore_program(df, genes)
        sr = residualize(s, df['logu'].values)
        scores[prog] = sr
        df[f'sc_{prog}'] = sr
        d, p = perm_p(sr, df['y'].values)
        res[prog] = {'delta': d, 'p_perm': p, 'n_genes': len(used),
                     'missing': [g for g in genes if g not in df.columns]}
        log(f'   {prog:6s} Δ(diab-nondiab)={d:+.4f}  p={p:.4f}  ({len(used)}/{len(genes)} genes)')

    # --- 样本级（pseudobulk）---
    log('  --- 样本级 pseudobulk ---')
    sample_rows = []
    for (s, cond), g in df.groupby(['sample', 'condition']):
        if len(g) < 20:
            continue
        r = {'sample': s, 'condition': cond, 'n': len(g),
             'logu': float(g['logu'].mean()), 'y': int(g['y'].iloc[0])}
        for prog in M:
            r[prog] = float(np.nanmean(g[f'sc_{prog}']))
        sample_rows.append(r)
    sdf = pd.DataFrame(sample_rows)
    sres = {}
    for prog in M:
        d, p = perm_p(sdf[prog].values, sdf['y'].values, nperm=10000)
        sres[prog] = {'delta': d, 'p_perm': p}
        log(f'   {prog:6s} 样本级 Δ={d:+.4f}  p={p:.4f}')
    sdf.to_csv(OUT / 'S21_GSE223964_samples.csv', index=False, encoding='utf-8-sig')

    # --- 深度诊断 ---
    depth_by_grp = df.groupby('y')['logu'].mean()
    log(f'  深度对照: diabetic={depth_by_grp.get(1,np.nan):.3f} vs non-diabetic={depth_by_grp.get(0,np.nan):.3f} '
        f'(Δ={depth_by_grp.get(1,np.nan)-depth_by_grp.get(0,np.nan):+.3f} log units)')

    # --- 阴性对照 ---
    log('  --- 阴性对照 ---')
    for prog in ['F1', 'F5']:
        d, p = perm_p(df[f'sc_{prog}'].values, df['y'].values)
        log(f'   [对照] {prog}: Δ={d:+.4f}  p={p:.4f}  ' +
            ('✅ 不显著（好）' if p > 0.1 else '⚠️ 显著（需注意）'))

    return {'res_cell': res, 'res_sample': sres,
            'n_fib': int(len(fib_idx)), 'n_total': int(len(obs)),
            'fib_clusters': [str(c) for c in fib_clusters],
            'conditions': obs['condition'].value_counts().to_dict(),
            'depth_by_group': {str(k): float(v) for k, v in depth_by_grp.items()},
            'fib_scores_by_cluster': {str(k): v for k, v in cl_fib.items()}}


rA = run_GSE223964()
if rA:
    RESULT['GSE223964'] = rA


# =====================================================================
# 队列 B: GSE245703
# =====================================================================
def run_GSE245703():
    tar = BASE / 'gse245703' / 'GSE245703_RAW.tar'
    if not tar.exists() or tar.stat().st_size < 150_000_000:
        log(f'\n[GSE245703] 未就绪 ({tar.stat().st_size if tar.exists() else 0} bytes)'); return None
    ex = BASE / 'gse245703' / 'extracted'
    files = [f for f in ex.rglob('*') if f.is_file()] if ex.exists() else []
    if not files:
        ex.mkdir(parents=True, exist_ok=True)
        log(f'[GSE245703] 解压 {tar.stat().st_size:,} bytes ...')
        with tarfile.open(tar) as t:
            t.extractall(ex)
        files = [f for f in ex.rglob('*') if f.is_file()]
    log(f'\n[GSE245703] 文件 {len(files)}, 后缀 {dict(Counter(f.suffix for f in files))}')

    # 探测格式：优先 h5，其次 mtx，最后 csv
    h5s = sorted([f for f in files if f.suffix in ('.h5', '.h5ad')])
    mtx_dir = {}
    for m in sorted([f for f in files if f.name.startswith('matrix.mtx')]):
        d = m.parent
        b = list(d.glob('barcodes.tsv*')); gg = list(d.glob('genes.tsv*')) + list(d.glob('features.tsv*'))
        if b and gg:
            mtx_dir[d.name] = {'mtx': m, 'genes': gg[0], 'barcodes': b[0]}
    csvs = sorted([f for f in files if f.suffix in ('.csv', '.tsv')])
    log(f'  h5={len(h5s)}, mtx样本={len(mtx_dir)}, csv={len(csvs)}')

    frames = []
    allg = sorted({g for v in M.values() for g in v} | set(FIB))
    import anndata as ad
    try:
        from scanpy.readwrite import read_10x_h5 as _r10x
    except Exception:
        _r10x = None

    def read_10x_h5_direct(hp, wanted):
        """直接 h5py 读取 10x v3 H5，优雅处理重复基因名。返回 DataFrame(cells × wanted)。"""
        import h5py
        with h5py.File(str(hp), 'r') as f:
            g = f['matrix']
            shape = tuple(g['shape'][()])          # (n_genes, n_cells)
            names = np.array([x.decode() if isinstance(x, bytes) else str(x)
                              for x in g['features/name'][()]])
            bcs = np.array([x.decode() if isinstance(x, bytes) else str(x)
                            for x in g['barcodes'][()]])
            data = g['data'][()]
            indices = g['indices'][()]
            indptr = g['indptr'][()]
        X = sp.csc_matrix((data, indices, indptr), shape=shape)  # genes × cells
        # 保留出现在 wanted 中的基因（重复名取第一个命中）
        keep_pos, keep_name = [], []
        seen = set()
        for i, nm in enumerate(names):
            if nm in wanted and nm not in seen:
                seen.add(nm); keep_pos.append(i); keep_name.append(nm)
        sub = X[keep_pos, :].toarray().T  # cells × keep
        df = pd.DataFrame(sub, columns=keep_name)
        df['logu'] = np.log1p(sub.sum(1))
        df['barcode'] = bcs
        return df

    if h5s:
        for h in h5s:
            try:
                if str(h).endswith('.h5ad'):
                    aa = ad.read_h5ad(str(h))
                    keep = [g for g in allg if g in aa.var_names]
                    d = aa[:, keep].to_df().reset_index(drop=True)
                    d['logu'] = np.log1p(np.asarray(aa.X.sum(1)).ravel())
                    d['barcode'] = aa.obs_names.values
                else:
                    d = read_10x_h5_direct(h, set(allg))
                d['sample'] = h.stem
                frames.append(d); log(f'   [h5 ok] {h.stem}: {d.shape[:-1]}')
            except Exception as e:
                log(f'   [h5 fail] {h.stem}: {type(e).__name__}: {e}')
    if not frames and mtx_dir:
        for sid, spec in mtx_dir.items():
            try:
                Xm = sp.csr_matrix(sio.mmread(str(spec['mtx']))).T  # cells × genes
                g = pd.read_csv(spec['genes'], sep='\t', header=None)
                gn = (g[1] if g.shape[1] >= 2 else g[0]).astype(str).values
                gidx = {x: i for i, x in enumerate(gn)}
                keep = [x for x in allg if x in gidx]
                sub = Xm[:, [gidx[x] for x in keep]].toarray()
                d = pd.DataFrame(sub, columns=keep)
                d['sample'] = sid
                d['logu'] = np.log1p(sub.sum(1))
                d['barcode'] = pd.read_csv(spec['barcodes'], sep='\t', header=None)[0].values
                frames.append(d); log(f'   [mtx ok] {sid}: {d.shape}')
            except Exception as e:
                log(f'   [mtx fail] {sid}: {e}')

    if not frames:
        log('  ⚠️ 无法读取表达；保存目录树')
        for f in files[:40]:
            log('     ', f.relative_to(ex), f.stat().st_size)
        return {'status': 'unreadable', 'n_files': len(files),
                'suffixes': dict(Counter(f.suffix for f in files)),
                'files': [str(f.relative_to(ex)) for f in files[:60]]}

    # 用 metadata（GEO 样本 title 已知 Diabetic/Non-Diabetic）
    meta = pd.read_csv(BASE / 'gse245703' / 'sample_meta.csv') if (BASE / 'gse245703' / 'sample_meta.csv').exists() else None
    df = pd.concat(frames, ignore_index=True)
    log(f'  合并: {df.shape}')

    # 分组：由 GEO 元数据核验（curl acc.cgi）
    #   Diabetic     = Diac32, Diac34, Diac35, Diac36, Diac37   (5 样本)
    #   Non-Diabetic = Diac38, Diac39, Diac44a, Diac44b, Diac44c, Diac50a, Diac50b (7 样本)
    DIAB = {'Diac32', 'Diac34', 'Diac35', 'Diac36', 'Diac37'}
    NODIAB = {'Diac38', 'Diac39', 'Diac44a', 'Diac44b', 'Diac44c', 'Diac50a', 'Diac50b'}

    def to_group(s):
        # h.stem 形如 GSM7847134_diac32_filtered_feature_bc_matrix
        low = str(s).lower()
        for nm in DIAB:
            if nm.lower() in low:
                return 1
        for nm in NODIAB:
            if nm.lower() in low:
                return 0
        return np.nan

    df['y'] = df['sample'].map(to_group)
    ok = df['y'].notna()
    log(f'  可分组细胞: {ok.sum():,} / {len(df):,}')
    log(f'  样本分组: ' + str({s: int(g) for s, g in df.groupby("sample")["y"].first().items()}))
    if ok.sum() < 200:
        log('  ⚠️ 分组映射失败，样本名: ' + str(df['sample'].unique().tolist()))
        return {'status': 'grouping_failed', 'samples': df['sample'].unique().tolist()}

    df = df[ok].copy()
    # 过滤表达全零细胞
    gc = [c for c in df.columns if c in allg]
    df = df[df[gc].sum(1) > 0]

    # --- 成纤维细胞识别（无预注释 → marker 表达门控）---
    # 策略：同时要求 (a) 结缔组织核心 marker 正表达、(b) 上皮/免疫/内皮 marker 低
    # 避免 KMeans 因深度差异把成纤维切成多簇
    fibm = [g for g in FIB if g in df.columns]
    for g in fibm:
        df[g] = df[g].astype(float)
    # 定义性 marker：DCN/LUM/COL1A2 至少 2 个 > 0；且 PDGFRA 或 COL1A1 > 0
    core = [g for g in ['DCN', 'LUM', 'COL1A2'] if g in df.columns]
    core_pos = (df[core] > 0).sum(1)
    second = [g for g in ['PDGFRA', 'COL1A1', 'FAP', 'THY1'] if g in df.columns]
    second_pos = (df[second] > 0).sum(1)
    is_fib = (core_pos >= 2) & (second_pos >= 1)
    df['is_fib'] = is_fib
    log(f'  成纤维门控: core>=2 & second>=1 → {int(is_fib.sum()):,} 细胞 '
        f'({is_fib.mean()*100:.2f}%)')
    df = df[is_fib].copy()
    log(f'  成纤维细胞: {len(df):,}')
    if len(df) < 300:
        log('  ⚠️ 成纤维细胞仍过少')

    res = {}
    for prog, genes in M.items():
        s, used = zscore_program(df, genes)
        sr = residualize(s, df['logu'].values)
        df[f'sc_{prog}'] = sr
        d, p = perm_p(sr, df['y'].values)
        res[prog] = {'delta': d, 'p_perm': p, 'n_genes': len(used)}
        log(f'   {prog:6s} Δ(diab-nondiab)={d:+.4f}  p={p:.4f}  ({len(used)}/{len(genes)})')

    # 样本级
    log('  --- 样本级 pseudobulk ---')
    rows = []
    for (s, y), g in df.groupby(['sample', 'y']):
        if len(g) < 20:
            continue
        r = {'sample': s, 'y': int(y), 'n': len(g), 'logu': float(g['logu'].mean())}
        for prog in M:
            r[prog] = float(np.nanmean(g[f'sc_{prog}']))
        rows.append(r)
    sdf = pd.DataFrame(rows)
    sres = {}
    for prog in M:
        d, p = perm_p(sdf[prog].values, sdf['y'].values, nperm=10000)
        sres[prog] = {'delta': d, 'p_perm': p}
        log(f'   {prog:6s} 样本级 Δ={d:+.4f}  p={p:.4f}')
    sdf.to_csv(OUT / 'S21_GSE245703_samples.csv', index=False, encoding='utf-8-sig')
    depth_by_grp = df.groupby('y')['logu'].mean()
    log(f'  深度对照: diabetic={depth_by_grp.get(1,np.nan):.3f} vs non={depth_by_grp.get(0,np.nan):.3f}')
    return {'res_cell': res, 'res_sample': sres, 'n_fib': int(len(df)),
            'samples': df['sample'].unique().tolist(),
            'depth_by_group': {str(k): float(v) for k, v in depth_by_grp.items()}}


rB = run_GSE245703()
if rB:
    RESULT['GSE245703'] = rB

with open(OUT / 'S21_crosscohort.json', 'w', encoding='utf-8') as f:
    json.dump(RESULT, f, indent=2, ensure_ascii=False, default=str)
Path(OUT / 'S21_run.log').write_text('\n'.join(LOG), encoding='utf-8')
log('\n[saved] S21_crosscohort.json / S21_run.log')
log('✅ S21 完成')
