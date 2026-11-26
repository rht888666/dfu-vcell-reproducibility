# -*- coding: utf-8 -*-
"""
S2: 从 Southard mean_pop 提取「TF 扰动谱矩阵」—— 虚拟细胞的扰动操作手册

输入: fibroblast_CRISPRa_mean_pop.h5ad (10916 guides x 4914 genes)
  - X: guide 对各基因的效应系数 β
  - obs.target_gene: TF 名
  - obs.active / masked_active: 是否判定为有效激活
  - obs.cluster: 表达程序编号

输出:
  A. 每个 TF 的平均扰动向量（在 active guides 上聚合）→ (n_TF, n_genes)
  B. 每个 TF 的扰动强度（norm）
  C. myo/universal/inflammatory 程序（cluster 4/31/32）的 top TF
  D. 与 DFU 相关的 TF 子集标记

严格性:
  - 只用 masked_active=True 的 guide（论文的严格版活性判定）
  - 若某 TF 无 active guide → 排除并记录
  - 报告每个 TF 的 guide 数（可靠性指标）
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

warnings.filterwarnings('ignore')
OUT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿/17_DFU_VirtualCell_10plus")
D = OUT / "data_southard"

print('[load] Southard mean_pop')
a = ad.read_h5ad(D / 'fibroblast_CRISPRa_mean_pop.h5ad')
print(f'  {a.shape}   guides x genes')
print(f'  target_gene: {a.obs["target_gene"].nunique()} unique')
print(f'  active: {int(a.obs["active"].sum())}   masked_active: {int(a.obs["masked_active"].sum())}')

X = np.asarray(a.X)
genes = a.var_names.astype(str)
obs = a.obs.copy()
print(f'  X shape {X.shape}  finite={np.isfinite(X).all()}')

meta = {'n_guides': int(X.shape[0]), 'n_genes': int(X.shape[1]),
        'n_tf': int(obs['target_gene'].nunique()),
        'n_active': int(obs['active'].sum()),
        'n_masked_active': int(obs['masked_active'].sum())}

# ---------------- A. 按 TF 聚合 ----------------
print('\n[1] 按 TF 聚合扰动向量（masked_active 优先）')
use = obs['masked_active'].astype(bool).values
if use.sum() < 100:
    print('  masked_active 太少，回退到 active')
    use = obs['active'].astype(bool).values
print(f'  使用 guide: {use.sum()} / {len(use)}')

Xn = X[use]
tg = obs['target_gene'].values[use]
cl = obs['cluster'].values[use]

tf_names = sorted(set(tg))
rows, tfs, nguide = [], [], []
for t in tf_names:
    mm = tg == t
    rows.append(Xn[mm].mean(0))
    tfs.append(t)
    nguide.append(int(mm.sum()))
M = np.array(rows)          # (n_TF, n_genes)
tfs = np.array(tfs)
nguide = np.array(nguide)
print(f'  扰动矩阵 M: {M.shape}  (TF x genes)')
print(f'  guides per TF: min={nguide.min()} median={np.median(nguide):.0f} max={nguide.max()}')

# 扰动强度
strength = np.linalg.norm(M, axis=1)
print(f'  扰动强度 ||β||: min={strength.min():.3f} median={np.median(strength):.3f} '
      f'max={strength.max():.3f}')
print(f'  Top 15 最强扰动 TF: {", ".join(tfs[np.argsort(-strength)[:15]])}')

# ---------------- B. myo/universal/inflammatory 程序 ----------------
print('\n[2] 表达程序（cluster）与 fibroblast state')
# 论文: 程序 4=myofibroblast, 31=universal, 32=inflammatory
PROG = {4: 'myofibroblast', 31: 'universal', 32: 'inflammatory'}
# cluster 是 float，可能经过 masked/stricter 重编号 → 检查
for cn, cnt in obs['cluster'].value_counts().head(12).items():
    print(f'    cluster {cn}: {cnt}')

prog_hits = {}
for cl_id, nm in PROG.items():
    if (obs['cluster'].values == cl_id).any():
        sub = obs[obs['cluster'].values == cl_id]
        prog_hits[nm] = {'n_guides': int(len(sub)),
                         'tfs': sorted(set(sub['target_gene']))[:40]}
        print(f'  程序 {cl_id} ({nm}): {len(sub)} guides, '
              f'{sub["target_gene"].nunique()} TFs')
        print(f'    TFs: {", ".join(sorted(set(sub["target_gene"]))[:30])}')
    else:
        print(f'  程序 {cl_id} ({nm}): 未找到')

# ---------------- C. DFU 相关 TF ----------------
print('\n[3] DFU 相关 TF 是否在矩阵中')
DFU_TF = {
    'hypoxia': ['HIF1A', 'EPAS1', 'ARNT', 'VHL'],
    'fibrosis': ['SMAD2', 'SMAD3', 'SMAD4', 'SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2'],
    'diabetes': ['FOXO1', 'FOXO3', 'PPARG', 'TCF7L2', 'HNF1A', 'HNF4A'],
    'southard': ['KLF2', 'KLF4', 'KLF5', 'PLAG1', 'PLAGL1', 'PLAGL2'],
    'inflammation': ['NFKB1', 'NFKB2', 'RELA', 'STAT1', 'STAT3', 'STAT6', 'CEBPB'],
    'master': ['MYC', 'TP53', 'JUN', 'FOS', 'ATF3', 'EGR1', 'SP1', 'YY1', 'NFE2L2'],
    'repair': ['GRHL2', 'GRHL3', 'TP63'],
}
dfu_cov = {}
for g, lst in DFU_TF.items():
    hit = [t for t in lst if t in set(tfs)]
    miss = [t for t in lst if t not in set(tfs)]
    dfu_cov[g] = {'hit': hit, 'miss': miss}
    print(f'  {g:<14} {len(hit)}/{len(lst)}  HIT={hit}'
          + (f'  MISS={miss}' if miss else ''))

# ---------------- D. 导出 ----------------
print('\n[4] 导出')
np.savez_compressed(D / 'tf_perturbation_matrix.npz',
                    M=M.astype(np.float32),
                    tfs=tfs,
                    genes=genes,
                    strength=strength.astype(np.float32),
                    n_guides=nguide.astype(np.int32))
print(f'  → tf_perturbation_matrix.npz  (M {M.shape})')

# 保存 gene-level 的 per-TF 效应表（可读）
np.save(D / 'tf_names.npy', tfs)
np.save(D / 'gene_names.npy', genes)

meta['prog_hits'] = {k: {'n_guides': v['n_guides'], 'n_tf': len(v['tfs'])}
                     for k, v in prog_hits.items()}
meta['prog_tfs'] = prog_hits
meta['dfu_coverage'] = dfu_cov
meta['top_strong_tfs'] = tfs[np.argsort(-strength)[:50]].tolist()
with open(OUT / 'S2_perturbation_extract.json', 'w', encoding='utf-8') as f:
    json.dump(meta, f, indent=2, ensure_ascii=False)
print('  → S2_perturbation_extract.json')
print(f'\n[done] {len(tfs)} TF x {len(genes)} genes')
