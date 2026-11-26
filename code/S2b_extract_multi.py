# -*- coding: utf-8 -*-
"""
S2b: 多版本 TF 扰动谱提取 —— 解决覆盖率问题

问题：masked_active 只有 661 guides（6%），只覆盖 213 TF，
      DFU 核心 TF（HIF1A/SMAD2/3/4/STAT3/NFKB1）全部缺失。

解决：提取三个版本的扰动矩阵
  V1 strict  : masked_active guides → 高置信度，但覆盖 213 TF
  V2 active  : active guides → 覆盖更多 TF
  V3 all     : 全部 guides（含未激活）→ 覆盖全部 1,837 TF
               （未激活 guide 的 β 接近 0，均值仍保留真实扰动信号）
  V4 strict_de: 按 de_genes 加权（效应大的 guide 权重高）

同时在 gene 空间也做：
  G_all  : 全部 4,914 基因
  G_ov   : 与 DFU 数据交集基因

关键：只有**与 DFU 数据基因空间有交集**的 TF 才能用于投影。
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

warnings.filterwarnings('ignore')
OUT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿/17_DFU_VirtualCell_10plus")
D = OUT / "data_southard"

a = ad.read_h5ad(D / 'fibroblast_CRISPRa_mean_pop.h5ad')
X = np.asarray(a.X).astype(np.float32)
genes = a.var_names.astype(str)
obs = a.obs.copy()
tg = obs['target_gene'].values
print(f'{a.shape}  targets={len(set(tg))}')

# DFU 基因空间
print('\n[0] 基因空间交集')
dfu = ad.read_h5ad(Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿/10_DFU_Dataset/dfu_outcome_dataset/GSE165816_sc.h5ad"),
                   backed='r')
dfu_genes = set(dfu.var_names.astype(str))
inter = [g for g in genes if g in dfu_genes]
print(f'  Southard genes: {len(genes)}   DFU genes: {len(dfu_genes)}')
print(f'  交集: {len(inter)}  ({100*len(inter)/len(genes):.1f}% of Southard)')
gi = {g: i for i, g in enumerate(genes)}
keep_idx = [gi[g] for g in inter]

meta = {'n_genes_southard': len(genes), 'n_genes_dfu': len(dfu_genes),
        'n_intersect': len(inter)}

def build(mask, tag, weight=None):
    tgm = tg[mask]
    Xm = X[mask]
    if weight is not None:
        w = weight[mask].astype(np.float32)
    tfs = sorted(set(tgm))
    rows, ns, nw = [], [], []
    for t in tfs:
        mm = tgm == t
        if weight is not None:
            ww = w[mm]
            if ww.sum() <= 0:
                rows.append(np.zeros(X.shape[1], dtype=np.float32)); ns.append(0); nw.append(0.0); continue
            rows.append((Xm[mm] * ww[:, None]).sum(0) / ww.sum())
            nw.append(float(ww.sum()))
        else:
            rows.append(Xm[mm].mean(0))
        ns.append(int(mm.sum()))
    M = np.array(rows, dtype=np.float32)
    print(f'  [{tag}] {M.shape}  TFs={len(tfs)}  '
          f'guides/TF median={np.median(ns):.0f}  '
          f'非零基因/TF median={int(np.median((np.abs(M)>1e-6).sum(1)))}')
    return M, np.array(tfs), np.array(ns), np.array(nw)

print('\n[1] 构建多版本扰动矩阵')
res = {}
M1, T1, N1, _ = build(obs['masked_active'].astype(bool).values, 'V1 strict')
M2, T2, N2, _ = build(obs['active'].astype(bool).values, 'V2 active')
M3, T3, N3, _ = build(np.ones(len(tg), bool), 'V3 all')
de = obs['de_genes'].astype(float).fillna(0).values
M4, T4, N4, W4 = build(obs['active'].astype(bool).values, 'V4 active+DE-weight', weight=np.maximum(de, 1))

# 只保留与 DFU 有交集的基因
for tag, M, T, N in [('V1', M1, T1, N1), ('V2', M2, T2, N2), ('V3', M3, T3, N3)]:
    Ms = M[:, keep_idx]
    np.savez_compressed(D / f'perturb_{tag}.npz', M=Ms, tfs=T, n_guides=N,
                        genes=np.array(inter))
    print(f'  保存 perturb_{tag}.npz  {Ms.shape}')
    res[tag] = {'n_tf': int(len(T)), 'n_genes': int(len(inter))}

np.savez_compressed(D / 'perturb_V4_dew.npz', M=M4[:, keep_idx], tfs=T4,
                    n_guides=N4, de_weight=W4, genes=np.array(inter))

# ---------------- 覆盖率报告 ----------------
print('\n[2] DFU 核心 TF 覆盖率（各版本）')
DFU_TF = ['HIF1A', 'EPAS1', 'ARNT', 'VHL', 'SMAD2', 'SMAD3', 'SMAD4', 'SNAI1', 'SNAI2',
          'TWIST1', 'ZEB1', 'ZEB2', 'FOXO1', 'FOXO3', 'PPARG', 'TCF7L2', 'HNF1A', 'HNF4A',
          'KLF2', 'KLF3', 'KLF4', 'KLF5', 'KLF6', 'KLF7', 'KLF9', 'KLF10', 'KLF11',
          'PLAG1', 'PLAGL1', 'PLAGL2', 'NFKB1', 'NFKB2', 'RELA', 'STAT1', 'STAT3',
          'STAT6', 'CEBPB', 'CEBPD', 'MYC', 'TP53', 'JUN', 'JUNB', 'FOS', 'ATF3',
          'EGR1', 'SP1', 'YY1', 'NFE2L2', 'GRHL2', 'GRHL3', 'TP63',
          'SMAD1', 'SMAD5', 'SMAD7', 'TGFBR1', 'TGFBR2', 'MMP1', 'MMP3', 'MMP9',
          'IL6', 'CXCL8', 'TNF', 'CXCL10', 'CCL2', 'VEGFA', 'COL1A1', 'COL3A1']
cov = {}
print(f'  {"TF":<10} {"V3(all)":>8} {"V2(active)":>11} {"V1(strict)":>11}')
h1, h2, h3 = set(T1), set(T2), set(T3)
for t in DFU_TF:
    cov[t] = {'V1': t in h1, 'V2': t in h2, 'V3': t in h3}
    if t in h3:
        print(f'  {t:<10} {"✓":>8} {"✓" if t in h2 else "-":>11} {"✓" if t in h1 else "-":>11}')
n3 = sum(v['V3'] for v in cov.values()); n2 = sum(v['V2'] for v in cov.values())
n1 = sum(v['V1'] for v in cov.values())
print(f'\n  汇总: V3={n3}/{len(DFU_TF)}  V2={n2}/{len(DFU_TF)}  V1={n1}/{len(DFU_TF)}')

meta['versions'] = res
meta['dfu_tf_coverage'] = cov
meta['coverage_summary'] = {'V3': n3, 'V2': n2, 'V1': n1, 'total': len(DFU_TF)}
with open(OUT / 'S2b_extract_multi.json', 'w', encoding='utf-8') as f:
    json.dump(meta, f, indent=2, ensure_ascii=False)
print('\n[saved] S2b_extract_multi.json')
