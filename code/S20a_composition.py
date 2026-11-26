# -*- coding: utf-8 -*-
"""
S20a: GSE241132 细胞组成的时间轨迹（不依赖表达矩阵，仅用 metadata）

目的：
  先看【成纤维细胞占比】以及【成纤维亚型】随时间的变化，
  为 S20 表达分析提供组成层面的背景。

注意：这是细胞【组成】分析，不是【状态】分析。状态分析见 S20。
"""
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
G = ROOT / "10_DFU_Dataset" / "gse241132"

df = pd.read_csv(G / 'cell_metadata.txt.gz', sep='\t')
print('shape:', df.shape)

# 统一派生
df['donor'] = df['Patient']
df['timepoint'] = df['Condition']
df['is_fib'] = df['newMainCellTypes'] == 'Fibroblast'

print('\n' + '=' * 78)
print('1. 各供体 × 时间点 的细胞总数')
print('=' * 78)
ct = df.pivot_table(index='donor', columns='timepoint', values='barcode', aggfunc='count')
order = ['Skin', 'Wound1', 'Wound7', 'Wound30']
ct = ct.reindex(columns=[c for c in order if c in ct.columns])
print(ct.to_string())
print('\n每时间点合计:')
print(ct.sum(0).to_string())

print('\n' + '=' * 78)
print('2. 【成纤维细胞占比】的时间轨迹（每供体内比例，避免样本量偏差）')
print('=' * 78)
rows = []
for (d, t), sub in df.groupby(['donor', 'timepoint']):
    rows.append({'donor': d, 'timepoint': t, 'n': len(sub),
                 'fib_n': int(sub['is_fib'].sum()),
                 'fib_pct': sub['is_fib'].mean()})
P = pd.DataFrame(rows)
pv = P.pivot_table(index='donor', columns='timepoint', values='fib_pct')
pv = pv.reindex(columns=[c for c in order if c in pv.columns])
print('成纤维占比:')
print((pv * 100).round(2).to_string())
print('\n均值:')
print((pv.mean(0) * 100).round(2).to_string())

print('\n' + '=' * 78)
print('3. 成纤维细胞【绝对数】的时间轨迹')
print('=' * 78)
pv2 = P.pivot_table(index='donor', columns='timepoint', values='fib_n')
pv2 = pv2.reindex(columns=[c for c in order if c in pv2.columns])
print(pv2.to_string())
print('\n均值:')
print(pv2.mean(0).round(1).to_string())

print('\n' + '=' * 78)
print('4. 各时间点的细胞类型构成（合并 3 供体）')
print('=' * 78)
comp = df.pivot_table(index='newMainCellTypes', columns='timepoint',
                      values='barcode', aggfunc='count').fillna(0)
comp = comp.reindex(columns=[c for c in order if c in comp.columns])
# 归一化为比例
comp_pct = comp / comp.sum(0) * 100
print(comp_pct.round(2).to_string())

print('\n' + '=' * 78)
print('5. 成纤维细胞【亚型】的时间轨迹（newCellTypes 中成纤维相关）')
print('=' * 78)
fib = df[df['is_fib']]
print('成纤维亚型:', sorted(fib['newCellTypes'].dropna().unique()))
sub_ct = fib.pivot_table(index='newCellTypes', columns='timepoint',
                         values='barcode', aggfunc='count').fillna(0)
sub_ct = sub_ct.reindex(columns=[c for c in order if c in sub_ct.columns])
sub_pct = sub_ct / sub_ct.sum(0) * 100
print('\n各亚型占比（%）：')
print(sub_pct.round(2).to_string())

# 保存
comp_pct.to_csv(OUT / 'S20a_celltype_time.csv', encoding='utf-8-sig')
sub_pct.to_csv(OUT / 'S20a_fib_subtype_time.csv', encoding='utf-8-sig')
P.to_csv(OUT / 'S20a_fib_pct.csv', index=False, encoding='utf-8-sig')
print('\n[saved] S20a_*.csv')

print('\n' + '=' * 78)
print('6. 深度（nCount_RNA）随时间的变化（重要：影响后续状态分析）')
print('=' * 78)
pv3 = df.pivot_table(index='donor', columns='timepoint', values='nCount_RNA', aggfunc='mean')
pv3 = pv3.reindex(columns=[c for c in order if c in pv3.columns])
print('平均 nCount_RNA:')
print(pv3.round(0).to_string())
print('\n各时间点均值:')
print(pv3.mean(0).round(0).to_string())
fibd = fib.pivot_table(index='donor', columns='timepoint', values='nCount_RNA', aggfunc='mean')
fibd = fibd.reindex(columns=[c for c in order if c in fibd.columns])
print('\n仅成纤维细胞的平均 nCount_RNA:')
print(fibd.round(0).to_string())
print('\n仅成纤维各时间点均值:')
print(fibd.mean(0).round(0).to_string())
