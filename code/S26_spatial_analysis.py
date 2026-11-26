# -*- coding: utf-8 -*-
"""
S26: GSE241124 空间转录组 —— F6/F7 程序的空间分布与时空演进
（"虚拟组织"维度：把 S20 的纵向状态转移投射到空间）

数据：GSE241124 Visium，15 切片（4 供体 × ~4 条件），22,915 spots
  - 每个 spot 有 F6/F7 程序得分（用 S20 相同 marker 集，spot 级 z-score）
  - metadata 提供 AnnoType（FB1/FB2/Wound_edges/Dermis1...）、Condition、Donor

核心问题：
  Q1. F6/F7 得分在【4 个时间条件】上如何变化？（空间版纵向轨迹，检验与 scRNA 一致）
  Q2. F6/F7 得分在【空间注释类别】上如何分布？（F6 是否富集于创缘/免疫区，F7 于真皮）
  Q3. FB1/FB2 的 spot 占比是否随条件变化？（独立验证 scRNA 的 FB-I D7 爆发）
  Q4. F6/F7 的"空间重心"是否随条件移动？（虚拟组织层的状态转移）
  Q5. 深度控制：spot 级 nCount_Spatial 校正，并做阴性对照 F1/F5

诚实原则：
  - spot 不是单细胞（~20 cells/spot），得分是"组织区域状态"，非细胞状态
  - 必须做深度校正（Visium spot 深度差异大）
  - F1/F5 阴性对照必须无趋势
"""
import json, gzip, warnings
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
G = ROOT / "10_DFU_Dataset" / "gse241124"
RNG = np.random.default_rng(20261006)

LOG = []
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOG.append(s)

log('=' * 84)
log('S26: GSE241124 空间转录组 —— F6/F7 空间分布与时空演进')
log('=' * 84)

# 程序 marker（与 S20 一致）
PROG = {
    'F6':  ['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
    'F7':  ['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
    'F7e': ['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
    'F1':  ['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
    'F5':  ['SCN7A','SOX10','PLP1','NGFR','S100B'],
}
ALL_MARKERS = sorted({g for v in PROG.values() for g in v})

# =========================================================
# 1. 读 metadata
# =========================================================
log('\n[1] 读 metadata')
meta = pd.read_csv(G / 'meta.txt.gz', sep='\t')
log(f'  {meta.shape[0]:,} spots × {meta.shape[1]} 列')
log(f'  条件: {sorted(meta["Condition"].unique())}')
log(f'  供体: {sorted(meta["Donor"].unique())}')
log(f'  注释类型 ({meta["AnnoType"].nunique()}): {sorted(meta["AnnoType"].unique())}')
log(f'  Condition × Donor 交叉:')
log(meta.groupby(['Donor','Condition']).size().unstack(fill_value=0).to_string())

# =========================================================
# 2. 读 h5，计算 spot 级程序得分
# =========================================================
log('\n[2] 读 15 个 h5 矩阵，计算 spot 级得分')

def read_h5_markers(hp, markers):
    """读 10x h5，只保留 marker 基因，返回 (DataFrame barcodes × markers, 总 counts)"""
    with h5py.File(str(hp), 'r') as f:
        M = f['matrix']
        shape = tuple(M['shape'][()])
        names = np.array([x.decode() if isinstance(x,bytes) else str(x) for x in M['features/name'][()]])
        bcs = np.array([x.decode() if isinstance(x,bytes) else str(x) for x in M['barcodes'][()]])
        data = M['data'][()]; indices = M['indices'][()]; indptr = M['indptr'][()]
    X = sp.csc_matrix((data, indices, indptr), shape=shape)
    total = np.asarray(X.sum(0)).ravel()
    keep_pos, keep_name, seen = [], [], set()
    for i, nm in enumerate(names):
        if nm in markers and nm not in seen:
            seen.add(nm); keep_pos.append(i); keep_name.append(nm)
    sub = X[keep_pos, :].toarray().T   # (n_spots, n_markers)
    return pd.DataFrame(sub, columns=keep_name), bcs, total

# 样本名映射：h5 文件名 -> orig.ident
fmap = {}
for f_ in sorted((G/'h5').glob('*.h5')):
    stem = f_.name.replace('_filtered_feature_bc_matrix.h5','')
    fmap[stem] = f_
log(f'  {len(fmap)} 个 h5 文件')

frames = []
for stem, hp in fmap.items():
    try:
        df, bcs, total = read_h5_markers(hp, set(ALL_MARKERS))
    except Exception as e:
        log(f'  ⚠️ {stem}: 读取失败 ({type(e).__name__}: {str(e)[:60]}) — 跳过')
        continue
    df['barcode'] = bcs
    df['total_counts'] = total
    df['h5_stem'] = stem
    log(f'  {stem:45s} {len(bcs):5d} spots, marker 命中 {df.shape[1]-3}/{len(ALL_MARKERS)}')
    frames.append(df)

SP = pd.concat(frames, ignore_index=True)
log(f'\n  合并: {len(SP):,} spots')

# =========================================================
# 3. 合并 metadata（按 barcode 后缀匹配）
# =========================================================
log('\n[3] 合并 metadata')
# metadata barcode 形如 Donor1_Skin_AAACAAGTATCTCCCA-1
# h5 barcode 形如 AAACAAGTATCTCCCA-1
meta['bc_suffix'] = meta['barcode'].astype(str).str.split('_').str[-1]
meta_key = meta.set_index('bc_suffix')

SP['bc_suffix'] = SP['barcode'].astype(str).str.split('_').str[-1]
# 因为不同切片 barcode 可能重复，需按 Sample_name 匹配
# 先尝试按 bc_suffix + orig.ident 匹配

# 建立 (orig.ident, bc_suffix) -> meta 行
meta_idx = {}
for _, r in meta.iterrows():
    meta_idx[(str(r['orig.ident']), str(r['bc_suffix']))] = r

# h5 stem -> orig.ident 的映射表（手工/自动）
# ⚠️ S31 审计修正：原表把 GSM7717016-19 错映射为 P17401_1001-1004（当时误信 GEO 标题
#    的 "Donor4-*" 与 P17401 编号冲突，实际二者不冲突：GSM7717016-19 是 Donor4，
#    orig.ident = P20063_109-112）。正确映射由 S31b 反向唯一性检验（16/16 样本 100% 覆盖）给出。
STEM2IDENT = {
    # --- 早期上传批次（2023-08），文件名无 Pxxx 编号，实为 Donor4 ---
    'GSM7717016_Skin':    'P20063_109',
    'GSM7717017_Wound1':  'P20063_110',
    'GSM7717018_Wound7':  'P20063_111',
    'GSM7717019_Wound30': 'P20063_112',
    # --- 主批次（2024-04）---
    'GSM8238462_P17401_1001': 'P17401_1001',
    'GSM8238463_P17401_1002': 'P17401_1002',
    'GSM8238464_P17401_1003': 'P17401_1003',
    'GSM8238465_P17401_1004': 'P17401_1004',
    'GSM8238466_P20063_101': 'P20063_101',
    'GSM8238467_P20063_102': 'P20063_102',
    'GSM8238468_P20063_103': 'P20063_103',
    'GSM8238469_P20063_104': 'P20063_104',
    'GSM8238470_P20063_105': 'P20063_105',
    'GSM8238471_P20063_106': 'P20063_106',
    'GSM8238472_P20063_107': 'P20063_107',
    'GSM8238473_P20063_108': 'P20063_108',
}
SP['orig_ident'] = SP['h5_stem'].map(STEM2IDENT)

n_match = 0
annos, conds, donors_, patients = [], [], [], []
for _, r in SP.iterrows():
    key = (r['orig_ident'], r['bc_suffix'])
    if key in meta_idx:
        mr = meta_idx[key]
        annos.append(mr['AnnoType']); conds.append(mr['Condition'])
        donors_.append(mr['Donor']); patients.append(mr['Patient'])
        n_match += 1
    else:
        annos.append(None); conds.append(None); donors_.append(None); patients.append(None)
SP['AnnoType'] = annos; SP['Condition'] = conds
SP['Donor'] = donors_; SP['Patient'] = patients
log(f'  匹配成功: {n_match:,}/{len(SP):,} ({n_match/len(SP):.1%})')

SPm = SP.dropna(subset=['Condition']).copy()
log(f'  可用 spots: {len(SPm):,}')

# =========================================================
# 4. 计算程序得分（深度校正）
# =========================================================
log('\n[4] 计算程序得分（log-normalize + z-score + 深度残差化）')
mk_present = [g for g in ALL_MARKERS if g in SPm.columns]
log(f'  实际可用 marker: {len(mk_present)}/{len(ALL_MARKERS)}')

# CPM-like normalize: log1p(count / total * 1e4)
tot = SPm['total_counts'].values.astype(float)
tot[tot == 0] = 1
for g in mk_present:
    SPm[g+'_n'] = np.log1p(SPm[g].values / tot * 1e4)

# z-score across all spots
for g in mk_present:
    v = SPm[g+'_n'].values
    SPm[g+'_z'] = (v - v.mean()) / (v.std() + 1e-9)

log_n = np.log1p(tot)

for pname, glist in PROG.items():
    cols = [g+'_z' for g in glist if g+'_z' in SPm.columns]
    if not cols: continue
    s = SPm[cols].mean(1).values
    # 深度残差化（用 log total counts）
    lc = log_n - log_n.mean()
    beta = (lc * (s - s.mean())).sum() / ((lc**2).sum() + 1e-12)
    SPm[pname] = s
    SPm[pname+'_res'] = s - beta * lc
    r_raw = spearmanr(s, tot)[0]
    r_res = spearmanr(SPm[pname+'_res'].values, tot)[0]
    log(f'  {pname:4s}: n_genes={len(cols)}  rho(raw,tot)={r_raw:+.4f}  rho(res,tot)={r_res:+.4f}')

# =========================================================
# Q1. 四条件下程序得分（空间版纵向）
# =========================================================
log('\n' + '=' * 84)
log('Q1. 四条件程序得分（空间版纵向轨迹）')
log('=' * 84)
COND = ['Skin','Wound1','Wound7','Wound30']
tbl = SPm.groupby('Condition')[['F6','F7','F7e','F1','F5']].mean().reindex([c for c in COND if c in SPm['Condition'].unique()])
log('\n  [原始得分]')
log(tbl.round(4).to_string())
tbl_r = SPm.groupby('Condition')[['F6_res','F7_res']].mean().reindex([c for c in COND if c in SPm['Condition'].unique()])
log('\n  [深度残差化]')
log(tbl_r.round(4).to_string())

# n spots
log('\n  spots/条件: ' + ', '.join(f'{c}={int((SPm["Condition"]==c).sum())}' for c in COND if c in SPm['Condition'].unique()))

# 趋势检验（供体内，条件作为序数）
ord_map = {'Skin':0,'Wound1':1,'Wound7':2,'Wound30':3}
SPm['cond_ord'] = SPm['Condition'].map(ord_map)

log('\n  [趋势 Spearman]（每供体内）')
for p in ['F6','F7','F7e','F1','F5','F6_res','F7_res']:
    if p not in SPm.columns: continue
    rhos = []
    for d in sorted(SPm['Donor'].dropna().unique()):
        m = SPm['Donor'] == d
        if m.sum() < 50: continue
        r, _ = spearmanr(SPm.loc[m,'cond_ord'], SPm.loc[m,p])
        rhos.append(r)
    rhos = np.array(rhos)
    log(f'    {p:8s} 供体 rho: ' + ', '.join(f'{x:+.3f}' for x in rhos) +
        f'   均值={rhos.mean():+.3f}  同号={sum(1 for x in rhos if x>0) if rhos.mean()>0 else sum(1 for x in rhos if x<0)}/{len(rhos)}')

# =========================================================
# Q2. 空间注释类别上的 F6/F7 分布
# =========================================================
log('\n' + '=' * 84)
log('Q2. 空间注释类别的 F6/F7 得分（空间定位）')
log('=' * 84)
anno_tbl = SPm.groupby('AnnoType').agg(
    n=('F6','size'), F6=('F6','mean'), F7=('F7','mean'),
    F7e=('F7e','mean'), F1=('F1','mean'), F5=('F5','mean')
).sort_values('n', ascending=False)
log('\n' + anno_tbl.round(4).to_string())

# 关键：F6/F7 最高的注释类别
log('\n  F6 最高的 5 个注释:')
log(anno_tbl.nlargest(5,'F6')[['n','F6']].round(4).to_string())
log('\n  F7 最高的 5 个注释:')
log(anno_tbl.nlargest(5,'F7')[['n','F7']].round(4).to_string())

# =========================================================
# Q3. FB1/FB2 spot 占比随条件变化
# =========================================================
log('\n' + '=' * 84)
log('Q3. FB1/FB2 spot 占比随条件变化（独立验证 scRNA FB-I 爆发）')
log('=' * 84)
for anno in ['FB1','FB2','Sweatgland_FB','Wound_edges','Dermis1','Papillary_dermis']:
    sub = SPm[SPm['AnnoType']==anno]
    if len(sub)==0: continue
    tot_by_cond = SPm.groupby('Condition').size()
    frac = sub.groupby('Condition').size() / tot_by_cond
    frac = frac.reindex([c for c in COND if c in frac.index])
    log(f'  {anno:18s} 占比: ' + ', '.join(f'{c}={frac[c]:.3%}' for c in frac.index))

# =========================================================
# Q4. F6/F7 的"空间重心"是否随条件移动
# =========================================================
log('\n' + '=' * 84)
log('Q4. 空间重心移动（需要坐标，检查是否有）')
log('=' * 84)
log('  ⚠️ 本 metadata 不含 spot 空间坐标（x/y），仅含 AnnoType —— 空间坐标在 images.zip 内')
log('  → Q4 用"注释组成"替代：用 AnnoType 分布刻画空间结构变化')

# =========================================================
# Q5. 阴性对照
# =========================================================
log('\n' + '=' * 84)
log('Q5. 阴性对照 F1/F5')
log('=' * 84)
for p in ['F1','F5']:
    rhos = []
    for d in sorted(SPm['Donor'].dropna().unique()):
        m = SPm['Donor'] == d
        if m.sum() < 50: continue
        r, _ = spearmanr(SPm.loc[m,'cond_ord'], SPm.loc[m,p])
        rhos.append(r)
    log(f'  {p}: 供体 rho 均值={np.mean(rhos):+.4f}  (应接近 0)')

# 保存
SPm.to_csv(OUT / "S26_spatial_scores.csv", index=False)
tbl.to_csv(OUT / "S26_condition_table.csv")
anno_tbl.to_csv(OUT / "S26_annotation_table.csv")
out = dict(
    n_spots=int(len(SPm)),
    conditions={c:int((SPm['Condition']==c).sum()) for c in COND if c in SPm['Condition'].unique()},
    condition_scores=tbl.round(5).to_dict(),
    annotation_scores=anno_tbl.round(5).to_dict(),
)
with open(OUT / "S26_spatial.json",'w',encoding='utf-8') as f:
    json.dump(out,f,indent=2,ensure_ascii=False,default=str)
log(f'\n[saved] S26_spatial_scores.csv + S26_spatial.json')
with open(OUT / "S26_run.log",'w',encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('\n✅ S26 完成')
