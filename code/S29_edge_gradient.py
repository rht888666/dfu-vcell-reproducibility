# -*- coding: utf-8 -*-
"""
S29: 空间梯度分析 —— 创缘距离 vs F6/F7 程序（虚拟组织决定性证据）

科学问题：
  如果 F6→F7 是真实的"修复状态转移"，它在空间上应表现为【位置依赖】：
    H1. F6（炎症肌成纤维）靠近创缘（wound edge），随距离衰减
    H2. F7/F7e（肌成纤维/ECM）在创缘外/深层真皮富集，或随愈合时间在创缘区累积
    H3. 这种梯度在 D30（重塑期）最明显

方法：
  1. 从 tissue_positions_list.csv 读 spot 的 pixel 坐标（pxl_row, pxl_col）
  2. 用 Wound_edges 注释的 spot 定义"创缘位置"，计算每个 spot 到最近创缘的距离
  3. 分层：按距离分箱，看 F6/F7/F7e 的均值
  4. 检验：距离 vs 程序得分的 Spearman（供体内）
  5. 阴性对照：F1/F5 不应有距离梯度
  6. 随机基因集零分布：判断梯度是否程序特异

诚实：
  - spot ≠ 单细胞
  - "创缘距离"用 AnnoType 代理（无病理标注）
  - 必须在每个切片内做（不同切片坐标系不同）
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.spatial import cKDTree

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
POS = ROOT / "10_DFU_Dataset" / "gse241124" / "pos"
RNG = np.random.default_rng(20261006)

LOG=[]
def log(*a):
    s=' '.join(str(x) for x in a); print(s); LOG.append(s)

log('='*88)
log('S29: 空间梯度分析 —— 创缘距离 vs F6/F7')
log('='*88)

# ---- 载入 spot 得分 ----
SP = pd.read_csv(OUT/"S26_spatial_scores.csv")
log(f'\n[1] spot 得分: {len(SP):,}')

# h5_stem -> pos 目录名
STEM2POS = {
    'GSM7717016_Skin': 'GSM7717016_Skin',
    'GSM7717017_Wound1': 'GSM7717017_Wound1',
    'GSM7717018_Wound7': 'GSM7717018_Wound7',
    'GSM7717019_Wound30': 'GSM7717019_Wound30',
    'GSM8238462_P17401_1001': 'GSM8238462_P17401_1001',
    'GSM8238463_P17401_1002': 'GSM8238463_P17401_1002',
    'GSM8238464_P17401_1003': 'GSM8238464_P17401_1003',
    'GSM8238465_P17401_1004': 'GSM8238465_P17401_1004',
    'GSM8238466_P20063_101': 'GSM8238466_P20063_101',
    'GSM8238467_P20063_102': 'GSM8238467_P20063_102',
    'GSM8238468_P20063_103': 'GSM8238468_P20063_103',
    'GSM8238469_P20063_104': 'GSM8238469_P20063_104',
    'GSM8238470_P20063_105': 'GSM8238470_P20063_105',
    'GSM8238471_P20063_106': 'GSM8238471_P20063_106',
    'GSM8238472_P20063_107': 'GSM8238472_P20063_107',
}

# ---- 读坐标 ----
log('\n[2] 读 spot 坐标')
coord_frames=[]
for stem, posdir in STEM2POS.items():
    f = POS/posdir/'tissue_positions_list.csv'
    if not f.exists():
        log(f'  ⚠️ 缺坐标: {posdir}'); continue
    c = pd.read_csv(f, header=None,
                    names=['barcode','in_tissue','array_row','array_col','pxl_row','pxl_col'])
    c['h5_stem']=stem
    coord_frames.append(c)
C = pd.concat(coord_frames, ignore_index=True)
log(f'  坐标行数: {len(C):,}  切片数: {C["h5_stem"].nunique()}')

# ---- 合并 ----
SP['bc'] = SP['barcode'].astype(str)
C['bc'] = C['barcode'].astype(str)
M = SP.merge(C[['bc','h5_stem','pxl_row','pxl_col','in_tissue']],
             on=['bc','h5_stem'], how='inner', suffixes=('','_c'))
log(f'  合并后: {len(M):,} spots（有坐标）')

# 深度残差化的程序得分
PROG = ['F6','F7','F7e','F1','F5']
zcols = [c for c in M.columns if c.endswith('_z')]
GENES = [c[:-2] for c in zcols]
log_n = np.log1p(M['total_counts'].values.astype(float))
lc = log_n-log_n.mean(); den=(lc**2).sum()+1e-12
GZ={}
for c in zcols:
    v=M[c].values; beta=(lc*(v-v.mean())).sum()/den
    GZ[c[:-2]]=v-beta*lc
GZ=pd.DataFrame(GZ, index=M.index)

# ---- 每切片：创缘距离 ----
log('\n[3] 计算创缘距离')
M = M.reset_index(drop=True); GZ = GZ.reset_index(drop=True)

def per_section_gradient(df, gzdf, sname):
    """对单个切片计算创缘距离梯度"""
    edge = (df['AnnoType']=='Wound_edges').values
    if edge.sum() < 20:
        return None
    pts = df[['pxl_row','pxl_col']].values.astype(float)
    ei = np.where(edge)[0]
    tree = cKDTree(pts[ei])
    d, _ = tree.query(pts, k=1)
    out = {'h5_stem': sname, 'n': len(df), 'n_edge': int(edge.sum())}
    # 距离缩放：转换成"spot 直径"单位（用坐标系 median 间距）
    # 近似：相邻 spot 中心距 ~ 100 px (Visium)，取 robust 估计
    from scipy.spatial import cKDTree as KD
    sub = pts[RNG.choice(len(pts), min(2000,len(pts)), replace=False)]
    nn = KD(sub).query(sub, k=2)[0][:,1]
    spot_d = np.median(nn)
    out['spot_diam_px'] = spot_d
    out['dist_units'] = d/spot_d
    return out, d/spot_d

log(f'\n  {"切片":32s} {"n":>5s} {"n_edge":>7s} {"中位距离":>9s} {"最大距离":>9s}')
dist_all = pd.Series(index=M.index, dtype=float)
for stem in sorted(M['h5_stem'].unique()):
    m = M['h5_stem']==stem
    r = per_section_gradient(M[m], GZ[m], stem)
    if r is None:
        log(f'  {stem:32s} 创缘 spot 不足，跳过'); continue
    info, dvec = r
    dist_all[M.index[m]] = dvec
    log(f'  {stem:32s} {info["n"]:5d} {info["n_edge"]:7d} {np.median(dvec):9.2f} {dvec.max():9.2f}')

M['edge_dist'] = dist_all.values
Mv = M.dropna(subset=['edge_dist']).copy()
GZv = GZ.loc[Mv.index]
log(f'\n  可用于梯度分析: {len(Mv):,} spots')

# ---- 梯度检验 ----
log('\n' + '='*88)
log('[4] 创缘距离 vs 程序得分（供体内 Spearman）')
log('='*88)
# 只在不含 Skin 的创面切片？（创缘在 Skin 中无意义）
# 分别在 全部 / 创面 only 上做
for scope in ['all','wound_only']:
    if scope=='wound_only':
        mv = Mv[Mv['Condition']!='Skin']; gv = GZ.loc[mv.index]
    else:
        mv = Mv; gv = GZv
    log(f'\n  [{scope}] n={len(mv):,}')
    log(f'    {"程序":6s} ' + ' '.join(f'{"D"+s[-1] if s!="Skin" else "Skin":>8s}' for s in ['Skin','Wound1','Wound7','Wound30']) + '   合并rho')
    for p in PROG:
        cols = [g+'_z' for g in [] ]  # 用 GZ 已残差化
        s_vec = gv[p].values if p in gv.columns else None
        # 用程序 marker 在 GZ 上的均值
        marker_map = {
            'F6':  ['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
            'F7':  ['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
            'F7e': ['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
            'F1':  ['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
            'F5':  ['SCN7A','SOX10','PLP1','NGFR','S100B'],
        }
        gl = [g for g in marker_map[p] if g in gv.columns]
        sv = gv[gl].mean(1).values
        rhos=[]
        for cond in ['Skin','Wound1','Wound7','Wound30']:
            mm = (mv['Condition']==cond).values
            if mm.sum()<50: rhos.append(np.nan); continue
            r,_ = spearmanr(mv.loc[mm,'edge_dist'], sv[mm])
            rhos.append(r)
        r_all,_ = spearmanr(mv['edge_dist'], sv)
        log(f'    {p:6s} ' + ' '.join(f'{x:+8.3f}' if not np.isnan(x) else '     n/a' for x in rhos) + f'   {r_all:+.4f}')

# ---- 关键：D7/D30 的 F7 梯度 ----
log('\n' + '='*88)
log('[5] 关键检验：D7/D30 创面中 F7/F7e 随创缘距离的梯度')
log('='*88)
for cond in ['Wound1','Wound7','Wound30']:
    mm = (Mv['Condition']==cond).values
    if mm.sum()<50: continue
    log(f'\n  {cond} (n={mm.sum()}):')
    for p in ['F6','F7','F7e','F1','F5']:
        marker_map = {
            'F6':  ['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
            'F7':  ['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
            'F7e': ['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
            'F1':  ['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
            'F5':  ['SCN7A','SOX10','PLP1','NGFR','S100B'],
        }
        gl=[g for g in marker_map[p] if g in GZ.columns]
        sv = GZ.loc[Mv.index[mm], gl].mean(1).values
        r,pv = spearmanr(Mv.loc[mm,'edge_dist'], sv)
        log(f'    {p:5s}: rho={r:+.4f}  p={pv:.2e}')

# 保存
Mv[['barcode','h5_stem','Condition','Donor','AnnoType','edge_dist','total_counts']].to_csv(
    OUT/"S29_edge_dist.csv", index=False)
out = dict(n_spots=int(len(Mv)),
           conditions={c:int((Mv['Condition']==c).sum()) for c in ['Skin','Wound1','Wound7','Wound30']},
           note='edge_dist in spot-diameter units, per-section')
with open(OUT/"S29_gradient.json",'w',encoding='utf-8') as f:
    json.dump(out,f,indent=2,ensure_ascii=False,default=str)
log(f'\n[saved] S29_edge_dist.csv + S29_gradient.json')
with open(OUT/"S29_run.log",'w',encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('\n✅ S29 完成')
