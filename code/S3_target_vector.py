# -*- coding: utf-8 -*-
"""
S3: 构建「DFU → 正常皮肤」的目标方向向量（不依赖 DE，不依赖结局标签）

设计动机（绕开 M10b 的零显著困境）：
  - 用 DE 定义目标 → M10b 证明患者级 FDR<0.05 = 0 基因  ❌
  - 用结局标签（healed/not_healed）→ n=14，M7/M8 推翻  ❌
  - 用"距 Healthy 欧氏距离" → M13/M14 证明深度敏感  ⚠️

本次采用【判别方向】方案：
  1. 在 wound_bed 内，用 4 组标签做【正则化多分类 / 二分类】，
     取"Healthy vs Disease"的判别方向 w（在基因空间）
  2. 该方向是【多基因加权】，比单基因 DE 稳健
  3. 用 LOPO 验证该方向的泛化性（真患者级）
  4. 输出三种目标向量，供 S4 投影：
     t_dir    : Healthy vs Disease 判别方向（逻辑回归系数）
     t_center : Healthy 患者中心 − Disease 患者中心（质心差）
     t_forearm: Forearm 皮肤中心 − DFU 创面中心（外部正常皮肤参照）

关键：所有向量都限制在 Southard ∩ DFU 的 4,827 交集基因上。
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
DATA = ROOT / "10_DFU_Dataset" / "dfu_outcome_dataset" / "GSE165816_sc.h5ad"
D = OUT / "data_southard"

# ---- 载入 Southard 交集基因（决定工作基因空间）----
pv3 = np.load(D / 'perturb_V3.npz', allow_pickle=True)
SGENES = pv3['genes'].astype(str)
print(f'Southard∩DFU 交集基因: {len(SGENES)}')

# ---- 载入 DFU 数据 ----
a = ad.read_h5ad(DATA)
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
genes = a.var_names.astype(str)
gi = {g: i for i, g in enumerate(genes)}

# 取交集基因的列索引
sel = [g for g in SGENES if g in gi]
SEL_IDX = np.array([gi[g] for g in sel])
print(f'实际可用（DFU 基因空间内）: {len(sel)}')

obs = a.obs.copy()
tc = obs['tissue_class'].astype(str).values
grp = obs['group'].astype(str).values
pat = obs['subject_hint'].astype(str).values

# ---- 成纤维细胞门控（与 M7-M14 一致）----
FIB = ['DCN', 'LUM', 'PDGFRA', 'COL1A1', 'COL1A2', 'FAP', 'THY1']
fib_cols = [gi[g] for g in FIB if g in gi]
sub = np.asarray(X[:, fib_cols].todense())
subz = (sub - sub.mean(0)) / (sub.std(0) + 1e-9)
fs = subz.mean(1)
fib = fs >= np.quantile(fs, .75)
print(f'成纤维门控细胞: {fib.sum():,}')

# ---- 构建患者级 pseudobulk（在交集基因上）----
# 三组：wound_bed (Healthy/DiabeticNoDFU/Healer/NonHealer) + forearm
def pseudobulk(mask_cells, label_fn):
    pats, rows, labs, tiss = [], [], [], []
    pat_sel = pat[mask_cells]
    grp_sel = grp[mask_cells]
    tc_sel = tc[mask_cells]
    Xs = X[mask_cells][:, SEL_IDX]
    for p in np.unique(pat_sel):
        mm = pat_sel == p
        if mm.sum() < 30:
            continue
        v = np.asarray(Xs[mm].mean(0)).ravel()
        if not np.all(np.isfinite(v)):
            continue
        rows.append(v)
        pats.append(p)
        labs.append(pd.Series(grp_sel[mm]).mode()[0])
        tiss.append(pd.Series(tc_sel[mm]).mode()[0])
    return np.array(rows), np.array(pats), np.array(labs), np.array(tiss)

# A: wound_bed 成纤维
wb = (tc == 'wound_bed') & fib
Pw, pw, lw, tw = pseudobulk(wb, None)
print(f'wound_bed pseudobulk: {Pw.shape}  '
      f'groups={dict(pd.Series(lw).value_counts())}')

# B: forearm（正常皮肤）成纤维
fb = (tc == 'non_wound_skin') & fib
Pf, pf, lf, tf = pseudobulk(fb, None)
print(f'forearm  pseudobulk: {Pf.shape}  '
      f'groups={dict(pd.Series(lf).value_counts())}')

res = {'n_intersect': int(len(sel)), 'genes': sel}

# ==================================================================
# 目标向量 1: t_center —— Healthy 质心 − Disease 质心
# ==================================================================
print('\n' + '=' * 78)
print('T1 质心差：Healthy 患者质心 − Disease 患者质心')
print('=' * 78)
h = lw == 'Healthy'
d = np.isin(lw, ['DiabeticNoDFU', 'Healer', 'NonHealer'])
t_center = Pw[h].mean(0) - Pw[d].mean(0)
t_center = t_center / (np.linalg.norm(t_center) + 1e-12)
print(f'  n_healthy={h.sum()}  n_disease={d.sum()}  '
      f'|t_center|={np.linalg.norm(Pw[h].mean(0)-Pw[d].mean(0)):.3f}')
res['t_center_norm'] = float(np.linalg.norm(Pw[h].mean(0) - Pw[d].mean(0)))

# ==================================================================
# 目标向量 2: t_forearm —— 正常前臂皮肤 − DFU 创面
# ==================================================================
print('\n' + '=' * 78)
print('T2 外部正常皮肤：Forearm 质心 − wound_bed 全体质心')
print('=' * 78)
t_fore = Pf.mean(0) - Pw.mean(0)
t_fore = t_fore / (np.linalg.norm(t_fore) + 1e-12)
print(f'  n_forearm={len(Pf)}  n_woundbed={len(Pw)}  '
      f'|raw|={np.linalg.norm(Pf.mean(0)-Pw.mean(0)):.3f}')
res['t_forearm_norm'] = float(np.linalg.norm(Pf.mean(0) - Pw.mean(0)))

# ==================================================================
# 目标向量 3: t_dir —— 判别方向（逻辑回归系数，LOPO 验证）
# ==================================================================
print('\n' + '=' * 78)
print('T3 判别方向：Healthy vs Disease 逻辑回归（LOPO 验证）')
print('=' * 78)

# 3a 在 SVD 空间做 LOPO 验证（证明方向真实泛化）
def lopo_dir(P, y, n_pc=15):
    Z = TruncatedSVD(n_pc, random_state=0).fit_transform(
        StandardScaler().fit_transform(P))
    sc = np.full(len(y), np.nan)
    for i in range(len(y)):
        tr = np.arange(len(y)) != i
        if len(np.unique(y[tr])) < 2:
            continue
        c = LogisticRegression(max_iter=3000, C=1.0).fit(Z[tr], y[tr])
        sc[i] = c.predict_proba(Z[i:i + 1])[0, 1]
    ok = np.isfinite(sc)
    auc = roc_auc_score(y[ok], sc[ok]) if len(np.unique(y[ok])) > 1 else np.nan
    # 置换
    rng = np.random.default_rng(7)
    null = []
    for _ in range(1000):
        yl = rng.permutation(y[ok])
        if len(np.unique(yl)) > 1:
            null.append(roc_auc_score(yl, sc[ok]))
    pp = float((np.abs(np.array(null) - .5) >= abs(auc - .5)).mean())
    return float(auc), pp, int(ok.sum())

y_hd = (lw != 'Healthy').astype(int)   # 1 = disease
auc_hd, p_hd, n_hd = lopo_dir(Pw, y_hd)
print(f'  Healthy vs Disease  LOPO AUC={auc_hd:.4f}  perm_p={p_hd:.4f}  n={n_hd}')
res['T3_lopo_healthy_vs_disease'] = {'auc': auc_hd, 'perm_p': p_hd, 'n': n_hd}

# 3b 在**基因空间**直接拟合（全数据），取系数作为方向
sc = StandardScaler().fit(Pw)
Pz = sc.transform(Pw)
clf = LogisticRegression(max_iter=5000, C=0.3, penalty='l2').fit(Pz, y_hd)
w = clf.coef_.ravel() / (np.linalg.norm(clf.coef_) + 1e-12)
# 换算回"原始表达单位"的方向（供投影）
t_dir = w / (sc.scale_ + 1e-12)
t_dir = t_dir / (np.linalg.norm(t_dir) + 1e-12)
# 训练集自评（乐观，仅供对照）
print(f'  基因空间拟合方向：|w| L2=1，非零系数={int((np.abs(w)>1e-6).sum())}'
      f'/ {len(w)}')
print(f'  训练集自评 AUC（乐观）={roc_auc_score(y_hd, clf.decision_function(Pz)):.4f}')
res['T3_train_auc_optimistic'] = float(roc_auc_score(y_hd, clf.decision_function(Pz)))

# ==================================================================
# 方向相似性检查（三个目标向量是否一致？）
# ==================================================================
print('\n' + '=' * 78)
print('三个目标向量的余弦相似度（一致性检查）')
print('=' * 78)
def cos(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))
pairs = {
    'center vs forearm': cos(t_center, t_fore),
    'center vs dir': cos(t_center, t_dir),
    'forearm vs dir': cos(t_fore, t_dir),
}
for k, v in pairs.items():
    print(f'  {k:<24} cos={v:+.4f}')
res['vector_consistency'] = pairs

# ==================================================================
# 保存
# ==================================================================
np.savez_compressed(D / 'target_vectors.npz',
                    genes=np.array(sel),
                    t_center=t_center.astype(np.float32),
                    t_forearm=t_fore.astype(np.float32),
                    t_dir=t_dir.astype(np.float32),
                    Pw=Pw.astype(np.float32), lw=lw, pw=pw,
                    Pf=Pf.astype(np.float32), lf=lf, pf=pf)
print('\n[saved] data_southard/target_vectors.npz')

res['groups_woundbed'] = {str(k): int(v) for k, v in
                          pd.Series(lw).value_counts().items()}
res['groups_forearm'] = {str(k): int(v) for k, v in
                         pd.Series(lf).value_counts().items()}
with open(OUT / 'S3_target_vector.json', 'w', encoding='utf-8') as f:
    json.dump(res, f, indent=2, ensure_ascii=False)
print('[saved] S3_target_vector.json')
