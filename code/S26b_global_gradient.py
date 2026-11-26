# -*- coding: utf-8 -*-
"""
S26b: 空间结果的深度/全局梯度裁决 —— F7 趋势是"程序特异"还是"全局伪影"？

S26 发现：
  🟢 F6 在 D1 峰、F7/F7e 在 D30 峰（与 scRNA 一致，独立技术）
  🟢 F7 富集于 FB1 spots (+0.613)；F6 富集于 Wound_edges/Immune
  🔴 阴性对照失败：F1 rho=+0.307, F5 rho=+0.140（也随条件上升）

必须回答：
  A. 全局梯度是否存在？—— 随机基因集的"条件趋势"零分布
  B. F7 的趋势是否【显著大于】随机基因集 / F1/F5？
  C. 在【双差】框架下：F7趋势 − 阴性对照趋势 是否显著？
  D. 空间特异性：F7 富集于 FB1 是否超过"任何高表达程序都富集于 FB1"的基线？

设计：
  - 构建 1000 个随机 5-7 基因集，各自算条件趋势 rho → 零分布
  - F7/F7e/F6 的 rho 在零分布中的位置（经验 p）
  - 双差法：逐 spot，F7 - F1（或 F7 - F5），看双差的条件趋势
  - 阳性-阴性对照配对法：对每个供体，比较 F7 rho 与 F1/F5 rho
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"
RNG = np.random.default_rng(20261006)

LOG = []
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOG.append(s)

log('=' * 84)
log('S26b: 空间结果的全局梯度裁决')
log('=' * 84)

SP = pd.read_csv(OUT / "S26_spatial_scores.csv")
log(f'\n[data] {len(SP):,} spots')
ord_map = {'Skin':0,'Wound1':1,'Wound7':2,'Wound30':3}
SP['cond_ord'] = SP['Condition'].map(ord_map)
donors = sorted(SP['Donor'].dropna().unique())
log(f'  供体: {donors}')

# 重建 gene-level z 矩阵（用于随机基因集）
zcols = [c for c in SP.columns if c.endswith('_z')]
GENES = [c[:-2] for c in zcols]
log(f'  可用基因 z 列: {len(GENES)} (marker 集)')

PROG = {
    'F6':  ['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
    'F7':  ['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
    'F7e': ['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
    'F1':  ['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
    'F5':  ['SCN7A','SOX10','PLP1','NGFR','S100B'],
}

def per_donor_rho(score_vec):
    """返回每供体的 (cond_ord vs score) Spearman rho 列表"""
    out = []
    for d in donors:
        m = (SP['Donor'] == d).values
        if m.sum() < 50: continue
        r, _ = spearmanr(SP.loc[m,'cond_ord'], score_vec[m])
        out.append(r)
    return np.array(out)

# 各程序
log('\n' + '=' * 84)
log('A. 各程序的条件趋势（供体内 rho）')
log('=' * 84)
prog_rhos = {}
for p, gl in PROG.items():
    cols = [g+'_z' for g in gl if g+'_z' in SP.columns]
    if not cols: continue
    # 深度残差化后的 score
    if p+'_res' in SP.columns:
        s = SP[p+'_res'].values
    else:
        s = SP[cols].mean(1).values
    rhos = per_donor_rho(s)
    prog_rhos[p] = rhos
    log(f'  {p:4s}: rho = ' + ', '.join(f'{x:+.3f}' for x in rhos) + f'   均值={rhos.mean():+.3f}')

# =========================================================
# B. 随机基因集零分布
# =========================================================
log('\n' + '=' * 84)
log('B. 随机基因集零分布（20 个 marker 基因随机抽 5-7 个）')
log('=' * 84)
# 深度残差化每个 marker 的 z
log_n = np.log1p(SP['total_counts'].values.astype(float))
lc = log_n - log_n.mean(); den = (lc**2).sum() + 1e-12
GZ = {}
for g in GENES:
    v = SP[g+'_z'].values
    beta = (lc*(v-v.mean())).sum()/den
    GZ[g] = v - beta*lc
GZ = pd.DataFrame(GZ)

n_rand = 2000
rand_rhos = []
sizes = [5,6,7]
for _ in range(n_rand):
    k = RNG.choice(sizes)
    gs = RNG.choice(list(GZ.columns), size=k, replace=False)
    s = GZ[gs].mean(1).values
    r = per_donor_rho(s)
    rand_rhos.append(r.mean())
rand_rhos = np.array(rand_rhos)
log(f'  零分布: mean={rand_rhos.mean():+.4f}, sd={rand_rhos.std():.4f}')
log(f'  分位: p50={np.percentile(rand_rhos,50):+.4f}, p90={np.percentile(rand_rhos,90):+.4f}, p95={np.percentile(rand_rhos,95):+.4f}, p99={np.percentile(rand_rhos,99):+.4f}')

log('\n  各程序 rho 在零分布中的经验 p (单侧):')
verdict = {}
for p, rhos in prog_rhos.items():
    obs = rhos.mean()
    p_emp = float((rand_rhos >= obs).sum() + 1) / (n_rand + 1)
    log(f'    {p:4s}: obs={obs:+.4f}  p_emp={p_emp:.4f}  {"✅ 超零分布" if p_emp<0.05 else "❌ 在零分布内"}')
    verdict[p] = dict(obs=float(obs), p_emp=p_emp)

# =========================================================
# C. 双差法：F7 - 阴性对照
# =========================================================
log('\n' + '=' * 84)
log('C. 双差法：F7 − 阴性对照（逐 spot）')
log('=' * 84)
for neg in ['F1','F5']:
    d = GZ[[g for g in PROG['F7'] if g in GZ.columns]].mean(1).values - GZ[[g for g in PROG[neg] if g in GZ.columns]].mean(1).values
    rhos = per_donor_rho(d)
    obs = rhos.mean()
    # 零分布：随机 F7-rand - rand2
    null = []
    for _ in range(2000):
        g1 = RNG.choice(list(GZ.columns), size=7, replace=False)
        g2 = RNG.choice(list(GZ.columns), size=5, replace=False)
        dd = GZ[g1].mean(1).values - GZ[g2].mean(1).values
        null.append(per_donor_rho(dd).mean())
    null = np.array(null)
    p = float((null >= obs).sum()+1)/(len(null)+1)
    log(f'  F7 − {neg}: rho=' + ', '.join(f'{x:+.3f}' for x in rhos) + f'  均值={obs:+.4f}')
    log(f'      零分布 sd={null.std():.4f}, p={p:.4f}  {"✅" if p<0.05 else "❌"}')

# =========================================================
# D. 配对检验：每供体 F7 rho vs F1/F5 rho
# =========================================================
log('\n' + '=' * 84)
log('D. 供体配对检验（F7 vs 阴性对照）')
log('=' * 84)
f7 = prog_rhos['F7']
for neg in ['F1','F5']:
    n_ = prog_rhos[neg]
    if len(f7)==len(n_):
        diff = f7 - n_
        log(f'  F7 − {neg} 逐供体差: ' + ', '.join(f'{x:+.3f}' for x in diff) + f'  均值={diff.mean():+.4f}')
        log(f'    同号 {sum(1 for x in diff if x>0)}/{len(diff)}')

# =========================================================
# E. 空间特异性：F7 富集于 FB1 是否显著
# =========================================================
log('\n' + '=' * 84)
log('E. 空间定位特异性：F7 富集于 FB1 的显著性')
log('=' * 84)
is_fb1 = (SP['AnnoType']=='FB1').values
is_fb2 = (SP['AnnoType']=='FB2').values
f7v = GZ[[g for g in PROG['F7'] if g in GZ.columns]].mean(1).values
log(f'  FB1 n={is_fb1.sum()}, FB2 n={is_fb2.sum()}')
log(f'  F7 在 FB1: {f7v[is_fb1].mean():+.4f}  在非FB1: {f7v[~is_fb1].mean():+.4f}  差={f7v[is_fb1].mean()-f7v[~is_fb1].mean():+.4f}')
log(f'  F7 在 FB2: {f7v[is_fb2].mean():+.4f}  在非FB2: {f7v[~is_fb2].mean():+.4f}  差={f7v[is_fb2].mean()-f7v[~is_fb2].mean():+.4f}')
# 置换：打乱 FB1 标签
obs = f7v[is_fb1].mean() - f7v[~is_fb1].mean()
cnt = 0
for _ in range(2000):
    perm = RNG.permutation(is_fb1)
    if (f7v[perm].mean() - f7v[~perm].mean()) >= obs: cnt += 1
log(f'  FB1 富集置换检验: obs={obs:+.4f}, p={(cnt+1)/2001:.4f}')

# 随机基因集在 FB1 的富集零分布（判断是否"任何基因都富集于FB1"）
rand_fb1 = []
for _ in range(2000):
    gs = RNG.choice(list(GZ.columns), size=7, replace=False)
    s = GZ[gs].mean(1).values
    rand_fb1.append(s[is_fb1].mean()-s[~is_fb1].mean())
rand_fb1 = np.array(rand_fb1)
p_fb1 = float((rand_fb1 >= obs).sum()+1)/2001
log(f'  随机基因集 FB1 富集零分布: sd={rand_fb1.std():.4f}, p(F7)={p_fb1:.4f}  {"✅ 特异" if p_fb1<0.05 else "❌ 非特异"}')

# 保存
out = dict(
    n_spots=int(len(SP)),
    donors=sorted(SP['Donor'].dropna().unique().tolist()),
    n_sections=int(SP['h5_stem'].nunique()) if 'h5_stem' in SP.columns else None,
    program_rhos={p: r.tolist() for p,r in prog_rhos.items()},
    null_rand=dict(mean=float(rand_rhos.mean()), sd=float(rand_rhos.std()),
                   p95=float(np.percentile(rand_rhos,95)), p99=float(np.percentile(rand_rhos,99))),
    empirical_p=verdict,
    fb1_enrichment=dict(obs=float(obs), p_perm=float((cnt+1)/2001),
                        rand_sd=float(rand_fb1.std()), p_specific=float(p_fb1)),
)
with open(OUT / "S26b_global_gradient.json",'w',encoding='utf-8') as f:
    json.dump(out,f,indent=2,ensure_ascii=False,default=str)
log(f'\n[saved] S26b_global_gradient.json')
with open(OUT / "S26b_run.log",'w',encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('\n✅ S26b 完成')
