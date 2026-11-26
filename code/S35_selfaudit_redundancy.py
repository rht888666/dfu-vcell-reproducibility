# -*- coding: utf-8 -*-
"""
S35 — 自审：空间分工是否为"同义反复"？

会诊最可能的攻击：**"F7 定位到 FB1（成纤维细胞）"是废话 —— F7 就是成纤维细胞 marker，
而 FB1 就是被标注为成纤维细胞的类别。这是循环论证。**

要证伪这个攻击，需要回答：
Q1. F7 在 FB1 的富集，是否显著高于 F7 在其他**成纤维细胞亚类**（FB2, Sweatgland_FB, Dermis1/2,
    Papillary_dermis）的富集？如果 F7 在所有成纤维亚类都高，则"定位"失去特异性。
Q2. 阴性对照程序 F1/F5 在 FB1 是否也高？（若也高 → 是"成纤维细胞通用"效应，非 F7 特异）
Q3. F6 在 Wound_edges 的富集，是否高于 F6 在其他上皮类（Suprabasal_epi, Basal_epi）？
    因为 Wound_edges 本身可能被标注为上皮/创缘。

输入: S26_annotation_table.csv, S26_condition_table.csv
输出: S35_selfaudit_redundancy.json / S35_selfaudit_stdout.txt
"""
import io, os, json
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    OUT.append(s)


ann = pd.read_csv(os.path.join(BASE, 'S26_annotation_table.csv'))
log('=' * 90)
log('S35 — SELF-AUDIT: is the spatial division of labour a tautology?')
log('=' * 90)
log('')
log('Annotation table: ' + str(len(ann)) + ' tissue categories')
log('')

FIBRO = ['FB1', 'FB2', 'Sweatgland_FB', 'Dermis1', 'Dermis2', 'Papillary_dermis']
EPI = ['Wound_edges', 'Suprabasal_epi', 'Basal_epi']
PROGS = ['F6', 'F7', 'F7e', 'F1', 'F5']

# ---------- Q1 / Q2: F7 在成纤维亚类间的特异性 ----------
log('--- Q1/Q2: F7 & F7e enrichment across FIBROBLAST subtypes ---')
sub = ann[ann['AnnoType'].isin(FIBRO)][['AnnoType', 'n'] + PROGS].copy()
sub = sub.sort_values('F7', ascending=False)
log('')
log('  category             n     F6      F7     F7e      F1      F5')
for _, r in sub.iterrows():
    log('  %-18s %5d  %+.3f  %+.3f  %+.3f  %+.3f  %+.3f' % (
        r['AnnoType'], r['n'], r['F6'], r['F7'], r['F7e'], r['F1'], r['F5']))
log('')

fb1 = ann[ann['AnnoType'] == 'FB1'].iloc[0]
others_fib = ann[ann['AnnoType'].isin([x for x in FIBRO if x != 'FB1'])]
log('FB1 vs mean(other fibroblast subtypes):')
q1 = {}
for p in PROGS:
    fb1_v = float(fb1[p])
    oth_v = float(others_fib[p].mean())
    delta = fb1_v - oth_v
    q1[p] = dict(fb1=fb1_v, other_fib_mean=oth_v, delta=delta)
    log('  %-4s  FB1=%+.3f  otherFib=%+.3f  delta=%+.3f' % (p, fb1_v, oth_v, delta))
log('')
log('  >> Q1 判据: 若 F7 的 delta 不显著大于 F1/F5 的 delta，则"FB1 定位"是同义反复。')

# 置换：把 FB1 的 spot 与其成分随机打散 —— 用类别级 bootstrap（按 n 加权）
rng = np.random.default_rng(20261006)
log('')
log('  -- Category-level permutation (shuffle program labels across categories, 5000x) --')
mats = ann[PROGS].values  # 17 x 5
weights = ann['n'].values.astype(float)
obs_f7 = q1['F7']['delta']
obs_f7e = q1['F7e']['delta']
obs_f1 = q1['F1']['delta']
obs_f5 = q1['F5']['delta']


def delta_for(col_idx):
    """给定程序列索引，算 FB1 - 其他成纤维类的加权均值差"""
    v = mats[:, col_idx]
    fb1_mask = (ann['AnnoType'] == 'FB1').values
    oth_mask = ann['AnnoType'].isin([x for x in FIBRO if x != 'FB1']).values
    return v[fb1_mask][0] - np.average(v[oth_mask], weights=weights[oth_mask])


null_f7, null_f7e, null_f1, null_f5 = [], [], [], []
for _ in range(5000):
    perm = rng.permutation(len(PROGS))
    null_f7.append(delta_for(perm[0]))
    null_f7e.append(delta_for(perm[1]))
    null_f1.append(delta_for(perm[2]))
    null_f5.append(delta_for(perm[3]))
null_f7 = np.array(null_f7); null_f7e = np.array(null_f7e)
null_f1 = np.array(null_f1); null_f5 = np.array(null_f5)

for nm, obs, nl in [('F7', obs_f7, null_f7), ('F7e', obs_f7e, null_f7e),
                    ('F1(neg)', obs_f1, null_f1), ('F5(neg)', obs_f5, null_f5)]:
    p = (np.sum(nl >= obs) + 1) / (len(nl) + 1)
    z = (obs - nl.mean()) / (nl.std() + 1e-12)
    log('  %-7s obs=%+.3f  null_mean=%+.3f sd=%.3f  p=%.4f  z=%+.2f' % (
        nm, obs, nl.mean(), nl.std(), p, z))

# ---------- Q3: F6 在创缘 vs 其他上皮 ----------
log('')
log('--- Q3: F6 enrichment across EPITHELIAL / edge categories ---')
ep = ann[ann['AnnoType'].isin(EPI)][['AnnoType', 'n'] + PROGS].sort_values('F6', ascending=False)
log('')
log('  category             n     F6      F7     F7e      F1      F5')
for _, r in ep.iterrows():
    log('  %-18s %5d  %+.3f  %+.3f  %+.3f  %+.3f  %+.3f' % (
        r['AnnoType'], r['n'], r['F6'], r['F7'], r['F7e'], r['F1'], r['F5']))
log('')
w = ann[ann['AnnoType'] == 'Wound_edges'].iloc[0]
oth_epi = ann[ann['AnnoType'].isin([x for x in EPI if x != 'Wound_edges'])]
q3 = {}
for p in PROGS:
    d = float(w[p]) - float(oth_epi[p].mean())
    q3[p] = dict(wound_edges=float(w[p]), other_epi_mean=float(oth_epi[p].mean()), delta=d)
    log('  %-4s  Edge=%+.3f  otherEpi=%+.3f  delta=%+.3f' % (p, w[p], oth_epi[p].mean(), d))

# ---------- 结论 ----------
log('')
log('=' * 90)
log('VERDICT')
log('=' * 90)
log('')
log('Q1 (F7 在成纤维亚类内部的特异性):')
log('   FB1 F7=%+.3f vs 其他成纤维均值 %+.3f  delta=%+.3f' % (
    q1['F7']['fb1'], q1['F7']['other_fib_mean'], q1['F7']['delta']))
if q1['F7']['delta'] > q1['F1']['delta'] and q1['F7']['delta'] > q1['F5']['delta']:
    log('   >> F7 的 FB1 delta 大于两个阴性对照 → 具有一定特异性 ✅')
else:
    log('   >> ⚠ F7 的 FB1 delta 未明显超过阴性对照 → 同义反复风险高')
log('')
log('Q2 (阴性对照在 FB1):')
log('   F1 delta=%+.3f  F5 delta=%+.3f' % (q1['F1']['delta'], q1['F5']['delta']))
log('')
log('Q3 (F6 在创缘 vs 其他上皮):')
log('   Wound_edges F6=%+.3f vs 其他上皮 %+.3f  delta=%+.3f' % (
    q3['F6']['wound_edges'], q3['F6']['other_epi_mean'], q3['F6']['delta']))
log('   F1 delta=%+.3f  F5 delta=%+.3f' % (q3['F1']['delta'], q3['F5']['delta']))

res = dict(q1_fb1_vs_otherFIB=q1, q3_edge_vs_otherEPI=q3,
           perm_q1=dict(F7=dict(obs=obs_f7, null_mean=float(null_f7.mean()),
                                null_sd=float(null_f7.std()),
                                p=float((np.sum(null_f7 >= obs_f7) + 1) / 5001)),
                        F7e=dict(obs=obs_f7e, null_mean=float(null_f7e.mean()),
                                 null_sd=float(null_f7e.std()),
                                 p=float((np.sum(null_f7e >= obs_f7e) + 1) / 5001)),
                        F1=dict(obs=obs_f1, p=float((np.sum(null_f1 >= obs_f1) + 1) / 5001)),
                        F5=dict(obs=obs_f5, p=float((np.sum(null_f5 >= obs_f5) + 1) / 5001))),
           annotation_table=ann.to_dict('records'))
with io.open(os.path.join(BASE, 'S35_selfaudit_redundancy.json'), 'w', encoding='utf-8') as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
with io.open(os.path.join(BASE, 'S35_selfaudit_stdout.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(OUT))
log('')
log('wrote S35_selfaudit_redundancy.json / S35_selfaudit_stdout.txt')
