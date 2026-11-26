#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
S31b: 空间切片身份判定 —— 用 nCount 联合指纹而非 barcode 交集

问题
----
S31 用 barcode 交集判定 h5 ↔ metadata 的对应，但 Visium barcode 集合在
不同玻片间高度重叠（同一 55µm 网格），导致覆盖率 0.72–0.94 的模糊匹配，
并产出 4 处"冲突"（其中若干个明显荒谬，如 GSM8238462 被匹配到 Donor2_Wound30）。

方法
----
用 (barcode, nCount_Spatial) 的【联合指纹】。
nCount 在单个 spot 上高度特异（整数计数），联合匹配几乎不可能是巧合。
对每个 h5：
  score(sample) = #{bc : nCount_h5[bc] == nCount_meta[bc] for bc in 交集}
取 argmax score。

同时报告：
  - 该样本自身的 spot 数 (n_meta)
  - 命中数 / 该样本总数 = 召回率
  - 命中数 / h5 总数 = 精确率
"""
import gzip, json, os
import numpy as np
import h5py

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(ROOT, '..', '10_DFU_Dataset', 'gse241124'))
H5 = os.path.join(DATA, 'h5')
META = os.path.join(DATA, 'GSE241124_spatialseq_metadata_acutewound.txt.gz')

lines = []
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); lines.append(s)

log('=' * 100)
log('S31b: 切片身份判定 —— (barcode, nCount) 联合指纹')
log('=' * 100)

# metadata
meta = {}   # sample -> {bc: nCount}
info = {}   # sample -> (donor, cond, oid)
with gzip.open(META, 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    ix = {k: i for i, k in enumerate(hdr)}
    for line in f:
        w = line.rstrip('\n').split('\t')
        if len(w) < len(hdr):
            continue
        sm = w[ix['Sample_name']]
        bc = w[ix['barcode']]
        bc = bc[len(sm) + 1:] if bc.startswith(sm + '_') else bc
        meta.setdefault(sm, {})[bc] = int(w[ix['nCount_Spatial']])
        info[sm] = (w[ix['Donor']], w[ix['Condition']], w[ix['orig.ident']])

log(f'\nmetadata: {len(meta)} 样本, {sum(len(v) for v in meta.values())} spots')

def read_h5(path):
    with h5py.File(path, 'r') as f:
        bc = [x.decode() if isinstance(x, bytes) else str(x) for x in f['matrix/barcodes'][()]]
        M = f['matrix']
        indptr = M['indptr'][()]
        data = M['data'][()]
    tot = np.add.reduceat(data, indptr[:-1]) if len(indptr) > 1 else np.array([])
    if len(indptr) > 1:
        tot = np.array([data[indptr[i]:indptr[i+1]].sum() for i in range(len(bc))])
    return dict(zip(bc, tot.astype(int)))

h5_files = sorted(f for f in os.listdir(H5) if f.endswith('.h5'))
log(f'h5: {len(h5_files)} 个\n')

log('-' * 100)
log(f'{"h5 stem":32s} {"n":>5s}  {"best match":18s} {"exact":>6s} {"prec":>6s} {"rec":>6s}  {"2nd best":18s} {"2nd":>5s}  donor/cond')
log('-' * 100)

results = []
for fn in h5_files:
    stem = fn.replace('_filtered_feature_bc_matrix.h5', '')
    hdat = read_h5(os.path.join(H5, fn))
    scored = []
    for sm, md in meta.items():
        common = set(hdat) & set(md)
        exact = sum(1 for b in common if hdat[b] == md[b])
        scored.append((exact, sm, len(common), len(md)))
    scored.sort(reverse=True)
    e1, s1, c1, n1 = scored[0]
    e2, s2, c2, n2 = scored[1]
    prec = e1 / len(hdat); rec = e1 / n1
    d, cond, oid = info[s1]
    # 模糊性：第2名是否接近
    ambiguous = ' ⚠️模糊' if e2 > 0.6 * e1 else ''
    log(f'{stem:32s} {len(hdat):5d}  {s1:18s} {e1:6d} {prec:6.3f} {rec:6.3f}  '
        f'{s2:18s} {e2:5d}  {d}/{cond}{ambiguous}')
    results.append({'h5': stem, 'n_h5': len(hdat), 'match': s1, 'exact': e1,
                    'prec': prec, 'rec': rec, 'second': s2, 'second_exact': e2,
                    'donor': d, 'cond': cond, 'oid': oid,
                    'exact_vs_second': e1 - e2})

log('\n' + '=' * 100)
log('[判定] 每个 h5 的真实身份（exact 唯一性检查）')
log('=' * 100)
log(f'{"h5 stem":32s} {"真实身份":18s} {"exact":>6s} {"次名":>6s} {"差":>5s} {"判定"}')
for r in results:
    gap = r['exact_vs_second']
    verd = '✅唯一' if gap > 0.5 * r['exact'] else '⚠️ 不唯一'
    log(f'{r["h5"]:32s} {r["match"]:18s} {r["exact"]:6d} {r["second_exact"]:6d} {gap:5d} {verd}')

log('\n' + '=' * 100)
log('[真实设计] donor × condition')
log('=' * 100)
design = {}
for r in results:
    design.setdefault(r['donor'], {})[r['cond']] = r['match']
conds = ['Skin', 'Wound1', 'Wound7', 'Wound30']
log(f'{"donor":10s} ' + ' '.join(f'{c:>16s}' for c in conds))
for dd in sorted(design):
    log(f'{dd:10s} ' + ' '.join(f'{design[dd].get(c,"—"):>16s}' for c in conds))

log('\n' + '=' * 100)
log('[GEO 标题一致性]')
GEO = {'GSM7717016': 'Donor4-Skin', 'GSM7717017': 'Donor4-Wound1',
       'GSM7717018': 'Donor4-Wound7', 'GSM7717019': 'Donor4-Wound30',
       'GSM8238462': 'Donor1-Skin', 'GSM8238463': 'Donor1-Wound1',
       'GSM8238464': 'Donor1-Wound7', 'GSM8238465': 'Donor1-Wound30',
       'GSM8238466': 'Donor2-Skin', 'GSM8238467': 'Donor2-Wound1',
       'GSM8238468': 'Donor2-Wound7', 'GSM8238469': 'Donor2-Wound30',
       'GSM8238470': 'Donor3-Skin', 'GSM8238471': 'Donor3-Wound1',
       'GSM8238472': 'Donor3-Wound7', 'GSM8238473': 'Donor3-Wound30'}
nconf = 0
for r in results:
    gsm = r['h5'].split('_')[0]
    gt = GEO.get(gsm, '?'); mt = f"{r['donor']}-{r['cond']}"
    ok = gt == mt
    if not ok: nconf += 1
    log(f'  {gsm}  GEO={gt:16s} 实测={mt:16s} {"✅" if ok else "❌冲突"}')
log(f'\n  冲突: {nconf}/16')

json.dump({'results': results, 'design': design, 'n_conflict': nconf},
          open(os.path.join(ROOT, 'S31b_identity.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
open(os.path.join(ROOT, 'S31b_stdout.txt'), 'w', encoding='utf-8').write('\n'.join(lines))
log('\n[saved] S31b_identity.json / S31b_stdout.txt')
log('\n✅ S31b 完成')
