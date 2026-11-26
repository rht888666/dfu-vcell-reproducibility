# -*- coding: utf-8 -*-
"""
S28: 跨队列证据的诚实重估 —— 哪些队列在"样本级 + 深度校正"下真正复现？

背景：
  S25 证明患者级信号可被深度混淆；S26b 证明必须用随机基因集零分布区分"全局位移"与"程序特异"。
  S21 的 GSE223964 细胞级"全部程序 p=0.0002"是全局位移的典型特征（连 F2/F3 都显著）。

本脚本汇总所有跨队列证据，统一到【样本级 + 深度校正 + 阴性对照】的判据下：
  判据（三重）：
    C1. 样本级 p < 0.05
    C2. 阴性对照 F1/F5 不显著（否则全局伪影）
    C3. 目标程序 F6/F7 显著

  只有 C1∧C2∧C3 全部满足，才判定"该队列复现"。

汇总对象：
  - GSE165816 (主队列)：S25
  - GSE231643 (S12)：已含 delta_res/p_res
  - GSE223964 (S21-A)
  - GSE245703 (S21-B)：已知被排除（F1 也显著）
  - GSE241132 (S20)：纵向（不同设计，单列）
  - GSE241124 (S26)：空间（不同设计，单列）
"""
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
ROOT = Path(r"C:/Users/Administrator/Desktop/虚拟细胞论著投稿")
OUT = ROOT / "17_DFU_VirtualCell_10plus"

LOG = []
def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOG.append(s)

log('=' * 88)
log('S28: 跨队列证据的诚实重估（样本级 + 深度校正 + 阴性对照）')
log('=' * 88)

rows = []

# ---- GSE165816 (主队列, S25) ----
p25 = json.load(open(OUT/"S25_depth_corrected.json", encoding='utf-8'))
r2 = p25['R2']
def get_p(p, key='resid_p'):
    return r2.get(p, {}).get(key, np.nan)
rows.append(dict(cohort='GSE165816', design='横断面 DFU (主队列)',
                 F6_p=get_p('F6_orig'), F7_p=get_p('F7_orig'), F7e_p=get_p('F7_ecm'),
                 F1_p=get_p('F1_orig'), F5_p=get_p('F5_orig'),
                 level='patient', depth='residualized', note='S25'))

# ---- GSE231643 (S12) ----
s12 = json.load(open(OUT/"S12_crosscohort_GSE231643.json", encoding='utf-8'))
def g12(p):
    v = s12.get(p, {})
    return v.get('p_res', v.get('p_raw', np.nan))
rows.append(dict(cohort='GSE231643', design='DFU 清创 (第2队列)',
                 F6_p=g12('F6'), F7_p=g12('F7'), F7e_p=g12('F7e'),
                 F1_p=g12('F1'), F5_p=g12('F5'),
                 level='patient', depth='residualized', note='S12'))

# ---- GSE223964 (S21-A) ----
s21 = json.load(open(OUT/"S21_crosscohort.json", encoding='utf-8'))
g = s21['GSE223964']
rs = g.get('res_sample', {})
def g21sample(p): return rs.get(p, {}).get('p_perm', np.nan)
rc = g.get('res_cell', {})
def g21cell(p): return rc.get(p, {}).get('p_perm', np.nan)
rows.append(dict(cohort='GSE223964', design='糖尿病/非糖尿病皮肤 (第3队列)',
                 F6_p=g21sample('F6'), F7_p=g21sample('F7'), F7e_p=g21sample('F7e'),
                 F1_p=g21sample('F1'), F5_p=g21sample('F5'),
                 level='sample', depth='residualized', note='S21-A'))
rows.append(dict(cohort='GSE223964 (cell)', design='同上, 细胞级(参考)',
                 F6_p=g21cell('F6'), F7_p=g21cell('F7'), F7e_p=g21cell('F7e'),
                 F1_p=g21cell('F1'), F5_p=g21cell('F5'),
                 level='cell', depth='residualized', note='S21-A 细胞级'))

# GSE245703 —— 从 S21_VERDICT 已知被排除
rows.append(dict(cohort='GSE245703', design='糖尿病创面 (第4队列)',
                 F6_p=np.nan, F7_p=np.nan, F7e_p=np.nan,
                 F1_p=0.0002, F5_p=np.nan,
                 level='cell', depth='residualized',
                 note='S21-B 已排除：F1 阴性对照也显著(p=0.0002) → 全局伪影'))

DF = pd.DataFrame(rows)
log('\n' + '=' * 88)
log('一、跨队列证据汇总表')
log('=' * 88)
log(f'\n  {"队列":22s} {"层级":8s} {"F6 p":>8s} {"F7 p":>8s} {"F7e p":>8s} {"F1 p":>8s} {"F5 p":>8s}')
for _, r in DF.iterrows():
    def fmt(x): return f'{x:.4f}' if pd.notna(x) else '  n/a '
    log(f'  {r["cohort"]:22s} {r["level"]:8s} {fmt(r["F6_p"]):>8s} {fmt(r["F7_p"]):>8s} {fmt(r["F7e_p"]):>8s} '
        f'{fmt(r["F1_p"]):>8s} {fmt(r["F5_p"]):>8s}')

# ---- 三重判据 ----
log('\n' + '=' * 88)
log('二、三重判据裁决（C1 样本级显著 ∧ C2 阴性对照不显著 ∧ C3 F6/F7 显著）')
log('=' * 88)
log('\n  判据: C1 F7 或 F7e p<0.05 | C2 F1 且 F5 p>0.05 | C3 F6 或 F7 p<0.05')
log('')
verdicts = []
for _, r in DF.iterrows():
    f1 = r['F1_p']; f5 = r['F5_p']
    c1 = (pd.notna(r['F7_p']) and r['F7_p']<0.05) or (pd.notna(r['F7e_p']) and r['F7e_p']<0.05)
    c2 = (pd.isna(f1) or f1>0.05) and (pd.isna(f5) or f5>0.05)
    c3 = (pd.notna(r['F6_p']) and r['F6_p']<0.05) or (pd.notna(r['F7_p']) and r['F7_p']<0.05)
    ok = c1 and c2 and c3
    tag = '✅ 复现' if ok else ('❌ 排除' if not c2 else '🟡 部分')
    if pd.isna(r['F7_p']) and pd.isna(r['F7e_p']):
        tag = '❌ 已排除' if (pd.notna(f1) and f1<0.05) else '⚪ 未评'
    log(f'  {r["cohort"]:22s} C1={c1!s:5s} C2={c2!s:5s} C3={c3!s:5s} -> {tag}   ({r["note"]})')
    verdicts.append(dict(cohort=r['cohort'], c1=bool(c1), c2=bool(c2), c3=bool(c3), verdict=tag))

# ---- 结论 ----
log('\n' + '=' * 88)
log('三、诚实结论')
log('=' * 88)
ok_list = [v['cohort'] for v in verdicts if v['verdict']=='✅ 复现']
part_list = [v['cohort'] for v in verdicts if v['verdict']=='🟡 部分']
fail_list = [v['cohort'] for v in verdicts if '排除' in v['verdict']]
log(f'\n  ✅ 通过三重判据: {ok_list}')
log(f'  🟡 部分（阴性对照可，但目标程序不显著）: {part_list}')
log(f'  ❌ 排除: {fail_list}')

log('\n  【关键诚实点】')
log('  - GSE223964 细胞级"全部程序 p=0.0002"是【全局位移】特征（连 F2/F3 都显著）')
log('  - 降级到样本级后，仅 F7 显著 (p=0.047)，F6 不显著 (p=0.268)')
log('  - → GSE223964 是"部分复现"：只能支持 F7，不支持 F6')
log('  - GSE245703 因 F1 阴性对照也显著而被排除')

out = dict(table=DF.fillna('NA').to_dict('records'), verdicts=verdicts,
           passed=ok_list, partial=part_list, excluded=fail_list)
with open(OUT/"S28_crosscohort_honest.json",'w',encoding='utf-8') as f:
    json.dump(out,f,indent=2,ensure_ascii=False,default=str)
log(f'\n[saved] S28_crosscohort_honest.json')
with open(OUT/"S28_run.log",'w',encoding='utf-8') as f:
    f.write('\n'.join(LOG))
log('\n✅ S28 完成')
