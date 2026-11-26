# -*- coding: utf-8 -*-
"""
S32b - Manuscript number audit v2 (R15).
不再硬编码"稿件写的值"，而是从 MANUSCRIPT_v3_integrated.md 里用正则抽取数字，
再与 JSON 产物比对。这样稿件每次修改后重跑都能自动回归。

用法: python S32b_manuscript_audit.py
"""
import json, os, re

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    OUT.append(s)


def J(name):
    with open(os.path.join(BASE, name), 'r', encoding='utf-8') as f:
        return json.load(f)


MS = open(os.path.join(BASE, 'MANUSCRIPT_v3_integrated.md'), 'r', encoding='utf-8').read()
R = []


def num(s):
    """把稿件里的数字字符串转为 float（处理 U+2212 减号）"""
    return float(s.replace('\u2212', '-').replace(',', ''))


def grab(pattern, group=1, flags=0):
    m = re.search(pattern, MS, flags)
    return num(m.group(group)) if m else None


def check(tag, ms, js, tol=1e-3, note=''):
    ok = (ms is not None) and abs(float(ms) - float(js)) <= tol
    log(('OK  ' if ok else 'BAD ') + '[' + tag + '] ms=' + str(ms) + '  json=' + str(js) +
        ('  ' + note if note else ''))
    R.append(dict(tag=tag, manuscript=str(ms), json=str(js), ok=bool(ok), note=note))


log('=' * 90)
log('S32b - MANUSCRIPT NUMBER AUDIT v2 (regex-extracted from .md)')
log('=' * 90)

# ---------------- S26b spatial
G = J('S26b_global_gradient.json')
ep = G['empirical_p']
log('')
log('--- S26b: spatial ---')
check('S2.3 spots 22,915', grab(r'(\d{2},\d{3}) spots'), G['n_spots'], 0)
check('S2.3 4 donors x 4 tp', grab(r'(\d) donors . (\d) timepoints, 16 sections', 1), 4, 0)
check('S2.3 16 sections', grab(r'(\d+) sections'), 16, 0)
_f7e = grab(r'donor-mean .{1,3} = \*\*\+([0-9]+\.[0-9]+)\*\*')
check('F7e temporal rho 0.470', _f7e, ep['F7e']['obs'], 5e-3)
check('F7e temporal p 0.002', grab(r'F7e accumulation exceeds a random-gene null across time.{0,40}?p.{0,3}=\s*0\.0([0-9]+)'), ep['F7e']['p_emp'], 5e-4) if False else None
# 直接抽 §2.3 段
check('F7e p (S2.3)', grab(r'donor-mean .{1,3} = \*\*\+0\.470\*\*, empirical \*\*p = ([0-9.]+)\*\*'), ep['F7e']['p_emp'], 5e-4)
check('F7 rho (S2.3)', grab(r'broader F7 programme is directionally identical \(4/4 donors, .{1,3} = \+([0-9.]+)\)'), ep['F7']['obs'], 5e-4)
check('F7 p (S2.3)', grab(r'falls within the null distribution \(p = ([0-9.]+)\)'), ep['F7']['p_emp'], 5e-4)
check('F1 p 0.687', grab(r'F1 and F5 behave correctly, remaining inside the null \(p = ([0-9.]+)'), ep['F1']['p_emp'], 5e-4)
check('F5 p 0.976', grab(r'remaining inside the null \(p = [0-9.]+, ([0-9.]+)\)'), ep['F5']['p_emp'], 5e-4)
check('null sd 0.125', grab(r'the null.s standard deviation is ([0-9.]+)'), G['null_rand']['sd'], 5e-4)
fb = G['fb1_enrichment']
check('FB1 delta 0.703', grab(r'mean F7 in FB1 spots = \+0\.652 versus .0\.051 elsewhere, . = \*\*\+([0-9.]+)\*\*'), fb['obs'], 5e-4)
check('FB1 p_perm 0.0005', grab(r'spot-level permutation \*\*p = ([0-9.]+)\*\*\)'), fb['p_perm'], 5e-4)
check('FB1 vs-random 0.025', grab(r'exceeds a random-gene-set null as well \(\*\*p = ([0-9.]+)\*\*\)'), fb['p_specific'], 5e-4)
check('Abstract F6 -0.280', grab(r'F6 edge gradient is .{1,3} = .([0-9]+\.?[0-9]*)'), abs(sum(G['program_rhos']['F6']) / 4), 5e-3)

# ---------------- S35 subtype-level (R15 第12次自我推翻核心)
A = J('S35_selfaudit_redundancy.json')
pm, q1 = A['perm_q1'], A['q1_fb1_vs_otherFIB']
log('')
log('--- S35: F7 vs F7e subtype-level localization ---')
check('S2.3 F7e subtype delta 0.895', grab(r'other fibroblast subtypes = \+([0-9]+\.[0-9]+)'), q1['F7e']['delta'], 5e-3)
check('S2.3 F7e subtype p 0.0002', grab(r'category-level (?:permutation )?\*{0,2}p\*{0,2} = ([0-9]+\.[0-9]+)'), pm['F7e']['p'], 1e-4)
check('S2.3 F7 subtype delta 0.758', grab(r'whereas \*\*F7 does not\*\* \(. = \+([0-9]+\.[0-9]+)'), q1['F7']['delta'], 5e-3)
check('S2.3 F7 subtype p 0.409', grab(r'\*\*p = ([0-9]+\.[0-9]+)\*\*\)\. F1 and F5 behave'), pm['F7']['p'], 5e-4)
_sm = [r for r in A['annotation_table'] if r['AnnoType'] == 'Smooth_muscle'][0]
check('S2.3 Smooth_muscle F7 0.510', grab(r'Smooth_muscle. category carries the second-highest F7 score of all 17 categories \(\+([0-9]+\.[0-9]+)\)'), _sm['F7'], 5e-4)
# F7e donor-mean rho +0.483 (Table 1)
t1 = re.search(r'\|\s*\*\*F7e\*\*\s*\|[^|]*\|[^|]*\|[^|]*\|\s*\*\*\+([0-9]+\.[0-9]+)\*\*', MS)
L1 = J('S20_trajectory.json')
_f7e_dm = sum(L1['trend']['F7e']['donor_rhos']) / len(L1['trend']['F7e']['donor_rhos'])
check('T1 F7e donor-mean 0.483', num(t1.group(1)) if t1 else None, _f7e_dm, 5e-4)

# ---------------- R15b: monotonicity guard (第13次自我推翻)
# 稿件曾声称 F6 "declines monotonically"；逐时点核查证伪。此处把真实时点序列锁进审计。
T = J('S20b_trajectory_tests.json')['timecourse']
tc_D0_F6 = grab(r'F6 \(pooled spot mean\) is ([+\-\u2212][0-9]+\.[0-9]+) at D0')
tc_D1_F6 = grab(r'rises to \+([0-9]+\.[0-9]+) at D1')
tc_D30_F6 = grab(r'and \u2212([0-9]+\.[0-9]+) \(D30\)')
tc_D30_F7 = grab(r'F7 moves [^;]*\u2192 \+([0-9]+\.[0-9]+)\.')
log('')
log('--- R15b: monotonicity guard (sampled timecourse) ---')
check('tc F6 D0 +0.093', tc_D0_F6, T['sc_F6']['Skin'], 5e-4)
check('tc F6 D1 +0.171', tc_D1_F6, T['sc_F6']['Wound1'], 5e-4)
check('tc F6 D30 -0.156', -tc_D30_F6 if tc_D30_F6 is not None else None, T['sc_F6']['Wound30'], 5e-4)
check('tc F7 D30 +0.261', tc_D30_F7, T['sc_F7']['Wound30'], 5e-4)
# 硬断言：真值序列本身不是单调的（若将来某版改动使序列变单调，此处会提示复核措辞）
_f6 = [T['sc_F6'][k] for k in ['Skin', 'Wound1', 'Wound7', 'Wound30']]
_f7 = [T['sc_F7'][k] for k in ['Skin', 'Wound1', 'Wound7', 'Wound30']]
_mono = all(_f6[i] >= _f6[i + 1] for i in range(3)) or all(_f7[i] <= _f7[i + 1] for i in range(3))
log(('OK  ' if not _mono else 'BAD ') + f'[monotonicity] F6={_f6} F7={_f7} globally_monotone={_mono} '
    f'({"must review wording" if _mono else "claim of monotonicity correctly withheld"})')
R.append(dict(tag='monotonicity', manuscript=str(_mono), json='not monotone', ok=bool(not _mono)))

# ---------------- S29b edge gradient
S = J('S29b_gradient_specificity.json')
pg, dd = S['program_gradients'], S['double_diff']
log('')
log('--- S29b: edge gradient (Table 3, D7 convention) ---')
# Table 3 行: | **F6** | **−0.160** | **0.00019** | ✅ exceeds null (z = −3.73) |
t3 = re.search(r'\|\s*\*\*F6\*\*\s*\|\s*\*\*\u2212([0-9.]+)\*\*\s*\|\s*\*\*([0-9.]+)\*\*\s*\|\s*\u2705 exceeds null \(z = \u2212([0-9.]+)\)', MS)
if t3:
    check('T3 F6 rho -0.160', -num(t3.group(1)), pg['Wound7_F6']['rho'], 5e-4)
    check('T3 F6 p 0.00019', num(t3.group(2)), pg['Wound7_F6']['p_emp'], 5e-5)
    check('T3 F6 z -3.73', -num(t3.group(3)), pg['Wound7_F6']['z'], 0.02)
else:
    check('T3 F6 row regex', None, 0, 0, 'regex miss')
for lbl, key in [('F7', 'Wound7_F7'), ('F7e', 'Wound7_F7e'), ('F1', 'Wound7_F1'), ('F5', 'Wound7_F5')]:
    m = re.search(r'\|\s*' + lbl + r'(?: \(neg\.\))?\s*\|\s*\+([0-9.]+)\s*\|\s*([0-9.]+)\s*\|', MS)
    if m:
        check('T3 ' + lbl + ' rho', num(m.group(1)), pg[key]['rho'], 5e-3)
        check('T3 ' + lbl + ' p', num(m.group(2)), pg[key]['p_emp'], 5e-3)
    else:
        check('T3 ' + lbl + ' regex', None, 0, 0, 'regex miss')
check('T3 F6-F5 @D7 rho', grab(r'F6\u2212F5 at D7: . = .([0-9.]+)'), abs(dd['F6 \u2212 F5']['Wound7']['rho']), 5e-4)
check('T3 F6-F5 @D7 p', grab(r'F6\u2212F5 at D7: .{0,20}?\*\*p = ([0-9.]+)\*\*'), dd['F6 \u2212 F5']['Wound7']['p'], 5e-4)
check('T3 F6-F1 @D7 rho', grab(r'F6\u2212F1 at D7: . = .([0-9.]+)'), abs(dd['F6 \u2212 F1']['Wound7']['rho']), 5e-4)
check('T3 F6-F1 @D7 p', grab(r'F6\u2212F1 at D7: .{0,20}?\*\*p = ([0-9.]+)\*\*'), dd['F6 \u2212 F1']['Wound7']['p'], 5e-4)

# ---------------- S20 longitudinal Table 1
T = J('S20_trajectory.json')
tr = T['trend']
log('')
log('--- S20: Table 1 longitudinal ---')
rows = re.findall(r'\|\s*(F6i|F7e|F7|F6|F1 \(neg\.\)|F5 \(neg\.\))\s*\|\s*([\u2212+]?[0-9.]+)\s*\|\s*([\u2212+]?[0-9.]+)\s*\|\s*([\u2212+]?[0-9.]+)\s*\|\s*\*?([\u2212+]?[0-9.]+)', MS)
for prog, a, b, c, dm in rows:
    key = prog.split(' ')[0]
    if key not in tr:
        continue
    dr = tr[key]['donor_rhos']
    check('T1 ' + key + ' PWH26', num(a), dr[0], 5e-4)
    check('T1 ' + key + ' PWH27', num(b), dr[1], 5e-4)
    check('T1 ' + key + ' PWH28', num(c), dr[2], 5e-4)
    check('T1 ' + key + ' mean', num(dm), sum(dr) / 3, 1e-3)
log('   (extracted ' + str(len(rows)) + ' rows)')
check('S2.2 F7 per-donor line', grab(r'F7 . = \+0\.25/\+0\.54/\+0\.55'), None, 1, 1, 'manual check') if False else None

# ---------------- S24b geometry
M = J('S24b_state_model.json')
log('')
log('--- S24b: geometry ---')
check('S2.2 R_bar 0.991', grab(r'circular concordance \*\*.{1,2} = ([0-9.]+)\*\*'), M['angle_R_bar'], 5e-4)
mj = M['dfu_projection']['per_group']
tbl_rows = re.findall(r'\|\s*(Healthy|Diabetic, no DFU|\*\*Healer\*\*|\*\*NonHealer\*\*)\s*\|\s*(\d+)\s*\|\s*\*?\*?([0-9.]+)', MS)
name_map = {'Healthy': 'Healthy', 'Diabetic, no DFU': 'DiabeticNoDFU',
            '**Healer**': 'Healer', '**NonHealer**': 'NonHealer'}
for lbl, n, p in tbl_rows:
    g = name_map[lbl]
    check('S2.5 ' + g + ' n', num(n), mj[g]['n'], 0)
    check('S2.5 ' + g + ' progress', num(p), mj[g]['progress_mean'], 5e-4)

# ---------------- S25b scheme matrix
W = J('S25b_scheme_comparison.json')['results']
log('')
log('--- S25b: Table 5 ---')
t5 = {
    'Raw score': 'raw', 'Depth-matched pairs': 'd_matched', 'Ratio F7/F6': 'e_ratio_F7F6',
    'Linear residualization': 'a_lin_resid', 'Rank residualization': 'b_rank_resid',
    'Stratified z (by depth tercile)': 'c_stratum_z'}
for label, key in t5.items():
    m = re.search(r'\|\s*\*{0,2}' + re.escape(label) + r'\*{0,2}\s*\|\s*\*{0,2}\+([0-9.]+)\*{0,2}\s*\|\s*\*{0,2}([0-9.]+)\*{0,2}\s*\|\s*\*{0,2}([\u2212+\-]?[0-9.]+|\u2014)\*{0,2}', MS)
    if m:
        d = W[key]
        check('T5 ' + key + ' obs', num(m.group(1)), d['obs'], 1.5e-3)
        check('T5 ' + key + ' p', num(m.group(2)), d['p'], 1.5e-3)
        if m.group(3) != '\u2014':  # em dash = 该方案无 rho
            _r = m.group(3).replace('\u2212', '-')
            check('T5 ' + key + ' rho', num(_r), d.get('umi_rho', 0.0), 1.5e-3)
    else:
        check('T5 ' + key + ' regex', None, 0, 0, 'regex miss')

# ---------------- S28 cross-cohort
C = J('S28_crosscohort_honest.json')
tbl = {r['cohort']: r for r in C['table']}
log('')
log('--- S28: Table 4 ---')
for coh in ['GSE165816', 'GSE231643', 'GSE223964']:
    m = re.search(r'\|\s*\*{0,2}' + coh + r'\*{0,2}\s*\|[^|]+\|[^|]+\|\s*\*{0,2}([0-9.]+)\*{0,2}\s*\|\s*\*{0,2}([0-9.]+)\*{0,2}\s*\|\s*\*{0,2}([0-9.]+)', MS)
    if m:
        r = tbl[coh]
        check('T4 ' + coh + ' F6', num(m.group(1)), r['F6_p'], 1.5e-3)
        check('T4 ' + coh + ' F7', num(m.group(2)), r['F7_p'], 1.5e-3)
        check('T4 ' + coh + ' F7e', num(m.group(3)), r['F7e_p'], 1.5e-3)
    else:
        check('T4 ' + coh + ' regex', None, 0, 0, 'regex miss')

# ---------------- 汇总
log('')
log('=' * 90)
n_ok = sum(1 for r in R if r['ok'])
n_bad = len(R) - n_ok
log('TOTAL ' + str(len(R)) + ' | OK ' + str(n_ok) + ' | BAD ' + str(n_bad))
if n_bad:
    log('')
    log('BAD LIST:')
    for r in R:
        if not r['ok']:
            log('   [' + r['tag'] + '] ms=' + r['manuscript'] + ' json=' + r['json'] + ' ' + r['note'])
with open(os.path.join(BASE, 'S32b_audit_report.json'), 'w', encoding='utf-8') as f:
    json.dump(dict(results=R, n_total=len(R), n_ok=n_ok, n_bad=n_bad), f, ensure_ascii=False, indent=2)
with open(os.path.join(BASE, 'S32b_audit_stdout.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(OUT))
log('')
log('wrote S32b_audit_report.json / S32b_audit_stdout.txt')
