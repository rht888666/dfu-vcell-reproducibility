# -*- coding: utf-8 -*-
"""
S32 - Manuscript number audit (R15).
核验 MANUSCRIPT_v3_integrated.md 中每一个数字是否与 JSON/CSV 产物一致。
硬约束：稿件中任何未被脚本生成的数字都必须被揪出来。

用法: python S32_manuscript_number_audit.py
"""
import json, os

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    OUT.append(s)


def J(name):
    with open(os.path.join(BASE, name), 'r', encoding='utf-8') as f:
        return json.load(f)


def check(tag, ms, js, tol=1e-3, note=''):
    try:
        ok = abs(float(ms) - float(js)) <= tol
    except Exception:
        ok = (str(ms) == str(js))
    log(('OK  ' if ok else 'BAD ') + '[' + tag + '] ms=' + str(ms) + '  json=' + str(js) +
        ('  ' + note if note else ''))
    return dict(tag=tag, manuscript=str(ms), json=str(js), ok=bool(ok), note=note)


R = []
log('=' * 90)
log('S32 - MANUSCRIPT v3 NUMBER AUDIT (R15)')
log('=' * 90)

# ============================================== S26b 空间
G = J('S26b_global_gradient.json')
log('')
log('--- S26b_global_gradient.json : spatial 4 donors / 16 sections ---')
log('n_spots=' + str(G['n_spots']) + '  donors=' + str(G['donors']) + '  n_sections=' + str(G['n_sections']))
R.append(check('S2.3 spots 22915', 22915, G['n_spots'], 0))
R.append(check('S2.3 donors 4', 4, len(G['donors']), 0))
R.append(check('S2.3 sections 16', 16, G['n_sections'], 0))

ep = G['empirical_p']
R.append(check('F7e temporal rho 0.470', 0.470, ep['F7e']['obs'], 5e-4))
R.append(check('F7e temporal p 0.002', 0.002, ep['F7e']['p_emp'], 5e-4))
R.append(check('F7 temporal rho 0.361', 0.361, ep['F7']['obs'], 5e-4))
R.append(check('F7 temporal p 0.238', 0.238, ep['F7']['p_emp'], 5e-4))
R.append(check('F1 temporal p 0.687', 0.687, ep['F1']['p_emp'], 5e-4))
R.append(check('F5 temporal p 0.976', 0.976, ep['F5']['p_emp'], 5e-4))
R.append(check('null sd 0.125', 0.125, G['null_rand']['sd'], 5e-4))

fb = G['fb1_enrichment']
R.append(check('FB1 delta 0.703', 0.703, fb['obs'], 5e-4))
R.append(check('FB1 p_perm 0.0005', 0.0005, fb['p_perm'], 5e-4))
R.append(check('FB1 vs-random p 0.025', 0.025, fb['p_specific'], 5e-4))

log('')
log('F6 donor-level edge rhos:')
_f6 = G['program_rhos']['F6']
for i, r in enumerate(_f6):
    log('   Donor' + str(i + 1) + ' = ' + format(r, '+.4f'))
log('   mean = ' + format(sum(_f6) / len(_f6), '+.4f'))
R.append(check('Abstract F6 rho -0.28 (donor mean)', -0.28, sum(_f6) / len(_f6), 5e-3))

# ============================================== S29b 创缘
S = J('S29b_gradient_specificity.json')
log('')
log('--- S29b_gradient_specificity.json : wound-edge distance gradient ---')
pg, dd = S['program_gradients'], S['double_diff']

log('F6 gradient by timepoint:')
for k in pg:
    if k.endswith('_F6'):
        v = pg[k]
        log('   ' + k + ': rho=' + format(v['rho'], '+.4f') + ' p=' + format(v['p_emp'], '.6f') +
            ' z=' + format(v['z'], '.3f'))

# 注意：Table 3 的 F6 行实际对应 Donor1(Skin 切片名) 的 rho = -0.1771，
# 该值与 program_rhos['F6'][0] 相同；而 S29b 按时点重组后 Skin 时点 rho=-0.1198。
# 稿件把两个口径混用了 —— 记为待修正项。
R.append(check('T3 F6 rho -0.177 vs donor1', -0.177, G['program_rhos']['F6'][0], 5e-4,
               'matches Donor1 slice group; NOT the S29b "Skin" timepoint (-0.1198)'))
R.append(check('T3 F6 rho -0.177 vs S29b Skin tp', -0.177, pg['Skin_F6']['rho'], 5e-4,
               'MISMATCH - manuscript mixes two conventions'))
R.append(check('T3 F6 p 0.0002 vs S29b Skin tp', 0.0002, pg['Skin_F6']['p_emp'], 5e-4,
               'Skin tp p=0.000284'))
R.append(check('T3 F6 z -3.7 vs S29b Skin tp z', -3.7, pg['Skin_F6']['z'], 0.05,
               'Skin z=-3.63 ; Wound7 z=-3.73'))

log('')
log('Table 3 rows - which timepoint does the manuscript use?')
log('   claimed: F7 +0.332/p0.058 ; F7e +0.229/p0.805 ; F1 +0.197/p0.423 ; F5 +0.270/p0.060')
for lbl in ['F7', 'F7e', 'F1', 'F5']:
    k = 'Skin_' + lbl
    log('   Skin_' + lbl + ': rho=' + format(pg[k]['rho'], '+.4f') + ' p=' + format(pg[k]['p_emp'], '.4f'))
for lbl in ['F7', 'F7e', 'F1', 'F5']:
    k = 'Wound7_' + lbl
    log('   Wound7_' + lbl + ': rho=' + format(pg[k]['rho'], '+.4f') + ' p=' + format(pg[k]['p_emp'], '.4f'))

log('')
R.append(check('S2.3 F6-F5 @D7 rho -0.379', -0.379, dd['F6 \u2212 F5']['Wound7']['rho'], 5e-4))
R.append(check('S2.3 F6-F5 @D7 p 0.005', 0.005, dd['F6 \u2212 F5']['Wound7']['p'], 5e-4))
R.append(check('S2.3 F6-F1 @D7 rho -0.305', -0.305, dd['F6 \u2212 F1']['Wound7']['rho'], 5e-4))
R.append(check('S2.3 F6-F1 @D7 p 0.029', 0.029, dd['F6 \u2212 F1']['Wound7']['p'], 5e-4))

# Table 3 其余行：按 S29b 的时点口径重新核验（稿件的值来自哪个时点？）
log('')
log('Table 3 residual rows vs S29b Wound7 convention:')
R.append(check('T3 F7 rho +0.332', 0.332, pg['Wound7_F7']['rho'], 5e-3, 'Wound7_F7=+0.351'))
R.append(check('T3 F7 p 0.058', 0.058, pg['Wound7_F7']['p_emp'], 5e-3, 'Wound7_F7 p=0.0583'))
R.append(check('T3 F7e rho +0.229', 0.229, pg['Wound7_F7e']['rho'], 5e-3, 'Wound7_F7e=+0.157'))
R.append(check('T3 F7e p 0.805', 0.805, pg['Wound7_F7e']['p_emp'], 5e-3, 'Wound7_F7e p=0.8046'))
R.append(check('T3 F1 rho +0.197', 0.197, pg['Wound7_F1']['rho'], 5e-3, 'Wound7_F1=+0.252'))
R.append(check('T3 F1 p 0.423', 0.423, pg['Wound7_F1']['p_emp'], 5e-3, 'Wound7_F1 p=0.4227'))
R.append(check('T3 F5 rho +0.270', 0.270, pg['Wound7_F5']['rho'], 5e-3, 'Wound7_F5=+0.350'))
R.append(check('T3 F5 p 0.060', 0.060, pg['Wound7_F5']['p_emp'], 5e-3, 'Wound7_F5 p=0.0603'))

# ============================================== S20 纵向
T = J('S20_trajectory.json')
log('')
log('--- S20_trajectory.json : longitudinal GSE241132 ---')
tr = T['trend']
for prog in ['F6', 'F6i', 'F7', 'F7e', 'F1', 'F5']:
    dr = tr[prog]['donor_rhos']
    log('   ' + prog.ljust(4) + ' donors=' + str([format(x, '+.3f') for x in dr]) +
        ' mean=' + format(sum(dr) / len(dr), '+.3f'))

R.append(check('T1 F6 PWH26 -0.430', -0.430, tr['F6']['donor_rhos'][0], 5e-4))
R.append(check('T1 F6 PWH27 -0.435', -0.435, tr['F6']['donor_rhos'][1], 5e-4))
R.append(check('T1 F6 PWH28 -0.475', -0.475, tr['F6']['donor_rhos'][2], 5e-4))
R.append(check('T1 F6 mean -0.447', -0.447, sum(tr['F6']['donor_rhos']) / 3, 5e-4))
R.append(check('T1 F6i mean -0.500', -0.500, sum(tr['F6i']['donor_rhos']) / 3, 5e-4))
R.append(check('T1 F7 mean +0.486', 0.486, sum(tr['F7']['donor_rhos']) / 3, 5e-4))
R.append(check('T1 F7e mean +0.521', 0.521, sum(tr['F7e']['donor_rhos']) / 3, 5e-4))

log('')
log('S2.2 claims F7 per-donor = +0.42/+0.53/+0.51')
log('   JSON F7 = ' + str([format(x, '+.3f') for x in tr['F7']['donor_rhos']]))
R.append(check('S2.2 F7 PWH26 +0.42', 0.42, tr['F7']['donor_rhos'][0], 5e-3, 'JSON=+0.248'))
R.append(check('S2.2 F7 PWH27 +0.53', 0.53, tr['F7']['donor_rhos'][1], 5e-3))
R.append(check('S2.2 F7 PWH28 +0.51', 0.51, tr['F7']['donor_rhos'][2], 5e-3))
log('S2.2 claims F6 per-donor = -0.43/-0.43/-0.47')
log('   JSON F6 = ' + str([format(x, '+.3f') for x in tr['F6']['donor_rhos']]))

# ============================================== S24b 几何
M = J('S24b_state_model.json')
log('')
log('--- S24b_state_model.json : geometry ---')
log('   angle_R_bar = ' + format(M['angle_R_bar'], '.4f'))
log('   normal_angle = ' + format(M['normal_angle'], '.2f'))
log('   donor_angles = ' + str([format(a, '.2f') for a in M['donor_angles']]))
log('   spread = ' + format(max(M['donor_angles']) - min(M['donor_angles']), '.2f'))
R.append(check('S2.2 R_bar 0.991', 0.991, M['angle_R_bar'], 5e-4))
R.append(check('S2.2 mean angle +119.1', 119.1, M['normal_angle'], 0.5, 'JSON=117.76'))

pj = M['dfu_projection']['per_group']
log('')
log('--- S24b dfu_projection (S2.5 table) ---')
for g, wn, wp in [('Healthy', 11, 0.282), ('DiabeticNoDFU', 6, 0.590),
                  ('Healer', 9, 0.640), ('NonHealer', 5, 0.299)]:
    got = pj[g]
    log('   ' + g + ': n=' + str(got['n']) + ' progress=' + format(got['progress_mean'], '.4f'))
    R.append(check('S2.5 ' + g + ' n', wn, got['n'], 0))
    R.append(check('S2.5 ' + g + ' progress', wp, got['progress_mean'], 5e-4))
hv = M['dfu_projection']['healer_vs_nonhealer']
R.append(check('S2.5 Healer-NonHealer obs 0.341', 0.341, hv['obs'], 5e-4))
R.append(check('S2.5 raw p 0.026', 0.026, hv['p_perm'], 5e-4))

# ============================================== S25b 五方案
W = J('S25b_scheme_comparison.json')['results']
log('')
log('--- S25b_scheme_comparison.json (Table 5) ---')
want5 = [('raw', 0.341, 0.026, 0.613), ('d_matched', 0.309, 0.025, None),
         ('e_ratio_F7F6', 0.178, 0.084, 0.504), ('a_lin_resid', 0.174, 0.143, -0.392),
         ('b_rank_resid', 0.336, 0.108, 0.000), ('c_stratum_z', 0.180, 0.193, 0.110)]
for key, obs, p, rho in want5:
    d = W[key]
    _ur = d.get('umi_rho', float('nan'))
    log('   ' + key + ': obs=' + format(d['obs'], '.4f') + ' p=' + format(d['p'], '.4f') +
        ' umi_rho=' + (format(_ur, '.4f') if _ur == _ur else 'n/a'))
    R.append(check('T5 ' + key + ' obs', obs, d['obs'], 1.5e-3))
    R.append(check('T5 ' + key + ' p', p, d['p'], 1.5e-3))
    if rho is not None and _ur == _ur:
        R.append(check('T5 ' + key + ' rho', rho, _ur, 1.5e-3))

# ============================================== S28 跨队列
C = J('S28_crosscohort_honest.json')
log('')
log('--- S28_crosscohort_honest.json (Table 4) ---')
tbl = {r['cohort']: r for r in C['table']}
for coh, f6, f7, f7e in [('GSE165816', 0.128, 0.056, 0.230),
                         ('GSE231643', 0.037, 0.208, 0.033),
                         ('GSE223964', 0.268, 0.047, 0.058)]:
    r = tbl[coh]
    log('   ' + coh + ': F6=' + format(r['F6_p'], '.4f') + ' F7=' + format(r['F7_p'], '.4f') +
        ' F7e=' + format(r['F7e_p'], '.4f'))
    R.append(check('T4 ' + coh + ' F6', f6, r['F6_p'], 1.5e-3))
    R.append(check('T4 ' + coh + ' F7', f7, r['F7_p'], 1.5e-3))
    R.append(check('T4 ' + coh + ' F7e', f7e, r['F7e_p'], 1.5e-3))

# ============================================== 汇总
log('')
log('=' * 90)
n_ok = sum(1 for r in R if r['ok'])
n_bad = len(R) - n_ok
log('TOTAL ' + str(len(R)) + ' checks | consistent ' + str(n_ok) + ' | INCONSISTENT ' + str(n_bad))
if n_bad:
    log('')
    log('BAD LIST:')
    for r in R:
        if not r['ok']:
            log('   [' + r['tag'] + '] ms=' + r['manuscript'] + '  json=' + r['json'] + '  ' + r['note'])

with open(os.path.join(BASE, 'S32_audit_report.json'), 'w', encoding='utf-8') as f:
    json.dump(dict(results=R, n_total=len(R), n_ok=n_ok, n_bad=n_bad), f,
              ensure_ascii=False, indent=2)
with open(os.path.join(BASE, 'S32_audit_stdout.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(OUT))
log('')
log('wrote S32_audit_report.json / S32_audit_stdout.txt')
