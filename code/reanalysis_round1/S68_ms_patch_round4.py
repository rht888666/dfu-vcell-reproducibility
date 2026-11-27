# -*- coding: utf-8 -*-
import sys, io, shutil, copy, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
from docx import Document
from docx.oxml.ns import qn

P = r'E:\虚拟细胞论著投稿\revision_work\round2\biomedicines-4701005_revised_clean_round2.docx'
d = Document(P)
log = []

def sub_text(p, pairs):
    t = p.text
    for a, b in pairs:
        if a in t:
            t = t.replace(a, b)
    if t != p.text:
        for r in list(p.runs):
            r._element.getparent().remove(r._element)
        p.add_run(t)
        return True
    return False

def sub_everywhere(doc, pairs, tag):
    n = 0
    for p in doc.paragraphs:
        if sub_text(p, pairs):
            n += 1
    for tb in doc.tables:
        for row in tb.rows:
            for c in row.cells:
                for p in c.paragraphs:
                    if sub_text(p, pairs):
                        n += 1
    log.append(f'{tag}: {n}')

# ---- global text pairs ----
pairs = [
    # duplicated parenthetical in 2.1
    ('(14 ulcer wound-bed, 19 non-ulcerated skin) (14 ulcer wound-bed, 19 non-ulcerated skin)',
     '(14 ulcer wound-bed, 19 non-ulcerated skin)'),
    # methods 2.3: drop normal-approximation clause
    (', tail probabilities are reported as normal approximations where empirical resolution is insufficient,', ''),
    # methods 2.3: precise 27-gene definition
    ('the 27-gene candidate panel (the union of the program panel genes; listed in the reproducibility package)',
     'the 27-gene candidate panel (the union of programs F6, F7, F7e, F1, and F5, excluding F6i; the pool contains the six F7e genes themselves)'),
    # heading six-scheme -> multi-scheme
    ('Six-Scheme Sequencing-Depth Audit', 'Multi-Scheme Sequencing-Depth Audit'),
    # methods 2.4: specify rank residualization
    ('Rank residualization was performed on the pooled fibroblast population of each dataset because within-library ranking is mathematically degenerate (it zeroes all library means); this choice is a limitation disclosed in Section 4.',
     'Rank residualization was performed on the pooled fibroblast population of each dataset because within-library ranking is mathematically degenerate (it zeroes all library means); this choice is a limitation disclosed in Section 4. Specifically, cell scores were converted to midranks within the pooled population, the ranks were residualized against log(UMI) by ordinary least squares fitted across all gate-selected cells and libraries, residuals were averaged per library, and the library difference was rescaled by the number of gated cells to yield scaled-rank units.'),
    # methods 2.4: dedupe progress definition - remove first vague mention
    (' A composite progress score projecting each library onto the acute-wound reference polyline is reported as a heuristic, cross-dataset descriptive measure only.', ''),
    ("the projection of each library's mean (F6, matrix) coordinate",
     "the projection of each library's mean (F6, F7) coordinate"),
    # curated object attribution
    ('in the curated object', 'in the deposited curated object [3]'),
    # table5 title garble
    ('Table 5. Multi-scheme sequencing-depth audit across 33 foot-tissue libraries of which 14 ulcer libraries in the legacy GSE165816 wound-bed analysis (primary healing contrast: 9 healing versus 5 non-healing units), collagen-free gating sensitivity, and cross-cohort triangulation.',
     'Table 5. Multi-scheme sequencing-depth audit in GSE165816 (primary contrast: 9 healing versus 5 non-healing ulcer wound-bed libraries), with collagen-free gating sensitivity and cross-cohort context.'),
    # gate-size explanation: correct mechanism
    ('Gate sizes differ (36,119 vs 34,891 cells) because the upper-quartile threshold is recomputed within each gate\'s own z-scored marker distribution and the collagen-free gate uses a different marker set (VIM replacing COL1A1/COL1A2).',
     'Both gates retain exactly the upper quartile of the same full cell object, so the gate sizes are identical before tissue selection; the reported counts (36,119 vs 34,891 cells) are the intersections with wound-bed libraries and differ because the two marker sets (the collagen-free gate replaces COL1A1/COL1A2 with VIM) rank different wound-bed cells into the top quartile.'),
    # table4 note: abolishes -> sign reversal
    ('removing ACTA2/TAGLN from F7 abolishes the mural signal in 4/4 donors',
     'removing ACTA2/TAGLN from F7 reverses the sign of the mural signal in 4/4 donors'),
    # discussion abolished -> reversed sign
    ('removing ACTA2/TAGLN from the broad panel abolished the smooth-muscle signal in all four spatial donors',
     'removing ACTA2/TAGLN from the broad panel reversed the sign of the smooth-muscle signal in all four spatial donors'),
    # conclusions
    ('implicated ACTA2/TAGLN as major carriers of mural-cell admixture and excess donor dispersion',
     'implicated ACTA2/TAGLN as carriers of mural-cell admixture and contributors to donor dispersion'),
    ('as the leading candidate regulator of the reparative matrix program',
     'as the leading candidate regulator of this matrix program'),
    # results 3.4 reparative softening + fibrotic caveat
    ('we asked which transcriptional regulators most plausibly drive the reparative matrix program',
     'we asked which transcriptional regulators most plausibly drive this matrix program'),
    ('These results prioritize EGR3 as the leading testable candidate while delimiting exactly what the in silico evidence does and does not support.',
     'These results prioritize EGR3 as the leading testable candidate while delimiting exactly what the in silico evidence does and does not support. Because the same matrix signature has been linked to disease-associated fibroblast states, this prioritization does not by itself imply a reparative effect.'),
    # Southard sentence strengthening
    ('Southard et al. associated collagen-high programs with disease-associated fibroblast states',
     'Southard et al. themselves listed EGR3 among the strong inducers of this collagen program, noted its reported role in fibrogenic responses, and identified EGR3 perturbation as a hit for a disease-associated SPARC+COL3A1+ inflammatory fibroblast state [16]'),
    # limitations: concrete truncation example + fibrotic already there; add acute-chronic extrapolation caveat at end
    ('Seventh, all analyses are computational reanalyses of public data; no new experimental validation was performed.',
     'Seventh, the acute-wound reference datasets cannot establish how fibroblast programs behave in chronic diabetic ulcers, and clinical covariates such as glycaemic control, perfusion, and neuropathy are unavailable in the public metadata; cross-sectional data cannot establish causality. Eighth, all analyses are computational reanalyses of public data; no new experimental validation was performed.'),
    # demonstrably merges -> concrete mechanism
    ('a sample-label truncation in the curated object demonstrably merges distinct subjects',
     'a sample-label truncation in the deposited object merges distinct subject identifiers that share a truncated prefix'),
    ('a sample-label truncation in the deposited curated object [3] merges distinct subjects',
     'a sample-label truncation in the deposited object merges distinct subject identifiers that share a truncated prefix'),
]
sub_everywhere(d, pairs, 'global-pairs')

# normalize truncation-example phrasing after the curated-object replacement
pat = re.compile(r'a sample-label truncation in the deposited curated object \[3\] (demonstrably |incorrectly )?merges distinct subjects')
for p in d.paragraphs:
    if pat.search(p.text):
        t = pat.sub('a sample-label truncation in the deposited object merges distinct subject identifiers that share a truncated prefix', p.text)
        for r in list(p.runs):
            r._element.getparent().remove(r._element)
        p.add_run(t)
        log.append('truncation phrasing normalized')

# ---- Table 2: insert F7e5 / F7noCT / F7epCT after F7e row (row index 4) ----
t2 = d.tables[2]
if not any('F7e5' in r.cells[0].text for r in t2.rows):
    new_rows = [
        ('F7e5 structural matrix vector', 'COL1A1, COL1A2, COL3A1, FN1, MMP2',
         'CRISPRa ranking vector; F7e excluding POSTN, which was not measured in the screen'),
        ('F7-ACTA2/TAGLN ablation panel', 'F7 minus ACTA2 and TAGLN',
         'Contractile-marker ablation; tests mural contribution in longitudinal and spatial analyses'),
        ('F7e+ACTA2/TAGLN addition panel', 'F7e plus ACTA2 and TAGLN',
         'Reverse-direction ablation control; tests robustness of the EGR3 prioritization'),
    ]
    def insert_row_after(table, after_idx, values):
        ref = table.rows[after_idx]._tr
        new_tr = copy.deepcopy(ref)
        for tc in new_tr.findall(qn('w:tc')):
            for p_ in tc.findall(qn('w:p')):
                for r_ in p_.findall(qn('w:r')):
                    p_.remove(r_)
        ref.addnext(new_tr)
        row = table.rows[after_idx + 1]
        for i, v in enumerate(values):
            if i < len(row.cells):
                row.cells[i].paragraphs[0].add_run(v)
    for vals in reversed(new_rows):
        insert_row_after(t2, 4, vals)
    log.append('table2: 3 rows inserted')

# Table 2 note: fix final sentence (P31) - read current text first
for p in d.paragraphs:
    if p.text.startswith('Note: Axis-M7 is a direction-signed'):
        t = p.text
        t = re.sub(r'EGR3 rankings[^.]*refer to Axis-M7\.?',
                   'EGR3 rankings are reported separately for Axis-A, Axis-M7, and F7e5. None directly tests the complete six-gene F7e panel because POSTN was not measured.',
                   t)
        # fallback: if sentence absent, append
        if 'reported separately for Axis-A' not in t:
            t = t.rstrip() + ' EGR3 rankings are reported separately for Axis-A, Axis-M7, and F7e5. None directly tests the complete six-gene F7e panel because POSTN was not measured.'
        for r in list(p.runs):
            r._element.getparent().remove(r._element)
        p.add_run(t)
        log.append('table2 note fixed')
        break

# ---- Table 4 fixes ----
t4 = d.tables[4]
for row in t4.rows:
    c0 = row.cells[1].text.strip()
    # temporal row: exceedance count phrasing
    if 'F7e (Core ECM' in c0:
        cell = row.cells[4]
        t = cell.text
        t = t.replace('p = 0.002 (empirical, 500 panel-matched null draws; observed',
                      'fewer than 1 in 500 null draws exceeded the observed correlation (empirical resolution floor; 500 panel-matched null draws; observed')
        if t != cell.text:
            for r in list(cell.paragraphs[0].runs):
                r._element.getparent().remove(r._element)
            cell.paragraphs[0].add_run(t)
            log.append('table4 temporal p reworded')
    # double-diff rows: add null mean
    if 'F6 - F5' in c0 or 'F6 - F1' in c0:
        cell = row.cells[3]
        t = cell.text
        if 'null mean' not in t.lower():
            t = t.replace('Pair null SD', 'Pair null mean ~0 (difference of symmetric draws from the same panel null), SD')
            for para in cell.paragraphs:
                for r in list(para.runs):
                    r._element.getparent().remove(r._element)
            cell.paragraphs[0].add_run(t)
            log.append(f'table4 null mean added: {c0[:12]}')
    # soften "specific" in interpretation
    interp = row.cells[5]
    if 'specific relative to' in interp.text:
        t = interp.text.replace('specific relative to', 'beyond the')
        for para in interp.paragraphs:
            for r in list(para.runs):
                r._element.getparent().remove(r._element)
        interp.paragraphs[0].add_run(t)
        log.append('table4 interpretation softened')

# Table 4 title fix (P43 caption paragraph)
for p in d.paragraphs:
    if p.text.startswith('Table 4.') and 'Smooth muscle versus FB1 admixture audit' in p.text:
        t = p.text.replace(', and Smooth muscle versus FB1 admixture audit', '')
        for r in list(p.runs):
            r._element.getparent().remove(r._element)
        p.add_run(t)
        log.append('table4 title fixed')

# ---- Table 5: scaled-rank effect sizes ----
t5 = d.tables[5]
for row in t5.rows:
    c0 = row.cells[0].text.strip()
    eff = row.cells[2]
    if 'Pooled Rank' in c0 and 'F7e p = 0.514' in eff.text and '+0.036' not in eff.text:
        t = eff.text.replace('F7e p = 0.514', 'F7e scaled-rank diff. = +0.036, p = 0.514')
        for para in eff.paragraphs:
            for r in list(para.runs):
                r._element.getparent().remove(r._element)
        eff.paragraphs[0].add_run(t)
        log.append('table5 primary rank eff sizes complete')
    if 'Collagen-Free' in c0 and 'pooled-rank p = 0.075' in eff.text and '+0.107' not in eff.text:
        t = eff.text.replace('pooled-rank p = 0.075', 'pooled-rank +0.107 (p = 0.075)')
        t = t.replace('pooled-rank p = 0.183', 'pooled-rank +0.096 (p = 0.183)')
        for para in eff.paragraphs:
            for r in list(para.runs):
                r._element.getparent().remove(r._element)
        eff.paragraphs[0].add_run(t)
        log.append('table5 collagen-free rank eff sizes complete')

# ---- reference [20] reinsertion in intro ----
for p in d.paragraphs[:15]:
    if ('virtual cell' in p.text.lower() or 'perturbation' in p.text.lower()) and 'deep learning' in p.text.lower() and '[20]' not in p.text:
        # cite after first mention of predictive/perturbation models
        m = re.search(r'(predictive[^.]*?models)', p.text)
        if m:
            t = p.text.replace(m.group(1), m.group(1) + ' [20]', 1)
            for r in list(p.runs):
                r._element.getparent().remove(r._element)
            p.add_run(t)
            log.append('ref20 cited')
            break

d.save(P)
shutil.copy(P, r'E:\虚拟细胞论著投稿\revision_work\round2\biomedicines-4701005_revised_redmarked_round2.docx')
for l in log:
    print(l)
print('DONE')
