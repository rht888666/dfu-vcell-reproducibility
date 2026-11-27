# -*- coding: utf-8 -*-
import sys, io, shutil, copy
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
from docx import Document
from docx.oxml.ns import qn

P = r'E:\虚拟细胞论著投稿\revision_work\round2\biomedicines-4701005_revised_clean_round2.docx'
shutil.copy(P, P.replace('.docx', '_pre67.docx'))
d = Document(P)

def sub_all(doc, pairs, name=''):
    n = 0
    def do(paragraphs):
        nonlocal n
        for p in paragraphs:
            t = p.text
            for a, b in pairs:
                if a in t:
                    t = t.replace(a, b)
            if t != p.text:
                n += 1
                for r in list(p.runs):
                    r._element.getparent().remove(r._element)
                p.add_run(t)
    do(doc.paragraphs)
    for tb in doc.tables:
        for row in tb.rows:
            for c in row.cells:
                do(c.paragraphs)
    print(f'{name}: {n} paras/cells changed')

pairs = [
    ('The upper-quartile reanalysis used the original frozen z-space as a legacy reference input.', ''),
    ('(14 ulcer wound-bed, 19 non-ulcerated skin), comprising 14 ulcer wound-bed and 19 non-ulcerated skin libraries',
     '(14 ulcer wound-bed, 19 non-ulcerated skin)'),
    ('identifier crosswalk, sample-label audit, and literature overlap check',
     'identifier crosswalk and sample-label audit'),
    ('high-depth admixture sensitivity audit, systematic smooth muscle versus FB1 admixture audit',
     'high-depth admixture sensitivity audit'),
    ('Smooth muscle versus FB1 admixture audit; spatial and temporal confidence structure (4 donors).',
     'Spatial and temporal confidence structure (4 donors).'),
    ('ACTA2/TAGLN removal abolished the pooled smooth-muscle association',
     'ACTA2/TAGLN removal reversed the pooled smooth-muscle association'),
    ('removal abolished the pooled smooth-muscle association',
     'removal reversed the pooled smooth-muscle association'),
    ('removal abolishes the pooled positive association',
     'removal reverses the sign of the pooled positive association'),
    ('This ablation abolishes the smooth-muscle contrast',
     'This ablation reverses the sign of the smooth-muscle contrast'),
    ('implicating ACTA2 and TAGLN as major carriers of excess donor dispersion',
     'implicating ACTA2 and TAGLN as contributors to donor dispersion'),
    (' When fewer than 500 draws are available, tail probabilities are reported as normal approximations using the empirical null mean and SD.', ''),
    ('where fewer than 500 draws are available, tail probabilities are reported as normal approximations using the empirical null mean and SD',
     'where the empirical resolution permits'),
    ('the 27-gene union of all five scored programs',
     'the 27-gene union of programs F6, F7, F7e, F1, and F5 (F6i excluded); the pool contains the six F7e genes themselves'),
    ('Rank residualization is reported separately as a scale-invariance stress test.',
     'Rank residualization is reported separately as a scale-invariance stress test: within each library, cell scores are converted to ranks, the ranks are residualized against log(UMI) by ordinary least squares pooled across libraries, residuals are averaged per library, and the library difference is rescaled by the number of gated cells.'),
    ('A cell-level ordinary least squares regression of program scores on log(UMI counts) was fitted',
     'A cell-level ordinary least squares regression of program scores on log(UMI counts) was fitted pooled across all gated cells and libraries'),
    ('the projection of the (F6, matrix) coordinate',
     'the projection of the (F6, F7) coordinate'),
    ('the library crosswalk in the curated object',
     'the library annotations in the GEO sample records'),
    ('Southard et al. also associated collagen-high programs with disease-associated fibroblast states.',
     'Southard et al. themselves listed EGR3 among the strong inducers of this collagen program, noted its reported role in fibrogenic responses, and identified EGR3 perturbation as a hit for a disease-associated SPARC+COL3A1+ inflammatory fibroblast state [16].'),
    ('IL6, and CCL2 to levels resembling the strongest observed matrix-axis shifts.',
     'IL6, and CCL2 to levels resembling the strongest observed matrix-axis shifts. Because the same matrix signature has been linked to disease-associated fibroblast states, this prioritization does not by itself imply a reparative effect.'),
    ('implicating EGR3 as the leading candidate regulator of the reparative matrix program',
     'implicating EGR3 as the leading candidate regulator of this matrix program'),
    ('in silico perturbation', 'CRISPRa perturbation screen'),
    ('spot-level distance correlations of donor mean programs; correlations are descriptive without inferential claims',
     'spot-level distance correlations of donor mean programs (n = 6,044 HE-Fibro spots in the 500 to 3500 um band); correlations are descriptive without inferential claims'),
]
sub_all(d, pairs, 'text-pairs')

# reference [20] reinsertion in intro
for p in d.paragraphs:
    if 'few of these perturbation-trained models have been evaluated' in p.text and '[20]' not in p.text:
        t = p.text.replace('few of these perturbation-trained models have been evaluated',
                           'few of these perturbation-trained models have been evaluated [20]')
        for r in list(p.runs):
            r._element.getparent().remove(r._element)
        p.add_run(t)
        print('ref20 cited in intro')

# Table 2: insert F7e5 / F7noCT / F7epCT rows after F7e row
t2 = d.tables[2]
print('table2 first cell:', t2.rows[0].cells[0].text.strip()[:40])
new_rows = [
    ('F7e5', 'F7e excluding POSTN', '5',
     'COL1A1, COL1A2, COL3A1, FN1, MMP2; used in the CRISPRa ranking because POSTN was not measured in the screen'),
    ('F7-ACTA2/TAGLN', 'F7 contractile-marker ablation', '5',
     'F7 minus ACTA2 and TAGLN; tests contractile-marker contribution'),
    ('F7e+ACTA2/TAGLN', 'F7e plus contractile markers', '8',
     'F7e plus ACTA2 and TAGLN; tests robustness of the EGR3 prioritization'),
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

f7e_idx = None
for i, row in enumerate(t2.rows):
    if row.cells[0].text.strip() == 'F7e':
        f7e_idx = i
        break
print('F7e row idx:', f7e_idx)
if f7e_idx is not None and not any(r.cells[0].text.strip() == 'F7e5' for r in t2.rows):
    for vals in reversed(new_rows):
        insert_row_after(t2, f7e_idx, vals)
    print('inserted 3 rows into table2')

# Table 2 note final sentence
for p in d.paragraphs:
    if 'EGR3 rankings refer to Axis-M7' in p.text:
        t = p.text.replace(
            'EGR3 rankings refer to Axis-M7',
            'EGR3 rankings are reported separately for Axis-A, Axis-M7, and F7e5. None directly tests the complete six-gene F7e panel because POSTN was not measured. EGR3 rankings')
        # fix resulting doubled phrase
        t = t.replace('EGR3 rankings for Axis-M7', 'those for Axis-M7')
        for r in list(p.runs):
            r._element.getparent().remove(r._element)
        p.add_run(t)
        print('table2 note fixed')

# Table 5: add scaled-rank effect sizes
t5 = d.tables[5]
rank_eff = {
    ('F7', 'Rank'): '+0.140',
    ('F7e', 'Rank'): '+0.036',
    ('F7 (collagen-free gate)', 'Rank'): '+0.107',
    ('F7e (collagen-free gate)', 'Rank'): '+0.096',
}
for row in t5.rows:
    prog = row.cells[0].text.strip()
    scheme = row.cells[1].text.strip()
    key = (prog, scheme)
    if key in rank_eff:
        eff = rank_eff[key]
        cur = row.cells[2].text.strip()
        if eff not in cur:
            row.cells[2].paragraphs[0].add_run(eff)
            print('table5 rank eff added:', prog, eff)

d.save(P)
shutil.copy(P, r'E:\虚拟细胞论著投稿\revision_work\round2\biomedicines-4701005_revised_redmarked_round2.docx')
print('saved')
