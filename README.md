# DFU-VCell reproducibility package (revision draft)

This package contains the analysis scripts, fixed-seed outputs, and processed tables used for the revised manuscript **“Virtual Cell State-Space Projection Prioritizes EGR3 as a Candidate Regulator in Diabetic Foot Ulcers.”**

## Contents

- `code/`: Python scripts S2/S3/S20-S36 and `04_insilico_screen_and_orthogonal_closure.py` used for longitudinal, spatial, clinical-depth, cross-cohort, CRISPRa-prioritization, audit, and figure-related analyses.
- `processed_results/`: machine-generated JSON outputs and audit verdicts used to verify manuscript numbers.
- `processed_tables/`: processed CSV tables, including `TF_CRISPRa_Screen_All.csv`, longitudinal donor summaries, clinical library summaries, spatial scores, edge distances, and transcription-factor rankings.
- `data/perturb_V3.npz`: compressed CRISPRa perturbation matrix used for the Axis-A/Axis-M7 reanalysis.

## Public input datasets

Raw public datasets are not redistributed here. They are available from the original repositories/accessions cited in the manuscript: GSE241132, GSE241124, GSE165816, GSE231643, GSE223964, GSE245703, and the Southard et al. CRISPRa resource.

## Notes for reproduction

1. Several legacy scripts contain the original local root path `C:/Users/Administrator/Desktop/虚拟细胞论著投稿`. Before rerunning locally, set that root to the local package/project path or replace it with the appropriate directory.
2. Processed outputs are provided so that the manuscript numbers can be checked without re-downloading all raw public data.
3. `S36_donor_stratified_spatial_sensitivity.py` was added during revision to calculate donor-stratified Day-7 wound-edge correlations without cross-donor pooling.
4. Random seeds are fixed inside the scripts where permutation procedures are used.

## Revision status

This archive is prepared for public deposit in Zenodo or GitHub at revision. The persistent DOI/URL must be inserted into the manuscript Data Availability Statement after upload.
