# DFU-VCell reproducibility package (revision draft)

This package contains the analysis scripts, fixed-seed outputs, and processed tables used for the revised manuscript **“Virtual Cell State-Space Projection Prioritizes EGR3 as a Candidate Regulator in Diabetic Foot Ulcers.”**

## Contents

- `code/`: Python scripts S2/S3/S20-S38 and `04_insilico_screen_and_orthogonal_closure.py` used for longitudinal, spatial, clinical-depth, cross-cohort, CRISPRa-prioritization, audit, and figure-related analyses.
- `processed_results/`: machine-generated JSON outputs and audit verdicts used to verify manuscript numbers.
- `processed_tables/`: processed CSV tables, including `TF_CRISPRa_Screen_All.csv`, longitudinal donor summaries, clinical library summaries, spatial scores, edge distances, and transcription-factor rankings.
- `data/perturb_V3.npz`: compressed CRISPRa perturbation matrix used for the Axis-A/Axis-M7 reanalysis.

## Public input datasets

Raw public datasets are not redistributed here. They are available from the original repositories/accessions cited in the manuscript: GSE241132, GSE241124, GSE165816, GSE231643, GSE223964, GSE245703, and the Southard et al. CRISPRa resource.

## Notes for reproduction

1. Several legacy scripts contain the original local root path `C:/Users/Administrator/Desktop/虚拟细胞论著投稿`. Before rerunning locally, set that root to the local package/project path or replace it with the appropriate directory.
2. Processed outputs are provided so that the manuscript numbers can be checked without re-downloading all raw public data.
3. `S36_donor_stratified_spatial_sensitivity.py` was added during revision to calculate donor-stratified Day-7 wound-edge correlations without cross-donor pooling.
4. `S37_gating_sensitivity_GSE165816.py` was added during revision to recompute the primary and collagen-free fibroblast-gating sensitivity analyses under the rank-residualization framework. Its machine-generated output is `processed_results/S37_gating_sensitivity_GSE165816.json`.
5. `S38_figure2_spatial_revision.py` regenerates revised manuscript Figure 2 from the corrected spatial provenance outputs. It uses all 22,915 annotated Visium spots for the temporal/compartment analyses and the auditable 6,044 Day-7 wound spots with edge-distance measurements for the Day-7 gradient panel. The generated PNG/PDF are stored in `figures/`, and key reproduction checks are summarized in `processed_results/S38_figure2_reproduction_summary.json`.
6. In GSE165816, “analytical units” refer to the legacy donor-hint aggregated units described in the manuscript; subject_hint is only a weak donor hint rather than verified patient identity.
7. Random seeds are fixed inside the scripts where permutation procedures are used.

## Revision status

This archive is publicly available at https://github.com/rht888666/dfu-vcell-reproducibility. The current revision release is https://github.com/rht888666/dfu-vcell-reproducibility/releases/tag/v1.1.1-revision. It supersedes the initial v1.1-revision archive by adding the Figure 2 regeneration script and outputs; the underlying analysis results are unchanged.
