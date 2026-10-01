# C10 Overnight Execution Plan — Second Replication Cohort + Manuscript Figure/Table Package

## Repository
All C10 work must be created under:

`/home/son/Research/GraphPathology`

C0–C9 are frozen and must not be modified.

This plan is intended to run unattended overnight.

Priority:
1. second fresh replication cohort;
2. publication-quality figures and tables;
3. exact numerical provenance;
4. rough manuscript draft;
5. concise morning summary.

The figures, tables, and traceability are more important than polished prose.

---

# Scientific story entering C10

The strongest currently supported publication story is:

> Sparse spatial transcriptomics measurements can reconstruct unmeasured spatial expression, with reconstruction quality governed by measurement coverage, spatial geometry, and distance to observed molecular anchors.

Do not make the paper primarily about:
- graph neural networks;
- H&E fusion;
- affine target harmonization;
- optimized/adaptive sampling.

Those directions were tested and were mixed or negative.

Supported/promising:
1. Sparse molecular anchors outperform H&E-only prediction.
2. More measured ST generally improves reconstruction.
3. Reconstruction degrades with distance from measured anchors.
4. Contiguous missing-region extrapolation is harder than open-tissue interpolation.
5. These findings have fresh external support, but the external evidence remains limited.

Not supported as headline claims:
1. H&E consistently improves over sparse-ST interpolation.
2. Exact graph topology is essential.
3. Sampling optimization reliably solves contiguous extrapolation.
4. Affine target harmonization produces a distinct interpolation gain in common units.

---

# C10 structure

```text
C10A — Second fresh cohort qualification
        ↓
C10B — Locked second-cohort replication
        ↓
C10C — Cross-cohort synthesis
        ↓
C10D — Figures and tables
        ↓
C10E — Rough manuscript draft
        ↓
C10F — Morning summary
```

If no suitable second cohort can be qualified, skip C10B and continue C10C–C10F using C0–C9. Clearly mark the missing replication as a limitation.

---

# C10A — Second fresh cohort qualification

Find a second independent cohort unused in C0–C9.

Automatically exclude:
- NCBI783
- NCBI785
- TENX95
- TENX99
- HER2ST A1–H1
- Wu cohort patients used in C9
- duplicates/reprocessed copies of prior samples

Preference order:
1. independent breast-cancer ST cohort with paired H&E;
2. independent breast-cancer cohort on a different platform;
3. another epithelial cancer with paired H&E + ST only if no adequate breast cohort exists.

Prefer:
- >=4 independent patients, ideally >=6;
- paired H&E;
- spatial coordinates;
- adequate tissue area;
- >=80% overlap with the frozen S2-50 panel;
- better H&E resolution than the C9 Wu cohort when possible;
- HEST-compatible data when practical.

Create:
- `C10/C10A_DATASET_QUALIFICATION.md`
- `C10/C10A_DATASET_MANIFEST.csv`

Manifest fields:
```text
dataset
study
patient_id
sample_id
disease
platform
num_locations
gene_overlap_S2_50
gene_overlap_S2_280
H&E resolution
coordinate availability
prior-use check
qualification status
reason
```

C10A status:
- `PASS_SECOND_REPLICATION_DATASET`
- `LIMITED_SECOND_REPLICATION_DATASET`
- `NO_SUITABLE_SECOND_REPLICATION_DATASET`
- `FAIL_TECHNICAL`

If no suitable cohort exists after a reasonable search, do not keep searching indefinitely. Record why and proceed with manuscript synthesis from C0–C9.

---

# C10B — Locked second-cohort replication

Run only if C10A qualifies a dataset.

This is validation, not method development.

## Primary method choice

The primary method is **simple sparse-ST interpolation**.

Do not center C10 on C8/C9 affine harmonization because C9 showed affine per-gene scaling commutes with weighted interpolation in common units.

Do not center C10 on H&E fusion because C7–C9 did not validate reproducible H&E gain.

### M0 — H&E-only baseline
Use the existing frozen H&E predictor.

### M1 — Sparse-ST interpolation
Use the frozen simple spatial interpolation method from C9.
This is the **primary method**.

### M2 — Optional affine-harmonized interpolation
Run only as an equivalence/invariance check if useful.
Do not present it as distinct if mathematically identical in common units.

### M3 — Frozen H&E + sparse-ST fusion
Secondary only.
Use the same frozen fusion rule from C9.
Do not tune its weight.
The paper does not require M3 to beat M1.

## Calibration/evaluation split

If cohort-specific normalization requires calibration:
- split patients before outcome evaluation;
- fit calibration/scaling on calibration patients only;
- evaluate on held-out patients.

Prefer approximately half calibration / half evaluation, with at least 2 evaluation patients.

If M1 requires no learned cross-patient calibration, evaluate all eligible patients under the frozen method and report patient-level results separately.

Do not invent unnecessary calibration.

## Gene panel
Primary: frozen S2-50 intersection.
Do not replace missing genes.
Prefer >=80% overlap.

Secondary: frozen S2-280 intersection if useful.
Do not call it whole transcriptome.

## Coverage
Use exactly:
```text
10%
20%
40%
```
Optional 5% only if trivial to execute.

## Regimes

### V1 — Open-tissue interpolation
Primary.
Use random distributed anchors.
At least 10 seeds per patient/coverage where feasible.

### V2 — Contiguous hidden-region extrapolation
Mandatory.
Use multiple predefined missing-region geometries where feasible.
Do not cherry-pick a favorable region.

### V3 — Simple spatially regular sampling
Secondary only.
Use one frozen coverage-maximizing policy.
Do not restart sampling optimization.

## Metrics
Compute only on hidden locations.

Per gene:
- Pearson
- Spearman
- MSE
- MAE

Per patient/panel:
- median gene-wise Pearson
- mean gene-wise Pearson
- median Spearman
- mean MSE
- mean MAE
- fraction Pearson >0, >0.1, >0.2, >0.3

Do not use spot-level p-values as if spots were independent patients.

## C10B primary questions
1. Does M1 beat H&E-only M0?
2. Does performance improve from 10% → 20% → 40%?
3. Does error increase with distance to nearest measured anchor?
4. Is V2 worse than V1?
5. Does M3 provide any reproducible incremental value over M1?

Question 5 is secondary and may remain negative.

## Distance-to-anchor
For each hidden location:
`d_i = distance to nearest measured ST anchor`

Use physical units when reliable.

Report:
- MSE vs distance;
- MAE vs distance;
- Pearson by distance bin where meaningful;
- counts per bin.

Do not force a single extrapolation horizon if unsupported.

## Negative controls
On one representative held-out patient / 20% coverage:

- NC1: permute anchor-expression assignments.
- NC2: permute anchor positions.
- NC3: permute H&E features/predictions for M3.

Use multiple seeds if practical.

## Leakage/integrity
For every run:
1. hidden expression never enters input construction;
2. hidden labels only used for final scoring;
3. save observed/hidden barcode lists;
4. save mask checksums;
5. run NaN/random hidden-label corruption invariance;
6. predictions must remain unchanged;
7. verify C0–C9 hashes unchanged.

Any violation = `FAIL_TECHNICAL`.

## C10B status
Use exactly one:
- `PASS_SECOND_REPLICATION`
- `PASS_LIMITED_SECOND_REPLICATION`
- `MIXED_SECOND_REPLICATION`
- `FAIL_SECOND_REPLICATION`
- `FAIL_TECHNICAL`

PASS requires:
1. M1 clearly beats H&E-only on multiple held-out patients;
2. benefit appears at multiple coverages;
3. coverage-response is coherent;
4. distance-to-anchor degradation is reproduced or clearly characterized;
5. V2 is evaluated honestly;
6. errors do not systematically contradict correlation;
7. technical audits pass.

H&E fusion does not need to pass.

---

# C10C — Cross-cohort synthesis

Create:
`C10/CROSS_COHORT_SUMMARY.csv`

Cover:
- development cohort;
- HER2ST zero-shot;
- HER2ST C8 calibrated/development result;
- Wu fresh C9 cohort;
- C10 second replication cohort if available.

Columns:
```text
cohort
study
patient
role
platform
primary_gene_count
coverage
regime
method
median_pearson
mean_mse
mean_mae
distance_metric_summary
notes
```

Do not compare medians from different gene sets as if directly equivalent.
Clearly mark 40-gene vs 50-gene analyses.
Distinguish development, transportability-development, and fresh external validation.

---

# Claim-evidence matrix

Create:
`C10/CLAIM_EVIDENCE_MATRIX.csv`

Columns:
```text
claim_id
claim_text
status
supporting_checkpoints
supporting_files
key_numbers
limitations
recommended_manuscript_wording
```

Include at least:

1. Sparse ST reconstructs hidden expression better than H&E alone.
2. Increasing measured coverage improves reconstruction.
3. Greater anchor distance generally reduces fidelity.
4. Contiguous extrapolation is harder than open-tissue interpolation.
5. H&E provides reproducible incremental value beyond sparse-ST interpolation.
6. Sampling optimization reliably improves reconstruction.
7. Graph topology is essential.

Expected statuses for 5–7 may be mixed/unsupported.

This file is critical for later human writing.

---

# C10D — Publication figure package

Create:
`C10/MANUSCRIPT/figures/`

For every main figure:
- PNG at minimum 300 dpi;
- PDF or SVG where practical;
- plotting script;
- exact source CSV/data used.

Do not use notebook screenshots.
Use consistent typography, axis labels, panel lettering, and metric definitions.
Do not manually alter values in plotting code.

## Figure 1 — Study overview
Simple technical schematic/sketch only; not the final artistic conceptual figure.

Show:
```text
partial measured ST + tissue coordinates (+ optional H&E baseline)
        ↓
sparse-ST reconstruction
        ↓
hidden expression
        ↓
coverage / distance / geometry evaluation
```

Show V1 distributed anchors vs V2 contiguous missing region.

Filename:
`Fig1_study_overview.*`

## Figure 2 — Development reconstruction
Use frozen C5 results.

Include:
- accuracy vs coverage;
- MSE vs coverage;
- H&E-only vs sparse-ST interpolation vs fusion where relevant;
- V1 vs V2.

Filename:
`Fig2_development_reconstruction.*`

## Figure 3 — Geometry and distance govern reliability
Use C5/C6 evidence.

Include:
- error vs nearest-anchor distance;
- V1 vs V2;
- representative spatial error maps;
- optional simple spatial-sampling result as secondary.

Headline is geometry/reliability, not sampling optimization.

Filename:
`Fig3_geometry_distance.*`

## Figure 4 — HER2ST transportability
Use C7/C8.

Include:
- H&E-only vs ST interpolation vs fusion;
- C7 zero-shot;
- calibration/scale issue;
- V1 vs V2;
- clearly show H&E fusion is mixed.

Do not present affine harmonization as a unique interpolation advantage.

Filename:
`Fig4_HER2ST_transportability.*`

## Figure 5 — Fresh Wu external validation
Use C9.

Include:
- M0 vs M1/M2/M3 at 20%;
- 10/20/40% coverage-response;
- patient-level effects;
- V1 vs V2;
- distance-to-anchor.

Emphasize:
- sparse-ST strongly beats H&E-only;
- H&E fusion does not add reproducible value;
- contiguous extrapolation is harder.

Filename:
`Fig5_fresh_external_Wu.*`

## Figure 6 — Second replication
If C10B runs successfully, create a main figure with:
- sparse ST vs H&E;
- coverage curve;
- V1 vs V2;
- distance curve;
- patient effects.

Filename:
`Fig6_second_replication.*`

If no cohort qualifies, do not fabricate this figure.

## Figure 7 — Cross-cohort synthesis
Create only if defensible.

Possible panels:
- patient-level sparse-ST effects across cohorts;
- coverage-response side by side;
- V1/V2 degradation across cohorts;
- dataset/platform summary.

Do not pool incomparable gene panels into a misleading single scalar.

Filename:
`Fig7_cross_cohort_synthesis.*`

---

# Supplementary figures

Create:
`C10/MANUSCRIPT/supp_figures/`

Include where data exist:
1. C1 morphology-only baselines and controls.
2. C2 graph-context comparison.
3. C3 topology-randomization controls.
4. C4 encoder comparison.
5. C5 leakage/corruption controls.
6. C6 sampling strategy comparisons.
7. C7 external patient-level distributions.
8. C8 calibration/fusion analyses.
9. C9 negative controls.
10. C10 QC/negative controls.

Negative/mixed results are part of the rigor and should not be hidden.

---

# C10D — Publication table package

Create:
`C10/MANUSCRIPT/tables/`

Save every table as:
- CSV;
- Markdown;
- optionally LaTeX `.tex`.

## Table 1 — Cohort summary
Columns:
```text
cohort
study
disease
platform
patient_count
evaluation_patient_count
location_count
H&E resolution
primary gene overlap
role in project
```

## Table 2 — Primary reconstruction results
For each cohort/coverage:
```text
cohort
coverage
regime
method
median Pearson
MSE
MAE
patient support count
```

## Table 3 — Coverage and geometry
Include:
```text
cohort
10% Pearson/MSE
20% Pearson/MSE
40% Pearson/MSE
V1 Pearson/MSE
V2 Pearson/MSE
V2-V1 degradation
```

## Table 4 — External replication summary
Fresh validation cohorts only:
```text
cohort
held-out patients
primary gene count
sparse-ST beats H&E?
coverage-response?
distance effect?
contiguous extrapolation?
H&E incremental gain?
validation status
```

Supplementary tables:
1. exact gene panels/intersections;
2. per-patient primary metrics;
3. per-gene metrics;
4. negative controls;
5. integrity/replay counts;
6. mask geometry summary;
7. distance-bin metrics;
8. full C0–C10 checkpoint timeline.

---

# Statistical presentation rules

1. Patient is the biological unit for patient-level inference.
2. Do not treat thousands of spots as independent n.
3. Distinguish patient variability, mask/random-seed variability, and model-seed variability.
4. Do not invent p-values.
5. If a valid frozen patient-level test already exists, report it exactly.
6. Otherwise use descriptive effects and patient support counts.
7. Do not add new outcome-driven statistical tests unless explicitly labeled exploratory.

---

# C10E — Rough manuscript draft

Create:
`C10/MANUSCRIPT/MANUSCRIPT_DRAFT.md`

This is a scientific sketch, not final prose.

Priority:
1. correct numbers;
2. correct claims;
3. correct figure/table references;
4. clear limitations;
5. source traceability;
6. prose polish last.

## Working title
Use a conservative title, e.g.:

> Reconstructing Spatial Transcriptomics from Sparse Molecular Measurements Across Independent Cohorts

Alternative:

> Sparse Spatial Transcriptomics Reconstruction Reveals Coverage- and Geometry-Dependent Limits of Molecular Interpolation

Do not put H&E, GNN, or harmonization in the title unless C10 unexpectedly supplies strong new support.

## Abstract
Draft only as placeholder.
Maximum ~250 words.
Mark:
`[DRAFT — REVISE AFTER HUMAN WRITING]`

## Introduction
4–6 paragraph sketch.
Do not invent citations.
Use `[REF]` placeholders.

Core motivation:
- spatial transcriptomics provides spatial molecular information but measurement density/coverage is limited;
- unmeasured locations motivate reconstruction;
- image-only prediction may not transport robustly;
- the practical question is how sparse molecular anchors, coverage, distance, and geometry govern reliability;
- this study systematically evaluates that question across cohorts.

## Results
Organize by scientific findings, not checkpoint numbers.

Suggested subsections:

1. Sparse molecular measurements reconstruct hidden spatial expression.
2. Reconstruction improves with measured coverage.
3. Distance and missing-region geometry govern reliability.
4. H&E-only prediction and image fusion show limited transportability.
5. Findings replicate across independent cohorts.

If C10 replication is mixed/negative, write that accurately.

## Methods
Draft more fully than the Introduction.

Include exact:
- datasets;
- patient split rules;
- target panels;
- preprocessing;
- H&E embeddings;
- masks;
- interpolation algorithms;
- fusion baseline;
- coverage regimes;
- distance analysis;
- contiguous extrapolation;
- negative controls;
- leakage checks;
- reproducibility/freeze protocol;
- hardware.

Use exact frozen config values.
If unavailable write `[VERIFY FROM CONFIG]`; do not guess.

## Discussion
Draft structured notes/prose covering:
- main sparse-ST result;
- external replication;
- coverage/distance dependence;
- contiguous extrapolation limitation;
- limited incremental H&E value;
- why simple interpolation may be strong;
- cross-platform differences;
- sample-size limitations;
- future work.

No unsupported mechanism claims.

---

# Do not fabricate references

Create:
`C10/MANUSCRIPT/CITATION_NEEDS.md`

List places where literature citations are required.

Examples:
- ST cost/coverage limitation;
- representative ST platforms;
- histology-to-expression work;
- CMRCNet;
- GHIST;
- MERGE / graph methods;
- spatial interpolation/statistics;
- cross-platform domain shift.

ChatGPT + user will later perform final literature search and citations.

---

# Figure legends

Create:
`C10/MANUSCRIPT/FIGURE_LEGENDS.md`

Every legend must state:
- cohort;
- patient/sample count;
- coverage;
- method;
- metric;
- whether values are patient medians, means, seed distributions, etc.;
- meaning of error bars;
- development vs external-validation role.

---

# Exact-results ledger

Create:
`C10/MANUSCRIPT/RESULTS_LEDGER.md`

For every important manuscript number record:

```text
Manuscript claim:
Exact number:
Metric definition:
Cohort:
Patients:
Gene panel:
Coverage:
Regime:
Method:
Source checkpoint:
Source file:
Source row/field:
```

This is one of the most important overnight outputs.

---

# Manuscript claim restrictions

Do not write:
- "whole transcriptome" for 50/280-gene panels;
- "clinical utility";
- "cost savings";
- "H&E improves reconstruction" as a general claim;
- "graph topology is essential";
- "optimal sampling";
- "causal";
- "single-cell" unless truly single-cell;
- "external validation" for development/transportability cohorts.

Preferred terminology:
- sparse spatial-transcriptomics reconstruction;
- hidden/unmeasured locations;
- measured-anchor coverage;
- spatial interpolation;
- contiguous-region extrapolation;
- independent cohort replication;
- patient-level held-out evaluation.

---

# C10F — Morning summary

Create:
`C10/C10_OVERNIGHT_SUMMARY.txt`

Target ~60–120 lines.

Use exactly:

```text
C10 OVERNIGHT SUMMARY

1. SECOND REPLICATION
Dataset:
Status:
Patients:
Evaluation patients:
Platform:
Primary gene overlap:

2. PRIMARY RESULT
At 20% V1:
H&E only:
Sparse-ST interpolation:
H&E + ST:
Patient support:

3. COVERAGE
10%:
20%:
40%:
Did coverage-response replicate? YES / PARTLY / NO

4. CONTIGUOUS EXTRAPOLATION
V1:
V2:
Did V2 remain harder? YES / PARTLY / NO

5. DISTANCE
Did error worsen with anchor distance? YES / PARTLY / NO

6. CROSS-COHORT CONCLUSION
Development:
HER2ST:
Wu:
Second replication:
Core sparse-ST claim: STRONG / MODERATE / MIXED / NOT SUPPORTED
H&E incremental claim: SUPPORTED / MIXED / NOT SUPPORTED

7. MANUSCRIPT PACKAGE
Draft:
Main figures:
Supplementary figures:
Main tables:
Supplementary tables:
Results ledger:
Claim-evidence matrix:
Citation-needs file:

8. PUBLICATION READINESS
Ready to begin human manuscript writing? YES / NO
Remaining major evidence gap:

9. FULL OUTPUT
C10 report:
Manuscript directory:
Git commit:
```

---

# Final Codex response requirement

When the overnight run finishes, Codex must:

1. return the full contents of `C10_OVERNIGHT_SUMMARY.txt`;
2. list paths to:
   - `MANUSCRIPT_DRAFT.md`
   - `FIGURE_LEGENDS.md`
   - `RESULTS_LEDGER.md`
   - `CLAIM_EVIDENCE_MATRIX.csv`
   - main figures directory
   - main tables directory
3. report final Git commit;
4. STOP.

Do not return only a GitHub link.

---

# Required directory structure

```text
C10/
├── C10_REPORT.md
├── C10_STATUS.json
├── C10_OVERNIGHT_SUMMARY.txt
├── C10A_DATASET_QUALIFICATION.md
├── C10A_DATASET_MANIFEST.csv
├── C10B/
│   ├── C10B_REPORT.md
│   ├── C10B_PANEL_METRICS.csv
│   ├── C10B_GENE_METRICS.csv
│   ├── C10B_COVERAGE_CURVES.csv
│   ├── C10B_DISTANCE_TO_ANCHOR.csv
│   ├── C10B_PATIENT_EFFECTS.csv
│   ├── C10B_NEGATIVE_CONTROLS.csv
│   ├── masks/
│   ├── predictions/
│   ├── logs/
│   └── figures/
├── CROSS_COHORT_SUMMARY.csv
├── CLAIM_EVIDENCE_MATRIX.csv
└── MANUSCRIPT/
    ├── MANUSCRIPT_DRAFT.md
    ├── FIGURE_LEGENDS.md
    ├── RESULTS_LEDGER.md
    ├── CITATION_NEEDS.md
    ├── figures/
    ├── supp_figures/
    ├── tables/
    └── supp_tables/
```

Append C10 events to:
`results/experiment_registry.csv`

Do not overwrite earlier rows.

---

# Publication-readiness recommendation

At end choose exactly one:

- `READY_FOR_HUMAN_MANUSCRIPT_WRITING`
- `READY_WITH_LIMITED_SECOND_REPLICATION`
- `NEED_ONE_MORE_REPLICATION_DATASET`
- `NEED_METHOD_REVIEW`
- `STOP_TECHNICAL`

## READY_FOR_HUMAN_MANUSCRIPT_WRITING
Use if:
- C9 fresh validation remains intact;
- C10 second cohort provides supportive replication;
- core sparse-ST claim is consistent;
- figures/tables are complete;
- no technical integrity issue remains.

H&E fusion may remain negative.

## READY_WITH_LIMITED_SECOND_REPLICATION
Use if C10 replication is supportive but cohort quality/patient count is limited.

## NEED_ONE_MORE_REPLICATION_DATASET
Use only if C10 cannot qualify a usable cohort or evidence remains too thin.

Do not add cohorts indefinitely if two independent validation studies already support the core claim.

---

# Explicit prohibitions

Do not:
- modify C0–C9;
- retune the primary method on C10 evaluation patients;
- redesign method based on C10 outcome;
- add new encoders;
- add new GNN architectures;
- restart sampling optimization;
- hide negative H&E results;
- fabricate references;
- fabricate missing experimental details;
- claim whole-transcriptome reconstruction;
- claim cost savings;
- claim clinical utility;
- change frozen results during manuscript generation;
- continue beyond C10 automatically.

---

# Final priority reminder

Overnight priority:

```text
1. rigorous second replication
2. publication-quality figures
3. complete tables
4. exact numerical provenance
5. manuscript sketch
6. prose polish
```

The user and ChatGPT will later rewrite the manuscript, add literature citations, refine the conceptual figure, and improve the narrative.

Codex should optimize for correct, traceable scientific evidence and publication-ready visual/data assets, not elegant prose.
