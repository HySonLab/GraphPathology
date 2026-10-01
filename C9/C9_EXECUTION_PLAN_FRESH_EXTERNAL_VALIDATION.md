# C9 Execution Plan — Fresh External Validation of Sparse-ST Reconstruction and Cross-Platform Calibration

## Repository

All C9 code, configs, manifests, calibration objects, masks, metrics, predictions, figures, logs, and reports must be created under:

`/home/son/Research/GraphPathology`

C0–C8 are frozen and must not be modified.

---

# Why C9 exists

C7 showed that the original sparse-ST reconstruction pipeline transferred only partially to HER2ST.

C8 then showed:

- `PASS_TARGET_HARMONIZATION`
- `MIXED_FUSION_CALIBRATION`
- overall `MIXED_CROSS_PLATFORM_CALIBRATION`

The strongest C8 result was target harmonization, not image fusion.

Therefore, C9 is the **fresh external-validation checkpoint**.

C9 must use a cohort that was not used anywhere in C0–C8.

HER2ST is no longer eligible as pristine validation because it was used for C8 transportability development.

---

# Publication claims under test

C9 should test the strongest claims that are currently defensible.

## Primary claim

> Sparse measured spatial transcriptomics can reconstruct unmeasured spatial expression in a fresh independent cohort.

## Secondary claim

> Cross-platform target harmonization improves transportability of sparse-ST reconstruction.

## Exploratory/secondary claim

> H&E may provide incremental value beyond sparse-ST interpolation in some settings, but this must be validated rather than assumed.

The paper must **not** depend on the H&E-fusion claim passing.

---

# C9 scientific objective

C9 asks:

1. Does sparse-ST reconstruction generalize to a fresh external cohort?
2. Does the C8 target-harmonization strategy remain useful on a new cohort?
3. Does a simple frozen H&E fusion provide reproducible gain beyond ST interpolation?
4. Do the C5/C7 coverage-response and distance-to-anchor behaviors replicate?
5. Does contiguous-region extrapolation remain systematically harder than interpolation?
6. Is the resulting evidence strong enough to support manuscript submission?

---

# C9 structure

C9 has three stages:

```text
C9A — Fresh external dataset qualification
        ↓ if suitable
C9B — Locked primary external validation
        ↓ optional if another fresh cohort is readily available
C9C — Independent replication cohort
```

C9 is a validation checkpoint, not a method-development checkpoint.

---

# Global rules

1. Keep C0–C8 frozen.
2. Do not use any C9 evaluation labels to redesign the method.
3. Do not add new image encoders.
4. Do not fine-tune ResNet18 or Phikon.
5. Do not add new graph architectures.
6. Do not optimize sampling policies.
7. Do not change gene panels after seeing C9 outcomes.
8. Do not change coverage levels after seeing C9 outcomes.
9. Do not change fusion weights after seeing C9 outcomes.
10. Patient identity is the biological split unit.
11. Hidden-location labels may only be accessed for final evaluation.
12. Any cohort-specific calibration must be fit on calibration patients only and then frozen.
13. C9 evaluation patients must remain untouched during calibration.
14. Report both correlation and error.
15. Do not claim success if correlation improves while MSE/MAE systematically worsen.
16. Do not treat spots as independent patients.
17. Stop after C9 and wait for review.

---

# Frozen method definitions carried into C9

C9 must evaluate a small fixed method set.

## M0 — H&E-only baseline

Use the frozen image-only predictor.

No C9 tuning.

## M1 — Sparse-ST interpolation only

This is the primary molecular baseline.

Use the same simple interpolation family frozen from C5/C7.

No C9 hyperparameter search beyond any strictly predeclared cohort-independent setting already frozen.

## M2 — Target-harmonized sparse-ST interpolation

Use the C8 target-harmonization strategy that passed C8A.

Primary choice:

**T1 — patient-disjoint per-gene z-score harmonization**

Important:

- any destination-cohort scaling statistics must be fit only on C9 calibration patients;
- evaluation patients must not contribute to those statistics;
- the exact calibration/evaluation split must be frozen before outcome evaluation.

## M3 — Conservative H&E + sparse-ST fusion

Use exactly one simple frozen fusion rule.

Recommended frozen rule:

**global image weight = 0.2**

because this was the dominant selected C8 fusion in most folds and is simple/interpretable.

Form:

\[
\hat{Y}
=
0.2\hat{Y}^{image}
+
0.8\hat{Y}^{ST}
\]

Do not tune the 0.2 weight on C9.

If the project repository already contains a more explicitly frozen final C8 fusion rule, use that instead and record it in `C9A_SPEC_LOCK.json`.

Do not choose between multiple C8 fusion rules using C9 outcomes.

---

# C9A — Fresh external dataset qualification

## Goal

Identify a truly independent cohort with:

- paired H&E;
- spatial transcriptomics;
- spatial coordinates;
- patient/sample identifiers;
- sufficient tissue coverage;
- sufficient overlap with the frozen gene panels;
- no sample used in C0–C8.

Prefer HEST-compatible datasets when practical.

---

# Automatic exclusions

Exclude:

- NCBI783
- NCBI785
- TENX95
- TENX99
- HER2ST A1–H1
- any duplicate/reprocessed form of the above
- any sample used in C0–C8 experiments, tuning, calibration, or diagnostics

Verify using metadata, source study, sample IDs, and hashes where available.

---

# Fresh-cohort preference order

## Tier 1 — Independent breast-cancer spatial-transcriptomics cohort

Preferred.

## Tier 2 — Independent breast-cancer cohort on a different platform

Strong transportability test.

## Tier 3 — Another epithelial cancer/tissue

Use only if no adequate breast-cancer cohort exists.

If using Tier 3, explicitly downgrade the interpretation from same-disease replication to broader transportability.

---

# Minimum cohort requirement

Preferred:

`>= 6 independent patients`

Minimum:

`>= 4 independent patients`

Reason:

C9 must separate calibration from evaluation patients.

If fewer than 4 independent patients are available, C9A should return:

`NO_SUITABLE_FRESH_EXTERNAL_COHORT`

unless a second independent cohort can be combined without violating study independence.

---

# Calibration/evaluation split

Before accessing C9 outcomes, freeze a patient-level split.

Preferred:

- 50% calibration patients
- 50% held-out evaluation patients

For 6–8 patients:
- use at least 3 calibration patients;
- use at least 3 evaluation patients.

For exactly 4 patients:
- 2 calibration;
- 2 evaluation.

The held-out evaluation patients form the primary C9 external-validation set.

Do not cross-fit across evaluation patients after outcome access.

---

# Primary gene panel

Primary:

**intersection of the original frozen S2-50 panel with the C9 assay**

Do not replace missing genes.

Recommended qualification threshold:

`>= 80% overlap`

If overlap is lower, C9A may still proceed only if:
- overlap remains scientifically meaningful;
- exact overlap is reported;
- status is marked limited.

Secondary:

- intersection with frozen S2-280 panel.

Do not call this whole-transcriptome validation.

---

# Dataset qualification manifest

Create:

`C9/C9A_EXTERNAL_DATASET_MANIFEST.csv`

Required fields:

```text
dataset_id
study
patient_id
sample_id
disease
platform
num_locations
num_assay_genes
overlap_S2_50
overlap_S2_280
image_resolution
spatial_coordinates_available
duplicate_check_pass
cohort_role
qualification_status
exclusion_reason
```

---

# C9A freeze

Before evaluating hidden labels, write and hash:

`C9/C9A_SPEC_LOCK.json`

Freeze:

- cohort;
- patient IDs;
- calibration patients;
- evaluation patients;
- gene intersections/order;
- M0/M1/M2/M3 definitions;
- image encoder;
- calibration transform;
- sampling strategy;
- coverage levels;
- mask regimes;
- random seeds;
- metric definitions;
- pass criteria.

Only after this lock exists may C9B evaluate held-out patients.

---

# C9A status

Assign one of:

- `PASS_FRESH_EXTERNAL_DATASET`
- `LIMITED_FRESH_EXTERNAL_DATASET`
- `NO_SUITABLE_FRESH_EXTERNAL_COHORT`
- `FAIL_TECHNICAL`

If no suitable cohort exists, STOP and return the reason.

---

# C9B — Locked primary external validation

## Primary question

> Does target-harmonized sparse-ST reconstruction generalize to held-out patients from a fresh external cohort?

---

# Sampling strategy

C9 is not a sampling-optimization checkpoint.

Use:

## S0 — Random/open-tissue sampling

Primary.

At least:

`10 random seeds`

per evaluation patient and coverage where feasible.

## S1 — Simple spatially regular sampling

Secondary.

Use a frozen label-free spatial coverage policy such as farthest-point/k-center sampling.

Do not use morphology-guided or adaptive sampling as the primary result because C6 did not pass sampling optimization.

---

# Coverage levels

Use:

```text
10%
20%
40%
```

Optional:

`5%`

only if the cohort has enough locations.

Do not add new coverage values after seeing outcomes.

---

# Mask regimes

## V1 — Open-tissue interpolation

Primary.

Measured anchors distributed across tissue.

## V2 — Contiguous hidden-region extrapolation

Mandatory.

A compact region is completely unmeasured.

## V3 — Local-to-global extrapolation

Optional if geometry permits.

Keep V1 and V2 separate in all headline results.

---

# Primary comparisons

At every budget compare:

```text
M0  H&E only
M1  raw/frozen ST interpolation
M2  target-harmonized ST interpolation
M3  target-harmonized H&E + ST fusion
```

The paper does not require M3 to beat M2.

The most important result is whether M2 provides robust sparse-ST reconstruction on fresh evaluation patients.

---

# Primary metrics

Compute only on hidden locations.

Per gene:

- Pearson
- Spearman
- MSE
- MAE

Per patient/panel:

- mean gene-wise Pearson
- median gene-wise Pearson
- mean/median Spearman
- mean/median MSE
- mean/median MAE
- fraction of genes with Pearson > 0
- fraction with Pearson > 0.1
- fraction with Pearson > 0.2
- fraction with Pearson > 0.3

Report each evaluation patient separately.

Do not report spot-level p-values as if spots were independent patients.

---

# Primary C9 effect sizes

## Sparse-ST value over H&E

\[
\Delta^{ST-H\&E}
=
Metric(M2)-Metric(M0)
\]

This is the primary publication effect.

## Harmonization value

\[
\Delta^{harm}
=
Metric(M2)-Metric(M1)
\]

This tests whether C8 target harmonization transports.

## Incremental H&E value

\[
\Delta^{image}
=
Metric(M3)-Metric(M2)
\]

This is secondary.

A negative/zero `Delta_image` does not invalidate the core paper if M2 validates strongly.

---

# Coverage-response replication

For each evaluation patient and method, generate:

\[
	ext{measured coverage}
ightarrow
	ext{reconstruction fidelity}
\]

Report:

- median Pearson vs coverage;
- MSE vs coverage.

Primary expectation:

more measured ST should generally improve reconstruction.

Do not require perfect monotonicity in every patient.

---

# Distance-to-anchor replication

For every hidden evaluation location:

\[
d_i
=
\min_{j\in\mathcal O}
\|p_i-p_j\|
\]

Use physical units if reliable.

Reproduce:

\[
	ext{error}(d)
\]

with frozen/predeclared bins.

Question:

> Does reconstruction degrade as hidden locations become farther from measured anchors?

Do not force a single extrapolation-horizon number unless the curve supports one.

---

# Contiguous-extrapolation validation

This is mandatory because C5–C8 consistently showed it is harder.

For V2 report separately:

- M1;
- M2;
- M3.

Quantify:

\[
	ext{V2 degradation}
=
Metric(V2)-Metric(V1)
\]

at matched budget.

Do not hide weak V2 results.

A strong paper may legitimately conclude:

> interpolation is reliable while contiguous extrapolation remains a limitation.

---

# C9 primary pass criteria

Assign one of:

- `PASS_FRESH_EXTERNAL_VALIDATION`
- `PASS_LIMITED_FRESH_EXTERNAL_VALIDATION`
- `MIXED_FRESH_EXTERNAL_VALIDATION`
- `FAIL_FRESH_EXTERNAL_VALIDATION`
- `FAIL_TECHNICAL`

## PASS_FRESH_EXTERNAL_VALIDATION

Require:

1. dataset qualification/integrity pass;
2. M2 substantially beats H&E-only M0 on held-out evaluation patients;
3. M2 improves over or is more stable than M1 in the fresh cohort;
4. effect is positive on more than one evaluation patient;
5. sparse-ST benefit is present at multiple coverages;
6. distance-to-anchor behavior is coherent;
7. V2 contiguous extrapolation is evaluated honestly;
8. correlation gains are not systematically contradicted by MSE/MAE;
9. no C9 evaluation patient influenced calibration.

H&E fusion M3 does **not** need to pass for overall C9 PASS.

## PASS_LIMITED_FRESH_EXTERNAL_VALIDATION

Use if:
- the main sparse-ST result replicates;
- but patient number, panel overlap, or platform limitations restrict strength.

## MIXED_FRESH_EXTERNAL_VALIDATION

Use if:
- some evaluation patients replicate;
- others do not;
- or M2 helps only at selected coverages/metrics.

Do not retune on C9 to rescue the result.

## FAIL_FRESH_EXTERNAL_VALIDATION

Use if M2 does not meaningfully beat H&E-only and/or fails to reconstruct hidden expression reliably after technical audits pass.

---

# C9C — Optional second fresh cohort

Run only if:
- C9B is complete and frozen;
- a second independent cohort is readily available;
- no method changes are made.

Purpose:

> replicate the frozen C9B conclusions in another study/platform.

Do not use C9C to redesign the method.

If run, report separately.

---

# Negative controls

At minimum on one held-out evaluation patient and representative budget:

## NC1 — Anchor-expression permutation

Permute observed molecular anchors among observed locations.

Expected:
M1/M2/M3 degrade.

## NC2 — Anchor-position permutation

Assign expression vectors to incorrect observed locations.

Expected:
spatial interpolation benefit decreases.

## NC3 — Image-location permutation

Permute H&E/image predictions across locations.

Expected:
any M3 incremental image benefit decreases.

Use multiple seeds where feasible.

---

# Leakage and calibration audit

For every C9 evaluation patient verify:

1. patient is absent from calibration-set scaler fitting;
2. patient is absent from fusion-weight fitting;
3. hidden labels are absent from all inputs;
4. mask selection uses no hidden labels;
5. predictions are invariant when hidden labels are replaced with NaN/random values.

Any failure:

`FAIL_TECHNICAL`

---

# Required figures

At minimum:

1. C9 external H&E + sampled anchors;
2. M0/M1/M2/M3 comparison;
3. coverage vs Pearson;
4. coverage vs MSE;
5. V1 interpolation vs V2 contiguous extrapolation;
6. error vs nearest-anchor distance;
7. per-patient effect sizes;
8. harmonization gain M2 vs M1;
9. incremental image gain M3 vs M2;
10. development/HER2ST/C9 comparison using comparable metrics.

Do not cherry-pick patients or genes.

---

# Development / C7 / C8 / C9 comparison

Create:

`C9/C9_TRANSPORTABILITY_SUMMARY.csv`

Fields:

```text
cohort
patient
role
panel
coverage
mask_regime
method
median_pearson
mse
mae
gain_vs_H&E
gain_vs_raw_interpolation
gain_vs_harmonized_interpolation
```

The purpose is to show the progression transparently:

```text
development
→ HER2ST zero-shot
→ HER2ST calibrated
→ fresh C9 validation
```

Do not hide domain shift.

---

# Required outputs

Create:

```text
C9/
├── C9_REPORT.md
├── C9_STATUS.json
├── C9_IMMEDIATE_SUMMARY.txt
├── C9_TRANSPORTABILITY_SUMMARY.csv
│
├── C9A/
│   ├── C9A_REPORT.md
│   ├── C9A_SPEC_LOCK.json
│   ├── C9A_EXTERNAL_DATASET_MANIFEST.csv
│   ├── C9A_GENE_INTERSECTION.csv
│   └── logs/
│
├── C9B/
│   ├── C9B_REPORT.md
│   ├── C9B_PANEL_METRICS.csv
│   ├── C9B_GENE_METRICS.csv
│   ├── C9B_COVERAGE_CURVES.csv
│   ├── C9B_DISTANCE_TO_ANCHOR.csv
│   ├── C9B_PATIENT_EFFECTS.csv
│   ├── C9B_NEGATIVE_CONTROLS.csv
│   ├── masks/
│   ├── predictions/
│   ├── logs/
│   └── figures/
│
└── C9C/     # optional
    ├── C9C_REPORT.md
    ├── metrics/
    ├── predictions/
    └── figures/
```

Append all runs to:

`results/experiment_registry.csv`

Do not overwrite prior rows.

---

# Reproducibility requirements

For every run save:

- dataset/study ID;
- patient/sample ID;
- cohort role;
- raw asset hashes where practical;
- gene intersection/order;
- calibration transform;
- calibration-patient IDs;
- encoder/cache hash;
- mask ID/checksum;
- observed barcode list;
- hidden barcode list;
- method config;
- fusion weight;
- random seed;
- prediction hash;
- runtime;
- peak RAM.

Replay saved predictions where practical.

Verify C0–C8 hashes remain unchanged.

---

# Mandatory immediate summary file

Codex must create:

`C9/C9_IMMEDIATE_SUMMARY.txt`

Keep it short: about 40–80 lines.

Use exactly:

```text
C9 IMMEDIATE SUMMARY

1. STATUS
C9A: ...
C9B: ...
C9C: ...
Overall: ...

2. FRESH EXTERNAL COHORT
Dataset/study:
Patients:
Calibration patients:
Evaluation patients:
Platform:
Primary gene overlap:
Why it is fresh/independent:

3. PRIMARY SPARSE-ST RESULT
Coverage:
H&E only:
Raw ST interpolation:
Harmonized ST interpolation:
H&E + harmonized ST:
Patient-level support:

4. TARGET HARMONIZATION
Did C8 harmonization replicate? YES / PARTLY / NO
Pearson effect:
MSE effect:

5. H&E INCREMENTAL VALUE
Did H&E improve beyond harmonized ST interpolation? YES / PARTLY / NO
Patient count improved:
Main effect size:

6. CONTIGUOUS EXTRAPOLATION
Result:
V1 vs V2 degradation:
Main limitation:

7. DISTANCE / COVERAGE
Coverage-response replicated? YES / PARTLY / NO
Distance-to-anchor behavior replicated? YES / PARTLY / NO

8. PUBLICATION READINESS
Is the core sparse-ST reconstruction claim externally validated? YES / PARTLY / NO
Is the H&E-fusion claim externally validated? YES / PARTLY / NO

9. RECOMMENDATION
Exactly one of:
PROCEED_TO_MANUSCRIPT
PROCEED_TO_SECOND_REPLICATION_COHORT
NEED_METHOD_REVIEW
STOP_TECHNICAL

10. FULL REPORT
Path:
Git commit:
```

---

# Codex final-response requirement

When C9 finishes, Codex must return the contents of:

`C9/C9_IMMEDIATE_SUMMARY.txt`

directly in its final response.

Do not return only a GitHub link.

The user should be able to see the C9 result immediately.

---

# Final recommendation

Choose exactly one:

- `PROCEED_TO_MANUSCRIPT`
- `PROCEED_TO_SECOND_REPLICATION_COHORT`
- `NEED_METHOD_REVIEW`
- `STOP_TECHNICAL`

## PROCEED_TO_MANUSCRIPT

Use if the fresh held-out C9 evaluation patients support the core sparse-ST reconstruction claim with acceptable cross-platform stability.

H&E fusion may remain secondary/mixed.

## PROCEED_TO_SECOND_REPLICATION_COHORT

Use if C9 is positive but limited by patient count, overlap, or platform.

## NEED_METHOD_REVIEW

Use if the fresh cohort fails to support the main sparse-ST reconstruction claim.

---

# Explicitly forbidden in C9

Do not:

- modify C0–C8;
- reuse HER2ST as the fresh validation cohort;
- tune on C9 held-out evaluation labels;
- change gene panels after outcome access;
- change coverage levels after outcome access;
- change fusion weights after outcome access;
- add new encoders;
- fine-tune encoders;
- revive multiresolution graph modeling;
- optimize sampling policies;
- run GNNExplainer;
- run pathway enrichment before C9 is frozen;
- claim whole-transcriptome reconstruction;
- claim clinical utility;
- claim monetary cost savings;
- proceed automatically to manuscript generation.

---

# End-of-checkpoint instruction

Before C9B outcome evaluation:

1. qualify a fresh cohort;
2. freeze calibration/evaluation patient split;
3. freeze primary gene intersection;
4. freeze M0–M3;
5. freeze fusion weight;
6. freeze coverage levels;
7. freeze mask regimes;
8. freeze random seeds;
9. freeze metric definitions and pass criteria;
10. write/hash `C9A_SPEC_LOCK.json`.

During execution:

11. never use evaluation-patient labels for calibration;
12. save all calibration objects;
13. run corruption-invariance checks;
14. preserve all prior hashes.

After execution:

15. write `C9A_REPORT.md`;
16. write `C9B_REPORT.md`;
17. optionally run/freeze C9C;
18. write `C9_REPORT.md`;
19. write `C9_STATUS.json`;
20. create `C9_IMMEDIATE_SUMMARY.txt`;
21. verify C0–C8 unchanged;
22. commit/push C9;
23. return `C9_IMMEDIATE_SUMMARY.txt` directly to the user;
24. STOP.

Do not proceed automatically to manuscript generation.
