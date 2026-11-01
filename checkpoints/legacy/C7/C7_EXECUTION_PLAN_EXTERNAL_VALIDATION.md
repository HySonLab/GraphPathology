# C7 Execution Plan — External Validation of Sparse Spatial Transcriptomics Reconstruction

## Repository
All C7 code, configs, dataset manifests, metrics, predictions, figures, logs, and reports must be created under:

`/home/son/Research/GraphPathology`

C0–C6 are frozen and must not be modified.

## Frozen context
Development patients used in C0–C6:
- `NCBI783`
- `NCBI785`
- `TENX95`
- `TENX99`

These four patients are development data only and must never be counted as external validation.

C5 status: `PASS_SPARSE_ST_RECONSTRUCTION`.

C6 status: `MIXED_SAMPLING_OPTIMIZATION`, frozen commit:
`25964a0ff1b796b73bfe761c838a7d1d6f9095e1`

C6 showed that sampling placement can help open-tissue interpolation, but the benefit did not reliably extend to contiguous missing regions. C6B adaptive sampling was therefore not run.

---

# Scientific objective

C7 asks:

> Do the main sparse-ST reconstruction findings from C5 generalize to an independent cohort that was never used for model development, target selection, hyperparameter selection, encoder selection, mask-policy design, or sampling optimization?

Primary claim under test:

`H&E everywhere + sparse measured ST -> reconstruction of unmeasured spatial expression beyond the development cohort.`

Secondary claim under test:

`Reconstruction reliability depends on measured coverage and distance from observed ST anchors.`

C7 is an external-validation checkpoint, not a new-method checkpoint.

---

# C7 structure

```text
C7A — External dataset qualification and freeze
        ↓ if suitable dataset exists
C7B — Locked external validation
        ↓ optional if data permit
C7C — Secondary transportability stress test
```

Do not modify the frozen C5/C6 method after seeing external outcomes.

---

# Global rules

1. Keep C0–C6 unchanged.
2. Do not use the original four patients as external validation.
3. Prefer an independent study/source cohort, not another slide from the same original study.
4. Prefer the same broad disease domain first (breast cancer / IDC where available).
5. Different disease/tissue may be used only as a secondary transportability test.
6. Do not introduce new image encoders.
7. Do not fine-tune ResNet18 or Phikon on external expression labels.
8. Do not redesign graph architecture.
9. Do not add multiresolution graph models.
10. Do not optimize sampling policy on external labels.
11. Do not change gene panels after inspecting external performance.
12. Do not change coverage levels after inspecting external performance.
13. Do not change model hyperparameters after external-label access.
14. External hidden labels may be accessed only for final evaluation after C7 settings are frozen.
15. Any external-domain adaptation must be a separate secondary analysis and must never replace the primary zero-shot result.
16. Report failures and domain shift explicitly.
17. Stop after C7 for review.

---

# C7A — External dataset qualification and freeze

## Goal
Identify at least one independent cohort with paired:
- H&E histology;
- spatial transcriptomics;
- spot/location coordinates;
- patient/sample identifiers;
- sufficient tissue area;
- sufficient overlap with frozen gene panels.

Prefer HEST-compatible data when practical because the preprocessing pipeline is already audited.

## Automatic exclusions
Exclude:
- NCBI783
- NCBI785
- TENX95
- TENX99
- exact duplicates or reprocessed copies of these samples
- any sample used in C0–C6 development

Check source study, identifiers, metadata, and hashes where available.

## External cohort preference order
1. Independent breast-cancer study/cohort.
2. Independent breast-cancer cohort with different platform/preparation.
3. Different cancer/tissue only as a secondary stress test.

## Dataset manifest
Create:

`C7/C7A_EXTERNAL_DATASET_MANIFEST.csv`

with fields:

```text
dataset_id
study_source
patient_id
sample_id
disease_tissue
platform
num_locations
num_assay_genes
num_overlap_S2_50
num_overlap_S2_250
num_overlap_S2_280
image_resolution
coordinate_available
patient_metadata_available
duplicate_check_pass
qualification_status
exclusion_reason
```

## Gene panels
### S2-50
Mandatory primary panel.

Use the exact frozen S2-50 genes where possible. If genes are missing externally:
- do not replace them;
- use the frozen-panel intersection;
- report exact retained genes and overlap fraction.

Freeze a minimum acceptable overlap before outcome evaluation. Recommended default: at least 80% of S2-50 unless assay limitations justify another threshold.

### S2-280
Secondary panel. Use the intersection with the frozen 280-gene shared assay panel. Do not call this whole-transcriptome validation.

### S2-250
Optional if overlap is high and computation is inexpensive.

## External preprocessing
- Use the frozen C5 preprocessing/model contract whenever technically compatible.
- Apply frozen encoder preprocessing to external images.
- Cache external embeddings separately.
- Verify barcode/location alignment.
- Any platform-specific deterministic preprocessing must be defined before outcome evaluation.
- Do not choose preprocessing alternatives using external test performance.
- If the external platform is incompatible with the frozen target representation, stop and document the incompatibility.

## Freeze before external evaluation
Create and hash:

`C7/C7A/C7A_SPEC_LOCK.json`

Freeze:
- external cohort(s)
- patient/sample IDs
- gene-panel intersections
- encoders
- reconstruction models
- mask regimes
- coverage levels
- random seeds
- preprocessing
- metrics
- pass criteria

Only after this lock is written may C7B evaluate external hidden labels.

## C7A status
Use one of:
- `PASS_EXTERNAL_DATASET_QUALIFICATION`
- `LIMITED_EXTERNAL_DATASET`
- `NO_SUITABLE_EXTERNAL_DATASET`
- `FAIL_TECHNICAL`

If no suitable external dataset exists, STOP and explain why.

---

# C7B — Locked external validation

## Primary question
Does the frozen C5 sparse-ST reconstruction approach improve prediction of hidden locations in an independent external cohort?

## Frozen method set
Evaluate at minimum:

### E0 — H&E only
Frozen image-only model contract.

### E1 — Sparse-ST interpolation only
Frozen/simple interpolation method from C5.

### E2 — H&E + sparse-ST fusion
Primary C5 method to validate externally.

### E3 — Learned propagation
Secondary only. Run only if the exact frozen C5 configuration is available and compatible. Do not retune externally.

## Sampling strategy
C7 is not another sampling-optimization checkpoint.

Use only simple, predeclared policies:
1. random sparse observation;
2. spatially regular / coverage-maximizing sampling.

Do not use adaptive sampling as the primary external result.

## Coverage levels
Use:
- 10%
- 20%
- 40%

Optionally include 5% if the external cohort has enough locations.

Do not add new coverage points after seeing outcomes.

## External mask regimes
Evaluate separately:

### V1 — Random/open-tissue interpolation
Measured anchors distributed across tissue.

### V2 — Contiguous hidden-region extrapolation
A compact tissue region is fully unmeasured.

### V3 — Local-to-global extrapolation
Optional if tissue size permits: a compact observed region is measured and the remainder is reconstructed.

Never combine V1/V2/V3 into one headline metric.

## Replication
For V1 random sampling, use at least 10 random seeds per patient/budget where feasible.

For contiguous masks, use multiple predefined geometries/positions where tissue size permits. Do not cherry-pick favorable missing regions.

---

# Leakage requirements

For every external run:
1. remove hidden expression before model-input construction;
2. hidden labels cannot enter interpolation/fusion/graph inputs;
3. hidden labels cannot influence sampling;
4. hidden labels cannot influence preprocessing or model selection;
5. hidden labels are used only for final metrics;
6. save observed/hidden barcode lists;
7. save mask checksums.

Run corruption-invariance tests:
- replace hidden external expression with NaN/random values;
- predictions must remain unchanged.

Any failure = `FAIL_TECHNICAL`.

---

# Primary metrics

Compute only on hidden locations.

Per gene:
- Pearson
- Spearman
- MSE
- MAE

Per panel:
- mean/median Pearson
- mean/median Spearman
- mean/median MSE
- mean/median MAE
- fraction Pearson > 0
- fraction Pearson > 0.1
- fraction Pearson > 0.2
- fraction Pearson > 0.3

Report separately by:
- external patient
- budget
- mask regime
- sampling strategy
- reconstruction model

Do not treat thousands of spots as independent patients.

---

# Primary effect sizes

## Fusion gain over H&E only

`Delta_fusion_image = Metric(E2) - Metric(E0)`

## Fusion gain over interpolation only

`Delta_fusion_interp = Metric(E2) - Metric(E1)`

For error metrics, use the sign convention so positive means improvement.

The main publication story is strongest if sparse molecular anchors improve over H&E alone and H&E adds value beyond interpolation alone.

---

# Coverage-response replication

Reproduce externally:

`measured coverage -> hidden-location reconstruction accuracy`

Plot separately for:
- V1 interpolation
- V2 contiguous extrapolation

The direction of the relationship matters, not only the best budget.

---

# Distance-to-anchor replication

For each hidden external location compute distance to the nearest observed ST anchor.

Reproduce:

`error(distance_to_anchor)`

using the same predeclared binning logic as C5 where possible.

Primary question:

> Does reconstruction degrade with increasing distance from measured anchors in the external cohort as it did in development?

Do not force a single extrapolation-horizon number if the external data do not support one.

---

# Contiguous extrapolation requirement

C5/C6 showed this is the hard case. C7 must report it explicitly.

Show interpolation and contiguous extrapolation side-by-side.

A method may pass external interpolation validation while remaining weak for contiguous extrapolation. State this directly if observed.

---

# External patient consistency

For each major comparison report:
- per-patient direction of effect;
- number/fraction of patients where E2 beats E0;
- number/fraction where E2 beats E1;
- median patient-level effect.

---

# C7B status

Use one of:
- `PASS_EXTERNAL_VALIDATION`
- `PASS_LIMITED_EXTERNAL_VALIDATION`
- `MIXED_EXTERNAL_VALIDATION`
- `FAIL_EXTERNAL_VALIDATION`
- `FAIL_TECHNICAL`

## PASS_EXTERNAL_VALIDATION
Require:
1. independent cohort passes integrity audits;
2. E2 improves over H&E-only on primary S2-50;
3. E2 improves over interpolation-only for a meaningful subset of budgets/patients;
4. benefit is not confined to one external patient;
5. coverage-response behavior is coherent;
6. distance-to-anchor behavior is reproduced or clearly characterized;
7. contiguous extrapolation is honestly evaluated;
8. correlation gains are not systematically contradicted by worse MSE.

## PASS_LIMITED_EXTERNAL_VALIDATION
Use if the main interpolation result reproduces but patient count, panel overlap, or extrapolation evidence is limited.

## MIXED_EXTERNAL_VALIDATION
Use if only some budgets/patients replicate or effects are strongly heterogeneous.

Do not tune the frozen method inside C7 to rescue a mixed result.

## FAIL_EXTERNAL_VALIDATION
Use if the frozen method does not meaningfully reproduce after technical integrity passes.

---

# C7C — Secondary transportability stress test

Optional and only after C7B is frozen.

Possible tests:
- different breast-cancer ST platform;
- different preparation/staining source;
- different cancer/tissue.

Use the same frozen C7B method.

Purpose: characterize transportability, not rescue C7B.

Report C7C separately.

---

# Optional external adaptation analysis

Only after the primary zero-shot C7B result is frozen.

If severe domain shift is observed, a secondary analysis may test minimal adaptation, for example:
- unlabeled external H&E feature centering;
- explicit small observed-ST calibration subset;
- lightweight fusion-layer refit.

Rules:
1. zero-shot remains primary;
2. adaptation budget is predefined;
3. calibration locations are separated from hidden evaluation locations;
4. no hidden evaluation labels are used;
5. report adaptation separately.

---

# Negative controls

On at least one representative external patient/budget run:

## NC1 — Anchor-expression permutation
Permute observed anchor expression among measured locations.

## NC2 — Anchor-position permutation
Assign measured expression vectors to incorrect measured positions.

## NC3 — Image-location permutation
Permute external H&E embeddings across locations.

Corresponding information-specific gains should degrade.

---

# Required figures

At minimum:
1. external H&E with sampled anchors;
2. ground-truth / prediction / error maps for a predeclared S2-50 gene subset;
3. coverage-vs-Pearson curves;
4. coverage-vs-MSE curves;
5. interpolation vs contiguous-extrapolation comparison;
6. error vs distance-to-anchor;
7. E0 vs E1 vs E2 comparison;
8. patient-level effect plot;
9. development vs external summary plot.

Do not cherry-pick genes based on external performance.

---

# Development-versus-external comparison

Create:

`C7/C7_DEVELOPMENT_EXTERNAL_COMPARISON.csv`

with fields such as:

```text
dataset
patient
panel
budget
mask_regime
model
median_pearson
mse
mae
fusion_gain_vs_image
fusion_gain_vs_interpolation
```

The goal is to make domain shift visible, not hide it.

---

# Reproducibility requirements

For every external run save:
- dataset/source ID;
- patient/sample ID;
- raw asset hashes where practical;
- gene intersection;
- preprocessing config;
- encoder-cache hash;
- mask ID/checksum;
- observed barcode list;
- hidden barcode list;
- model config/checkpoint hash;
- seed;
- predictions;
- runtime;
- peak RAM.

Replay predictions where practical.

Verify:
- corruption invariance;
- no hidden-label input;
- no external-label model selection;
- frozen C0–C6 hashes unchanged.

---

# Required C7 files

```text
C7/
├── C7_REPORT.md
├── C7_STATUS.json
├── C7_IMMEDIATE_SUMMARY.txt
├── C7_DEVELOPMENT_EXTERNAL_COMPARISON.csv
├── C7A/
│   ├── C7A_REPORT.md
│   ├── C7A_SPEC_LOCK.json
│   ├── C7A_EXTERNAL_DATASET_MANIFEST.csv
│   ├── C7A_GENE_INTERSECTION.csv
│   └── logs/
├── C7B/
│   ├── C7B_REPORT.md
│   ├── C7B_PANEL_METRICS.csv
│   ├── C7B_GENE_METRICS.csv
│   ├── C7B_COVERAGE_CURVES.csv
│   ├── C7B_DISTANCE_TO_ANCHOR.csv
│   ├── C7B_PATIENT_EFFECTS.csv
│   ├── C7B_NEGATIVE_CONTROLS.csv
│   ├── masks/
│   ├── predictions/
│   ├── logs/
│   └── figures/
└── C7C/                 # optional
    ├── C7C_REPORT.md
    ├── metrics/
    ├── predictions/
    └── figures/
```

Append all runs to:

`results/experiment_registry.csv`

Do not overwrite prior rows.

---

# REQUIRED: immediate-return text file

Codex must create:

`C7/C7_IMMEDIATE_SUMMARY.txt`

This file is mandatory even if C7 stops early.

Target length: approximately 40–80 lines.

It must contain these sections:

```text
C7 IMMEDIATE SUMMARY

1. STATUS
C7A: ...
C7B: ...
C7C: ...
Overall: ...

2. EXTERNAL DATASET
Dataset/study:
Patients:
Platform:
Gene overlap:
Why this qualifies as external:

3. PRIMARY RESULT
S2-50:
H&E only:
ST interpolation only:
H&E + sparse ST:
Best coverage:
Contiguous extrapolation:

4. GENERALIZATION CHECK
Did the main C5 finding replicate? YES / PARTLY / NO
Did coverage-response replicate? YES / PARTLY / NO
Did distance-to-anchor behavior replicate? YES / PARTLY / NO

5. LIMITATIONS
- ...
- ...

6. RECOMMENDATION
Exactly one of:
PROCEED_TO_MANUSCRIPT_BENCHMARKING
PROCEED_TO_BIOLOGICAL_INTERPRETATION
NEED_EXTERNAL_DATASET_EXPANSION
NEED_METHOD_REVIEW
STOP_TECHNICAL

7. FULL REPORT
Path to C7_REPORT.md:
Git commit:
```

## Codex final-response requirement

When C7 finishes, Codex must **return the contents of `C7_IMMEDIATE_SUMMARY.txt` directly in its final message**, or attach/provide the text file directly if supported.

Do not return only a GitHub link.

Do not require the user to open `C7_REPORT.md` first.

The user should be able to see the outcome immediately from this text file.

---

# C7 overall status

Use one of:
- `PASS_EXTERNAL_VALIDATION`
- `PASS_LIMITED_EXTERNAL_VALIDATION`
- `MIXED_EXTERNAL_VALIDATION`
- `FAIL_EXTERNAL_VALIDATION`
- `NO_SUITABLE_EXTERNAL_DATASET`
- `FAIL_TECHNICAL`

---

# Final report questions

`C7_REPORT.md` must answer:
1. Was a genuinely independent external cohort used?
2. Did sparse measured ST improve reconstruction over H&E alone?
3. Did H&E add value beyond ST interpolation alone?
4. Did the result reproduce across external patients?
5. Did coverage-response behavior reproduce?
6. Did distance-to-anchor behavior reproduce?
7. How much weaker was contiguous extrapolation than interpolation?
8. How large was development-to-external domain shift?
9. Does current evidence support the core publication story?
10. What should the next checkpoint be?

---

# Final recommendation

Choose exactly one:
- `PROCEED_TO_MANUSCRIPT_BENCHMARKING`
- `PROCEED_TO_BIOLOGICAL_INTERPRETATION`
- `NEED_EXTERNAL_DATASET_EXPANSION`
- `NEED_METHOD_REVIEW`
- `STOP_TECHNICAL`

---

# Explicitly forbidden in C7

Do not:
- modify C0–C6;
- tune hyperparameters on external hidden labels;
- change gene panels after external outcomes are known;
- add new image encoders;
- fine-tune encoders on external labels;
- optimize sampling on external outcomes;
- revive multiresolution graph modeling;
- use external test residuals for acquisition;
- hide weak contiguous-extrapolation results;
- call S2-280 whole transcriptome;
- make clinical-utility or monetary-savings claims;
- treat external spots as independent patients;
- proceed automatically beyond C7.

---

# End-of-checkpoint instruction

Before external outcome evaluation:
1. identify and audit the external cohort;
2. freeze `C7A_SPEC_LOCK.json`;
3. freeze gene intersections;
4. freeze models;
5. freeze coverage levels;
6. freeze mask regimes;
7. freeze seeds and metrics;
8. verify no external hidden labels were used for selection.

After execution:
9. complete corruption/leakage audits;
10. verify frozen C0–C6 hashes;
11. write C7A/C7B/(optional C7C) reports;
12. write `C7_REPORT.md`;
13. write `C7_STATUS.json`;
14. create `C7_IMMEDIATE_SUMMARY.txt`;
15. commit/push C7;
16. **return `C7_IMMEDIATE_SUMMARY.txt` directly to the user in the final Codex response**;
17. STOP.

Do not proceed automatically to manuscript writing, biological interpretation, or any later checkpoint.
