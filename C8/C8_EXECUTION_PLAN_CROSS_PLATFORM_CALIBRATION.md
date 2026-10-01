# C8 Execution Plan — Cross-Platform Target Harmonization and Fusion Calibration

## Repository
All C8 code and artifacts must be created under:

`/home/son/Research/GraphPathology`

C0–C7 are frozen and must not be modified.

## Frozen context

### C7
Status: `MIXED_EXTERNAL_VALIDATION`

External cohort:
- HER2ST (Andersson et al. 2021)
- 8 independent breast-cancer patients
- 3,455 locations
- 100-µm capture spots
- 40/50 overlap for the primary panel
- 226/280 overlap for the secondary panel

C7 established:
- sparse ST anchors generalize externally;
- H&E-only prediction is weak externally;
- frozen H&E + sparse-ST fusion improves over H&E alone;
- but fusion does not consistently improve over ST interpolation alone;
- at 20% coverage, fusion worsened MSE versus ST interpolation in all 8 patients;
- contiguous extrapolation remained substantially harder;
- cross-platform count scale, capture geometry, and image distribution are important domain-shift factors.

C7 recommendation: `NEED_METHOD_REVIEW`.

Important consequence:
HER2ST may be used in C8 for method development/calibration, but it must NOT be reused later as the pristine final external-validation cohort.

---

# Scientific objective

C8 asks:

> Can cross-platform target harmonization and fusion calibration make H&E contribute robustly beyond sparse-ST interpolation?

The key failure to address is:

`H&E + sparse ST` does not yet consistently beat `ST interpolation only` under platform shift.

C8 focuses on:
1. target representation/harmonization;
2. image-versus-ST fusion calibration;
3. reliability-aware fusion weighting;
4. simultaneous improvement in correlation and error.

C8 is not a new encoder checkpoint and not a graph-architecture checkpoint.

---

# C8 structure

```text
C8A — Target harmonization study
        ↓ only if promising
C8B — Fusion calibration study
        ↓ only if fusion beats interpolation
C8C — Internal transportability stress test
```

C8 uses HER2ST as a method-development transportability cohort.
Do not call C8 final external validation.

---

# Global rules

1. Keep C0–C7 unchanged.
2. Preserve C7 zero-shot results exactly.
3. Do not add new image encoders.
4. Do not fine-tune ResNet18 or Phikon.
5. Do not revive multiresolution H&E-only graph modeling.
6. Do not optimize sampling policies.
7. Patient identity is the biological split unit.
8. Any learned scaling/calibration must be fit only on allowed training/calibration patients.
9. Never use held-out-patient labels to fit target transforms, fusion weights, gates, or thresholds.
10. Report both correlation and error.
11. A method cannot pass by improving Pearson while systematically worsening MSE.
12. Stop after C8.

---

# Dataset roles

Use:
- original C0–C6 cohort as source-development data;
- HER2ST as transportability/calibration data.

Within HER2ST use patient-disjoint folds.

Recommended:
`4-fold patient-disjoint cross-validation`
with 2 held-out patients per fold.

Do not randomly split spots across patients.

---

# Primary gene panel

Primary:
- frozen 40-gene HER2ST intersection of original S2-50.

Secondary:
- frozen intersection with S2-280 if useful.

Do not replace missing genes.
Do not create a new panel based on C8 outcome.

---

# C8A — Target harmonization

## Question

Which target representation preserves spatial biological variation while reducing assay-specific scale mismatch?

Compare a small predefined set only.

### T0 — C7 frozen target representation
Reference baseline.

### T1 — Per-gene z-score
Fit mean/std on training patients only and apply unchanged to held-out patients.

### T2 — Per-gene robust scaling
Fit median/IQR on training patients only and apply unchanged.

### T3 — Rank/quantile spatial target
Optional but recommended.
Convert each gene's within-sample spatial expression to a rank/quantile-like target.
If used, clearly state that the task becomes relative spatial expression rather than absolute assay counts.

### T4 — Deterministic cross-platform log-normalized relative-expression target
Optional if technically justified.
No outcome-guided tuning of normalization constants.

For each transform, evaluate:
- H&E-only;
- ST interpolation-only;
- simple H&E + ST fusion.

Report:
- median/mean gene-wise Pearson;
- MSE;
- MAE;
- per-patient metrics;
- 10%, 20%, 40% coverage;
- contiguous extrapolation at 20%.

Also record target variance and cross-patient scale consistency.

---

# C8A status

Assign one of:
- `PASS_TARGET_HARMONIZATION`
- `MIXED_TARGET_HARMONIZATION`
- `NO_TARGET_HARMONIZATION_GAIN`
- `FAIL_TECHNICAL`

PASS only if at least one transform:
1. improves transportability vs T0;
2. improves/preserves Pearson while reducing error instability;
3. behaves consistently across multiple held-out HER2ST patients;
4. uses no held-out-patient fitting;
5. is not catastrophic under contiguous extrapolation.

Only PASS proceeds to C8B.

---

# C8B — Fusion calibration

## Question

Can H&E contribute conditionally rather than through one globally frozen fusion weight?

Primary baseline:
### F0 — ST interpolation only

Any fusion method must beat F0 to justify the H&E contribution.

### F1 — Global scalar fusion
`Y_hat = lambda * Y_hat_image + (1-lambda) * Y_hat_ST`

Choose lambda using training/calibration patients only.
Use a small frozen grid, e.g. 0.0, 0.1, ..., 1.0.

### F2 — Distance-aware fusion
Let image/ST weight depend on distance to nearest measured ST anchor.

Use only simple interpretable forms:
- clipped linear;
- monotonic sigmoid;
- small binned lookup.

Fit only on training/calibration patients.

### F3 — Interpolation-confidence-aware fusion
Allow fusion weight to depend on:
- nearest-anchor distance;
- kNN anchor distances;
- local anchor density;
- interpolation disagreement/variance.

Allowed gating model:
- linear/logistic;
- tiny MLP with strict parameter cap.

No hidden expression input.

### F4 — Image-confidence-aware fusion
Optional.
Use only label-free uncertainty such as ensemble/seed variance or distance to source embedding distribution.

### F5 — Gene-specific fusion weights
Optional secondary experiment.
Fit per-gene weights on training/calibration patients only with strong regularization.

Do not make this primary unless clearly justified.

---

# Fusion selection objective

Do not optimize Pearson alone.

Recommended frozen rule:
1. maximize median Pearson;
2. among models within a small predeclared tolerance of best Pearson, choose lowest MSE.

Any alternative composite objective must be frozen before held-out evaluation.

---

# C8B pass condition

Relative to ST interpolation only, require:
1. median Pearson improves;
2. MSE does not materially worsen;
3. majority of held-out HER2ST patients improve;
4. benefit appears at more than one coverage;
5. benefit is not restricted to one mask geometry.

Assign:
- `PASS_FUSION_CALIBRATION`
- `MIXED_FUSION_CALIBRATION`
- `NO_FUSION_ADVANTAGE`
- `FAIL_TECHNICAL`

---

# C8C — Internal transportability stress test

Run only if C8B passes.

Evaluate the frozen calibrated method across:
- all patient-disjoint HER2ST folds;
- 10%, 20%, 40% coverage;
- random/open-tissue interpolation;
- contiguous hidden-region extrapolation.

Report:
- patient-level effects;
- budget-level effects;
- gene-level effects;
- distance-to-anchor behavior.

---

# Development-to-HER2ST comparison

Create:

`C8/C8_DEVELOPMENT_TRANSPORTABILITY.csv`

Fields:
```text
dataset
patient
target_transform
fusion_method
coverage
mask_regime
median_pearson
mse
mae
gain_vs_ST_interp
gain_vs_H&E
```

Show residual domain shift explicitly.

---

# Negative controls

At minimum:

## NC1 — Random fusion weights
Replace calibrated weights with matched random weights.

## NC2 — Distance permutation
For distance-aware fusion, permute distance features across locations.

## NC3 — Image permutation
Permute image predictions/embeddings across locations while keeping ST interpolation fixed.

The corresponding fusion advantage should diminish.

---

# Leakage/calibration audit

For each outer fold:
1. define training/calibration patients;
2. fit target transform using allowed patients only;
3. fit fusion calibration using allowed patients only;
4. freeze;
5. evaluate held-out patients;
6. save all scaler/gating/weight objects and hashes.

Run corruption-invariance:
- replace held-out hidden expression with NaN/random values;
- predictions must remain unchanged.

Any failure => `FAIL_TECHNICAL`.

---

# Required outputs

```text
C8/
├── C8_REPORT.md
├── C8_STATUS.json
├── C8_IMMEDIATE_SUMMARY.txt
├── C8_DEVELOPMENT_TRANSPORTABILITY.csv
├── C8A/
│   ├── C8A_REPORT.md
│   ├── C8A_SPEC_LOCK.json
│   ├── C8A_TARGET_TRANSFORMS.csv
│   ├── C8A_PANEL_METRICS.csv
│   ├── C8A_PATIENT_METRICS.csv
│   ├── scalers/
│   ├── predictions/
│   ├── logs/
│   └── figures/
├── C8B/
│   ├── C8B_REPORT.md
│   ├── C8B_SPEC_LOCK.json
│   ├── C8B_FUSION_MODELS.csv
│   ├── C8B_PANEL_METRICS.csv
│   ├── C8B_GENE_METRICS.csv
│   ├── C8B_PATIENT_EFFECTS.csv
│   ├── calibration/
│   ├── predictions/
│   ├── logs/
│   └── figures/
└── C8C/   # only if C8B passes
    ├── C8C_REPORT.md
    ├── C8C_FOLD_METRICS.csv
    ├── predictions/
    └── figures/
```

Append all runs to:
`results/experiment_registry.csv`

Do not overwrite prior rows.

---

# Required figures

At minimum:
1. C7 frozen baseline vs C8 target transforms;
2. per-patient Pearson/MSE by target transform;
3. ST interpolation vs each fusion strategy;
4. fusion gain vs nearest-anchor distance;
5. 10/20/40% coverage-response curves;
6. interpolation vs contiguous extrapolation;
7. patient-level fusion effect;
8. gene-level fusion gain;
9. development vs HER2ST domain-gap summary.

Do not cherry-pick patients or genes.

---

# C8 overall status

Assign one of:
- `PASS_CROSS_PLATFORM_CALIBRATION`
- `MIXED_CROSS_PLATFORM_CALIBRATION`
- `NO_CROSS_PLATFORM_FUSION_ADVANTAGE`
- `FAIL_TECHNICAL`

PASS requires:
1. C8A passes;
2. C8B passes;
3. H&E+sparse ST beats ST interpolation on held-out HER2ST patients;
4. benefit exists at multiple coverages;
5. MSE is not systematically harmed;
6. contiguous extrapolation is reported honestly;
7. integrity audits pass.

---

# Mandatory immediate summary file

Codex must create:

`C8/C8_IMMEDIATE_SUMMARY.txt`

Keep it to about 40–80 lines.

Use exactly:

```text
C8 IMMEDIATE SUMMARY

1. STATUS
C8A: ...
C8B: ...
C8C: ...
Overall: ...

2. TARGET HARMONIZATION
Best transform:
Why selected:
Primary Pearson:
Primary MSE:
Improvement vs C7 target:

3. FUSION CALIBRATION
Best fusion:
ST interpolation Pearson / MSE:
Calibrated fusion Pearson / MSE:
Patient improvement count:
Coverage levels with benefit:

4. CONTIGUOUS EXTRAPOLATION
Result:
Main limitation:

5. DOMAIN-SHIFT RECOVERY
C7 frozen external result:
C8 calibrated result:
Did H&E add value beyond ST interpolation? YES / PARTLY / NO

6. LIMITATIONS
- ...
- ...

7. RECOMMENDATION
Exactly one of:
PROCEED_TO_FRESH_EXTERNAL_VALIDATION
NEED_MORE_TRANSPORTABILITY_DATA
NEED_METHOD_REVIEW
STOP_TECHNICAL

8. FULL REPORT
Path:
Git commit:
```

## Codex final-response requirement

When C8 finishes, Codex must return the contents of `C8_IMMEDIATE_SUMMARY.txt` directly in its final response.

Do not return only a GitHub link.

---

# Final recommendation

Choose exactly one:
- `PROCEED_TO_FRESH_EXTERNAL_VALIDATION`
- `NEED_MORE_TRANSPORTABILITY_DATA`
- `NEED_METHOD_REVIEW`
- `STOP_TECHNICAL`

Use `PROCEED_TO_FRESH_EXTERNAL_VALIDATION` only if the calibrated method consistently beats ST interpolation on held-out HER2ST patients.

The next validation cohort must be independent and unused in C0–C8.

HER2ST cannot be reused as the final validation cohort.

---

# Explicitly forbidden

Do not:
- modify C0–C7;
- overwrite C7 zero-shot results;
- use HER2ST held-out patients to fit scalers or fusion weights;
- choose target transforms from held-out outcomes;
- add new image encoders;
- fine-tune encoders;
- revive multiresolution graph modeling;
- optimize sampling policies;
- run GNNExplainer;
- run pathway enrichment;
- claim final external validation;
- claim clinical utility;
- claim monetary savings;
- call 226 genes whole transcriptome;
- proceed automatically to the next checkpoint.

---

# End-of-checkpoint instruction

Before C8 evaluation:
1. freeze the 40-gene primary panel;
2. freeze patient-disjoint folds;
3. freeze target transforms;
4. freeze fusion candidates;
5. freeze 10/20/40% coverages;
6. freeze mask regimes;
7. freeze selection objective and pass criteria;
8. write/hash `C8A_SPEC_LOCK.json`.

During execution:
9. never fit preprocessing/fusion calibration using held-out patient labels;
10. save all calibration objects and hashes;
11. run corruption-invariance tests.

After execution:
12. write `C8A_REPORT.md`;
13. if C8A passes, run C8B;
14. if C8B passes, optionally run C8C;
15. write `C8_REPORT.md`;
16. write `C8_STATUS.json`;
17. create `C8_IMMEDIATE_SUMMARY.txt`;
18. verify C0–C7 hashes unchanged;
19. commit/push C8;
20. return `C8_IMMEDIATE_SUMMARY.txt` directly to the user;
21. STOP.

Do not proceed automatically to fresh external validation.
