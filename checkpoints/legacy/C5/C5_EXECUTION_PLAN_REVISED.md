# C5 Execution Plan — Sparse Spatial Transcriptomics Reconstruction

## Repository

All C5 code and artifacts must be created under:

`/home/son/Research/GraphPathology`

C0–C4 are frozen and must not be modified.

## Frozen context

- **C0:** `PASS_DEVELOPMENT_ONLY`; 4 HEST IDC patients, 35,536 paired locations, 280 harmonized genes, frozen 50/250/280-gene panels, frozen masks, frozen splits/preprocessing.
- **C1:** `PASS_DEVELOPMENT`; local H&E morphology contains measurable expression signal.
- **C2:** `MIXED_CONTEXT_SIGNAL`.
- **C3:** `MIXED_SPATIAL_SPECIFICITY`; true topology did not reproducibly beat matched controls.
- **C4:** `C4A PASS` only for S2-50 with Phikon; `C4B MIXED_SPATIAL_RECHECK`; larger panels did not improve consistently.

Therefore, C5 must **not** attempt to rescue the H&E-only graph-topology story. The primary task is sparse spatial-transcriptomics reconstruction.

---

# Scientific objective

C5 asks:

> If H&E is available everywhere but spatial transcriptomics is measured only at a subset of tissue locations, how accurately can expression at the remaining locations be reconstructed?

C5 must quantify:

1. reconstruction accuracy versus measured-ST coverage;
2. reconstruction accuracy versus distance to the nearest measured ST anchor;
3. H&E-only versus sparse-ST-only versus H&E+sparse-ST;
4. simple interpolation versus learned propagation;
5. interpolation versus true contiguous-region extrapolation;
6. gene-specific measurement requirements.

The prediction target is the **entire hidden location**, not individual hidden genes within an otherwise measured location.

---

# Core task definition

For every tissue/sample, partition locations into:

- `O`: observed ST locations
- `H`: hidden ST locations

For each hidden location `i in H`:

- H&E features remain available;
- coordinates remain available;
- the **entire gene-expression vector must be hidden**;
- no component of hidden expression may be used in model input, scaling, interpolation anchors, retrieval banks, contrastive positives, or auxiliary targets.

Predict:

`Y_hat_i in R^G`

for hidden locations only.

---

# Gene-panel settings

Run:

- **S2-50**
- **S2-250**
- **S2-280** = full shared assay panel

S1 may be used for debugging/visualization only.

Do not call S2-280 whole-transcriptome reconstruction.

---

# Image representations

Carry forward both frozen representations:

- ResNet18
- Phikon

Do not introduce new encoders.

Do not fine-tune either encoder.

Do not choose one based on C4 test results and discard the other.

---

# Mask regimes

Use the 216 frozen C0 masks where compatible, but classify them by geometry and realized coverage.

Evaluate at least:

## M1 — Random sparse observation

Measured locations are distributed across tissue.

Purpose: interpolation-dominated setting.

## M2 — Contiguous hidden region

A compact square/circular region is completely unmeasured.

Purpose: extrapolation into a missing tissue region.

## M3 — Contiguous observed region

Only a compact region is measured and the rest must be reconstructed.

Purpose: local-to-global extrapolation and strongest cost-stress setting.

If C0 masks do not cover all regimes, create separately versioned C5 masks. Do not modify C0.

For every mask save:

`mask_id, patient_id, mask_type, nominal_coverage, realized_coverage, num_observed, num_hidden, geometry, seed, checksum`

---

# Coverage levels

Where feasible evaluate approximately:

- 5%
- 10%
- 20%
- 40%
- 60%

measured-ST coverage.

Always record **realized coverage** exactly.

Do not merge interpolation and extrapolation masks into a single summary.

---

# Leakage rules — mandatory

For every mask:

1. remove hidden expression before constructing any molecular input;
2. hidden locations must not appear in interpolation anchors, retrieval banks, contrastive molecular pairs, expression-neighbor inputs, expression-derived scaling fits, or auxiliary molecular targets;
3. test hidden labels may be accessed only for final metric computation;
4. save observed/hidden barcode lists and mask checksum with each prediction;
5. run hidden-label corruption-invariance tests by replacing hidden labels with NaN/random values and verifying predictions are unchanged.

Any failed corruption-invariance test = `FAIL_TECHNICAL`.

---

# Baseline hierarchy

## B0 — Training mean

Predict training-set mean expression.

## B1 — H&E only

Use the frozen local predictor with:

- ResNet18
- Phikon

No measured ST anchors.

## B2 — Sparse ST interpolation only

Use observed ST anchors, but no H&E.

### B2a — kNN interpolation
Use a small predefined validation-only `k` grid.

### B2b — inverse-distance weighting
Use a small predefined validation-only exponent grid.

Optional Gaussian-kernel interpolation is allowed if frozen before test evaluation.

Do not begin with complex geostatistical methods.

## B3 — H&E + simple sparse-ST fusion

Combine image prediction and ST interpolation, e.g.:

`Y_hat = lambda * Y_hat_image + (1-lambda) * Y_hat_interp`

Use a small validation-only lambda grid.

Also allow one small learned fusion model using inputs such as:

- image embedding or image-only prediction;
- local anchor-expression summary;
- distance to nearest anchor;
- local anchor count.

Keep it small.

## B4 — Learned sparse-ST propagation

Only after B2/B3 are established.

Every location has image features.

Observed locations additionally have measured expression.

Hidden locations have no expression input.

Recommended node input:

`[image_embedding, observed_indicator, observed_indicator * expression]`

The observation indicator is mandatory.

Keep the graph/message-passing architecture compact.

Do not introduce multiresolution graphs.

---

# Central comparisons

The primary comparisons are:

- **B1 vs B2:** morphology-only versus sparse molecular interpolation
- **B2 vs B3:** does H&E add information beyond sparse ST?
- **B3 vs B4:** does learned propagation beat simple fusion?

Do not claim graph value merely because B4 beats B1.

---

# Primary metrics

Compute metrics **only on hidden locations**.

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
- fraction Pearson > 0, >0.1, >0.2, >0.3

Report separately by:

- patient
- mask
- coverage
- geometry
- encoder
- model

Do not treat all spots as independent biological replicates.

---

# Coverage-response analysis

For every model and mask regime, plot:

- median Pearson versus realized measured coverage
- MSE versus realized measured coverage

Keep separate curves for:

- random interpolation
- contiguous hidden-region extrapolation
- contiguous observed-region local-to-global extrapolation

This analysis is mandatory.

---

# Distance-to-anchor analysis

For each hidden location compute:

`d_i = distance to nearest observed ST anchor`

Use physical units when reliable; otherwise label coordinate units explicitly.

Use predefined distance bins.

For each bin report:

- MSE
- MAE
- Pearson where meaningful
- number of hidden locations

Generate error-versus-distance curves.

This analysis is mandatory.

---

# Effective extrapolation horizon

Define a development-only **effective extrapolation horizon** only if a fidelity criterion can be frozen before final test evaluation.

If a robust threshold is not justified, report the continuous error-versus-distance curve and do not force a single horizon value.

Do not interpret this as a biological diffusion length.

---

# Gene-specific measurement requirement

For each gene, evaluate performance across coverage levels.

Optionally define:

`C*_g = minimum measured-ST coverage required to reach a predefined fidelity threshold`

only if:

- the threshold is frozen before test evaluation;
- the result is stable across masks/seeds;
- patient sensitivity is acceptable.

Otherwise report the full gene-by-coverage matrix.

No pathway enrichment in C5.

---

# Interpolation versus extrapolation

Treat as separate tasks.

A model that succeeds with randomly distributed anchors but fails when an entire region is unmeasured must not be described as generally reconstructing unmeasured tissue.

---

# H&E contribution

For each mask/coverage setting calculate the gain of B3 over B2.

This is the main measure of whether morphology adds value after sparse transcriptomics is already available.

Compute separately for ResNet18 and Phikon.

---

# Learned-propagation contribution

Calculate B4 minus B3.

A learned propagation gain matters only if it is reproducible across masks/patients and improves both correlation and acceptable error behavior.

If B4 fails to beat B3, do not increase graph complexity.

---

# Four-patient conditional sensitivity analysis

Use the same interpretation as previous checkpoints:

**conditional sensitivity analysis — not unbiased external validation**

For each held-out patient:

- train model parameters on the other three;
- fit fold-specific scalers on fold-training patients only;
- evaluate masks only within the held-out patient;
- keep primary-selected configurations fixed.

At minimum run:

- S2-50
- S2-280
- 10%, 20%, 40% coverage
- one interpolation regime
- one contiguous extrapolation regime

---

# Negative controls

For a representative S2-50 / ~20% coverage setting, run at least 3 seeds of:

## NC1 — Anchor-expression permutation
Permute observed expression vectors among observed anchor locations.

## NC2 — Anchor-position permutation
Keep observed expression vectors but assign them to different observed positions.

## NC3 — Image-location permutation
For image-using models, permute image features among locations while keeping sparse ST anchors fixed.

Expected outcome: the corresponding information-specific gain should collapse.

---

# Optional cross-modal auxiliary learning

Do **not** make this the primary C5 model.

Only if B3/B4 are promising, optionally test:

`L = L_hidden_expression + lambda_c L_contrastive + lambda_m L_masked_reconstruction`

Strict rules:

- hidden test expression cannot enter contrastive/retrieval target banks;
- auxiliary molecular targets must come only from permitted observed/training locations;
- compare directly against the same model without auxiliary losses;
- label as exploratory.

Do not call it CMRCNet unless reproducing the official method faithfully.

---

# C5 status

Assign one of:

- `PASS_SPARSE_ST_RECONSTRUCTION`
- `MIXED_SPARSE_ST_RECONSTRUCTION`
- `FAIL_SPARSE_ST_RECONSTRUCTION`
- `FAIL_TECHNICAL`

## PASS

Require all of:

1. sparse measured ST improves hidden-location prediction over H&E-only models;
2. H&E+sparse-ST improves over sparse-ST interpolation for a meaningful subset of settings;
3. gains reproduce across multiple masks and more than one patient;
4. contiguous-region extrapolation shows nontrivial reconstruction ability;
5. corruption-invariance/leakage audits pass;
6. coverage-response behavior is coherent;
7. distance-to-anchor degradation is quantified.

Learned graph propagation is **not required** for PASS.

## MIXED

Use if:

- sparse anchors help but H&E adds little;
- H&E helps only in some patients/masks;
- interpolation succeeds but extrapolation is weak;
- correlation and error disagree substantially.

STOP for review.

## FAIL

Use if sparse measured ST cannot reliably improve hidden-location reconstruction beyond simple baselines after integrity checks pass.

Do not add architectural complexity merely to force a positive result.

---

# Required outputs

Create:

```text
C5/
├── C5_REPORT.md
├── C5_STATUS.json
├── C5_SPEC_LOCK.json
├── C5_MASK_INDEX.csv
├── C5_PANEL_METRICS.csv
├── C5_GENE_METRICS.csv
├── C5_COVERAGE_CURVES.csv
├── C5_DISTANCE_TO_ANCHOR.csv
├── C5_PATIENT_SENSITIVITY.csv
├── C5_NEGATIVE_CONTROLS.csv
├── C5_RESOURCE_PROFILE.csv
├── masks/
├── configs/
├── models/
├── predictions/
├── logs/
└── figures/
```

Append all runs to:

`results/experiment_registry.csv`

Do not overwrite prior rows.

---

# Required figures

At minimum:

1. measured coverage vs median Pearson;
2. measured coverage vs MSE;
3. separate interpolation and extrapolation curves;
4. error vs distance to nearest measured anchor;
5. B1 vs B2 vs B3 vs B4;
6. ResNet18 vs Phikon under sparse-ST fusion;
7. representative spatial reconstruction maps:
   - observation mask
   - ground truth
   - prediction
   - absolute error
8. gene-by-coverage heatmap;
9. patient-resolved summary.

Use deterministic/predeclared examples. Do not cherry-pick favorable maps.

---

# Reproducibility requirements

For every run save:

- source commit
- mask ID/checksum
- observed and hidden barcode lists
- encoder identity
- feature-cache hash
- model config
- scaler hashes
- seed
- checkpoint hash
- predictions
- runtime
- peak RAM

Replay saved models/predictions where practical.

Verify:

- hidden-label corruption invariance
- barcode-order consistency
- no hidden expression in model inputs
- C0–C4 frozen hashes unchanged

---

# C5 report questions

`C5_REPORT.md` must explicitly answer:

1. How much do sparse ST anchors improve over H&E alone?
2. Does H&E improve over sparse-ST interpolation alone?
3. Does learned propagation beat simple fusion?
4. How does performance change with measured-ST coverage?
5. How does error change with distance to the nearest observed anchor?
6. How different are interpolation and contiguous extrapolation?
7. Which genes appear to need dense versus sparse measurements?
8. Are results robust across patients and mask geometries?
9. Is a later sampling-optimization checkpoint justified?

Final recommendation must be exactly one of:

- `PROCEED_TO_SAMPLING_OPTIMIZATION`
- `PROCEED_TO_EXTERNAL_VALIDATION`
- `NEED_SPARSE_ST_REVIEW`
- `STOP_TECHNICAL`

---

# Explicitly forbidden in C5

Do not:

- modify C0–C4;
- reveal hidden-location expression to the model;
- mask individual genes instead of whole locations for the main task;
- use hidden locations in molecular contrastive/retrieval targets;
- introduce new image encoders;
- fine-tune ResNet18 or Phikon;
- introduce multiresolution H&E-only graphs;
- run pathway enrichment;
- run GNNExplainer;
- claim external validation;
- claim whole-transcriptome reconstruction;
- claim experimental cost savings from four patients;
- use primary test labels to choose coverage/fidelity thresholds.

---

# Decision tree

```text
C4 complete
   ↓
C5 sparse-ST reconstruction
   │
   ├── FAIL_TECHNICAL
   │     → repair only
   │
   ├── FAIL
   │     → stop methodological escalation
   │
   ├── MIXED
   │     → review interpolation vs extrapolation
   │
   └── PASS
         ↓
   Does performance depend meaningfully
   on coverage/location choice?
         │
         ├── YES → C6 sampling optimization
         └── NO  → consider external validation
```

---

# End-of-checkpoint instruction

Before execution:

1. freeze `C5_SPEC_LOCK.json`;
2. freeze mask categories, coverage bins, model candidate sets, metric definitions, and any fidelity criteria;
3. verify no primary-test hidden labels were used for selection.

After execution:

4. write `C5_REPORT.md`;
5. write `C5_STATUS.json`;
6. complete leakage and corruption-invariance audits;
7. verify all frozen C0–C4 hashes;
8. state the exact next-checkpoint recommendation;
9. STOP.

Do not proceed automatically to C6 or external validation.
