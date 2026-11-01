# C6 Execution Plan — Sampling Optimization for Sparse Spatial Transcriptomics

## Repository

All C6 code, configs, policies, traces, metrics, figures, and reports must be created under:

`/home/son/Research/GraphPathology`

C0–C5 are frozen and must not be modified.

---

# Frozen prior checkpoint context

## C5

Status: `PASS_SPARSE_ST_RECONSTRUCTION`

Frozen source/report commit:

`e8b27529c4563785597ea9075c14d32121b602db`

C5 established that sparse measured spatial transcriptomics can support reconstruction of hidden locations in this four-patient development cohort.

Important C5 findings:

- sparse-ST reconstruction passed for S2-50 and S2-280;
- all leakage/corruption checks passed;
- contiguous extrapolation is harder than interpolation;
- learned-propagation gains depend on patient and mask geometry;
- C0–C4 remained unchanged.

C6 must build directly on the frozen C5 task and ask a new question:

> Given a fixed measurement budget, where should spatial transcriptomics be measured to maximize reconstruction fidelity over the unmeasured tissue?

C6 is still a **development-only** checkpoint.

No population-level or cost-saving claim is allowed from four patients.

---

# C6 scientific objective

C6 evaluates whether **measurement placement** matters after the total measured-ST budget is fixed.

The core idea is:

\[
\boxed{
\text{fixed assay budget}
\rightarrow
\text{choose measurement locations}
\rightarrow
\text{reconstruct the rest}
}
\]

C6 should answer:

1. Does sampling strategy matter at fixed coverage?
2. Which non-adaptive strategy is strongest?
3. Does H&E-guided sampling improve over random or regular spatial coverage?
4. If strategy matters, does adaptive sequential sampling improve further?
5. Can better sampling reduce errors in regions far from measured anchors?
6. Does optimized sampling help the difficult contiguous-extrapolation regime?

---

# C6 structure

C6 has two gated stages:

```text
C6A — Fixed-budget sampling strategy comparison
        ↓ only if sampling placement matters
C6B — Adaptive sequential sampling
```

Do not run C6B unless C6A shows a meaningful strategy effect.

---

# Global C6 rules

1. Keep C0–C5 frozen.
2. Use the exact C5 target panels, normalization, and hidden-location reconstruction definition.
3. Do not introduce new image encoders.
4. Do not fine-tune ResNet18 or Phikon.
5. Do not introduce multiresolution H&E-only graph models.
6. Do not alter the definition of observed versus hidden locations.
7. Hidden labels must never be used to choose sampling locations.
8. All sampling-policy hyperparameters must be selected using training/validation data only.
9. Primary test hidden labels may be accessed only after the full sampling trajectory is frozen.
10. Sampling strategies may use:
    - H&E/image embeddings;
    - spatial coordinates;
    - currently observed ST values;
    - model predictions;
    - model uncertainty estimated without hidden labels.
11. Sampling strategies may **not** use:
    - hidden gene expression;
    - test error maps;
    - test residuals;
    - test-derived gene rankings;
    - oracle knowledge of future prediction error.
12. Every sampling trajectory must be saved exactly and replayable.
13. Every strategy comparison must use the same budget.
14. C6 must distinguish pointwise interpolation from contiguous extrapolation.
15. Stop after C6 and wait for review.

---

# Reconstruction model used in C6

C6 is about **sampling policy**, not inventing a new reconstruction architecture.

Use the frozen C5 reconstruction pipeline.

## Primary reconstruction model

Use the C5 primary-selected model/configuration according to C5's own frozen selection procedure.

Do not choose a new model because it happened to perform best on the C5 primary test patient.

## Required baseline reconstruction

Also retain one simple C5 baseline:

- sparse-ST interpolation only, or
- the strongest simple C5 fusion baseline selected without test leakage.

This allows separation of:

\[
\text{sampling-policy gain}
\]

from:

\[
\text{reconstruction-model complexity}.
\]

If C5 has multiple frozen image representations, keep the C5 predeclared primary representation(s); do not start a new encoder search in C6.

---

# Gene-panel settings

Prioritize:

- **S2-50**
- **S2-280**

S2-250 is optional if computationally inexpensive.

Do not call S2-280 whole-transcriptome reconstruction.

---

# Measurement budgets

Evaluate fixed realized measurement fractions:

```text
5%
10%
20%
40%
```

Include 60% only if useful for saturation analysis.

For each patient and strategy record:

- nominal budget;
- realized budget;
- number of measured locations;
- number of hidden locations.

All strategies compared at a given budget must use the same number of measured locations, up to unavoidable geometric constraints.

---

# C6A — Fixed-budget sampling strategy comparison

## Scientific question

> At the same measurement budget, does the choice of measured locations materially affect reconstruction of the remaining tissue?

At minimum compare the following strategies.

---

## A0 — Random sampling

Uniformly sample measurement locations from eligible tissue locations.

Use at least:

`10 random seeds`

per patient and budget if computationally feasible.

This is the primary stochastic baseline.

Report the full distribution, not only the mean.

---

## A1 — Spatially regular / coverage-maximizing sampling

Use a deterministic spatial strategy designed to distribute anchors across the tissue.

Possible implementations:

- farthest-point sampling in physical coordinate space;
- grid-like selection;
- k-center coverage in coordinate space.

Prefer a method that:

- minimizes large unmeasured spatial gaps;
- does not use expression labels;
- is deterministic once initialized.

Record the exact algorithm.

---

## A2 — Morphology-diverse sampling

Select locations to cover diverse H&E morphology.

Use the frozen image embedding space.

Acceptable approaches:

- k-medoids on H&E embeddings;
- farthest-point sampling in H&E embedding space;
- diversity maximization.

Do not use expression to form clusters.

Purpose:

> test whether measuring morphologically diverse tissue improves molecular reconstruction at fixed budget.

---

## A3 — Joint spatial + morphology coverage

Combine physical-space and H&E-space diversity.

Example:

\[
D(i,j)
=
\alpha D_{\text{spatial}}(i,j)
+
(1-\alpha)D_{\text{morph}}(i,j)
\]

Use only a small predefined validation-only alpha grid.

Freeze alpha before test evaluation.

Purpose:

> avoid selecting many morphologically distinct points from the same spatial region or many spatially dispersed points with redundant morphology.

---

## A4 — Uncertainty-guided one-shot sampling

Optional in C6A.

Use model uncertainty estimated **without hidden labels** to select a fixed set before measuring those locations.

Examples:

- ensemble predictive variance;
- seed variance;
- distance-to-training-distribution;
- calibrated uncertainty available from the frozen C5 model.

Do not derive uncertainty from test residuals.

If the C5 model does not provide a defensible uncertainty estimate, skip A4 rather than inventing one.

---

# Sampling-policy fairness

At each patient/budget:

- same eligible tissue locations;
- same measured-location count;
- same reconstruction model;
- same gene panel;
- same downstream evaluation set definition;
- no hidden-label access.

For random strategies, report mean ± SD and the empirical distribution.

For deterministic strategies, report the single frozen trajectory.

---

# C6A geometry diagnostics

For every selected measurement set compute:

- mean distance from each hidden point to nearest measured anchor;
- median nearest-anchor distance;
- 90th percentile nearest-anchor distance;
- maximum nearest-anchor distance;
- spatial coverage radius;
- number of large unmeasured connected regions if definable;
- morphology-embedding coverage/diversity score.

These diagnostics are explanatory, not success metrics by themselves.

---

# C6A reconstruction metrics

Evaluate only on hidden locations.

Per panel:

- median gene-wise Pearson;
- mean gene-wise Pearson;
- MSE;
- MAE;
- fraction of genes above Pearson thresholds used in C5.

Report separately by:

- patient;
- budget;
- sampling strategy;
- reconstruction model;
- mask/extrapolation regime where applicable.

---

# Primary C6A effect size

For strategy \(s\) at budget \(b\):

\[
\Delta^{sample}_{s,b}
=
Metric(s,b)
-
Metric(random,b)
\]

Use random-sampling mean as the primary reference.

For error metrics, use the corresponding signed convention so that positive gain always means improvement.

Report both:

- correlation gain;
- MSE change.

A strategy is not considered superior if correlation improves while MSE systematically worsens.

---

# Coverage-efficiency summary

For every strategy, construct:

\[
A_s(b)
=
\text{reconstruction performance at budget } b.
\]

Report the area under the budget-performance curve:

**AUC-coverage**

for descriptive comparison.

Use the same budget points for all strategies.

Also report raw curves; do not reduce the result to AUC alone.

---

# Contiguous-extrapolation stress test

C5 showed contiguous extrapolation is harder.

C6A must explicitly test whether optimized placement reduces this weakness.

At minimum define one evaluation where:

- a contiguous tissue region remains unmeasured;
- sampling occurs outside or around that region;
- all strategies have the same total measurement budget.

Compare reconstruction inside the held-out contiguous region.

Important:

> Do not let the sampling strategy place anchors inside a region defined as unavailable for the extrapolation stress test.

This tests whether intelligent placement around a missing region improves extrapolation.

---

# Distance-to-anchor reanalysis

Repeat the C5 distance-to-anchor analysis under each sampling strategy.

For hidden location \(i\):

\[
d_i
=
\min_{j \in \mathcal O}
\|p_i-p_j\|.
\]

Plot:

\[
\text{error}(d)
\]

for each strategy.

Key question:

> Does optimized sampling improve reconstruction primarily by reducing anchor distance, or does it improve prediction even at matched distance?

Where feasible, compare strategies within matched distance bins.

---

# C6A patient sensitivity

Use all four patients in the same conditional development framework used in prior checkpoints.

Sampling-policy hyperparameters must be frozen from training/validation patients.

For a held-out patient, the policy may use that patient's:

- H&E;
- coordinates;
- currently measured ST anchors;

but may not use hidden ST values.

Do not describe this as unbiased external validation.

---

# C6A negative controls

At minimum include:

## NC1 — Random strategy distribution

Use enough random seeds to establish the natural variance of random sampling.

## NC2 — Morphology-permuted sampling

For the morphology-guided strategy, permute image embeddings across locations before selecting anchors.

If morphology-guided sampling truly matters, its advantage should decrease.

## NC3 — Spatial-coordinate permutation

For the spatial-coverage strategy, permute coordinates before selection while preserving reconstruction evaluation on the true coordinates.

Expected: coverage advantage should collapse.

Run controls on representative budgets such as 10% and 20%.

---

# C6A gate

Assign one of:

- `PASS_SAMPLING_EFFECT`
- `MIXED_SAMPLING_EFFECT`
- `NO_MEANINGFUL_SAMPLING_EFFECT`
- `FAIL_TECHNICAL`

## PASS_SAMPLING_EFFECT

Use if at least one non-random strategy:

1. reproducibly improves hidden-location reconstruction over random sampling at matched budget;
2. shows benefit on more than one patient;
3. does not systematically worsen MSE;
4. remains beneficial in at least one difficult contiguous-extrapolation setting;
5. passes leakage and policy-replay audits.

If PASS, proceed to C6B.

## MIXED_SAMPLING_EFFECT

Use if:

- strategy effects exist but are small, patient-dependent, budget-dependent, or geometry-dependent.

STOP for review.

Do not automatically run adaptive sampling.

## NO_MEANINGFUL_SAMPLING_EFFECT

Use if random/regular strategies are essentially equivalent within development variability.

STOP C6.

This does **not** invalidate the C5 sparse-ST reconstruction result.

Recommend external validation of the simplest defensible strategy.

---

# C6B — Adaptive sequential sampling

Run only if:

`C6A == PASS_SAMPLING_EFFECT`

## Scientific question

> Can sequentially choosing the next ST measurements based on current observations improve reconstruction efficiency beyond the best one-shot strategy?

---

# Adaptive protocol

Start from a small initial measured set.

Recommended:

```text
initial budget = 5%
```

Then increase sequentially:

```text
5% → 10% → 20% → 40%
```

At each step:

1. reconstruct the unmeasured tissue;
2. compute allowed acquisition scores using only currently available information;
3. select new locations;
4. reveal their ST expression;
5. update/retrain the reconstruction model according to a frozen protocol;
6. continue to the next budget.

Save every selected location in order.

---

# C6B adaptive strategies

At minimum compare:

## B0 — Random sequential acquisition

Randomly add new measured locations at each step.

Use multiple seeds.

## B1 — Best C6A non-adaptive strategy applied incrementally

Example:

- spatial coverage;
- morphology diversity;
- joint spatial+morphology strategy.

## B2 — Uncertainty-guided adaptive sampling

Select locations with highest current predictive uncertainty, optionally with a diversity constraint.

Allowed uncertainty signals must be label-free for unmeasured locations.

## B3 — Hybrid uncertainty + diversity

Optional.

Combine uncertainty with spatial/morphological diversity to prevent clustering all new measurements in one region.

Keep the scoring rule simple and validation-frozen.

---

# Adaptive update rule

Do not allow different strategies different model-update budgets.

For every acquisition step, use the same frozen update protocol, e.g.:

- refit only the lightweight fusion layer; or
- retrain the frozen C5 reconstruction model from scratch using currently observed anchors; or
- a fixed number of update epochs.

Choose one protocol on validation and freeze it.

Do not tune update effort separately by strategy.

---

# C6B metrics

At every acquisition budget record:

- median Pearson;
- MSE;
- MAE;
- distance-to-anchor statistics;
- runtime;
- selected-location trace.

Compute:

- AUC-coverage;
- performance at 10%, 20%, and 40%;
- improvement over random sequential acquisition;
- improvement over the best C6A one-shot strategy.

---

# Early-efficiency metric

Optionally define:

\[
B^\star
=
\text{minimum coverage required to reach a predeclared fidelity target}.
\]

Only report \(B^\star\) if the fidelity threshold was frozen before final test evaluation.

Otherwise report continuous budget-performance curves.

Do not claim monetary savings.

---

# C6B status

Assign one of:

- `PASS_ADAPTIVE_SAMPLING`
- `MIXED_ADAPTIVE_SAMPLING`
- `NO_ADAPTIVE_ADVANTAGE`
- `FAIL_TECHNICAL`

## PASS_ADAPTIVE_SAMPLING

Use only if adaptive sampling:

- beats random sequential acquisition;
- beats or complements the strongest C6A strategy;
- is reproducible across multiple patients;
- shows coherent benefit at low or moderate budgets;
- passes trace/replay/leakage audits.

Adaptive sampling is **not required** for the overall publication direction to remain viable.

---

# Reproducibility and leakage audit

For every policy run save:

- source commit;
- patient ID;
- panel;
- budget;
- sampling strategy;
- strategy config;
- strategy seed;
- selected barcode list in exact order;
- observed-label reveal sequence;
- reconstruction model config;
- scaler hashes;
- prediction artifact hash;
- metrics;
- runtime.

For adaptive runs save a full trace:

```text
step
current_budget
candidate_scores
selected_locations
newly_revealed_locations
model_checkpoint
prediction_hash
```

Verify:

1. hidden labels cannot affect selection scores;
2. replacing unrevealed labels with NaN/random values leaves the next selected locations unchanged;
3. all sampling traces replay exactly;
4. C0–C5 frozen hashes remain unchanged.

A failed selection corruption-invariance test is `FAIL_TECHNICAL`.

---

# Required outputs

Create:

```text
C6/
├── C6_REPORT.md
├── C6_STATUS.json
├── C6A/
│   ├── C6A_REPORT.md
│   ├── C6A_SPEC_LOCK.json
│   ├── C6A_POLICY_MANIFEST.csv
│   ├── C6A_SELECTIONS.csv
│   ├── C6A_PANEL_METRICS.csv
│   ├── C6A_DISTANCE_METRICS.csv
│   ├── C6A_GEOMETRY_DIAGNOSTICS.csv
│   ├── C6A_PATIENT_SENSITIVITY.csv
│   ├── C6A_NEGATIVE_CONTROLS.csv
│   ├── policies/
│   ├── predictions/
│   ├── logs/
│   └── figures/
│
└── C6B/                         # only if C6A passes
    ├── C6B_REPORT.md
    ├── C6B_SPEC_LOCK.json
    ├── C6B_ADAPTIVE_TRACES.jsonl
    ├── C6B_PANEL_METRICS.csv
    ├── C6B_BUDGET_CURVES.csv
    ├── C6B_PATIENT_SENSITIVITY.csv
    ├── policies/
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

1. reconstruction performance vs measurement budget for all C6A strategies;
2. MSE vs measurement budget;
3. random-sampling distribution vs deterministic strategies;
4. selected measurement locations overlaid on H&E for representative patients;
5. nearest-anchor distance distributions by strategy;
6. error vs distance-to-anchor by strategy;
7. contiguous-extrapolation performance by strategy;
8. morphology-diversity/coverage diagnostics;
9. patient-resolved strategy comparison;
10. if C6B runs: adaptive trajectories and AUC-coverage comparison.

Use deterministic/predeclared examples where possible.

Do not cherry-pick only successful slides or genes.

---

# C6 overall status

Assign one of:

- `PASS_SAMPLING_OPTIMIZATION`
- `MIXED_SAMPLING_OPTIMIZATION`
- `NO_SAMPLING_ADVANTAGE`
- `FAIL_TECHNICAL`

A valid `PASS_SAMPLING_OPTIMIZATION` requires at minimum `PASS_SAMPLING_EFFECT` in C6A.

C6B adaptive success is optional.

---

# Final C6 report questions

`C6_REPORT.md` must explicitly answer:

1. Does measurement placement matter at fixed ST coverage?
2. Which strategy performs best without hidden-label access?
3. Does H&E-guided morphology diversity improve over random/grid sampling?
4. Does optimized placement reduce nearest-anchor distances and contiguous-extrapolation error?
5. Is any gain preserved across multiple patients?
6. If C6B ran, does adaptive acquisition beat the best one-shot strategy?
7. At what budgets are strategy differences largest?
8. Is the result strong enough to justify external validation?

Final recommendation must be exactly one of:

- `PROCEED_TO_EXTERNAL_VALIDATION`
- `NEED_SAMPLING_REVIEW`
- `STOP_TECHNICAL`

---

# Explicitly forbidden in C6

Do not:

- modify C0–C5;
- use hidden expression to choose sampling locations;
- derive acquisition scores from hidden test residuals;
- introduce new image encoders;
- fine-tune ResNet18 or Phikon;
- introduce multiresolution H&E-only graph models;
- run GNNExplainer;
- run pathway enrichment;
- claim clinical utility;
- claim real monetary savings;
- claim unbiased generalization;
- optimize policies directly on the primary test patient;
- proceed automatically to external validation.

---

# Decision tree

```text
C5 = PASS_SPARSE_ST_RECONSTRUCTION
        ↓
C6A: Does sampling placement matter?
        │
        ├── FAIL_TECHNICAL
        │     → repair only
        │
        ├── NO MEANINGFUL EFFECT
        │     → stop C6
        │     → external validation with simple policy
        │
        ├── MIXED
        │     → stop for review
        │
        └── PASS
              ↓
C6B: Does adaptive sampling add further value?
              │
              ├── NO
              │     → keep best C6A strategy
              │
              ├── MIXED
              │     → stop for review
              │
              └── PASS
                    → retain adaptive strategy
        ↓
Next major checkpoint:
EXTERNAL VALIDATION
```

---

# End-of-checkpoint instruction

Before execution:

1. freeze `C6A_SPEC_LOCK.json`;
2. freeze budgets, sampling policies, policy hyperparameter grids, random seeds, metrics, and extrapolation stress-test geometry;
3. verify no test hidden labels were used for policy selection.

After C6A:

4. write `C6A_REPORT.md`;
5. assign C6A status;
6. if C6A is not `PASS_SAMPLING_EFFECT`, STOP.

Only if C6A passes:

7. freeze `C6B_SPEC_LOCK.json`;
8. execute adaptive sampling;
9. write `C6B_REPORT.md`;
10. assign C6B status.

Then:

11. write `C6_REPORT.md`;
12. write `C6_STATUS.json`;
13. complete selection corruption-invariance and replay audits;
14. verify all frozen C0–C5 hashes;
15. state the exact recommendation for external validation;
16. STOP.

Do not proceed automatically to the next checkpoint.
