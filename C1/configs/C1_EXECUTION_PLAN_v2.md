# C1 Execution Plan — Morphology-Only Spatial Transcriptome Prediction

## Repository
All code, configs, logs, metrics, plots, and reports for this checkpoint must be created under:

`/home/son/Research/GraphPathology`

Do not modify or overwrite frozen C0 artifacts.

---

## C0 input status

C0 status: `PASS_DEVELOPMENT_ONLY`

Frozen source commit:

`2b54db3379500fec4603a61751145c929db9a72a`

Use the exact frozen C0 data, splits, masks, feature caches, target preprocessing, and hashes recorded in:

- `C0_SPEC_LOCK.json`
- `C0_DATA_MANIFEST.json` if present
- `C0_SPLITS.json` if present
- `C0_MASK_INDEX.json`
- `C0_MASKS/`

Do not silently regenerate or replace any C0 artifact.

C0 contains:

- 35,536 paired spatial locations
- 4 HEST IDC patients
- 280 harmonized biological genes
- 50-gene panel
- 250-gene panel
- frozen single-gene target: top training-only log-variance gene
- frozen ImageNet ResNet18 image embeddings
- 216 frozen observation masks

Known limitations:

- only 4 patients
- development-only cohort
- no 1000/3000-gene panel
- generic ImageNet ResNet18 features
- no validated cell/nucleus annotations
- simulated square/circular masks
- no cell-level registration claim

C1 must not make publication-level claims.

---

# C1 objective

Determine whether local H&E morphology alone contains recoverable signal for spatial gene-expression prediction before introducing graph message passing, multiresolution modeling, sparse-ST extrapolation, or pathology-specific foundation embeddings.

The central question is:

> Can frozen local image features predict spatial gene expression on held-out patients?

C1 is a morphology-only checkpoint.

Do **not** use spatial graph edges.

Do **not** use neighboring transcriptomic measurements.

Do **not** use the 216 sparse-observation masks for training.

Do **not** run S3 partial-ST completion in C1.

---

# Prediction settings

Run all supported settings.

## S1 — Single-gene spatial field

Use exactly the single gene frozen by C0:

> the top training-only log-expression-variance gene

Predict one expression value at every target spatial location.

Output spatial maps for held-out data.

---

## S2-50 — 50-gene joint prediction

Predict:

\[
f(x_i) \rightarrow \hat{y}_i \in \mathbb{R}^{50}
\]

Use the exact frozen C0 50-gene panel.

---

## S2-250 — 250-gene joint prediction

Predict:

\[
f(x_i) \rightarrow \hat{y}_i \in \mathbb{R}^{250}
\]

Use the exact frozen C0 250-gene panel.

---

## S2-280 — Full shared assay panel

Predict:

\[
f(x_i) \rightarrow \hat{y}_i \in \mathbb{R}^{280}
\]

Use all 280 harmonized biological targets.

Label this experiment:

**Full shared assay panel**

Do not call it whole-transcriptome prediction.

---

# Models

Use the frozen ResNet18 embeddings from C0.

No pathology foundation model is introduced in C1.

Run the following baselines.

## B0 — Training-set mean

For each gene:

\[
\hat{y}_g = \mathrm{mean}_{train}(y_g)
\]

This is the null baseline.

---

## B1 — Ridge regression

Train a multi-output Ridge model from frozen image embeddings to gene expression.

Requirements:

- training data only for fitting
- validation data only for selecting regularization strength
- test data untouched until final evaluation
- use a small, predefined alpha grid
- no large hyperparameter sweep

Recommended alpha grid:

`[1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100]`

Select alpha using validation-set median gene-wise Pearson correlation.

Tie-break using validation MSE.

---

## B2 — Small MLP

Use a deliberately small network.

Recommended default:

```text
Input embedding
→ Linear(d, 256)
→ GELU
→ Dropout(0.1)
→ Linear(256, output_dim)
```

Allowed tuning should remain minimal.

Candidate hidden dimensions:

`[128, 256]`

Candidate dropout:

`[0.0, 0.1]`

Use early stopping on validation loss.

Do not conduct an architecture search.

---

# Data-split rules

The frozen C0 train/validation/test split is the primary development split.

Patient separation must remain exact.

For the **primary frozen development split**:

- **Training patients** may be used to fit model parameters and any learned preprocessing/scaling.
- **Validation patients** may be used to select hyperparameters, choose among predefined candidate configurations, and drive early stopping.
- **Test patients must remain completely untouched until final evaluation.**

No sample, patch, barcode, or expression target from the **test patient(s)** may influence:

- normalization or learned scaling
- gene ranking
- model fitting
- hyperparameter selection
- architecture/configuration choice
- early stopping
- checkpoint selection
- feature scaling
- target scaling

Any learned feature scaling must be fit on the training patients only and then applied unchanged to validation and test patients.

Any learned target scaling/standardization must be fit on the training patients only and then applied unchanged to validation and test patients. Invert it before reporting metrics in the original normalized target space.

Validation data may guide model selection, but it must never be merged into training before the final test evaluation unless a separately specified post-selection refit protocol is explicitly frozen in advance. C1 does not use such a refit by default.

---

# Primary evaluation

For every gene \(g\), compute on held-out test data:

- Pearson correlation
- Spearman correlation
- MSE
- MAE
- target variance

Store one row per gene.

For S2-50, S2-250, and S2-280 report:

- mean Pearson
- median Pearson
- mean Spearman
- median Spearman
- mean MSE
- median MSE
- mean MAE
- median MAE
- fraction of genes with Pearson > 0
- fraction of genes with Pearson > 0.1
- fraction of genes with Pearson > 0.2
- fraction of genes with Pearson > 0.3

Do not report only one aggregate correlation.

---

# Patient-resolved evaluation

Report metrics separately for each held-out patient whenever the split permits.

Never treat 35,536 locations as 35,536 independent biological replicates.

The patient is the primary biological unit.

For pooled spot-level metrics, clearly label them as descriptive.

---

# Gene-gene structure preservation

For S2-50, S2-250, and S2-280, evaluate whether the predicted transcriptomic profiles preserve co-expression structure.

For measured test expression:

\[
C = \mathrm{corr}(Y)
\]

For predicted test expression:

\[
\hat{C} = \mathrm{corr}(\hat{Y})
\]

Compute:

\[
R_{coexpression}
=
\mathrm{corr}
\left(
\mathrm{vec}_{upper}(C),
\mathrm{vec}_{upper}(\hat{C})
\right)
\]

Report the result for each panel.

For S2-50, generate:

1. ground-truth gene-gene correlation heatmap
2. predicted gene-gene correlation heatmap
3. difference heatmap

Use the same gene ordering for all three.

---

# Spatial visualization

For S1, generate for the held-out test patient:

1. H&E reference
2. measured expression map
3. predicted expression map
4. absolute-error map

Use identical coordinates and spatial extent.

For S2, create the same four-panel visualization for a deterministic set of genes.

Do not cherry-pick genes using final test performance.

Use the following deterministic selection based only on training data:

- highest training variance gene
- median training variance gene
- lowest nonzero training variance gene

For S2-50, include these at minimum.

---

# Negative controls

Negative controls are mandatory.

## NC1 — Training-label permutation

Randomly permute expression targets among training locations.

Preserve validation/test sets unchanged.

Retrain B1 and B2.

Expected outcome:

- held-out correlations should collapse toward chance
- performance should be substantially worse than the correctly paired model

Run at least 3 permutation seeds for S2-50.

---

## NC2 — Training-feature permutation

Randomly permute image embeddings among training locations while preserving targets.

Validation/test sets remain unchanged.

Retrain B1 and B2.

Expected outcome:

- held-out correlations should collapse substantially

Run at least 3 permutation seeds for S2-50.

---

# Secondary four-patient robustness audit

In addition to the frozen primary split, run a secondary leave-one-patient-out (LOPO) **conditional sensitivity analysis**.

For each of the four patients:

```text
held-out patient = evaluation patient
remaining patients = fold-training pool
```

This analysis is **not an unbiased validation estimate**. The gene panels and primary hyperparameters were selected under the frozen primary C0/C1 development design. Therefore, when an original primary-training patient becomes the LOPO held-out patient, the frozen panel itself may already reflect information from that patient.

Accordingly, LOPO is used only to answer:

> Given the already frozen target panels and primary-selected model configuration, is the observed morphology-expression signal catastrophically dependent on one patient?

Do not describe LOPO as external validation, independent validation, or unbiased generalization performance.

## LOPO preprocessing rule — mandatory

Each LOPO fold must create its own fold-specific learned preprocessing state using **only the remaining fold-training patients**.

For each fold, independently fit and save:

- feature scaler, if used
- target scaler/standardizer, if used
- any learned normalization statistics required by the model

Do **not** reuse learned C0/C1 primary-split scaling statistics inside LOPO when the associated fitting population includes the LOPO held-out patient.

Store these objects under a fold-specific C1 path, for example:

```text
C1/lopo/<heldout_patient_id>/
    feature_scaler.*
    target_scaler.*
    config.json
    metrics.*
```

The frozen expression normalization/data representation from C0 remains unchanged; the rule above applies to any **learned C1 model-input or target scaling** layered on top of those frozen C0 values.

## LOPO model-selection rule

Do not retune architecture or hyperparameters separately per LOPO held-out patient.

Use:

- the B1 alpha selected by the primary train/validation development experiment
- the B2 architecture and training configuration selected by the primary train/validation development experiment

Retrain model parameters from scratch on the fold-training patients using the fold-specific scaling objects.

If early stopping is used in LOPO, do not use the held-out patient for it. Either:

1. use a fully predefined epoch count taken from the primary experiment; or
2. derive an internal stopping rule using fold-training patients only.

Prefer option 1 for this four-patient sensitivity audit to avoid unstable nested tuning.

Report all four held-out patient results separately.

Label the entire analysis:

**Conditional four-patient robustness sensitivity analysis — not unbiased validation**

Do not merge these estimates with the primary frozen test result.

A fully unbiased LOPO study would require fold-specific target-panel selection/ranking and a separately designed nested model-selection protocol. That is outside the scope of C1 and, if needed later, must be created as a separately versioned experiment rather than silently substituted here.

---

# Seeds

Use fixed seeds.

For B1:

- deterministic where possible

For B2 main experiments:

`[0, 1, 2]`

For negative controls:

`[0, 1, 2]`

Record every seed in the experiment registry.

---

# Resource profiling

For B1 and B2 record:

- training wall time
- validation/inference wall time
- peak resident RAM
- model parameter count
- input feature dimension
- output dimension
- number of training locations
- CPU model/thread count if available

C1 must remain CPU-only.

---

# Required outputs

Create:

```text
C1/
├── C1_REPORT.md
├── C1_METRICS_SUMMARY.csv
├── C1_GENE_METRICS.csv
├── C1_PATIENT_METRICS.csv
├── C1_COEXPRESSION_METRICS.csv
├── C1_NEGATIVE_CONTROLS.csv
├── C1_ROBUSTNESS_LOPO.csv
├── C1_RESOURCE_PROFILE.csv
├── C1_SPEC_LOCK.json
├── configs/
├── logs/
├── models/
├── predictions/
└── figures/
```

Also append every run to the global experiment registry:

`results/experiment_registry.csv`

Do not overwrite prior rows.

---

# C1_SPEC_LOCK.json

At completion, freeze at least:

- source commit
- exact C0 hashes used
- train/validation/test patient IDs
- feature cache hash
- target cache hash
- gene-panel lists and order
- single-gene target
- model configurations
- selected Ridge alpha
- selected MLP configuration
- seeds
- preprocessing parameters
- primary-split learned feature/target scalers, if any
- fold-specific LOPO scaler hashes/paths, if any
- library versions
- CPU environment

C1 completion must be reproducible from this lock.

---

# Required summary table

`C1_REPORT.md` must contain one primary table with rows:

```text
B0
B1
B2
```

and columns for:

```text
S1 Pearson
S2-50 median Pearson
S2-250 median Pearson
S2-280 median Pearson
S2-50 coexpression preservation
S2-250 coexpression preservation
S2-280 coexpression preservation
training time
peak RAM
```

For B2, report mean ± standard deviation across seeds.

---

# Required gene-level analysis

For each multi-gene setting, produce:

1. histogram of gene-wise Pearson correlations
2. target variance vs. Pearson scatter
3. top 10 most predictable genes
4. bottom 10 genes
5. number and fraction of genes exceeding Pearson thresholds

The ranking is descriptive only.

Do not make biological claims in C1.

---


## Interpretation rule for the robustness audit

The primary frozen train/validation/test result is the main C1 development estimate.

The four-patient LOPO analysis is secondary and conditional on the previously frozen gene panels and primary-selected hyperparameters. It may support statements about **patient dependence/stability**, but it must not be used as an unbiased cross-validated estimate of model performance.


# C1 acceptance criteria

C1 status may be one of:

- `PASS_DEVELOPMENT`
- `PASS_WEAK_SIGNAL`
- `FAIL_SIGNAL`
- `FAIL_TECHNICAL`

## PASS_DEVELOPMENT

Assign only if all are true:

1. B1 and/or B2 consistently outperform B0.
2. Correct image-expression pairing substantially outperforms both negative controls.
3. A nontrivial subset of genes has positive held-out predictability.
4. Predictive signal is not confined entirely to one patient in the LOPO audit.
5. Multi-gene predictions preserve nonzero gene-gene structure.
6. Spatial maps show meaningful nonconstant variation rather than only regression to the global mean.
7. No leakage or integrity failure is detected.

No fixed Pearson threshold is required.

---

## PASS_WEAK_SIGNAL

Use if:

- correct models beat B0 and controls,
- but the signal is small, unstable across patients, or weak for most genes.

If this occurs, the recommended next step is:

> introduce a frozen pathology-specific image encoder as a separately versioned feature checkpoint before judging the graph hypothesis.

Do not reject the project solely because generic ResNet18 is weak.

---

## FAIL_SIGNAL

Use only if:

- correctly paired models do not meaningfully outperform B0 and negative controls,
- after technical integrity checks pass.

Do not proceed directly to graph modeling.

Audit image representation and spatial alignment first.

---

## FAIL_TECHNICAL

Use for:

- leakage
- broken cache integrity
- target mismatch
- split violation
- failed reproducibility
- NaN/shape/indexing errors
- inconsistent gene ordering

Repair the technical issue before interpreting performance.

---

# Explicitly forbidden in C1

Do not:

- modify C0
- add spatial edges
- train a GNN
- use neighboring expression values
- run S3 partial-ST extrapolation
- use C0 sparse masks as model inputs
- fine-tune an image encoder
- replace ResNet18 embeddings
- perform biological pathway enrichment
- perform explainer analysis
- claim external validation
- claim whole-transcriptome reconstruction
- make publication-level claims from four patients

---

# End-of-checkpoint instruction

After all C1 experiments and audits are complete:

1. write `C1_REPORT.md`
2. write and hash `C1_SPEC_LOCK.json`
3. report the final C1 status
4. list all important limitations
5. list the exact recommendation for C2
6. STOP

Do not execute C2 automatically.

The next decision must be based on C1 results.

---

# Scientific interpretation boundary

C1 answers only:

> Does local H&E morphology, represented by the frozen generic ResNet18 features, contain measurable signal for predicting spatial expression in this development cohort?

C1 does **not** answer:

- whether pathology foundation embeddings are optimal
- whether graph context improves prediction
- whether multiresolution modeling is beneficial
- whether sparse spatial transcriptomics can be extrapolated
- whether explanations are biologically faithful
- whether findings generalize externally
- whether the method is publication-ready

Those questions belong to later checkpoints.
