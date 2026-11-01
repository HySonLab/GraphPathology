# C4 Execution Plan — Pathology Representation and Cross-Modal Benchmarking

## Repository
All C4 code and artifacts must be created under:

`/home/son/Research/GraphPathology`

C0–C3 are frozen and must not be modified.

## Frozen context

- **C0:** `PASS_DEVELOPMENT_ONLY`, commit `2b54db3379500fec4603a61751145c929db9a72a`
- **C1:** `PASS_DEVELOPMENT`; generic ImageNet ResNet18 features contain measurable local morphology→expression signal.
- **C2:** `MIXED_CONTEXT_SIGNAL`, commit `62805742791b3636a895782bbd3400594eab0295`
- **C3:** `MIXED_SPATIAL_SPECIFICITY`, commit `5ba150951084ea6da003b345e66dec1a966097a5`; C3B was not run because true topology did not reproducibly beat matched controls.

C4 must determine whether the current limitation is mainly caused by the generic image representation.

---

# Scientific objective

C4 asks:

1. **Does a frozen pathology-specific encoder materially improve spatial gene-expression prediction over frozen ResNet18 features?**
2. **If it does, does correct spatial topology become more informative with the stronger representation?**

C4 has two gated stages:

```text
C4A — Pathology representation benchmark
        ↓ only if clear PASS
C4B — Minimal spatial-specificity recheck
```

Do not introduce multiresolution modeling in C4.

---

# Global rules

1. Keep C0 patient splits, target panels, target normalization, and gene ordering frozen.
2. Keep C1–C3 outputs unchanged.
3. In C4A, the primary experimental variable is the image representation.
4. All image encoders must remain frozen; no fine-tuning.
5. Compute embeddings once and cache them.
6. Do not overwrite the C0 ResNet18 cache.
7. Validation data may guide model/hyperparameter selection. Test data must remain untouched until final evaluation.
8. Use the same downstream model families, tuning budget, seeds, and metrics across encoders.
9. Any learned scaling must be fit on training patients only.
10. CPU-first downstream training is preferred. GPU may be used for frozen feature extraction if necessary; record hardware and wall time.
11. No pathway analysis, explainability, clinical claims, S3 partial-ST reconstruction, or publication-level claims in C4.
12. Stop after C4 and wait for review.

---

# C4A — Pathology representation benchmark

## A. Encoders

### E0 — Frozen ResNet18 reference
Use the exact C0 ResNet18 embeddings.

### E1 — Primary pathology foundation encoder
Use one strong frozen pathology encoder.

Preferred order:
1. UNI2-h if technically practical;
2. otherwise UNI;
3. otherwise another strong public pathology foundation encoder.

Freeze before test evaluation:
- exact model/release/version;
- source repository;
- preprocessing;
- patch size;
- magnification/resolution assumptions;
- embedding dimension;
- device/precision/batch size;
- model/cache hashes.

### E2 — Optional second pathology encoder
Use only if practical and scientifically useful, preferably from a different representation family (e.g. CONCH or Virchow2).

Do not run a broad encoder sweep.

Preferred maximum:
`ResNet18 + 1 mandatory pathology encoder + 1 optional pathology encoder`.

---

## B. Patch compatibility audit

Before extraction, verify for each pathology encoder:
- expected patch size;
- normalization;
- color preprocessing;
- magnification assumptions;
- field of view;
- compatibility with frozen C0 patch geometry.

If C0 patches are incompatible, do not alter C0. Create a separately versioned **C4-only patch representation** while preserving the same patient/sample/barcode/location identities. Re-audit spatial alignment and record the new patch hashes.

---

## C. Embedding cache

Create separate caches:

```text
C4/embeddings/
    resnet18_reference/
    pathology_encoder_1/
    pathology_encoder_2/   # optional
```

For every location save:
- patient/sample/barcode ID;
- spatial coordinates;
- embedding;
- source patch identifier/hash;
- encoder/preprocessing hash.

Verify exact barcode order against C0.

---

## D. Prediction settings

Run:
- **S1:** frozen single-gene task
- **S2-50**
- **S2-250**
- **S2-280:** full shared assay panel

Do not call S2-280 whole-transcriptome prediction.

---

## E. Downstream models

Use the same downstream model families for every encoder.

### D0 — Ridge
Use a small predefined alpha grid selected on validation only.

Because embedding distributions differ, a new C4 validation-only grid is allowed, but it must be frozen before final test evaluation.

### D1 — Small MLP
Use the same compact candidate family for all encoders:

```text
hidden width: [128, 256]
dropout: [0.0, 0.1]
```

Use validation data for selection and early stopping.

Do not give one encoder a larger downstream model or larger tuning budget.

---

## F. Primary metrics

For each gene:
- Pearson
- Spearman
- MSE
- MAE
- target variance
- prediction variance

For each panel:
- mean/median Pearson
- mean/median Spearman
- mean/median MSE
- mean/median MAE
- fraction of genes with Pearson > 0, >0.1, >0.2, >0.3

Retain gene-gene correlation preservation only as a secondary structural diagnostic.

---

## G. Encoder gain

For every gene \(g\):

\[
\Delta_g^{encoder}
=
r_g^{pathology}
-
r_g^{ResNet18}
\]

and

\[
\Delta MSE_g^{encoder}
=
MSE_g^{ResNet18}
-
MSE_g^{pathology}.
\]

Report:
- mean/median encoder gain;
- fraction of genes improved in Pearson;
- fraction improved in both Pearson and MSE;
- gene-level gain table.

---

## H. Negative controls

For the strongest pathology encoder, at minimum on S2-50:
- training-label permutation;
- training feature-location permutation;
- at least 3 seeds each.

Validation/test remain unchanged.

Purpose: ensure apparent gains disappear when image-expression pairing is broken.

---

## I. Four-patient conditional sensitivity analysis

Compare the strongest pathology encoder against ResNet18 across all four held-out-patient sensitivity folds.

This remains **conditional sensitivity analysis**, not unbiased cross-validation.

For each fold:
- fit feature and target scaling using fold-training patients only;
- use frozen panels and the primary-selected downstream configuration;
- do not retune architecture using the held-out patient.

Prioritize S2-50 and S2-280. Include S2-250 if inexpensive.

---

## J. Resource/representation diagnostics

For every encoder record:
- embedding dimension;
- extraction wall time;
- extraction device;
- cache size;
- downstream training time;
- peak RAM;
- parameter count of downstream model.

Optional diagnostics:
- PCA/UMAP colored by patient/source study;
- embedding-distance summaries.

Do not interpret latent dimensions biologically.

---

# C4A status

Assign one of:

- `PASS_PATHOLOGY_REPRESENTATION`
- `MIXED_PATHOLOGY_REPRESENTATION`
- `FAIL_PATHOLOGY_REPRESENTATION`
- `FAIL_TECHNICAL`

## PASS
Use only if the pathology encoder:
1. improves at least one major multi-gene panel over ResNet18 on the primary test patient;
2. improves a meaningful subset of genes;
3. does not systematically trade higher correlation for clearly worse MSE;
4. shows benefit on more than one conditional held-out patient;
5. passes integrity and negative controls.

If PASS, C4B may run.

## MIXED
Use if gains are small, patient-dependent, or inconsistent between correlation and error.

STOP for review unless there is a strong technical reason to perform one minimal C4B recheck.

## FAIL
Use if pathology-specific encoders do not materially improve the frozen ResNet18 reference.

STOP. Do not keep adding larger encoders to force a positive result.

---

# C4B — Minimal spatial-specificity recheck

Run only if C4A clearly passes.

## Scientific question

> Does true spatial topology become more informative when node features are pathology-specific?

This is not a repeat of the full C3 stress test.

## Model set

Use only:
1. pathology-encoder local MLP;
2. neighbor mean-pooling + MLP;
3. true spatial GraphSAGE;
4. strongest matched C3 randomized control, preferably the distance-bin-preserving control;
5. optional self-only GraphSAGE.

## Geometry

Use the frozen best C2/C3 single-scale graph geometry.

Do not introduce:
- new radii;
- multiresolution graphs;
- gene-adaptive fusion.

The intended experimental change from C3 is the **node representation**.

## Tasks

Prioritize:
- S2-50
- S2-280

S2-250 is optional if needed.

## Metrics

For each gene:

\[
\Delta_g^{spatial,path}
=
r_g^{true\ spatial,path}
-
r_g^{local,path}
\]

and

\[
\Delta_g^{matched,path}
=
r_g^{true\ spatial,path}
-
r_g^{distance\ matched,path}.
\]

Also report corresponding MSE changes.

Compare directly with frozen C3 ResNet18 results.

---

# C4B status

Assign one of:
- `PASS_SPATIAL_RECHECK`
- `MIXED_SPATIAL_RECHECK`
- `FAIL_SPATIAL_RECHECK`
- `NOT_RUN`
- `FAIL_TECHNICAL`

## PASS
Require:
- true topology beats local morphology and matched random topology for a reproducible subset of genes;
- gains are not confined to one patient;
- MSE behavior is acceptable.

Only a PASS can justify reconsidering multiresolution graph modeling later.

## MIXED
STOP for review.

## FAIL
If matched controls remain comparable to or better than true topology, deprioritize multiresolution graph modeling for this cohort and move next to sparse-ST reconstruction.

---

# Cross-modal benchmark compatibility

Create:

`C4/CROSS_MODAL_COMPATIBILITY.md`

Assess whether a faithful CMRCNet-style comparison is feasible.

Document:
- dataset compatibility;
- gene-panel compatibility;
- required image resolution;
- official preprocessing compatibility;
- required dependencies;
- expected CPU/GPU requirements;
- whether official code can be reproduced without distorting the method.

If exact reproduction is feasible and inexpensive, it may be run as an auxiliary benchmark.

Otherwise defer full CMRCNet reproduction to a later publication-benchmark checkpoint.

Do not implement an approximate method and label it CMRCNet.

---

# Required outputs

```text
C4/
├── C4_REPORT.md
├── C4_STATUS.json
├── CROSS_MODAL_COMPATIBILITY.md
├── C4A/
│   ├── C4A_REPORT.md
│   ├── C4A_SPEC_LOCK.json
│   ├── C4A_ENCODER_MANIFEST.csv
│   ├── C4A_PANEL_METRICS.csv
│   ├── C4A_GENE_METRICS.csv
│   ├── C4A_ENCODER_GAIN.csv
│   ├── C4A_PATIENT_SENSITIVITY.csv
│   ├── C4A_NEGATIVE_CONTROLS.csv
│   ├── C4A_RESOURCE_PROFILE.csv
│   ├── embeddings/
│   ├── configs/
│   ├── models/
│   ├── predictions/
│   ├── logs/
│   └── figures/
└── C4B/                       # only if executed
    ├── C4B_REPORT.md
    ├── C4B_SPEC_LOCK.json
    ├── C4B_PANEL_METRICS.csv
    ├── C4B_GENE_METRICS.csv
    ├── C4B_PATIENT_SENSITIVITY.csv
    ├── configs/
    ├── graphs/
    ├── models/
    ├── predictions/
    ├── logs/
    └── figures/
```

Append all extraction/training/evaluation events to:

`results/experiment_registry.csv`

Do not overwrite existing rows.

---

# Reproducibility requirements

For each encoder:
- exact name/version;
- source;
- preprocessing config;
- model hash if available;
- embedding cache hash;
- patch/barcode order verification;
- extraction hardware/time.

For each downstream model:
- config;
- scaler hashes;
- seed;
- checkpoint;
- predictions;
- exact replay verification.

For C4B:
- graph hash;
- control/randomization seed;
- graph-statistic audit.

Verify that all frozen C0–C3 hashes remain unchanged.

---

# C4 report questions

`C4_REPORT.md` must start with:

```text
C4A status: ...
C4B status: NOT_RUN / ...
Overall C4 status: ...
```

Then answer explicitly:

1. Does a pathology-specific encoder materially improve H&E→expression prediction over ResNet18?
2. Which panels/genes benefit most?
3. Are gains reproducible across conditional held-out patients?
4. If C4B ran, does true tissue topology become more informative?
5. What should happen next?

The final recommendation must be exactly one of:

- `PROCEED_TO_SPARSE_ST_RECONSTRUCTION`
- `RECONSIDER_MULTISCALE_GRAPH`
- `NEED_REPRESENTATION_REVIEW`
- `STOP_TECHNICAL`

`RECONSIDER_MULTISCALE_GRAPH` requires `PASS_SPATIAL_RECHECK`.

---

# Explicitly forbidden in C4

Do not:
- modify C0–C3;
- fine-tune image encoders;
- add multiresolution graphs;
- add gene-adaptive scale fusion;
- run GNNExplainer;
- run pathway enrichment;
- use measured neighboring gene expression as input;
- execute S3 partial-ST reconstruction;
- use test data for encoder/model selection;
- claim external validation;
- claim whole-transcriptome reconstruction;
- run an uncontrolled encoder sweep.

---

# Decision tree

```text
C3 = MIXED_SPATIAL_SPECIFICITY
        ↓
C4A: Better pathology representation?
        │
        ├── FAIL
        │     → STOP C4
        │     → recommend sparse-ST reconstruction
        │
        ├── MIXED
        │     → STOP FOR REVIEW
        │
        └── PASS
              ↓
C4B: Does true topology now beat matched controls?
              │
              ├── FAIL
              │     → deprioritize multiscale graph
              │     → recommend sparse-ST reconstruction
              │
              ├── MIXED
              │     → STOP FOR REVIEW
              │
              └── PASS
                    ↓
          multiscale graph may be reconsidered
          in a later checkpoint
```

---

# End-of-checkpoint instruction

After C4A:
1. write `C4A_REPORT.md`;
2. freeze and hash `C4A_SPEC_LOCK.json`;
3. assign C4A status;
4. if C4A is not a clear PASS, STOP for review.

Only if C4A clearly passes:
5. execute the minimal C4B recheck;
6. write `C4B_REPORT.md`;
7. freeze and hash `C4B_SPEC_LOCK.json`;
8. assign C4B status.

Then:
9. write `C4_REPORT.md`;
10. write `C4_STATUS.json`;
11. write `CROSS_MODAL_COMPATIBILITY.md`;
12. verify frozen C0–C3 hashes;
13. state the exact next-checkpoint recommendation;
14. STOP.

Do not proceed automatically to any later checkpoint.
