# C3 Execution Plan — Spatial Specificity and Scale Identification

## Repository

All code, configs, logs, metrics, figures, and reports for C3 must be created under:

`/home/son/Research/GraphPathology`

C0, C1, and C2 are frozen and must remain unchanged.

---

## Frozen prior checkpoints

### C0

Status: `PASS_DEVELOPMENT_ONLY`

Source commit:

`2b54db3379500fec4603a61751145c929db9a72a`

### C1

Status: `PASS_DEVELOPMENT`

C1 established measurable local morphology-expression signal using frozen generic ImageNet ResNet18 embeddings.

### C2

Status: `MIXED_CONTEXT_SIGNAL`

Frozen GitHub commit:

`62805742791b3636a895782bbd3400594eab0295`

Primary-test median Pearson:

| Gene panel | Frozen local MLP | GraphSAGE |
|---|---:|---:|
| 50 | 0.684 | 0.678 |
| 250 | 0.511 | 0.526 |
| 280 | 0.459 | 0.485 |

C2 also found:

- six genes meeting descriptive contextual-gain criteria;
- context gains that depend on patient;
- five of seven qualifying gene-panel results with worse average LOPO MSE;
- GraphSAGE retaining strong performance when neighbors were shuffled.

Therefore, C2 does **not** yet establish that the observed graph gain is caused by biologically meaningful spatial tissue topology.

C3 must resolve this ambiguity before any gene-adaptive multiresolution architecture is introduced.

---

# C3 scientific objective

C3 asks two sequential questions:

1. **Does correct tissue spatial topology carry information beyond local morphology, model capacity, and generic neighborhood aggregation?**
2. **If true spatial information matters, do different genes benefit from different physical context scales?**

C3 is divided into two gated stages:

```text
C3A — Spatial specificity
        ↓ only if PASS
C3B — Physical-scale identification
```

Do not execute C3B unless C3A passes its predefined gate.

Do not implement gene-adaptive multiresolution fusion in C3.

That architecture belongs to a later checkpoint only if C3B demonstrates genuine scale-dependent behavior.

---

# Global C3 rules

1. Keep all C0 preprocessing, patient splits, gene panels, normalization, and cached ResNet18 embeddings unchanged.
2. Keep the frozen C1 local-MLP results as the primary local baseline.
3. Reuse the best C2 single-scale GraphSAGE configuration as the principal graph architecture unless a comparison requires a tightly matched alternative.
4. Do not introduce a pathology foundation encoder in C3.
5. Do not use sparse spatial-transcriptomics measurements as model inputs.
6. Do not run S3 partial-ST extrapolation.
7. Do not perform pathway enrichment.
8. Do not perform GNNExplainer or biological interpretation.
9. Do not tune on the primary test patient.
10. Validation patient(s) may guide predefined hyperparameter/model selection; test patient(s) remain untouched until final evaluation.
11. All learned feature/target scaling must follow the C1/C2 split discipline.
12. C3 must distinguish correlation improvement from error/calibration improvement.
13. All randomized graph controls must preserve the relevant graph statistic they are intended to control.
14. Every randomization must be reproducible from stored seeds and graph hashes.
15. Stop immediately after C3A if the spatial-specificity gate fails.

---

# C3A — Spatial specificity

## Scientific question

> Does the **correct spatial arrangement of tissue locations** improve spatial gene-expression prediction beyond what can be explained by generic neighborhood aggregation, model capacity, or regularization?

The primary quantity is:

\[
\Delta_g^{spatial}
=
r_g^{true}
-
r_g^{matched\ control}
\]

where \(r_g\) is held-out gene-wise Pearson correlation.

Do not define spatial evidence solely as:

\[
r_g^{graph} - r_g^{MLP}.
\]

A graph can outperform an MLP even when its neighbors are biologically meaningless.

---

# C3A prediction settings

Run:

- S1 — frozen single gene
- S2-50
- S2-250
- S2-280

The 280-gene setting remains:

**Full shared assay panel**

Do not call it whole-transcriptome prediction.

---

# C3A model/control matrix

Use the same frozen ResNet18 node features.

## A0 — Frozen local MLP

Use the exact frozen C1 result.

Do not rerun or retune unless needed only for deterministic artifact loading.

---

## A1 — Self-only GraphSAGE

Use the same GraphSAGE architecture as C2 but remove neighbor information.

Conceptually:

\[
h_i^{next} = f(h_i)
\]

with no aggregated neighbor contribution.

Purpose:

- controls for graph-model parameterization;
- tests whether the architecture itself acts as a stronger MLP.

Parameter count should remain as close as practical to the true GraphSAGE model.

---

## A2 — Neighbor mean pooling + MLP

For each node:

\[
\bar{x}_{N(i)}
=
\frac{1}{|N(i)|}
\sum_{j \in N(i)} x_j
\]

Predict using either:

\[
[x_i \Vert \bar{x}_{N(i)}]
\]

or the exact C2 mean-pooling formulation if already frozen.

Purpose:

- controls for simple local smoothing;
- determines whether learned message passing is needed.

---

## A3 — True spatial GraphSAGE

Use the best frozen C2 GraphSAGE geometry/configuration.

This is the reference spatial model.

No retuning based on C3 test results.

---

## A4 — Fully shuffled-neighbor GraphSAGE

Randomly reassign neighbors while preserving node count and approximate degree distribution.

Requirements:

- same number of nodes;
- same feature matrix;
- same labels;
- approximately same degree distribution;
- graph must remain valid for the implementation;
- no use of target values in randomization.

Use at least 5 graph-randomization seeds.

Purpose:

- tests whether arbitrary message passing reproduces the gain.

---

## A5 — Degree-preserving random graph

Construct a degree-preserving rewired version of the true spatial graph.

Prefer edge swaps or another reproducible rewiring procedure.

Preserve:

- node identities;
- node degrees exactly or as closely as technically possible.

Destroy:

- biological spatial adjacency.

Use at least 5 randomization seeds.

This is stronger than naive neighbor shuffling.

---

## A6 — Distance-bin-preserving shuffled graph

This control is important.

Preserve the empirical distribution of physical edge lengths while breaking the exact tissue topology.

Procedure:

1. compute every true edge length in physical units if available;
2. partition lengths into predefined bins;
3. reassign edges so the resulting graph approximately preserves:
   - degree distribution;
   - edge-length-bin distribution;
   - node count;
4. prevent duplicate/self edges unless the original graph allows them.

Use at least 5 randomization seeds.

Purpose:

> distinguish exact spatial organization from merely connecting nodes at similar physical distances.

---

## A7 — Feature-nearest graph

Optional but recommended if computationally inexpensive.

Construct neighbors using ResNet18 embedding similarity rather than physical proximity.

Purpose:

- tests whether contextual benefit comes from grouping morphologically similar patches instead of spatial adjacency.

Do not mix spatial and feature edges in C3A.

---

# Graph controls audit

For every graph variant report:

- number of nodes;
- number of edges;
- mean degree;
- degree distribution summary;
- connected components;
- mean physical edge length where meaningful;
- edge-length distribution;
- graph hash;
- randomization seed.

Create a graph-control audit table.

For randomized graphs, verify numerically that the intended preserved statistics were actually preserved.

---

# C3A metrics

For every model and gene compute:

- Pearson
- Spearman
- MSE
- MAE
- prediction variance
- target variance

For every panel report:

- median gene-wise Pearson
- mean gene-wise Pearson
- median Spearman
- mean MSE
- median MSE
- fraction Pearson > 0
- fraction Pearson > 0.1
- fraction Pearson > 0.2
- fraction Pearson > 0.3

Also retain the C1/C2 coexpression metric as secondary only.

Do not use coexpression preservation as evidence of spatial specificity.

---

# Gene-specific spatial gain

For every gene \(g\), compute at minimum:

\[
\Delta_g^{MLP}
=
r_g^{true\ spatial}
-
r_g^{local\ MLP}
\]

and

\[
\Delta_g^{shuffle}
=
r_g^{true\ spatial}
-
\mathbb{E}[r_g^{shuffled}]
\]

and

\[
\Delta_g^{degree}
=
r_g^{true\ spatial}
-
\mathbb{E}[r_g^{degree-preserved}]
\]

and

\[
\Delta_g^{distance}
=
r_g^{true\ spatial}
-
\mathbb{E}[r_g^{distance-bin-control}].
\]

These quantities must be stored gene by gene.

The strongest evidence of true spatial specificity is a positive gain against the **distance-matched topology control**, not merely against the MLP.

---

# C3A patient dependence

Repeat the already defined conditional four-patient sensitivity framework.

Important:

- preserve the interpretation as a conditional sensitivity analysis;
- fit fold-specific feature/target scaling using fold-training patients only;
- keep frozen C3A model configurations;
- do not retune architecture separately for each held-out patient;
- do not call it unbiased cross-validation.

For computational practicality, the full randomized-control matrix may first be run on:

- S2-50;
- S2-280;

then expanded to S2-250 only if needed for interpretation.

At minimum, true spatial versus the strongest matched randomized control must be evaluated across all four conditional held-out patients.

---

# C3A randomization uncertainty

For A4, A5, and A6 use at least:

`[0, 1, 2, 3, 4]`

random graph seeds.

Report:

\[
r_g^{control}
=
\text{mean} \pm \text{SD}
\]

across graph randomizations.

For panel-level results, show the full randomization distribution rather than only one random control graph.

---

# C3A descriptive spatial-specificity criteria

Do not use a universal biological significance threshold.

For a gene to be labeled a **candidate spatial-context gene** in C3A, require all of:

1. true-spatial GraphSAGE has positive held-out Pearson;
2. true-spatial performance exceeds frozen local MLP;
3. true-spatial performance exceeds the mean degree-preserving randomized control;
4. true-spatial performance exceeds the mean distance-bin-preserving control;
5. MSE does not materially deteriorate relative to the local baseline;
6. effect direction is not catastrophically inconsistent across the conditional patient-sensitivity analysis.

These are descriptive development criteria, not population-level statistical significance.

Record exact thresholds/tolerances in `C3A_SPEC_LOCK.json` before final test evaluation.

Do not choose them after inspecting C3A test genes.

---

# C3A required figures

Generate:

1. panel-level comparison of all graph controls;
2. true-spatial minus shuffled performance distributions;
3. true-spatial minus degree-preserved performance distributions;
4. true-spatial minus distance-matched performance distributions;
5. per-gene \(\Delta_g^{spatial}\) plot;
6. patient-resolved true-vs-random graph comparison;
7. correlation gain versus MSE change;
8. graph edge-length/degree diagnostic plots.

For candidate spatial-context genes, create spatial expression maps only as descriptive diagnostics.

Do not perform biological interpretation.

---

# C3A gate

C3A status must be one of:

- `PASS_SPATIAL_SPECIFICITY`
- `MIXED_SPATIAL_SPECIFICITY`
- `FAIL_SPATIAL_SPECIFICITY`
- `FAIL_TECHNICAL`

## PASS_SPATIAL_SPECIFICITY

Use only if:

- true spatial topology reproducibly outperforms the strongest matched randomized topology for a meaningful subset of genes;
- the effect is not explained solely by mean pooling or self-only GraphSAGE;
- at least some gains are supported by both correlation and error behavior;
- the signal is not entirely confined to one patient;
- graph-control audits pass.

If PASS, proceed to C3B.

---

## MIXED_SPATIAL_SPECIFICITY

Use if:

- some genes show apparent true-topology benefit;
- but effects are small, patient-dependent, calibration-inconsistent, or not clearly separated from matched random graphs.

If MIXED, STOP.

Do not proceed automatically to C3B.

Review whether stronger pathology image features or a different graph construction should be tested in a separately versioned checkpoint.

---

## FAIL_SPATIAL_SPECIFICITY

Use if:

- true spatial GraphSAGE does not meaningfully outperform matched randomized graphs;
- or shuffled/rewired graphs perform comparably or better across most targets.

If FAIL, STOP.

Do not construct a multiresolution GNN.

Interpret the C2 graph improvement as insufficient evidence that correct tissue topology adds information.

---

# C3B — Physical-scale identification

Execute only if:

`C3A == PASS_SPATIAL_SPECIFICITY`

## Scientific question

> For genes with genuine spatial-context dependence, does predictive performance depend on the physical radius of histological context?

C3B should test the scale hypothesis **before** building any gene-adaptive multiresolution architecture.

---

# Physical-scale design

Use physical distances in microns whenever the metadata support reliable conversion.

Do not define scale only in graph hops.

Create a small predefined radius set covering local to broader context.

The exact radii must be chosen from the assay geometry and frozen before test evaluation.

Example conceptual levels:

```text
R1 — immediate/local context
R2 — small microenvironment
R3 — intermediate microenvironment
R4 — broader tissue context
```

Do not use these example values as numeric radii without deriving the real physical conversion from the dataset.

Save the final physical radii in:

`C3B_SPEC_LOCK.json`

---

# Radius-controlled graphs

For every radius \(r_s\), construct a graph:

\[
G^{(s)} = (V, E^{(s)})
\]

with edges defined by physical distance.

Where graph density becomes excessive, use a frozen maximum-neighbor cap selected on validation only.

Record:

- radius in microns;
- average degree;
- edge count;
- connected components;
- edge-length distribution.

---

# C3B models

Use a fixed compact GraphSAGE architecture across scales.

Do not change hidden size, depth, or regularization per radius unless a small predefined validation-only tuning rule is explicitly frozen.

The main variable should be **context scale**, not model capacity.

Also evaluate simple radius-based neighbor mean pooling at each scale.

This distinguishes:

\[
\text{larger context itself}
\]

from

\[
\text{learned message passing at that context}.
\]

---

# C3B scale-response curve

For every gene:

\[
r_g(s)
=
\text{held-out Pearson at scale } s.
\]

Also record:

\[
MSE_g(s).
\]

Define:

\[
s_g^\star
=
\arg\max_s r_g(s)
\]

only if:

- performance differences are nontrivial;
- the selected scale is reasonably stable across seeds;
- MSE behavior is not contradictory.

Do not call \(s_g^\star\) a biological interaction radius.

Call it:

**preferred predictive morphological context scale**

or

**effective predictive context scale**.

---

# Context-scale effect size

For each gene define:

\[
\Delta_g^{scale}
=
\max_s r_g(s)
-
r_g(R1).
\]

Also compute:

\[
\Delta MSE_g^{scale}
=
MSE_g(R1)
-
\min_s MSE_g(s).
\]

Genes with higher correlation but substantially worse MSE should be flagged as calibration-inconsistent.

---

# C3B scale stability

Evaluate:

- training-seed stability;
- conditional patient sensitivity;
- neighboring-radius stability.

A gene whose optimal radius jumps arbitrarily among all scales across seeds should not be described as scale-specific.

Define a predefined stability metric before final test use, for example:

- fraction of seeds selecting the same or adjacent scale;
- variance of the effective selected radius.

---

# Candidate morphology-scale genes

C3B may define a descriptive set of genes satisfying:

1. C3A spatial-specificity criteria;
2. positive \(\Delta_g^{scale}\);
3. stable preferred scale;
4. acceptable MSE behavior;
5. non-catastrophic patient inconsistency.

Do not perform pathway enrichment yet.

This set will become the hypothesis-generating input for a later multiresolution checkpoint.

---

# Relationship to local predictability

For each gene compare:

\[
P_g^{local}
\]

from C1 with:

\[
\Delta_g^{spatial}
\]

from C3A and:

\[
s_g^\star
\]

from C3B.

Analyze whether genes that are poorly predicted locally benefit more from wider context.

This is descriptive only.

---

# Spatial autocorrelation analysis

Recommended but secondary.

For measured gene expression, calculate a frozen spatial-autocorrelation statistic such as Moran's I where appropriate.

Then examine:

\[
\text{Moran's I}
\quad \text{vs.} \quad
\Delta_g^{spatial}
\]

and:

\[
\text{Moran's I}
\quad \text{vs.} \quad
s_g^\star.
\]

Do not treat correlation between these quantities as mechanistic evidence.

---

# C3B required figures

At minimum:

1. median panel performance versus physical radius;
2. per-gene scale-response heatmap;
3. examples of local-preference genes;
4. examples of broader-context genes;
5. distribution of preferred scales;
6. \(\Delta_g^{scale}\) distribution;
7. local predictability versus contextual gain;
8. spatial autocorrelation versus contextual gain, if computed;
9. patient-resolved scale-response summaries.

Do not annotate pathways or biological mechanisms in C3.

---

# C3B acceptance status

C3B status must be one of:

- `PASS_SCALE_STRUCTURE`
- `MIXED_SCALE_STRUCTURE`
- `FAIL_SCALE_STRUCTURE`
- `FAIL_TECHNICAL`

## PASS_SCALE_STRUCTURE

Use only if:

- multiple genes show reproducible scale-dependent prediction;
- scale effects survive comparison with simple radius mean pooling;
- preferred scales are not random across seeds;
- the phenomenon is not entirely driven by one patient;
- correlation gains are not systematically contradicted by worse error.

A PASS supports moving to a later checkpoint for a gene-adaptive multiresolution architecture.

---

## MIXED_SCALE_STRUCTURE

Use if some scale structure exists but is weak, unstable, or strongly patient-dependent.

STOP for review.

Do not automatically build a gene-adaptive hierarchy.

---

## FAIL_SCALE_STRUCTURE

Use if context radius has little reproducible effect after controlling for graph density/model capacity.

STOP.

Do not construct a multiresolution model merely because it was part of the original plan.

---

# Required C3 files

Create:

```text
C3/
├── C3_REPORT.md
├── C3_STATUS.json
│
├── C3A/
│   ├── C3A_REPORT.md
│   ├── C3A_SPEC_LOCK.json
│   ├── C3A_GRAPH_AUDIT.csv
│   ├── C3A_PANEL_METRICS.csv
│   ├── C3A_GENE_METRICS.csv
│   ├── C3A_RANDOMIZATION_METRICS.csv
│   ├── C3A_PATIENT_SENSITIVITY.csv
│   ├── configs/
│   ├── graphs/
│   ├── predictions/
│   ├── logs/
│   └── figures/
│
└── C3B/                      # create only if C3A passes
    ├── C3B_REPORT.md
    ├── C3B_SPEC_LOCK.json
    ├── C3B_GRAPH_AUDIT.csv
    ├── C3B_PANEL_METRICS.csv
    ├── C3B_GENE_SCALE_METRICS.csv
    ├── C3B_PATIENT_SENSITIVITY.csv
    ├── configs/
    ├── graphs/
    ├── predictions/
    ├── logs/
    └── figures/
```

Append every fit/evaluation/randomization event to:

`results/experiment_registry.csv`

Do not overwrite existing rows.

---

# C3A experiment registry additions

Add fields where missing:

```text
checkpoint
graph_control_type
graph_randomization_seed
graph_hash
num_nodes
num_edges
mean_degree
num_components
mean_edge_length_um
median_edge_length_um
degree_preserved
distance_distribution_preserved
```

---

# C3B experiment registry additions

Add:

```text
physical_radius_um
neighbor_cap
mean_graph_degree
preferred_scale_candidate
```

---

# Reproducibility requirements

For every stored model:

- save config;
- save seed;
- save preprocessing/scaler hashes;
- save graph hash;
- save predictions;
- verify prediction replay from checkpoint;
- record source commit.

For every randomized graph:

- save seed;
- save graph hash;
- save preservation diagnostics.

All C0/C1/C2 frozen hashes must remain unchanged.

---

# C3 report structure

`C3_REPORT.md` must begin with:

```text
C3A status: ...
C3B status: NOT_RUN / ...
Overall C3 status: ...
```

Then answer these questions explicitly.

## Question 1

Does true tissue topology outperform matched randomized topology?

Answer using:

- correlation;
- MSE;
- patient sensitivity;
- graph-control distributions.

## Question 2

Is any graph benefit more than simple neighborhood averaging?

## Question 3

Which genes show the strongest descriptive true-spatial benefit?

Do not attach biological interpretation.

## Question 4

If C3B ran, do genes show reproducible physical-scale preferences?

## Question 5

Is the evidence sufficient to justify a gene-adaptive multiresolution model?

This answer must be one of:

- `YES`
- `NOT YET`
- `NO`

with technical justification.

---

# Explicitly forbidden in C3

Do not:

- modify C0/C1/C2;
- change image encoder;
- use pathology foundation embeddings;
- use measured neighboring gene expression as input;
- run S3 sparse-ST completion;
- introduce cross-modal CMRCNet-style training;
- use gene-expression anchors;
- run GNNExplainer;
- claim causal morphology-expression relations;
- perform pathway enrichment;
- make clinical claims;
- use test data for scale/radius selection;
- introduce a gene-adaptive multiresolution architecture before C3B passes;
- interpret attention weights as explanations.

---

# Decision tree

```text
C2 = MIXED_CONTEXT_SIGNAL
        ↓
C3A: Does correct spatial topology matter?
        │
        ├── FAIL
        │     → STOP
        │     → Do not build multiresolution GNN
        │
        ├── MIXED
        │     → STOP FOR REVIEW
        │
        └── PASS
              ↓
C3B: Do genes have reproducible physical context scales?
              │
              ├── FAIL
              │     → STOP
              │
              ├── MIXED
              │     → STOP FOR REVIEW
              │
              └── PASS
                    ↓
       Recommend later gene-adaptive
       multiresolution checkpoint
```

---

# Scientific interpretation boundary

C3 may establish:

> correct spatial tissue topology contributes predictive information beyond local morphology and matched graph controls;

and, if C3B passes:

> different genes exhibit different effective predictive morphological context scales.

C3 does **not** establish:

- causal tissue interactions;
- biological pathway mechanisms;
- faithful explainability;
- external generalization;
- clinical utility;
- sparse-ST reconstruction performance;
- superiority of pathology foundation embeddings;
- publication-level evidence.

Those questions belong to later checkpoints.

---

# End-of-checkpoint instruction

After C3A:

1. write `C3A_REPORT.md`;
2. freeze and hash `C3A_SPEC_LOCK.json`;
3. assign the C3A status;
4. if status is not `PASS_SPATIAL_SPECIFICITY`, STOP immediately.

Only if C3A passes:

5. freeze C3B physical-scale definitions before final test evaluation;
6. execute C3B;
7. write `C3B_REPORT.md`;
8. freeze and hash `C3B_SPEC_LOCK.json`;
9. assign C3B status;
10. write the combined `C3_REPORT.md`;
11. state whether a later gene-adaptive multiresolution model is justified;
12. STOP.

Do not proceed automatically to any later checkpoint.
