C6 overall status: **MIXED_SAMPLING_OPTIMIZATION**

Final recommendation: **NEED_SAMPLING_REVIEW**

# C6 — Sampling optimization for sparse spatial transcriptomics

C6A returned **MIXED_SAMPLING_EFFECT** under the frozen rules. **C6B was not run** because C6A did not pass. C5 remains frozen at PASS_SPARSE_ST_RECONSTRUCTION. **Spatially regular placement improves open-tissue interpolation, but the benefit does not extend reliably to the fixed unavailable region.** Morphology-only diversity does not establish an open-sampling advantage, and joint diversity is patient/geometry dependent. The strongest sampling effects occur with simple reconstruction; the ResNet18 propagation model is less sensitive to placement. No conclusion here invalidates C5 sparse-ST reconstruction.

## 1. Does measurement placement matter at fixed ST coverage?

The table summarizes primary-test effects over both panels, encoders and the four matched budgets, separately by geometry. Each row has 16 combinations; positive values mean improvement over the ten-seed random mean. Exact budget-, patient- and model-specific results appear in the [C6A report](C6A/C6A_REPORT.md).

| Geometry | Policy | Qualifying combinations / 16 | Mean Pearson gain | Mean relative MSE reduction |
| --- | --- | --- | --- | --- |
| contiguous | A1 | 0 | -0.020 | -2.6% |
| contiguous | A2 | 1 | +0.000 | -0.5% |
| contiguous | A3 | 2 | +0.004 | -0.1% |
| open | A1 | 8 | +0.016 | +4.6% |
| open | A2 | 0 | -0.036 | -10.7% |
| open | A3 | 2 | +0.001 | -0.7% |

For the simple interpolation baseline at 20% open coverage, A1 improves median Pearson by about 0.028/0.027 and reduces MSE by 11.7%/7.1% on S2-50/S2-280, respectively. These pointwise interpolation gains do not satisfy the separate contiguous-region requirement.

The frozen PASS requires portable benefit in both open and contiguous-unavailable sampling, acceptable MSE behavior, and complete integrity audits. Its conjunction was not satisfied. Open scoring membership varies by policy; contiguous scoring uses the same withheld region for every policy and budget.

## 2. Which strategy performs best without hidden-label access?

**A1 spatial FPS** is the strongest nonrandom strategy selected by validation for both encoders. That selection is fixed; held-out winners are not substituted afterward. Random sampling remains a ten-seed distribution, with deterministic A1/A2/A3 trajectories compared at exactly the same number of measured locations. C5 validation chooses ResNet18 B4 and Phikon B3 reconstruction, alongside the simple C5 molecular interpolation baseline. Encoder-branch comparisons therefore also involve different reconstruction families. Phikon S2-280 open B3 has lambda=0 and equals molecular interpolation, so its sampling gains do not require H&E prediction input.

## 3. Does H&E-guided morphology diversity improve over random/grid sampling?

Morphology-only A2 and joint A3 are evaluated against both random A0 and spatial A1. The [morphology comparisons and full gate table](C6A/C6A_REPORT.md) report positive and negative outcomes. The joint-policy grid selected alpha=0.75 for both encoders, yet A1 had the best validation MSE. A2 has no qualifying open-sampling setting for either panel/encoder in the primary or conditional analyses. Its image-permutation control improves, rather than removes, its average open-sampling performance. Joint A3 shows isolated gains but fails the full cross-geometry gate. These findings concern the fixed eight-component diversity heuristic; they do not negate the H&E fusion signal found in C5.

## 4. Does optimized placement reduce nearest-anchor distances and contiguous-extrapolation error?

Distance reduction and reconstruction improvement are evaluated separately. At primary-patient 20% coverage, the ResNet18 sampling trajectories have these 90th-percentile anchor distances:

| Scoring regime | Random mean distance (µm) | Spatial FPS distance (µm) |
| --- | --- | --- |
| open | 200.0 | 141.4 |
| contiguous | 4063.1 | 4001.2 |

Spatial FPS targets geometric coverage. Contiguous targets remain permanently unavailable, so distributing anchors elsewhere cannot remove the unmeasured region. Distance-bin errors and [matched-bin contrasts](C6A/C6A_MATCHED_DISTANCE.csv) quantify remaining error; the fixed gate did not establish the required portable combined benefit. No scalar extrapolation horizon is inferred.

## 5. Is any gain preserved across multiple patients?

The [patient-resolved effects](C6A/C6A_EFFECTS.csv) and C6A report distinguish qualifying settings from the stronger joint gate. The following counts show patients with at least one qualifying setting; the full gate also requires primary support and acceptable aggregate error behavior.

| Panel | Encoder | Policy | Open patients / 4 | Contiguous patients / 4 |
| --- | --- | --- | --- | --- |
| S2-280 | Phikon | A1 | 4 | 0 |
| S2-280 | Phikon | A2 | 0 | 0 |
| S2-280 | Phikon | A3 | 3 | 1 |
| S2-280 | ResNet18 | A1 | 1 | 0 |
| S2-280 | ResNet18 | A2 | 0 | 0 |
| S2-280 | ResNet18 | A3 | 0 | 0 |
| S2-50 | Phikon | A1 | 4 | 0 |
| S2-50 | Phikon | A2 | 0 | 0 |
| S2-50 | Phikon | A3 | 4 | 0 |
| S2-50 | ResNet18 | A1 | 2 | 0 |
| S2-50 | ResNet18 | A2 | 0 | 0 |
| S2-50 | ResNet18 | A3 | 2 | 0 |

The final status reflects all four conditional patients, both encoders, both panels and both geometries. These are overlapping conditional development folds, not unbiased external validation. One unavailable region per patient limits extrapolation-geometry generalization.

## 6. Does adaptive acquisition beat the best one-shot strategy?

**Not assessed.** C6B's required C6A PASS gate was not met. No adaptive trajectory, uncertainty model or model-update search was run.

## 7. At what budgets are strategy differences largest?

The largest average absolute Pearson differences on the primary patient occur at **contiguous: 5% (mean absolute Pearson difference 0.021); open: 5% (mean absolute Pearson difference 0.023)**. This averages the paired policy-minus-random differences over both encoders, panels and nonrandom policies within each geometry; magnitude includes harmful differences and is not an optimal-budget claim.

All 5/10/20/40% curves, random-seed distributions and budget-specific effect ranges are in the [C6A report](C6A/C6A_REPORT.md). [AUC-coverage](C6A/C6A_AUC_COVERAGE.csv) summarizes the same four points descriptively; it does not replace the raw curves or justify a minimum required budget. No test-derived fidelity cutoff or cost-saving estimate is used.

## 8. Is the result strong enough to justify external validation?

**NEED_SAMPLING_REVIEW**

Sampling effects require review before committing to an optimized policy for external validation. C5 remains a viable sparse-reconstruction result.

Technical audit: **PASS**. All 296 test sampling trajectories and 464 stratified prediction/metric runs replayed exactly; 592 selection and 8,896 reconstruction corruption comparisons passed. C0–C5 and the prior registry prefix are unchanged. Code, configs, traces, metrics, figures, reports and resource provenance are saved under C6. No external validation was executed. **Stop after C6 for review.**
