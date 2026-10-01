# C9B: locked held-out validation



Ten random seeds per patient/coverage, with 10/20/40% measured locations. Tables average seeds within geometry, then geometries within patient; headline values are patient medians, not pooled spots. Correlations exclude undefined constant-gene values; threshold fractions count such genes as not passing. Errors use identical C0 standardized log1p-count units across methods. Pearson gains below are paired by patient; subtracting cohort medians is not substituted for the paired effect.

At 20% V1, all three patients beat frozen H&E in correlation, MSE and MAE; each does so at all three budgets. Median paired Pearson gain is **+0.4042** and median paired relative MSE reduction is **96.10%**. These very large error reductions are partly relative to a badly transported H&E baseline; they do not establish accurate absolute molecular abundance across assays. The image-resolution and count-distribution shifts limit how strongly M2-M0 can be interpreted as an algorithmic advantage.

| method   |   median_pearson |   mean_mse |   mean_mae |
|:---------|-----------------:|-----------:|-----------:|
| M0       |           0.0590 |     2.9155 |     1.5578 |
| M1       |           0.3590 |     0.0987 |     0.2179 |
| M2       |           0.3590 |     0.0987 |     0.2179 |
| M3       |           0.3349 |     0.1002 |     0.2321 |

All held-out patients at 20% V1:

| patient_id   | method   |   median_pearson |   mean_spearman |   mean_mse |   mean_mae |
|:-------------|:---------|-----------------:|----------------:|-----------:|-----------:|
| 1160920F     | M0       |          -0.0452 |         -0.0429 |     0.8460 |     0.7358 |
| 1160920F     | M1       |           0.3590 |          0.3286 |     0.0817 |     0.1951 |
| 1160920F     | M2       |           0.3590 |          0.3381 |     0.0817 |     0.1951 |
| 1160920F     | M3       |           0.3349 |          0.3164 |     0.0949 |     0.2129 |
| CID4290      | M0       |           0.1614 |          0.1467 |     2.9155 |     1.5578 |
| CID4290      | M1       |           0.3253 |          0.3002 |     0.0987 |     0.2179 |
| CID4290      | M2       |           0.3253 |          0.3021 |     0.0987 |     0.2179 |
| CID4290      | M3       |           0.3329 |          0.3104 |     0.1002 |     0.2321 |
| CID44971     | M0       |           0.0590 |          0.0605 |     2.9491 |     1.5710 |
| CID44971     | M1       |           0.5017 |          0.4097 |     0.1150 |     0.2437 |
| CID44971     | M2       |           0.5017 |          0.4140 |     0.1150 |     0.2437 |
| CID44971     | M3       |           0.4957 |          0.4052 |     0.1180 |     0.2553 |

## What harmonization can and cannot establish

The prospective lock explicitly records that T1 is affine per gene. For weights summing to one, inverse_transform(sum(w*transform(y))) = sum(w*y). Consequently M1/M2 are equal in all 414 cases, with maximum absolute difference **2.66e-15**, entirely numerical roundoff. Their Pearson and common-unit error effects are zero. **Raw Spearman values can differ** because machine-precision perturbations break tied ranks differently; the table above preserves the locked raw metric. These rank differences are numerical artifacts, not harmonization gains. `C9B_NUMERICAL_RANK_SENSITIVITY.csv` quantifies them without replacing predictions or metrics. This is not replication of C8's target-calibrated image/fusion benefit. C9 has no matched unharmonized-fusion arm, so M3-M2 cannot isolate that calibration effect either.

Full PASS required M2 to improve upon or be more stable than M1. That clause is scientifically unidentifiable with these specified methods, and was not waived after observing outcomes. Overall **MIXED** therefore coexists with positive core sparse-ST evidence. Native z-score errors are not reported as if they were gains in common-unit errors.

## H&E fusion

M3 uses exactly 0.2 calibrated image plus 0.8 harmonized ST, with no weight fitting. At 20% V1 only 1/3 patients improves Pearson, while 0/3 improves MSE or MAE. The median paired Pearson effect is **-0.0060**, median relative MSE increase **2.62%**. Image fusion fails the predeclared secondary criterion. The absence of benefit applies to these frozen inputs and this cohort; it does not establish that native high-resolution histology never helps.

| patient_id   |   delta_pearson |   relative_mse_reduction |   mae_reduction |
|:-------------|----------------:|-------------------------:|----------------:|
| 1160920F     |         -0.0242 |                  -0.1616 |         -0.0178 |
| CID4290      |          0.0076 |                  -0.0150 |         -0.0143 |
| CID44971     |         -0.0060 |                  -0.0262 |         -0.0115 |

## Coverage, geometry and distance

|   budget | method   |   median_pearson |   mean_mse |   mean_mae |
|---------:|:---------|-----------------:|-----------:|-----------:|
|   0.1000 | M0       |           0.0616 |     2.9196 |     1.5592 |
|   0.1000 | M1       |           0.3127 |     0.1080 |     0.2276 |
|   0.1000 | M2       |           0.3127 |     0.1080 |     0.2276 |
|   0.1000 | M3       |           0.2979 |     0.1073 |     0.2395 |
|   0.2000 | M0       |           0.0590 |     2.9155 |     1.5578 |
|   0.2000 | M1       |           0.3590 |     0.0987 |     0.2179 |
|   0.2000 | M2       |           0.3590 |     0.0987 |     0.2179 |
|   0.2000 | M3       |           0.3349 |     0.1002 |     0.2321 |
|   0.4000 | M0       |           0.0582 |     2.9057 |     1.5543 |
|   0.4000 | M1       |           0.3934 |     0.0899 |     0.2084 |
|   0.4000 | M2       |           0.3934 |     0.0899 |     0.2084 |
|   0.4000 | M3       |           0.3881 |     0.0934 |     0.2249 |

M2 correlation increases and MSE decreases from 10% to 40% in all 3 patients. Spatial FPS is a separately reported, label-free secondary policy, not an optimized alternative. The regular-policy 20% comparison is below; it was not selected using outcomes. All regular-policy budgets remain in `C9B/C9B_PATIENT_METRICS.csv`.

| patient_id   | regime   | method   |   median_pearson |   mean_mse |   mean_mae |
|:-------------|:---------|:---------|-----------------:|-----------:|-----------:|
| 1160920F     | V1       | M0       |          -0.0421 |     0.8452 |     0.7351 |
| 1160920F     | V1       | M1       |           0.3882 |     0.0808 |     0.1938 |
| 1160920F     | V1       | M2       |           0.3882 |     0.0808 |     0.1938 |
| 1160920F     | V1       | M3       |           0.3451 |     0.0946 |     0.2129 |
| 1160920F     | V2       | M0       |          -0.0469 |     0.8307 |     0.7231 |
| 1160920F     | V2       | M1       |           0.1804 |     0.1022 |     0.2196 |
| 1160920F     | V2       | M2       |           0.1804 |     0.1022 |     0.2196 |
| 1160920F     | V2       | M3       |           0.1526 |     0.1188 |     0.2360 |
| CID4290      | V1       | M0       |           0.1656 |     2.9335 |     1.5647 |
| CID4290      | V1       | M1       |           0.3556 |     0.0947 |     0.2129 |
| CID4290      | V1       | M2       |           0.3556 |     0.0947 |     0.2129 |
| CID4290      | V1       | M3       |           0.3477 |     0.0968 |     0.2286 |
| CID4290      | V2       | M0       |           0.1729 |     2.9558 |     1.5720 |
| CID4290      | V2       | M1       |           0.0311 |     0.1313 |     0.2547 |
| CID4290      | V2       | M2       |           0.0311 |     0.1313 |     0.2547 |
| CID4290      | V2       | M3       |           0.0846 |     0.1262 |     0.2610 |
| CID44971     | V1       | M0       |           0.0587 |     2.9930 |     1.5869 |
| CID44971     | V1       | M1       |           0.5212 |     0.1109 |     0.2386 |
| CID44971     | V1       | M2       |           0.5212 |     0.1109 |     0.2386 |
| CID44971     | V1       | M3       |           0.5164 |     0.1154 |     0.2521 |
| CID44971     | V2       | M0       |           0.0162 |     3.0503 |     1.6043 |
| CID44971     | V2       | M1       |           0.1020 |     0.1746 |     0.2994 |
| CID44971     | V2       | M2       |           0.1020 |     0.1746 |     0.2994 |
| CID44971     | V2       | M3       |           0.0926 |     0.1647 |     0.2996 |

At 20% V2, M2 median patient Pearson falls to **0.0815**, with MSE **0.1366** and MAE **0.2593**. Every patient worsens versus V1. At CID4290, V2 M2 Pearson is below M0 despite much lower error; this is disclosed rather than hidden. V2 pools no spots across patients, and all three predeclared geometry results are available in `C9B_GEOMETRY_METRICS.csv`. V2 differs both in missing-region shape and the frozen interpolation family (mean16 vs inverse-distance16 in V1), so this is protocol-level degradation, not an isolated causal estimate of mask geometry.

| patient_id   | method   |   median_pearson |   mean_mse |   mean_mae |
|:-------------|:---------|-----------------:|-----------:|-----------:|
| 1160920F     | M0       |          -0.0469 |     0.8307 |     0.7231 |
| 1160920F     | M1       |           0.1764 |     0.1028 |     0.2207 |
| 1160920F     | M2       |           0.1764 |     0.1028 |     0.2207 |
| 1160920F     | M3       |           0.1569 |     0.1170 |     0.2348 |
| CID4290      | M0       |           0.1729 |     2.9558 |     1.5720 |
| CID4290      | M1       |           0.0360 |     0.1366 |     0.2593 |
| CID4290      | M2       |           0.0360 |     0.1366 |     0.2593 |
| CID4290      | M3       |           0.0855 |     0.1290 |     0.2635 |
| CID44971     | M0       |           0.0162 |     3.0503 |     1.6043 |
| CID44971     | M1       |           0.0815 |     0.1767 |     0.2988 |
| CID44971     | M2       |           0.0815 |     0.1767 |     0.2988 |
| CID44971     | M3       |           0.0845 |     0.1675 |     0.3004 |

Paired V2-minus-V1 M2 degradation:

| patient_id   |   delta_pearson_V2_minus_V1 |   mse_increase_V2_minus_V1 |   relative_mse_increase |   mae_increase |
|:-------------|----------------------------:|---------------------------:|------------------------:|---------------:|
| 1160920F     |                     -0.1827 |                     0.0211 |                  0.2578 |         0.0256 |
| CID4290      |                     -0.2893 |                     0.0379 |                  0.3841 |         0.0414 |
| CID44971     |                     -0.4202 |                     0.0617 |                  0.5361 |         0.0551 |

M2 MSE rises across occupied distance bins in each patient's 20% V1 curve. V2 broadly worsens with distance but is nonmonotonic, with gene/tissue-composition and geometric differences. Empty bins are omitted; binwise means average available runs/geometries and do not imply independent biological replicates. No unique extrapolation horizon or spot-level p-value is claimed.

## Controls and integrity

Expression-anchor and position-assignment permutations reduce M2 correlation and worsen error on the predeclared 1160920F / 20% / seed 0 case in both regimes, across 3 permutation seeds. NC1/NC2 are related assignment disruptions, not independent biological evidence. Image permutation slightly improves M3 here, consistent with weak/misaligned image contribution rather than supportive morphology dependence. M0 remains unchanged under molecular controls and M1/M2 unchanged under image permutation.

| regime   | control   | method   |   delta_pearson |   mse_reduction |   mae_reduction |
|:---------|:----------|:---------|----------------:|----------------:|----------------:|
| V1       | NC1       | M2       |         -0.3640 |         -0.0479 |         -0.0556 |
| V1       | NC1       | M3       |         -0.3479 |         -0.0424 |         -0.0429 |
| V1       | NC2       | M2       |         -0.3601 |         -0.0477 |         -0.0551 |
| V1       | NC2       | M3       |         -0.3464 |         -0.0423 |         -0.0429 |
| V1       | NC3       | M2       |          0.0000 |          0.0000 |          0.0000 |
| V1       | NC3       | M3       |          0.0062 |          0.0002 |          0.0005 |
| V2       | NC1       | M2       |         -0.0878 |         -0.0217 |         -0.0343 |
| V2       | NC1       | M3       |         -0.0833 |         -0.0326 |         -0.0334 |
| V2       | NC2       | M2       |         -0.1023 |         -0.0287 |         -0.0390 |
| V2       | NC2       | M3       |         -0.0952 |         -0.0408 |         -0.0387 |
| V2       | NC3       | M2       |          0.0000 |          0.0000 |          0.0000 |
| V2       | NC3       | M3       |          0.0026 |          0.0004 |          0.0002 |

Technical audit **PASS**: 132 mask trajectories replayed, 1,656 method predictions replayed exactly, 6,624 per-gene metric-array replays, 3,312 hidden-label NaN/random invariance checks, exact calibration-scaler rebuild using allowed patients only, explicit forbidden evaluation-patient calibration access rejected, and source/model/image hash checks. All 414 run records contain hashes, barcodes, masks, calibration IDs, runtime and peak RAM. All C0-C8 hashes and prior registry bytes are unchanged. The registry contains 132 mask events, 1 calibration event and 414 evaluation events.

## Gene-specific evidence and transportability

All per-gene Pearson/Spearman/MSE/MAE values, patient-level gains, and all gene support counts are saved. At 20% V1, **50/50 genes** improve Pearson, MSE, MAE over M0 in at least 2/3 patients; **0/50** do so for M3 over M2. These are descriptive results, not a newly selected panel or multiple-testing-adjusted biological discovery.

`C9_TRANSPORTABILITY_SUMMARY.csv` uses the same C7/C8 common 40 genes and common error units for development, HER2ST and fresh C9. The C8A T0 rows replay C7 exactly; C8A T1 and C8B rows are transportability development, not fresh independent validation. The main C9 results use 50 genes and must not be directly compared to 40-gene median results. Absolute errors remain sensitive to platform capture/depth despite common mathematical units; no clinical utility, cost saving or whole-transcriptome claim is made.

## Recommendation and stopping point

**PROCEED_TO_SECOND_REPLICATION_COHORT**. A second fresh cohort with better-resolution paired H&E would strengthen the core sparse-ST claim and clarify its limits. The frozen four-method C9 comparison cannot prove an affine-harmonization benefit over identical raw interpolation. Review that estimand before designing another harmonization test; do not reinterpret C9 or retune on its held-out patients. No second cohort was readily qualified after C9B, so optional C9C was not run. Stop for review; do not generate the manuscript automatically.

