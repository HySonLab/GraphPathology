# C8A — Target harmonization

**PASS_TARGET_HARMONIZATION**. All four nested training-only selections chose T1. HER2ST is development transportability data, not final external validation.

T1 uses separately fitted source and destination mean/SD; frozen source image outputs are mapped through source statistics, and HER2ST targets/anchors through allowed training-patient statistics. This is an assay-level affine calibration, not a refitted encoder/head. Native MSE cannot be compared across transforms; all selection/gates use predictions inverse-mapped to the original C7 units.

At 20% V1, selected fusion vs T0 has median paired Pearson gain +0.01372 and relative common-unit MSE reduction 79.99%; 7/8 patients satisfy joint support. Common-MSE patient IQR falls from 0.06440 to 0.02519. All three V1 budgets have positive median correlation/error gains. V2 median gain is +0.00861 and 8/8 patients pass the error safeguard.

This passes the harmonization gate but does not establish a fusion advantage over ST interpolation; C8B is required. Exact fold/scaler objects were sealed before held-out evaluation in FIT_LOCK.json.

20% coverage, medians across eight patient summaries; native mean_mse and common_mse are different units:

| transform   | model   | regime   |   median_pearson |   mean_mse |   common_mse |   common_mae |
|:------------|:--------|:---------|-----------------:|-----------:|-------------:|-------------:|
| T0          | E0      | V1       |          0.14319 |    4.36650 |      4.36650 |      1.93756 |
| T0          | E0      | V2       |          0.11105 |    4.46020 |      4.46020 |      1.95830 |
| T0          | E1      | V1       |          0.19431 |    0.06797 |      0.06797 |      0.17298 |
| T0          | E1      | V2       |          0.06379 |    0.07751 |      0.07751 |      0.19623 |
| T0          | E2      | V1       |          0.25463 |    0.34809 |      0.34809 |      0.51525 |
| T0          | E2      | V2       |          0.13160 |    1.20468 |      1.20468 |      1.00309 |
| T1          | E0      | V1       |          0.14319 |    1.62996 |      0.23699 |      0.35568 |
| T1          | E0      | V2       |          0.11105 |    1.69841 |      0.24668 |      0.36082 |
| T1          | E1      | V1       |          0.19431 |    0.61722 |      0.06797 |      0.17298 |
| T1          | E1      | V2       |          0.06379 |    0.66610 |      0.07751 |      0.19623 |
| T1          | E2      | V1       |          0.26187 |    0.63949 |      0.07388 |      0.19683 |
| T1          | E2      | V2       |          0.15476 |    0.85746 |      0.12519 |      0.26408 |
| T2          | E0      | V1       |          0.14319 |    0.82677 |      0.24038 |      0.33988 |
| T2          | E0      | V2       |          0.11105 |    0.85458 |      0.25130 |      0.33858 |
| T2          | E1      | V1       |          0.19431 |    0.46510 |      0.06797 |      0.17298 |
| T2          | E1      | V2       |          0.06379 |    0.47571 |      0.07751 |      0.19623 |
| T2          | E2      | V1       |          0.26056 |    0.45176 |      0.07341 |      0.18668 |
| T2          | E2      | V2       |          0.15927 |    0.55824 |      0.11923 |      0.24273 |

T0 reproduces C7 without modification. T2 uses a training-SD fallback for zero IQR; all fallback indicators are in C8A_TARGET_TRANSFORMS.csv. T3/T4 were omitted prospectively.
