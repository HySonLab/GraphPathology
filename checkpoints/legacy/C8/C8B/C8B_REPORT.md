# C8B — Fusion calibration

**MIXED_FUSION_CALIBRATION**; technical audit **PASS**. C8C was not run.

All four folds used training-selected T1. F1 with image weight 0.2 was selected in folds 0–2; the five-parameter F3 gate was selected in fold 3. No method family was selected from outer patient outcomes.

At 20% V1, median paired ΔPearson is +0.00856 and median paired relative MSE reduction is 0.995%, but only 3/8 patients meet joint support. MSE harm >2% occurs in 9/24 primary patient-budget settings, and MAE worsens in all 8 at 20%. Only 1/3 V2 geometries supports both metrics. The frozen pass gate fails and its mixed fallback applies.

| regime   |   budget | model    |   median_pearson |   common_mse |   common_mae |
|:---------|---------:|:---------|-----------------:|-------------:|-------------:|
| V1       |  0.10000 | F0       |          0.15255 |      0.07245 |      0.17718 |
| V1       |  0.20000 | F0       |          0.19431 |      0.06797 |      0.17298 |
| V1       |  0.40000 | F0       |          0.24135 |      0.06523 |      0.17052 |
| V2       |  0.20000 | F0       |          0.06379 |      0.07751 |      0.19623 |
| V1       |  0.10000 | F1       |          0.19430 |      0.07467 |      0.19331 |
| V1       |  0.20000 | F1       |          0.23764 |      0.07004 |      0.18709 |
| V1       |  0.40000 | F1       |          0.27947 |      0.06715 |      0.18243 |
| V2       |  0.20000 | F1       |          0.13967 |      0.08086 |      0.21031 |
| V1       |  0.10000 | F2       |          0.19430 |      0.07749 |      0.19828 |
| V1       |  0.20000 | F2       |          0.23764 |      0.07289 |      0.19607 |
| V1       |  0.40000 | F2       |          0.27947 |      0.06950 |      0.19356 |
| V2       |  0.20000 | F2       |          0.13550 |      0.08364 |      0.21551 |
| V1       |  0.10000 | F3       |          0.17268 |      0.07153 |      0.18317 |
| V1       |  0.20000 | F3       |          0.21506 |      0.06710 |      0.17807 |
| V1       |  0.40000 | F3       |          0.26219 |      0.06464 |      0.17455 |
| V2       |  0.20000 | F3       |          0.10855 |      0.07790 |      0.20229 |
| V1       |  0.10000 | SELECTED |          0.17268 |      0.07467 |      0.19331 |
| V1       |  0.20000 | SELECTED |          0.21506 |      0.07004 |      0.18709 |
| V1       |  0.40000 | SELECTED |          0.26219 |      0.06712 |      0.18243 |
| V2       |  0.20000 | SELECTED |          0.10855 |      0.08086 |      0.21031 |

See [full C8 report](../C8_REPORT.md) for every patient, controls, all forty genes, units, fold protocol, locked gates and limitations. Common-unit errors are in the original C7 target units. Native transform errors are saved separately. [Audit](../logs/INTEGRITY_AUDIT.json) passed; parent artifacts and C7 predictions remain unchanged.
