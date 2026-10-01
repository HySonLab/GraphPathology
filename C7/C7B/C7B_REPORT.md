# C7B — Locked external evaluation

**MIXED_EXTERNAL_VALIDATION**, technical audit **PASS**. Recommendation: **NEED_METHOD_REVIEW**.

The complete methods, all ten checkpoint answers, gate, patient effects, negative controls, domain-shift limitations and figure captions are in [C7_REPORT.md](../C7_REPORT.md). This stage completed 8,448 main and 108 control evaluations; no tuning followed external-label access.

At 20% random primary coverage, fusion improves over H&E in all eight patients (median paired Δr +0.0603) but improves correlation over interpolation in only three and MSE in none. The separate cohort medians must not be substituted for paired effects.

| regime   |   nominal_budget | model   |   median_pearson |   mean_mse |   mean_mae |
|:---------|-----------------:|:--------|-----------------:|-----------:|-----------:|
| V1       |           0.1000 | E0      |           0.1393 |     4.3683 |     1.9379 |
| V1       |           0.1000 | E1      |           0.1525 |     0.0724 |     0.1772 |
| V1       |           0.1000 | E2      |           0.2426 |     0.3582 |     0.5216 |
| V1       |           0.1000 | E3      |           0.2356 |     0.5514 |     0.6437 |
| V1       |           0.2000 | E0      |           0.1432 |     4.3665 |     1.9376 |
| V1       |           0.2000 | E1      |           0.1943 |     0.0680 |     0.1730 |
| V1       |           0.2000 | E2      |           0.2546 |     0.3481 |     0.5153 |
| V1       |           0.2000 | E3      |           0.2507 |     0.5643 |     0.6537 |
| V1       |           0.4000 | E0      |           0.1671 |     4.3621 |     1.9393 |
| V1       |           0.4000 | E1      |           0.2414 |     0.0652 |     0.1705 |
| V1       |           0.4000 | E2      |           0.2718 |     0.3425 |     0.5122 |
| V1       |           0.4000 | E3      |           0.2605 |     0.5806 |     0.6651 |
| V2       |           0.1000 | E0      |           0.1111 |     4.4602 |     1.9583 |
| V2       |           0.1000 | E1      |           0.0482 |     0.0767 |     0.1954 |
| V2       |           0.1000 | E2      |           0.1245 |     1.2102 |     1.0086 |
| V2       |           0.1000 | E3      |           0.1288 |     1.3722 |     1.0571 |
| V2       |           0.2000 | E0      |           0.1111 |     4.4602 |     1.9583 |
| V2       |           0.2000 | E1      |           0.0638 |     0.0775 |     0.1962 |
| V2       |           0.2000 | E2      |           0.1316 |     1.2047 |     1.0031 |
| V2       |           0.2000 | E3      |           0.1316 |     1.3704 |     1.0586 |
| V2       |           0.4000 | E0      |           0.1111 |     4.4602 |     1.9583 |
| V2       |           0.4000 | E1      |           0.0915 |     0.0773 |     0.1961 |
| V2       |           0.4000 | E2      |           0.1360 |     1.1982 |     0.9993 |
| V2       |           0.4000 | E3      |           0.1379 |     1.3743 |     1.0626 |

[C7B_PANEL_METRICS.csv](C7B_PANEL_METRICS.csv) contains all per-run metrics and Pearson fractions. The compressed gene table reproduces the local plain CSV. [PATIENT_EFFECTS](C7B_PATIENT_EFFECTS.csv), [COVERAGE_CURVES](C7B_COVERAGE_CURVES.csv), [DISTANCE_TO_ANCHOR](C7B_DISTANCE_TO_ANCHOR.csv) and [NEGATIVE_CONTROLS](C7B_NEGATIVE_CONTROLS.csv) retain every declared setting.

Independent replay and frozen-parent verification passed: [INTEGRITY_AUDIT.json](logs/INTEGRITY_AUDIT.json).
