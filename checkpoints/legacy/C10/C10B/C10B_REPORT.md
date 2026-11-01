# C10B — Locked second replication

**PASS_LIMITED_SECOND_REPLICATION**
**Recommendation: READY_WITH_LIMITED_SECOND_REPLICATION**

The fresh GSE210616 cohort contributes 22 TNBC patients and 30,003 retained tissue locations, using the first numbered section for every published patient. All 50 original primary genes are present. No C10 patient was used for training, calibration or selection. Public H&E inputs contain only 24–29 native pixels per 112-µm patch; this limits conclusions about the image baseline and fusion.

## Primary result: 20% distributed coverage

| method   |   median_pearson |   mean_mse |   mean_mae |   patient_support_vs_M0_all_metrics |
|:---------|-----------------:|-----------:|-----------:|------------------------------------:|
| M0       |         0.074917 |  1.96734   |   1.21201  |                                   0 |
| M1       |         0.415862 |  0.0922563 |   0.203274 |                                  22 |
| M3       |         0.406872 |  0.102208  |   0.22876  |                                  22 |

M0 is H&E only, M1 is sparse-ST interpolation, and M3 is fixed fusion. Each value is a median across patient summaries; random seeds are averaged within patients. All 22 patients met the prespecified meaningful ST criterion at all three budgets: Pearson gain at least 0.02, relative MSE reduction at least 2%, and non-worsening MAE. All 22 improved correlation and MSE between 10% and 40% V1 coverage.

Fusion improved 20% V1 correlation in only 3/22 patients and improved MSE or MAE in none. The median paired correlation change was −0.010830 and the median relative MSE increase was 10.7078%. This is negative evidence for the tested fixed fusion rule, not proof that morphology can never help. M3 deliberately reuses the frozen Wu calibration object; it was not refit on C10.

## Coverage, geometry and distance

At 10%, 20% and 40% V1 coverage, ST median patient correlations were 0.380374, 0.415862 and 0.465266; MSE values were 0.100076, 0.092256 and 0.085249. MAE was 0.210553, 0.203274 and 0.196547.

At 20% V2, ST correlation/MSE/MAE were 0.170747/0.149598/0.262427. Every patient worsened in all three metrics relative to V1. The median paired correlation change was −0.266923 and relative MSE increase was 52.8206%. This contrast includes the scoring population and the frozen interpolation family as well as missing-region geometry. C5 block masks are a different development protocol.

Distance replication is partial. Farthest occupied-bin MSE exceeded nearest-bin MSE in 15/22 V1 and 20/22 V2 patients. Strict increases across every occupied bin occurred in only 6/22 and 13/22 patients. Sparse remote bins, tissue composition and changing patient support constrain pooled curves. No universal distance cutoff is claimed. Complete counts and eligible bin correlations are retained in C10B_DISTANCE_TO_ANCHOR.csv and Supplementary Table 7.

Secondary spatial FPS effects are reported in MANUSCRIPT/data/C10_FPS_effects.csv; they do not reopen sampling optimization. Anchor-expression and anchor-position permutations degraded reconstruction on the prespecified control patient, whereas image permutation produced much smaller changes. Matched control values are in C10B_NEGATIVE_CONTROLS.csv. Descriptive per-gene gains are in C10B_GENE_PATIENT_GAINS.csv and C10B_GENE_GAIN_SUMMARY.csv; they are not biological interpretations or outcome-selected gene panels.

## Integrity

The audit reproduced 8,766 method predictions and 35,064 per-gene metric arrays exactly, replayed 968 label-free mask trajectories, and verified 17,532 recorded hidden-label corruption comparisons. Raw inputs, image caches, parent models, scalers, barcodes and the pre-C10 registry prefix matched their hashes. The protocol and mask commits precede target access. Affine M1/M2 equivalence was checked at every inference, and duplicate M2 metric rows were omitted.

C0–C9 were verified read-only through the full chained manifests; their experiments were neither rerun nor retuned. Before/after verification logs are under C10A/logs and C10B/logs. Original registry bytes are preserved. Qualification, downloads, physical-scale checks, resource records and exact run identities are retained with the scientific locks.

## Manuscript handoff

The package contains seven main and ten supplementary figures, each as a 300-dpi PNG and SVG, with source CSVs and figure-specific plotting entrypoints. Four main and eight supplementary tables are available as CSV and Markdown; the large per-gene CSV is compressed. The draft, legends, citation worklist, source index, numerical ledger and claim-evidence matrix retain the negative findings and separate 40-gene HER2ST transportability development from 50-gene fresh validation.

The manuscript is ready for human writing under the stated limitations. It does not establish general H&E benefit, essential graph topology, optimal sampling, whole-transcriptome recovery, cost savings, clinical utility or causal biology. Native-resolution external H&E remains a gap for stronger morphology claims, but does not block the limited sparse-ST manuscript. Stop after C10.
