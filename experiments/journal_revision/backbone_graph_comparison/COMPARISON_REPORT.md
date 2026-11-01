# Matched backbone comparison report

Status: COMPLETE. ResNet18 selected heads/results: REUSED_VERIFIED; Phikon selected heads/results: COMPLETE. All 72 source-validation candidates are preserved; 36 validation-selected predictors are evaluated. This exploratory follow-up uses previously inspected patients, not fresh blind validation.

Primary condition: historical distributed 20%, exact masks 0/1/2; three training seeds averaged within patient. External predictions use the fixed C9 T1 contract once. Pearson, MSE and MAE are interpreted jointly.

## Question A — graph context versus a matched local image head

Primary external interpretation: image-only gains depend on the encoder and endpoint. ResNet18 GCN improves the three Bassiouni median endpoints, while GraphSAGE trades lower errors for lower correlation. Phikon GCN has small error reductions without a clear patient-level correlation gain; Phikon GraphSAGE worsens median errors. Wu Phikon image-only correlation gains also accompany worse errors.

With identical anchors, GCN improves all three endpoints over matched local reconstruction in all 22 Bassiouni patients and all three Wu patients in both backbones. Bassiouni relative MSE reductions are 17.84% for ResNet18 and 1.70% for Phikon. These are gains versus local-assisted prediction, not versus IDW: both GCN variants worsen all three Bassiouni cohort-median endpoints against IDW. Phikon GCN in Wu slightly improves median correlation/MSE against IDW but worsens MAE in all three patients. GraphSAGE and the controls remain reported below. The separate common-valid sensitivity leaves both primary GCN correlation gains unchanged.

### ResNet18

| Cohort | Comparison | Median Δr | Median relative MSE (%) | Median ΔMAE | Improved patients (r/MSE/MAE) |
|---|---|---:|---:|---:|---|
| Development | GCN − local | -0.0280 | -1.26 | -0.0012 | 0/1; 1/1; 1/1 |
| Development | GraphSAGE − local | -0.0041 | -2.80 | -0.0069 | 0/1; 1/1; 1/1 |
| Wu | GCN − local | +0.0621 | -10.81 | -0.0087 | 3/3; 2/3; 2/3 |
| Wu | GraphSAGE − local | +0.0225 | -4.88 | -0.0098 | 2/3; 3/3; 3/3 |
| Bassiouni | GCN − local | +0.0079 [-0.0088, +0.0196] | -12.97 [-16.64, -8.90] | -0.0150 [-0.0274, -0.0039] | 14/22; 22/22; 16/22 |
| Bassiouni | GraphSAGE − local | -0.0089 [-0.0228, +0.0035] | -4.92 [-7.82, -2.22] | -0.0059 [-0.0148, -0.0021] | 8/22; 19/22; 17/22 |

Opposite Pearson/error directions represent a tradeoff. Mean and shuffled controls are reported separately below; these primary contrasts alone do not establish topology specificity.

### Phikon

| Cohort | Comparison | Median Δr | Median relative MSE (%) | Median ΔMAE | Improved patients (r/MSE/MAE) |
|---|---|---:|---:|---:|---|
| Development | GCN − local | -0.0165 | +5.02 | +0.0256 | 0/1; 0/1; 0/1 |
| Development | GraphSAGE − local | -0.0068 | -1.48 | -0.0107 | 0/1; 1/1; 1/1 |
| Wu | GCN − local | +0.0020 | +0.40 | +0.0016 | 2/3; 1/3; 0/3 |
| Wu | GraphSAGE − local | +0.0036 | +9.39 | +0.0092 | 2/3; 0/3; 0/3 |
| Bassiouni | GCN − local | +0.0007 [-0.0041, +0.0071] | -0.43 [-0.61, -0.16] | -0.0002 [-0.0007, +0.0009] | 11/22; 16/22; 12/22 |
| Bassiouni | GraphSAGE − local | -0.0016 [-0.0118, +0.0051] | +8.21 [+5.69, +8.91] | +0.0023 [-0.0021, +0.0139] | 9/22; 1/22; 9/22 |

Opposite Pearson/error directions represent a tradeoff. Mean and shuffled controls are reported separately below; these primary contrasts alone do not establish topology specificity.

## Question B — observed-residual reconstruction, then molecular support

All contrasts below use the same anchors and IDW weights within a run. Gains versus local and versus IDW answer separate questions.

### ResNet18

| Cohort | Candidate / reference | Median Δr | Median relative MSE (%) | Median ΔMAE | Improved patients (r/MSE/MAE) |
|---|---|---:|---:|---:|---|
| Development | GCN / local | -0.0108 | -2.21 | -0.0075 | 0/1; 1/1; 1/1 |
| Development | GraphSAGE / local | +0.0045 | -8.90 | -0.0211 | 1/1; 1/1; 1/1 |
| Development | local / IDW | +0.0286 | -4.37 | -0.0013 | 1/1; 1/1; 1/1 |
| Development | GCN / IDW | +0.0178 | -6.48 | -0.0089 | 1/1; 1/1; 1/1 |
| Development | GraphSAGE / IDW | +0.0331 | -12.88 | -0.0224 | 1/1; 1/1; 1/1 |
| Wu | GCN / local | +0.0455 | -14.29 | -0.0234 | 3/3; 3/3; 3/3 |
| Wu | GraphSAGE / local | +0.0357 | -11.05 | -0.0180 | 3/3; 3/3; 3/3 |
| Wu | local / IDW | -0.0732 | +17.85 | +0.0316 | 0/3; 0/3; 0/3 |
| Wu | GCN / IDW | -0.0317 | +4.15 | +0.0082 | 1/3; 1/3; 0/3 |
| Wu | GraphSAGE / IDW | -0.0503 | +6.19 | +0.0137 | 1/3; 0/3; 0/3 |
| Bassiouni | GCN / local | +0.0481 [+0.0408, +0.0555] | -17.84 [-18.64, -13.89] | -0.0295 [-0.0342, -0.0249] | 22/22; 22/22; 22/22 |
| Bassiouni | GraphSAGE / local | +0.0368 [+0.0329, +0.0460] | -13.27 [-15.02, -11.91] | -0.0216 [-0.0259, -0.0199] | 22/22; 22/22; 22/22 |
| Bassiouni | local / IDW | -0.0816 [-0.0974, -0.0613] | +26.01 [+21.13, +33.41] | +0.0427 [+0.0359, +0.0490] | 0/22; 0/22; 0/22 |
| Bassiouni | GCN / IDW | -0.0297 [-0.0353, -0.0243] | +4.64 [+3.53, +7.05] | +0.0124 [+0.0110, +0.0169] | 0/22; 1/22; 0/22 |
| Bassiouni | GraphSAGE / IDW | -0.0428 [-0.0568, -0.0277] | +7.92 [+6.11, +14.10] | +0.0182 [+0.0170, +0.0259] | 0/22; 0/22; 0/22 |

### Phikon

| Cohort | Candidate / reference | Median Δr | Median relative MSE (%) | Median ΔMAE | Improved patients (r/MSE/MAE) |
|---|---|---:|---:|---:|---|
| Development | GCN / local | -0.0211 | +13.32 | +0.0218 | 0/1; 0/1; 0/1 |
| Development | GraphSAGE / local | -0.0073 | +3.56 | +0.0052 | 0/1; 0/1; 0/1 |
| Development | local / IDW | +0.0548 | -22.54 | -0.0458 | 1/1; 1/1; 1/1 |
| Development | GCN / IDW | +0.0338 | -12.23 | -0.0240 | 1/1; 1/1; 1/1 |
| Development | GraphSAGE / IDW | +0.0475 | -19.78 | -0.0407 | 1/1; 1/1; 1/1 |
| Wu | GCN / local | +0.0045 | -0.46 | -0.0015 | 3/3; 3/3; 3/3 |
| Wu | GraphSAGE / local | +0.0006 | -0.09 | -0.0005 | 3/3; 2/3; 3/3 |
| Wu | local / IDW | -0.0042 | +0.49 | +0.0026 | 1/3; 1/3; 0/3 |
| Wu | GCN / IDW | +0.0004 | -0.04 | +0.0010 | 2/3; 2/3; 0/3 |
| Wu | GraphSAGE / IDW | -0.0040 | +0.49 | +0.0019 | 1/3; 1/3; 0/3 |
| Bassiouni | GCN / local | +0.0075 [+0.0054, +0.0089] | -1.70 [-2.08, -1.33] | -0.0045 [-0.0057, -0.0039] | 22/22; 22/22; 22/22 |
| Bassiouni | GraphSAGE / local | +0.0018 [-0.0000, +0.0056] | -0.67 [-0.97, -0.52] | -0.0019 [-0.0024, -0.0016] | 15/22; 20/22; 22/22 |
| Bassiouni | local / IDW | -0.0156 [-0.0272, -0.0128] | +1.91 [+1.39, +2.82] | +0.0070 [+0.0059, +0.0090] | 0/22; 0/22; 0/22 |
| Bassiouni | GCN / IDW | -0.0097 [-0.0199, -0.0065] | +0.37 [+0.07, +0.59] | +0.0028 [+0.0020, +0.0034] | 1/22; 4/22; 0/22 |
| Bassiouni | GraphSAGE / IDW | -0.0147 [-0.0228, -0.0096] | +1.27 [+0.90, +1.92] | +0.0051 [+0.0042, +0.0072] | 0/22; 0/22; 0/22 |

## Wu individual primary effects

The Wu medians above are descriptive summaries of only three patients, without population intervals. All three patient effects are given here; training and mask seeds were averaged within each patient. Positive Pearson and negative errors favor the candidate. Bassiouni table brackets above are descriptive 95% patient-bootstrap intervals.

| Backbone | Patient | Route | Candidate / reference | Δr | Relative MSE (%) | ΔMAE |
|---|---|---|---|---:|---:|---:|
| ResNet18 | 1160920F | Residual | GCN / local | +0.0415 | -10.66 | -0.0176 |
| ResNet18 | CID4290 | Residual | GCN / local | +0.0455 | -15.43 | -0.0262 |
| ResNet18 | CID44971 | Residual | GCN / local | +0.0633 | -14.29 | -0.0234 |
| ResNet18 | 1160920F | Residual | GCN / IDW | -0.0317 | +4.15 | +0.0088 |
| ResNet18 | CID4290 | Residual | GCN / IDW | +0.0212 | -0.33 | +0.0075 |
| ResNet18 | CID44971 | Residual | GCN / IDW | -0.0317 | +5.24 | +0.0082 |
| ResNet18 | 1160920F | Image | GCN / local | +0.0239 | -32.48 | -0.0892 |
| ResNet18 | CID4290 | Image | GCN / local | +0.0874 | +3.73 | +0.0152 |
| ResNet18 | CID44971 | Image | GCN / local | +0.0621 | -10.81 | -0.0087 |
| ResNet18 | 1160920F | Residual | GraphSAGE / local | +0.0228 | -8.91 | -0.0140 |
| ResNet18 | CID4290 | Residual | GraphSAGE / local | +0.0357 | -12.29 | -0.0198 |
| ResNet18 | CID44971 | Residual | GraphSAGE / local | +0.0411 | -11.05 | -0.0180 |
| ResNet18 | 1160920F | Residual | GraphSAGE / IDW | -0.0503 | +6.19 | +0.0123 |
| ResNet18 | CID4290 | Residual | GraphSAGE / IDW | +0.0114 | +3.37 | +0.0139 |
| ResNet18 | CID44971 | Residual | GraphSAGE / IDW | -0.0540 | +9.22 | +0.0137 |
| ResNet18 | 1160920F | Image | GraphSAGE / local | -0.0217 | -4.88 | -0.0098 |
| ResNet18 | CID4290 | Image | GraphSAGE / local | +0.0555 | -2.88 | -0.0034 |
| ResNet18 | CID44971 | Image | GraphSAGE / local | +0.0225 | -7.16 | -0.0123 |
| ResNet18 | 1160920F | Residual | local / IDW | -0.0732 | +16.58 | +0.0264 |
| ResNet18 | CID4290 | Residual | local / IDW | -0.0243 | +17.85 | +0.0337 |
| ResNet18 | CID44971 | Residual | local / IDW | -0.0951 | +22.79 | +0.0316 |
| Phikon | 1160920F | Residual | GCN / local | +0.0047 | -0.46 | -0.0015 |
| Phikon | CID4290 | Residual | GCN / local | +0.0005 | -0.28 | -0.0025 |
| Phikon | CID44971 | Residual | GCN / local | +0.0045 | -0.53 | -0.0014 |
| Phikon | 1160920F | Residual | GCN / IDW | -0.0044 | +0.17 | +0.0010 |
| Phikon | CID4290 | Residual | GCN / IDW | +0.0045 | -0.52 | +0.0013 |
| Phikon | CID44971 | Residual | GCN / IDW | +0.0004 | -0.04 | +0.0007 |
| Phikon | 1160920F | Image | GCN / local | -0.0032 | +0.40 | +0.0016 |
| Phikon | CID4290 | Image | GCN / local | +0.0020 | +0.45 | +0.0021 |
| Phikon | CID44971 | Image | GCN / local | +0.0153 | -0.19 | +0.0005 |
| Phikon | 1160920F | Residual | GraphSAGE / local | +0.0037 | -0.15 | -0.0007 |
| Phikon | CID4290 | Residual | GraphSAGE / local | +0.0006 | -0.09 | -0.0005 |
| Phikon | CID44971 | Residual | GraphSAGE / local | +0.0002 | +0.04 | -0.0003 |
| Phikon | 1160920F | Residual | GraphSAGE / IDW | -0.0055 | +0.49 | +0.0019 |
| Phikon | CID4290 | Residual | GraphSAGE / IDW | +0.0046 | -0.33 | +0.0032 |
| Phikon | CID44971 | Residual | GraphSAGE / IDW | -0.0040 | +0.53 | +0.0018 |
| Phikon | 1160920F | Image | GraphSAGE / local | -0.0088 | +9.39 | +0.0107 |
| Phikon | CID4290 | Image | GraphSAGE / local | +0.0163 | +11.67 | +0.0092 |
| Phikon | CID44971 | Image | GraphSAGE / local | +0.0036 | +6.01 | +0.0019 |
| Phikon | 1160920F | Residual | local / IDW | -0.0092 | +0.64 | +0.0026 |
| Phikon | CID4290 | Residual | local / IDW | +0.0040 | -0.24 | +0.0037 |
| Phikon | CID44971 | Residual | local / IDW | -0.0042 | +0.49 | +0.0021 |

## Aggregation and topology controls

These Bassiouni primary-condition contrasts assess specificity separately from the graph-versus-local question. Positive Pearson and negative errors favor the candidate. Mixed directions are tradeoffs.

| Backbone | Route | Candidate / reference | Median Δr | Relative MSE (%) | ΔMAE | Improved patients (r/MSE/MAE) |
|---|---|---|---:|---:|---:|---|
| ResNet18 | Image | GCN / mean | -0.0068 [-0.0175, +0.0005] | +0.38 [-1.63, +3.27] | -0.0042 [-0.0080, +0.0026] | 7/22; 11/22; 15/22 |
| ResNet18 | Image | GCN / GCN shuffled | +0.0382 [+0.0185, +0.0664] | -1.37 [-5.46, +3.55] | +0.0004 [-0.0092, +0.0118] | 18/22; 14/22; 11/22 |
| ResNet18 | Image | GraphSAGE / mean | -0.0135 [-0.0617, +0.0029] | +6.60 [-0.43, +15.98] | -0.0032 [-0.0139, +0.0215] | 8/22; 7/22; 13/22 |
| ResNet18 | Image | GraphSAGE / GraphSAGE shuffled | +0.0041 [-0.0078, +0.0121] | -7.69 [-10.01, -5.40] | -0.0096 [-0.0170, -0.0065] | 14/22; 22/22; 20/22 |
| ResNet18 | Residual | GCN / mean | -0.0019 [-0.0037, -0.0006] | +0.31 [+0.26, +0.71] | +0.0005 [+0.0004, +0.0009] | 2/22; 0/22; 0/22 |
| ResNet18 | Residual | GCN / GCN shuffled | +0.0019 [-0.0060, +0.0058] | +1.13 [-0.23, +1.81] | +0.0011 [-0.0007, +0.0022] | 12/22; 8/22; 7/22 |
| ResNet18 | Residual | GraphSAGE / mean | -0.0132 [-0.0217, -0.0083] | +5.30 [+2.84, +6.07] | +0.0080 [+0.0056, +0.0094] | 1/22; 0/22; 0/22 |
| ResNet18 | Residual | GraphSAGE / GraphSAGE shuffled | +0.0237 [+0.0178, +0.0290] | -6.84 [-8.66, -5.63] | -0.0105 [-0.0140, -0.0096] | 22/22; 22/22; 22/22 |
| Phikon | Image | GCN / mean | -0.0013 [-0.0022, +0.0003] | +0.04 [+0.01, +0.07] | +0.0001 [+0.0000, +0.0001] | 8/22; 5/22; 6/22 |
| Phikon | Image | GCN / GCN shuffled | +0.0165 [+0.0065, +0.0331] | +3.37 [+0.65, +4.86] | -0.0018 [-0.0070, +0.0055] | 17/22; 6/22; 12/22 |
| Phikon | Image | GraphSAGE / mean | -0.0082 [-0.0166, -0.0008] | +8.32 [+5.68, +10.00] | +0.0027 [-0.0020, +0.0142] | 6/22; 1/22; 10/22 |
| Phikon | Image | GraphSAGE / GraphSAGE shuffled | +0.0112 [+0.0070, +0.0214] | -0.34 [-1.38, +0.95] | +0.0013 [-0.0011, +0.0047] | 18/22; 12/22; 10/22 |
| Phikon | Residual | GCN / mean | -0.0004 [-0.0009, -0.0000] | +0.09 [+0.07, +0.15] | +0.0002 [+0.0001, +0.0004] | 7/22; 1/22; 0/22 |
| Phikon | Residual | GCN / GCN shuffled | +0.0009 [-0.0003, +0.0027] | -0.27 [-0.43, -0.17] | -0.0005 [-0.0009, -0.0003] | 14/22; 19/22; 21/22 |
| Phikon | Residual | GraphSAGE / mean | -0.0070 [-0.0078, -0.0032] | +1.20 [+0.84, +1.69] | +0.0030 [+0.0021, +0.0043] | 1/22; 0/22; 0/22 |
| Phikon | Residual | GraphSAGE / GraphSAGE shuffled | +0.0021 [+0.0014, +0.0040] | -0.61 [-0.80, -0.37] | -0.0012 [-0.0019, -0.0009] | 21/22; 22/22; 22/22 |

## Common-valid-gene sensitivity

Primary definitions remain unchanged. The following separately labelled sensitivity intersects the finite-correlation genes of each paired run before patient averaging. Gene counts are patient means over mask/training seeds, not extra biological replicates.

| Backbone | Route | Candidate / reference | Primary Δr | Common-valid Δr | Mean valid genes, patient range |
|---|---|---|---:|---:|---|
| ResNet18 | Image | GCN / local | +0.0079 | +0.0079 | 46.67–50.00 |
| ResNet18 | Image | GraphSAGE / local | -0.0089 | -0.0089 | 46.67–50.00 |
| ResNet18 | Residual | GCN / local | +0.0481 | +0.0481 | 46.67–50.00 |
| ResNet18 | Residual | GCN / IDW | -0.0297 | -0.0202 | 44.67–50.00 |
| ResNet18 | Residual | GraphSAGE / local | +0.0368 | +0.0368 | 46.67–50.00 |
| ResNet18 | Residual | GraphSAGE / IDW | -0.0428 | -0.0303 | 44.67–50.00 |
| Phikon | Image | GCN / local | +0.0007 | +0.0007 | 46.67–50.00 |
| Phikon | Image | GraphSAGE / local | -0.0016 | -0.0016 | 46.67–50.00 |
| Phikon | Residual | GCN / local | +0.0075 | +0.0075 | 46.67–50.00 |
| Phikon | Residual | GCN / IDW | -0.0097 | -0.0010 | 44.67–50.00 |
| Phikon | Residual | GraphSAGE / local | +0.0018 | +0.0018 | 46.67–50.00 |
| Phikon | Residual | GraphSAGE / IDW | -0.0147 | -0.0060 | 44.67–50.00 |

## Controls, uncertainty and complete supporting results

`PAIRED_SUMMARY.csv` includes all prespecified spatial/mean/shuffled, molecular, route and budget contrasts with descriptive 10,000-patient-resampling intervals in Bassiouni. `PAIRED_PATIENT_EFFECTS.csv` retains all three individual Wu effects and the single development patient; cohorts are never pooled. `DEFINED_GENE_COUNTS.csv` and the separately labelled common-valid-gene sensitivity preserve undefined-correlation information. Full per-run/per-gene data or checksum-verified prior pointers remain available. Fixed fusion, 10/40% coverage, identical-H placement and historical central-V2 are secondary conditions, not outcome-selected alternatives.

## Provenance and limits

The model manifest records every seed, state hash, selected/run epoch, update count and actual observed cost. Each configuration received the same two learning-rate trials, selected on mean TENX95 MSE across all three seeds before new selected outcomes. The validation-guided refinement and initial early-overfitting diagnostic remain documented. Head capacities are approximately matched within each backbone, not between encoders. Graph operators are fixed and encoders are frozen. Cross-backbone differences also include pretraining, feature dimension and effective field of view. The prior 48 advanced fits and original 30 minibatch fits remain intact, with their distinct update schedules. Original primary image-head binaries remain missing. All new full exports/states are author-local; no public release or journal submission occurred.
