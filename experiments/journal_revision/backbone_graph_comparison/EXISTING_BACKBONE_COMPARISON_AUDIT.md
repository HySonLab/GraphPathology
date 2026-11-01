# Existing backbone comparison audit

Source snapshot: ab0811134bb1d62719206ab11a43bd7c4dd5d813. Contract audit precedes new matched outcome inspection. Exact raw records and checksums are in REUSED_FIT_MANIFEST.json; source protocols and original outputs remain immutable.

## Usable matched ResNet18 configurations

Every row below uses frozen torchvision ResNet18 IMAGENET1K_V1 (official f37072fd weights,512-D), official resize256/center-crop224 normalization, identical registered112um-field RGB patches, source-only scaler fit on NCBI783/NCBI785, C0-normalized50-gene targets without another head target scaler, AdamW0.001/decay0.0001, pooled512-row minibatches/max80/patience10/min-delta0.0001. TENX95 alone selects minimum MSE. All external heads receive the frozen C9 T1 map once. The exact source/development/external masks0/1/2 and all three routes are available in extension evaluation bundles. Those are distinct from original ten-mask headline analyses. Checkpoint, prediction, support and export hashes are checked before reuse; no arrays will be invented.

| Configuration | Seed | Parameters | Selected/run epochs | Checkpoint | Reuse |
|---|---:|---:|---:|---|---|
| local / local | 0 | 72114 | 12/22 | `checkpoints/journal_revision/hest_learned/journal_v2_local_local_seed0.pt` | REUSED_VERIFIED |
| MEAN / spatial | 0 | 72114 | 12/22 | `checkpoints/journal_revision/hest_learned/journal_v2_mean_spatial_seed0.pt` | REUSED_VERIFIED |
| GCN / spatial | 0 | 72114 | 12/22 | `checkpoints/journal_revision/hest_learned/journal_v2_gcn_spatial_seed0.pt` | REUSED_VERIFIED |
| GCN / shuffled | 0 | 72114 | 1/11 | `checkpoints/journal_revision/hest_learned/journal_v2_gcn_shuffled_seed0.pt` | REUSED_VERIFIED |
| SAGE / spatial | 0 | 72075 | 7/17 | `checkpoints/journal_revision/hest_learned/journal_v2_sage_spatial_seed0.pt` | REUSED_VERIFIED |
| SAGE / shuffled | 0 | 72075 | 9/19 | `checkpoints/journal_revision/hest_learned/journal_v2_sage_shuffled_seed0.pt` | REUSED_VERIFIED |
| local / local | 1 | 72114 | 14/24 | `checkpoints/journal_revision/hest_learned/journal_v2_local_local_seed1.pt` | REUSED_VERIFIED |
| MEAN / spatial | 1 | 72114 | 16/26 | `checkpoints/journal_revision/hest_learned/journal_v2_mean_spatial_seed1.pt` | REUSED_VERIFIED |
| GCN / spatial | 1 | 72114 | 14/24 | `checkpoints/journal_revision/hest_learned/journal_v2_gcn_spatial_seed1.pt` | REUSED_VERIFIED |
| GCN / shuffled | 1 | 72114 | 9/19 | `checkpoints/journal_revision/hest_learned/journal_v2_gcn_shuffled_seed1.pt` | REUSED_VERIFIED |
| SAGE / spatial | 1 | 72075 | 10/20 | `checkpoints/journal_revision/hest_learned/journal_v2_sage_spatial_seed1.pt` | REUSED_VERIFIED |
| SAGE / shuffled | 1 | 72075 | 18/28 | `checkpoints/journal_revision/hest_learned/journal_v2_sage_shuffled_seed1.pt` | REUSED_VERIFIED |
| local / local | 2 | 72114 | 17/27 | `checkpoints/journal_revision/hest_learned/journal_v2_local_local_seed2.pt` | REUSED_VERIFIED |
| MEAN / spatial | 2 | 72114 | 14/24 | `checkpoints/journal_revision/hest_learned/journal_v2_mean_spatial_seed2.pt` | REUSED_VERIFIED |
| GCN / spatial | 2 | 72114 | 14/24 | `checkpoints/journal_revision/hest_learned/journal_v2_gcn_spatial_seed2.pt` | REUSED_VERIFIED |
| GCN / shuffled | 2 | 72114 | 15/25 | `checkpoints/journal_revision/hest_learned/journal_v2_gcn_shuffled_seed2.pt` | REUSED_VERIFIED |
| SAGE / spatial | 2 | 72075 | 7/17 | `checkpoints/journal_revision/hest_learned/journal_v2_sage_spatial_seed2.pt` | REUSED_VERIFIED |
| SAGE / shuffled | 2 | 72075 | 7/17 | `checkpoints/journal_revision/hest_learned/journal_v2_sage_shuffled_seed2.pt` | REUSED_VERIFIED |

## Other existing configurations and limits

- Original three selected ResNet18 width256 image heads: historical primary ensemble,144178 parameters each; binaries remain missing. Exported original results remain preserved, but these are not substitutes for the matched width128 minibatch comparison.
- Later full-section ResNet18 local/mean/GCN/SAGE fits: all48 fits and advanced-model controls remain retained. One update per training section per epoch differs from pooled minibatches; these are not optimization-matched to the18 reused states.
- Historical C4 Phikon local MLP/ridge and C6 spatial sensitivity: final-LayerNorm CLS768-D, pinned owkin/phikon057cc0295895c2df3dd7681a89680da6015cbefe, official224px ViTImageProcessor/ImageNet normalization. The local selected width256 and an additional head-target scaler used a validation-Pearson selection grid; spatial widths128/67 were inherited. Development/conditional holdout exports, preprocessing and extraction registries exist, but exact embedding arrays and original primary model binaries are absent. These are valid historical development evidence, not the new matched contract or external predictions. Their registered FOV was approximately112um; no native-resolution recovery is claimed.
- New Phikon matched configurations: local,MEAN,GCN,SAGE,GCN-shuffled,SAGE-shuffled × seeds0/1/2. Source-only scaler is separately fit on the new768-D embeddings, direct common C0 targets, same pooled optimization/selection opportunities as reused ResNet. Local/MEAN/GCN128 (104882 parameters), SAGE66 (104792). Inputs and weights are hash-verified; these are NEW_FITS_REQUIRED until completed, never recovered historical heads.

## Common comparison contract

Fit NCBI783/NCBI785, select TENX95, report TENX99 separately, Wu1160920F/CID4290/CID44971 and the locked22 Bassiouni patients separately. No external expression determines tuning or calibration. Predictions/routes: image-only; fixed fusion0.2g+0.8IDW separately; observed-residual g+IDW(y_O-g_O). Use identical 10/20/40% historical neural masks and B2 supports across all comparators; primary20% distributed. The molecular IDW baseline always remains in anchor-assisted comparisons. All per-run/per-gene/prediction artifacts are identified in the study manifests. Encoder changes include pretraining, dimension and transforms; only within-backbone contrasts isolate added context. Fixed aggregation+MLP heads are not end-to-end encoder training or new message-passing architectures.
