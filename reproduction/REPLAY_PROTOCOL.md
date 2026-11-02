# Locked post-publication replay protocol

Locked before any new training or fresh outcome inspection on this branch. This is an informed post-publication reproducibility exercise, not blind validation. All existing paper/results files are read-only references. No new manuscript build or experimental search is permitted.

## Recovery and budget

The three original primary states are missing. Only three replacement fits are authorized by this protocol: seeds 0, 1 and 2, CPU float32, four Torch/BLAS threads, at most 80 epochs each; stop after 10 epochs without validation-MSE improvement of at least 0.0001. No retries, architecture/learning-rate search or external-outcome tuning. Available journal states are reused, not retrained. A failed fit is retained and reported. Resource budget: three hours CPU wall time for fitting and verification; no external computing purchase.

Fit only NCBI783 then NCBI785 in original barcode order; stop only on TENX95. TENX99 is development evaluation and Wu/Bassiouni are evaluation only. The original ordered 50 genes, C0 log1p normalization, persisted 512-dimensional feature scaler and 50-dimensional head-target scaler govern fitting. Audit source files, barcodes, targets, coordinate/encoder provenance and scaler consistency before fitting. Cached frozen ResNet18 features may be reused only after their recorded hashes/provenance pass. Independently recompute source target normalization from original counts, avoiding a premature float32 rounding of head targets. Record any historical input reconstruction differences.

The head is Linear(512,256), exact GELU, Dropout(0.1), Linear(256,50). AdamW lr0.001, decay0.0001, default betas/epsilon, all trainable parameters; shuffle with torch.randperm each epoch, batches512 retaining partial batch, MSEmean, no scheduler. Use the maintained historical-contract core.train_mlp. Save the validation-selected state, never force historical epochs12/11/12. Preserve all histories. Inverse head-target scaling per seed, average in float64, then cast ensemble to float32. Use archived C9 calibration once inside external.predict; no refit or Bassiouni calibration.

## Tier A: strict deterministic replay

Recompute molecular interpolation and the twenty-estimator V1 control grid, original V2 geometry checks and identical-hidden-target B2 designs with the existing journal_anchors runner on the three Wu and22 Bassiouni sections, budgets10/20/40%, ten frozen mask seeds. Write only to a new reproduction run directory. Match historical random-mask molecular primary Pearson/MSE/MAE with locked rtol=atol=1e-10; exact NaN/random-hidden-target input invariance. Spearman is explicitly secondary and may differ with floating-point near ties; do not claim exact ranks. Recompute B2 fixed-H central-placement paired patient effects against original derived exports with the same primary tolerances. No changed tolerance after seeing outcomes.

## Tier B: new primary fits

Evaluate all25 fixed external sections, ten original V1 random trajectories,20% anchors, image-only (M0) and fixed calibrated20% image fusion (M3); also retain individual seed outcomes. Compare each run to the corresponding original C9/C10 frozen rows. Report per-patient means across masks, then cohort medians, plus seed variation. Practical agreement requires ALL: absolute cohort median Pearson difference <=0.05; absolute relative cohort MSE difference <=15%; absolute relative MAE difference <=15%. These are exploratory practical criteria, not equivalence tests. Report negative findings and all failures. TENX99 whole-section source-unit image performance is a separate development check, not an external selection signal. Do not tune against any of these outcomes.

## Tier C: recovered states

Check all30 minibatch,48 advanced and72 candidate matched state references (132 distinct best states), their recorded hashes/configs/selection links, safe weights_only loading and bounded finite real-feature forwards. Replay an actual measured patient with original masks through each family; distinguish subset replay from whole-cohort reproduction. Reuse recovered selected RN/Phikon heads to replay the22-patient, three-mask,20%-anchor, residual GCN-versus-local comparisons only if feasible within verification budget; otherwise mark those cohort targets PARTIAL rather than claiming them reproduced from CSV regeneration.

## Reports, preservation and releases

Status labels are EXACT_REPLAY, APPROXIMATE_RETRAIN, PARTIAL and NOT_REPRODUCED, scoped per finding. Comparisons retain paper/reference values separately from regenerated values with signed, absolute and relative differences and aggregation. Preserve frozen file hashes and git trees before/after; no original identity is assigned to new weights. Inventory missing historical supporting states without needlessly retraining them. Publish only rights-cleared, verified new/recovered artifacts; hold uncertain derivative redistribution in an explicitly explained draft. No raw patient images/counts or third-party pretrained encoders in assets. Upload a draft, verify size/hash by downloading actual GitHub assets into a clean checkout, run safe load and bounded real inference, then publish permitted assets and verify anonymous URLs. New code/docs PR stays unmerged; an immutable tag may reference its tested commit, disclosed explicitly.
