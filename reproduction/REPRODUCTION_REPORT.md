# GraphPathology artifact recovery and reproduction

The three original selected primary checkpoints remain missing. Three separately identified replacement heads were trained once, using the historical configuration and source-only splits. The external ensemble results satisfy every advance-declared agreement criterion; full-precision cohort metrics match the frozen records within floating-point precision. No manuscript, supplementary material, figure, table, bibliography, PDF or historical result was edited or rebuilt in place.

## What was recovered

Five author-local prepared archives were recovered and verified. All 30 minibatch, 48 advanced and 72 candidate matched state references are available. Eighteen matched ResNet references reuse minibatch states: there are 132 distinct best states, plus 48 optimizer/resume files. All were safely loaded and all best states performed bounded real-section inference. The original 3 primary heads and 1,379 supporting historical model/PCA entries remain explicitly missing. Only the primary ensemble was necessary to retrain. A 27-archive model-member search found no registered missing original; accessible Documents/Downloads, mounted-media locations, git model-object names, release/tag/CI listings were also checked. CI listings and their producing workflow definitions were inspected; their archive payloads were not exhaustively downloaded. No older drive was mounted and no private storage credentials were requested.

See INVENTORY.csv, RECOVERY_REPORT.md, ARCHIVE_SEARCH.json and RECOVERED_BUNDLES.json. Original expected hashes are preserved; replacement identities are never substituted into original doctor/infer routes.

## Training and evaluation

Protocol commit 1b19cf9d9 precedes all new fits. Source-code audit corrected the protocol's initial target-precision prose before fitting: the historical C0 transformation casts to float32 before the separate head-target scaler. NCBI783/NCBI785 provide 7,015 fitting locations; TENX95 provides 7,760 stopping locations. Source features, ordered genes/barcodes, raw data identities, 224-pixel RGB geometry and persisted scalers passed audit. Reconstructed feature-scaler statistics match exactly; the target-scaler standard-deviation difference is at most 1.54e-14. The persisted original scalers were used, not refitted replacements.

The head is 512→256→GELU→Dropout0.1→50, AdamW lr0.001/decay0.0001, batch512, max80, patience10/min_delta0.0001, CPU float32, deterministic algorithms, four threads. Seeds0/1/2 selected epochs12/11/12 independently. Validation median Pearson was0.6858333517/0.6797878971/0.6861210879. No external/development outcome selected weights, hyperparameters, calibration or stopping. Each head prediction is inverse target-scaled; the float64 ensemble mean is cast to float32. Archived C9 calibration is applied once to external fusion, without fitting Bassiouni.

Only one three-seed training run occurred. Fresh identities live under checkpoints/reproduced/primary_resnet18/primary-20261101-v1, with origin RETRAINED_NOT_ORIGINAL. Tensor-level agreement in a rerun is not recovery of an original checkpoint/container identity.

## Results by scope

| Finding | Actual rerun | Status |
|---|---:|---|
|20% distributed molecular correlation, Wu |0.3590441 (paper0.359)|EXACT_REPLAY against full-precision exports|
|20% distributed molecular correlation, Bassiouni |0.4158616 (paper0.416)|EXACT_REPLAY against full-precision exports|
|Bassiouni central fixed-H placement median Δr |-0.1540804 (paper-0.1541)|EXACT_REPLAY|
|Bassiouni placement relative MSE |+44.9864% (paper+45.0%)|EXACT_REPLAY|
|ResNet18 GCN-vs-local residual Δr / MSE |+0.0481 /-17.84%|EXACT_REPLAY of recovered weights|
|Phikon GCN-vs-local residual Δr / MSE |+0.0075 /-1.70%|EXACT_REPLAY locally; distribution review pending|

Tier A recomputed 22,500 runs across all 25 sections, all 20 V1 estimators, V2 geometries,10/20/40% budgets and B2 identical-H placements. All 3,375,000 primary per-gene Pearson/MSE/MAE values matched frozen full-precision references at rtol=atol=1e-10; 45,000 NaN/random-all-unobserved corruption checks passed. Spearman remains explicitly near-tie sensitive, with no exact secondary-rank claim.

Tier B evaluated image-only M0 and fixed calibrated fusion M3 on all 25 external sections, ten original V1 mask trajectories and 20% anchors, retaining each individual seed and the ensemble. All 12 ensemble cohort criteria passed (|Δr|≤0.05 and absolute relative MSE/MAE difference≤15%). Per-patient, per-mask, seed and cohort differences are in RESULTS_COMPARISON.csv and verification/tier-b-v2_*; the ensemble's largest cohort median-correlation delta is 5.56e-17. These are newly fitted models, labeled APPROXIMATE_RETRAIN even where exported metrics match. The weak image-only performance and the fusion disadvantage relative to molecular interpolation remain present. TENX99 is a separate development-only evaluation. Per-patient sample variance/std/range over the three training seeds is retained in verification/SEED_VARIATION.csv; seeds were not selected or treated as independent patients. Other historical neural protocols/LOPO/PCA branches were not retrained or claimed reproduced.

Tier C hash/config/selected-state/safe-load/real-forward checks cover 150 logical references and 132 unique best states on the actual 554-location BASS_P19 section. 114 references also match existing numerical exports at 1e-10; 36 unselected validation-only trial references have no historical external export and are marked PARTIAL/safe-forward-only. No such missing reference was fabricated. All 48 optimizer/resume files load safely. The separate 22-patient selected GCN/local replay computed 792 residual-route rows from real predictions, both backbones, three training seeds and three original neural masks; all pass strict archive comparisons.

The first completed recovered-head replay used four inference threads everywhere. Its float32 reduction differences reached 6.09e-8 in a panel correlation. Original extension records specify one inference thread, while later Phikon scoring used four. One corrective replay restored those immutable runtime settings; the initial outcomes remain as sensitivity evidence. No thresholds, data, models or weights changed. EXECUTION_ISSUES.md also records two pre-outcome input/reference-path assertion failures and their regression tests.

Displayed paper values are rounded. Exactness refers to full-precision archive checks, not an assertion that a rounded printed number agrees to 1e-10. This is informed post-publication reproducibility, not new blind validation or a scientific equivalence test.

## Resources and verification

Available hardware: 12 logical CPU cores, approximately 30 GiB RAM; CUDA unavailable. Fitting took 14.92s including input audit, with peak sampled RSS 583 MiB. Full deterministic replay took 521.22s; corrected matched replay 67.39s; corrected artifact verification 12.51s. Detailed per-stage times, hashes, environment and histories are retained in COMPUTE_SUMMARY.json and verification/. All research inputs were recovered locally, so no new research dataset was downloaded. License files and the official encoder used in the clean-checkout test are separately logged. No paid computing service was initiated.

The original 142-test baseline is preserved, with new tests for replacement identity, scaler/ensemble ordering, source-split protection, protected output paths, checkpoint tampering, preserved external float64 precision, original reference/runtime restoration and safe asset restoration. All 29 publication targets and 71 regression checks passed in a separate scratch output directory; no manuscript PDF was rebuilt. Preservation of 83,136 frozen scientific/history files and original Git trees is verified separately in PRESERVATION_REPORT.json.

## Distribution and limitations

The public-eligible package contains three replacement primary heads and 96 recovered ResNet18-based heads, with 48 additional resume files, scalers/configs/indexes, permitted derived numerical results and SHA256 manifests. Third-party pretrained encoders and raw patient data are excluded. Thirty-six verified Phikon heads and their derivative result archives remain local; a metadata-only draft records the specific licensing blocker described in LICENSE_AUDIT.md. This does not block the verified ResNet artifacts.

The research tag points to the tested code commit on the separate, unmerged PR branch, not an alteration of main. Actual URLs, upload/download checksums and the clean-checkout verification receipt are recorded in docs/RELEASES.md and RELEASE_VERIFICATION.json. Release publication does not imply that the PR has been reviewed or merged. No missing original checkpoint is represented as recovered.
