# Post-v1 bounded Transformer / EGNN extension

Continue draft PR #3; do not merge before author review. The local outcome-blind protocol is [PROTOCOL.md](../experiments/journal_revision/transformer_egnn/PROTOCOL.md), committed as `2e03fae636d6f6527d25db58f4cbab72f7c1a286` before fitting. This is `post_v1_exploratory`, not formal preregistration. All previous A–E results, 30 journal-v2 fits, original manuscript files and C0–C10 records remain unchanged.

The extension completed 48 fits (16 configurations × seeds 0/1/2), plus evaluation of those and the 30 frozen heads in 26 previously inspected patients: TENX99, Wu3 and Bassiouni22. All 2,028 model–patient bundles and 218,772 route/reference rows are retained. HEST training is restricted to NCBI783/NCBI785; TENX95 selects stopping. No evaluation outcome selected architecture, coefficients or seed. External images use the existing C9 T1 map without destination refitting. Training uses four CPU threads per process; the large external evaluation was distributed across CPU processes with one thread each, recorded per bundle. The standalone inference benchmark uses four threads on the same validation section after the experimental workers finished. No GPU or paid computing service was used.

The primary neural result is negative. At 20% distributed Bassiouni coverage, GPS full and latent EGNN residual routes changed median paired patient Pearson by −0.0504/−0.0480 and relative MSE by +13.4%/+13.0% versus IDW. All 22 patients deteriorated in Pearson/MSE/MAE for both primary models. Component ablations did not establish the hypothesized gains. This supports neither a stronger neural novelty claim nor a general conclusion that native-resolution geometry-aware methods cannot help. The complete source-only/low-resolution experiment is retained without selecting a winner.

The recovered 22-patient matched spatial baseline study is the strongest new evidence. All 88 GEO source files match historical SHA256 values. Selective HTTP206 TAR ranges transferred 291,569,893 bytes, under a 1GiB cap; no 37.7GB whole-series transfer occurred. All 19,800 A/B1/B2 runs and 39,600 exact corruption checks completed. At 20% identical-H/equal-count B2 coverage, every patient deteriorated under region exclusion for both fixed estimators. IDW median paired Δr was −0.1541 and MSE increased 45.0%; uniform Δr was −0.1502 and MSE increased 39.9%. Patient bootstrap intervals are descriptive (10,000 resamples, seed 20261101), with no spot/seed/region pseudoreplication. Wu effects remain individual; HEST has one evaluation patient.

Read the [scientific report](../experiments/journal_revision/transformer_egnn/SCIENTIFIC_REPORT.md), [experiment matrix](../experiments/journal_revision/transformer_egnn/summary/EXPERIMENT_MATRIX.csv), [all paired patient effects](../experiments/journal_revision/transformer_egnn/summary/PAIRED_PATIENT_EFFECTS.csv) and [validation report](../experiments/journal_revision/transformer_egnn/VALIDATION_REPORT.json). Full image/fixed/residual, topology/distance, seed, gene and occupied-distance-bin outputs are retained. Undefined correlations and raw Spearman near-tie sensitivity are explicit. The original selected primary heads remain missing; these 78 post-v1 states are not replacements.

## Reproduction

Imports/help/tests do not download or train. From an editable repository installation with the recorded environment and existing verified HEST source caches:

```bash
python -m graphpathology.extension preflight
python -m graphpathology.extension acquire-bassiouni
python -m graphpathology.extension prepare-wu
python -m graphpathology.extension prepare-bassiouni
python -m graphpathology.extension train-all --reproduction-output outputs/fresh-training-replay
python -m graphpathology.extension evaluate
python -m graphpathology.extension diagnostics
python -m graphpathology.extension report
python -m graphpathology.extension benchmark
python -m graphpathology.extension verify
```

For a fresh single-fit replay, use `train --variant gps_full --seed 0 --reproduction-output outputs/fresh-single-fit`. New replay directories link only the verified source caches/scaler and write new states/histories under ignored outputs; they never overwrite the committed fit records or original states. Without this option, `train` resumes an incomplete original run or verifies/skips an installed complete state. Restore the source feature scaler from the author bundle or original preparation pipeline before replay. For per-patient evaluation, use `evaluate --patient BASS_P01 --threads 1`; complete bundles are hash-verified and skipped under per-case locks. The recovered-cohort baseline runner uses the unchanged `journal_anchors.run_patient(root, geo_data_root, 'bassiouni', patient, new_output)` specification, followed by `summarize_all(new_output)`. Exact designated primary replay must pass before each patient's new baseline scoring. `extension_verify.audit_bassiouni` checks all 2,640 historical random-reference runs and 396,000 primary per-gene values. Full cohort image diagnostics use `extension_diagnostics.recovered_patient`; all recorded patient IDs come from the exact selected-source manifest.

Run the matched baseline cohort before the recovered-cohort neural comparisons, as recorded in [DATA_RECOVERY.md](../experiments/journal_revision/transformer_egnn/DATA_RECOVERY.md). Do not refit any scaler on evaluation patients. The nominal112um crop is followed by the original Resize256/CenterCrop224 encoder transform, retaining approximately0.875 of its native linear field. External native image resolution remains limited. Conditional E(2)/permutation diagnostics do not establish full CNN/image-pipeline equivariance. Actual rewiring acceptance counts delimit weak/infeasible degree–length controls. Both distance-fixed and latent models are E(2)-compatible, and the final latent-coordinate update has no downstream scalar loss.

## Author bundles and public checkout checks

Raw counts/tissue images/model binaries are not committed. High-volume per-gene result files are preserved outside Git in a 523MB author bundle, alongside the approximately101MB bundle of 48 new best/resume/optimizer states and unchanged source feature scaler. [RELEASE_READY_MANIFEST.json](../experiments/journal_revision/transformer_egnn/RELEASE_READY_MANIFEST.json) records every path, checksum and author-local output directory. No public release upload was performed. The author can review these concrete bundles before approving a licensed release/DOI destination. The older30-state bundle remains distinct.

Extract the full-gene ZIP into the repository root to restore exact experiment-relative paths. For a checkout without author-local weights/exports, use `verify --committed-only`: it verifies committed panels/bins/manifests and explicitly reports skipped, uninstalled author-local files. It does not claim those missing binaries were checked. Full `verify` requires the local complete bundles/caches. `bundle --output /chosen/author/directory` recreates release-ready artifacts without uploading.

The new journal version is under [manuscript/journal_revision/transformer_egnn](../manuscript/journal_revision/transformer_egnn/). Earlier journal files remain historical and byte-preserved. Rebuild the committed version without research downloads/training:

```bash
python -m graphpathology.extension paper --build-output outputs/journal-extension-rebuild
```

The clean source package includes the new drawing code, all needed figure-source tables and existing29 original publication targets. Its `scripts/regenerate_extension_figures.py` redraws the extension from exported tables without raw research data or checkpoints. Full fitting/evaluation still requires the repository and public inputs. New PDFs, source ZIP and release-ready bundles are author-review artifacts; PR #3 remains draft and unmerged.
