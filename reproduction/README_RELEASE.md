# GraphPathology research artifacts v1

This package separates missing original weights from verified post-v1 and newly retrained models. The immutable research-artifacts-v1 tag identifies the tested code on an unmerged review branch. Main and every manuscript/supplement/figure/table remain unchanged.

- Three primary ResNet18-head replacements are RETRAINED_NOT_ORIGINAL, seeds0/1/2, validation-selected epochs12/11/12. Original checkpoint paths/hashes remain missing.
- The exploratory head asset contains96 distinct recovered ResNet18-based best states:30 minibatch,48 advanced and18 additional matched-rate trial heads. Eighteen matched selected states reuse minibatch binaries; indexes distinguish references from unique states. The48 resume files are not extra fits.
- Thirty-six Phikon heads and associated derivative result archives are held locally pending the custom-license review. The Phikon draft has metadata only. No Phikon encoder/head binary is in this public package.
- Numerical bundles contain derived per-gene/panel metrics and new replay records, not raw patient tissue/counts or pretrained encoders. Frozen source exports retain their bytes and namespaces. New reruns use reproduction/runs/.

Read REPRODUCTION_REPORT.md for exact versus approximate versus partial scopes, RESULTS_COMPARISON.csv for numerical differences, and LICENSE_AUDIT.md for artifact-specific rights. The new training took place only once; no evaluation-driven tuning was performed. Tiny thread-count sensitivity outcomes are retained. The molecular/geometry full replay and the selected22-patient matched comparison pass strict full-precision checks; new primary ensemble metrics satisfy all practical agreement criteria.

## Restore and load

1. Clone HySonLab/GraphPathology at the research-artifacts-v1 tag, on a fresh directory, and install the documented Python dependencies with `python -m pip install -e '.[vision]'`. Use requirements-v1.txt to pin the recorded package versions; install the +cpu Torch packages from the official CPU index. The tag commit is specified in graphpathology-model-and-results-manifest-v1.json; it need not be main.
2. Download the assets attached to this exact release. In their download directory run `sha256sum -c SHA256SUMS.txt`. The manifest itself is covered by that checksum index; the index cannot hash itself. Compare the separately fetched GitHub asset digests/verification receipt as well. Checksums establish integrity, not a cryptographic author signature.
3. Run from the checkout:

```bash
python -m graphpathology.reproduction_release restore \
  --assets /absolute/path/to/downloads \
  --manifest /absolute/path/to/downloads/graphpathology-model-and-results-manifest-v1.json
```

The restore command verifies every outer/member hash and preflights all paths before writing; it refuses conflicting existing files and leaves identical files untouched. It cannot write manuscript paths. A model can be loaded only with safe `torch.load(...,weights_only=True)` and its indexed configuration/feature scaler. Do not load a primary replacement as an original v1 checkpoint.

For aligned frozen ResNet18 embeddings (N×512, no hidden expression), the explicit primary route is:

```bash
python -m graphpathology.reproduction_primary predict \
  --checkpoints checkpoints/reproduced/primary_resnet18/primary-20261101-v1 \
  --features /absolute/path/to/real_aligned_embeddings.npz \
  --output reproduction/runs/your-new-check/primary_predictions.npz
```

The NPZ needs an embeddings array in the original official preprocessing contract and gene/barcode alignment must be independently verified. The checkpoint folder includes the historical feature and head-target scalers, original50-gene order, source checksums, selected epochs and complete training histories. External C9 transport is separate and applied once by the documented evaluation route; no external refitting is allowed.

To retrain after reconstructing the registered HEST inputs, use the locked config and a NEW run identity:

```bash
python -m graphpathology.reproduction_primary train \
  --config configs/reproduction/primary-v1.json \
  --checkpoints checkpoints/reproduced/primary_resnet18/your-new-run \
  --output reproduction/runs/your-new-training-run
```

The original doctor --require-models still fails for missing historical primary weights; this is intentional. Opt-in reproduction succeeds independently. Download the official ResNet18 IMAGENET1K_V1 encoder separately from https://download.pytorch.org/models/resnet18-f37072fd.pth (SHA256 f37072fd47e89c5e827621c5baffa7500819f7896bbacec160b1a16c560e07ec) when extracting image features. No encoder is redistributed here. HEST source acquisitions and transforms are documented in docs/JOURNAL_REVISION.md and the immutable cache provenance; source samples remain NCBI783/NCBI785, validation TENX95, development-only TENX99.

The accompanying RELEASE_VERIFICATION.json records actual GitHub download/checksum/clean-checkout/real-image inference checks. Public URLs and the metadata-only Phikon draft are recorded in docs/RELEASES.md. The PR remains unmerged for author review.
