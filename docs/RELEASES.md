# Verified research artifacts v1

[Public research release](https://github.com/HySonLab/GraphPathology/releases/tag/research-artifacts-v1) —  15 assets; 706,647,760 bytes. All asset sizes/SHA256 values match local staging and GitHub digests; all 15 URLs responded without credentials. The primary-head ZIP, README and verification receipt were additionally downloaded anonymously and rehashed.

Code tag: research-artifacts-v1 → 44918ef8f9e9876234ceb02a8b4a815d8b4072d4. This is the tested commit on the separate, unmerged [PR #4](https://github.com/HySonLab/GraphPathology/pull/4), not a change to main. Later PR commits add download/publication receipts and the auditable clean-checkout driver; the packaged inference/training code is unchanged.

Public artifacts contain 99 distinct best states: three RETRAINED_NOT_ORIGINAL primary replacements and 96 recovered ResNet18-based exploratory states. Forty-eight resume files are additional files, not extra fitted models. No raw patient data, pretrained encoder or Phikon head is included. The original primary identities remain missing, and the original doctor still fails as expected.

## Downloads

- [graphpathology-retrained-primary-resnet18-heads-v1.zip](https://github.com/HySonLab/GraphPathology/releases/download/research-artifacts-v1/graphpathology-retrained-primary-resnet18-heads-v1.zip) — 1,635,176 bytes.
- [graphpathology-exploratory-neural-heads-v1.zip](https://github.com/HySonLab/GraphPathology/releases/download/research-artifacts-v1/graphpathology-exploratory-neural-heads-v1.zip) — 114,543,063 bytes.
- [graphpathology-source-data-results-v1.zip](https://github.com/HySonLab/GraphPathology/releases/download/research-artifacts-v1/graphpathology-source-data-results-v1.zip) — 587,712,241 bytes.
- [graphpathology-model-and-results-manifest-v1.json](https://github.com/HySonLab/GraphPathology/releases/download/research-artifacts-v1/graphpathology-model-and-results-manifest-v1.json) — 738,362 bytes.
- [SHA256SUMS.txt](https://github.com/HySonLab/GraphPathology/releases/download/research-artifacts-v1/SHA256SUMS.txt) — 1,230 bytes.
- [RELEASE_VERIFICATION.json](https://github.com/HySonLab/GraphPathology/releases/download/research-artifacts-v1/RELEASE_VERIFICATION.json) — 24,786 bytes.

The complete 15-asset inventory, stable URLs, sizes and SHA256 hashes— including the checksum index and verification receipt themselves— is in reproduction/RELEASE_ASSETS.json. SHA256SUMS.txt covers the initial archives/metadata/manifest; it cannot hash itself and predates the later receipt. The receipt’s own hash is verified against GitHub’s digest and recorded in RELEASE_ASSETS.json. Checksums are integrity checks, not author signatures.

## Actual clean-checkout test

All 14 initial assets (706,622,974 bytes) were downloaded from the actual authenticated draft asset IDs into a separate directory. A fresh GitHub clone of the exact tag restored the downloaded archives, preflighting every member and refusing existing conflicts. All 147 public PT files loaded with weights_only=True; all 96 recovered public heads performed real 554-location BASS_P19 inference. The three replacement heads additionally ran on eight real HEST image patches using a fresh official 46,830,571-byte ResNet18 download; its SHA256 matched and extracted features were byte-identical to the recovered reference. All 83,136 frozen files and original scientific Git trees passed in that checkout; tracked files remained unchanged.

The same pinned Python environment on the same host was reused; this is not a claim of independent hardware/environment or newly acquired evaluation data. Exact inputs, file hashes, outputs and scope are in reproduction/CLEANROOM_VERIFICATION.json and RELEASE_VERIFICATION.json. The auditable driver is reproduction/CLEANROOM_PROBE.py; point PYTHONPATH at the tagged checkout and supply the documented registered local inputs. Dataset/encoder acquisition is separate from artifact redistribution.

## Restore

Follow reproduction/README_RELEASE.md and requirements-v1.txt. From a fresh editable checkout, download the exact release, check SHA256SUMS.txt, then run:

```bash
python -m graphpathology.reproduction_release restore \
  --assets /absolute/path/to/downloads \
  --manifest /absolute/path/to/downloads/graphpathology-model-and-results-manifest-v1.json
```

Use the explicit reproduction_primary route for replacement heads; do not rename them into historical original paths. Read the original ordered genes/scalers/configuration from the supplied manifest. Fresh run identities are required for any rerun.

## Phikon draft and unresolved rights

[Metadata-only Phikon rights-review draft](https://github.com/HySonLab/GraphPathology/releases/tag/untagged-84e85d779dcbe207b3ab) (authenticated maintainer access; intentionally unpublished). Its only asset is PHIKON_RIGHTS_REVIEW.md. Thirty-six verified head binaries and their derivative numerical archives remain local; none was uploaded to this draft or the public release. Owkin’s custom nonprofit/transfer/same-terms Results provisions and HEST’s separate share-alike conditions require clarified compatible onward-distribution rights. See reproduction/LICENSE_AUDIT.md for exact pinned sources and hashes. No blanket software or manuscript license was changed.

## Scientific preservation

No manuscript, supplementary material, figure, table, bibliography, PDF, historical result or reported number changed. Original trees are byte-identical; main remains 67738be6f0aaace07b55afac44db3d3d23c41c9b. The PR stays unmerged. Exact deterministic replay, approximate new fitting, scoped artifact checks, seed variation, source/runtime sensitivity and failures are separately documented in reproduction/REPRODUCTION_REPORT.md, RESULTS_COMPARISON.csv, INVENTORY.csv and verification/. No paid computing service was used.
