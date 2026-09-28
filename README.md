# GraphPathology

Reproducible, checkpointed research on spatial transcriptome reconstruction from histopathology.
This repository implements **C0: data, alignment, observation masks, and leakage audit** from the supplied September 2026 execution plan. C1-C8 have not been executed; there are no prediction-performance or biological-discovery results yet.

## Current development cohort

The revision-pinned public [HEST-Benchmark](https://huggingface.co/datasets/MahmoodLab/hest-bench) IDC task supplies four paired breast cancer H&E/Xenium slides. Patient identities and estimated microns per pixel come from pinned [HEST metadata](https://github.com/mahmoodlab/HEST/blob/3ddb5eaf5bd2a8133e0c0e8015816489a3d99dc3/assets/HEST_v1_1_0.csv).

| Split | Slides | HEST patient IDs |
|---|---|---|
| Train | NCBI783, NCBI785 | patient 5, patient 4 |
| Validation | TENX95 | patient 2 |
| Test | TENX99 | patient 1 |

This development split is intentionally small and differs from the official HEST cross-validation protocol. It is insufficient for publication-level patient statistics. Cohort/platform and source differences are confounded with split; evaluate this at C1 before expanding claims.

Patch barcodes are joined explicitly to expression barcodes. Patch coordinate order is row/column, while AnnData spatial coordinates are x/y. Expression locations without cached patches are excluded and counted. Negative controls and blank codewords are excluded from biological targets by fixed identifier prefixes; complete downloaded raw counts remain unchanged. Biological panel sizes also differ between slides; genes are intersected and reordered by identifier before training-only variance ranking. This harmonization uses assay gene identifiers, never held-out expression values.

## Reproduce

Run from `/home/son/Research/GraphPathology`:

```bash
python -m venv --system-site-packages .venv
.venv/bin/pip install -e '.[test]'
.venv/bin/pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu
.venv/bin/python scripts/download_c0.py
.venv/bin/python -m pytest -q
.venv/bin/python scripts/cache_features.py --threads 4
.venv/bin/python scripts/alignment_figures.py
.venv/bin/python scripts/patch_audit.py
.venv/bin/python scripts/mask_figures.py
.venv/bin/python scripts/verify_c0.py --include-data
```

`configs/environment.txt` records exact versions from the executed environment. Approximately 6 GB of disk is needed for downloaded assets, plus environment and derived caches. Downloads require network access; benchmark assets are ungated. No authentication or credential files are included in source.

The checked-in C0 freeze is immutable: `scripts/run_c0.py` refuses to overwrite an existing `C0_SPEC_LOCK.json`. It generated the committed report, normalized targets, panels, masks, manifests and lock after the visual audit. To create a revised scientific checkpoint, use a separately versioned configuration/output directory and review that revision; do not silently replace the current freeze. Verification can run without data; `--include-data` also hashes the downloaded assets and feature/target caches.

## Contracts

- `src/graphpathology/data.py`: strict reader, barcode alignment, physical coordinates, patient-disjoint splits, control-feature exclusion.
- `src/graphpathology/preprocessing.py`: training-only mean/scale of log1p counts; training-only variance-ranked panels. This is not a formal HVG algorithm. No test library-size normalization is performed; depth effects remain a C1 diagnostic.
- `src/graphpathology/masks.py`: fixed random, square block holdout, and circular capture-window masks at 5/10/20/40/60/80% nominal observed coverage, three seeds each. The whole target vector is hidden at a masked location.
- `scripts/cache_features.py`: frozen ResNet18 ImageNet embeddings, using official preprocessing; separate auditable RGB/texture statistics. These are a generic vision reference, not pathology foundation embeddings or validated cell/nuclear features.
- `scripts/run_c0.py`: checksums, statistics, frozen experiment specification, and hidden-target invariance audits.
- `tests/`: patient leakage, swapped coordinates, barcode reorder/subsets, negative controls, training-only preprocessing, missingness flags, geometric masks and hidden-label corruption.

The model-facing S3 interface is `CountTransform.observed_inputs(counts, observed_mask)`. It indexes observed rows **before** transformation, returns zeros at hidden rows plus an explicit observation flag, and remains unchanged when every hidden label is replaced by NaN. Full normalized target caches are for supervised training/evaluation only and must never be passed directly to an S3 forward call.

All mask files include ordered barcodes, a boolean `observed` vector, and physical distance to the nearest anchor. Geometry metadata records actual spot coverage, which can differ from nominal coverage because all boundary ties are retained. No platform-specific capture footprint or tissue-area fraction is claimed. S3 evaluation must use only `~observed` locations.

## Reports and next checkpoint

Read [C0_REPORT.md](reports/C0/C0_REPORT.md), [C0_SPEC_LOCK.json](reports/C0/C0_SPEC_LOCK.json), [C0_DATA_MANIFEST.json](reports/C0/C0_DATA_MANIFEST.json), and the figures in `reports/C0/alignment/`. The source commit and exact file hashes connect the data preparation run to its implementation. A subsequent commit stores generated results.

The plan requests review after every checkpoint. C1 will fit mean, ridge and small MLP H&E-only baselines only after C0 review. Graphs, gene-adaptive fusion, explanation claims, and external validation belong to later checkpoints. The current assay cannot supply the requested 1,000/3,000 biological gene panels; broad transcriptome experiments require another cohort.

## Attribution and data use

HEST: Jaume et al., *HEST-1k: A Dataset for Spatial Transcriptomics and Histology Image Analysis*, NeurIPS 2024. [Official implementation and terms](https://github.com/mahmoodlab/HEST).
The upstream HEST release specifies CC BY-NC-SA 4.0. Downloaded datasets and derived images retain applicable upstream terms; local raw/feature assets are excluded from Git. See the source metadata for original study links and sample-specific provenance. No license is assigned here to upstream model weights. ResNet18 weights are obtained from the official PyTorch distribution.
