# Inputs, models and acquisition

The separately authorized post-v1 research workflow is documented in [JOURNAL_REVISION.md](JOURNAL_REVISION.md). Its new fits/results and blocked raw-data cases are distinct from the original-state recovery/replay limitations below.

No dataset downloads or fitted-state replacements were performed. Acquisition is explicit; ordinary tests/plots use fixtures or existing CSV exports. Official source locations below are also identified in the supplied manuscript and archived acquisition records.

| Source | Official acquisition | Included roles and required files |
|---|---|---|
| HEST | https://huggingface.co/datasets/MahmoodLab/hest-bench and https://github.com/mahmoodlab/HEST | IDC NCBI783/NCBI785 fit, TENX95 validation, TENX99 development evaluation; H5AD counts/barcodes/genes/coordinates and cached RGB patch H5 |
| HER2ST | https://zenodo.org/records/4751624 | A1–H1, one section each; counts, barcodes/positions, tissue images; 40-gene intersection; reused calibration-development cohort |
| Wu local Visium | https://zenodo.org/records/4739739 | 1142243F/CID4465/CID4535 calibration; 1160920F/CID4290/CID44971 evaluation; feature/barcode/matrix, positions, images and scale factors |
| Bassiouni | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616 | First numbered section for BASS_P01–P09, P11–P15, P17–P24; source filtered feature-barcode H5, tissue positions, scalefactor JSON and image |

Do not use all sections from a study interchangeably. Exact selected original filenames, URLs, checksums and qualification settings remain in `checkpoints/legacy/C7/raw/sources.json`, C9 acquisition/download logs and C10/C10A/logs/SOURCE_SAMPLE_SELECTION.json and DOWNLOADS.json. These are reference records; historical acquisition scripts are not maintained commands and must not be run in place. Dataset and weight licenses must be checked at those sources; this revision grants no redistribution rights.

## Join and transformation contract

Join by barcode, never row position. Patch H5 coordinates are row/column; H5AD spatial coordinates are x/y. The reader rejects duplicate barcodes, wrong coordinate orientation and missing alignment. Harmonize gene identifiers, then retain exact ordered panels; no destination-expression panel selection. `data/raw/<cohort>/` and `data/cache/<cohort>/` are suggested configurable new locations, excluded from Git. Original caches have canonical missing paths in the manifest.

Targets are raw-count log1p followed by C0 development-training mean/population scale. No library-size or cross-assay correction is introduced. Preserve `legacy/reports/C0/C0_NORMALIZATION.json`, ordered panels, feature and head-target scalers. The selected head scalers are under `legacy/C4/C4A/preprocessing/resnet18_reference/primary/S2-50/` (inspect manifest for exact path). Source-to-Wu calibration is `legacy/C9/C9A/T1_SCALER.json`. Bassiouni reuses that calibration without fitting.

Nominal 112-um external crops are rounded at native pixels, resized with LANCZOS to 224, then officially resized to 256, center-cropped to 224, ImageNet normalized. Wu has only 15/25 native pixels per crop, Bassiouni 24–29; resizing adds no missing morphology. Visium joins use source in-tissue flags and the locked 100-um hexagonal-lattice image-coordinate registration.

## Required original primary heads

Install at `checkpoints/models/` followed by the original path below. No release URL was supplied. Obtain originals from the author; ImageNet ResNet18 weights are a separate backbone and cannot replace these learned heads.

| Original path | SHA256 |
|---|---|
| C4/C4A/models/C4A_primary_resnet18_reference_primary_S2-50_D1_h256d0.1_s0.pt | 4aaa11ab8998627dd4cac9a2b8da59e9fc8e43182e31ffee85f9204a0500b91e |
| C4/C4A/models/C4A_primary_resnet18_reference_primary_S2-50_D1_h256d0.1_s1.pt | 6c86019b731b67549b175d8b39a2ff4a3ced30d5a612ee2d8cdd4ee57f3e67e6 |
| C4/C4A/models/C4A_primary_resnet18_reference_primary_S2-50_D1_h256d0.1_s2.pt | 51ebf67a05e986586efe4564bc797cac6e31dd7352308b0fba7cf0ff0e9c9e68 |

All three are missing. This blocks `doctor --require-models` and actual `infer --config configs/replay.example.json`, as well as original image-prediction replay. The compact manifest also identifies supporting ridge, MLP, sparse residual and PCA state paths; none is installed. Feature caches and predictions are not learned states. Missing models cause a specific nonzero error with the expected checksum; no random fallback or automatic training is permitted.

The replay NPZ contract requires aligned `xy_um`, `embeddings` (N,512), `targets` (N,50), `barcodes`, ordered `genes`, boolean `observed` and `scored`. All arrays share one row order. Targets must already use unchanged development units. Use one section per file; V2 observed anchors must be outside the locked region, with scored rows equal to that region. The CLI verifies shapes/disjointness/gene order and original model/scaler hashes; the caller must supply the correct locked cohort/mask identity. This is a bounded embedding-based route, not an implemented raw-cohort ingestion pipeline. `max_locations` is capped at 10,000 to prevent a cohort job.

## Current publication metadata (not released fitted weights)

The PREFLIGHT_FIXED author sources and all numerical exports are now under `manuscript/revised/`. `methods_data/image_feature_scaler.csv` contains 512 ordered feature centers/scales; `target_and_calibration_constants.csv` contains the 50-gene development normalization, additional head-target transform and source/Wu moments. Publication regression compares these values with the original archived NPZ/JSON objects, at rtol=atol=1e-12. `primary_model_details.json` records architecture/optimization settings, and `repository_documentation/MODEL_FILES_TO_RELEASE.csv` records all three required original head identities. These metadata/scaler files are available; the original fitted .pt states remain unavailable. Metadata does not constitute a model release.

The canonical package is hash-covered in `manuscript/PACKAGE_PROVENANCE.json`. Builds and redraws copy it to a new output tree. The author package's replacement synthetic graph is for illustration only; it does not change empirical graph selection or fitted parameters. See `REPRODUCING.md` for the now-working publication/document/source-package workflows.
