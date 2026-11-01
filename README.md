# GraphPathology

**Coverage and tissue geometry govern the reconstruction of spatial transcriptomics from sparse molecular measurements**  
Truong-Son Hy, University of Alabama at Birmingham.

This study asks how measured-location coverage and missing-region protocols affect reconstruction of a fixed expression panel. Image-only prediction, same-section molecular interpolation and image–molecular fusion receive different information. Distributed molecular anchors help reconstruction; contiguous regional completion is weaker. The tested fixed histology contribution did not reproducibly improve external interpolation. External public images are downsampled, and V1/V2 differ in weighting and scoring support as well as geometry. These results do not establish universal GNN superiority, transcriptome-wide recovery or clinical utility.

| Cohort | Role | Panel |
|---|---|---|
| HEST NCBI783 / NCBI785 | Fitting | 50 ordered genes |
| HEST TENX95 / TENX99 | Validation / primary development evaluation | 50 |
| HER2ST A1–H1 | Transportability, subsequently calibration development on the same patients | 40-gene intersection |
| Wu 1142243F, CID4465, CID4535 | Calibration | 50 |
| Wu 1160920F, CID4290, CID44971 | Fresh evaluation | 50 |
| Bassiouni, first section for each of 22 patients | Fresh replication, no refitting | 50 |

## Install and CPU quick start

Tested on Python 3.12 with CPU torch 2.8.0 / torchvision 0.23.0. See [exact validation versions](configs/validation-environment.txt); the original execution record remains in the archive.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m pip install torch==2.8.0 torchvision==0.23.0 --index-url https://download.pytorch.org/whl/cpu
.venv/bin/python -m graphpathology smoke --output outputs/quickstart
.venv/bin/python -m pytest -q
```

If your OS lacks ensurepip, install its Python venv package, or use the tested bootstrap described in [REPRODUCING](docs/REPRODUCING.md). The smoke command runs two synthetic CPU routes and two head-training steps; it does not reproduce the original fitted predictions. Each output path must be new. Imports/help/tests do not fetch datasets or model weights. Basic asset inspection uses `python -m graphpathology doctor`.

## Repository layout

- `src/graphpathology/`: maintained numerical contracts and explicit CLI.
- `checkpoints/legacy/`: original C1–C10, `reports/C0`, shared source/configuration, append-only registry and old manuscript, preserving relative paths and bytes.
- `checkpoints/models/`: installation home for checksum-verified original fitted states. **No original fitted states were present in the checkout.**
- `checkpoints/manifest.csv`: present and referenced missing assets, identities and canonical paths; [migration map](checkpoints/path_map.csv).
- `manuscript/`: supplied revised PDFs and revision brief. The revised LaTeX/source export package was not supplied or found on GitHub.
- `docs/`, `configs/`, `tests/`: reproduction guide, portable example and bounded regression suite.
- `outputs/`: new ignored test/reproduction outputs.

## Publication and reproducibility status

The current article has **6 main figures, 12 supplementary figures, 2 main tables and 9 supplementary tables**. All 29 are listed in [REPRODUCTION_TARGETS.csv](docs/REPRODUCTION_TARGETS.csv). Legacy display numbers are not the current numbers: old `Fig5_source` supplies current Figure 4.

```bash
python -m graphpathology publication --dry-run
python -m graphpathology publication --available-only --output outputs/publication
python -m graphpathology build-documents --dry-run
python -m graphpathology verify-archive --historical --report outputs/archive-audit.json
```

The first validation returns nonzero for missing revised sources. `--available-only` generated 26 labeled export views using archived numerical tables and maintained plotting functions. Figures 1/2/6 need missing revised schematic/map inputs. **Neither revised LaTeX document could be built, and exact revised-display equivalence has not been established.** Original PDFs remain available for reading. No numerical values were estimated from plots. The old package verifier/audit remains historical evidence, not the validator for this article.

[Full dependency order and workflows A–E](docs/REPRODUCING.md), [data/model acquisition](docs/DATA_AND_MODELS.md), [methods map](docs/MANUSCRIPT_MAP.csv), [migration](docs/MIGRATION.md), [test report](docs/TEST_REPORT.md), [remaining work](docs/REMAINING_WORK.md).

## Data use, citation and licensing

Please cite the manuscript title and author above and cite the source studies listed in the manuscript. No DOI, journal acceptance or publication year is invented here. HEST, HER2ST, Wu and Bassiouni data retain their own terms; consult their official records before redistribution. The starting repository contains no root software LICENSE, so this revision does not invent a license grant. Upstream model weights retain their terms. Raw datasets, upstream weights and new scientific outcomes are not included as part of this revision.
