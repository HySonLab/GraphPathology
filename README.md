# GraphPathology

**GraphPathology enables geometry-aware reconstruction of spatial transcriptomics from sparse molecular measurements**
Truong-Son Hy, University of Alabama at Birmingham.

**Latest manuscript:** [Combined paper + SI](manuscript/journal_revision/transformer_egnn/GraphPathology_Combined_Manuscript_Transformer_EGNN.pdf) · [Article only](manuscript/journal_revision/transformer_egnn/GraphPathology_Article_Only_Transformer_EGNN.pdf) · [Supplement](manuscript/journal_revision/transformer_egnn/GraphPathology_Supplementary_Transformer_EGNN.pdf) · [Editable sources and build guide](manuscript/journal_revision/transformer_egnn/README.md)

GraphPathology investigates how spatial measurement geometry, histology representations and graph learning affect reconstruction of spatial transcriptomics from sparse molecular observations. Controlled experiments in 22 patients emphasize the importance of where molecular measurements are collected, alongside how predictions are computed.

The study compares image-only prediction, molecular interpolation and image–molecular fusion. Matched comparisons evaluate graph and local predictors with the same encoder and molecular support. GCN-assisted reconstruction improves over local-assisted prediction across the 22-patient cohort with both tested encoders, while anchor-only inverse-distance interpolation remains stronger on cohort-median endpoints. Mixed and negative findings, neighborhood averaging and shuffled controls qualify the interpretation of graph topology. The tested fixed image contribution did not reproducibly improve external molecular interpolation.

Distributed and contiguous observation settings differ in weighting and scoring support as well as geometry. Previously inspected cohorts are not fresh blind validation. A potential fivefold reduction concerns selective acquisition and proportional variable location costs, rather than demonstrated total assay savings.

## Data and study design

| Cohort | Role | Panel |
|---|---|---|
| HEST NCBI783 / NCBI785 | Fitting | 50 ordered genes |
| HEST TENX95 / TENX99 | Validation / primary development evaluation | 50 |
| HER2ST A1–H1 | Transportability, then calibration development on the same patients | 40-gene intersection |
| Wu 1142243F, CID4465, CID4535 | Calibration | 50 |
| Wu 1160920F, CID4290, CID44971 | Initially held-out evaluation; subsequently inspected | 50 |
| Bassiouni, first section for each of 22 patients | Initially held-out evaluation without refitting; subsequently inspected | 50 |

See the [experimental workflow](docs/JOURNAL_REVISION.md), [matched-backbone comparison](docs/MATCHED_BACKBONE_COMPARISON.md) and [graph architecture experiments](docs/GRAPH_TRANSFORMER_EGNN.md) for protocols, baselines, controls, selection procedures and limitations.

## Install and CPU quick start

Tested with Python 3.12, torch 2.8.0+cpu and torchvision 0.23.0+cpu. [Validation versions](configs/validation-environment.txt) are distinct from the original research environment.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r configs/validation-environment.txt
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m graphpathology smoke --output outputs/quickstart
.venv/bin/python -m pytest -q
```

The smoke route uses tiny synthetic heads and two backward steps to verify the software. Each output directory must be new. Imports, help and tests do not download datasets or model weights. If venv creation lacks ensurepip, see [REPRODUCING.md](docs/REPRODUCING.md).

## Build the latest manuscript

The maintained [manuscript sources and build guide](manuscript/journal_revision/transformer_egnn/README.md) include the article, supplementary information, vector figures and tables. Build all three documents from committed sources without downloading research data or running experiments:

```bash
python -m graphpathology.backbone_paper build --output outputs/current-manuscript
```

The output directory must be new. TeX Live, BibTeX and Poppler are required. Builds check figure and table numbering, references, continuous page numbering and TeX diagnostics.

## Models and reproducibility

Verified redistributable model heads, configurations, checksums and supporting numerical artifacts are available in the [research artifacts release](https://github.com/HySonLab/GraphPathology/releases/tag/research-artifacts-v1). The [release and reproduction guide](docs/RELEASES.md) distinguishes recovered checkpoints, newly retrained replacements and artifacts withheld from publication, with verification and restoration instructions. Phikon model weights remain private pending clarification of licensing and redistribution permissions.

The three original selected primary head checkpoints remain missing; retrained replacements have separate identities and never silently substitute for them. To inspect the original checkpoint requirements and bounded inference configuration:

```bash
python -m graphpathology doctor --require-models
python -m graphpathology infer --config configs/replay.example.json --dry-run
```

These commands report required paths and checksums. The maintained inference route starts from aligned embeddings and normalized targets. The [experimental workflow](docs/JOURNAL_REVISION.md) documents public-data acquisition, frozen feature extraction, training and evaluation. Reproduction reports retain numerical differences, negative findings and incomplete checks transparently.

## Repository organization

- `src/graphpathology/`, `tests/`, `configs/`: numerical package, CLI, tests and configurations.
- `manuscript/journal_revision/transformer_egnn/`: latest manuscript sources, figures and PDFs.
- `manuscript/revised/`, `manuscript/history/`, `docs/history/`: preserved source packages, earlier documents and validation records; [package provenance](manuscript/PACKAGE_PROVENANCE.json).
- `reproduction/`: recovery inventories, reproduction reports, verification records and release manifests.
- `checkpoints/`: model inventories, installation locations and preserved analysis archives.
- `outputs/`: ignored builds, regression traces and fixtures.

Archived analysis integrity can be checked without running archived research scripts:

```bash
python -m graphpathology verify-archive --historical --report outputs/archive-audit.json
```

## Citation and terms

Use the manuscript title and author above; [CITATION.cff](CITATION.cff) supplies citation metadata. HEST, HER2ST, Wu, Bassiouni and upstream model weights retain their own terms. See the [release guide](docs/RELEASES.md) for artifact-specific distribution conditions. Raw tissue data and upstream pretrained encoder weights are not redistributed in this repository. The repository has no root software LICENSE; inclusion here does not establish a software license grant.
