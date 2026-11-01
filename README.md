# GraphPathology

**GraphPathology enables geometry-aware reconstruction of spatial transcriptomics from sparse molecular measurements**
Truong-Son Hy, University of Alabama at Birmingham.

[arXiv v1 combined paper + SI](manuscript/GraphPathology_Combined_Manuscript.pdf) · [Article only](manuscript/GraphPathology_Article_Only.pdf) · [Supplement](manuscript/GraphPathology_Supplementary.pdf)

GraphPathology separates contextual image representations, reconstruction from observed same-section molecular anchors, and image–molecular fusion. Graph neural predictors are development comparisons; the external headline results come from molecular interpolation. The tested fixed image contribution did not reproducibly improve external interpolation. V1/V2 differ in weighting and scoring support as well as geometry. External images are downsampled. A potential fivefold reduction concerns selective acquisition and proportional variable location costs, not proven total assay savings or learned long-range correlations.

| Cohort | Role | Panel |
|---|---|---|
| HEST NCBI783 / NCBI785 | Fitting | 50 ordered genes |
| HEST TENX95 / TENX99 | Validation / primary development evaluation | 50 |
| HER2ST A1–H1 | Transportability, then calibration development on the same patients | 40-gene intersection |
| Wu 1142243F, CID4465, CID4535 | Calibration | 50 |
| Wu 1160920F, CID4290, CID44971 | Held-out at v1; now inspected | 50 |
| Bassiouni, first section for each of 22 patients | Held-out at v1, no refitting; now inspected | 50 |

## Post-v1 Scientific Reports journal revision

[Journal experiments and reproduction](docs/JOURNAL_REVISION.md) · [Separate journal sources/PDFs](manuscript/journal_revision/README.md) · [Preregistered protocol](experiments/journal_revision/PROTOCOL_V1.md). New real-data matched-information and identical-target analyses, patient uncertainty and controlled neural training are separate from the arXiv v1 package. Wu comparisons and original 22-patient Bassiouni uncertainty are completed; new Bassiouni raw-data comparisons are blocked by source transfer failure. See the journal status for precise completed/blocked/skipped work. Previously inspected cohorts are not a new blind validation.

## Install and CPU quick start

Tested with Python 3.12, torch 2.8.0+cpu and torchvision 0.23.0+cpu. [Validation versions](configs/validation-environment.txt) are distinct from the original research environment.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r configs/validation-environment.txt
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m graphpathology smoke --output outputs/quickstart
.venv/bin/python -m pytest -q
```

The smoke route uses tiny synthetic heads and two backward steps; it is a software fixture, not original fitted-prediction replay. Each output directory must be new. Imports/help/tests do not download datasets or original weights. If venv creation lacks ensurepip, see the tested bootstrap in [REPRODUCING.md](docs/REPRODUCING.md).

## Paper and figures from frozen exports

The canonical author-supplied **PREFLIGHT_FIXED reconstructed** sources live in [manuscript/revised](manuscript/revised/README.md). They produce a 52-page combined manuscript, a 29-page article and a 23-page supplement. This is not a byte-identical recovery of the earlier 54-page source. All **29 display targets** (18 figures and 11 tables) are mapped in [REPRODUCTION_TARGETS.csv](docs/REPRODUCTION_TARGETS.csv).

```bash
python -m graphpathology publication --dry-run
python -m graphpathology publication --output outputs/publication
python -m graphpathology build-documents --publication outputs/publication --output outputs/paper
python -m graphpathology arxiv --output outputs/arxiv
```

Regeneration runs the exact author renderers in isolated copies, compares every figure's numerical artist arrays/masks/colors/limits/order against the unmodified author baseline, and checks printed table cells, definitions, graphs and constants. It preserves scientific export bytes, including the full gene CSV.gz. Builds reject unresolved citations/references, duplicate destinations and TeX warnings. The combined document places the supplement after References with continuous pages and correct S numbering. The arXiv command creates a single-root source ZIP and compiles a clean extraction, with numeric ancillary files under `anc/`. No original model binaries, raw tissue data or research reruns are needed for these workflows.

Install a TeX distribution (`pdflatex`, `bibtex8`/`bibtex`, Latin Modern, latex-extra packages) and Poppler (`pdfinfo`, `pdftotext`, `pdftoppm`) for document builds and inspection. Plotting uses the Python environment above. [Actual results and tolerances](docs/TEST_REPORT.md), [per-target status](docs/PUBLICATION_STATUS.json), [numerical regression](docs/PUBLICATION_VALIDATION.json).

## Original predictions and retraining are separate workflows

```bash
python -m graphpathology doctor --require-models
python -m graphpathology infer --config configs/replay.example.json --dry-run
```

The original fitted states are still missing, including all three selected primary heads. These commands report the required paths and checksums and never fall back to new/random weights. The maintained inference route starts from aligned embeddings and normalized targets; it is not a complete raw-cohort research orchestrator. The separate explicit [journal research workflow](docs/JOURNAL_REVISION.md) now supports public-data acquisition, frozen feature extraction, controlled new training and evaluation. New journal fits have independent identities and cannot substitute for replay of original fits. The original bounded inference route retains its existing contract.

## Organization and provenance

- `src/graphpathology/`, `tests/`, `configs/`: maintained numerical package, CLI and bounded tests.
- `manuscript/revised/`: current immutable author sources, vector figures, exports and metadata. [ZIP/file provenance](manuscript/PACKAGE_PROVENANCE.json).
- `manuscript/history/`, `docs/history/`: explicitly superseded earlier PDFs and partial-validation reports.
- `checkpoints/legacy/`: byte-preserved C1–C10, `reports/C0`, shared source/configuration and historical registry.
- `checkpoints/models/`: canonical installation home for original states; [inventory](checkpoints/manifest.csv) records unavailable assets honestly.
- `outputs/`: new ignored builds, regression traces and fixtures.

```bash
python -m graphpathology verify-archive --historical --report outputs/archive-audit.json
```

This checks archived bytes and declared registry prefixes without running historical replay scripts. [Methods map](docs/MANUSCRIPT_MAP.csv), [migration map](checkpoints/path_map.csv), [migration notes](docs/MIGRATION.md), [remaining limitations](docs/REMAINING_WORK.md).

## Citation and terms

Use the manuscript title and author above; [CITATION.cff](CITATION.cff) supplies citation metadata. No journal acceptance or DOI is asserted. HEST, HER2ST, Wu, Bassiouni and upstream model weights retain their own terms. The original repository had no root software LICENSE; this revision does not invent a license grant or redistribute raw tissue data or upstream weights.

### Post-v1 bounded Transformer / EGNN extension

The separate [extension documentation](docs/GRAPH_TRANSFORMER_EGNN.md) and [scientific report](experiments/journal_revision/transformer_egnn/SCIENTIFIC_REPORT.md) describe 48 new controlled fits and 78-model evaluation in 26 previously inspected patients. Exact GEO range recovery completed the 22-patient matched spatial baseline study. The primary neural models did not improve same-information IDW; all negative findings, controls and patient effects are retained. The versioned journal manuscript is under [manuscript/journal_revision/transformer_egnn](manuscript/journal_revision/transformer_egnn/). Earlier arXiv, journal and A–E evidence remains unchanged. Full gene exports/checkpoints have hashed author-local bundles and a release-ready manifest; raw datasets/model binaries are excluded from Git. Draft PR #3 requires author review before merge.
