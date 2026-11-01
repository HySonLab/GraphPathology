# Reproduction workflows

Use Python 3.12 and the installed versions in `configs/validation-environment.txt`. The historical executed environment remains unchanged at `checkpoints/legacy/configs/environment.txt`. The actual earlier host bootstrap, when ensurepip was missing, was `python3.12 -m venv work/venv`, then `python3.12 -m pip --python work/venv/bin/python install pip`. Normal systems can use the README venv commands.

Editable installs locate the repository automatically. Outside a checkout or with a non-editable installation, pass `python -m graphpathology --root /absolute/path/to/GraphPathology ...`. Explicit output paths are relative to the caller's working directory and must be new; canonical sources and archives are protected. Replay config paths resolve relative to the config file.

## A. Build existing reviewed figure assets

```bash
python -m graphpathology build-documents --dry-run
python -m graphpathology build-documents --output outputs/paper-reviewed-assets
```

Inputs: the canonical `manuscript/revised/` package, checked against all author-provided file hashes. Outputs: `main.pdf`, `article_only.pdf`, `supplementary.pdf`, isolated `project/` build tree and `BUILD_REPORT.json`. Combined output is 52 pages; article is 29; supplement is 23. Supplement starts after References on combined page 30, with continuous pages, 18 figure captions, 11 table captions and 29 bibliography entries. Undefined references/citations, duplicate PDF destinations and TeX warnings/underfull/overfull boxes cause failure. No-shell-escape is used.

Dependencies: TeX Live with pdfLaTeX, bibtex8 or BibTeX, Latin Modern, geometry, ragged2e, microtype, natbib, caption, placeins, lineno, hyperref and other packages in the supplied preamble. On Ubuntu the tested CI installation uses `texlive-latex-extra texlive-fonts-recommended lmodern`; Poppler provides pdfinfo/pdftotext/pdftoppm. No model/data downloads are involved.

## B. Regenerate and numerically validate all publication displays

```bash
python -m graphpathology publication --dry-run
python -m graphpathology publication --output outputs/publication
python -m graphpathology build-documents --publication outputs/publication --output outputs/paper
```

The atomic publication command regenerates all six main and twelve supplementary figures with the exact supplied author scripts. It operates only on new copies. All eleven table templates/digital sources are included and validated, rather than fabricated as new result tables. It produces `publication/`, a separately rendered unmodified-author `baseline/`, `REGRESSION.json` and 29 validated rows in `TARGET_STATUS.json`. Plot arrays, masks/NaNs, colors, color limits and gene/patient ordering are compared with explicit tolerances; every scientific export is byte-checked. The `--available-only` argument is retained for compatibility but there are no missing publication targets in this package. Single-target partial rendering is intentionally rejected; use the author functions in a scratch copy for figure-specific editing.

Dependency order: frozen CSV/gzip results and methods constants -> original selections/aggregation in author renderers -> numerical artist traces -> independent printed-cell, patient-median, pairing, graph and constant checks -> new figure assets -> document build. No fitted predictions are recomputed. Main Figure 4 uses `Fig5_source.csv`, Figure 5 uses `Fig6_source.csv`; legacy filenames are not display identifiers. Full SFig6_source has 4,448 rows. Figure 6's cap uses the entire map export, including undisplayed capture-window values. Figure 2 is a replacement synthetic illustration with 36 locations and 94/22/98/85 edges, not a change to patient graphs or research outcomes.

All 35 scientific export hashes remain unchanged after redraw. The complete compressed per-gene export is preserved; text definitions do not substitute for it. Printed rounding is checked exactly at the supplied decimal precision. Floating arrays/constants use rtol=atol=1e-12 with identical NaN/mask/order patterns. Synthetic coordinates/distances are exported to six decimals, so the independent distance check allows 2e-6 um; edge identity is exact. PDF-byte equality is not required across fonts/plotting versions.

## C. arXiv single-root source ZIP

```bash
python -m graphpathology arxiv --dry-run
python -m graphpathology arxiv --output outputs/arxiv
```

Output: `GraphPathology_arxiv_source.zip` plus `ARXIV_REPORT.json`. The dependency closure begins at main.tex and includes shared content, caption sources, figure PDFs, bibliography/style and main.bbl. It excludes extra document roots, compiled manuscript PDFs, logs and author notes. Numerical exports and methods metadata are ancillary files under `anc/` and are not TeX dependencies. A clean extraction is compiled using the included main.bbl without running BibTeX, then numbering and reference structure are checked.

## D. Bounded software checks

```bash
python -m pytest -q
python -m graphpathology smoke --output outputs/fixture
python -m graphpathology verify-archive --historical --report outputs/audit.json
```

Default tests are offline and CPU compatible after environment setup. The fixture uses 20 locations, 50 genes, synthetic local/graph heads and two backward steps. It exercises the development graph route and final external image/anchor/calibration/fusion route, including serialization. This is not original-state replay. Archive verification hashes 74,253 files and declared registry prefixes without executing historical analysis scripts. Full publication/TeX integration is a separate bounded CI job.

## E. Original fitted predictions and future retraining

```bash
python -m graphpathology doctor --require-models
python -m graphpathology infer --config configs/replay.example.json --dry-run
python -m graphpathology infer --config configs/replay.example.json
```

Currently blocked by missing original heads and aligned real section/cache inputs. The JSON is an example schema, not an included dataset. Required upstream order: acquire exact selected cohort/sections -> checksum verification -> barcode/coordinate/gene joins -> unchanged development normalization -> frozen crop/encoder -> original feature scaler -> three original heads -> inverse target scaling/equal average -> frozen source/Wu calibration -> locked observed anchors -> V1/V2 interpolation -> 0.8/0.2 fusion -> hidden-only metrics -> seed/geometry/patient summaries. Bassiouni is never refitted; final external inference has no graph head.

The maintained bounded route accepts one aligned section of up to 10,000 locations. It does not ingest a raw cohort. See `DATA_AND_MODELS.md` for sources, licenses, sample identities, arrays, original checksums and acquisition records. Full extraction, training and cohort evaluation require a separately implemented/validated orchestrator; historical campaign scripts must not be restarted in place. The tested `core.train_mlp` component is not a full cohort-training command. Retraining new heads cannot fill the missing-original-weight gap. No full research jobs are authorized or executed by these publication workflows.
