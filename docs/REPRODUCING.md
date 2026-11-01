# Reproduction workflows and dependency order

Use the maintained CLI from an editable installation. All scientific inputs are explicit, repository-root paths resolve through `--root`, and replay config paths resolve relative to the config file. From outside the checkout, use `python -m graphpathology --root /absolute/path/to/GraphPathology ...`. All outputs must use a new directory; archive/source/manuscript overwrites are rejected.

## Setup

Use the README install commands. In this host, system Python 3.12 had no ensurepip. The actual clean bootstrap was `python3.12 -m venv work/venv` (failed after creating the interpreter), then `python3.12 -m pip --python work/venv/bin/python install pip`, followed by normal installation using that new interpreter. No system environment was modified. Exact validated versions are in configs/validation-environment.txt; original versions remain in checkpoints/legacy/configs/environment.txt. TeX/Poppler are optional external build/view tools; no current TeX source was available to establish a complete TeX dependency list. Tests are offline after package installation and CPU compatible. No universal runtime/GPU/disk requirement is claimed.

## A. Build supplied documents from existing figures

`python -m graphpathology build-documents --dry-run`

Blocked: only the two revised PDFs were supplied. Revised TeX, bibliography, figure files and build.sh are missing both locally and in the public checkout. The command reports this with exit 2. Existing PDFs are preserved under manuscript/. Building from a PDF is not rebuilding its LaTeX source. Install the revised package under manuscript/revised/, inspect its script, then implement/validate the actual build order and label/citation checks. Do not run the old C10 verifier as the current publication validator.

## B. Regenerate figure/table views from numerical exports

```bash
python -m graphpathology publication --dry-run
python -m graphpathology publication --available-only --output outputs/publication
python -m graphpathology publication --target Figure4 --output outputs/wu
```

The complete all-target validation reports three missing targets and returns 2. The explicitly partial command generates 15 figure views and 11 table views, with 29-item TARGET_STATUS.json. It does not claim a complete revised publication build. Table S8 uses a new editable JSON schema transcribed from the supplied PDF and checked against production parameter-count formulas. Output folders hold PNG/SVG/PDF figure views, editable CSV/HTML tables, preserved complete compressed gene export, VALUE_TRACES.json and plotted-value JSON traces. Table previews are labeled first 100 rows; full CSVs are retained. Figures use saved numerical exports; no model replay or targets from plots.

Dependency order: original exported metrics / selected protocol -> patient summaries -> within-cohort plotting/table filtering -> value traces and numerical checks -> new figure/table views. Main Figure 4 uses old Fig5 values, Figure 5 uses old Fig6. Figure S1 uses old HER2ST figure 4; Figure S2 uses old synthesis figure 7; S3–S12 map to old S1–S10. S7 is explicitly restricted to the revised single random-mask/replicate-0 control condition. S5 retains both historical SFig3 inputs. Full C6 input rows are retained for current S8 (historical SFig6), including undisplayed rows. Main Figure 6 is blocked by the missing full map export; do not reconstruct error values from images or compute its color cap using only the two displayed conditions. Figure 2 needs exact 36-node coordinate literals for 106/55/113/84 edges.

No revised renderer/array baseline was supplied. Source row/schema checks and table-number comparisons are recorded, but exact revised layout/plotted-array equivalence is still unverified. See REPRODUCTION_TARGETS.csv and PUBLICATION_VALIDATION.json.

## C. Bounded software checks

```bash
python -m pytest -q
python -m graphpathology smoke --output outputs/fixture
python -m graphpathology verify-archive --historical --report outputs/audit.json
python -m graphpathology doctor
```

The fixture uses 20 locations, 50 genes, seeded synthetic local/graph heads and two backward steps. It tests development coordinates/embeddings -> aggregation -> head -> metrics, and external local image ensemble -> observed anchors -> calibration/interpolation -> fusion -> metrics. Save/load predictions are identical. This is software testing, not replay of fitted scientific results. Archive verification hashes files and recorded registry prefixes only, never runs historical replay scripts.

## D. Bounded original-state replay

```bash
python -m graphpathology doctor --require-models
python -m graphpathology infer --config configs/replay.example.json --dry-run
python -m graphpathology infer --config configs/replay.example.json
```

Currently blocked by original heads and aligned section/cache inputs. The example is a schema, not a working dataset. Follow DATA_AND_MODELS.md. Order: obtain original data -> verify checksums and selected patient/section -> barcode/coordinate/gene joins -> unchanged count normalization -> frozen image crop/encoder -> original feature scaling -> original three heads -> inverse head-target scaling and equal mean -> frozen source/Wu calibration -> locked observed anchors/V1 or V2 -> molecular interpolation -> 0.8/0.2 fusion -> hidden-only metrics -> seed/geometry/patient aggregation. No default missing-weight fallback exists. The final external route has no graph head. Full-cohort replay is not an implemented maintained command and was not executed.

## E. Future training/full evaluation

This revision does not restart the C0–C10 campaign. Historical scripts/configs record its dependency order, but their old generic imports/path assumptions and write side effects make them unsuitable as maintained reproduction commands. The maintained core.train_mlp is tested on tiny synthetic inputs, not a fully migrated cohort training orchestrator. A future orchestrator must explicitly configure new raw/cache/output roots, preserve fitting-only transforms and protocol locks, respect the order above, then add real-data integrations before claiming executable full research reproduction. Training new heads never fills the missing-original-state gap.
