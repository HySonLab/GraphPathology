# C5: sparse spatial-transcriptomics reconstruction

See [C5_REPORT.md](C5_REPORT.md) for results and [C5_STATUS.json](C5_STATUS.json) for the checkpoint decision. C0–C4 remain immutable. No C6 experiment is part of this checkpoint.

The predictor receives H&E features and coordinates everywhere and measured expression only at observed locations. All genes at each hidden location are withheld together. Results distinguish random anchors, contiguous hidden regions, and contiguous observed windows; measured coverage is the realized fraction of locations.

## Execution and reproducibility

Run from the repository root with its existing `.venv`. The original workflow was:

```bash
.venv/bin/python C5/prepare.py
.venv/bin/python C5/lock.py
.venv/bin/python -m pytest C5/tests
.venv/bin/python C5/run.py select
.venv/bin/python C5/run.py final
.venv/bin/python C5/analyze.py
.venv/bin/python C5/verify.py
MPLCONFIGDIR=/tmp/c5-mpl XDG_CACHE_HOME=/tmp/c5-cache .venv/bin/python C5/figures.py
.venv/bin/python C5/report.py
.venv/bin/python C5/freeze.py
```

The specification was frozen before validation selection, and selection was locked before primary hidden-label evaluation. Do not rerun mutating stages in the frozen output directory. The read-only integrity check is:

```bash
.venv/bin/python C5/freeze.py --verify
```

`common.py`, `engine.py`, `run.py`, `analyze.py`, and the protocol/mask configurations are checked against the pretest source lock. `verify.py` checks every prediction hash and barcode partition, training-only image scalers, training-anchor means, checkpoint provenance, all saved hidden-label invariance results, and an independent stratified prediction/metric replay. It also verifies the complete frozen parent chain. Model selection uses validation MSE; primary results never choose models, thresholds, map examples, or mask geometries.

## Artifacts

- `C5_PANEL_METRICS.csv` and `C5_GENE_METRICS.csv.gz`: hidden-only correlation and error metrics, retaining each mask and patient.
- `C5_COVERAGE_CURVES.csv`, `C5_DISTANCE_TO_ANCHOR.csv`, and `C5_GENE_COVERAGE.csv.gz`: coverage and distance behavior without forcing a fidelity horizon or minimum gene coverage.
- `C5_SETTING_GAINS.csv`, `C5_GENE_GAINS.csv.gz`, and `C5_PATIENT_SENSITIVITY.csv`: paired B2–B1, B3–B2, and B4–B3 comparisons.
- `C5_NEGATIVE_CONTROLS.csv`, `C5_RESOURCE_PROFILE.csv`, `logs/`, and `C5_REGISTRY_EXTENSIONS.csv`: controls, resource use and lossless provenance alongside the append-only global registry.
- `C5_SPEC_LOCK.json`, `SELECTION_LOCK.json`, and `C5_MANIFEST.json`: pretest decisions, validation choices and final artifact hashes.

Large prediction/model/cache files and uncompressed gene tables remain local and are listed in the final manifest. Their gzip metric counterparts are committed. Reproduction requires the existing frozen C0–C4 data/model artifacts and Python environment. This four-patient conditional sensitivity analysis is not external validation.
