# C7 external validation

Final scientific outcome: **MIXED_EXTERNAL_VALIDATION**, technical integrity **PASS**, recommendation **NEED_METHOD_REVIEW**. See [C7_IMMEDIATE_SUMMARY.txt](C7_IMMEDIATE_SUMMARY.txt) and [C7_REPORT.md](C7_REPORT.md). Stop for review; this repository checkpoint does not authorize adaptation or a later experiment.

## Scope and provenance

C0–C6 are immutable dependencies. The active cohort is HER2ST, eight independent source patients, one section each, A1–H1. [Zenodo source](https://zenodo.org/records/4751624), [study](https://doi.org/10.1038/s41467-021-26271-2), [official source code](https://github.com/almaan/her2st). Public archive passwords embedded in the reader are published in that source README; they are not private credentials. Source archive MD5 and local SHA-256 records are retained. Respect source data terms when reusing images.

GSE311609 and Xenium Prime were screened and excluded before evaluating outcomes. Their image/metadata acquisition and exploratory registration scripts are retained for provenance, not alternate evaluation pipelines. Target acquisition for those excluded candidates is explicitly disabled. Exploratory registration scores are not certified registration accuracy. No candidate molecular values informed qualification or model settings.

## Environment

Commands were run from `/home/son/Research/GraphPathology` using `.venv/bin/python` (Python 3.13, existing PyTorch/scientific environment), four CPU threads and deterministic inference. Exact package versions/platform are in `C7A/logs/ENVIRONMENT.json`. Additional acquisition/inspection dependencies are isolated in ignored `C7/deps`; no base environment upgrade was performed. Only pyzipper and pycryptodomex are needed from this set for the admitted cohort's encrypted archives; OpenSlide/OpenCV/imagecodecs were for excluded-candidate registration.

```bash
.venv/bin/python -m pip install --target C7/deps --no-deps pyzipper pycryptodomex
```

For figure rendering use `MPLCONFIGDIR=/tmp/c7-mpl XDG_CACHE_HOME=/tmp/c7-cache`.

## Read-only verification of the completed checkpoint

```bash
.venv/bin/python C7/freeze.py --verify
```

This verifies the complete C0–C6 chain, the C7 spec, all C7 manifest artifacts, and the append-only registry prefix. It requires the local large artifacts; a Git-only clone is intentionally insufficient for full prediction verification. Do not rerun generators in the frozen checkpoint: audit logs, scientific sources and outputs are themselves hashed.

`verify.py` implements the independent scientific audit (352 mask replays, every prediction hash/identity, 620 stratified prediction/metric replays and corruption invariance). Its completed result is `C7B/logs/INTEGRITY_AUDIT.json`; invoking it again writes that log, so run it only in an isolated audit copy if necessary. `presentation.py` also independently reconstructs all 8,556 saved distances/bin occupancies. The pre-outcome contract tests are in `tests/test_contract.py` and their passed output is retained in C7A logs.

## Historical execution sequence

The order below documents the execution; it is not a command to overwrite the frozen checkpoint.

1. Read public cohort metadata and exclude development/related samples. `acquire_her2.py` reconstructs the three source archive downloads from captured URLs and verifies MD5 without parsing expression values.
2. `her2_metadata.py`: gene headers/barcode identities only, supplied coordinates, physical scale, H&E extraction. `prepare.py`: deterministic image crops and draft protocol. `extract.py`: exact frozen ResNet18 and C4 image-only heads; no external labels. `qualification_figures.py` and `qualify.py`: image/coordinate QA, manifest and qualification report.
3. Finalize the recorded protocol before outcomes, including removal of its draft authorization field; `lock.py` hashes all scientific sources, exact protocol, frozen models, scalers, raw assets and external image caches. Spec lock is the authoritative evaluation authorization. The draft protocol generator must not replace the checked-in finalized `configs/protocol.json`.
4. Commit qualification/spec (`c06e92d`). `policy.py` freezes 352 trajectories with two selection-corruption checks each. Commit masks (`455b9a9`) **before** expression access.
5. `execute.py`: integer/finite count checks, exact C0 normalization, 8,448 main and 108 control runs, observed-only inputs, all-run hidden-label corruption checks, append-only registry. It is resumable by run ID but not a way to retune this checkpoint.
6. `analyze.py`: frozen patient-level summaries and gate. `verify.py`: independent replays and parent checks. `comparison.py`: re-summarize existing C5 gene metrics on identical retained genes, no parent model rerun.
7. `presentation.py`: descriptive tables, paired gene gains, distance verification and 13 figures. `domain_diagnostics.py`: descriptive scale/sparsity summaries after evaluation, no calibration. `report.py`: reports/status/60-line immediate summary. These presentation additions do not alter the locked decision rule.
8. Commit results and source, record results commit, freeze final hashes with `freeze.py`, commit manifest, verify and push. No C7C, adaptation or later checkpoint.

A prospective re-execution must use a separate copy/output root and preserve these locks; do not delete locks or change flags to rerun in place. Bitwise cache-file identity may require the archived original NPZ/HDF5 bytes even when numerical reconstruction is identical.

## Artifact distribution

Git includes source/configs, source metadata, qualification diagnostics, all mask NPZs and barcode traces, per-run JSON metrics/provenance, appended registry rows, compressed gene metrics, descriptive CSVs, PNG figures and reports. `C7B/C7B_GENE_METRICS.csv.gz` decompresses exactly to the required local `C7B/C7B_GENE_METRICS.csv`.

Large raw assets (`raw/` binary files), image/target/embedding caches (`cache/`), and prediction arrays (`C7B/predictions/`) remain local and are hashed by `C7_MANIFEST.json`. Dataset download URLs/hashes, frozen model paths, and reconstruction code are provided; a public repository alone does not include those large artifacts. `C7_MANIFEST.json` records the results/source commit and all local hashes; it necessarily excludes itself and its digest.

Per-run JSON files record sample/patient, panel, model contract hash, image/cache/target hashes, mask checksum, observed/hidden counts, seeds, metrics, runtime and peak RSS. NPZs retain predictions, exact observed/hidden barcodes, gene metrics/variance and nearest-anchor distance. All unobserved expression is removed before constructing input payloads, including non-scored unobserved tissue outside V2 regions.

The summary's Git commit identifies the completed results snapshot. The later manifest commit seals the checkpoint; neither hash can self-reference the file that creates it.
