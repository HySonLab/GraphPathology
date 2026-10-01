# C8 cross-platform calibration

Outcome: **MIXED_CROSS_PLATFORM_CALIBRATION**. C8A passed target harmonization; C8B was mixed, so C8C did not run. Recommendation: **NEED_METHOD_REVIEW**. Start with [C8_IMMEDIATE_SUMMARY.txt](C8_IMMEDIATE_SUMMARY.txt) or [C8_REPORT.md](C8_REPORT.md).

HER2ST is development transportability data in C8. Neither this cross-validation nor the previously examined C7 cohort is pristine final external validation. No later experiment is authorized by this checkpoint.

## Read-only verification

From `/home/son/Research/GraphPathology`:

```bash
.venv/bin/python C8/freeze.py --verify
```

This checks the full frozen C0–C7 chain, the C8 scientific spec, all C8 artifacts, and the append-only registry prefix. Full verification requires the local prediction arrays and inner calibration caches; Git intentionally excludes these larger arrays while retaining their hashes. It does not write to frozen outputs.

`verify.py` implements the scientific replay: all prediction identities and C7 T0 matches, explicit held-out-label exclusion while rebuilding every scaler, inner selection replay, exact F3 refits, stratified prediction/metric replay, and corruption invariance. Its completed result is `logs/INTEGRITY_AUDIT.json`. Running it again overwrites that audit log, so use an isolated audit copy rather than the frozen checkout.

## Scientific implementation

- `contract.py`: patient-scoped target access; equal-patient source/destination moments; observed-only inputs; frozen interpolation/image output reuse; common-unit errors; small fusion families and calibration objective.
- `run.py`: inner fitting, outer scaler/model locks, held-out predictions, corruption checks, negative controls, append-only registry. Imports parent modules read-only. Never use a generic `import run` from an audit script: parent directories contain files with the same module name; use an explicit file import as in `verify.py`.
- `analyze.py`: frozen A/B/C scientific gates and patient summaries. Its source was hashed before fitting/evaluation and never modified afterward.
- `presentation.py`: descriptive comparisons, all-gene effects, negative controls and nine figures. No model selection.
- `report.py`: the completed mixed-result reports/status/61-line summary; `C8_RESULTS_COMMIT.txt` identifies the completed results snapshot.
- `lock.py`: original prospective scientific lock. `freeze.py`: final immutable artifact manifest.

T1/T2 harmonization uses separately fitted source and destination scales. The source is the original C0 training patients only; destination statistics are fit on allowed HER2ST training patients. Frozen source image outputs are standardized with the source map, while observed anchors/targets use the destination map. Prediction errors are evaluated both natively and after inversion to original C7 units. Only common-unit errors enter selection and pass criteria; numerical shrinkage of the target scale cannot satisfy the gate.

All C8A inner and outer scaler objects are saved under `C8A/scalers`. All C8B families and selected models are under `C8B/calibration`. Each prediction NPZ saves native/common predictions, gene metrics, exact observed/hidden barcode lists and distances. Per-run JSONs retain scaler/model/mask/prediction hashes and the two hidden-expression corruption checks. The `runtime_seconds` field measures the audit/metric/serialization part of `save_outer`, not complete end-to-end inference time; it is not used to make speed claims.

## Historical execution order

This documents the completed run. Do not invoke generators in place on the frozen checkpoint.

1. Verify C7 and read the revised plan. Fix the 40-gene panel, four patient-disjoint folds, T0/T1/T2, original masks, fusion families/objective, quantitative gates and control definitions in `configs/protocol.json`.
2. Run the synthetic contract tests, then `lock.py`; commit scientific source/spec (`b75b764`). Tests used synthetic target values; no C8 outcomes informed the lock.
3. `run.py fit_A`: each outer training set has six patients; inner leave-one-patient-out scalers have five. All four selected transforms and 84 scalers are locked before outer A access. Commit (`5a97385`).
4. `run.py evaluate_A`, then `analyze.py A`. 4,320 outer A runs, all three transforms and E0/E1/E2. The frozen A gate passes; write C8A report.
5. `run.py fit_B`: use the selected-transform inner records only; 48 candidates per fold, including F0. Save parameters, inner scoring tables, fit hashes and B spec lock; commit (`d29ed78`).
6. `run.py evaluate_B`, then `analyze.py B`. 2,400 clean outer B evaluations plus 576 controls. The frozen B gate is mixed, so no `evaluate_C` or `analyze.py C` execution and no C8C directory.
7. `verify.py`, `presentation.py`, `report.py`; inspect figures, verify compressed/plain gene-table equivalence, commit results, record results SHA, freeze final manifest, commit it, verify and push.

The data and source dependencies already exist under frozen C0–C7. No new downloads or packages were required. Exact Python/package versions and fixed thread settings are in `logs/ENVIRONMENT.json`. Figure commands used `MPLCONFIGDIR=/tmp/c8-mpl XDG_CACHE_HOME=/tmp/c8-cache`.

## Stored outputs

The required uncompressed `C8B/C8B_GENE_METRICS.csv` exists locally; its byte-equivalent gzip is committed. Large `*/predictions/` and `C8A/inner_cache/` arrays are local and hashed. Code, configs, all scaler/model objects, masks by immutable C7 reference, per-run metadata, metric CSVs, figures, reports and registry events are committed. The final manifest includes local arrays and records its source commit; it cannot include its own hash recursively.

Only C8 and new registry rows are changed. C7 reference predictions, source weights, original normalizers, panels, masks and all parent artifacts remain frozen. The repository does not claim that a Git-only clone contains every source dataset or cached prediction array.
