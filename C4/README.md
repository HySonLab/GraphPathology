# C4: frozen pathology representation benchmark

This checkpoint follows [the revised user plan](C4_EXECUTION_PLAN.md). It preserves C0–C3 and tests one accessible frozen pathology encoder, Phikon, against the exact cached ResNet18 representation. UNI2-h and UNI were inaccessible without gated credentials; the access audit records that limitation. No second encoder or biological interpretation is included.

The immutable [scientific specification](C4A/C4A_SPEC_LOCK.json) records the source hashes and decision rules before test extraction. [Protocol](C4A/configs/protocol.json) details matched tuning, train-only scaling, three MLP seeds, deterministic Ridge, two negative controls, and four conditional patient folds. The existing C0 target normalization remains frozen, so sensitivity folds are conditional rather than an untouched external validation.

## Execution

Prerequisites are the frozen C0–C3 data/artifacts and the recorded base environment. For a separate reproduction, prepare C4-only directories and packages with `.venv/bin/python C4/bootstrap.py --install-isolated-dependencies`; this never upgrades the base environment.

From the repository root, using the unchanged `.venv`:

```bash
.venv/bin/python C4/download_encoder.py
.venv/bin/python C4/extract.py trainval
.venv/bin/python C4/run.py select
.venv/bin/python C4/extract.py test
.venv/bin/python C4/run.py final
.venv/bin/python C4/analyze.py
.venv/bin/python C4/verify.py
.venv/bin/python C4/figures.py
.venv/bin/python C4/report.py
.venv/bin/python C4/freeze.py
.venv/bin/python C4/freeze.py --verify
```

These commands document execution, not permission to rerun frozen experiments. `select`, specification creation and final freeze refuse to overwrite locks. Completed fits and extraction caches are reused only after applicable hash checks. The isolated extraction packages under `C4/deps` are recorded in [the environment record](C4A/configs/execution_environment.json); the old environment is unchanged. Restore matching dependencies and source weights before a fresh reproduction in a separate output directory.

C4B is conditional on a clear C4A pass. In this execution the predeclared S2-50 gate passed and authorized B; the larger primary panels did not pass. The B specification and source were frozen before spatial training/evaluation. Its separate workflow is:

```bash
.venv/bin/python C4/b_lock.py
.venv/bin/python C4/b_spatial.py primary
.venv/bin/python C4/b_spatial.py final
.venv/bin/python C4/b_analyze.py
.venv/bin/python C4/b_verify.py
MPLCONFIGDIR=/tmp/c4-mpl .venv/bin/python C4/b_report.py
```

These also document the completed execution; do not rewrite existing locks or frozen outputs. Final `C4/freeze.py` runs only after both stages finish. Do not launch it for mixed, failed or technical-failure results. The final report and status file record whether this gate was reached. Stop after C4 for user review.

## Artifacts

- `C4_REPORT.md`, `C4_STATUS.json`: conclusions and next-checkpoint recommendation.
- `C4A/C4A_REPORT.md`, `C4A/*.csv`, `C4A/figures/`: detailed results and figures.
- `C4A/logs/`: run/evaluation records, controls, replay audit, timings and immutable registry events.
- `CROSS_MODAL_COMPATIBILITY.md`: official CMRCNet compatibility assessment; no approximate reproduction scores.
- `C4_MANIFEST.json`: hashes of local artifacts, model weights, embeddings, checkpoints and predictions, plus preserved global-registry prefix.

Large weights, embeddings, models and predictions remain local and are excluded from Git. Source, configuration, hashes and metric tables are committed. The final freeze verifies C0–C3 recursively without rerunning or tuning them.
