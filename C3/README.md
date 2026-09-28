# C3: spatial specificity, then gated scale identification

The revised execution plan is preserved in `C3_EXECUTION_PLAN_UPDATED.md`. C0–C2 remain frozen. C3A tests topology specificity using immutable A0/A2/A3 references, a parameter-matched self-only control, and five graph randomizations crossed with three training seeds for each required randomized control. C3B may be created only after `PASS_SPATIAL_SPECIFICITY`; a later multiresolution architecture is outside C3.

`C3A/C3A_SPEC_LOCK.json` fixes configuration, graph audits, descriptive gene criteria, and the gate before new test evaluation. `C3A/ANALYSIS_IMPLEMENTATION_LOCK.json` separately records the pre-test gate implementation. The final `C3_MANIFEST.json` hashes outputs without rewriting those pre-test records.

Execution from the repository root (do not regenerate a frozen checkpoint):

```bash
.venv/bin/python C2/freeze.py --verify
.venv/bin/python -m pytest -q C3/C3A/tests
.venv/bin/python C3/C3A/lock.py
.venv/bin/python C3/C3A/run.py prepare
.venv/bin/python C3/C3A/run.py train_primary
.venv/bin/python C3/C3A/run.py final_evaluation
.venv/bin/python C3/C3A/analyze.py
.venv/bin/python C3/C3A/verify.py
.venv/bin/python C3/C3A/figures.py
.venv/bin/python C3/C3A/report.py
.venv/bin/python C3/freeze.py
```

After freezing, audit without rewriting any result:

```bash
.venv/bin/python C3/freeze.py --verify
```

The runner can resume completed stages by loading their recorded artifacts. Graph preparation and model execution should run sequentially. Every new graph/fit/evaluation has a lossless JSON registry event and an append-only event in `results/experiment_registry.csv`. The original CSV header is frozen by C1/C2; new graph fields are supplied by `C3A/C3A_REGISTRY_EXTENSIONS.csv`, keyed by event and run ID. Multi-slide fits carry per-slide graph references.

Checkpoints and prediction arrays remain locally in `C3A/models` and `C3A/predictions`, are included in final hashes, and are excluded from Git. Reproduction requires the same frozen C0 cache and local C1/C2 artifacts. No old model is refitted or replayed.

See [C3_REPORT.md](C3_REPORT.md) and [C3_STATUS.json](C3_STATUS.json) for the final gate decision. The full shared assay panel contains 280 genes.
