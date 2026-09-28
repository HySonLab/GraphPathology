# C1: local morphology-only prediction

The frozen execution specification is `configs/C1_EXECUTION_PLAN_v2.md`; the finite implementation choices are in `configs/protocol.json`. C0 is read-only.

Run from the repository root, using the existing CPU environment:

```bash
.venv/bin/python scripts/verify_c0.py --include-data
.venv/bin/python -m pytest -q C1/tests
.venv/bin/python C1/run.py select
.venv/bin/python C1/run.py evaluate
.venv/bin/python C1/run.py lopo
```

Stages refuse duplicate run IDs. `select` loads only primary training and validation data, runs the predefined tuning grid, fits final primary models and negative controls, and writes `SELECTION_LOCK.json`. `evaluate` checks that fitting source has not changed and evaluates the frozen primary test patient. `lopo` uses fixed primary-selected parameters and fold-specific scalers, without evaluation-patient early stopping.

All four target settings are modeled separately. MLP configurations are selected using seed 0 and validation MSE; final primary MLPs use seeds 0, 1, 2. LOPO epochs are the rounded median primary best epoch per panel. Both negative controls use the already selected model configuration and permute entire rows, retaining gene-gene structure within each target vector.

C0 normalized targets are the reporting space. C1 feature and target scaling is fit on the appropriate training population and inverted for target metrics. Within each LOPO fold, restandardization of C0 targets algebraically removes the original affine shift/scale from model inputs; gene panels and primary-selected hyperparameters still make LOPO conditional. LOPO error units remain the frozen C0 units.

A constant predictor has undefined Pearson/Spearman and co-expression correlations. Such values are stored as blank CSV/JSON null values, never fabricated as zero. B0 comparison uses error and prediction variance. Threshold fractions use all panel genes as denominator; counts of defined correlations are reported separately.

`results/experiment_registry.csv` is append-only with separate fit/evaluation events. Candidate fits are included. Reported peak RSS is sampled process RSS (including resident input arrays), not isolated model allocation. CPU thread counts are fixed at four.

Model checkpoints and prediction arrays are retained locally under `C1/models` and `C1/predictions`, excluded from Git to avoid large binary history. Their checksums are recorded in run logs and the final specification. Source, configurations, metrics, figures, scalers, logs and reports are versioned. Reproduction on a different machine requires reconstructing the verified C0 caches first and running in a fresh separately versioned output location rather than overwriting this freeze.
