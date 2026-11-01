# Archived stages and fitted model states

`legacy/` preserves every tracked file from starting commit `03506defbe205cabd13865780d6a1224b348a781`, including the earlier README/environment and original shared source. Original C0 is `legacy/reports/C0`, not a root C0 directory. All experimental stages are now under this directory.

`models/` is the canonical installation home for original fitted states, including supporting heads, ridge/PCA objects. None were present initially. Keep each historical subpath below models/ to avoid basename collisions; retain filenames. Feature caches, masks, predictions and scalers remain distinct asset types.

`manifest.csv` inventories all present files and known missing references. `path_map.csv` records the one-to-one move with pre-move SHA256 and byte count. `historical_expectations.csv` preserves declared digests and registry-prefix lengths separately from the migration snapshot. Missing does not mean downloadable; obtain original states from the author unless an actual release URL is documented. Multiple historical expected digests are retained in the expectations file.

```bash
python -m graphpathology verify-archive --historical --report outputs/audit.json
python -m graphpathology resolve C9/C9A/T1_SCALER.json
python -m graphpathology doctor --require-models
```

Verification is read-only. No historical lock/hash is rewritten. Do not execute historical verify.py, freeze.py, run.py or verify_package.py as ordinary tests: these may run large analyses, append registries or write audits. Maintained code imports no archived runtime modules. See ../docs/MIGRATION.md.
