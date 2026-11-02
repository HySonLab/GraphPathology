# Executed and repeatable commands

Run in a checkout with the registered inputs acquired and hash-verified as documented in the original data guides. Restore missing public head/result assets first. Phikon distribution remains blocked; the complete author-local integrity command requires those recovered files.

```bash
python -m pytest -q
python -m graphpathology.reproduction_primary train \
  --config configs/reproduction/primary-v1.json \
  --checkpoints checkpoints/reproduced/primary_resnet18/primary-20261101-v1 \
  --output reproduction/runs/primary-20261101-v1
python -m graphpathology.reproduction_replay tier-a \
  --output reproduction/runs/tier-a-v1
python -m graphpathology.reproduction_replay tier-b \
  --checkpoints checkpoints/reproduced/primary_resnet18/primary-20261101-v1 \
  --output reproduction/runs/tier-b-v2
python -m graphpathology.reproduction_artifacts verify \
  --output reproduction/runs/tier-c-v3
python -m graphpathology.reproduction_artifacts matched \
  --output reproduction/runs/matched-v2
python -m graphpathology.reproduction_reporting
python -m graphpathology publication --output work/reproduction/publication
python -m graphpathology.reproduction_release preserve \
  --output work/reproduction/your-fresh-preservation-report.json
python -m graphpathology.reproduction_release build \
  --output /absolute/path/to/a-new-release-staging-directory
```

Fresh run/output identities are mandatory: do not rerun training into an existing folder or overwrite a completed output. Numerical rerun commands use the same fixed design; the three-seed fitting budget was exhausted after the single completed training run. Set OPENBLAS_NUM_THREADS=4, OMP_NUM_THREADS=4 for computation; recorded per-model Torch thread overrides are read from immutable reference metadata. Tests ran with one BLAS/OpenMP thread. Exact hardware, packages and stage durations appear in ENVIRONMENT.json and COMPUTE_SUMMARY.json. The trained primary code is pinned in its own manifest at commit10b6827f2; later interface validation does not alter its fitting function. Final packaging points to the tested release code commit.

The original doctor --require-models intentionally exits1; ORIGINAL_IDENTITY_CHECK.json verifies three absent original identities. No PDF build command was run. Publication figure regression regenerated separate scratch outputs and left every frozen file unchanged. Initial failed assertions and the four-thread sensitivity run remain documented in EXECUTION_ISSUES.md; they did not trigger additional fits or changed thresholds.
