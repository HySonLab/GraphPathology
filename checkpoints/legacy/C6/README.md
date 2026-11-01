# C6: fixed-budget spatial transcriptomics sampling

Start with [C6_REPORT.md](C6_REPORT.md), [C6_STATUS.json](C6_STATUS.json), and the detailed [C6A report](C6A/C6A_REPORT.md). C0–C5 remain frozen. The revised plan requires an audited C6A pass before adaptive C6B; external validation is not part of this checkpoint.

## Reproduction

Use the existing repository `.venv` and the local frozen C0–C5 artifacts. The original sequence was:

```bash
.venv/bin/python -m pytest C6/tests -q
.venv/bin/python C6/prepare.py
.venv/bin/python C6/execute.py select
.venv/bin/python C6/execute.py trajectories
.venv/bin/python C6/execute.py final
.venv/bin/python C6/analyze.py
.venv/bin/python C6/verify.py
MPLCONFIGDIR=/tmp/c6-mpl XDG_CACHE_HOME=/tmp/c6-cache .venv/bin/python C6/figures.py
.venv/bin/python C6/report.py
```

The source/spec lock precedes validation selection. The full trajectory lock precedes every C6 test-expression reveal. Models, preprocessing, gene panels and target units are inherited from C5; no reconstruction weights are refitted. C5 validation MSE determines the primary reconstruction family; a molecular-only interpolation baseline is also retained. Policies compare identical anchor counts. Open sampling and a fixed unavailable-region stress test remain separate.

Do not rerun mutating stages in a frozen checkpoint directory. After finalization, the read-only artifact-chain check is:

```bash
.venv/bin/python C6/freeze.py --verify
```

`support.py` provides the observed-only molecular boundary; `policy.py` implements exact deterministic FPS and random trajectories; `execute.py` performs validation selection, trajectory freezing and reconstruction; `analyze.py` contains the pretest decision rules. `verify.py` independently replays every test policy and stratified predictions and verifies the parent chain. `figures.py` and `report.py` generate presentation artifacts from saved outcomes.

## Results and traces

C6A retains per-run panel and gene metrics, the complete random distribution, signed sampling gains, AUC-coverage, geometry diagnostics, matched-distance comparisons, negative controls, runtime and resource records. `policies/*.npz` saves ordered selected indices/barcodes, acquisition scores, eligible indices and the unavailable-region mask. A budget selects the exact prefix specified by the frozen count rule. Every evaluation records its one-shot observed-label reveal prefix and hidden-only predictions.

Large prediction arrays and full uncompressed gene/selection CSVs remain local and are included in the final hash manifest; gzip copies of the tables are committed. Saved traces, configurations, source, reports and figures are committed. Four-patient conditional development is not unbiased external validation, and no monetary savings are inferred.
