# C2: isolated spatial-context comparison

C0 and C1 are read-only. C1 is neither retrained nor replayed. Baselines are read from its frozen per-gene metrics, summaries and conditional-LOPO results. C1 primary/fold preprocessing objects are reused without fitting new scaling.

The exact finite design is in `configs/protocol.json`. MEAN is a self-inclusive spatial mean followed by the C1-selected small MLP. GCN is one symmetric-normalized graph convolution followed by the same activation, dropout and direct readout. GraphSAGE uses a single concatenated self/neighbor-mean convolution, with width chosen to match C1 parameter count within 1%. Because embeddings are frozen and aggregation is linear, aggregating before the learned linear map is mathematically equivalent and saves CPU work. These are deliberately shallow models; one-layer SAGE is also expressible as an MLP on concatenated local and neighborhood features.

Every graph is built within one slide using physical coordinates only. There are four predefined geometries, with no multiscale fusion. Search is validation-only. Three seeds assess main results; a vertex-permuted graph control tests proximity specificity. No coordinates or neighboring expression are supplied to a learned predictor beyond the graph operator. No graph may connect patients or slides.

Run the stages from the repository root:

```bash
.venv/bin/python C1/freeze.py --verify  # read-only integrity check, no C1 inference
.venv/bin/python -m pytest -q C2/tests
.venv/bin/python C2/run.py select
.venv/bin/python C2/run.py evaluate
.venv/bin/python C2/run.py lopo
```

The selection stage writes its lock before test access. Validation per-gene C1 predictions were not stored, so the validation gene shortlist uses graph-versus-MEAN gains only; no new C1 inference is used to invent that missing baseline. Test and conditional-LOPO contextual gains are compared directly to the existing C1 local B2 results by patient/panel/seed. Graph sensitivity is reported on validation, not used to retune against test data.

Conditional LOPO retains primary-selected graphs/configurations and fixed primary-derived epoch counts. Its targets and hyperparameters were not selected independently per fold; it is patient-dependence sensitivity, not unbiased or external validation.

Large models and predictions are stored locally and hashed, but excluded from Git. Code, scalers referenced by immutable hashes, graph files, metrics, figures, logs and reports are versioned. Run IDs are unique and the global registry is append-only. Stop after C2; C3 requires review.

Final reporting and read-only integrity verification:

```bash
.venv/bin/python C2/summarize.py
.venv/bin/python C2/graph_figures.py
.venv/bin/python C2/verify.py
.venv/bin/python C2/build_report.py
.venv/bin/python C2/freeze.py
.venv/bin/python C2/freeze.py --verify
```

After freezing, use only the final `--verify` command to audit the checkpoint; do not regenerate frozen outputs. See [C2_REPORT.md](C2_REPORT.md) for the completed `MIXED_CONTEXT_SIGNAL` decision, limitations, and review gate before C3.
