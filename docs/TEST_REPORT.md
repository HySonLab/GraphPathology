# Actual validation report

Validation date: 2026-10-31 (America/Chicago). Starting revision: 03506defbe205cabd13865780d6a1224b348a781. Dedicated branch: codex/publication-reproducibility. Python 3.12; torch 2.8.0+cpu and torchvision 0.23.0+cpu, no downloaded encoder weights. Exact installed versions: configs/validation-environment.txt. Historical execution versions are preserved separately at checkpoints/legacy/configs/environment.txt. Tested versions do not establish bitwise original-model replay across environments.

## Executed

| Command / check | Result | Measured duration |
|---|---|---|
| Clean Python environment + editable install of .[test], CPU torch/torchvision | Installed; no system package modifications | Not separately measured |
| `python -m pytest -q` | **84 passed, 0 failed, 0 skipped** | **5.39 s**, after relocation to the delivered repository |
| `python -m graphpathology smoke --output 'outputs/quick start with spaces'` | 20 locations,50 genes,two synthetic training steps; both routes; serialization equality | Not separately measured |
| `python -m graphpathology verify-archive --historical` | 74,253 pre-move SHA256 identities unchanged; no migration damage; available historical digests match;10 registry-prefix scopes retained | 13.15 s |
| `python -m graphpathology publication --available-only --output outputs/publication` | 26 export views generated;three targets explicitly blocked;29-item report | Not separately measured |
| `python scripts/validate_export_views.py --output outputs/publication` | 66 checks passed;printed Table2 cells,external curves/paired effects,source bytes/rows/nulls,traces | Not separately measured |
| `python -m graphpathology publication --dry-run` | Expected exit2: missing revised inputs | Not separately measured |
| `python -m graphpathology build-documents --dry-run` | Expected exit2: revised LaTeX/build source absent | Not separately measured |
| `python -m graphpathology doctor --require-models` | Expected exit1: original fitted heads absent | Not separately measured |
| `git push --dry-run origin HEAD:refs/heads/codex/publication-reproducibility` | Failed authentication: could not read Username, terminal prompts disabled | Not separately measured |

Default tests cover numerical/data contracts, fixture hidden-label NaN/random corruption (including unscored V2 hidden rows),nested trajectories,graph formulas/ties/locality,parameter counts,frozen encoders,dropout,three-head inverse scaling,calibration moments/direction,gate population variance,34-column residual/intercept convention,CSV control typing,metrics/aggregation,missing-state errors,output safety,and CLI help from outside the checkout. No test collection imports archived runtime chains. No research registry is appended.

## Historical test disposition

All original test files remain byte-identical in checkpoints/legacy. Root C0 tests are retained in the supported suite. C1 numerical tests are promoted to test_core;C2 to test_graph;C3 rewiring to test_topology_controls;C4 CSV typing to test_csvio;C5 to test_residual;C6 to test_sampling;C9 to test_external. C3/C4 config/capacity and C7 label-free sampling checks are covered in test_additional_contracts/test_workflow. C8 patient access/selection/gate features are tested on explicit fixtures;C10 no-new-calibration is checked against its saved protocol and C9 scaler.

The original C7 real-scene feature/metadata checks,C8 real-scene calibration test,and C10 original predictor comparison cannot run without ignored research caches/states. Their files/assertions are preserved;they are not default passes or silently weakened/skipped tests. The maintained synthetic routes test software behavior,not equality to unavailable original patient predictions. Full legacy suite execution remains unvalidated and is not claimed.

## Publication limits and visual checks

15 PNG/SVG/PDF figure export views and11 editable table export views generated. Contact-sheet visual review found readable axes and no clipped data;the revised documents are not rebuilt or exactly layout-matched. The synthesis view was changed to five rows/two columns to match the current SI organization;plot values remain traced. Table previews are labeled;complete CSV/gzip data remain retained. Only selected numerical comparisons against the supplied manuscript could be checked without its revised source package. Full caption/citation/label resolution,36-node schematic equivalence,FOXA1 full-value color limit,and both LaTeX builds are blocked. See PUBLICATION_VALIDATION.json and PUBLICATION_STATUS.json.

No cohort training,feature-extraction campaign,calibration search,full-cohort inference,large evaluation or historical verify.py replay was run. No random/zero parameters were installed as original models. All source data/model identities were preserved. GitHub Actions workflow is added but has not run remotely.
