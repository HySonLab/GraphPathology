# Codex handoff — GraphPathology publication synchronization

**Input:** this newly reconstructed full LaTeX ZIP. **Target:** public GitHub repository `HySonLab/GraphPathology`.

## Mandatory preflight

1. Extract the ZIP and verify its manifest. Confirm that `main.tex`, `article_only.tex`, `supplementary.tex`, `article_content.tex`, `supplementary_content.tex`, `figures/Figure2.pdf`, `scripts/draw_graph_architecture.py`, and `repository_documentation/REPRODUCTION_TARGETS.csv` exist. Reject a ZIP containing only the earlier 5-figure article.
2. Confirm title starts **GraphPathology enables geometry-aware...**, precisely 6 main figures, 2 main tables, 12 supplementary figures, 9 supplementary tables, and 29 BibTeX entries. The current combined PDF is 52 pages: this package is a reconstruction, not the previously claimed 54-page byte-identical source.
3. Read `README.md`, the current GitHub repository README, and the repository's package/tests before editing.

## Publication integration

1. Place or integrate all canonical LaTeX sources, vector figures, source tables, bibliography and figure-generation scripts. Avoid duplicate publication sources with competing titles.
2. Ensure the repository's publication build regenerates `main.pdf` (combined), `article_only.pdf` and `supplementary.pdf` in an isolated output directory. The main output must append Supplementary Information after references, with correct S-numbering.
3. Update `README.md` with concise quickstart and separate workflows: (a) rebuild paper/figures from frozen exports, (b) replay original predictions if fitted weights available, (c) retrain from raw input. Never present these as identical or assume all external assets have been released.
4. Use `repository_documentation/REPRODUCTION_TARGETS.csv` as the checklist for all 29 display items. Verify that figure/table paths exist and no narrative claims are inferred solely from a drawing script.
5. Protect Figure 2 during regeneration: `scripts/redraw_development_figures.py` now targets Figures 3 and 6; `scripts/redraw_external_figures.py` targets Figures 4 and 5; `scripts/draw_architecture.py` targets Figure 1 and `scripts/draw_graph_architecture.py` targets Figure 2. Regenerate into temporary output folders first. Verify the revised synthetic Figure 2 has 94/22/98/85 edges from its archived CSVs; the earlier 106/55/113/84 counts belong to a different synthetic example, not an empirical inconsistency.
6. Reconcile original C0..C10 stage paths and checkpoints with the maintained `checkpoints/` organization. Repair imports, CLI paths, doc references and manifests; preserve frozen original artifacts and hashes. Do not overwrite original result exports.
7. Update citation, paper link and code-availability documentation. Ensure README no longer says development stages are unexecuted.

## Bounded tests; do not rerun the scientific campaign

- Clean editable LaTeX build with no missing files, undefined citations, or references.
- Verify all 29 publication items are represented in the sources and numerical mapping.
- Run `python scripts/validate_publication_sources.py` as a bounded preflight before integration.
- Synthetic graph construction checks: deterministic seed, 36 coordinates, four graph rules, undirected edges, radius/cap consistency, invariant to node ordering.
- Small fixture CLI tests: install/import, sparse masking/hidden-label invariance, adjacency, model forward/backward, serialization, calibration, numeric interpolation, metric aggregation.
- Figure-generation smoke test on supplied exports in temporary output directory. Do not modify frozen exports while testing.
- Explicitly report the unavailable original fitted head binaries (if still absent) rather than claiming fully reproduced predictions.
- Provide final commit/PR summary, tests run, regenerated output paths and remaining unavailable assets.

**Do not rerun HEST/HER2ST/Wu/Bassiouni complete training or external evaluation; do not alter numerical results.**
