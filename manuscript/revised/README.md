# GraphPathology — corrected reconstructed publication source

**Paper:** GraphPathology enables geometry-aware reconstruction of spatial transcriptomics from sparse molecular measurements

This is a newly **reconstructed source package**, not a byte-identical recovery of the previously discussed 54-page archive. It was rebuilt from the available earlier 110-file LaTeX project plus the author-requested structural and scientific edits. Its current verified combined PDF has **52 pages**. Codex should integrate these actual sources, not rely on the earlier page-count statement or older archive metadata. Original exported numerical results were retained; no research experiments were rerun.

## Required contents

- `main.tex` — **combined article and supplementary information** (shared content via `article_content.tex` and `supplementary_content.tex`). The supplement begins after references.
- `article_only.tex` — main article only.
- `supplementary.tex` — standalone Supplementary Information.
- `article_content.tex` and `supplementary_content.tex` — canonical editable text sources shared by the above wrappers.
- `references.bib`, `main.bbl`, `article_only.bbl` — 29-entry bibliography and generated files.
- `figures/` — six vector main figures and six caption sources; Figure 2 is the geometry-aware graph illustration (synthetic points, four topology rules).
- `supplementary_figures/` — 12 supplementary vector figures.
- `tables/` — three table sources; Table 3 is **included as Supplementary Table S9**, not as a third main display item.
- `source_data/` and `methods_data/` — archived numerical exports and model-method settings; `Figure2_synthetic_*.csv` are newly generated geometry illustrations only.
- `scripts/draw_graph_architecture.py` — regenerate the new Figure 2 only (supports `--output-dir` and `--data-dir` for isolated testing). `scripts/redraw_development_figures.py` regenerates Figures **3 and 6**; `scripts/redraw_external_figures.py` regenerates Figures **4 and 5**; `scripts/draw_architecture.py` regenerates Figure **1**. All development and external scripts now have corrected publication figure numbers. Prefer passing `--output-dir` to a temporary directory while testing.
- `repository_documentation/REPRODUCTION_TARGETS.csv` — 29-item target map for publication regeneration.
- `Codex_Publication_Sync_Instructions.md` — specific integration and bounded testing requirements.

## Publication preflight and figure provenance

Run `python scripts/validate_publication_sources.py` to verify the archived export checksums, the required figure/table assets, publication figure-number mappings, and the four synthetic graph edge counts. Run redraw scripts with `--output-dir` pointed at a temporary directory to avoid changing reviewed figure files.

**Figure 2 provenance:** the current 36-location, 6×6 jittered *synthetic* example has **94 / 22 / 98 / 85** retained undirected edges (symmetrized kNN / radius 150 µm / radius 250 µm / capped Delaunay). The earlier **106 / 55 / 113 / 84** counts came from a different illustration. They are not empirical model results and are not expected to match; the current figure, caption, edge CSV and drawing code agree.

## Build

Install LaTeX (`pdflatex`, `bibtex8` or `bibtex`) plus the included `naturemag.bst` bibliography style. Run:

```bash
bash build.sh
```

Outputs: `main.pdf` (combined), `article_only.pdf`, `supplementary.pdf`, and `GraphPathology_Combined_Manuscript.pdf` (copy of combined).

For arXiv, use **`main.tex`** as the compilation entrypoint. No experiments, model weights, raw images or network are needed to typeset the paper from this source.

## Important limitations

This source was reconstructed after the original later archive was not present in the available files. It reflects the latest *requested* framing and incorporates the available numerical outputs, but it must be reviewed by the author against any previously saved 54-page PDF. GraphPathology's graph-based morphology contextualization is a development study; the external replication benefits arise from spatially weighted molecular anchors. The theoretical fivefold reduction concerns *measured locations/variable location-dependent costs*, not validated total assay expenditure. This package does not contain original fitted image-head weights; see `repository_documentation/MODEL_FILES_TO_RELEASE.csv`.
