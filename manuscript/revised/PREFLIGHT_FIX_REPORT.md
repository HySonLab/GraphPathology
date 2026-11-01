# Publication preflight corrections (2026-11-01)

This ZIP is a revision of the reconstructed 52-page GraphPathology package. It is **not** the missing earlier 54-page archive.

## Findings addressed

1. **Figure output mappings repaired:**
   - `scripts/draw_architecture.py` → Figure 1.
   - `scripts/draw_graph_architecture.py` → Figure 2.
   - `scripts/redraw_development_figures.py` → Figures 3 and 6.
   - `scripts/redraw_external_figures.py` → Figures 4 and 5.
   - Supplement redraws → Figures S1–S12.
   Redraw smoke tests used separate temporary output directories; none modified manuscript figures or scientific data.
2. **Synthetic graph counts clarified:** Figure 2 has 36 synthetic locations; the correct retained undirected-edge counts for the current script are 94 / 22 / 98 / 85 (kNN 6 capped 250 µm / radius 150 µm / radius 250 µm / Delaunay capped 250 µm). These values are validated against the Figure 2 synthetic edge CSV. The earlier 106 / 55 / 113 / 84 values came from a *different* synthetic demonstration, not the real tissue data. No empirical result changed.
3. **Underfull box warning resolved:** the TESLA-related-work sentence was rephrased for line breaking, with no change to its scientific meaning. The refreshed LaTeX builds have no `Underfull \hbox`, `Overfull \hbox`, or `LaTeX Warning` entries.

## Automated validation

`python scripts/validate_publication_sources.py` passed: **35 exact export SHA-256 comparisons**, 29 expected publication targets, six main figures, twelve supplemental figures, 36 synthetic locations, four graph edge lists, edge-cap checks, and old-script numbering protection.

`bash build.sh` succeeded: 52-page combined PDF (`main.pdf`), 29-page article (`article_only.pdf`) and 23-page supplement (`supplementary.pdf`). All references/citations resolved; 29 bibliography items.

## Caveats

The package is still a reconstructed publication source, not byte-identical to the unavailable earlier 54-page version. It has not rerun research experiments and does not include the original external image model head binaries.
