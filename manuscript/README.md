# Current publication sources

`revised/` is the author-supplied `Scientific_Reports_Revised_LaTeX_Package_PREFLIGHT_FIXED.zip`, preserved file-for-file. `PACKAGE_PROVENANCE.json` records the ZIP SHA256 and every supplied file SHA256. Current title: **GraphPathology enables geometry-aware reconstruction of spatial transcriptomics from sparse molecular measurements**.

This is explicitly a reconstructed 52-page combined source, not byte-identical source for the previous 54-page PDF. The current synthetic 36-location Figure 2 has 94/22/98/85 edges; the prior 106/55/113/84 illustration is superseded, not a changed empirical result. Historical PDFs and the earlier manual Table S8 schema are clearly separated in `history/`.

Use the maintained repository CLI from any working directory (pass `--root` for a non-editable installation). It copies sources and writes only to new output trees:

```bash
python -m graphpathology publication --output outputs/publication
python -m graphpathology build-documents --publication outputs/publication --output outputs/paper
python -m graphpathology arxiv --output outputs/arxiv
```

The author `build.sh` is preserved as source evidence. The maintained wrapper handles isolated paths, no-shell-escape builds, diagnostics, numbering checks and optional regenerated figures. Never run redraw scripts against canonical sources or the historical archive. The author files inside the ZIP are publication inputs, not independent authorization to run research campaigns.

Current review PDFs in this directory were generated from the numerically validated figures and shared sources. Their document structure is recorded in `../docs/DOCUMENT_BUILD_REPORT.json`. Compiled manuscript PDFs are intentionally excluded from the arXiv source ZIP.
