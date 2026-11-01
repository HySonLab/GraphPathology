# Publication integration handoff

Base commit: `96efca22d383dd37732bc62bc2c2b7842d92d915`. Branch: `codex/final-publication-sync`.

Integrated the author-supplied PREFLIGHT_FIXED package under `manuscript/revised/`, with its exact ZIP/file hashes in `manuscript/PACKAGE_PROVENANCE.json`. The package is explicitly a reconstructed 52-page source, not the missing earlier 54-page source. Earlier PDFs and partial reports are retained in clearly marked history directories; all original checkpoint/archive bytes remain unchanged.

The maintained publication command now runs exact author renderers in isolated output copies and compares all 18 numerical figure traces with an independently rendered unmodified-author baseline. Table sources and all 35 scientific export hashes are preserved. Independent checks verify printed cells, gene/patient order, cohort aggregation, paired effects, full-map color limits, graph topology/permutation invariance and original preprocessing/calibration constants. All 29 display targets have explicit validated statuses; no missing targets are represented as passes.

Replaced the former always-raising build gate with isolated pdfLaTeX/BibTeX builds of combined main.pdf, article_only.pdf and supplementary.pdf. Added final-log, citation/reference, S-numbering and continuous-page checks. Added a tested one-root arXiv source command, including main.bbl, all input/figure dependencies and ancillary numeric files under anc/. A clean extraction compiles successfully without BibTeX. New builds/refreshed review PDFs use the current manuscript title. Source, migration and citation documentation have been updated.

Results: **98 tests passed**; **71 publication regression checks passed**; **29 targets validated**; document pages **52/29/23**; arXiv clean-extraction build **52 pages**; all **74,253** historical file hashes unchanged. Current detailed reports are `TEST_REPORT.md`, `PUBLICATION_VALIDATION.json`, `PUBLICATION_STATUS.json`, `DOCUMENT_BUILD_REPORT.json` and `ARXIV_REPORT.json`.

Run the README commands from an editable install. Full validated output trees are under local `outputs/publication-sync-final`, `outputs/publication-sync-documents` and `outputs/publication-sync-arxiv`; generated regression/work trees are not committed. Review PDFs are tracked under manuscript/. Current author figures/exports are tracked as input sources, not new scientific outcomes.

Original fitted states remain unavailable. Their exact identities and affected doctor/infer commands are documented in DATA_AND_MODELS.md and the checkpoint manifest. Publication reproduction does not require them. Full raw-cohort orchestration and original prediction replay remain explicitly incomplete and are not claimed by this publication-focused revision. No research campaign was restarted, and no original weights or numerical results were replaced.

This PR is for review, not automatic merging. Its GitHub checks distinguish bounded CPU validation from the separate publication/document workflow.
