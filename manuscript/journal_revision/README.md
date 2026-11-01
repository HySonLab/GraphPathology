# Scientific Reports journal revision — post-arXiv v1

This is a separate working journal manuscript based on executed new analyses. The original arXiv v1 package remains byte-preserved in `../revised/`, and original root manuscript PDFs remain unchanged. New comparisons reuse previously inspected cohorts; no new blind validation or journal acceptance is asserted.

- [Combined manuscript + SI](GraphPathology_Combined_Manuscript_Journal_Revision.pdf)
- [Article only](GraphPathology_Article_Only_Journal_Revision.pdf)
- [Supplementary information](GraphPathology_Supplementary_Journal_Revision.pdf)
- [Executed protocols, status and commands](../../docs/JOURNAL_REVISION.md)
- [New display-target ledger](JOURNAL_TARGETS.csv)

The package retains the original 29 display targets and adds nine separately sourced journal targets: two main figures, one main table, two real-data supplementary maps and four supplementary tables. Original v1 metadata notes are prefixed `V1_`; they describe inherited displays, not validation of the new targets. The new build report is `JOURNAL_BUILD_REPORT.json`.

Wu matched-information, operator-factorial and identical-target analyses are completed. Original-export patient uncertainty includes all 22 Bassiouni and all three Wu patients. Thirty new local/mean/GCN/GraphSAGE development heads and all control/seed evaluations are complete. Their new weights never substitute for the missing original states. New Bassiouni predictions and BASS_P01 maps remain BLOCKED because official transfers returned errors/empty bodies; no predictions were fabricated. Optional extensions explicitly skipped in the preregistered protocol are not presented as completed experiments.

`source_data/journal/` contains actual new numeric exports; raw full tissue images, complete raw count matrices and upstream encoder weights are not redistributed. Real Wu maps and their derived numeric sources are attributed to Wu et al., Zenodo 4739739, CC-BY4.0. Inherited v1 source exports remain unchanged. New model training used pinned HEST inputs under their CC-BY-NC-SA4.0 research terms; full provenance is under `experiments/journal_revision/`.

For an offline source build, run `bash build.sh` in this directory. For strict numbering/warning checks and regeneration from committed journal exports, run from the checkout root:

```bash
python -m graphpathology.journal paper --build-output outputs/new_journal_build
```

This builds three PDFs in an isolated staging copy and updates only journal outputs. It neither downloads research inputs nor trains. TeX Live, BibTeX and Poppler are required. The source `main.tex` places SI after completed References; all numbering and citations are checked. Source tables/plots retain all declared estimators and seeds, including negative and constant-field examples. Raw Spearman near-tie sensitivity and float32 self-control precision are explicitly diagnosed; they do not weaken the original molecular replay tolerance or exact hidden-label checks.

Verified PDFs: combined 63 pages, article 35 pages, SI 28 pages; 22 figures, 16 tables and 29 references, with no TeX diagnostics, unresolved marks or page-number gaps. Updated pages were visually reviewed.
