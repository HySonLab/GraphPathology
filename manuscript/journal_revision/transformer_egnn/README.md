# Scientific Reports editorial manuscript

This is the current author-review manuscript. Its scientific narrative centers on controlled measurement placement in 22 patients, with the complete original reconstruction results and negative graph comparisons retained. The abstract is 172 words and names GNNs only at a high level. Detailed development methods are consolidated in Supplementary Methods S6.

See `EDITORIAL_CHANGELOG.md`, `CLAIM_EVIDENCE_AUDIT.md` and `AUTHOR_RELEASE_CHECKLIST.md`. All 40 figures, 22 tables and 31 references remain. The [pre-edit manuscript snapshot](https://github.com/HySonLab/GraphPathology/tree/d1a0d1bb4e9d802b92186c010b1f86f3efdff2d6/manuscript/journal_revision/transformer_egnn) is preserved in Git history. Original arXiv sources, earlier journal files, experimental protocols, numerical exports and fitted states are unchanged.

Build the current source without research downloads or training:

```bash
python -m graphpathology.extension paper --build-output outputs/editorial-manuscript
python -m graphpathology.editorial check --output outputs/editorial-integrity.json
```

A checkout without author-local full gene files/states uses `check --committed-only`; skipped files are explicitly counted. For a clean standalone ZIP that is tested without BibTeX:

```bash
python -m graphpathology.editorial package --output outputs/GraphPathology_Scientific_Reports_Submission_Source.zip --check-directory work/source-check
```

The source ZIP contains `.bbl` files and all figure assets; `bash build.sh` typesets all three documents. Use `bash build.sh --bibtex` after editing references. `scripts/regenerate_extension_figures.py` redraws from existing summary tables, then applies the current readable route/cohort labels. No raw research inputs or checkpoints are needed for document/figure regeneration. Do not invoke the historical research assembly to overwrite this editorial source.

The full per-gene and model bundles remain author-local, with release-ready manifests and no claimed public deposit. Three original primary image-head binaries remain missing; publication generation is distinct from full original image-model replay. AI-assisted computational tools are disclosed in Methods, and the author's Editorial Board disclosure remains in Competing interests. PR #3 stays draft/unmerged for author review; no journal submission has occurred.
