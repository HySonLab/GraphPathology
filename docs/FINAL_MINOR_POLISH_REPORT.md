# Final minor editorial polish

Status: COMPLETE. Editorial only; no experiments, model fitting, data/metric changes or manuscript restructuring.

## Exact manuscript source changes

All paths below are relative to `manuscript/journal_revision/transformer_egnn/`:

- `article_content.tex`: natural Phikon encoder/preprocessing introduction; removed the unrelated ResNet18 diagram/cohort sentence from Figure 6; replaced the existing future-acquisition sentence with one prospective direction limited to the predefined expression panel.
- `supplementary_content.tex`: S6 caption distinguishes the original frozen ResNet18 evaluation from the later exploratory Phikon evaluation in Wu/Bassiouni; the S6 plot remains a one-patient development comparison. Changed “final sparse pipeline” to “original sparse pipeline” in the historical representation-method description.
- `figures/Figure1.tex`: caption title identifies the original study/architecture.

Documentation changes: `README.md` preserves the exact model identifier/revision in reproducibility prose; `EDITORIAL_CHANGELOG.md` records this small revision. This report is new. Generated changes are the three canonical PDFs and their `EXTENSION_BUILD_REPORT.json`; figure assets, bibliography, table contents, numbering map, code and protocols are unchanged.

## Preservation and checks

- Abstract is byte-identical to the prior 178-word version. SHA-256: `dc66615cc33b432baf989ba39257d81e2145d097b0b93bedbeab03d879d54981`.
- All numerical tokens in the three edited TeX sources are unchanged after accounting for the explicitly removed immutable version hash and misplaced ResNet18 diagram reference. All remaining Phikon preprocessing/extraction details are byte-identical. The pinned `owkin/phikon` revision remains `057cc0295895c2df3dd7681a89680da6015cbefe` in unchanged configuration and provenance records.
- Positive GCN-versus-local findings, molecular-reference findings, limitations, author-directed AI disclosure and Editorial Board statement are unchanged.
- 142 offline tests pass. Original publication regression passes 29 display targets and 71 numerical checks.
- 8,887 protected prior scientific files are unchanged; the 74,253-file historical audit has zero changed/missing files. All 10 protected numerical table bodies, 72 fitted states and 936 full per-gene result arrays pass existing integrity checks.
- Combined/article/supplement builds have 94/42/52 pages, 40 figures, 24 tables and 40 cited references, with zero final TeX diagnostics and correct continuous/S numbering.
- Rendered and visually inspected affected combined pages 7,14,27,29–30,51 and 89, plus nearby page/float transitions 6,8,13,15,28,50 and 52. No clipping, caption/label errors or stale external-evaluation statement remains.
- The new author-review source ZIP also compiles all three documents in a clean extraction without BibTeX, datasets or checkpoints. Earlier scientific bundles and their manifests remain unchanged.

Baseline commit: `5696498aaf1588878b8d067c7e056cad1d8b5825`. Updated existing draft PR #3; it must remain open and unmerged for author review. No release or journal/arXiv submission was performed. Live GitHub checks provide the CI status.
