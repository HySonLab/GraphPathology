# Publication synchronization validation

Validated on 2026-10-31 in the local host, using Python 3.12 and the installed versions in `configs/validation-environment.txt`; graphpathology editable version 0.3.0. Base: merged commit `96efca22d383dd37732bc62bc2c2b7842d92d915`. Branch: `codex/final-publication-sync`. The original research versions remain separate, unchanged archive records.

| Executed command/check | Actual result | Duration |
|---|---|---|
| Editable install, `pip install --no-deps --no-build-isolation -e .` in the existing isolated validation venv | graphpathology 0.3.0 installed | Not separately timed |
| `python -m pytest -q` | **98 passed, 0 failed, 0 skipped; no warnings** | **5.71 s** |
| `python -m graphpathology publication --output outputs/publication-sync-final` | **29 targets validated; 71 regression checks** | **42.07 s** |
| `python -m graphpathology build-documents --publication outputs/publication-sync --output outputs/publication-sync-documents` | Combined/article/supplement **52/29/23 pages**, no final TeX diagnostics | **9.92 s** |
| `python -m graphpathology arxiv --output outputs/publication-sync-arxiv` | Single-root ZIP; clean extraction compiled without BibTeX; **52 pages** | **5.53 s** |
| `python -m graphpathology verify-archive --historical --report outputs/publication-sync-archive.json` | **74,253 historical files unchanged**; available declarations/prefix scopes match; no new issues | **13.52 s** |

Publication timings refer to the final full renderer/trace run, including the printed inline-S9/source equivalence and original-head-identity checks. An early validation-only pandas deprecation was fixed; the final test and regression runs are clean. The earlier 84-pass/26-view/66-check report is retained in `docs/history/INITIAL_TEST_REPORT.md` and is not current validation evidence.

## Numerical evidence and limits

`PUBLICATION_VALIDATION.json` lists all checks and tolerances. All 122 author-package file SHA256 values match the input ZIP; all 35 hash-listed scientific exports remain byte-identical after rendering. All 18 figures are regenerated with the unmodified author scripts in independent output and baseline copies. Their actual line/scatter/image/patch arrays, masks/NaNs, scalar colors, color limits, labels and ordering agree at rtol=atol=1e-12. This verifies the integration against the current author baseline; it does not independently validate the original scientific experiments.

Independent checks protect all printed Table 2 metric cells and cohort medians, all 50 Table S1 gene ranks, all 150 printed Table S2 metrics and 25-patient ordering, Table S9 coverage cells/aggregation, Table 1 counts, Table S5 recorded audit cells, Table S7 support counts, Table S8 dimensions/counts, and definition-table schemas. The full per-gene CSV.gz remains intact. Figure 6's shared 98th-percentile cap is checked against the full export, including undisplayed capture-window rows, along with exact selected map points/errors. All 4,448 sampling rows are retained.

All 50-gene normalization, head-target and source/Wu calibration constants and the 512-feature scaler agree with archived NPZ/JSON originals. Synthetic Figure 2 matches the maintained production graph functions, exact undirected edge identities and three node permutations. Its rounded coordinate/distance exports justify a 2e-6-um distance tolerance; edge counts are exactly 94/22/98/85. Those are replacement illustration counts, not new patient results.

## Document and visual checks

All three document forms build from isolated source copies and regenerated figures. Final logs have no warnings, undefined references/citations, underfull/overfull boxes or duplicate destinations. Combined PDF has all 18 figure captions, 11 table captions and 29 bibliography entries. Supplement begins after References on page 30; page numbering is continuous through 52. `DOCUMENT_BUILD_REPORT.json` and `ARXIV_REPORT.json` record checks. The arXiv ZIP contains 34 TeX/figure/bibliography dependencies with only main.tex as a document root; numerical ancillary files are under anc/ and not TeX dependencies. Extra roots, compiled manuscript PDFs, logs and author notes are excluded.

All 52 combined pages were rendered and reviewed in contact sheets; the article/supplement share the same content sources. No clipped or overlapping content was found in that review. The author's whitespace/float layout is preserved, including its section-divider pages. This is not byte-identical recovery of the earlier 54-page document.

## Separate evidence

The base-main GitHub CI run `41262294332` passed before this revision. It is evidence for the merged initial revision, not for these new changes. A bounded CPU workflow and a separate full publication/TeX CI workflow accompany this PR; their actual status is reported on GitHub. Initial local records and authentication failures remain historical notes under docs/history, not current blockers.

No original fitted-head replay was possible. The three original heads and supporting fitted states remain missing; no research training, feature-extraction campaign, calibration search, full-cohort inference or historical replay was run. Default tests remain offline/CPU fixtures, with publication/TeX integration exercised separately. No historical tests were deleted or weakened to obtain a passing result; the previous asset-dependent limitations documented in the initial report remain.
