# Current Scientific Reports matched-backbone manuscript

This author-review version centers controlled measurement placement in 22 previously evaluated patients. Its expanded Introduction and separate Related Work explain the distinction between image context and observed molecular support. The abstract is 178 words and uses the general GNN description without an architecture inventory. The original molecular headline correlations remain correctly attributed.

The matched follow-up separates graph versus local prediction within ResNet18 and Phikon from reconstruction with identical molecular anchors and from anchor-only IDW. GCN-assisted reconstruction improves all three endpoints over local-assisted prediction in all 22 Bassiouni patients with both encoders, while IDW remains stronger on the cohort-median endpoints. Image-only comparisons and other heads include mixed and negative findings. Neighborhood averaging and shuffled controls qualify topology-specific interpretation. All earlier empirical evidence, including the nonlocal and latent-coordinate negative findings, remains present.

The [comparison protocol and reproduction guide](../../../docs/MATCHED_BACKBONE_COMPARISON.md) describe actual source/validation-guided selection and resource revisions. This is exploratory follow-up on inspected patients, not fresh blind validation. The original arXiv v1 and earlier journal packages remain preserved. The pre-edit current-source snapshot is [commit ab0811134](https://github.com/HySonLab/GraphPathology/tree/ab0811134bb1d62719206ab11a43bd7c4dd5d813/manuscript/journal_revision/transformer_egnn).

`CURRENT_DISPLAY_MAP.csv` is the current numbering and page map. `DISPLAY_BASELINE.json` and `DISPLAY_CHANGES.json` document deliberate changes. Historical V1, JOURNAL and EXTENSION target ledgers remain historical regression records. The conceptual information-boundary schematic is archived outside this active package; no empirical display was removed. The current manuscript has 40 figures, 24 tables and 40 cited references. The verified PDFs have 94 combined pages, 42 article pages and 52 supplement pages, with no final TeX diagnostics. The build report records the checks.

Build maintained sources without research downloads or model execution:

```bash
python -m graphpathology.backbone_paper build --output outputs/matched-manuscript
python -m graphpathology.backbone_selected_verify --output outputs/matched-integrity.json
python -m graphpathology.editorial check --output outputs/prior-integrity.json
```

A public checkout uses `--committed-only` for integrity checks; absent author-local states and full gene arrays are explicitly counted, not replayed. Source output and extraction directories must be new:

```bash
python -m graphpathology.backbone_paper package --output outputs/matched-source.zip --check-directory outputs/matched-source-check
```

The standalone source ZIP contains fresh `.bbl` files, vector assets, editable tables, source summaries and the exact pure plotting/table functions. `bash build.sh` builds all three PDFs without BibTeX; use `bash build.sh --bibtex` after editing references. `python scripts/regenerate_backbone_figures.py` regenerates the new comparison displays from `source_data/backbone/` with NumPy, pandas and Matplotlib, without model weights or research downloads. `scripts/regenerate_extension_figures.py` maintains earlier displays and removes the retired schematic from the active path. Historical research assembly cannot overwrite this maintained opening or restore the schematic.

See `EDITORIAL_CHANGELOG.md`, `CLAIM_EVIDENCE_AUDIT.md` and `AUTHOR_RELEASE_CHECKLIST.md`. Full gene exports and head-state bundles remain author-local; no release or DOI is claimed. The original three primary image-head binaries remain missing. Sole-human authorship, the approved distinction between author direction and AI assistance, and the Editorial Board disclosure remain intact. PR #3 stays draft and unmerged; no submission occurred.

The article describes the frozen Phikon encoder in journal prose. Exact reproducibility identity remains `owkin/phikon`, revision `057cc0295895c2df3dd7681a89680da6015cbefe`, in the unchanged protocol configuration and provenance manifests. The original frozen ResNet18 evaluation and later matched-backbone Phikon evaluation have separate caption scopes.
