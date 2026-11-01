# Reference and claim audit

Prepared against the attached 83-page manuscript and PR #3 commit `ab0811134bb1d62719206ab11a43bd7c4dd5d813` on 1 November 2026. This file separates published-literature descriptions from results reported in the user's manuscript. It is not a new empirical validation.

## Eight new bibliography entries

### `han2022vig`
Vision GNN: An Image is Worth Graph of Nodes (2022).

General-domain image graphs: patches as nodes, graph aggregation and feature transformations; evaluated for recognition and detection. Does not establish molecular reconstruction performance.

Primary source: https://papers.nips.cc/paper_files/paper/2022/hash/3743e69c8e47eb2e6d3afaea80e439fb-Abstract-Conference.html

### `wang2018nonlocal`
Non-Local Neural Networks (2018).

Nonlocal visual feature aggregation; not claimed to be the same implementation as the bounded landmark model in GraphPathology.

Primary source: https://openaccess.thecvf.com/content_cvpr_2018/html/Wang_Non-Local_Neural_Networks_CVPR_2018_paper.html

### `pati2022hact`
Hierarchical graph representations in digital pathology (2022).

HACT-Net: hierarchical cell/tissue graph representations for histopathological classification; not spatial-expression prediction.

Primary source: https://doi.org/10.1016/j.media.2021.102264

### `chen2021patchgcn`
Whole Slide Images are 2D Point Clouds: Context-Aware Survival Prediction Using Patch-Based Graph Convolutional Networks (2021).

Patch-based spatial histology graphs for survival prediction. Correct venue is MICCAI 2021, not ICCV.

Primary source: https://doi.org/10.1007/978-3-030-87237-3_33

### `hu2021spagcn`
SpaGCN: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network (2021).

Expression, spatial location and histology for spatial domains and spatially variable genes.

Primary source: https://www.nature.com/articles/s41592-021-01255-8

### `dong2022stagate`
Deciphering spatial domains from spatially resolved transcriptomics with an adaptive graph attention auto-encoder (2022).

Adaptive graph-attention autoencoder for spatial representations/domains.

Primary source: https://www.nature.com/articles/s41467-022-29439-6

### `long2023graphst`
Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST (2023).

Graph-based contrastive learning for spatial clustering, integration and deconvolution.

Primary source: https://www.nature.com/articles/s41467-023-36796-3

### `song2023gntd`
GNTD: reconstructing spatial transcriptomes with graph-guided neural tensor decomposition informed by spatial and functional relations (2023).

Direct transcriptome reconstruction precedent: neural tensor decomposition with spatial and gene-functional graph regularization. Do not imply all earlier graph-ST work is clustering only.

Primary source: https://www.nature.com/articles/s41467-023-44017-0

## Other primary sources used for close-method comparisons

- GraphGPS: https://arxiv.org/abs/2205.12454 — local message passing, attention and structural/positional encodings; cite as existing architecture principles, not an exact reproduction claim.
- EGNN: https://proceedings.mlr.press/v139/satorras21a.html — Euclidean-equivariant graph networks; GraphPathology's full image pipeline is not thereby invariant.
- Hist2ST: https://doi.org/10.1093/bib/bbac297 — contextual image-to-expression modeling using Transformer and graph components.
- THItoGene: https://doi.org/10.1093/bib/bbad464 — local/global/graph-attention histology representation.
- FmH2ST: https://doi.org/10.1093/nar/gkaf865 — foundation features and spatial/semantic graph modeling.
- TESLA: https://pmc.ncbi.nlm.nih.gov/articles/PMC10246692/ — includes masked-spot recovery as well as enhanced spatial mapping.
- iStar: https://pmc.ncbi.nlm.nih.gov/articles/PMC11260191/ — high-resolution histology/ST integration.
- PanoSpace: https://doi.org/10.1038/s43588-025-00938-y — combines low-resolution ST, high-resolution histology and matched single-cell information.
- Sparser2Sparse: https://pubmed.ncbi.nlm.nih.gov/41915529/ — sparse-to-sparse within-sample imputation; DOI 10.1109/JBHI.2026.3679252.
- SCR2-ST: https://proceedings.mlr.press/v315/zhu26a.html — active sampling and prediction with single-cell priors.
- Spatial validation: https://doi.org/10.1111/ecog.02881 — structured validation for dependent observations.
- Cross-study histology-expression benchmark: https://doi.org/10.1038/s41467-025-56618-y.
- TISSUE: https://doi.org/10.1038/s41592-024-02184-y — calibrated expression-prediction uncertainty and downstream analyses, not merely patient-bootstrap intervals.

These are literature comparisons, not claims that GraphPathology ran new head-to-head experiments against all cited systems. The full list of the 34 works cited in the replacement sections appears in READING_COPY.md and the compiled preview; existing manuscript references elsewhere must not be removed arbitrarily.

## Manuscript-derived evidence to preserve

| Claim in the replacement | Source in supplied 83-page PDF | Interpretation constraint |
|---|---|---|
| Original 25 evaluation patients: Wu 3 and Bassiouni 22 | pp. 2–3, Results 2.1 / Table 1 | Follow-up controls reuse patients; do not increase biological n |
| 50-gene panel, complete-vector withholding | pp. 21, 25–27, Methods 4.1, 4.6–4.7 | Not whole-transcriptome recovery or prediction of an unassayed panel |
| 20% distributed r = 0.359 and 0.416 | p. 7, Table 2 | Cohort medians of patient summaries; original masks/regimes |
| Placement Δr = -0.1541, relative MSE +45.0%, all22 deteriorate | pp. 7,10, Results 2.3 / Table 3 | Same hidden H, anchor count and estimator; a specified central geometry |
| Uniform placement agrees in direction | p. 10, Table 3 | A second fixed estimator, not independent additional patients |
| Twenty molecular estimators, three budgets | pp. 11–12, Results 2.4 / Figure 8 | Complete declared grid, not external-outcome-selected winner |
| No consistent external gain from tested GNN predictors | pp. 12–16, Results 2.5 / Table 4 | Source-trained frozen features, low-resolution external images, specific controls |
| Two graph roles and prediction routes | pp. 22–27, Methods 4.3–4.7 | Molecular support is not a learned neural attention operator |
| Figure10 schematic unnecessary/overflowing | p.14, caption 'Graph predictors and information boundaries' | Remove current conceptual display, not quantitative experiments |
| AI assistance and author responsibility | pp.29–30 | Preserve approved disclosure |
| Public/author-local artifact availability | p.30 | Do not claim author-local weights/exports have been publicly released |

## Editorial safeguards

The contribution statement describes the framework and controlled design without a first-ever claim. Every Related Work paragraph closes with a task-, input-, representation- or evaluation-specific distinction. Negative findings and the scope of cost claims remain unchanged. New literature is used for context, not to fill gaps in the study's evidence.

## Repository integration verification

Merged all eight additions into the active bibliography and retained every original key. The existing `he2016` reference supplies ResNet; a duplicate proposed ResNet key was rejected before writing. Added `dosovitskiy2021vit` for the original ViT paper and retained `filiot2023` for Phikon. DOI/title/author/publisher metadata are captured in `REFERENCE_METADATA_AUDIT.json`; primary NeurIPS, Nature, author/institutional full-text and arXiv sources were checked. Automated direct access to CVF, OpenReview and some publisher pages was limited; the original author arXiv versions and publisher-registered metadata were used, without account login or access-term acceptance. No literature model is claimed as a newly executed benchmark.

Primary sources additionally checked: https://arxiv.org/abs/1711.07971 (Non-local); https://arxiv.org/abs/2010.11929 (ViT); https://arxiv.org/abs/2107.13048 (Patch-GCN); https://www.nature.com/articles/s41467-023-44017-0 (GNTD); https://www.nature.com/articles/s41467-023-36796-3 (GraphST). The Vision GNN bibliography URL uses its DOI instead of the long proceedings hash URL solely for typesetting. CMRCNet remains cited within the image-prediction paragraph before its closing task distinction.
