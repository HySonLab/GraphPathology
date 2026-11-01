SCIENTIFIC SOURCE DATA

These files contain the existing numerical results used by the main article
and Supplementary Information. They were prepared for readable publication
from the supplied research exports. No model was fitted and no prediction,
metric, numerical precision, or scientific data row was changed during this
editorial preparation.

TABLES
SuppTable1_gene_panels.csv: gene-panel inclusion, assay availability and order.
SuppTable2_per_patient.csv: patient-level results across the analytical cohorts.
SuppTable3_per_gene.csv.gz: full gene-level results for fresh external analyses.
SuppTable4_negative_controls.csv: development and external control results.
SuppTable5_integrity.csv: numerical prediction/metric replay and leakage checks.
SuppTable6_masks.csv: identifiers and design metadata for observation masks.
SuppTable7_distance_bins.csv: distance-resolved errors and support counts.

ADDITIONAL RESULTS
Bassiouni_FPS_effects.csv compares spatial farthest-point with random sampling.
Bassiouni_distance_trends.csv contains occupied-bin error-trend summaries.
Bassiouni_support_statistics.csv contains patient support counts and paired effects.
development_patient.csv contains the primary development-patient summaries.

FIGURE SOURCES
Fig*_source.csv, Fig3_map_source.csv and SFig*_source.csv retain their source
basenames; these numbers are export identifiers, not necessarily the revised
article's display numbers. Use the figure/source mapping supplied with the
revised package and the figure legends when selecting records. Several files
include additional panels, budgets or development settings beyond a displayed
subset. Their rows have been retained.

ANALYSIS ROLES
The analysis column uses descriptive names: HEST preprocessing, HEST local
predictors, HEST graph context, HEST graph controls, HEST representation
benchmark, HEST sparse reconstruction, HEST sampling, HER2ST transportability,
HER2ST calibration, Wu evaluation and Bassiouni replication. Not every table
contains every analysis. HER2ST transportability and calibration reuse the same
eight patients. The Wu study contributes three calibration patients and three
held-out patients; Bassiouni contributes 22 held-out patients. Repeated masks,
model seeds and reused patients are not additional biological replicates.

EDITORIAL COLUMN CHANGES
Internal execution-stage fields are named analysis and use the scientific
roles above. Run/mask/policy identifiers have consistent descriptive prefixes;
their identity mappings are bijective. Source paths, file digests, internal
row pointers, execution-only status fields and the experiment-status timeline
are omitted. Exact scientific numeric strings, biological gene identifiers,
patient/sample identifiers and the data-row order are unchanged. A patient or
sample label is never interpreted as an internal execution code.

Reconstruction outcomes and Boolean scientific integrity checks remain in
the release. A blank cell means that the original source did not provide that
field for the row; it is not an imputed zero. Constant-vector correlations
remain undefined. Missingness patterns and the differences between mean and
median aggregation should be interpreted using the article Methods and each
figure/table legend.

Panel and method labels retain their experimental meanings. The primary
development and fresh external analyses use 50 genes; HER2ST uses the shared
40-gene intersection. Do not pool scores from different panels. Molecular
interpolation and image-only prediction receive different information at
inference. Regional completion and distributed interpolation differ in mask
geometry, scored locations and interpolation weights.
