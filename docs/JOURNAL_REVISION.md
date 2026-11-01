# Separate Scientific Reports journal experiments

The original arXiv v1 package at `manuscript/revised/`, its three published PDFs, 29 reproduction targets and archived C0–C10 evidence remain unchanged. New journal sources and PDFs live at `manuscript/journal_revision/`; all new small numerical results, protocols and provenance are in `experiments/journal_revision/`. These analyses use previously inspected cohorts and are not another blind external validation.

The main protocol was committed before new comparator outcomes: `2e70b30722b56317e82005999ed90d8faeb2537a`. Learned implementation details, including the explicitly exploratory residual-anchor family, were committed before new neural outcomes: `7e3cc167fd066e8c9a104cfc696f6cbd4e343682`. Architecture widths, hyperparameters, controls, seeds and stopping rules were fixed; no external winner was selected from the estimator grid or the neural runs.

## Empirical status

| Analysis | Status and scope |
|---|---|
| A: 20 matched-information molecular estimators | COMPLETE for three Wu patients, all ten original random trajectories and 10/20/40% budgets; new Bassiouni predictions BLOCKED by raw-source transfer failure |
| B1: two operators crossed with original V1/V2 | COMPLETE for Wu, every original random mask and scored set, all budgets and three V2 geometries; Bassiouni BLOCKED |
| B2: identical H, equal anchor count | COMPLETE for Wu, one preregistered central region, two operators, ten seeds, all budgets; Bassiouni BLOCKED |
| C: patient uncertainty | COMPLETE original-export reanalysis, 22 Bassiouni and three Wu patients; bootstrap patients only (10,000, seed20261101); no primary Wu bootstrap claim |
| D: predeclared gene maps | COMPLETE actual Wu1160920F V1/central-V2, ERBB2/EPCAM/FOXA1/POSTN; BASS_P01 BLOCKED; no synthetic replacement |
| E: new local/mean/GCN/GraphSAGE heads | Real HEST inputs and official encoder acquired; see `results/hest_learned/EVALUATION_STATUS.json` for executed completion and exact all-seed outputs |
| Harmonic extension, FPS B2, novel/multihop architectures | SKIPPED in protocol V1, not failed experiments |
| Historical selected heads | Still missing; new journal fits never substitute for their original identities |

Wu produced 2,700 evaluated runs and 5,400 all-unobserved-label NaN/random invariance checks. V1 operator rows are reused in B1 rather than rerun or counted twice. At 20% identical-target coverage, region-exclusion versus distributed IDW reduced patient correlation by 0.106–0.245 and increased MSE by 20.1–45.4%. Uniform interpolation also worsened in each patient. All individual effects and all candidates are retained; no best external method is promoted.

The original Bassiouni 20% V1 export reanalysis gives fusion-minus-anchor median delta r=-0.010830 (95% patient percentile interval -0.018738 to -0.008136) and relative MSE +10.7078% (6.6170–14.5048%). These are descriptive post-v1 patient intervals, not p-values, simultaneous bands or new predictions.

## Acquisition and terms

`DATA_STATUS.json`, `ACQUISITION_PREFLIGHT.json`, `MIRROR_ATTEMPTS.json` and the ignored acquisition manifests record source URLs, checksums, attempts and license checks. Wu's three official Zenodo archives match published MD5 checksums and have new SHA256 records (CC-BY4.0). All nine pinned HEST objects match the original full SHA256 values (CC-BY-NC-SA4.0; academic research). Official ResNet18 ImageNet1K_V1 weights match the advertised SHA256 prefix and have a recorded full digest. HER2ST metadata/license were verified (CC-BY4.0); new calibration was not part of this protocol.

GEO returned HTTP403/500, incomplete bodies, and zero-byte direct anonymous FTP transfers despite nonzero SIZE responses. Empty transfers fail the original full SHA256 requirements. No replacement patients, normalized substitutes or fabricated predictions were used. The 37.7GB whole-series bundle was not downloaded because it exceeds the predeclared initial20GiB quota. Recovered original exports allow C to proceed independently. To complete the blocked cases later, provide the exact selected raw files with the hashes in the original C10 download manifest, then rerun the same frozen protocol in a new result directory.

Raw data, upstream weights, processed embeddings and new checkpoint binaries remain in ignored caches. Numeric derived maps are actual Wu data, attributed to Zenodo4739739 under CC-BY4.0; they are not complete raw tissue images or count matrices. No new blanket software license or journal acceptance is asserted.

## Reproduction

Install the existing validation environment and package; the exact executed package versions are in `experiments/journal_revision/ENVIRONMENT.json`. Each command is explicit: imports, help and tests do not fetch data or train. Use new output directories for new analyses; original evidence paths are rejected. Downloads use two workers and resumable partial files with digest verification. Independent acquisition jobs initially overlapped (up to four aggregate workers); that resource deviation is recorded and changed no analytical choices. Neural operations use deterministic CPU float32 and four torch/BLAS threads; no paid compute was initiated.

```bash
python -m graphpathology.journal acquire --cohorts wu bassiouni --dry-run
python -m graphpathology.journal acquire --cohorts wu bassiouni
python -m graphpathology.journal original-uncertainty --output experiments/journal_revision/results/my_original_reanalysis
python -m graphpathology.journal anchors --cohort wu --output experiments/journal_revision/results/my_wu_anchors
python -m graphpathology.journal anchors --cohort bassiouni --output experiments/journal_revision/results/my_bassiouni_anchors

python -m graphpathology.journal --data data/journal_revision/hest_training acquire --cohorts hest resnet18
python -m graphpathology.journal --data data/journal_revision/hest_training learned-prepare
python -m graphpathology.journal --data data/journal_revision/hest_training learned-train \
  --output experiments/journal_revision/results/hest_learned \
  --checkpoints "$(pwd)/checkpoints/journal_revision/hest_learned"
python -m graphpathology.journal --data data/journal_revision/hest_training learned-evaluate \
  --output experiments/journal_revision/results/hest_learned \
  --checkpoints "$(pwd)/checkpoints/journal_revision/hest_learned"
python -m graphpathology.journal paper --build-output outputs/my_journal_build
```

The paper command reads the committed numerical exports and compiles isolated journal sources; it does not train. Actual map exports are included. Existing trained/feature outputs resume only with matching stored digests; checkpoint files are intentionally not committed. A complete neural rerun therefore downloads approximately5.6GB of HEST inputs plus a46.8MB encoder, recomputes frozen features, and trains new heads. Publication PDFs and small result tables are reviewable without those binaries.

## Numerical limits and provenance

Read `NUMERICAL_DIAGNOSIS.md`: source lattice recalibration differed by <=5.46e-12um and changed nearest-neighbor ties. After raw coordinate/barcode checks, using the exact preserved label-free mask coordinates restores every primary per-gene Pearson/MSE/MAE value in the three initial real replay cases within rtol=atol=1e-10. A subsequent all-reference audit matched all 360 original random-reference runs and 54,000 gene-level primary values; maximum absolute Pearson difference was 3.46e-14, below the declared tolerance. Raw Spearman is not numerically identical because floating-point near ties alter ranks; it remains a secondary flagged export and supports no improvement claim. No historical values were changed or tolerance quietly relaxed.

Each completed patient/model records protocol/code commits, actual code-file hashes, input and mask hashes, environment, seeds, timings, memory and output digests. The initial empty original-export analysis is explicitly INVALID; only corrected `original_v1_reanalysis_r2` is used. Confidence intervals never resample genes/spots as independent units. New controlled neural analyses use one previously used evaluation patient; training seeds and masks are computational variability, not biological replication. Fivefold fewer observed locations remains conditional variable-cost arithmetic, not demonstrated total assay savings.
