# Sparse Spatial Transcriptomics Reconstruction Reveals Coverage- and Geometry-Dependent Limits of Molecular Interpolation

[DRAFT — REVISE AFTER HUMAN WRITING]

## Abstract

Sparse spatial measurements may support reconstruction between measured locations, but their reliability depends on coverage and missing-region geometry. We compared frozen image-only prediction, spatial interpolation and fixed image–molecular fusion across development and external breast cancer datasets. Fresh validation included three held-out Wu patients and 22 independently sourced Bassiouni patients, using the same 50-gene panel. At 20% distributed coverage, median patient Pearson correlations were 0.359 and 0.416 for sparse-ST interpolation, versus 0.059 and 0.075 for image-only prediction. In the second replication, all patients met the prespecified joint correlation/error criterion at every tested coverage. Increasing distributed coverage improved reconstruction, whereas contiguous extrapolation reduced correlation and increased errors in every second-cohort patient. Distance-related error trends varied among patients and were not universally monotonic. Fixed image fusion did not provide reproducible incremental benefit. Public external images were downsampled, limiting conclusions about the potential of morphology. These results support coverage- and geometry-aware assessment of sparse molecular reconstruction, without establishing whole-transcriptome recovery, optimal sampling or clinical utility.

## Introduction

Spatial transcriptomics measures molecular variation together with tissue location, but incomplete sampling creates a distinction between measurements and reconstructed values. Reconstruction should therefore be evaluated under explicit observation masks and with molecular information withheld from unmeasured locations. [REF: spatial transcriptomics and missing-data evaluation]

Tissue images offer a complementary source of information. A useful benchmark asks whether a frozen image predictor is improved by sparse local molecular measurements and whether images add information after spatial interpolation already uses those measurements. Correlation alone is insufficient because predictions can have spatially coherent rankings while remaining poorly calibrated in absolute target units. [REF: image-to-expression benchmarks and calibration]

Coverage and geometry address different scientific questions. Distributed anchors surround many hidden locations; a contiguous missing region can require extrapolation across a molecular boundary. Nearest-anchor distance may describe part of this difficulty, but it does not capture tissue heterogeneity or establish a causal reconstruction horizon. [REF: spatial interpolation and spatial validation]

Cross-platform evaluation further changes the target distribution, capture support and image resolution. Development, inspected transportability datasets and genuinely held-out validation must remain distinguishable. Reusing a cohort for calibration-development can be informative without making it a fresh validation cohort. [REF: external validation and dataset shift]

Here we assemble a frozen sequence of breast tissue experiments around a conservative claim: sparse measured expression can improve hidden-location reconstruction, subject to coverage and geometry limitations. We retain negative graph, sampling and image-fusion findings, and report patient-level replication with exact numerical provenance.

## Results

### Distributed molecular anchors improved reconstruction over the frozen image baseline

On the development primary patient, at 20% distributed coverage, mean-over-mask median gene Pearson correlations were 0.695 for H&E, 0.792 for selected ST interpolation and 0.822 for fixed fusion. Corresponding mean gene MSE values were 0.911, 0.268 and 0.279. Thus the development fusion correlation increase did not imply lower error than interpolation. These are one-patient development findings, not an external population estimate (Figure 2; Tables 1–2; ledger entries keyed Development/0.2/V1).

Fresh Wu validation retained all original primary genes. At 20% V1, sparse ST achieved r/MSE/MAE of 0.359044/0.098712/0.217868, compared with 0.059037/2.915478/1.557786 for H&E. Fixed fusion yielded 0.334853/0.100191/0.232128 (Figure 5). C9 remains frozen as MIXED_FRESH_EXTERNAL_VALIDATION: the sparse-ST comparison was supportive, but affine harmonization could not provide an identifiable interpolation advantage and the public image inputs limited H&E inference.

The second replication used one prospectively selected first section for each of 22 Bassiouni patients, totaling 30,003 source-filtered tissue locations. No second-cohort patient was used for training, calibration or selection. At 20% V1, sparse ST achieved 0.415862/0.092256/0.203274, versus 0.074917/1.967339/1.212010 for H&E and 0.406872/0.102208/0.228760 for fixed fusion. All 22 patients met the prespecified sparse-ST criterion at all three budgets: a Pearson gain of at least 0.02, relative MSE reduction of at least 2%, and non-worsening MAE. The replication is PASS_LIMITED_SECOND_REPLICATION because public H&E inputs contained only 24–29 native pixels across the unchanged 112-µm patch (Figure 6; Tables 1–2,4).

### Coverage benefits were strongest for distributed measurements

Across Bassiouni patients, median patient ST correlations at 10%, 20% and 40% V1 coverage were 0.380, 0.416 and 0.465; MSE values were 0.100, 0.092 and 0.085. Every patient improved both correlation and MSE from 10% to 40%. The Wu coverage-response was directionally consistent (Figures 5–7; Table 3). These tested budgets do not establish an optimal acquisition level or a monetary saving.

### Contiguous extrapolation and distance exposed reliability limits

At 20% coverage, Bassiouni median patient ST correlation fell from 0.416 under V1 to 0.171 under V2. MSE rose from 0.092 to 0.150; all 22 patients worsened in correlation, MSE and MAE. The median paired relative MSE increase was 52.8%. This is a prespecified protocol contrast: the scoring region and frozen interpolation family both differ, so it cannot isolate geometry as a causal factor (Figures 3,5–6; Table 3).

Distance effects were heterogeneous. In the second cohort at 20%, the farthest occupied bin had greater MSE than the nearest in 15/22 V1 patients and 20/22 V2 patients. Strictly increasing occupied-bin curves occurred in only 6/22 and 13/22 patients, respectively. Remote bins can contain few locations or different subsets of patients; complete counts and correlation eligibility are retained (Supplementary Table 7). We therefore describe distance dependence without asserting a universal monotonic curve or reliability radius.

### Image fusion and affine calibration did not establish a portable extra benefit

HER2ST retained a common 40-gene subset and serves as transportability development. Calibration reduced scale-related image/fusion error relative to zero-shot transport, but an invertible per-gene affine transform commutes with normalized linear interpolation after returning to common units. A distinct harmonization benefit cannot be assigned to that interpolation identity (Figure 4; Supplementary Figure 8).

In Bassiouni 20% V1, fusion increased correlation in only 3/22 patients and improved neither MSE nor MAE in any patient. The median paired correlation change was −0.0108, with a median relative MSE increase of 10.7%. Wu likewise showed no reproducible joint benefit. This limits the tested fixed fusion rule, especially when transporting a Wu-derived scaler to the second cohort; it does not exclude benefits from higher-resolution images or separately justified future methods (Figures 5–7).

### Controls and frozen provenance support a limited reconstruction claim

Anchor-vector and anchor-position permutations degraded reconstruction on the prospectively designated control patient; image permutation had a substantially smaller effect (Supplementary Figures 9–10). Earlier context/topology and sampling experiments remain mixed and do not justify claims that graph topology is essential or morphology-guided sampling is reliably optimal (Supplementary Figures 1–8; claim matrix).

The second replication comprised 2,904 clean and 18 control bundles, with 8,766 method predictions. All 17,532 hidden-label NaN/random corruption comparisons were invariant. Independent replay reproduced all predictions and 35,064 per-gene metric arrays exactly, and reconstructed all 968 label-free mask trajectories (Supplementary Tables 5–6). C0–C9 artifacts are checked read-only through their chained manifests.

## Methods

### Cohorts, identities and analytical roles

Table 1 and Supplementary Table 1 give source studies, patient counts, resolution and exact gene intersections. Development retained the frozen C0 split: NCBI783/NCBI785 training, TENX95 validation, TENX99 primary test. Conditional leave-one-patient-out analyses in earlier checkpoints are sensitivity analyses after development choices, not pristine cross-validation. HER2ST uses eight patient-specific first sections; C7 zero-shot and subsequent C8 calibration-development are distinct analytical roles on the same inspected cohort. Forty-gene HER2ST summaries are not pooled with fifty-gene primary summaries.

Wu local Visium patients 1142243F, CID4465 and CID4535 formed the frozen calibration set; 1160920F, CID4290 and CID44971 were held out. The second fresh study is GSE210616, Bassiouni et al., Cancer Research (2023), doi:10.1158/0008-5472.CAN-22-2682. All published patient IDs with the first numbered section were selected from source metadata before target access. This yields BASS_P01–P09, P11–P15 and P17–P24; gaps are source numbering, not exclusions based on outcomes. The original publication's reference/validation groups are source-study roles, not calibration roles here. Anonymized source identity, repository identifier searches and available image hashes support freshness; they cannot prove patient identity across unrelated anonymized studies.

### Inputs and preprocessing

C10 downloaded the filtered feature-barcode H5, tissue image, tissue positions and scale-factor JSON for each selected section. Barcode joins and source in-tissue flags determine inclusion; no expression-derived quality selection, pathology annotations or deconvolution results are inputs. Target counts must be finite, nonnegative integers. Targets use log1p counts followed by unchanged C0 training-gene centering/scaling. No destination library-size or capture-area correction is fitted. Exact symbol intersection preserves the original 50-gene order; optional 280-panel overlap is documented without secondary model evaluation.

Coordinates use a least-squares mapping of the source Visium hexagonal 100-µm lattice into image coordinates. Residual and anisotropy thresholds are in the locked qualification code. Physical agreement tests coordinate consistency and does not independently establish cellular registration. All image overlays were visually checked before outcomes. Source resolution, downloaded checksums, metadata hashes, and patch hashes are retained. External downsampling is not corrected by superresolution.

### Frozen image predictor and molecular reconstruction

The image branch is the unchanged ResNet18 IMAGENET1K_V1 representation with frozen C4 D1 primary MLP heads, averaged over training seeds 0,1,2. The heads use hidden width 256 and dropout 0.1; complete checkpoint and preprocessing identities are in C10A/SPEC_LOCK.json and per-patient EXTRACTION.json. Patches span nominally 112 µm, rounded at native resolution, resized by LANCZOS to 224 RGB pixels and passed through the official saved encoder transform. Training is not repeated in C10. Detailed earlier training schedules and checkpoint choices remain in C4's frozen configuration; this draft makes no new training-hyperparameter claims.

M0 is the frozen image-only prediction. M1 uses the exact C9 molecular rule: up to 16 nearest observed anchors, normalized inverse-distance-squared weights in V1 and an unweighted nearest-anchor mean in V2. No external target fitting is required. M3 combines 0.8 M1 with 0.2 calibrated image prediction. Calibration is the already-frozen C9 T1 scaler, with source moments from NCBI783/NCBI785 and destination moments from the three Wu calibration patients; C10 does not fit a new scaler or weight. Source/destination moments are patient-weighted, as recorded in the frozen object.

M2, the affine-transform/interpolate/inverse-transform equivalent of M1, is checked numerically at every C10 prediction but omitted from duplicate metric tables. C9 raw Spearman differences can arise from rank ties affected by floating-point roundoff despite prediction equivalence; these are not biological harmonization gains. MSE and MAE are compared in common C0 normalized log1p units, not in incomparable native transformed scales.

### Masks, coverage, geometry and secondary sampling

Before evaluating expression, all trajectories were frozen. Budgets are 10%,20%,40% of the full retained tissue, using floor(N×budget+0.5) anchors. Random sampling uses seeds 0–9 with nested prefixes. Secondary spatial farthest-point sampling is deterministic and label-free. It is stored as a sampling strategy alongside V1/V2; the plan's V3 label is not treated as a third missing-region geometry. FPS selection starts at the eligible point farthest from the eligible-coordinate centroid.

V1 measures distributed anchors and scores all unmeasured locations. V2 excludes compact regions containing the nearest 40% of locations around centers at the 25th,50th,75th percentiles of the first coordinate principal-component rank. Anchors lie outside the forbidden region, and scoring uses that region only. All other unobserved labels are also removed from inputs. Three V2 geometries are averaged equally. C5's earlier block and capture-window masks differ and are labeled separately. Differences between V1/V2 include support, scoring composition and selected interpolation family.

### Metrics and aggregation

For each gene, Pearson and Spearman correlations, MSE and MAE are computed only on hidden locations. Panel summaries include gene means/medians and fractions above Pearson thresholds 0,0.1,0.2,0.3. Undefined correlations are excluded from correlation summaries and do not pass fraction thresholds. Within each patient and geometry, random seeds are averaged; geometry means are averaged equally; headlines are medians across patients. Mask seeds and spots are not independent biological replicates. We report descriptive patient effects without new p-values or post-hoc significance tests. Median paired differences need not equal differences between cohort medians.

Nearest-anchor distance bins are [0,150),[150,300),[300,600),[600,1200),[1200,2400),[2400,infinity) µm. Error is reported for occupied bins; bin correlation requires at least 20 locations. Missing bins are not filled with zero. Supplementary Table 7 retains counts for each run. Curves describe different location subsets and do not estimate an adjusted causal distance effect.

### Controls, prospective decision rule and software

At 20% on BASS_P01, random-mask seed0, V1 geometry0 and V2 geometry1, NC1 permutes anchor expression-vector assignments, NC2 permutes anchor-position assignments and NC3 permutes image predictions, using control seeds0,1,2. NC1/NC2 are related association nulls, not independent mechanistic tests. Every prediction is repeated with all unobserved labels replaced by NaNs and by independent random numbers. Saved barcodes, masks, predictions, source hashes and exact metric replays support the audit.

The prospective sparse-ST criterion requires Δ median Pearson≥0.02, relative mean MSE reduction≥0.02 and non-worsening MAE. The cohort gate requires at least ceil(2N/3) patients meeting this criterion at at least two budgets, positive median benefits at every budget, and the same patient-support threshold for correlation/MSE improvement between 10% and 40%; distance, V2 and controls must be fully reported and technical checks pass. For N=22 the threshold is15. Resolution-limited qualification prevents an unrestricted PASS. The separate fusion gate is frozen in configs/protocol.json and is not required for sparse-ST support.

C10 uses Python, NumPy, SciPy, pandas, h5py, PyTorch/torchvision and Matplotlib, with CPU execution and four computational threads. Exact runtime versions, processor description, per-run runtime/RAM and command logs are in C10B/logs/ENVIRONMENT.json and logs/runs. No new encoder, graph model, biological interpretation or sampling optimization is introduced. All source code, figures, source data, table scripts and manifest hashes are stored in the repository; bulky raw arrays and prediction archives are retained locally and hashed.

## Discussion

The central observation is that distributed molecular anchors support reconstruction beyond the frozen image baseline across two fresh studies. The second study increases patient diversity, but it remains a targeted-panel, breast-tissue result under specified input and mask protocols. It should motivate coverage-aware reporting rather than claims of recovering all molecular features.

Missing-region geometry is a major practical limitation. The contiguous protocol worsened outcomes throughout the second cohort, while distance curves were neither universally smooth nor monotonic. A global correlation cannot identify where reconstructed values remain unreliable; spatial error maps and per-bin counts complement the primary metrics.

Image-fusion findings are negative for the frozen method. Their interpretation is constrained by downsampled public images, capture-platform differences and the deliberate reuse of a calibration object from another cohort. These limitations prevent broad conclusions about morphology, yet do not justify hiding the observed increases in error. Affine interpolation equivalence similarly limits the interpretation of calibration as a separate molecular reconstruction method.

Sampling and graph-control results argue against adding unsupported architectural or acquisition claims. Some spatial placement effects appeared in development, but robust gains across contiguous regimes and patients were not established. The present study does not demonstrate essential topology, a universally optimal strategy, monetary savings, clinical utility, causal biological mechanisms or single-cell resolution.

The next human-writing step is to verify literature citations, dataset permissions/acknowledgments and the exact scope of each claim, then refine prose and the conceptual schematic. Higher-resolution external H&E would strengthen future morphology comparisons, but additional cohorts are not required to begin writing this deliberately limited sparse-ST report. This execution stops after C10.

## Data and code availability (draft)

Repository: https://github.com/HySonLab/GraphPathology . Dataset accessions, manifests and source URLs are retained in C0/C7/C9/C10 qualification records. Results are traceable through RESULTS_LEDGER.md, SOURCE_INDEX.json, exact figure CSVs and chained frozen manifests. Human authors must add final repository release/DOI, licenses, author contributions and acknowledgments before submission.
