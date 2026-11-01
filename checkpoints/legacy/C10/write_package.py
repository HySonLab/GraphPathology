"""Render traceable prose, ledgers, claims and figure/table documentation."""
from synthesis import *
import shutil

def main():
 t2=pd.read_csv(M/'tables/Table2_reconstruction.csv');t1=pd.read_csv(M/'tables/Table1_cohorts.csv');t3=pd.read_csv(M/'tables/Table3_coverage_geometry.csv');stats=pd.read_csv(M/'data/C10_support_statistics.csv');audit=json.loads((B/'logs/INTEGRITY_AUDIT.json').read_text());gate=json.loads((B/'logs/SCIENTIFIC_GATE.json').read_text());ledger=[]
 def add(claim,value,metric,cohort,patients,panel,coverage,regime,method,path,row,field,upstream=''):
  ident=f'L{len(ledger)+1:04d}';ledger.append({'ledger_id':ident,'claim':claim,'value':value,'metric':metric,'cohort':cohort,'patient_scope':patients,'gene_panel':panel,'coverage':coverage,'regime':regime,'method':method,'source_checkpoint':path.split('/')[0],'source_file':path,'source_row_or_field':f'row {row}; {field}' if row else field,'upstream':upstream,'source_sha256':sha(ROOT/path)});return ident
 ids={}
 for i,r in t2.iterrows():
  for met in ['median_pearson','mean_mse','mean_mae','patient_count','patient_support_vs_M0_all_metrics']:
   ids[(r.cohort,r.coverage,r.regime,r.method,met)]=add('Descriptive reconstruction summary',r[met],met,r.cohort,f'{r.patient_count} patients',str(r.primary_gene_count),r.coverage,r.regime,r.method,'C10/MANUSCRIPT/tables/Table2_reconstruction.csv',i+2,met,f'{r.source_file}; rows {r.source_rows}')
 for i,r in stats.iterrows():add(r.statistic,r.value,r.unit,'Bassiouni second replication','22 evaluation patients','50',.2,'V1 or V2 as named','M1/M3','C10/MANUSCRIPT/data/C10_support_statistics.csv',i+2,'value',f'{r.source_file}; rows {r.source_rows}; {r.source_field}')
 for i,r in t1.iterrows():
  for met in ['patient_count','evaluation_patient_count','location_count','evaluation_locations','H&E resolution','primary_gene_overlap']:
   add('Cohort description',r[met],met,r.cohort,'Table1 scope','50 or40 as specified','NA','NA','NA','C10/MANUSCRIPT/tables/Table1_cohorts.csv',i+2,met,f'{r.source_file}; rows {r.source_rows}')
 for i,r in t3.iterrows():
  for col in t3.columns:
   if 'Pearson' in col or 'pearson' in col or 'mse' in col:add('Coverage/geometry summary',r[col],col,r.cohort,'All listed evaluation patients',r.primary_gene_count,'10/20/40% as named','V1/V2','M1','C10/MANUSCRIPT/tables/Table3_coverage_geometry.csv',i+2,col,f'{r.source_file}; rows {r.source_rows}')
 for field,val in audit.items():add('C10 integrity audit',str(val),field,'Bassiouni second replication','22','50','all','all','all','C10/C10B/logs/INTEGRITY_AUDIT.json',None,field)
 for field in ['required_patients','coverage_response_patients','patient_meaningful_budget_counts']:add('Prospective replication gate',str(gate[field]),field,'Bassiouni second replication','22','50','10/20/40%','V1','M1 versus M0','C10/C10B/logs/SCIENTIFIC_GATE.json',None,field)
 for field,val in cfg().items():
  if field not in ['inputs']:add('Frozen scientific method',str(val),'protocol','C10','22','50','protocol','protocol','protocol','C10/configs/protocol.json',None,field)
 old=json.loads((ROOT/'C9/C9_STATUS.json').read_text())
 for field,val in old['technical_audit'].items():add('C9 audit',str(val),field,'Wu','3 heldout','50','all','all','all','C9/C9_STATUS.json',None,'technical_audit.'+field)
 def val(cohort,method,metric='median_pearson',budget=.2,reg='V1',digits=3):
  q=t2[(t2.cohort==cohort)&(t2.coverage==budget)&(t2.regime==reg)&(t2.method==method)];assert len(q)==1;return f'{q.iloc[0][metric]:.{digits}f}'
 C='Bassiouni second replication';W='Wu fresh validation';D='Development'
 def trip(cohort,method,budget=.2,reg='V1'):return '/'.join(val(cohort,method,met,budget,reg,6) for met in ['median_pearson','mean_mse','mean_mae'])
 claimrows=[
 ('CL1','Sparse ST reconstructs hidden expression better than frozen H&E alone.','SUPPORTED_WITH_SCOPE','C5;C7;C9;C10','C10/MANUSCRIPT/tables/Table2_reconstruction.csv',f'20%V1: Wu ST r={val(W,"M1")}, H&E={val(W,"M0")}; Bassiouni ST={val(C,"M1")}, H&E={val(C,"M0")};22/22joint benefit atall3budgets','Open coverage; targeted panels; downsampled external H&E; not universal under contiguous missingness.','Distributed sparse measurements consistently improved reconstruction over the frozen image baseline in the two fresh validation studies.'),
 ('CL2','Increasing measured coverage improves reconstruction.','SUPPORTED_V1','C5;C7;C9;C10','C10/MANUSCRIPT/tables/Table3_coverage_geometry.csv',f'Bassiouni V1 r: {val(C,"M1",budget=.1)}, {val(C,"M1")}, {val(C,"M1",budget=.4)};22/22r andMSE improve40vs10','V2error not monotonic; no cost or optimal-budget conclusion.','Greater distributed coverage improved correlation and error across fresh validation patients.'),
 ('CL3','Greater anchor distance generally reduces fidelity.','PARTLY_SUPPORTED','C5;C6;C9;C10','C10/C10B/C10B_DISTANCE_TRENDS.csv;C10/MANUSCRIPT/supp_tables/SuppTable7_distance_bins.csv','C10farthest-bin MSE worse:15/22V1,20/22V2;strictincrease:6/22V1,13/22V2','Sparse remote bins and composition differences; distance observational; no universal cutoff.','Errors were often larger farther from anchors, with heterogeneous and nonmonotonic patient-level curves.'),
 ('CL4','Contiguous extrapolation is harder than open interpolation.','SUPPORTED_PROTOCOL_CONTRAST','C5;C6;C9;C10','C10/MANUSCRIPT/tables/Table3_coverage_geometry.csv;C10/C10B/C10B_V2_DEGRADATION.csv',f'C10ST20% r V1={val(C,"M1")},V2={val(C,"M1",reg="V2")};all22worse r/MSE/MAE','Scoring population and frozen interpolation family change; C5geometry differs; not a geometry-only causal effect.','The prespecified contiguous-extrapolation protocol reduced fidelity relative to distributed interpolation.'),
 ('CL5','H&E provides reproducible incremental value beyond sparse ST.','NOT_SUPPORTED_IN_FRESH_VALIDATION','C5;C8;C9;C10','C10/MANUSCRIPT/data/C10_support_statistics.csv;C9/C9B/C9B_PATIENT_EFFECTS.csv','At20%V1 fusion r improves1/3Wu,3/22Bassiouni;MSE/MAE improve0patients inbothstudies','Frozen fusion, transported calibration and limited image resolution constrain inference; not evidence that H&E can never help.','The tested fixed fusion did not provide reproducible incremental benefit across correlation and error metrics.'),
 ('CL6','Sampling optimization reliably improves reconstruction.','MIXED_NOT_ESTABLISHED','C6;C9;C10','C6/C6A/C6A_REPORT.md;C10/MANUSCRIPT/data/C10_FPS_effects.csv','C6MIXED_SAMPLING_EFFECT; external FPS secondary only','Some interpolation-only open benefits; no robust morphology-guided or contiguous advantage; no adaptive claims.','Simple spatial placement sometimes helped, but robust sampling optimization was not established.'),
 ('CL7','Graph topology is essential.','NOT_ESTABLISHED','C2;C3;C4','C3/C3A/C3A_REPORT.md;C3/C3_STATUS.json','MIXED_CONTEXT_SIGNAL;MIXED_SPATIAL_SPECIFICITY;MIXED_SPATIAL_RECHECK','Randomized topologies retained substantial performance; not tested as mechanism inC10.','Earlier graph controls did not establish a topology-specific requirement for reconstruction.')]
 claims=pd.DataFrame(claimrows,columns=['claim_id','claim_text','status','supporting_checkpoints','supporting_files','key_numbers','limitations','recommended_manuscript_wording']);claims.to_csv(OUT/'CLAIM_EVIDENCE_MATRIX.csv',index=False)
 # Narrative uses only numbers in tables, fixed protocol or audit fields registered above.
 draft=f'''# Sparse Spatial Transcriptomics Reconstruction Reveals Coverage- and Geometry-Dependent Limits of Molecular Interpolation

[DRAFT — REVISE AFTER HUMAN WRITING]

## Abstract

Sparse spatial measurements may support reconstruction between measured locations, but their reliability depends on coverage and missing-region geometry. We compared frozen image-only prediction, spatial interpolation and fixed image–molecular fusion across development and external breast cancer datasets. Fresh validation included three held-out Wu patients and 22 independently sourced Bassiouni patients, using the same 50-gene panel. At 20% distributed coverage, median patient Pearson correlations were {val(W,'M1')} and {val(C,'M1')} for sparse-ST interpolation, versus {val(W,'M0')} and {val(C,'M0')} for image-only prediction. In the second replication, all patients met the prespecified joint correlation/error criterion at every tested coverage. Increasing distributed coverage improved reconstruction, whereas contiguous extrapolation reduced correlation and increased errors in every second-cohort patient. Distance-related error trends varied among patients and were not universally monotonic. Fixed image fusion did not provide reproducible incremental benefit. Public external images were downsampled, limiting conclusions about the potential of morphology. These results support coverage- and geometry-aware assessment of sparse molecular reconstruction, without establishing whole-transcriptome recovery, optimal sampling or clinical utility.

## Introduction

Spatial transcriptomics measures molecular variation together with tissue location, but incomplete sampling creates a distinction between measurements and reconstructed values. Reconstruction should therefore be evaluated under explicit observation masks and with molecular information withheld from unmeasured locations. [REF: spatial transcriptomics and missing-data evaluation]

Tissue images offer a complementary source of information. A useful benchmark asks whether a frozen image predictor is improved by sparse local molecular measurements and whether images add information after spatial interpolation already uses those measurements. Correlation alone is insufficient because predictions can have spatially coherent rankings while remaining poorly calibrated in absolute target units. [REF: image-to-expression benchmarks and calibration]

Coverage and geometry address different scientific questions. Distributed anchors surround many hidden locations; a contiguous missing region can require extrapolation across a molecular boundary. Nearest-anchor distance may describe part of this difficulty, but it does not capture tissue heterogeneity or establish a causal reconstruction horizon. [REF: spatial interpolation and spatial validation]

Cross-platform evaluation further changes the target distribution, capture support and image resolution. Development, inspected transportability datasets and genuinely held-out validation must remain distinguishable. Reusing a cohort for calibration-development can be informative without making it a fresh validation cohort. [REF: external validation and dataset shift]

Here we assemble a frozen sequence of breast tissue experiments around a conservative claim: sparse measured expression can improve hidden-location reconstruction, subject to coverage and geometry limitations. We retain negative graph, sampling and image-fusion findings, and report patient-level replication with exact numerical provenance.

## Results

### Distributed molecular anchors improved reconstruction over the frozen image baseline

On the development primary patient, at 20% distributed coverage, mean-over-mask median gene Pearson correlations were {val(D,'M0')} for H&E, {val(D,'M1')} for selected ST interpolation and {val(D,'M3')} for fixed fusion. Corresponding mean gene MSE values were {val(D,'M0','mean_mse')}, {val(D,'M1','mean_mse')} and {val(D,'M3','mean_mse')}. Thus the development fusion correlation increase did not imply lower error than interpolation. These are one-patient development findings, not an external population estimate (Figure 2; Tables 1–2; ledger entries keyed Development/0.2/V1).

Fresh Wu validation retained all original primary genes. At 20% V1, sparse ST achieved r/MSE/MAE of {trip(W,'M1')}, compared with {trip(W,'M0')} for H&E. Fixed fusion yielded {trip(W,'M3')} (Figure 5). C9 remains frozen as MIXED_FRESH_EXTERNAL_VALIDATION: the sparse-ST comparison was supportive, but affine harmonization could not provide an identifiable interpolation advantage and the public image inputs limited H&E inference.

The second replication used one prospectively selected first section for each of 22 Bassiouni patients, totaling 30,003 source-filtered tissue locations. No second-cohort patient was used for training, calibration or selection. At 20% V1, sparse ST achieved {trip(C,'M1')}, versus {trip(C,'M0')} for H&E and {trip(C,'M3')} for fixed fusion. All 22 patients met the prespecified sparse-ST criterion at all three budgets: a Pearson gain of at least 0.02, relative MSE reduction of at least 2%, and non-worsening MAE. The replication is PASS_LIMITED_SECOND_REPLICATION because public H&E inputs contained only 24–29 native pixels across the unchanged 112-µm patch (Figure 6; Tables 1–2,4).

### Coverage benefits were strongest for distributed measurements

Across Bassiouni patients, median patient ST correlations at 10%, 20% and 40% V1 coverage were {val(C,'M1',budget=.1)}, {val(C,'M1')} and {val(C,'M1',budget=.4)}; MSE values were {val(C,'M1','mean_mse',.1)}, {val(C,'M1','mean_mse')} and {val(C,'M1','mean_mse',.4)}. Every patient improved both correlation and MSE from 10% to 40%. The Wu coverage-response was directionally consistent (Figures 5–7; Table 3). These tested budgets do not establish an optimal acquisition level or a monetary saving.

### Contiguous extrapolation and distance exposed reliability limits

At 20% coverage, Bassiouni median patient ST correlation fell from {val(C,'M1')} under V1 to {val(C,'M1',reg='V2')} under V2. MSE rose from {val(C,'M1','mean_mse')} to {val(C,'M1','mean_mse',reg='V2')}; all 22 patients worsened in correlation, MSE and MAE. The median paired relative MSE increase was 52.8%. This is a prespecified protocol contrast: the scoring region and frozen interpolation family both differ, so it cannot isolate geometry as a causal factor (Figures 3,5–6; Table 3).

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
'''
 (M/'MANUSCRIPT_DRAFT.md').write_text(draft)
 # Every numeric summary cell has a machine-readable ledger; method numbers point to frozen config.
 ld=pd.DataFrame(ledger);ld.to_csv(M/'RESULTS_LEDGER.csv',index=False);(M/'RESULTS_LEDGER.md').write_text('# Results ledger\n\nCSV row numbers are one-based physical records including the header (first data row=2). JSON references name exact fields. Source SHA-256 values identify the exact bytes used. Table2 cells point to explicit patient rows in CROSS_COHORT_SUMMARY.csv, which point to original checkpoint rows. All manuscript scientific numbers come from these entries; formatting/figure numbering and bibliographic years are not experimental results. Values are not rounded in the CSV.\n\n'+ld.to_markdown(index=False)+'\n')
 (M/'CITATION_NEEDS.md').write_text('''# Citation needs and author checks

This is a citation worklist, not a completed literature review. No reference metadata have been invented. Replace [REF] placeholders after human verification.

- Spatial transcriptomics measurement principles, sparse reconstruction and spatially blocked validation.
- Image-to-expression benchmarks; correlation versus absolute-error calibration.
- External validation, dataset shift and patient-level statistical units.
- Xenium development data: source studies/pages in reports/C0/C0_DATASET_STATISTICS.csv; HEST patch/physical-scale provenance in the C0 manifest.
- HER2ST: source DOI https://doi.org/10.1038/s41467-021-26271-2 and authors' deposit https://zenodo.org/records/4751624 . Verify full citation, licensing and acknowledgments.
- Wu local Visium: https://doi.org/10.1038/s41588-021-00911-1 and https://zenodo.org/records/4739739 . Distinguish local patients from the publication's separate HER2ST reanalysis.
- Bassiouni second replication: https://doi.org/10.1158/0008-5472.CAN-22-2682 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616 . Verify author list/title and source-study reference/validation roles; C10 uses all22patients as fresh evaluation.
- ResNet18, ImageNet, torchvision transforms and PyTorch; cite exact implementations/versions in ENVIRONMENT.json.
- Phikon and graph architectures only for historical supplementary comparisons, not the primary C10 method.
- Nearest-neighbor interpolation, inverse-distance weighting and farthest-point sampling as appropriate; do not claim algorithmic novelty without literature review.
- Authors must verify ethical/data-use statements, original consent/approvals, dataset licenses, funding, competing interests and acknowledgments. No clinical validation or prospective collection is claimed here.
''')
 # Human-readable replication report shares exact table values with the draft.
 primary=t2[(t2.cohort==C)&(t2.coverage==.2)&(t2.regime=='V1')][['method','median_pearson','mean_mse','mean_mae','patient_support_vs_M0_all_metrics']]
 report=f'''# C10 second replication and manuscript package

**{gate['status_pending_technical_audit']}**. Recommendation: **{gate['recommendation']}**.

The new GSE210616 study contributes22 fresh TNBC patients and30,003 tissue locations, one first section perpatient. All50primary genes are present. NoC10patient was used for fitting/calibration/selection. Source image resolution is limited to24–29nativepixels per112µmpatch, so this is limited replication rather than unrestricted validation of the image baseline.

## Primary20%V1 result

{primary.to_markdown(index=False)}

Values are median patient summaries after averaging10random seeds. M0=H&E,M1=ST,M3=fixedfusion. All22patients meet the meaningful sparse-ST joint criterion atall3budgets; all22improve r/MSE from10%to40%. Fusion improves20%V1Pearson in3/22but MSE/MAE in0/22. Median paired fusionΔr=−0.010830; relativeMSEharm=10.7078%.

## Coverage, geometry and distance

ST r/MSE/MAE at10%: {trip(C,'M1',.1)};20%: {trip(C,'M1')};40%: {trip(C,'M1',.4)}.

At20%V2: {trip(C,'M1',reg='V2')}. Everypatient worsens r/MSE/MAE versusV1. Median pairedΔr=−0.266923; relativeMSEincrease52.8206%. The V1/V2comparison changes scoring region and frozen interpolationfamily, not onlygeometry. C5block is a differentdevelopmentprotocol.

Distance replication is PARTLY: farthestoccupied-binMSE exceedsnearest in15/22V1and20/22V2; strictmonotonicity only6/22and13/22. Sparsebins and varyingpatientcomposition limitpooledcurves; nofixedhorizonisclaimed. Completecounts/correlations/errors areinC10B_DISTANCE_TO_ANCHOR.csv.

The frozen secondary spatialFPScomparison is in MANUSCRIPT/data/C10_FPS_effects.csv. It is descriptive and does not reopen samplingoptimization. Anchorassociationpermutations degradeM1/M3; imagepermutationhasmuchsmallerchanges. Allmatchedcontrols areinC10B_NEGATIVE_CONTROLS.csv.

## Integrity and scope

Independent audit: {audit['replayed_method_predictions']} exact methodpredictionreplays, {audit['replayed_metric_arrays']} exact pergenemetricarrays, {audit['hidden_corruption_checks']} recordedhiddenNaN/randominvariancecomparisons, {audit['label_free_mask_replays']} label-freemaskreplays. Sourceimages,rawcounts,models,scalers,barcodesandpriorregistryprefixverified. Protocol/source/image lock and masklock precede targetaccess. AffineM1/M2equivalence checked, duplicateM2rowsomitted. NoC0–C9experimentisrerunorretuned; parentverificationisread-only.

C10A_DATASET_QUALIFICATION.md documents metadata-only selection and alternatives screened. C10B/logs contains rawdownloads,qualification,locks,runrecords,hardware and audits. FrozenC9calibration is deliberately transported for secondaryM3, not fitted onthiscohort. H&E limitations and negativefusionfindings are retained.

## Manuscript package

Sevenmainfigures andtensupplementaryfigures, each300dpiPNGandSVG; exactsourceCSVs andperfigureplotentrypoints. Fourmain andeightsupplementarytables withCSV/Markdown (largepergeneCSVcompressed). A roughdraft,completefigurelegends,citationworklist,sourcehashindex,resultsledger andseven-claim evidence matrix are underMANUSCRIPT andC10root. Allgraphicsaredescriptive; nofabricatedpvalues. Forty-geneHER2ST andfifty-geneprimaryanalyses remainseparate. See C10_OVERNIGHT_SUMMARY.txt for the reviewhandoff.

Stop afterC10. Humanrevisionandcitationcheckingarethenextstep; noautomaticC11orpublication.
'''
 (OUT/'C10_REPORT.md').write_text(report);(B/'C10B_REPORT.md').write_text(report.replace('# C10 second replication and manuscript package','# C10B — locked second replication',1))
 # Required concise summary; final result commit is filled once the package is committed.
 summary=f'''C10 OVERNIGHT SUMMARY

1. SECOND REPLICATION
Dataset: GSE210616, Bassiouni et al. (2023)
Status: PASS_LIMITED_SECOND_REPLICATION
Patients: 22; one prospectively selected first section each
Evaluation patients: 22; no C10 calibration or training
Platform: Visium fresh frozen, primary TNBC
Primary gene overlap: 50/50
Evaluation locations: 30,003
Image limitation: 24–29 native pixels per 112-µm patch

2. PRIMARY RESULT
At 20% V1 (median patient Pearson / mean MSE / mean MAE):
H&E only: {trip(C,'M0')}
Sparse-ST interpolation: {trip(C,'M1')}
H&E + ST: {trip(C,'M3')}
Patient support: 22/22 meet joint meaningful ST benefit at all three budgets
Fusion support: 3/22 improve Pearson; 0/22 improve MSE or MAE

3. COVERAGE
10%: ST {trip(C,'M1',.1)}
20%: ST {trip(C,'M1')}
40%: ST {trip(C,'M1',.4)}
Did coverage-response replicate? YES
Support: 22/22 improve Pearson and MSE from 10% to 40% V1

4. CONTIGUOUS EXTRAPOLATION
V1: ST {trip(C,'M1')}
V2: ST {trip(C,'M1',reg='V2')}
Did V2 remain harder? YES
Support: 22/22 worsen Pearson, MSE and MAE at 20%
Median paired Pearson change: -0.266923
Median paired relative MSE increase: 52.82%
Contrast includes scoring geometry and frozen interpolation family

5. DISTANCE
Did error worsen with anchor distance? PARTLY
Farthest-bin MSE worse: 15/22 V1; 20/22 V2
Strictly increasing occupied-bin MSE: 6/22 V1; 13/22 V2
Sparse bins and patient composition preclude a universal distance cutoff

6. CROSS-COHORT CONCLUSION
Development: strong distributed sparse-ST benefit; contiguous limitations
HER2ST: supportive transportability evidence; calibration/fusion mixed; 40 genes
Wu: fresh sparse-ST support in 3 held-out patients; image inputs limited
Second replication: supportive in 22 fresh patients; image inputs limited
Core sparse-ST claim: STRONG
Scope: targeted-panel distributed reconstruction under frozen protocols
H&E incremental claim: NOT SUPPORTED
No general claim that morphology can never help

7. MANUSCRIPT PACKAGE
Draft: C10/MANUSCRIPT/MANUSCRIPT_DRAFT.md
Main figures: 7; C10/MANUSCRIPT/figures/; 300-dpi PNG + SVG
Supplementary figures: 10; C10/MANUSCRIPT/supp_figures/
Main tables: 4; C10/MANUSCRIPT/tables/
Supplementary tables: 8; C10/MANUSCRIPT/supp_tables/
Results ledger: C10/MANUSCRIPT/RESULTS_LEDGER.md (+ CSV)
Claim-evidence matrix: C10/CLAIM_EVIDENCE_MATRIX.csv
Citation-needs file: C10/MANUSCRIPT/CITATION_NEEDS.md
Figure legends: C10/MANUSCRIPT/FIGURE_LEGENDS.md

8. PUBLICATION READINESS
Ready to begin human manuscript writing? YES
Recommendation: READY_WITH_LIMITED_SECOND_REPLICATION
Remaining major evidence gap: native-resolution external H&E for morphology claims
That gap does not block writing the limited sparse-ST manuscript
Technical integrity: 8,766 exact prediction replays; 35,064 metric-array replays
Hidden-label corruption checks: 17,532; all invariant
Mask trajectories independently replayed: 968
C0–C9: preserved; read-only chained verification
Human work: citations, prose, final conceptual artwork and submission statements

9. FULL OUTPUT
C10 report: C10/C10_REPORT.md
Manuscript directory: C10/MANUSCRIPT/
Git commit: RESULTS_COMMIT_PENDING
Execution stopped after C10
'''
 (OUT/'C10_OVERNIGHT_SUMMARY.txt').write_text(summary)
 print('WROTE DRAFT',len(draft.split()),'words;',len(ledger),'ledgerentries',flush=True)
if __name__=='__main__':main()
