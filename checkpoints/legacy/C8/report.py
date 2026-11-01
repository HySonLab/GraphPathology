"""Render C8 scientific reports without changing fitted objects or gate logic."""
from contract import *
def table(df):return df.to_markdown(index=False,floatfmt='.5f')
def main():
 sa=json.loads((A/'C8A_STATUS.json').read_text());sb=json.loads((B/'C8B_STATUS.json').read_text());audit=json.loads((OUT/'logs/INTEGRITY_AUDIT.json').read_text());assert audit['status']=='PASS';assert sa['status']=='PASS_TARGET_HARMONIZATION' and sb['status']=='MIXED_FUSION_CALIBRATION';assert not (OUT/'C8C').exists()
 status={'C8A_status':sa['status'],'C8B_status':sb['status'],'C8C_status':'NOT_RUN_C8B_DID_NOT_PASS','C8_status':'MIXED_CROSS_PLATFORM_CALIBRATION','recommendation':'NEED_METHOD_REVIEW','technical_audit':audit,'C0_C7_unchanged':True,'C7_zero_shot_preserved_exactly':True,'HER2ST_role':'method-development transportability; not final external validation','patients':8,'outer_folds':4,'primary_genes':40,'secondary_panel':'not run, optional','A_gate':sa,'B_gate':sb,'source_spec_commit':'b75b764f725e9da394f8efd230bc6f589d38dcc7','source_spec_sha256':sha(A/'C8A_SPEC_LOCK.json'),'main_outer_evaluations':6720,'negative_controls':576,'inner_A_evaluations':3888,'target_fit_registry_events':12,'fusion_candidate_fit_events':192,'registry_events':11388,'no_encoder_or_graph_changes':True,'stopped_after_C8':True};dump(OUT/'C8_STATUS.json',status)
 a=pd.read_csv(A/'C8A_COHORT_MEDIANS.csv');b=pd.read_csv(B/'C8B_COHORT_MEDIANS.csv');p=pd.read_csv(B/'C8B_PATIENT_EFFECTS.csv');scales=pd.read_csv(A/'C8A_TARGET_SCALE_CONSISTENCY.csv');models=pd.read_csv(B/'C8B_FUSION_MODELS.csv');nc=pd.read_csv(B/'C8B_NEGATIVE_CONTROL_PATIENTS.csv');genes=pd.read_csv(B/'C8B_GENE_GAIN_SUMMARY.csv')
 models_display=[]
 for fold in cfg()['folds']:
  obj=json.loads((B/'calibration'/f'{fold}.json').read_text());models_display.append({'fold':fold,'held_out':','.join(obj['held_out_patients']),'target':'T1','F1_lambda':obj['family_models']['F1']['lambda'],'F2_weights':str(obj['family_models']['F2']['weights']),'selected_family':obj['selected_family']})
 nc_summary=nc.groupby(['model','regime','control'],as_index=False)[['Pearson_change','common_MSE_change']].median();nc_summary.to_csv(B/'C8B_NEGATIVE_CONTROL_SUMMARY.csv',index=False)
 gene_summary=genes.groupby('regime').apply(lambda g:pd.Series({'positive_median_Pearson_genes':int((g.pearson_gain>0).sum()),'positive_median_MSE_genes':int((g.common_mse_reduction>0).sum()),'both_positive_genes':int(((g.pearson_gain>0)&(g.common_mse_reduction>0)).sum())}),include_groups=False).reset_index()
 fractions=[c for c in pd.read_csv(B/'C8B_PANEL_METRICS.csv',nrows=0).columns if c.startswith('fraction_pearson')];rf=pd.read_csv(B/'C8B_PANEL_METRICS.csv');rf=rf[rf.control.isna()];rf.groupby(['fold','sample_id','model','regime','budget'],as_index=False)[fractions+['n_defined_pearson']].mean().to_csv(B/'C8B_PATIENT_CORRELATION_FRACTIONS.csv',index=False)
 report='''# C8 — Cross-platform target harmonization and fusion calibration

**Overall: MIXED_CROSS_PLATFORM_CALIBRATION. Recommendation: NEED_METHOD_REVIEW.** C8A passed target harmonization; C8B produced mixed fusion calibration and did not meet its frozen pass criteria. **C8C was not run**, as required by the stage gate. C0–C7 remain unchanged, including the exact C7 zero-shot predictions. No later checkpoint or fresh external validation was started.

The central result is a substantial repair of cross-platform scale error, with only modest and heterogeneous evidence that H&E adds value beyond molecular interpolation. At 20% V1 coverage, training-selected calibrated fusion has median patient Pearson **0.21505**, common-unit MSE **0.07004**, and MAE **0.18709**, compared with ST interpolation **0.19431 / 0.06797 / 0.17298**. The median paired Pearson gain is **+0.00856**, and the median paired relative MSE reduction is **0.995%**, but only **3/8** patients jointly improve Pearson with at most 1% MSE harm. **MAE worsens in all eight patients.** These paired effects, not differences of independently computed cohort medians, determine the scientific conclusion.

## Dataset roles, folds and unchanged inputs

HER2ST is now **method-development transportability data**, not a pristine final external cohort. Its outcomes were examined in C7, and the C8 plan explicitly permits calibration development on it. C8's patient-disjoint outer cross-validation evaluates the predefined calibration procedure within this cohort; it cannot erase the earlier development knowledge or substitute for new external validation.

The eight patients and one section each (A1–H1, 3,455 locations) come from the [HER2ST study](https://doi.org/10.1038/s41467-021-26271-2) and [source deposit](https://zenodo.org/records/4751624). Source acquisition/provenance remains frozen under C7. The primary panel is exactly C7's common **40/50 genes**; the optional 226-gene panel was not run so this checkpoint could focus on the calibration question. No missing gene was replaced and no gene was selected using C8 performance.

Four outer folds were fixed by patient IDs before C8 fitting: A/B, C/D, E/F and G/H. Each has six training patients and two held-out patients. Within each six-patient training set, leave-one-patient-out records use scalers fit on the other five patients. The original C0 training patients **NCBI783 and NCBI785** supply source target moments; source validation/test patients TENX95/TENX99 do not fit a C8 transform. Only saved source-test predictions are re-summarized later for descriptive comparison.

Image embeddings, encoder weights, image heads, original C7 panel, interpolation rules and observation masks are read-only. There is no encoder fine-tuning, head refit, learned graph, multiresolution model or sampling search. Sampling is the predeclared **random C7 trajectories only**: ten outer seeds at 10/20/40% V1, three fixed V2 missing-region geometries at 20%. Inner fitting uses seeds 0/1/2 only. V2 excludes all anchors inside each fixed 40% region; every other unobserved label is also withheld from model inputs. V1 target membership changes with coverage; V2 targets remain identical across models for each geometry.

## What harmonization means here

Let `y` and frozen image output `p` be expressed in C7's original C0-standardized log-count units. For each gene, fit source center/scale `(a_s,b_s)` on the original source-training patients and destination center/scale `(a_t,b_t)` on the allowed HER2ST training patients. Then:

```text
image_native  = (p - a_s) / b_s
anchor_native = (y_observed - a_t) / b_t
prediction_common_C7 = a_t + b_t * prediction_native
```

This applies **separate assay-specific scales**, so it genuinely changes the image prediction after inversion. Applying one identical affine transform to targets and predictions would only relabel errors and could not fix the C7 mismatch. All transform/fusion selection and error gates therefore use **inverse-mapped errors in the original C7 units**, alongside Pearson. Native transform MSE/MAE are saved separately but must not be compared numerically across transforms.

- **T0:** exact C7 representation and predictions; identity source/destination maps.
- **T1:** per-gene mean/SD estimated with equal total weight per training patient in each domain. Scales below 1e-8 use 1. Source moments are recomputed from the original two source-training patients under this declared patient-weighting rule; no C0 file is edited.
- **T2:** patient-weighted pooled median/IQR in each domain, using a stable-sort left inverse empirical CDF. IQR below 1e-8 falls back to that domain's training SD, then 1 if needed. Destination fallback occurred for 15, 20, 12 and 13 of 40 genes in folds 0–3; source fallback was zero. The fallback is predefined, never selected from held-out performance.

These positive per-gene affine maps preserve within-sample gene ordering and correlation under simple coordinate conversion. They do not establish assay-equivalent absolute RNA abundance or eliminate between-patient variation. T3 within-sample ranks were omitted because computing full held-out ranks would use hidden expression; T4 was omitted to keep the candidate set small. No held-out patient estimates its own scaler, even from its observed anchors. Anchors are used only for the allowed reconstruction inputs.

## Selection and lock chronology

[C8A_SPEC_LOCK.json](C8A/C8A_SPEC_LOCK.json), SHA-256 `3d087940666ec8c0662fbd241d0d66f234f0a4401ba22c2364dcb501c477fcc0`, seals the complete source implementation, folds, transforms, fusion grids/forms, mask choices, selection rule, gates, original inputs and parent manifest. It was committed as `b75b764` before C8 fitting. All A fold scaler objects and training-only selections were sealed in [FIT_LOCK.json](C8A/FIT_LOCK.json), committed as `5a97385`, before outer A evaluation. After A passed, B calibration objects were sealed in [C8B_SPEC_LOCK.json](C8B/C8B_SPEC_LOCK.json), committed as `d29ed78`, before outer B evaluation. No locked scientific source was altered after external C8 predictions were evaluated.

For selection, average mask seeds/geometries within each inner patient/regime/budget, then the four conditions (V1 at three budgets; V2 at 20%) equally, then the six training patients equally. Maximize mean patient median-gene Pearson; among candidates within **0.005** of the best, choose the lowest **common-unit MSE**, with deterministic ties. A chooses a transform using fixed E2 fusion. B chooses parameters within each fusion family and then a family, including F0, using the same inner-record rule. The selected-per-fold procedure is the primary algorithm; no family was promoted because it looked better on outer patients.

## C8A: target harmonization passed

All four inner selections chose **T1**. With fixed C7 fusion weights (25% image in V1; 50% in V2), primary V1 at 20% improves median paired Pearson by **+0.01372** over T0 and reduces common-unit MSE by **79.99%** (median patient-relative reduction). Seven of eight patients meet the joint harmonization support criterion. Patient common-MSE IQR falls **0.06440 → 0.02519**. All three V1 budgets have nonnegative median Pearson gain and lower MSE. V2's median Pearson gain is **+0.00861**, and all eight patients pass its error safety bound.

20% coverage, medians across eight held-out patient summaries:

'''+table(a[a.budget==.2][['transform','model','regime','median_pearson','mean_mse','common_mse','common_mae']])+'''

Native `mean_mse` uses each transform's units; `common_mse` and `common_mae` use unchanged C7 units. T1 reduces H&E-only common MSE substantially without changing its gene-wise correlation, because scale calibration cannot create new spatial information in a frozen image predictor. Molecular interpolation remains effectively identical in common units under these affine transformations, apart from numerical precision.

The held-out target variance is not made uniform: at 20% V1, mean-gene native variance has median/range **0.07438 / 0.03672–0.17203** for T0, **0.63652 / 0.19308–2.59243** for T1, and **0.45019 / 0.12149–2.20531** for T2. T1 held-out mean levels range −0.5074 to +0.9851 in its native units. Training moment alignment therefore leaves substantial patient-specific distribution shift. Exact values are in [C8A_TARGET_SCALE_CONSISTENCY.csv](C8A/C8A_TARGET_SCALE_CONSISTENCY.csv).

The A gate required a non-T0 choice in ≥3 folds, preserved/improved Pearson with ≥10% median MSE reduction at 20%, joint support in ≥6/8 patients, non-increased error IQR, support at ≥2 V1 budgets, and non-catastrophic V2 behavior. All requirements passed. This gate concerns harmonization versus T0, not superiority over ST interpolation.

## C8B: fusion calibration is mixed

F0 is ST interpolation. F1 uses one global scalar image weight per fold, shared across regimes/budgets, from 0.0–1.0 in 0.1 steps. F2 uses nondecreasing image weights in nearest-anchor bins <300, 300–600 and ≥600 µm, with 35 ordered triples from {0,0.1,0.25,0.5,1}. F3 is a five-parameter sigmoid linear gate using log nearest/mean-kNN anchor distances, anchor count within 600 µm and mean gene anchor-disagreement variance. Feature scaling and gate fitting use only the inner out-of-fold training records. Its frozen optimizer minimizes equal-patient/equal-condition common-unit error normalized by the corresponding interpolation error, with fixed regularization. All four optimizers converged and were independently refit exactly.

Selected fold objects:

'''+table(pd.DataFrame(models_display))+'''

F2 selected **constant** weights in every fold, so no distance dependence was retained by that candidate family. The primary procedure selects F1 with λ=0.2 for A/B, C/D and E/F, and F3 for G/H. The complete F3 coefficients and feature normalizers are in the hashed calibration JSONs. F4/F5 were not run.

All candidate-family outer summaries are shown for transparency, but only `SELECTED` determines the gate:

'''+table(b[['regime','budget','model','median_pearson','common_mse','common_mae']])+'''

At 20% V1, selected fusion improves Pearson in **5/8** patients and MSE in **5/8**, but the improving sets differ: only D/G/H improve both. Patient E has **47.90% worse MSE**, and MAE worsens in **8/8**. The median paired relative MSE reduction (+0.995%) coexists with a larger separate cohort-median MSE (0.07004 versus 0.06797); these are distinct summaries, neither should replace the other.

Paired patient effects at 20%:

'''+table(p[p.budget==.2][['sample_id','regime','SELECTED_median_pearson','F0_median_pearson','SELECTED_common_mse','F0_common_mse','gain_pearson','relative_mse_reduction','mae_reduction']])+'''

Frozen V1 budget gates:

'''+table(pd.DataFrame(sb['budget_support']))+'''

PASS required median ΔPearson ≥0.01, nonnegative median relative MSE reduction and joint patient support ≥5/8 at at least two budgets; at most 25% of V1 patient–budget settings could have >2% MSE harm. **No budget meets all conditions**, and **9/24 settings (37.5%)** exceed the error-harm limit. At 40%, seven patients improve correlation but only three improve MSE; the median relative MSE change is adverse. The mixed fallback is satisfied because multiple patients have genuine paired correlation/error gains. The exact locked analyzer therefore returns **MIXED_FUSION_CALIBRATION**.

F3 considered across all folds has lower aggregate MSE than the selected procedure in some settings, but selecting it now would use held-out outcomes. It remains a secondary candidate comparison; no post-hoc replacement, grid expansion or tolerance change was made.

## Contiguous extrapolation and distance

At 20% V2, selected fusion has median Pearson **0.10855**, common MSE **0.08086**, MAE **0.21031**, versus interpolation **0.06379 / 0.07751 / 0.19623**. Median paired ΔPearson is **+0.02787** and median relative MSE reduction **0.833%**. Six patients improve correlation, four improve MSE, and only one improves MAE. This is weaker absolute correlation than V1 despite scale calibration. Patient E again has substantial harm (48.72% MSE increase).

By fixed missing-region geometry:

'''+table(pd.DataFrame(sb['geometry_support']))+'''

Only **1/3** geometries has both nonnegative median correlation and error improvement, below the required 2/3. C8C's additional V2 10/40% stress runs are therefore **not executed**. This is an intentional stage-gated stop, not missing execution of an authorized unconditional stage.

Distance summaries retain the predeclared C7 bins and equal-patient aggregation after weighting repeated mask incidences within each patient/bin. At 20% V1, median paired MSE reductions are approximately 0.000616 at 150–300 µm, 0.000049 at 300–600 µm and 0.001851 at 600–1200 µm. V2 is adverse at 150–300 µm (−0.000387), positive but small farther away. The V1 <150 µm bin has one patient and six repeated incidences; the farthest V2 bin has six patients. These descriptive, occupancy-confounded results do not validate a distance horizon or establish reliability gating as the source of gains.

## Gene-level evidence

At 20%, average matched-mask gene effects within patient first, then summarize across patients. All 40 genes are retained, with their original order and no outcome-based selection:

'''+table(gene_summary)+'''

Thus 27/40 V1 genes and 29/40 V2 genes have positive median paired correlation **and** common-MSE gains. These are descriptive, not multiplicity-controlled discoveries, and do not override adverse patient-level errors or the frozen gate. Undefined constant-gene correlations remain missing; panel fractions use all 40 genes in the denominator. White heatmap cells indicate undefined correlations, not zero gains. No pathway enrichment, cell-type inference or biological interpretation was performed.

## Negative controls

All eight held-out patients, 20% coverage, seed 0, V1 geometry 0 and V2 geometry 1; F1/F2/F3/selected methods; NC1/NC2/NC3; three control seeds: **576 runs**. NC1 randomly permutes the fitted location-weight vector, preserving its empirical distribution; NC2 permutes distance-derived features across locations; NC3 permutes image predictions while retaining ST and gate inputs.

NC1 is a matched empirical random-weight control. For a constant scalar it is necessarily an identity operation, and cannot test whether that scalar value is optimal. F1 and all selected F2 lookup tables are constant; their NC1/NC2 invariance is expected. Consequently these controls provide limited evidence about global scalar calibration. They must not be reported as a successful degradation test for a quantity that is unchanged.

Median patient changes (each patient first averaged across three control seeds; control minus clean):

'''+table(nc_summary)+'''

Image permutation reduces selected-method Pearson by **0.01925 in V1** and **0.04229 in V2**, with higher MSE, supporting a location-specific H&E signal. Confidence-gate controls are mixed: F3 distance permutation can improve correlation and slightly improve V2 MSE. The reliability features therefore do not have uniform negative-control support. These controls do not rescue the failed main gate.

## Domain-shift recovery and remaining limitations

For identical 40 genes and 20% V1, the saved C5 source-test fusion Pearson was **0.84091**; C7 zero-shot fusion was **0.25463** (eight-patient median), C8A harmonized fixed fusion **0.26187**, and C8B selected calibrated fusion **0.21505**. Corresponding common-unit MSE values are **0.27352, 0.34809, 0.07388, 0.07004**. Calibration repairs much of the count-scale error but does not recover source-level spatial correlation. Its more conservative image contribution reduces absolute median correlation relative to fixed C7/C8A fusion while producing much lower error. This is not a uniformly better predictor on every metric.

[Development transportability table](C8_DEVELOPMENT_TRANSPORTABILITY.csv) includes patient, transform, method, coverage, regime, Pearson/MSE/MAE and gains over ST/H&E. Source-development rows reuse saved C5 predictions; no parent model rerun. C5 has one primary test patient, whereas HER2ST has eight; C5 contiguous geometry differs from C7/C8 fixed regions. The contrast is descriptive, not a matched cohort experiment. Assay capture area/depth and gene sparsity still differ, so common mathematical units do not establish common measurement calibration.

This remains one small, already examined cohort with overlapping CV training sets and 40 retained genes. A and B selection reuse inner records but outer-patient labels are excluded throughout; reported inner optima are not unbiased test estimates. Neither additional folds nor thousands of spots create additional independent biological patients. No confidence interval treats spots as independent. The known HER2ST cohort cannot later serve as pristine final validation for a revised method. There are no clinical utility, cost savings or whole-transcriptome claims.

## Integrity and reproducibility

The independent audit passed **7,296 prediction-file checks**, including every barcode/mask/scaler/model hash and distance identity; **1,440 T0 predictions exactly match frozen C7**. It rebuilt all **84** outer/inner scaler objects while explicitly forbidding outer-held-out label access, replayed all eight transform/family selections, and independently refit all four F3 gates exactly. **1,248** stratified predictions and gene metrics replayed exactly; **2,496** extra hidden-label corruption comparisons passed, in addition to **14,592** checks recorded across all outer runs. Both stage-fit locks preceded their outer-evaluation access records. There are **11,388** appended registry events, with the prior prefix unchanged.

Counts: 3,888 inner A evaluations; 12 outer scaler-fit events; 4,320 outer A evaluations; 192 fusion candidate-fit events; 2,400 outer B evaluations (including the selected algorithm); 576 negative controls. No C8C outputs are created. Scaler/gating objects, locks, tables, per-run provenance, figures and compressed gene metrics are committed. Large predictions and inner-cache arrays remain local under C8 and are covered by the final manifest. The required plain gene CSV is present locally and has an identical decompressed gzip counterpart in Git. Read-only final verification checks the complete C0–C7 parent chain and C8 artifact hashes.

## Recommendation and stop

**NEED_METHOD_REVIEW.** The next decision should address the residual patient-specific error/MAE harm and whether H&E's modest contribution is reliable enough to justify another method revision. Any revision needs a separately declared training/calibration protocol and subsequent evaluation on a fresh cohort unused in C0–C8. This checkpoint does **not** authorize that next work. C8 is complete at a mixed result and stops for review.

## Figures

'''
 figures=[('C8A/figures/target_transforms.png','Frozen C7 versus target transforms'),('C8A/figures/patient_transforms.png','Target transforms for every patient'),('C8A/figures/target_coverage.png','Target harmonization coverage response'),('C8B/figures/fusion_strategies.png','Fusion candidates and interpolation versus contiguous regions'),('C8B/figures/coverage_response.png','Fusion coverage response; V2 only at permitted 20%'),('C8B/figures/patient_fusion_effects.png','Paired patient effects'),('C8B/figures/gene_fusion_effects.png','All 40 genes; white cells have undefined correlation'),('C8B/figures/fusion_distance_gain.png','Paired fusion gain by distance'),('C8B/figures/domain_gap.png','Source-to-HER2ST residual domain shift')]
 for path,title in figures:report+=f'### {title}\n\n![{title}]({path})\n\n'
 (OUT/'C8_REPORT.md').write_text(report.rstrip()+'\n')
 (B/'C8B_REPORT.md').write_text('# C8B — Fusion calibration\n\n**MIXED_FUSION_CALIBRATION**; technical audit **PASS**. C8C was not run.\n\nAll four folds used training-selected T1. F1 with image weight 0.2 was selected in folds 0–2; the five-parameter F3 gate was selected in fold 3. No method family was selected from outer patient outcomes.\n\nAt 20% V1, median paired ΔPearson is +0.00856 and median paired relative MSE reduction is 0.995%, but only 3/8 patients meet joint support. MSE harm >2% occurs in 9/24 primary patient-budget settings, and MAE worsens in all 8 at 20%. Only 1/3 V2 geometries supports both metrics. The frozen pass gate fails and its mixed fallback applies.\n\n'+table(b[['regime','budget','model','median_pearson','common_mse','common_mae']])+'\n\nSee [full C8 report](../C8_REPORT.md) for every patient, controls, all forty genes, units, fold protocol, locked gates and limitations. Common-unit errors are in the original C7 target units. Native transform errors are saved separately. [Audit](../logs/INTEGRITY_AUDIT.json) passed; parent artifacts and C7 predictions remain unchanged.\n')
 summary='''C8 IMMEDIATE SUMMARY

1. STATUS
C8A: PASS_TARGET_HARMONIZATION
C8B: MIXED_FUSION_CALIBRATION
C8C: NOT RUN — C8B did not pass.
Overall: MIXED_CROSS_PLATFORM_CALIBRATION; technical audit PASS.
C0–C7 unchanged; 1,440 reference predictions match frozen C7 exactly.

2. TARGET HARMONIZATION
Best transform: T1, source/destination training-patient gene-wise z-scores.
Why selected: chosen by inner patient-disjoint selection in all four folds.
Primary Pearson: 0.26187, fixed fusion at 20% V1 coverage.
Primary MSE: 0.07388 in unchanged C7 units.
Improvement vs C7 target: paired ΔPearson +0.01372; MSE reduction 79.99%.
Seven of eight patients meet the joint harmonization support rule.
No held-out patient fits a scaler; native-scale errors do not decide the gate.

3. FUSION CALIBRATION
Best fusion: training-selected F1 (image weight 0.2) in folds 0–2;
F3 confidence gate in fold 3. No outer-outcome method selection.
ST interpolation Pearson / MSE: 0.19431 / 0.06797 at 20% V1.
Calibrated fusion Pearson / MSE: 0.21505 / 0.07004 at 20% V1.
Patient improvement count: 5/8 Pearson, 5/8 MSE, only 3/8 jointly.
Median paired gain: ΔPearson +0.00856; relative MSE reduction 0.995%.
Cohort medians and medians of paired effects are different summaries.
Coverage levels with benefit: modest paired gains at 10% and 20%;
none of 10/20/40% meets every frozen budget-level pass condition.
MSE harm >2% occurs in 9/24 patient–budget settings; allowed maximum was 25%.
MAE worsens in all eight patients at 20% V1; patient E MSE worsens ~48%.

4. CONTIGUOUS EXTRAPOLATION
Result: 20% V2 fusion Pearson / MSE = 0.10855 / 0.08086.
ST interpolation: 0.06379 / 0.07751.
Median paired ΔPearson +0.02787; relative MSE reduction 0.833%.
Main limitation: only 1/3 fixed geometries supports both metrics.
C8C's additional stress runs were correctly stopped by the failed B gate.

5. DOMAIN-SHIFT RECOVERY
C7 frozen external result: fusion Pearson 0.25463; MSE 0.34809.
C8 calibrated result: Pearson 0.21505; MSE 0.07004 (20% V1).
Scale error improves substantially; spatial correlation remains far below source.
Did H&E add value beyond ST interpolation? PARTLY, not robustly enough to pass.
27/40 genes improve median paired Pearson and MSE in V1 at 20%.

6. LIMITATIONS
- HER2ST is development transportability data, not final external validation.
- Eight patients, four folds, one cohort, forty frozen genes.
- Residual patient-specific scale/variance shift and systematic MAE harm.
- Distance-gate selection collapsed to constant weights in every fold.
- Weight/distance permutations are null for constant gates; other controls mixed.
- 7,296 prediction files checked; 1,248 exact replays; leakage checks passed.

7. RECOMMENDATION
NEED_METHOD_REVIEW
Do not reuse HER2ST as pristine final validation for a revised method.
Stopped after C8; no fresh external cohort or later checkpoint started.

8. FULL REPORT
Path: /home/son/Research/GraphPathology/C8/C8_REPORT.md
Git commit: d29ed78 (sealed calibration objects); results commit recorded in history.
'''
 ref=OUT/'C8_RESULTS_COMMIT.txt'
 if ref.exists():summary=summary.replace('d29ed78 (sealed calibration objects); results commit recorded in history',ref.read_text().strip()+' (completed results snapshot)')
 assert 40<=len(summary.splitlines())<=80;(OUT/'C8_IMMEDIATE_SUMMARY.txt').write_text(summary);print('Reports written;',len(summary.splitlines()),'summary lines')
if __name__=='__main__':main()
