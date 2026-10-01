"""Write C9 evidence and immediate summary from locked evaluation outputs."""
from contract import *

def md(df):return df.to_markdown(index=False,floatfmt='.4f')
def main():
 p=pd.read_csv(B/'C9B_PATIENT_METRICS.csv');e=pd.read_csv(B/'C9B_PATIENT_EFFECTS.csv');g=json.loads((B/'logs/GATE_RESULT.json').read_text());audit=json.loads((B/'logs/INTEGRITY_AUDIT.json').read_text());assert audit['status']=='PASS'
 r=p[(p.strategy=='random')&(p.regime=='V1')];v=p[(p.strategy=='random')&(p.regime=='V2')];headline=r[r.budget==.2].groupby('method')[['median_pearson','mean_mse','mean_mae']].median();vhead=v[v.budget==.2].groupby('method')[['median_pearson','mean_mse','mean_mae']].median();coverage=r.groupby(['budget','method'])[['median_pearson','mean_mse','mean_mae']].median().reset_index();coverage.to_csv(B/'C9B_COHORT_MEDIANS.csv',index=False)
 a=r[(r.method=='M2')&(r.budget==.2)].set_index('patient_id');b=v[(v.method=='M2')&(v.budget==.2)].set_index('patient_id');degradation=pd.DataFrame({'delta_pearson_V2_minus_V1':b.median_pearson-a.median_pearson,'mse_increase_V2_minus_V1':b.mean_mse-a.mean_mse,'relative_mse_increase':b.mean_mse/a.mean_mse-1,'mae_increase':b.mean_mae-a.mean_mae});degradation.to_csv(B/'C9B_V2_DEGRADATION.csv')
 gene=pd.read_csv(B/'C9B_GENE_METRICS.csv');clean=gene[gene.control.isna()];gp=clean.groupby(['patient_id','regime','strategy','budget','method','gene','geometry'])[['pearson','spearman','mse','mae']].mean().groupby(['patient_id','regime','strategy','budget','method','gene']).mean().reset_index();gp.to_csv(B/'C9B_GENE_PATIENT_METRICS.csv',index=False)
 ga=[]
 for key,grp in gp.groupby(['patient_id','regime','strategy','budget','gene']):
  d=grp.set_index('method')
  for x,y,name in [('M2','M0','sparse_ST'),('M2','M1','harmonization'),('M3','M2','image')]:ga.append({**dict(zip(['patient_id','regime','strategy','budget','gene'],key)),'comparison':name,'delta_pearson':d.loc[x,'pearson']-d.loc[y,'pearson'],'mse_reduction':d.loc[y,'mse']-d.loc[x,'mse'],'mae_reduction':d.loc[y,'mae']-d.loc[x,'mae']})
 ga=pd.DataFrame(ga);ga.to_csv(B/'C9B_GENE_PATIENT_GAINS.csv',index=False);gs=ga[(ga.strategy=='random')&(ga.regime=='V1')&(ga.budget==.2)].copy();gs['joint_positive']=(gs.delta_pearson>1e-12)&(gs.mse_reduction>1e-12)&(gs.mae_reduction>1e-12);gene_support=gs.groupby(['comparison','gene']).agg(patients_joint_positive=('joint_positive','sum'),median_delta_pearson=('delta_pearson','median'),median_mse_reduction=('mse_reduction','median')).reset_index();gene_support.to_csv(B/'C9B_GENE_GAIN_SUMMARY.csv',index=False)
 nc=pd.read_csv(B/'C9B_NEGATIVE_CONTROLS.csv');ncs=nc.groupby(['regime','control','method'])[['delta_pearson','mse_reduction','mae_reduction']].mean().reset_index();ncs.to_csv(B/'C9B_NEGATIVE_CONTROL_SUMMARY.csv',index=False)
 effects=e[(e.strategy=='random')&(e.regime=='V1')&(e.budget==.2)];st=effects[effects.comparison=='sparse_ST'];im=effects[effects.comparison=='image'];q=json.loads((A/'logs/METADATA_QUALIFICATION.json').read_text());status={'checkpoint':'C9','C9A':'LIMITED_FRESH_EXTERNAL_DATASET','C9B':g['C9B_status'],'C9C':'NOT_RUN_NO_SECOND_COHORT_QUALIFIED','overall':g['C9B_status'],'core_sparse_ST_validated_with_limitations':g['core_sparse_support'],'target_harmonization_replication':'NOT_IDENTIFIABLE_M1_M2_AFFINE_INVARIANCE','H&E_fusion_validated':g['fusion_secondary_pass'],'recommendation':g['recommendation'],'patients':SAMPLES,'calibration_patients':CAL,'evaluation_patients':EVAL,'primary_gene_overlap':'50/50','image_resolution_limitation':'15/25nativepixels per112um patch; no source fullresWSI supplied','technical_audit':audit,'C0_C8_unchanged':True,'parent_verification_log':'C9/C9B/logs/PARENT_VERIFICATION.txt','spec_sha256':sha(A/'C9A_SPEC_LOCK.json'),'no_retuning':True,'spearman_roundoff_tie_limitation_disclosed':True,'stopped_after_C9':True};dump(OUT/'C9_STATUS.json',status)
 h=headline
 summary=f'''C9 IMMEDIATE SUMMARY

1. STATUS
C9A: LIMITED_FRESH_EXTERNAL_DATASET
C9B: {g['C9B_status']}
C9C: NOT_RUN; no second fresh cohort qualified.
Overall: {g['C9B_status']}; technical audits PASS.

2. FRESH EXTERNAL COHORT
Dataset/study: Wu et al. 2021 local Visium cohort; Zenodo 4739739.
Patients: 6 independent patients; 3 held out (8,481 locations).
Calibration patients: 1142243F, CID4465, CID4535.
Evaluation patients: 1160920F, CID4290, CID44971.
Platform: Visium; supplied H&E has only 15/25 pixels per 112 um patch.
Primary gene overlap: 50/50; original frozen order.
Why it is fresh/independent: Distinct study/patients/assets; no C0-C8 use.

3. PRIMARY SPARSE-ST RESULT
Coverage: 10%, 20%, 40%; headline 20% V1, medians across patients.
H&E only: r={h.loc['M0','median_pearson']:.3f}; MSE={h.loc['M0','mean_mse']:.4f}; MAE={h.loc['M0','mean_mae']:.4f}.
Raw ST interpolation: r={h.loc['M1','median_pearson']:.3f}; MSE={h.loc['M1','mean_mse']:.4f}; MAE={h.loc['M1','mean_mae']:.4f}.
Harmonized ST interpolation: Same as raw interpolation in common units.
H&E + harmonized ST: r={h.loc['M3','median_pearson']:.3f}; MSE={h.loc['M3','mean_mse']:.4f}; MAE={h.loc['M3','mean_mae']:.4f}.
Patient-level support: All 3 patients beat H&E on correlation/MSE/MAE at all 3 budgets.

4. TARGET HARMONIZATION
Did C8 harmonization replicate? NO; benefit is unidentifiable here.
Pearson effect: Zero; affine harmonization commutes with interpolation.
MSE effect: Zero in common units; no native-scale improvement claimed.

5. H&E INCREMENTAL VALUE
Did H&E improve beyond harmonized ST interpolation? NO.
Patient count improved: 1/3 Pearson; 0/3 MSE or MAE at 20% V1.
Main effect size: Median paired delta-r={im.delta_pearson.median():+.4f}; MSE {(-im.relative_mse_reduction.median())*100:.2f}% worse.

6. CONTIGUOUS EXTRAPOLATION
Result: M2 median r={vhead.loc['M2','median_pearson']:.3f}, MSE={vhead.loc['M2','mean_mse']:.4f} at 20% V2.
V1 vs V2 degradation: All 3 patients worsen; paired median delta-r={degradation.delta_pearson_V2_minus_V1.median():+.3f}.
Main limitation: Weak contiguous extrapolation; sparse-ST/H&E comparison also resolution-limited.

7. DISTANCE / COVERAGE
Coverage-response replicated? YES; correlation rises and MSE falls in all 3 patients.
Distance-to-anchor behavior replicated? YES in V1; V2 trends nonmonotonic.

8. PUBLICATION READINESS
Is the core sparse-ST reconstruction claim externally validated? YES, limited to this 3-patient test.
Is the H&E-fusion claim externally validated? NO.

9. RECOMMENDATION
{g['recommendation']}

10. FULL REPORT
Path: C9/C9_REPORT.md (repository); C0-C8 hashes unchanged.
Git commit: RESULTS_COMMIT_PENDING
'''
 result_file=OUT/'C9_RESULTS_COMMIT.txt'
 if result_file.exists():
  ref=result_file.read_text().strip();summary=summary.replace('RESULTS_COMMIT_PENDING',f'[{ref}](https://github.com/HySonLab/GraphPathology/commit/{ref})').replace('Path: C9/C9_REPORT.md (repository); C0-C8 hashes unchanged.',f'Path: [Full C9 report](https://github.com/HySonLab/GraphPathology/blob/{ref}/C9/C9_REPORT.md); C0-C8 unchanged.')
 (OUT/'C9_IMMEDIATE_SUMMARY.txt').write_text(summary)
 report=f'''# C9: fresh external validation

**{g['C9B_status']}**, with **positive but limited external evidence for sparse-ST interpolation**, no identifiable M2-M1 harmonization gain, and no reproducible H&E-fusion benefit. Recommendation: **{g['recommendation']}**. No method was retuned; C9C and manuscript generation were not performed.

The six-patient [Wu local Visium cohort](https://zenodo.org/records/4739739) is independent of development and HER2ST. Three patients calibrated T1 and three were held out throughout fitting. Qualification is **LIMITED_FRESH_EXTERNAL_DATASET** because the supplied public H&E has only 15/25 native pixels per 112 um patch. The primary panel contains all 50 original genes. The optional 278/280 intersection is metadata only. See [qualification](C9A/C9A_REPORT.md), [locked protocol](configs/protocol.json), and [specification hash](C9A/C9A_SPEC_LOCK.sha256).

## Primary held-out result

Ten random seeds per patient/coverage, with 10/20/40% measured locations. Tables average seeds within geometry, then geometries within patient; headline values are patient medians, not pooled spots. Correlations exclude undefined constant-gene values; threshold fractions count such genes as not passing. Errors use identical C0 standardized log1p-count units across methods. Pearson gains below are paired by patient; subtracting cohort medians is not substituted for the paired effect.

At 20% V1, all three patients beat frozen H&E in correlation, MSE and MAE; each does so at all three budgets. Median paired Pearson gain is **{st.delta_pearson.median():+.4f}** and median paired relative MSE reduction is **{st.relative_mse_reduction.median()*100:.2f}%**. These very large error reductions are partly relative to a badly transported H&E baseline; they do not establish accurate absolute molecular abundance across assays. The image-resolution and count-distribution shifts limit how strongly M2-M0 can be interpreted as an algorithmic advantage.

{md(headline.reset_index())}

All held-out patients at 20% V1:

{md(r[r.budget==.2][['patient_id','method','median_pearson','mean_spearman','mean_mse','mean_mae']])}

## What harmonization can and cannot establish

The prospective lock explicitly records that T1 is affine per gene. For weights summing to one, inverse_transform(sum(w*transform(y))) = sum(w*y). Consequently M1/M2 are equal in all 414 cases, with maximum absolute difference **{g['max_M1_M2_prediction_difference']:.3g}**, entirely numerical roundoff. Their Pearson and common-unit error effects are zero. **Raw Spearman values can differ** because machine-precision perturbations break tied ranks differently; the table above preserves the locked raw metric. These rank differences are numerical artifacts, not harmonization gains. `C9B_NUMERICAL_RANK_SENSITIVITY.csv` quantifies them without replacing predictions or metrics. This is not replication of C8's target-calibrated image/fusion benefit. C9 has no matched unharmonized-fusion arm, so M3-M2 cannot isolate that calibration effect either.

Full PASS required M2 to improve upon or be more stable than M1. That clause is scientifically unidentifiable with these specified methods, and was not waived after observing outcomes. Overall **MIXED** therefore coexists with positive core sparse-ST evidence. Native z-score errors are not reported as if they were gains in common-unit errors.

## H&E fusion

M3 uses exactly 0.2 calibrated image plus 0.8 harmonized ST, with no weight fitting. At 20% V1 only 1/3 patients improves Pearson, while 0/3 improves MSE or MAE. The median paired Pearson effect is **{im.delta_pearson.median():+.4f}**, median relative MSE increase **{-im.relative_mse_reduction.median()*100:.2f}%**. Image fusion fails the predeclared secondary criterion. The absence of benefit applies to these frozen inputs and this cohort; it does not establish that native high-resolution histology never helps.

{md(effects[effects.comparison=='image'][['patient_id','delta_pearson','relative_mse_reduction','mae_reduction']])}

## Coverage, geometry and distance

{md(coverage)}

M2 correlation increases and MSE decreases from 10% to 40% in all 3 patients. Spatial FPS is a separately reported, label-free secondary policy, not an optimized alternative. The regular-policy 20% comparison is below; it was not selected using outcomes. All regular-policy budgets remain in `C9B/C9B_PATIENT_METRICS.csv`.

{md(p[(p.strategy=='spatial_FPS')&(p.budget==.2)][['patient_id','regime','method','median_pearson','mean_mse','mean_mae']])}

At 20% V2, M2 median patient Pearson falls to **{vhead.loc['M2','median_pearson']:.4f}**, with MSE **{vhead.loc['M2','mean_mse']:.4f}** and MAE **{vhead.loc['M2','mean_mae']:.4f}**. Every patient worsens versus V1. At CID4290, V2 M2 Pearson is below M0 despite much lower error; this is disclosed rather than hidden. V2 pools no spots across patients, and all three predeclared geometry results are available in `C9B_GEOMETRY_METRICS.csv`. V2 differs both in missing-region shape and the frozen interpolation family (mean16 vs inverse-distance16 in V1), so this is protocol-level degradation, not an isolated causal estimate of mask geometry.

{md(v[v.budget==.2][['patient_id','method','median_pearson','mean_mse','mean_mae']])}

Paired V2-minus-V1 M2 degradation:

{md(degradation.reset_index())}

M2 MSE rises across occupied distance bins in each patient's 20% V1 curve. V2 broadly worsens with distance but is nonmonotonic, with gene/tissue-composition and geometric differences. Empty bins are omitted; binwise means average available runs/geometries and do not imply independent biological replicates. No unique extrapolation horizon or spot-level p-value is claimed.

## Controls and integrity

Expression-anchor and position-assignment permutations reduce M2 correlation and worsen error on the predeclared 1160920F / 20% / seed 0 case in both regimes, across 3 permutation seeds. NC1/NC2 are related assignment disruptions, not independent biological evidence. Image permutation slightly improves M3 here, consistent with weak/misaligned image contribution rather than supportive morphology dependence. M0 remains unchanged under molecular controls and M1/M2 unchanged under image permutation.

{md(ncs[ncs.method.isin(['M2','M3'])])}

Technical audit **PASS**: 132 mask trajectories replayed, 1,656 method predictions replayed exactly, 6,624 per-gene metric-array replays, 3,312 hidden-label NaN/random invariance checks, exact calibration-scaler rebuild using allowed patients only, explicit forbidden evaluation-patient calibration access rejected, and source/model/image hash checks. All 414 run records contain hashes, barcodes, masks, calibration IDs, runtime and peak RAM. All C0-C8 hashes and prior registry bytes are unchanged. The registry contains 132 mask events, 1 calibration event and 414 evaluation events.

## Gene-specific evidence and transportability

All per-gene Pearson/Spearman/MSE/MAE values, patient-level gains, and all gene support counts are saved. At 20% V1, **{int(((gene_support.comparison=='sparse_ST')&(gene_support.patients_joint_positive>=2)).sum())}/50 genes** improve Pearson, MSE, MAE over M0 in at least 2/3 patients; **{int(((gene_support.comparison=='image')&(gene_support.patients_joint_positive>=2)).sum())}/50** do so for M3 over M2. These are descriptive results, not a newly selected panel or multiple-testing-adjusted biological discovery.

`C9_TRANSPORTABILITY_SUMMARY.csv` uses the same C7/C8 common 40 genes and common error units for development, HER2ST and fresh C9. The C8A T0 rows replay C7 exactly; C8A T1 and C8B rows are transportability development, not fresh independent validation. The main C9 results use 50 genes and must not be directly compared to 40-gene median results. Absolute errors remain sensitive to platform capture/depth despite common mathematical units; no clinical utility, cost saving or whole-transcriptome claim is made.

## Recommendation and stopping point

**{g['recommendation']}**. A second fresh cohort with better-resolution paired H&E would strengthen the core sparse-ST claim and clarify its limits. The frozen four-method C9 comparison cannot prove an affine-harmonization benefit over identical raw interpolation. Review that estimand before designing another harmonization test; do not reinterpret C9 or retune on its held-out patients. No second cohort was readily qualified after C9B, so optional C9C was not run. Stop for review; do not generate the manuscript automatically.

## Reproducibility and figures

All source/configuration/masks/scalers/metrics/logs/reports/figures are under C9. Large public raw archives, image caches and predictions remain local under C9 and are SHA-256 covered, with acquisition/source code committed. `README.md` lists stage commands. Parent verification and final manifest cover the complete checkpoint.

'''
 # Normalize a few compact labels in prose without changing scientific files.
 report=report.replace('H ER2ST','HER2ST')
 for pth in sorted((B/'figures').glob('[0-9][0-9]_*.png')):report+=f'![{pth.stem}](C9B/figures/{pth.name})\n\n'
 (OUT/'C9_REPORT.md').write_text(report.rstrip()+'\n')
 (B/'C9B_REPORT.md').write_text('# C9B: locked held-out validation\n\n'+report.split('## Primary held-out result',1)[1].split('## Reproducibility and figures',1)[0].rstrip()+'\n')
 print(summary,end='')
if __name__=='__main__':main()
