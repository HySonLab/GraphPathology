"""Generate the C1 checkpoint report directly from recorded metrics."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C1'

def fmt(values,sd=False,digits=3):
 a=np.asarray(values,dtype=float);a=a[np.isfinite(a)]
 if not len(a):return 'undefined'
 base=f'{a.mean():.{digits}f}'
 return base+(f' ± {a.std(ddof=1):.{digits}f}' if sd and len(a)>1 else '')

def table(headers,rows):
 return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(x) for x in row)+' |' for row in rows])

def main():
 p=pd.read_csv(OUT/'C1_METRICS_SUMMARY.csv');nc=pd.read_csv(OUT/'C1_NEGATIVE_CONTROLS.csv');lo=pd.read_csv(OUT/'C1_ROBUSTNESS_LOPO.csv')
 resources=pd.read_csv(OUT/'C1_RESOURCE_PROFILE.csv');rank=pd.read_csv(OUT/'C1_GENE_RANKINGS.csv')
 lock=json.loads((OUT/'SELECTION_LOCK.json').read_text());audit=json.loads((OUT/'logs/INTEGRITY_AUDIT.json').read_text())
 assert audit['status']=='PASS'
 def values(df,model,panel,col):return df[(df.model==model)&(df.panel==panel)][col]
 error_ok=all(values(p,m,panel,'mean_mse').max()<values(p,'B0',panel,'mean_mse').iloc[0] for m in ['B1','B2'] for panel in lock['panels'])
 controls_ok=all(values(p,m,'S2-50','median_pearson').min()>nc[nc.model==m].median_pearson.max() for m in ['B1','B2'])
 lopo_ok=all(lo[(lo.model=='B2')&(lo.panel=='S2-50')&(lo.sample_id==s)].median_pearson.mean()>0 for s in lo.sample_id.unique())
 subset_ok=all(values(p,'B2',panel,'fraction_pearson_gt_0.2').mean()>.5 for panel in ['S2-50','S2-250','S2-280'])
 # These descriptive checks substantiate the plan's qualitative criteria; they are not new statistical significance thresholds.
 assert error_ok and controls_ok and lopo_ok and subset_ok
 status='PASS_DEVELOPMENT'
 rows=[]
 for model in ['B0','B1','B2']:
  row=[model]
  row += [fmt(values(p,model,panel,'median_pearson'),model=='B2') for panel in ['S1','S2-50','S2-250','S2-280']]
  row += [fmt(values(p,model,panel,'coexpression_pearson'),model=='B2') for panel in ['S2-50','S2-250','S2-280']]
  r=resources[(resources.phase=='primary')&(resources.model==model)]
  row += [fmt(r.groupby('seed').train_seconds.sum(),model=='B2',4)+' s',fmt(r.groupby('seed').peak_rss_mb.max(),model=='B2',1)+' MiB']
  rows.append(row)
 primary_table=table(['Model','S1 Pearson','S2-50 median Pearson','S2-250 median Pearson','S2-280 median Pearson','S2-50 coexpression','S2-250 coexpression','S2-280 coexpression','Training time','Peak RSS'],rows)
 mse_rows=[]
 for panel in lock['panels']:
  baseline=values(p,'B0',panel,'mean_mse').iloc[0]
  mse_rows.append([panel,fmt([baseline]),fmt(values(p,'B1',panel,'mean_mse')),fmt(values(p,'B2',panel,'mean_mse'),True),f'{100*(1-values(p,"B2",panel,"mean_mse").mean()/baseline):.1f}%'])
 nc_rows=[]
 for model in ['B1','B2']:
  for phase in ['primary','NC1','NC2']:
   frame=p if phase=='primary' else nc[nc.phase==phase]
   nc_rows.append([model,phase,fmt(values(frame,model,'S2-50','median_pearson'),True),fmt(values(frame,model,'S2-50','mean_mse'),True),fmt(values(frame,model,'S2-50','coexpression_pearson'),True)])
 lopo_rows=[]
 for sample in ['NCBI783','NCBI785','TENX95','TENX99']:
  f=lo[(lo.sample_id==sample)&(lo.panel=='S2-50')]
  lopo_rows.append([sample,fmt(f[f.model=='B1'].median_pearson),fmt(f[f.model=='B2'].median_pearson,True),fmt(f[f.model=='B0'].mean_mse),fmt(f[f.model=='B1'].mean_mse),fmt(f[f.model=='B2'].mean_mse,True)])
 config_rows=[]
 for panel,s in lock['selected'].items():config_rows.append([panel,s['ridge']['alpha'],s['mlp']['hidden'],s['mlp']['dropout'],s['lopo_fixed_epochs']])
 thresholds=[]
 for panel in ['S2-50','S2-250','S2-280']:
  for model in ['B1','B2']:
   thresholds.append([panel,model]+[fmt(values(p,model,panel,f'n_pearson_gt_{t}'),model=='B2',1)+' / '+fmt(100*values(p,model,panel,f'fraction_pearson_gt_{t}'),model=='B2',1)+'%' for t in ['0','0.1','0.2','0.3']])
 top_bottom=[]
 for panel in ['S2-50','S2-250','S2-280']:
  for group in ['top10','bottom10']:
   f=rank[(rank.panel==panel)&(rank.model=='B2')&(rank['list']==group)].sort_values('rank')
   top_bottom.append(f'- **{panel}, B2 {group}:** '+', '.join(f'{r.gene} ({r.pearson:.3f})' for r in f.itertuples())+'.')
 text=f'''# C1 report: morphology-only spatial expression prediction

**Status: {status}.** The frozen generic ResNet18 features contain measurable local morphology-expression signal in this four-patient development cohort. This is not external validation or a publication-level conclusion.

C1 v2 was executed with the frozen C0 primary split: NCBI783 and NCBI785 for training (7,015 locations), TENX95 for validation (7,760), and TENX99 for final evaluation (20,761). The full cohort has 35,536 paired locations and 280 shared biological targets. The S1 target is **FOXA1**. S2-280 is the **full shared assay panel**, not the whole transcriptome.

## Primary results

All correlation summaries describe locations within the single primary test patient, TENX99. B2 values are mean ± sample standard deviation over seeds 0, 1, 2; this is training-seed variability, not patient-level uncertainty. B0 predictions are constant, so its correlations and coexpression statistic are mathematically undefined.

{primary_table}

Training time is the sum over the four primary settings for each seed; it excludes tuning, controls, validation, inference and artifact writing. Peak RSS is the maximum sampled process RSS across those four fits, including resident data; it is not isolated model memory. CPU: Intel Core i7-14700, four computation threads. Timing can be very short because all image embeddings were already cached in C0. Per-run training/validation/inference times, parameter counts and dimensions are in [C1_RESOURCE_PROFILE.csv](C1_RESOURCE_PROFILE.csv).

### Prediction error

MSE is measured after inverting C1 target scaling, in the original frozen C0 normalized target units.

{table(['Setting','B0 MSE','B1 MSE','B2 MSE','B2 relative MSE reduction vs B0'],mse_rows)}

Every primary B1/B2 setting and all primary B2 seeds improve mean MSE over B0. The gene-level CSV also contains Spearman, MSE, MAE, target variance and prediction variance. Summary CSVs include means, medians and every required Pearson-threshold fraction.

## Negative controls

NC1 permutes whole expression rows among training locations. NC2 permutes whole embedding rows. Both preserve validation/test data and retain within-row gene structure. Each uses seeds 0, 1, 2 with primary-selected hyperparameters; B2 early stopping still uses the unchanged validation patient. No hyperparameter grid was expanded after test access.

{table(['Model','Condition','S2-50 median Pearson','MSE','Coexpression preservation'],nc_rows)}

The correctly paired model exceeds every individual negative-control run in median gene Pearson for both B1 and B2. Negative-control correlations are variable across three seeds and do not equal zero exactly. This is descriptive evidence, not a permutation p-value.

**Coexpression preservation alone is not evidence of spatial prediction.** Permuted Ridge models reach approximately 0.81 coexpression correlation, exceeding the paired Ridge value, while location-wise gene prediction largely collapses. Row permutation preserves target-vector covariance, which a multivariate fitted model can retain without faithful spatial prediction. The coexpression requirement is met as a descriptive structural result, but it cannot establish biological fidelity.

![Negative controls](figures/negative_controls.png)

## Conditional four-patient robustness sensitivity analysis — not unbiased validation

Each fold trains on the other three patients, fits and saves its own feature and target scalers, and uses the already selected primary hyperparameters. B2 trains for a fixed epoch count from the primary runs, never early-stopping on the held-out fold patient. The frozen gene panels and primary-selected configuration make these estimates conditional. They are kept separate from the primary result.

{table(['Held-out sample','B1 median Pearson','B2 median Pearson','B0 MSE','B1 MSE','B2 MSE'],lopo_rows)}

This table shows S2-50. [C1_ROBUSTNESS_LOPO.csv](C1_ROBUSTNESS_LOPO.csv) contains all four settings and every seed/fold. Positive multi-gene signal is not confined to one patient. However, correlation and calibration differ: for S1/FOXA1 on NCBI785, both model families have positive correlation but worse mean MSE than B0; Ridge also slightly worsens S1 MSE on NCBI783. C1 therefore does not establish uniform benefit for every patient/gene.

![Conditional LOPO](figures/conditional_lopo.png)

## Gene-level predictability

Counts and fractions below use all genes in each panel as the denominator. B2 shows mean ± SD of seed-specific counts/fractions. Undefined B0 correlations are not treated as measured zero correlations.

{table(['Panel','Model','Pearson > 0: count / fraction','> 0.1','> 0.2','> 0.3'],thresholds)}

The following rankings use mean Pearson across B2 seeds on the primary test patient and are **descriptive only**. They were not used to select models, panels or map genes. Equivalent B1 lists and seed SDs are in [C1_GENE_RANKINGS.csv](C1_GENE_RANKINGS.csv).

{chr(10).join(top_bottom)}

Histograms and target-variance/Pearson scatterplots are in `figures/*_gene_predictability.png`. Predictability is heterogeneous and some genes have near-zero or negative correlations. No pathway or biological explanation claims are made.

## Spatial-map audit

Map genes were frozen before primary test access using training log-count variance within each panel: S1 FOXA1; S2-50 FOXA1/LUM/CCDC80; S2-250 FOXA1/MRC1/KLRF1; S2-280 FOXA1/IGSF6/UCP1. Seed 0 is the predetermined visualization seed. Every map uses identical coordinate limits and extent for H&E, measured, predicted and absolute-error panels; measured/predicted expression share color limits.

The inspected FOXA1 and LUM maps show nonconstant spatial variation, including broad correspondence, but understate some regional expression contrasts. The low-variance UCP1 example shows that weak sparse signals are not recovered well. Error maps expose regional failures rather than hiding them. Thumbnail H&E does not validate cell-level registration. Coexpression heatmaps use the same frozen ordering for measured, predicted and difference matrices.

![S1 FOXA1](figures/S1_FOXA1_B2_spatial.png)

## Frozen model choices and preprocessing

{table(['Setting','Ridge alpha','MLP hidden width','Dropout','Fixed LOPO epochs'],config_rows)}

Ridge selected alpha 100, the upper boundary of the predefined grid, in all settings. The grid was not extended after test evaluation. MLP tuning comprised only the four specified width/dropout combinations with seed 0, chosen by validation MSE. AdamW learning rate 0.001, weight decay 0.0001, batch size 512, maximum 80 epochs, patience 10 and minimum improvement 0.0001 were fixed beforehand. Primary early stopping selects the best improving validation checkpoint.

C1 scales features and targets using training-only means/standard deviations and inverts target scaling for reporting. LOPO restandardizes the unchanged C0 values using each fold-training pool. This cancels the original affine C0 shift/scale in model input targets apart from stored floating-point precision, but reporting units and frozen panel selection remain tied to C0. No C0 file, graph edge, neighbor expression, sparse mask or encoder weight was changed or used as a new C1 model input.

## Integrity and reproducibility

- 156 completed fits: 44 tuning candidates, 20 primary fits, 12 negative-control fits and 80 LOPO fits.
- 112 evaluated models; 15,125 gene-level metric rows; 268 append-only registry fit/evaluation events.
- All 20 scaler pairs recomputed from their declared training populations and verified.
- All 112 stored prediction arrays reproduced exactly by reloading model checkpoints (maximum absolute difference 0).
- All 261 frozen C0 file hashes still match. Primary test access occurred after the selection lock was written.
- Eight C1 unit tests and the 28 C0 tests pass. An initial single-target Ridge shape error was fixed before any successful candidate fit; its failure log and regression test are retained.

See [logs/INTEGRITY_AUDIT.json](logs/INTEGRITY_AUDIT.json), [C1_SPEC_LOCK.json](C1_SPEC_LOCK.json), and `C1_SPEC_LOCK.sha256`. Fitting source commit: `{lock['source_commit']}`. Reporting source and all source/artifact hashes are recorded in the final lock. Large model/prediction files remain locally in this repository with recorded hashes; source, metrics, figures, scalers, logs and report are Git-versioned.

## Acceptance decision and limitations

`PASS_DEVELOPMENT` is supported by primary error improvements, substantially better paired-versus-permuted gene prediction, a nontrivial predictable gene subset, signal across all conditional LOPO patients, nonzero measured/predicted coexpression association, nonconstant spatial maps, and successful integrity/replay audits. These are checkpoint judgments; no new significance test or universal PCC cutoff is asserted.

Important limitations:

- Only four patients, with one primary test patient; locations are correlated observations, not independent biological replicates. No spot-level significance tests or population confidence intervals are reported.
- Conditional LOPO reuses frozen panels and primary-selected configurations. It is not unbiased cross-validation, independent validation, or external validation.
- The primary split also separates source studies, and assay/preparation/stain differences can confound transfer. Metadata patient grouping is inherited from HEST.
- Frozen generic ImageNet ResNet18 embeddings and simple CPU baselines were used. Encoder optimality, cell/nuclear interpretation and pathology-specific representations remain untested.
- The targeted assay supports 280 shared biological genes only; no 1,000/3,000-gene or whole-transcriptome result exists.
- C0 log-count normalization retains sequencing-depth and composition effects; correlations do not establish a mechanistic morphology-expression relation.
- Gene–gene structure can remain high under permutation controls. It is not a standalone success criterion for spatial reconstruction.
- Predictability varies across genes/patients; maps smooth or attenuate some regions, and positive correlation does not guarantee calibrated error improvement.
- Three seeds and a small fixed validation search give limited sensitivity information. Alpha selection at the grid boundary is documented, not retrospectively retuned.
- H&E thumbnail review cannot certify cell-level alignment. No biological discovery, faithful explanation, clinical utility or publication-level generalization is claimed.

## Exact recommendation for C2

**Proceed, after review, to a separately frozen C2 single-scale context comparison.** Keep the C0 embeddings, primary splits, target panels, normalization and C1 baselines unchanged. Compare local-neighborhood mean pooling plus MLP against compact single-scale GCN/GraphSAGE baselines, with graph geometry and hyperparameters chosen on validation patients only. Evaluate both correlation and error, report gene-specific contextual gains and patient dependence, and distinguish graph benefit from simple spatial averaging. Treat the results as development evidence; do not introduce S3, multiresolution fusion or biological explanations at this step.

C2 has **not** been executed. Stop at this checkpoint.
'''
 (OUT/'C1_REPORT.md').write_text(text)
 decision={'status':status,'primary_errors_better_than_mean':error_ok,'paired_pcc_exceeds_each_control':controls_ok,
           'positive_signal_all_four_conditional_lopo_patients':lopo_ok,'substantial_positive_gene_subset':subset_ok,
           'coexpression_caution':'Permutation controls retain coexpression; not independent proof of spatial fidelity',
           'visual_audit':'Representative maps show nonconstant broad variation with regional attenuation; sparse low-variance failure retained',
           'recommendation':'Review and proceed to separate C2 single-scale graph versus local averaging comparison; no C2 executed'}
 (OUT/'logs/C1_DECISION.json').write_text(json.dumps(decision,indent=2)+'\n')
 print(status)
if __name__=='__main__':main()
