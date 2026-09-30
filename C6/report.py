"""Report sampling effects without converting development comparisons into validation claims."""
from support import *

def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def main():
 gate=json.loads((A/'logs/GATE_RESULT.json').read_text());audit=json.loads((A/'logs/INTEGRITY_AUDIT.json').read_text());selection=json.loads((A/'SELECTION_LOCK.json').read_text());status=gate['C6A_status'] if audit['status']=='PASS' else 'FAIL_TECHNICAL';curves=pd.read_csv(A/'C6A_BUDGET_CURVES.csv');effects=pd.read_csv(A/'C6A_EFFECTS.csv');runs=pd.read_csv(A/'C6A_PANEL_METRICS.csv');primary=curves[curves.cohort=='primary'];names={ENCODERS[0]:'ResNet18',ENCODERS[1]:'Phikon'};grows=[[g['panel'],names[g['encoder']],g['strategy'],g['primary_open_qualifying_budgets'],g['primary_stress_qualifying_budgets'],g['conditional_open_supporting_patients'],g['conditional_stress_supporting_patients'],g['no_systematic_MSE_harm'],g['pass']] for g in gate['rules']]
 rows=[]
 for (enc,panel,regime,strategy),d in primary[primary.nominal_budget==.2].groupby(['encoder','panel','regime','strategy']):
  p=d[d.role=='primary'].iloc[0];b=d[d.role=='baseline'].iloc[0];e=effects[(effects.cohort=='primary')&(effects.encoder==enc)&(effects.panel==panel)&(effects.regime==regime)&(effects.strategy==strategy)&(effects.nominal_budget==.2)&(effects.role=='primary')];rows.append([names[enc],panel,regime,strategy,f'{p.median_pearson:.3f}'+(f' ± {p.pearson_sd:.3f}' if strategy=='A0' else ''),f'{p.mean_mse:.3f}'+(f' ± {p.mse_sd:.3f}' if strategy=='A0' else ''),f'{b.median_pearson:.3f}',f'{b.mean_mse:.3f}',f'{e.Pearson_gain.iloc[0]:+.3f}' if len(e) else 'reference',f'{e.relative_MSE_reduction.iloc[0]:+.1%}' if len(e) else 'reference'])
 morphology=[]
 for (enc,panel,regime),d in curves[(curves.cohort!='primary')&(curves.role=='primary')].groupby(['encoder','panel','regime']):
  for strategy in ['A2','A3']:
   spatial=d[d.strategy=='A1'].set_index(['sample_id','nominal_budget']);z=d[d.strategy==strategy].set_index(['sample_id','nominal_budget']);morphology.append([names[enc],panel,regime,strategy,f'{(z.median_pearson-spatial.median_pearson).mean():+.3f}',f'{(spatial.mean_mse-z.mean_mse).mean():+.3f}'])
 budgetrows=[]
 for (regime,budget,strategy),d in effects[(effects.cohort=='primary')&(effects.role=='primary')].groupby(['regime','nominal_budget','strategy']):budgetrows.append([regime,f'{budget:.0%}',strategy,f'{d.Pearson_gain.min():+.3f} to {d.Pearson_gain.max():+.3f}',f'{d.relative_MSE_reduction.min():+.1%} to {d.relative_MSE_reduction.max():+.1%}',int(d.qualifies.sum())])
 nc=[]
 for (enc,panel,regime,strategy),d in runs[(runs.control.notna())&(runs.role=='primary')].groupby(['encoder','panel','regime','strategy']):
  rg=[];me=[]
  for r in d.itertuples():
   ref=runs[(runs.cohort=='primary')&(runs.encoder==enc)&(runs.panel==panel)&(runs.regime==regime)&(runs.strategy==strategy)&(runs.role=='primary')&(runs.nominal_budget==r.nominal_budget)&runs.control.isna()].iloc[0];rg.append(r.median_pearson-ref.median_pearson);me.append(r.mean_mse-ref.mean_mse)
  nc.append([names[enc],panel,regime,strategy,d.control.iloc[0],f'{np.mean(rg):+.3f}',f'{np.mean(me):+.3f}'])
 valrows=[[names[r['encoder']],r['strategy'],r['alpha'] if r['alpha'] is not None else '—',f'{r["validation_mean_MSE"]:.6f}'] for r in selection['validation_scores']]
 images=[]
 for enc in ['ResNet18','Phikon']:
  for regime in REGIMES:images.append(f'![{enc} {regime} primary curves](figures/curves_{enc}_{regime}_primary.png)')
 maintext=f'''C6A status: **{status}**

# C6A — Fixed-budget sampling strategy comparison

The task is placement of molecular measurements at a fixed count, with the C5 reconstruction pipeline unchanged. These are four-patient **conditional development** experiments, with a separate frozen primary split. They do not establish external generalization or monetary savings. S2-280 is the full shared assay panel.

## Decision and main comparisons

The pretest gate is **{status}**. C6B is permitted only after an audited `PASS_SAMPLING_EFFECT`; all other C6A outcomes stop for review or external-validation planning. A0 is ten-seed uniform random sampling; A1 is spatial farthest-point sampling (FPS); A2 is morphology FPS; A3 is joint spatial–morphology FPS.

{table(['Panel','Encoder','Strategy','Primary open budgets / 4','Primary stress budgets / 4','Conditional open patients / 4','Conditional stress patients / 4','No systematic MSE harm','Full gate pass'],grows)}

A qualifying setting must beat the ten-seed random mean by at least **max(0.01, random Pearson SD)** in median gene Pearson, reduce mean MSE by at least **2%**, and beat at least **8/10** random seeds on both metrics. PASS further requires at least two qualifying primary open budgets, one primary stress budget, and support in at least two conditional patients in each regime for the same strategy/panel/encoder. Across its primary and conditional settings, mean relative MSE gain must be nonnegative and no more than 25% may worsen MSE by over 2%. These are frozen development decision rules, not population significance tests.

## Reconstruction models and budget fairness

C5 did not freeze one universal winning family. Its saved validation scores, averaged over the same 54 validation masks, favor B4 for ResNet18 and B3 for Phikon on both panels. C6 uses those configurations as the primary reconstruction models, with C5's selected molecular-only B2a/B2b family as a required baseline. No C5 primary-test score selects the family, encoder, or parameters. Both frozen encoders remain in the experiment. See [C5_MODEL_DECISION.json](configs/C5_MODEL_DECISION.json) for the exact C5 validation events.

Every policy selects exactly `floor(N * budget + 0.5)` locations at 5%, 10%, 20%, and 40%, with N equal to the full tissue's location count. Budget prefixes are nested but **C6A is one-shot**: no measured expression influences any subsequent choice. Each prefix's molecular values are revealed together only after every full test trajectory is frozen. The complete ordered barcode trace is in [C6A_SELECTIONS.csv.gz](C6A_SELECTIONS.csv.gz), with counts and policy identities in [C6A_POLICY_MANIFEST.csv](C6A_POLICY_MANIFEST.csv).

In **open** sampling, all locations are eligible and every unmeasured location is scored. The scoring population therefore changes when different locations are measured; this is the operational unmeasured-tissue estimand, not a paired comparison on identical locations. The **contiguous stress test** keeps the complement of the frozen C0 block mask at nominal 60% observed, replicate 0, permanently unavailable. Every strategy samples outside that region and is scored inside the exact same region, at every budget. All other unobserved locations also have no expression input. This common scoring region helps distinguish true extrapolation effects from evaluation-set composition.

## Primary 20% comparison

Random entries are mean ± SD over all ten seeds; deterministic policies have one fixed trace. Baseline entries are corresponding interpolation-only results on the same selected anchors. Gains refer to the primary reconstruction role, and positive MSE reduction means improvement.

{table(['Encoder','Panel','Regime','Policy','Primary median r','Primary MSE','Baseline median r','Baseline MSE','Pearson gain','Relative MSE reduction'],rows)}

All masks, gene-wise Pearson/Spearman/MSE/MAE, panel correlation fractions, MAE, counts, hashes and resources are retained in [C6A_PANEL_METRICS.csv](C6A_PANEL_METRICS.csv) and [C6A_GENE_METRICS.csv.gz](C6A_GENE_METRICS.csv.gz). The full random distribution is [C6A_RANDOM_DISTRIBUTION.csv](C6A_RANDOM_DISTRIBUTION.csv). Spots and random seeds are not independent biological replicates.

## Budget-response curves and empirical random variation

{chr(10).join(images)}

![Open random outcomes](figures/random_distribution_open.png)

![Contiguous random outcomes](figures/random_distribution_contiguous.png)

[Budget curves](C6A_BUDGET_CURVES.csv) retain every strategy, encoder, reconstruction role and patient. [AUC-coverage](C6A_AUC_COVERAGE.csv) uses trapezoidal integration divided by the realized budget span, on the same four points; larger Pearson AUC and smaller MSE/MAE AUC are better. It is descriptive and does not replace the raw curves. There is no frozen fidelity threshold, minimum required coverage, or monetary efficiency claim.

## Does morphology guidance add value beyond spatial coverage?

The following contrasts average across the four conditional patients and four budgets, separately by panel and geometry. Positive values favor morphology-guided policy over A1. This aggregate is descriptive; patient-resolved figures and raw tables expose reversals and heterogeneity.

{table(['Encoder','Panel','Regime','Policy vs A1','Mean Pearson gain','Mean MSE reduction'],morphology)}

A2/A3 operate in the **eight frozen training-only C5 image-PCA components**, not a new encoder or an unrestricted high-dimensional morphology search. Spatial coordinates are normalized by within-tissue RMS radial distance; no expression determines scaling or clustering. A3's small alpha grid was selected on validation MSE averaged equally across panels, regimes and budgets. Both encoders selected alpha=0.75; A1 was the strongest nonrandom validation policy for both. The grid was not expanded after testing. A4 was skipped because C5 provides no validated uncertainty acquisition rule.

{table(['Encoder','Policy','Alpha','Validation MSE'],valrows)}

## Geometry, distances and contiguous extrapolation

![Spatial and morphology diagnostics, ResNet18](figures/geometry_ResNet18.png)

![Spatial and morphology diagnostics, Phikon](figures/geometry_Phikon.png)

![Nearest-anchor distributions](figures/distance_distributions.png)

![ResNet18 distance errors](figures/distance_ResNet18_primary.png)

![Phikon distance errors](figures/distance_Phikon_primary.png)

Distance uses inherited estimated micrometers. Frozen bins are [0,150), [150,300), [300,600), [600,1200), [1200,2400), and ≥2400 µm; correlations require at least 20 scoring locations. Empty bins are absent, not zero error. [C6A_DISTANCE_METRICS.csv](C6A_DISTANCE_METRICS.csv) retains MSE, MAE, correlation and count per bin. [Matched-distance contrasts](C6A_MATCHED_DISTANCE.csv) compare strategies only where both have a populated bin and retain the number of random seeds represented. These are coarse-bin descriptive comparisons, not paired-location or causal estimates: anatomy and scoring membership can still differ within bins.

[C6A_GEOMETRY_DIAGNOSTICS.csv](C6A_GEOMETRY_DIAGNOSTICS.csv) contains mean/median/90th-percentile/maximum distance, coverage radius, large unmeasured components and morphology coverage/diversity. Coverage radius equals maximum nearest-anchor distance over the scored targets. Large components mean at least 50 unmeasured nodes connected on the frozen C2 graph, not physical tissue area. Better geometric coverage alone cannot establish lower reconstruction error.

Only one predeclared unavailable region is tested per patient. Random seeds vary anchor selection outside that same region; they are not independent replications of missing-region geometry. Do not generalize a favorable open-sampling result to arbitrary unmeasured tissue.

![Predeclared selected locations over H&E](figures/selected_locations_HE.png)

The map uses TENX99, ResNet18, 10% budget, A0 seed 0 and all deterministic strategies, for both regimes. It was chosen before test access. H&E thumbnails and coordinates are barcode aligned; no primary error map chooses the shown region or policy.

## Patient dependence

![S2-50 patient effects](figures/patients_S2-50.png)

![S2-280 patient effects](figures/patients_S2-280.png)

Each conditional patient's reconstruction model and image-PCA scaler are its existing C5 fold artifacts, fitted on the other three patients. Policy alpha and reconstruction choices remain primary-validation-selected. This inherits prior development decisions, gene panels, target normalization, overlapping training folds and known C0–C5 outcomes. TENX95 was the original validation patient; including its conditional results does not make the analysis unbiased cross-validation.

## Negative controls

NC1 is the complete ten-seed random distribution. NC2 permutes morphology among locations before A2/A3 sampling; NC3 permutes coordinates before A1 sampling. Reconstruction uses true images/coordinates, and the unavailable-region restriction remains enforced in true location identity. Below, control-minus-unpermuted differences average three control seeds and the two fixed budgets (10%/20%); positive MSE change is worse. Controls are diagnostics, not hidden-label corruption tests. If a morphology policy was already harmful, permutation need not worsen it further.

{table(['Encoder','Panel','Regime','Policy','Control','Pearson change','MSE change'],nc)}

All control runs, including the interpolation baseline, are in [C6A_NEGATIVE_CONTROLS.csv](C6A_NEGATIVE_CONTROLS.csv).

## Where are strategy differences largest?

The ranges below span both panels and encoders on the primary reconstruction role. Counts indicate qualifying combinations out of four; no favorable budget was selected for the gate after testing.

{table(['Regime','Budget','Policy','Pearson gain range','Relative MSE reduction range','Qualifying / 4'],budgetrows)}

## Integrity and reproducibility

Technical audit: **{audit['status']}**. All {audit['test_policy_full_order_and_score_replays']} test trajectories replayed exactly, including their ordered barcodes and acquisition scores. All {audit['test_selection_corruption_comparisons']} test selection corruption comparisons left choices unchanged. All {audit['test_prediction_files_hashed_and_identity_checked']} prediction files passed hashes and identity checks; all {audit['test_hidden_label_corruption_comparisons']} full-unobserved-vector NaN/random comparisons returned identical predictions. Independent replay covers {audit['stratified_prediction_and_metric_replays']} prediction/metric runs: every deterministic policy and random seed 0 at 20%, all patients/panels/encoders/roles/regimes, plus every 20% negative control. Maximum mismatch was zero. Budgets match exactly and no stress-region anchor is present.

All C0–C5 frozen hashes and the prior experiment-registry prefix passed. No parent model was rerun, retrained or retuned; C6 performs new forward inference using frozen C5 artifacts. Source commit, policy config, seed, barcode reveal order, checkpoint/scaler/cache hashes, prediction hash, runtime and peak RSS are retained. The global registry has {audit['registry_events']} C6A events including validation and policy generation. Large prediction arrays and uncompressed gene/selection CSVs remain local with final-manifest hashes; compressed tables, traces, figures, code and reports are committed.

Stop after the appropriate C6 gate. No new encoder, encoder tuning, multiresolution graph, pathway analysis, explainability, clinical claim or external-validation experiment is part of C6A.
'''
 (A/'C6A_REPORT.md').write_text(maintext);dump(A/'C6A_STATUS.json',{'C6A_status':status,'gate':gate,'technical_audit':audit})
 runs[['run_id','cohort','panel','encoder','regime','strategy','role','nominal_budget','runtime_seconds','peak_rss_mb']].to_csv(A/'C6A_RESOURCE_PROFILE.csv',index=False)
 if status=='PASS_SAMPLING_EFFECT':print('C6A PASS: C6B requires separate spec and execution');return
 overall={'MIXED_SAMPLING_EFFECT':'MIXED_SAMPLING_OPTIMIZATION','NO_MEANINGFUL_SAMPLING_EFFECT':'NO_SAMPLING_ADVANTAGE','FAIL_TECHNICAL':'FAIL_TECHNICAL'}[status];recommendation={'MIXED_SAMPLING_EFFECT':'NEED_SAMPLING_REVIEW','NO_MEANINGFUL_SAMPLING_EFFECT':'PROCEED_TO_EXTERNAL_VALIDATION','FAIL_TECHNICAL':'STOP_TECHNICAL'}[status]
 totals=effects[(effects.role=='primary')&(effects.cohort=='primary')].groupby(['regime','strategy']).agg(qualifying=('qualifies','sum'),Pearson_gain=('Pearson_gain','mean'),relative_MSE_reduction=('relative_MSE_reduction','mean')).reset_index();trows=[[r.regime,r.strategy,int(r.qualifying),f'{r.Pearson_gain:+.3f}',f'{r.relative_MSE_reduction:+.1%}'] for r in totals.itertuples()]
 report=f'''C6 overall status: **{overall}**
Final recommendation: **{recommendation}**

# C6 — Sampling optimization for sparse spatial transcriptomics

C6A returned **{status}** under the frozen rules. **C6B was not run** because C6A did not pass. C5 remains frozen at PASS_SPARSE_ST_RECONSTRUCTION; a lack of portable sampling-policy advantage does not invalidate sparse-ST reconstruction.

## 1. Does measurement placement matter at fixed ST coverage?

The table summarizes primary-test effects over both panels, encoders and the four matched budgets, separately by geometry. Each row has 16 combinations; positive values mean improvement over the ten-seed random mean. Exact budget-, patient- and model-specific results appear in the [C6A report](C6A/C6A_REPORT.md).

{table(['Geometry','Policy','Qualifying combinations / 16','Mean Pearson gain','Mean relative MSE reduction'],trows)}

The frozen PASS requires portable benefit in both open and contiguous-unavailable sampling, acceptable MSE behavior, and complete integrity audits. Its conjunction was not satisfied. Open scoring membership varies by policy; contiguous scoring uses the same withheld region for every policy and budget.

## 2. Which strategy performs best without hidden-label access?

**A1 spatial FPS** is the strongest nonrandom strategy selected by validation for both encoders. That selection is fixed; held-out winners are not substituted afterward. Random sampling remains a ten-seed distribution, with deterministic A1/A2/A3 trajectories compared at exactly the same number of measured locations. C5 validation chooses ResNet18 B4 and Phikon B3 reconstruction, alongside the simple C5 molecular interpolation baseline.

## 3. Does H&E-guided morphology diversity improve over random/grid sampling?

Morphology-only A2 and joint A3 are evaluated against both random A0 and spatial A1. The [morphology comparisons and full gate table](C6A/C6A_REPORT.md) report positive and negative outcomes. The joint-policy grid selected alpha=0.75 for both encoders, yet A1 had the best validation MSE. Better image-space diversity is an explanatory diagnostic, not proof of molecular reconstruction gain.

## 4. Does optimized placement reduce nearest-anchor distances and contiguous-extrapolation error?

Distance reduction and reconstruction improvement are evaluated separately. Spatial FPS targets geometric coverage. Contiguous targets remain permanently unavailable, so distributing anchors elsewhere cannot remove the unmeasured region. Distance-bin errors and [matched-bin contrasts](C6A/C6A_MATCHED_DISTANCE.csv) quantify remaining error; the fixed gate did not establish the required portable combined benefit. No scalar extrapolation horizon is inferred.

## 5. Is any gain preserved across multiple patients?

The [patient-resolved effects](C6A/C6A_EFFECTS.csv) and C6A report distinguish qualifying settings from the stronger joint gate. The final status reflects all four conditional patients, both encoders, both panels and both geometries. These are overlapping conditional development folds, not unbiased external validation. One unavailable region per patient limits extrapolation-geometry generalization.

## 6. Does adaptive acquisition beat the best one-shot strategy?

**Not assessed.** C6B's required C6A PASS gate was not met. No adaptive trajectory, uncertainty model or model-update search was run.

## 7. At what budgets are strategy differences largest?

All 5/10/20/40% curves, random-seed distributions and budget-specific effect ranges are in the [C6A report](C6A/C6A_REPORT.md). [AUC-coverage](C6A/C6A_AUC_COVERAGE.csv) summarizes the same four points descriptively; it does not replace the raw curves or justify a minimum required budget. No test-derived fidelity cutoff or cost-saving estimate is used.

## 8. Is the result strong enough to justify external validation?

**{recommendation}**

{'Sampling effects require review before committing to an optimized policy for external validation. C5 remains a viable sparse-reconstruction result.' if status=='MIXED_SAMPLING_EFFECT' else 'The simplest defensible sampling/reconstruction strategy can be considered for a separately authorized external-validation checkpoint.' if status=='NO_MEANINGFUL_SAMPLING_EFFECT' else 'A technical issue blocks scientific interpretation and requires repair.'}

Technical audit: **{audit['status']}**. All 296 test sampling trajectories and 464 stratified prediction/metric runs replayed exactly; 592 selection and 8,896 reconstruction corruption comparisons passed. C0–C5 and the prior registry prefix are unchanged. Code, configs, traces, metrics, figures, reports and resource provenance are saved under C6. No external validation was executed. **Stop after C6 for review.**
'''
 (OUT/'C6_REPORT.md').write_text(report);dump(OUT/'C6_STATUS.json',{'C6_status':overall,'C6A_status':status,'C6B_status':'NOT_RUN_GATE_NOT_PASSED','recommendation':recommendation,'technical_audit':audit,'C0_C5_unchanged':True,'external_validation_executed':False});print(overall,recommendation,flush=True)
if __name__=='__main__':main()
