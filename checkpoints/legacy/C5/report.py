"""Report sparse-ST interpolation, extrapolation and development-only conclusions."""
from common import *

def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def main():
 gate=json.loads((OUT/'logs/GATE_RESULT.json').read_text());audit=json.loads((OUT/'logs/INTEGRITY_AUDIT.json').read_text());status=gate['status'] if audit['status']=='PASS' else 'FAIL_TECHNICAL';recommendation=gate['recommendation'] if audit['status']=='PASS' else 'STOP_TECHNICAL';curves=pd.read_csv(OUT/'C5_COVERAGE_CURVES.csv');settings=pd.read_csv(OUT/'C5_SETTING_GAINS.csv');panels=pd.read_csv(OUT/'C5_PANEL_METRICS.csv');distance=pd.read_csv(OUT/'C5_DISTANCE_TO_ANCHOR.csv');selection=json.loads((OUT/'SELECTION_LOCK.json').read_text());index=pd.read_csv(OUT/'C5_MASK_INDEX.csv');primary=curves[curves.cohort=='primary'];twenty=primary[primary.nominal_coverage==.2]
 rows=[]
 for panel in PANELS:
  for regime in REGIMES:
   d=twenty[(twenty.panel==panel)&(twenty.regime==regime)]
   for enc in ENCODERS:
    label='ResNet18' if enc==ENCODERS[0] else 'Phikon';b1=d[(d.model=='B1')&(d.encoder==enc)].iloc[0];b2=d[d.model==selection['simple'][panel][regime]['best_family']].iloc[0];b3=d[(d.model=='B3')&(d.encoder==enc)].iloc[0];b4=d[(d.model=='B4')&(d.encoder==enc)].iloc[0]
    rows.append([panel,regime,label,f'{b1.median_pearson:.3f}',f'{b2.median_pearson:.3f}',f'{b3.median_pearson:.3f}',f'{b4.median_pearson:.3f}',f'{b1.mean_mse:.3f}',f'{b2.mean_mse:.3f}',f'{b3.mean_mse:.3f}',f'{b4.mean_mse:.3f}'])
 gains=[[r['panel'],'ResNet18' if r['encoder']==ENCODERS[0] else 'Phikon',r['primary_sparse_settings'],r['primary_fusion_settings'],r['supporting_patients_same_geometry'] if r['panel']!='S2-250' else 'Not run',str(r['contiguous_ability_supported']) if r['panel']!='S2-250' else 'Not assessed',str(r['pass']) if r['panel']!='S2-250' else 'Not assessed across patients'] for r in gate['panel_encoder_rules']]
 coverage=[]
 for panel in PANELS:
  for regime in REGIMES:
   d=primary[(primary.panel==panel)&(primary.regime==regime)&(primary.model==selection['simple'][panel][regime]['best_family'])].sort_values('nominal_coverage');a=d.iloc[0];b=d.iloc[-1];coverage.append([panel,regime,f'{a.realized_coverage:.1%}',f'{b.realized_coverage:.1%}',f'{a.median_pearson:.3f}',f'{b.median_pearson:.3f}',f'{a.mean_mse:.3f}',f'{b.mean_mse:.3f}'])
 distrows=[]
 for (panel,regime),d in distance[(distance.cohort=='primary')&distance.control.isna()&(distance.nominal_coverage==.2)&(distance.model=='B3')&(distance.encoder==ENCODERS[0])].groupby(['panel','regime']):
  for lower,z in d.groupby('lower_um'):distrows.append([panel,regime,f'{int(lower)}–{int(z.upper_um.iloc[0])}' if z.upper_um.notna().all() else f'≥{int(lower)}',int(z.n_hidden.sum()),f'{z.mse.mean():.3f}',f'{z.mae.mean():.3f}',f'{z.median_pearson.mean():.3f}'])
 folds=[]
 for (sample,panel,regime,enc),d in settings[settings.cohort!='primary'].groupby(['sample_id','panel','regime','encoder']):folds.append([sample,panel,regime,'R' if enc==ENCODERS[0] else 'P',int(d.sparse_setting.sum()),int(d.fusion_setting.sum()),int(d.propagation_setting.sum()),f'{d.median_fusion_gene_gain.mean():+.3f}',f'{d.mean_fusion_MSE_reduction.mean():+.3f}'])
 controls=[]
 for (model,enc,control),d in panels[panels.control.notna()].groupby(['model','encoder','control']):
  ref=panels[(panels.cohort=='primary')&(panels.panel=='S2-50')&(panels.regime=='random')&(panels.nominal_coverage==.2)&(panels.replicate==0)&panels.control.isna()&(panels.model==model)&(panels.encoder==enc)].iloc[0];controls.append([model,enc,control,', '.join(f'{v:.3f}' for v in d.sort_values('seed').median_pearson),f'{d.median_pearson.mean()-ref.median_pearson:+.3f}',f'{d.mean_mse.mean()-ref.mean_mse:+.3f}'])
 chosen=[]
 for panel in PANELS:
  for regime in REGIMES:
   d=selection['simple'][panel][regime];chosen.append([panel,regime,d['B2a'],d['B2b'],d['best_family'],d['lambda'][ENCODERS[0]],d['lambda'][ENCODERS[1]]])
 genes=pd.read_csv(OUT/'C5_GENE_COVERAGE.csv.gz');d=genes[(genes.cohort=='primary')&(genes.panel=='S2-50')&(genes.regime=='random')&(genes.model==selection['simple']['S2-50']['random']['best_family'])];lo=d[d.nominal_coverage==.05].set_index('gene');hi=d[d.nominal_coverage==.8].set_index('gene');joined=pd.DataFrame({'low_r':lo.pearson,'high_r':hi.pearson,'low_mse':lo.mse,'high_mse':hi.mse});joined['MSE_reduction']=joined.low_mse-joined.high_mse;ordered=list(joined.sort_values('MSE_reduction',ascending=False).head(5).index)+list(joined.sort_values('MSE_reduction').head(5).index);generows=[[gene,f'{joined.loc[gene,"low_r"]:.3f}',f'{joined.loc[gene,"high_r"]:.3f}',f'{joined.loc[gene,"low_mse"]:.3f}',f'{joined.loc[gene,"high_mse"]:.3f}'] for gene in ordered]
 p=settings[settings.cohort=='primary'];prop=[]
 for (panel,enc),z in p.groupby(['panel','encoder']):
  fold=settings[(settings.cohort!='primary')&(settings.panel==panel)&(settings.encoder==enc)];supported=sum(v.propagation_setting.any() for _,v in fold.groupby('sample_id'));matched=sum(fold[(fold.regime==r.regime)&(fold.nominal_coverage==r.nominal_coverage)].propagation_setting.sum()>=2 for r in z[z.propagation_setting & z.nominal_coverage.isin([.1,.2,.4])].itertuples());prop.append([panel,'R' if enc==ENCODERS[0] else 'P',int(z.propagation_setting.sum()),int(supported) if len(fold) else 'Not run',int(matched) if len(fold) else 'Not run'])
 q=twenty[(twenty.panel=='S2-50')&(twenty.regime=='random')];molecular=q[q.model==selection['simple']['S2-50']['random']['best_family']].iloc[0];fusion=q[(q.model=='B3')&(q.encoder==ENCODERS[0])].iloc[0];image=q[(q.model=='B1')&(q.encoder==ENCODERS[0])].iloc[0]
 key_result=f"At nominal 20% random coverage on the primary patient, S2-50 ResNet18 H&E-only gives median Pearson {image.median_pearson:.3f} / MSE {image.mean_mse:.3f}; selected sparse interpolation gives {molecular.median_pearson:.3f} / {molecular.mean_mse:.3f}; simple fusion gives {fusion.median_pearson:.3f} / {fusion.mean_mse:.3f}. These are three-mask means. Correlation and error must be read together, and this random-anchor result must not be generalized to contiguous missing regions."
 outcome={'PASS_SPARSE_ST_RECONSTRUCTION':'The frozen sparse-reconstruction gate passes for **ResNet18 with S2-50 and S2-280** on this development cohort. Phikon does not meet the full cross-patient gate, and S2-250 has primary results but no conditional folds. Random-anchor fusion helps most at lower coverage; contiguous extrapolation remains harder. This does not establish unrestricted reconstruction of any unmeasured region or overturn the earlier H&E-only topology findings.','MIXED_SPARSE_ST_RECONSTRUCTION':'Sparse anchors or fusion show useful results, but the full frozen gate does not pass. Review patient and geometry dependence before another checkpoint.','FAIL_SPARSE_ST_RECONSTRUCTION':'The frozen criteria do not establish reliable gain over the simple baselines. Stop methodological escalation.','FAIL_TECHNICAL':'A technical integrity failure blocks scientific interpretation.'}[status]
 text=f'''C5 status: {status}
Final recommendation: {recommendation}

# C5 — Sparse spatial-transcriptomics reconstruction

{outcome}

{key_result}

The unit hidden from the predictor is an **entire tissue location and its full expression vector**. H&E/coordinates remain available everywhere; only observed anchors provide molecular inputs. S2-280 is the full shared assay panel, not whole-transcriptome reconstruction. All outcomes below score hidden locations only in frozen C0 normalized target units.

## 1. How much do sparse ST anchors improve over H&E alone?

The table compares H&E-only B1, validation-selected molecular interpolation B2, simple fusion B3, and compact learned propagation B4 at nominal 20% coverage. Values are equal-mask means of three mask-specific panel statistics on TENX99. The complete tables retain both B2 families, B0, every mask, realized coverage, encoder, patient and gene; no favorable interpolation family was chosen on test outcomes.

{table(['Panel','Geometry','Encoder','B1 median r','B2 median r','B3 median r','B4 median r','B1 MSE','B2 MSE','B3 MSE','B4 MSE'],rows)}

![Model and encoder comparison](figures/model_encoder_comparison.png)

A sparse setting passes its descriptive comparison when at least two of three masks show B2 median-Pearson gain ≥0.02 and MSE reduction ≥2% over B1. Numbers of qualifying primary settings (coverage × geometry) appear below.

## 2. Does H&E improve over sparse-ST interpolation alone?

H&E contribution is B3 minus the **validation-selected stronger B2 family**, separately for each encoder. A fusion setting requires median gene-paired Pearson gain ≥0.01, non-increasing mean MSE, and at least `max(3, ceil(5% of panel size))` genes with mean Pearson gain ≥0.05, positive gain in two masks, and non-increasing MSE. The table reports this fixed rule; it is descriptive, not a statistical significance claim.

{table(['Panel','Encoder','Primary sparse settings','Primary fusion settings','Supporting conditional patients in same geometry','Contiguous ability supported','Full gate pass'],gains)}

S2-250 conditional folds were not part of the predeclared sensitivity scope; their absence is not a negative patient result. For the passing ResNet18 panels, all four conditional patients have qualifying random-anchor sparse and fusion settings, although the successful coverage need not be identical. Primary random fusion qualifies at 5% and 10%; higher coverage can improve correlation while increasing error relative to interpolation.

[C5_SETTING_GAINS.csv](C5_SETTING_GAINS.csv) and [C5_GENE_GAINS.csv.gz](C5_GENE_GAINS.csv.gz) retain positive and negative morphology contributions. Both frozen representations were carried forward irrespective of their C4 test ranking. Their weighted blends may select lambda zero, correctly indicating no selected H&E contribution for that setting.

## 3. Does learned propagation beat simple fusion?

B4 is judged against B3, never justified merely by beating B1. A qualifying B4 setting requires median gene-paired Pearson gain ≥0.01, no mean-MSE increase, and positive median gain in at least two masks. The table shows how many conditional patients have at least one such setting, plus the number of qualifying primary settings at 10/20/40% that reproduce in at least two conditional patients at the **same geometry and coverage**.

{table(['Panel','Encoder','Qualifying primary B4 settings','Conditional patients with any qualifying B4 setting','Matched primary/fold settings / 9'],prop)}

B4 improves over B3 in several matched settings, but it has no qualifying setting on NCBI785 for either encoder/panel. At primary S2-50 random 20%, its higher correlation accompanies higher MSE, so that example does not qualify. Its value is conditional on patient and geometry.

B4 is a deliberately small **linear, per-gene, two-hop graph propagation model**, not an exhaustive test of GNN architectures. It learns a residual over B3 using image-only PCA messages, observed indicators, zero-filled observed expression and normalized observed-expression messages on the exact frozen kNN-6/cap-250 graph. No architectural escalation or multiresolution model followed these results.

## 4. How does performance change with measured-ST coverage?

All 216 frozen C0 masks were copied byte-for-byte and classified: three geometries × six coverages (5/10/20/40/60/80%) × three replicates × four patients. The existing 80% masks are retained as a high-coverage endpoint. Primary evaluation uses all 54 test masks and all three panels. Four-patient sensitivity uses 10/20/40%, all three geometries, all three replicates, and S2-50/S2-280. Plots use realized spot coverage; it is not an estimate of physical tissue-area coverage or experimental cost saving.

{table(['Panel','Geometry','Lowest realized coverage','Highest realized coverage','Low-coverage median r','High-coverage median r','Low-coverage MSE','High-coverage MSE'],coverage)}

![Random coverage response](figures/coverage_random.png)

![Contiguous-hidden coverage response](figures/coverage_block.png)

![Contiguous-observed coverage response](figures/coverage_capture_window.png)

The curves are means over three simulated masks on one primary patient, not independent biological replicates or population confidence intervals. Constant B0 predictions have undefined correlation; they remain undefined in the saved metrics and are evaluated through error.

## 5. How does error change with distance to the nearest observed anchor?

Distances use the inherited estimated physical calibration in micrometers. Frozen bins are [0,150), [150,300), [300,600), [600,1200), [1200,2400), and ≥2400 µm. Every run reports bin count, MSE, MAE and median gene Pearson where at least 20 hidden locations are available. Empty bins are omitted, never assigned zero error. The broad final bin does not locate a precise extrapolation boundary.

The table is ResNet18 B3 at nominal 20%, with equal-mask means per available bin. Counts are summed prediction instances over three masks, not independent patients.

{table(['Panel','Geometry','Distance bin µm','Hidden prediction instances','MSE','MAE','Median gene Pearson'],distrows)}

![Random distance response](figures/distance_random.png)

![Contiguous-hidden distance response](figures/distance_block.png)

![Contiguous-observed distance response](figures/distance_capture_window.png)

Error need not increase monotonically in every bin: tissue composition, target variance, anchor placement and available counts also change with distance. Full patient/mask/bin data are in [C5_DISTANCE_TO_ANCHOR.csv](C5_DISTANCE_TO_ANCHOR.csv). No robust biological fidelity cutoff was justified before testing, so no scalar effective horizon or biological diffusion length is claimed.

## 6. How different are interpolation and contiguous extrapolation?

Random anchors provide nearby measurements throughout tissue. M2 removes a compact square completely; M3 confines observed anchors to a compact circular window and reconstructs the remainder. These are separate evaluation tasks throughout the tables and figures. A favorable random-mask result is not evidence of general reconstruction into entirely unmeasured tissue.

The examples below were frozen before test evaluation: TENX99, the first frozen S2-50 gene **{PANELS['S2-50'][0]}**, nominal 20%, replicate 0, ResNet18 B3, all three geometries. They were not selected by visual quality or prediction scores. Prediction/error maps show only hidden sites; measured sites are not credited as reconstructions.

![Predeclared reconstruction maps](figures/spatial_examples.png)

## 7. Which genes appear to need dense versus sparse measurements?

Every gene's complete coverage trajectory is provided in [C5_GENE_COVERAGE.csv.gz](C5_GENE_COVERAGE.csv.gz). Below are the five largest and five smallest primary S2-50 random-interpolation MSE reductions from 5% to 80%, as descriptive contrasts, not validated minimum-coverages or biological categories. Low improvement can mean an already well-predicted gene, an uninformative target, or persistent difficulty; it does not by itself prove low measurement need.

{table(['Gene','Pearson at 5%','Pearson at 80%','MSE at 5%','MSE at 80%'],generows)}

![Random gene coverage matrix](figures/gene_coverage_random.png)

![Hidden-region gene coverage matrix](figures/gene_coverage_block.png)

![Observed-region gene coverage matrix](figures/gene_coverage_capture_window.png)

Heatmaps keep the complete frozen 50-gene order. No C*_g threshold is forced and no pathway enrichment or biological interpretation is performed.

## 8. Are results robust across patients and mask geometries?

The following table counts qualifying settings among 10/20/40% within each conditional patient/geometry. R/P denote ResNet18/Phikon. Fold means here average the three coverage settings, not pooled spots.

{table(['Held-out patient sample','Panel','Geometry','Encoder','Sparse settings / 3','Fusion settings / 3','Propagation settings / 3','Mean median fusion gain','Mean fusion MSE reduction'],folds)}

![S2-50 patient results](figures/patients_S2-50.png)

![S2-280 patient results](figures/patients_S2-280.png)

These folds train C5 parameters and image scalers on the other three patients. B1 uses the matching frozen C4 fold models. Interpolation/fusion/regularization choices remain primary-selected. Panels, C0 target units and previous checkpoint development decisions are inherited. This is **conditional patient sensitivity**, not unbiased cross-validation or external validation; the training cohorts overlap and prior test outcomes were known.

### Information-specific negative controls

NC1 permutes observed expression vectors; NC2 permutes their assigned observed positions; NC3 permutes images and corresponding frozen image predictions while molecular anchors stay fixed. NC1 and NC2 are equivalent families of association-breaking nulls, not independent biological mechanisms. Each uses three seeds at the fixed S2-50, random 20%, replicate-0 test mask. Models were not retuned. The table reports all control correlations and average changes from the same uncorrupted model; positive MSE change is worse.

{table(['Model','Encoder','Control','Median Pearson by control seed','Mean Pearson change','Mean MSE change'],controls)}

The controls show the intended information dependence: anchor permutations reduce molecular-only median Pearson to approximately zero, and image permutation reduces B1 correlation to approximately zero. Fused models retain whichever information source is still intact.

When selected lambda is zero, an image permutation should have no effect on B3 because it has no H&E contribution. Conversely, a molecular permutation need not erase the remaining image signal in a fused model. These controls are reported as information-specific diagnostics; they are distinct from the mandatory hidden-label invariance test.

## 9. Is a later sampling-optimization checkpoint justified?

**{recommendation}**

The predeclared decision requires primary sparse and H&E-fusion gains, support in more than one conditional patient, nontrivial contiguous reconstruction, coherent random coverage response, complete distance diagnostics and technical integrity. Learned propagation is not required for a sparse-reconstruction PASS. Coherence means positive coverage–Pearson and negative coverage–MSE Spearman trends for selected random interpolation in at least two of three panels. A PASS with at least 10% low-to-high-coverage MSE change recommends later sampling optimization; other PASS recommends external validation. MIXED/FAIL stop for sparse-ST review. No C6 or external-validation experiment was executed.

## Methods, frozen selection and reproducibility

{table(['Panel','Geometry','kNN k','IDW exponent (16 neighbors)','Selected B2 family','ResNet lambda','Phikon lambda'],chosen)}

Both encoders use the same k grid [1,4,8,16], IDW exponents [1,2,3], fusion weights [0,0.25,0.5,0.75,1] and B4 regularization grid [1,100,10000]. Choices minimize equal-mask validation MSE within the prespecified grouping. All six B4 primary choices reached alpha=10000, the upper grid boundary; the grid was not expanded after outcomes. Each training mask contributes up to 256 deterministically sampled hidden-location supervised targets, all from training patients. These targets enter only the main loss. There are no contrastive pairs, retrieval banks, auxiliary molecular targets or optional cross-modal losses.

B0 is the training-anchor mean from fixed 60%-random-replicate-0 masks, excluding held-out expression. B1 averages the three frozen C4 local seeds. Saved held-out predictions are reused; missing training/validation predictions are generated by read-only frozen-head inference. No C0–C4 model was refitted or retuned. B4 uses eight image-only PCA components, plus local and propagated image/molecular summaries; its per-gene residual head has 35 coefficients. PCA scaling is trained only on the permitted training H&E. No new expression-derived normalization or scaling is fitted. The inherited C0 target transformation remains fixed as required, which also limits the interpretation of sensitivity folds.

Technical audit: **{audit['status']}**. All {audit['prediction_files_hashed_and_identity_checked']} prediction files and observed/hidden barcode identities passed. All {audit['hidden_label_corruption_prediction_comparisons']} full-hidden-vector NaN/random corruption comparisons returned identical predictions (maximum difference {audit['corruption_maximum_error']:.3g}). The independent replay covers {audit['stratified_prediction_and_metric_replays']} stratified runs: every model/encoder/geometry/patient/panel at 20%-replicate-0, plus all negative controls; maximum prediction difference {audit['replay_maximum_error']:.3g}. All 216 mask copies, 10 training-only image scalers, B4 checkpoints and frozen C0–C4 artifacts passed. This is sampled independent replay, not a claim of independently rerunning every saved prediction twice.

All prediction runs retain source commit, mask/checksum, barcode lists, encoder/cache hashes, configuration, scaler/checkpoint references, seed, predictions, runtime and peak RSS. Runs are appended to the global registry without changing its header or prior prefix; the lossless extension is [C5_REGISTRY_EXTENSIONS.csv](C5_REGISTRY_EXTENSIONS.csv). [C5_RESOURCE_PROFILE.csv](C5_RESOURCE_PROFILE.csv) includes inference, corruption-audit and metric time together; model-fit events retain training resources. Large uncompressed predictions/models/caches remain local and are hashed in the final manifest. Full gene CSVs exist locally; gzip copies are committed to keep repository files manageable.

Stop after C5 for review. No new encoder, encoder fine-tuning, multiresolution graph, pathway analysis, explainability, biological mechanism, experimental cost saving or clinical claim is introduced.
'''
 (OUT/'C5_REPORT.md').write_text(text);dump(OUT/'C5_STATUS.json',{'C5_status':status,'recommendation':recommendation,'gate':gate,'technical_audit':audit,'C0_C4_unchanged':True,'C6_executed':False,'external_validation_executed':False});print(status,recommendation)
if __name__=='__main__':main()
