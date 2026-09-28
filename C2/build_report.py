"""Create C2's evidence report without retuning or changing any experiment."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C2'

def fmt(v,d=3):
 a=np.asarray(v,dtype=float);a=a[np.isfinite(a)]
 if not len(a):return 'undefined'
 return f'{a.mean():.{d}f}'+(f' ± {a.std(ddof=1):.{d}f}' if len(a)>1 else '')
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])

def main():
 p=pd.read_csv(OUT/'C2_METRICS_SUMMARY.csv');g=pd.read_csv(OUT/'C2_CONTEXTUAL_GAIN_SUMMARY.csv');c=pd.read_csv(OUT/'C2_GRAPH_CONTROLS.csv');candidate=pd.read_csv(OUT/'C2_GENE_CANDIDATE_AUDIT.csv');supported=candidate[candidate.supported_candidate]
 local=pd.read_csv(ROOT/'C1/C1_METRICS_SUMMARY.csv');local=local[local.model=='B2'];runs=pd.read_csv(OUT/'C2_RESOURCE_PROFILE.csv')
 lock=json.loads((OUT/'SELECTION_LOCK.json').read_text());audit=json.loads((OUT/'logs/INTEGRITY_AUDIT.json').read_text());assert audit['status']=='PASS'
 primary_rows=[]
 for panel in ['S1','S2-50','S2-250','S2-280']:
  f=local[local.panel==panel];primary_rows.append([panel,'Frozen local B2',fmt(f.median_pearson),'0 (reference)',fmt(f.median_spearman),fmt(f.mean_mse),fmt(f.mean_mae)])
  for model in ['MEAN','GCN','SAGE']:
   f=p[(p.panel==panel)&(p.model==model)];h=g[(g.phase=='primary')&(g.panel==panel)&(g.model==model)]
   primary_rows.append([panel,model,fmt(f.median_pearson),fmt(h.median_delta_local),fmt(f.median_spearman),fmt(f.mean_mse),fmt(f.mean_mae)])
 gain_rows=[]
 for panel in ['S2-50','S2-250','S2-280']:
  for model in ['MEAN','GCN','SAGE']:
   h=g[(g.phase=='primary')&(g.panel==panel)&(g.model==model)]
   gain_rows.append([panel,model,fmt(h.mean_delta_local),fmt(h.median_delta_vs_mean),fmt(h.mean_mse_change),fmt(100*h.fraction_genes_positive,1)+'%',fmt(100*h.fraction_genes_gain_ge_005,1)+'%'])
 candidate_rows=[]
 for r in supported.sort_values(['panel','gene']).itertuples():candidate_rows.append([r.panel,r.gene,f'{r.primary_delta_local:+.3f}',f'{r.primary_delta_mean:+.3f}',f'{r.primary_mse_change:+.3f}',f'{r.lopo_positive_both_patients}/4',f'{r.lopo_mean_delta_local:+.3f}',f'{r.lopo_mean_mse_change:+.3f}'])
 patient_rows=[]
 for panel in ['S2-50','S2-250','S2-280']:
  for sample in ['NCBI783','NCBI785','TENX95','TENX99']:
   h=g[(g.phase=='lopo')&(g.panel==panel)&(g.model=='SAGE')&(g.sample_id==sample)]
   patient_rows.append([panel,sample,fmt(h.median_delta_local),fmt(h.median_delta_vs_mean),fmt(h.mean_mse_change)])
 control_rows=[]
 for model in ['MEAN','GCN','SAGE']:
  for phase,frame in [('Real spatial graph',p[p.panel=='S2-50']),('Vertex-permuted graph',c)]:
   f=frame[frame.model==model];control_rows.append([model,phase,fmt(f.median_pearson),fmt(f.mean_mse),fmt(f.mean_mae)])
 config_rows=[]
 for panel,models in lock['selected'].items():
  for model,choice in models.items():config_rows.append([panel,model,choice['graph']['name'],choice['params']['hidden'],choice['params']['dropout'],choice['lopo_epochs']])
 resources=[]
 for model in ['MEAN','GCN','SAGE']:
  r=runs[(runs.phase=='primary')&(runs.model==model)];resources.append([model,fmt(r.groupby('seed').train_seconds.sum()),fmt(r.groupby('seed').prepare_seconds.sum()),fmt(r.groupby('seed').peak_rss_mb.max(),1),fmt(r.groupby('seed').inference_seconds.sum())])
 unique=', '.join(sorted(supported.gene.unique()));error_degrade=int((supported.lopo_mean_mse_change>0).sum())
 status='MIXED_CONTEXT_SIGNAL'
 report=f'''# C2 report: does single-scale spatial context add information?

**Scientific outcome: {status}. Technical audit: PASS.** C2 is complete. C0 and accepted C1 remain unchanged; no C1 model was retrained, replayed or retuned. C3 has not been executed.

GraphSAGE gives modest primary-test median-correlation gains on the 250- and 280-gene settings, but not on S2-50. Per-gene gains are heterogeneous and patient-dependent. The frozen descriptive shortlist/support rule identifies **{len(supported)} gene-panel results covering {supported.gene.nunique()} distinct genes** ({unique}). These are candidates for further context-specificity checks, not established biological findings. Error disagreement and the shuffled-GraphSAGE control prevent a blanket claim that GNNs add reliable spatial information beyond local morphology and neighborhood averaging.

## Controlled design

The frozen C0 cohort, ResNet18 embeddings, split, panels and target normalization are unchanged. The primary split remains NCBI783/NCBI785 training, TENX95 validation, TENX99 test. All four settings were run: S1/FOXA1, S2-50, S2-250 and S2-280 (full shared assay panel).

The comparison is frozen local C1 B2 -> spatial MEAN + MLP -> one-layer GCN / mean-GraphSAGE. MEAN and GCN retain the corresponding C1 hidden width, dropout, direct output head, optimizer and training settings. GraphSAGE concatenates self and neighbor-mean features and reduces hidden width to match parameter count; the largest mismatch is 0.282%. Thus parameter budgets are controlled, but GraphSAGE's local/self pathway and width still differ architecturally from the frozen local baseline.

- MEAN: row-normalize A+I, aggregate image embeddings, then the C1-selected MLP.
- GCN: symmetric degree normalization of A+I, a learned linear graph convolution, GELU, frozen dropout, direct readout.
- SAGE: concatenate the local embedding with neighbor-only mean, learned linear map, GELU, frozen dropout, direct readout; isolated nodes receive a zero neighbor vector.

Precomputing aggregation of frozen embeddings is exactly equivalent to a single linear message-passing layer before the activation; a numerical equivalence test verifies this. These are shallow one-layer models, not multi-hop or multiresolution networks. The normalization/root-feature conventions follow the standard [GCN operator](https://github.com/pyg-team/pytorch_geometric/blob/master/torch_geometric/nn/conv/gcn_conv.py) and [mean GraphSAGE operator](https://pytorch-geometric.readthedocs.io/en/2.8.0/generated/torch_geometric.nn.conv.SAGEConv.html); implementation uses SciPy sparse operations and PyTorch directly.

Four predefined, undirected, slide-local graphs were compared: capped kNN k=6 (250 µm edge limit), 150 µm radius, 250 µm radius, and Delaunay with edges longer than 250 µm removed. Graphs use coordinates only. They never contain targets or cross-slide edges. Each experiment uses one graph scale. C0 physical calibration is estimated, not a biological interaction radius.

For each panel/model, graph choice uses seed-0 validation median gene Pearson, with validation MSE and graph name as tie-breaks. The architecture is not retuned. Final primary models use seeds 0, 1, 2 and the unchanged C1 training/early-stopping settings. The model-selection lock and analysis-rule lock were written before C2 test access. Since C1's test result was already known, this is a sequential development checkpoint, not a fresh untouched confirmatory study.

## Primary results: correlation and error

Values are mean ± sample SD over three seeds on **one test patient**. They are not population confidence intervals. MSE/MAE are in frozen C0 normalized units after inverting the C1 target scaler.

{table(['Setting','Model','Median gene Pearson','Median per-gene Δ vs local','Median gene Spearman','Mean MSE','Mean MAE'],primary_rows)}

**The two correlation columns answer different questions.** The main scientific quantity is Δ_g = r_g(model) − r_g(frozen local), paired by gene/patient/seed. The median of these differences is not generally the difference of median correlations. For example, S2-50 GraphSAGE has a lower aggregate median Pearson than local B2, yet its median paired gene gain is positive; gene ordering changes between predictors. Both summaries are retained.

S1 benefits on the primary test slide, especially from simple averaging. For the broader settings, most incremental changes are small, with both positive and negative gene effects. GCN does not show a consistent advantage over MEAN.

![Primary comparison](figures/primary_comparison.png)

{table(['Panel','Model','Mean per-gene Δ vs local','Median per-gene Δ vs MEAN','Mean MSE change vs local','Genes with Δ > 0','Genes with Δ ≥ 0.05'],gain_rows)}

All per-gene Pearson, Spearman, MSE, MAE and changes are available in [C2_GENE_CONTEXTUAL_GAINS.csv](C2_GENE_CONTEXTUAL_GAINS.csv). Fractions are calculated per seed over the full panel before aggregation; no negative-gain genes were excluded.

## Gene-specific candidate audit

C1 did not retain validation per-gene predictions. To honor the instruction not to rerun C1, the pre-test validation shortlist compares graph models with MEAN only. It requires mean validation gain ≥0.05 over MEAN, a positive gain in every seed, and non-increasing validation MSE. That creates 11 S2-50, 140 S2-250 and 139 S2-280 SAGE gene-panel candidates; no GCN gene meets this shortlist rule.

A candidate is then flagged as supported by the **descriptive** rule frozen before test access if: primary mean gain over local is ≥0.05; each primary seed improves over local; primary mean gain over MEAN is positive; primary mean MSE does not increase; and at least three of four conditional LOPO patients have positive mean gains over both local and MEAN. This rule is not a significance test. The large validation shortlist is not itself evidence of gain beyond local morphology.

{table(['Panel','Gene','Primary Δ vs local','Primary Δ vs MEAN','Primary MSE change','LOPO positive vs both','LOPO mean Δ vs local','LOPO mean MSE change'],candidate_rows)}

These {len(supported)} rows represent only {supported.gene.nunique()} unique genes because LRRC15 occurs in two panels. **{error_degrade} of the {len(supported)} rows worsen average LOPO MSE**, despite correlation support. In this audit, CCDC80 and LRRC15 in S2-280 also improve average LOPO MSE; this still does not prove spatial specificity on those panels. No pathway or biological interpretation is assigned to any gene.

The complete shortlist, rejected candidates, fold-level gains and error changes are in [C2_GENE_CANDIDATE_AUDIT.csv](C2_GENE_CANDIDATE_AUDIT.csv) and [C2_GENE_PATIENT_DEPENDENCE.csv](C2_GENE_PATIENT_DEPENDENCE.csv).

![Gene-level correlation/error changes](figures/S2-280_gene_gains.png)

## Patient dependence

This is a **conditional four-patient robustness sensitivity analysis, not unbiased or external validation**. Each C2 fold reads the exact frozen C1 fold-specific scalers and compares with the exact frozen C1 LOPO baseline for the same patient, panel and seed. Models are fitted only on the other three patients. Graph/configuration choice stays fixed from primary validation, and epochs are fixed to the rounded median primary best epoch; evaluation patients do not drive early stopping.

The following table highlights GraphSAGE. Full MEAN/GCN/SAGE results, including S1 and per-seed values, are in [C2_ROBUSTNESS_LOPO.csv](C2_ROBUSTNESS_LOPO.csv) and [C2_CONTEXTUAL_GAIN_SUMMARY.csv](C2_CONTEXTUAL_GAIN_SUMMARY.csv).

{table(['Panel','Evaluation sample','Median per-gene Δ vs local','Median per-gene Δ vs MEAN','Mean MSE change vs local'],patient_rows)}

NCBI783 shows negative median contextual gains for all three model families across the multi-gene panels. GraphSAGE's broader-panel gains are more favorable on the TENX patients than the NCBI patients, and several folds worsen MSE. Four patients cannot separate patient effects from source/staining/preparation differences. Reused panels, hyperparameters and overlapping fold-training cohorts further limit independence.

![Patient dependence](figures/S2-280_patient_dependence.png)

## Spatial-specificity control

The S2-50 control permutes graph vertex labels independently by slide and seed. It preserves the graph's topology and degree distribution while breaking the link between edges and physical proximity. All features and expression targets remain correctly paired at their own nodes. These controls were specified before test access and use the selected graph geometry; no target values enter the permutation.

{table(['Model','Graph condition','Median gene Pearson','Mean MSE','Mean MAE'],control_rows)}

Mean pooling and GCN deteriorate markedly under permutation: they rely on real neighborhoods to retain predictive morphology. This does **not** by itself show that they outperform local prediction. GraphSAGE retains strong performance with shuffled neighbors, consistent with its ability to use the unaltered local-feature branch. Therefore its S2-50 performance cannot be attributed uniquely to physical neighbors. The shuffled control was not run on S2-250/S2-280, so broader-panel and shortlisted-gene specificity remains unresolved; results from the 50-gene control cannot establish it.

![Shuffled-graph controls](figures/shuffled_graph_controls.png)

## Graph sensitivity and selected configurations

{table(['Panel','Model','Selected graph','Hidden width','Dropout','Fixed LOPO epochs'],config_rows)}

Sensitivity across all four graph choices is recorded on **validation only** in [C2_GRAPH_SENSITIVITY.csv](C2_GRAPH_SENSITIVITY.csv). Test scores were not used to choose a graph or extend the search. Delaunay is selected for most comparisons; S2-280 selects radius150 for MEAN/GCN and capped kNN6 for SAGE. Thus on S2-280, the main comparison includes independently validation-selected geometry as well as model family; this is not a common-geometry architectural ablation.

Fifteen graph files were needed: all four choices on the primary training/validation slides, plus only selected choices on the primary test slide. No unused radius250 test graph was constructed. Graph statistics include edge lengths, degree distribution, isolated nodes and connected components. Pruning removes long edges but does not guarantee every edge avoids a small tissue gap. Delaunay choices on near-grid coordinates can also depend on triangulation tie-breaking; cached graph hashes fix the executed version.

![Validation sensitivity](figures/validation_graph_sensitivity.png)

![Graph construction](figures/graph_construction_audit.png)

## Resources and audit

CPU-only execution used four threads on an Intel Core i7-14700. The following summarizes primary fits across the four settings, per seed. Training time excludes graph/feature preparation, which is shown separately; training RSS includes prepared arrays and resident data. Inference time includes model loading and feature aggregation. No encoder extraction was repeated.

{table(['Model','Training seconds','Preparation seconds','Peak training RSS MiB','Inference seconds'],resources)}

[C2_RESOURCE_PROFILE.csv](C2_RESOURCE_PROFILE.csv) records every fit and inference, including parameter counts and graph-preparation time. [C2_GRAPH_STATISTICS.csv](C2_GRAPH_STATISTICS.csv) records one-time graph-construction costs.

- 237 fits: 48 graph-search candidates, 36 primary fits, 9 shuffled-graph controls and 144 conditional LOPO fits.
- 189 evaluations and 26,595 gene-metric rows; 426 C2 fit/evaluation events appended without changing the C1 registry prefix.
- All 15 graph files reconstructed exactly from the frozen coordinates and configurations.
- All 189 C2 saved predictions replayed within the frozen numerical tolerance; maximum absolute replay error: {audit['max_replay_error']:.3g}.
- Ten C2 unit tests pass, including one-layer GCN equivalence, physical bounds, graph symmetry, determinism, disconnected-region locality and parameter matching.
- C0 and all frozen C1 artifacts pass their hash audit. No C1 training or inference was invoked.
- The first audit script expected an unused sixteenth graph. Its expectation was corrected to the actual validation/selected-test schedule; no data, graph, model or result was regenerated. The initial audit log is retained.

`C2_SPEC_LOCK.json` records source, inputs, selected graphs/configurations, scaler references, library versions, result/model/graph hashes and the append-only registry prefix. Model checkpoints/prediction arrays remain locally under C2 with hashes; source, metrics, graphs, figures, logs and reports are Git-versioned.

## Decision and recommendation before C3

**C2 is technically complete, with mixed and subset-specific evidence. Do not claim a general graph advantage.** The modest broader-panel improvements and six candidate genes warrant focused follow-up, but patient-dependent losses, correlation/error disagreement, and the GraphSAGE shuffled-neighbor result leave spatial attribution unresolved.

Before committing to a C3 multiresolution claim, review this C2 freeze. The most informative next authorized experiment would be a separately versioned control study on S2-250/S2-280 and the fixed candidate list: compare true neighbors with shuffled neighbors and a self-only GraphSAGE control that retains the same width/parameterization and training protocol. Keep graph search and all tuning on validation. Advance candidate-driven multiresolution experiments only if their benefits survive those controls and the error tradeoffs are acceptable.

Important limits remain: four patients; one primary test patient; correlated locations; conditional overlapping LOPO folds; generic ImageNet embeddings; 280 targeted assay genes; estimated spatial calibration; no cell-level alignment certification; one-layer graph models; training-seed variability rather than patient confidence intervals; primary-test knowledge inherited from C1; and no external, biological, mechanistic or clinical validation. No pathology encoder, S3, GAT, multiresolution fusion, explanation or biological interpretation was introduced.

**STOP: C3 has not started. Await review.**
'''
 (OUT/'C2_REPORT.md').write_text(report)
 decision={'status':status,'technical_status':'PASS','supported_gene_panel_results':len(supported),'unique_candidate_genes':sorted(supported.gene.unique()),
           'candidate_rows_with_worse_mean_lopo_mse':error_degrade,'conclusion':'Modest broader-panel gains and a descriptive subset, but no general spatial-GNN advantage established',
           'C3_recommendation':'Review; resolve broader-panel shuffled-neighbor and matched self-only controls in a separate version before a multiresolution claim','C3_executed':False}
 (OUT/'logs/C2_DECISION.json').write_text(json.dumps(decision,indent=2)+'\n');print(status)
if __name__=='__main__':main()
