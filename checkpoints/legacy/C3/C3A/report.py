"""Report the locked C3A gate, including null results and control limitations."""
from common import *

def table(headers,rows):
 return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def pm(s):return f'{s.mean():.3f} ± {s.std(ddof=1):.3f}'

def main():
 p=pd.read_csv(OUT/'C3A_PANEL_METRICS.csv');g=pd.read_csv(OUT/'C3A_GENE_SPATIAL_GAINS.csv');c=pd.read_csv(OUT/'C3A_CANDIDATE_AUDIT.csv');a=pd.read_csv(OUT/'C3A_GRAPH_AUDIT.csv');gate=json.loads((OUT/'logs/GATE_RESULT.json').read_text());audit=json.loads((OUT/'logs/INTEGRITY_AUDIT.json').read_text());status=gate['status'] if audit['status']=='PASS' else 'FAIL_TECHNICAL'
 primary=p[p.phase=='primary'];pg=g[g.phase=='primary'];fold=g[g.phase=='lopo'];candidates=c[c.candidate];robust=c[c.robust]
 rows=[]
 for panel in PANELS:
  for control in ['A0','A1','A2','A3','A4','A5','A6']:
   d=primary[(primary.panel==panel)&(primary.control==control)]
   if control in ['A4','A5','A6']:d=d.groupby('graph_randomization_seed')[['median_pearson','mean_mse','mean_mae','median_spearman']].mean()
   rows.append([panel,control,pm(d.median_pearson),pm(d.median_spearman),pm(d.mean_mse),pm(d.mean_mae)])
 gains=[]
 for panel in PANELS:
  d=pg[pg.panel==panel];gains.append([panel]+[f'{d[f"delta_{control}"].median():+.4f}' for control in ['A0','A1','A2','A4','A5','A6']]+[f'{d.mse_change_A6.mean():+.4f}'])
 patients=[]
 for panel in PANELS:
  for sample in IDS:
   d=fold[(fold.panel==panel)&(fold.sample_id==sample)];patients.append([panel,sample,f'{d.delta_A6.median():+.4f}',f'{d.mse_change_A6.mean():+.4f}',f'{(d.delta_A6>0).mean():.1%}'])
 top=c.sort_values('delta_A6',ascending=False).head(12)
 gene_rows=[[r.panel,r.gene,f'{r.delta_A0:+.3f}',f'{r.delta_A5:+.3f}',f'{r.delta_A6:+.3f}',f'{r.mse_change_A6:+.3f}',f'{r.lopo_positive_A6_patients}/4',str(r.candidate),str(r.robust)] for r in top.itertuples()]
 random=a[a.graph_control_type.isin(['A4','A5','A6'])];distance=random[random.graph_control_type=='A6'];graphrows=[]
 for control in ['A4','A5','A6']:
  d=random[random.graph_control_type==control];graphrows.append([control,len(d),f'{d.changed_edge_fraction.min():.1%}–{d.changed_edge_fraction.max():.1%}',str(bool(d.degree_distribution_preserved.all())),str(bool(d.degree_preserved.all())),str(bool(d.distance_distribution_preserved.all())),f'{d.edge_length_ks_rounded_0_001_um.max():.3f}'])
 candidate_text=', '.join(sorted(candidates.gene.unique())) or 'None'
 candidate_detail=(f"Their primary distance-control gains range from {candidates.delta_A6.min():+.4f} to {candidates.delta_A6.max():+.4f}. Counts by panel are "+', '.join(f'{panel}: {len(candidates[candidates.panel==panel])}' for panel in PANELS)+". These small effects are sensitive to the predeclared development tolerances; the 0.02 threshold is not a biological boundary." if len(candidates) else 'No genes meet the descriptive candidate criteria.')
 justification='YES' if status=='PASS_SPATIAL_SPECIFICITY' else ('NOT YET' if status=='MIXED_SPATIAL_SPECIFICITY' else 'NO')
 # A C3A pass authorizes a scale experiment, not a multiresolution architecture.
 architecture='NOT YET' if status in ['PASS_SPATIAL_SPECIFICITY','MIXED_SPATIAL_SPECIFICITY'] else 'NO'
 q1=('Yes, for the robust subset under this development protocol; this is not external or population-level validation.' if status=='PASS_SPATIAL_SPECIFICITY' else 'A reproducible spatial-specificity advantage was not established under the predeclared gate. Positive individual-gene differences must be read alongside the matched-control distributions, error changes and patient sensitivity below.')
 q2=('There is a supported subset under the frozen rule.' if len(robust) else 'Not reproducibly under the frozen rule: no gene-panel result clears the full control, uncertainty, error and patient-consistency requirements.')
 gate_text=('The predefined spatial-specificity gate passed. C3B is the next authorized stage; this report alone does not justify a multiresolution architecture.' if status=='PASS_SPATIAL_SPECIFICITY' else 'The predefined spatial-specificity gate did not pass. Stop here: C3B was not created or executed. No multiresolution model is justified by this checkpoint.')
 text=f'''C3A status: {status}
C3B status: NOT_RUN
Overall C3 status: {status}

# C3A — Spatial specificity

{gate_text}

Technical audit: **{audit['status']}**. All frozen C0/C1/C2 hashes are unchanged. No C1/C2 model was refitted, retuned or replayed. Their saved predictions and exact metrics were imported as immutable references.

## Question 1 — Does true tissue topology outperform matched randomized topology?

**{q1}**

The table below reports the primary-patient results for all controls. A0 is frozen local MLP; A1 is self-only SAGE; A2 is frozen mean pooling + MLP; A3 is frozen true-spatial SAGE; A4 permutes vertex labels; A5 rewires edges with exact node degrees; A6 rewires edges with exact node degrees and physical-length-bin counts. Every random control uses **five graph seeds crossed with three training seeds**.

{table(['Panel','Control','Median gene Pearson','Median Spearman','Mean MSE','Mean MAE'],rows)}

For A0–A3, ± is sample SD across three training seeds. For A4–A6, each graph seed first averages the three training seeds, then ± is sample SD over five graph seeds. These SDs describe different sources of development variability; neither is a patient confidence interval. [C3A_PANEL_METRICS.csv](C3A_PANEL_METRICS.csv) retains all 15 outcomes for each randomized primary comparison, including mean/median Pearson, Spearman, MSE, MAE, threshold fractions, and secondary coexpression. The 280-gene setting is the **full shared assay panel**, not whole-transcriptome prediction.

The central endpoint is the paired per-gene difference, not subtraction of panel medians. The following medians average training seeds and, for controls, graph seeds within each gene before summarizing genes:

{table(['Panel','Δ vs local A0','Δ vs self A1','Δ vs mean A2','Δ vs shuffle A4','Δ vs degree A5','Δ vs distance A6','Mean MSE change vs A6'],gains)}

All signed gains and error changes remain in [C3A_GENE_SPATIAL_GAINS.csv](C3A_GENE_SPATIAL_GAINS.csv); negative effects are retained. [C3A_RANDOMIZATION_METRICS.csv](C3A_RANDOMIZATION_METRICS.csv) contains every graph-seed/gene result, and [C3A_RANDOMIZATION_UNCERTAINTY.csv](C3A_RANDOMIZATION_UNCERTAINTY.csv) reports the corresponding mean and sample SD. Per-training-seed contrasts are separate in [C3A_SEED_CONTRASTS.csv](C3A_SEED_CONTRASTS.csv).

![All primary controls](figures/panel_controls.png)

![True minus shuffled](figures/gain_distribution_shuffle.png)

![True minus degree-matched](figures/gain_distribution_degree.png)

![True minus distance-matched](figures/gain_distribution_distance.png)

### Patient dependence and calibration

The same four conditional held-out patients were evaluated for A1 and A6 on **all four panels**, with all five A6 graph seeds and three training seeds. A0/A2/A3 retain the exact frozen C1/C2 conditional-fold results. A4/A5 are primary-only; they were not expanded to the folds. The strongest spatially matched control, A6, was specified before outcomes rather than chosen post hoc from test performance.

{table(['Panel','Held-out sample','Median gene Δ vs A6','Mean MSE change vs A6','Genes with positive Δ vs A6'],patients)}

These are **conditional patient-sensitivity analyses, not unbiased cross-validation**. Fold feature/target scalers are the frozen C1 objects and were checked numerically against the training patients only. Every new fold model uses the frozen C2 SAGE fixed epoch count, architecture and optimizer; it never early-stops on its held-out patient. Primary-selected panels/configurations and overlapping fold training cohorts limit independence. The primary patient is already part of the prior development history.

![Patient-resolved comparison](figures/patient_comparison.png)

![Correlation versus error](figures/correlation_vs_error.png)

## Question 2 — Is any benefit more than simple averaging or self-only prediction?

**{q2}**

A pass requires a meaningful subset that exceeds **each** of self-only SAGE, frozen mean pooling, fully shuffled, degree-preserving and distance-matched controls, not just frozen local MLP. There are **{len(robust)} robust gene-panel rows covering {robust.gene.nunique()} unique genes** under that rule. Inspect the A1/A2 gains above alongside the graph-null gains; gains relative to A0 alone do not demonstrate spatial specificity.

A1 feeds `[x_i, x_i]` into the exact C2 SAGE dimensions/hidden width/dropout. Both branches have trainable weights and only local information; nominal and active parameter counts match A3 exactly. The duplicated branches remain functionally redundant, so identical parameter count does not imply identical function-space dimension or regularization. A2 uses the exact frozen C2 self-inclusive mean formulation and selected geometry. In S2-280 its radius150 graph differs from A3's capped kNN6, so that comparison includes validation-selected graph geometry. SAGE itself is one linear neighbor aggregation before GELU and readout; it is also expressible as an MLP on local-plus-neighbor features, not evidence that a deep message-passing network is necessary.

## Question 3 — Which genes have the strongest descriptive true-spatial benefit?

The frozen candidate criteria identify **{len(candidates)} gene-panel rows and {candidates.gene.nunique()} distinct genes**: {candidate_text}. These are descriptive development candidates, not population-significant or biological findings. {candidate_detail}

The table ranks the strongest primary mean A6 correlation gains, whether or not they qualify, to expose rejection by error or patient criteria:

{table(['Panel','Gene','Δ local','Δ degree','Δ distance','MSE change vs distance','Positive patient directions','Candidate','Robust'],gene_rows)}

[C3A_CANDIDATE_AUDIT.csv](C3A_CANDIDATE_AUDIT.csv) records every gene, each criterion, every rejection and both candidate/robust flags. The six C2 candidate identities were retained in the specification for traceability, but C3A does not assume they are spatially specific. Candidate spatial maps, when candidates exist, are saved as `figures/map_<panel>_<gene>.png`; maps share a color scale between measured, local, true-graph and distance-control predictions and carry no mechanistic interpretation. Prediction maps average the available training/graph seeds; these diagnostic ensemble maps are not used to recompute the reported per-run metrics.

![All gene spatial gains](figures/gene_spatial_gains.png)

### Exact frozen decision rule

[C3A_SPEC_LOCK.json](C3A_SPEC_LOCK.json) was written and hashed before new C3 test evaluation; the analysis implementation was separately hashed before test access. Candidate: positive true Pearson; positive mean gain over A0/A5/A6; primary MSE no more than 2% worse than A0; positive A6 correlation direction in at least 3/4 conditional patients; no patient A6 gain below −0.05; and no patient MSE more than 25% worse than local. These are explicit pragmatic development tolerances, not universal biological cutoffs.

Robust: candidate plus primary mean gain ≥0.02 over **each A1/A2/A4/A5/A6**; positive A6 gain in all three training seeds; mean A6 gain greater than one graph-randomization SD; no primary MSE increase versus A0 or A6; and no mean conditional-fold MSE increase versus A6. PASS requires at least three distinct robust genes and one multi-gene panel with at least `max(3, ceil(2% of panel size))` robust genes. If no pass, MIXED requires at least one candidate or at least three distinct genes with primary A6 gain ≥0.02 and non-increasing A6 MSE. Otherwise FAIL. Any failed technical/control audit overrides these with FAIL_TECHNICAL. The thresholds were not revised after inspecting C3 predictions.

## Graph controls audit and limits

{table(['Control','Graphs','Edges changed','Degree histogram exact','Per-node degree exact','Length-bin counts exact','Max length KS'],graphrows)}

All 120 random graphs use fixed seeds derived from sample, geometry, control and randomization seed. Undirected graphs have no self/duplicate edges and never cross slides. A4 preserves graph topology under vertex relabeling and the degree histogram. A5 performs double-edge swaps while preserving every node degree. A6 additionally preserves counts in `[0,175)`, `[175,250)`, `[250,∞)` µm. Each proposal is label-free; no target expression enters graph construction or prediction inputs.

The frozen C0 spatial calibration is estimated. Coordinates form an approximately 100 µm lattice, with immediate cardinal/diagonal distances of 100/141 µm. A coordinate-only pilot found that bins separating those distances allowed less than 2% edge turnover under exact degree constraints; the broader immediate-neighborhood bin was chosen before C3 test evaluation. **A6 preserves binned lengths, not every exact length**. Within-bin empirical length KS (rounded to 0.001 µm to suppress floating-point lattice artifacts) and mean-length shifts are audited against predeclared limits of 0.20 and 10%; edge turnover must exceed 10%. A4/A5 must replace at least 90% of edges.

A6 changes only {distance.changed_edge_fraction.min():.1%}–{distance.changed_edge_fraction.max():.1%} of true edges. This is an intentionally constrained null and still retains most true adjacency. It tests whether exact adjacency matters beyond nearby connections of similar lengths; it cannot establish absence of all spatial information. Swap chains are reproducible but are not claimed to sample uniformly from every admissible graph or prove full mixing. Fixed maximum attempts are reported even when the target accepted-swap count is not reached. Stronger randomization would relax the matching constraints and answer a different question.

[C3A_GRAPH_AUDIT.csv](C3A_GRAPH_AUDIT.csv) includes every graph's nodes, edges, degree summaries, components, physical-distance distribution, overlap, algorithm, attempts, accepted swaps, seed, source and content hash. Per-graph JSON files retain detailed provenance. A7, the optional feature-nearest graph, was omitted to focus on the required matched controls.

![Graph diagnostics](figures/graph_diagnostics.png)

![S2-280 graph diagnostics](figures/graph_diagnostics_knn6.png)

The degree histograms overlap exactly because their distributions are preserved.

## Reproducibility and execution

- {audit['new_models_replayed']} new fits and saved-checkpoint prediction replays; maximum absolute replay difference {audit['maximum_absolute_prediction_error']:.3g}.
- {audit['evaluations']} evaluated/reused model outcomes: 192 new primary fits, 288 new conditional fits, and 180 immutable reference outcomes.
- All {audit['graphs_reconstructed']} stored graph variants reconstructed exactly; 120 are randomized controls.
- {audit['frozen_scaler_training_only_checks']} scaler mean/scale checks against the correct training patients; no refitting used for prediction.
- A test/audit module initially resolved the frozen C2 runner because both runner files share a name. Its import was made explicit; the retained initial test log documents this harness error. The training runner, predictions, and pre-test scientific specification did not change.
- Eight unit tests cover randomization constraints/determinism, label-free graph inputs, self-only locality, capacity matching and the balanced seed schedule.
- Frozen C0/C1/C2 hash audits pass; zero C1/C2 model runs. New preprocessing, encoder extraction, S3 inputs, explainability, enrichment and biological interpretation were not introduced.
- {audit['registry_events']} events are appended to the global registry. Its original header/prefix remains byte-identical because prior checkpoint hashes cover it. [C3A_REGISTRY_EXTENSIONS.csv](C3A_REGISTRY_EXTENSIONS.csv) supplies the new schema fields keyed by `(event, run_id)`; multi-slide fit events also retain their per-slide graph references.
- The pre-test specification is immutable. The separate final manifest hashes sources, outputs, graph files, local model/prediction artifacts and the registry prefix without rewriting the pre-test lock.

Models and prediction arrays remain under C3A locally with hashes and are excluded from Git because of size. Code, graph files, configs, logs, metrics, figures, and reports are versioned. Scripts reference the frozen C1/C2 implementations read-only. Primary controls use unchanged C1 optimizer/early-stopping settings; no C3 hyperparameter search was performed.

## Question 4 — Do genes show reproducible physical-scale preferences?

C3B is **NOT_RUN** at this stage. No physical-scale preferences are claimed and no C3B directory is created unless the C3A gate passes.

## Question 5 — Is a gene-adaptive multiresolution model justified?

**{architecture}.** {gate_text} A C3A pass would authorize only the C3B scale-identification experiment; a later gene-adaptive architecture also requires a C3B pass. This checkpoint supports no causal, pathway, external-generalization, clinical or whole-transcriptome claim.
'''
 (OUT/'C3A_REPORT.md').write_text(text);(ROOT/'C3/C3_REPORT.md').write_text(text.replace('](C3A_','](C3A/C3A_').replace('](figures/','](C3A/figures/'))
 dump(ROOT/'C3/C3_STATUS.json',{'C3A_status':status,'C3B_status':'NOT_RUN','overall_C3_status':status,'multiresolution_justified':architecture,'C3B_authorized_by_gate':status=='PASS_SPATIAL_SPECIFICITY','technical_audit':audit,'gate':gate,'stop_for_review':status!='PASS_SPATIAL_SPECIFICITY'})
 print(status)
if __name__=='__main__':main()
