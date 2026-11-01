"""C4B scientific diagnostics and final gated C4 report."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=OUT/'C4B'
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def fmt(v):return f'{v.mean():.3f} ± {v.std(ddof=1):.3f}'
def save(fig,name):fig.tight_layout();fig.savefig(B/'figures'/name,dpi=170,bbox_inches='tight');plt.close(fig)
def main():
 gate=json.loads((B/'logs/GATE_RESULT.json').read_text());audit=json.loads((B/'logs/INTEGRITY_AUDIT.json').read_text());aaudit=json.loads((A/'logs/INTEGRITY_AUDIT.json').read_text());status=gate['status'] if audit['status']=='PASS' else 'FAIL_TECHNICAL';p=pd.read_csv(B/'C4B_PANEL_METRICS.csv');g=pd.read_csv(B/'C4B_GENE_METRICS.csv');gain=pd.read_csv(B/'C4B_GENE_GAIN.csv');support=pd.read_csv(B/'C4B_GENE_SUPPORT.csv');comparison=pd.read_csv(B/'C4B_VS_FROZEN_C3.csv');panels=['S2-50','S2-280'];controls=['LOCAL','MEAN','TRUE','MATCHED'];rows=[]
 for panel in panels:
  for control in controls:
   d=p[(p.phase=='primary')&(p.panel==panel)&(p.control==control)]
   if control=='MATCHED':d=d.groupby('graph_seed')[['median_pearson','median_spearman','mean_mse','mean_mae']].mean()
   rows.append([panel,control,fmt(d.median_pearson),fmt(d.median_spearman),fmt(d.mean_mse),fmt(d.mean_mae)])
 contrast=[]
 for (phase,panel),d in comparison.groupby(['phase','panel']):
  # Each patient has the same gene panel. Summarize fold medians equally across patients.
  grouped=d.groupby('sample_id')[['C3_delta_spatial','delta_spatial','C3_delta_matched','delta_matched']].median().mean()
  contrast.append([phase,panel,*[f'{grouped[c]:+.3f}' for c in ['C3_delta_spatial','delta_spatial','C3_delta_matched','delta_matched']]])
 folds=[]
 for (panel,sample),d in gain[gain.phase=='lopo'].groupby(['panel','sample_id']):folds.append([panel,sample,f'{d.delta_spatial.median():+.3f}',f'{d.delta_matched.median():+.3f}',f'{d.delta_MSE_spatial.mean():+.3f}',f'{d.delta_MSE_matched.mean():+.3f}',f'{d.mse_TRUE.mean()/d.mse_LOCAL.mean()-1:+.1%}',f'{d.mse_TRUE.mean()/d.mse_MATCHED.mean()-1:+.1%}'])
 gates=[[r['panel'],r['robust_genes'],r['required_genes'],r['primary_candidates'],f"{r['median_spatial_gene_gain']:+.3f}",f"{r['median_matched_gene_gain']:+.3f}",f"{r['mean_spatial_MSE_reduction']:+.3f}",f"{r['mean_matched_MSE_reduction']:+.3f}"] for r in gate['panels']]
 best=support.sort_values('delta_matched',ascending=False).head(15);genes=[[r.panel,r.gene,f'{r.delta_spatial:+.3f}',f'{r.delta_mean:+.3f}',f'{r.delta_matched:+.3f}',f'{r.delta_MSE_matched:+.3f}',r.supporting_patients,str(r.robust_gene)] for r in best.itertuples()]
 graph_rows=[]
 for panel in panels:
  true=g[(g.phase=='primary')&(g.panel==panel)&(g.control=='TRUE')].groupby('gene')[['pearson','mse']].mean()
  for gs,d in g[(g.phase=='primary')&(g.panel==panel)&(g.control=='MATCHED')].groupby('graph_seed'):
   d=d.groupby('gene')[['pearson','mse']].mean();graph_rows.append([panel,int(gs),f'{(true.pearson-d.pearson).median():+.4f}',f'{(d.mse-true.mse).mean():+.4f}'])
 # Scientific figures: primary endpoints, per-gene contrasts, patient dependence.
 fig,axes=plt.subplots(2,2,figsize=(12,8))
 for j,panel in enumerate(panels):
  for ax,metric in zip(axes[:,j],['median_pearson','mean_mse']):
   for i,c in enumerate(controls):
    d=p[(p.phase=='primary')&(p.panel==panel)&(p.control==c)]
    if c=='MATCHED':d=d.groupby('graph_seed')[[metric]].mean()
    v=d[metric];ax.bar(i,v.mean(),color=['#888888','#82ad92','#267b60','#8a75b4'][i],alpha=.75);ax.scatter(i+np.linspace(-.14,.14,len(v)),v,s=22,c='black',zorder=3)
   ax.set_xticks(range(4),['Local MLP','Mean + MLP','True SAGE','Matched SAGE']);ax.set_title(panel);ax.set_ylabel('Median gene Pearson' if metric=='median_pearson' else 'Mean gene MSE')
 fig.suptitle('Phikon spatial recheck: training seeds; matched dots average three seeds per frozen graph');save(fig,'primary_comparison.png')
 fig,axes=plt.subplots(1,2,figsize=(12,5))
 for ax,panel in zip(axes,panels):
  d=support[support.panel==panel];ax.scatter(d.delta_matched,d.delta_MSE_matched,c=np.where(d.robust_gene,'#d89b30','#267b60'),s=20,alpha=.7);ax.axhline(0,color='black',lw=.8);ax.axvline(0,color='black',lw=.8);ax.set_title(panel);ax.set_xlabel('True − distance-matched gene Pearson');ax.set_ylabel('Matched MSE − true MSE')
 fig.suptitle('Primary per-gene specificity; upper-right improves both endpoints; gold = robust rule');save(fig,'matched_gene_gain.png')
 fig,axes=plt.subplots(1,2,figsize=(12,4))
 for ax,metric,title in zip(axes,['delta_matched','delta_MSE_matched'],['Median true − matched gene Pearson','Mean matched MSE − true MSE']):
  d=gain[gain.phase=='lopo'].pivot_table(index='sample_id',columns='panel',values=metric,aggfunc='median' if metric=='delta_matched' else 'mean').reindex(index=IDS,columns=panels);v=d.to_numpy();bound=max(abs(v).max(),1e-6);im=ax.imshow(v,vmin=-bound,vmax=bound,cmap='RdBu',aspect='auto');ax.set_yticks(range(4),IDS);ax.set_xticks(range(2),panels);ax.set_title(title)
  for i in range(4):
   for j in range(2):ax.text(j,i,f'{v[i,j]:+.3f}',ha='center',va='center',color='white' if abs(v[i,j])/bound>.5 else 'black')
  fig.colorbar(im,ax=ax)
 fig.suptitle('Conditional patient sensitivity; positive favors true topology');save(fig,'patient_sensitivity.png')
 recommendation={'PASS_SPATIAL_RECHECK':'RECONSIDER_MULTISCALE_GRAPH','MIXED_SPATIAL_RECHECK':'NEED_REPRESENTATION_REVIEW','FAIL_SPATIAL_RECHECK':'PROCEED_TO_SPARSE_ST_RECONSTRUCTION','FAIL_TECHNICAL':'STOP_TECHNICAL'}[status]
 verdict={'PASS_SPATIAL_RECHECK':'The predeclared reproducible-subset criterion passes. This permits later reconsideration of multiresolution graphs, but no such experiment was run here.','MIXED_SPATIAL_RECHECK':'Some primary genes meet the candidate rule, but the full reproducible spatial-specificity criterion does not pass. Stop for review.','FAIL_SPATIAL_RECHECK':'The predeclared reproducible spatial-specificity criterion does not pass, and the primary candidate fallback is not reached. Deprioritize multiresolution graphs for this cohort.','FAIL_TECHNICAL':'Technical integrity failed; scientific interpretation is blocked.'}[status]
 report=f'''C4B status: {status}

# Minimal spatial-specificity recheck

{verdict}

C4A authorized this stage through **S2-50 only**: +0.026 median paired gene Pearson gain, 15 strong genes, approximately 0.9% lower primary mean MSE, and two conditional patients improving both endpoints. The 250/280-gene primary results worsened; this is not a broad representation improvement. S2-280 remains the prespecified stress test here.

## Matched spatial experiment

The node features are the exact C4A Phikon caches. LOCAL reuses the 30 saved C4A MLP outcomes for these two panels. MEAN uses the exact C2 self-inclusive neighbor-mean operator. TRUE uses the frozen C2 one-hop mean GraphSAGE operator, concatenating local and neighbor means before the learned MLP. MATCHED uses the identical architecture/configuration and exact C3 A6 adjacency matrices. All nonlocal models use three training seeds; MATCHED crosses those with all five frozen graph seeds. No encoder, radius, graph, architecture or hyperparameter search was added.

Frozen geometry: S2-50 uses capped Delaunay for both MEAN and SAGE. S2-280 uses the frozen radius-150 mean baseline and kNN-6 capped-at-250 SAGE geometry. No new graph construction is introduced. SAGE hidden widths/dropout remain 67/0.1 and 156/0.1 for 50/280 genes; MEAN widths/dropout remain 128/0.1 and 256/0.1. Fixed sensitivity epochs remain C2's 7/13 for SAGE and 14/17 for MEAN. Primary validation-MSE early stopping remains the original C1 protocol.

Spatial hyperparameters and schedules match C3, but input dimension changes parameter count. The local comparator is C4A's newly selected matched-budget MLP, which can differ from C3's frozen C1 local configuration. Thus the direct C3 contrast is informative development evidence, not a claim that every downstream parameter or local optimization setting was identical. Feature and target scalers are reused from the audited C4A training-only fits. Measured neighboring expression never enters the inputs.

The C4B source and decision rules were locked before 42 primary fits, and all 42 model records were locked before B test prediction. C4A and older test results were already known. This remains a sequential development experiment. Four-patient results use primary-frozen choices and fold-training scalers without held-out early stopping; they are conditional sensitivity, not unbiased cross-validation.

## Primary correlation and error

{table(['Panel','Model','Median gene Pearson','Median Spearman','Mean MSE','Mean MAE'],rows)}

For LOCAL/MEAN/TRUE, values are mean ± sample SD over three training seeds. For MATCHED, average the three training seeds within each graph first, then report mean ± sample SD over five graph realizations. All 15 individual matched outcomes remain in the metric tables. No population significance or independent-patient uncertainty is implied. Errors remain in frozen C0 normalized units.

![Primary endpoints](figures/primary_comparison.png)

## Gene-specific gains and frozen decision

{table(['Panel','Robust genes','Required','Primary candidates','Median spatial gain','Median matched gain','Mean spatial MSE reduction','Mean matched MSE reduction'],gates)}

A robust gene must gain at least 0.02 mean primary Pearson over LOCAL, MEAN and MATCHED, improve against LOCAL and MATCHED in every training seed, exceed one matched-graph SD, and have non-increasing primary MSE versus LOCAL/MATCHED. At least two conditional patients must have positive LOCAL/MATCHED correlation gain with MSE worsening at most 2% against both; mean fold-relative MSE worsening versus MATCHED must be at most 5%. A panel needs `max(3, ceil(2% of panel size))` robust genes. The MIXED fallback requires at least three primary genes with matched gain ≥0.02, positive local gain and no MSE worsening versus either. These descriptive criteria were frozen before B outcomes.

The largest primary matched gains below are shown regardless of final robustness. All genes, including negative gains, are retained in [C4B_GENE_GAIN.csv](C4B_GENE_GAIN.csv) and [C4B_GENE_SUPPORT.csv](C4B_GENE_SUPPORT.csv).

{table(['Panel','Gene','Gain vs local','Gain vs mean','Gain vs matched','Matched MSE reduction','Supporting patients / 4','Robust'],genes)}

![Matched gene-specific gains](figures/matched_gene_gain.png)

## Patient dependence

{table(['Panel','Held-out sample','Median local gain','Median matched gain','Mean local MSE reduction','Mean matched MSE reduction','Relative MSE increase vs local','Relative MSE increase vs matched'],folds)}

![Patient sensitivity](figures/patient_sensitivity.png)

## Graph sensitivity

All 56 referenced graph files were copied exactly from frozen C2/C3 artifacts; 40 are matched graphs (two SAGE geometries × four slides × five seeds). Their per-node degrees and global physical-distance-bin counts remain exact. Each matched graph has at least 10% edge turnover, length-distribution KS ≤0.2, and absolute relative mean edge-length change ≤10%, using the frozen C3 audit. This constrained randomization is not claimed uniform.

Each row below averages training seeds before calculating the gene-paired primary contrast against one graph realization. Detailed graph statistics/hashes are under `graphs/`.

{table(['Panel','Frozen graph seed','Median true-minus-matched gene Pearson','Mean matched-minus-true MSE'],graph_rows)}

## Direct comparison with frozen C3 ResNet18 results

{table(['Phase','Panel','C3 median spatial gain','C4B median spatial gain','C3 median matched gain','C4B median matched gain'],contrast)}

Sensitivity entries are equal-weight means of four patient-specific gene medians. Primary entries refer to TENX99. [C4B_VS_FROZEN_C3.csv](C4B_VS_FROZEN_C3.csv) preserves every gene/patient contrast and corresponding MSE changes. C3 results were read from frozen metrics only; no C0–C3 model was rerun or retuned.

## Integrity and recommendation

Audit: **{audit['status']}**. All {audit['saved_models_replayed']} new models and held-out predictions replayed, plus {audit['primary_validation_replays']} primary validation outputs; maximum prediction difference {audit['maximum_prediction_error']:.3g}. Reused local prediction files are byte-identical; metric CSV round-tripping differs by at most {audit['reused_metric_CSV_roundtrip_max_error']:.3g}. All 30 reused local outcomes, 56 graph references, all frozen C4A artifacts and all C0–C3 hashes passed. The registry has {audit['registry_events']} C4B events appended without changing its prior header/prefix. Per-fit parameter counts, timings and peak RAM are in [C4B_RESOURCE_PROFILE.csv](C4B_RESOURCE_PROFILE.csv).

A CSV parser correction was necessary: pandas inferred the per-run control label TRUE as Boolean. The typed reader adapter restores text labels before executing the unchanged locked analysis source. Saved metrics/predictions and scientific decision rules were not changed; the correction and regression test are logged.

Final recommendation: **{recommendation}**. Stop after C4 for user review. No multiresolution, partial-ST reconstruction, interpretation, explainability or clinical claims were executed.
'''
 (B/'C4B_REPORT.md').write_text(report)
 ap=pd.read_csv(A/'C4A_PANEL_METRICS.csv');arows=[]
 for panel in ['S1','S2-50','S2-250','S2-280']:
  d=ap[(ap.phase=='primary')&(ap.model=='D1')].groupby(['panel','encoder'])[['median_pearson','mean_mse']].mean();r=d.loc[(panel,'resnet18_reference')];h=d.loc[(panel,'pathology_encoder_1')];arows.append([panel,f'{r.median_pearson:.3f}',f'{h.median_pearson:.3f}',f'{r.mean_mse:.3f}',f'{h.mean_mse:.3f}'])
 encgenes=pd.read_csv(A/'C4A_GENE_SUPPORT.csv');encgenes=encgenes[(encgenes.phase=='primary')&(encgenes.model=='D1')&(encgenes.panel=='S2-50')&encgenes.strong_gene].sort_values('delta_encoder',ascending=False).head(5);leading_genes=', '.join(f"{r.gene} ({r.delta_encoder:+.3f})" for r in encgenes.itertuples())
 root=f'''C4A status: PASS_PATHOLOGY_REPRESENTATION
C4B status: {status}
Overall C4 status: {status}

# C4 — Representation and minimal spatial recheck

## 1. Does a pathology-specific encoder materially improve H&E→expression prediction?

**A limited pass on S2-50, not a broad improvement across panels.** The single public frozen Phikon fallback improves the 50-gene MLP median correlation from 0.690 to 0.710, with approximately 0.9% lower MSE. Median gene-paired gain is +0.026; 15/50 genes satisfy the strong-gain rule, 62% improve Pearson and 36% improve both Pearson and MSE. The locked criterion passes on this panel and two conditional patients support it. S1 also improves, while the 250- and 280-gene primary outcomes worsen substantially.

{table(['Panel','ResNet MLP median Pearson','Phikon MLP median Pearson','ResNet MLP mean MSE','Phikon MLP mean MSE'],arows)}

These are means over three training seeds on the primary patient. Full Ridge/MLP metrics, individual seeds, gene gains and controls are in [C4A_REPORT.md](C4A/C4A_REPORT.md). UNI2-h and UNI were gated without configured access, so Phikon was the single planned public fallback. The C4A report is the frozen snapshot written before B; the final B status is shown above. Both encoders received the same compact downstream grid; no extra encoder was added. The conclusion applies to this encoder pipeline and cohort, not all pathology foundation models.

## 2. Which panels and genes benefit most?

S1 and S2-50 benefit most on the primary test patient. Fifteen S2-50 genes meet the encoder strong-gain criterion. The largest primary Pearson gains among these are **{leading_genes}**. Their identities, full per-gene correlation/error changes and fold support are in [C4A_GENE_SUPPORT.csv](C4A/C4A_GENE_SUPPORT.csv) and the detailed C4A report. Spatial gene-specific gains and failure cases are in [C4B_GENE_SUPPORT.csv](C4B/C4B_GENE_SUPPORT.csv). No biological significance, pathway or clinical interpretation is assigned to these identities.

{table(['Panel','Spatial robust genes','Required','Primary candidates','Median spatial gain','Median matched gain','Mean spatial MSE reduction','Mean matched MSE reduction'],gates)}

## 3. Are gains reproducible across conditional held-out patients?

C4A's S2-50 gain meets the joint correlation/error condition in two of four patients, with important error deterioration in NCBI785. It is patient dependent. All four sensitivity folds preserve frozen panels/normalization and primary-selected configurations, fit C4 scalers using fold-training patients only, and avoid held-out early stopping. These are conditional analyses, not external validation or unbiased cross-validation. Full stage-specific patient tables are retained. The following table reports **C4B spatial contrasts**, not encoder gains; positive reductions favor true topology.

{table(['Panel','Held-out sample','Median local gain','Median matched gain','Mean local MSE reduction','Mean matched MSE reduction','Relative MSE increase vs local','Relative MSE increase vs matched'],folds)}

## 4. Does true tissue topology become more informative?

**{status}.** {verdict} The recheck used only LOCAL, neighbor MEAN, true SAGE and the strongest frozen distance-matched C3 control on S2-50/S2-280. It used three training seeds and five matched-graph seeds, exact C2/C3 single-scale geometry, and no new graph search.

{table(['Phase','Panel','C3 median spatial gain','C4B median spatial gain','C3 median matched gain','C4B median matched gain'],contrast)}

See [C4B_REPORT.md](C4B/C4B_REPORT.md) for graph sensitivity, gene-level criteria, correlation/error metrics and replay evidence. C3 rows are exact frozen metrics, not rerun models. Spatial model settings remain frozen C2/C3 choices; the local C4A probe can have different selected settings than the C3 local baseline, and input dimensions change parameter counts.

![Spatial primary results](C4B/figures/primary_comparison.png)

## 5. What should happen next?

**{recommendation}**

Stop here for review. No later checkpoint was executed. A full CMRCNet reproduction is deferred because the official preprocessing/panels/training and CUDA assumptions do not support a cheap faithful comparison in this frozen CPU setup; see [CROSS_MODAL_COMPATIBILITY.md](CROSS_MODAL_COMPATIBILITY.md). No approximate implementation is labeled CMRCNet.

## Limits and reproducibility

C0–C3 remain unchanged. C4A replayed 244 stored models with zero prediction discrepancy, checked all 35,536 RGB patch hashes, eight caches and 80 training-only scaler fits, and exactly replayed 64 sampled encoder rows. C4B replayed 210 new models with maximum prediction discrepancy {audit['maximum_prediction_error']:.3g}, verified 30 reused local outcomes and 56 graph files, and preserved the complete frozen C4A artifact set. The registry contains 420 C4A plus 506 C4B events appended to the existing header/prefix. Code, configurations, reports, metrics and hashes are versioned; large weights/caches/checkpoints/predictions remain local and are hashed in the final manifest.

A logged CSV-reader correction preserved TRUE as a text control label; the locked analysis formulas and thresholds remained unchanged.

The compared encoder pipelines also differ in architecture, feature dimension and effective crop field of view: Phikon uses the full approximately 112 µm source tile; frozen ResNet18 retained approximately 98 µm after its official crop. Public-pretraining slide overlap has not been independently excluded. All Ridge selections reached the maximum predeclared alpha; Phikon MLPs often selected epoch 1. Negative controls pass their frozen average-over-seeds criterion, but some individual label-permutation seeds retain positive correlation (up to approximately 0.42); the individual results are disclosed and are not evidence of null signal in every run. Prior checkpoint/test outcomes were known, thresholds are descriptive, and no statistical, external-validation, causal or clinical claim follows.
'''
 (OUT/'C4_REPORT.md').write_text(root)
 dump(OUT/'C4_STATUS.json',{'C4A_status':'PASS_PATHOLOGY_REPRESENTATION','C4B_status':status,'overall_C4_status':status,'recommendation':recommendation,'C4A_gate':json.loads((A/'logs/GATE_RESULT.json').read_text()),'C4B_gate':gate,'technical_audit':{'status':'PASS' if aaudit['status']==audit['status']=='PASS' else 'FAIL','C4A':aaudit,'C4B':audit},'C0_C3_unchanged':True,'later_checkpoints_executed':False});print(status,recommendation)
if __name__=='__main__':main()
