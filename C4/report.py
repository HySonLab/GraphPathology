"""Report C4A without turning a single-encoder development study into a general claim."""
from common import *

def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def fmt(s):return f'{s.mean():.3f}'+(f' ± {s.std(ddof=1):.3f}' if len(s)>1 else ' (deterministic)')
def main():
 gate=json.loads((A/'logs/GATE_RESULT.json').read_text());audit=json.loads((A/'logs/INTEGRITY_AUDIT.json').read_text());status=gate['status'] if audit['status']=='PASS' else 'FAIL_TECHNICAL';p=pd.read_csv(A/'C4A_PANEL_METRICS.csv');g=pd.read_csv(A/'C4A_GENE_SUPPORT.csv');primary=p[p.phase=='primary'];gain=pd.read_csv(A/'C4A_GAIN_SUMMARY.csv');extract=[json.loads(f.read_text()) for f in (A/'logs/extraction').glob('*.json')]
 rows=[]
 for panel in PANELS:
  for model in ['D0','D1']:
   for encoder in ['resnet18_reference','pathology_encoder_1']:
    d=primary[(primary.panel==panel)&(primary.model==model)&(primary.encoder==encoder)];rows.append([panel,model,'ResNet18' if encoder=='resnet18_reference' else 'Phikon',fmt(d.median_pearson),fmt(d.median_spearman),fmt(d.mean_mse),fmt(d.mean_mae)])
 gainrows=[[r.panel,f'{r.median_gene_gain:+.3f}',f'{r.mean_gene_gain:+.3f}',f'{r.fraction_Pearson_improved:.1%}',f'{r.fraction_both_improved:.1%}',r.strong_genes,r.conditional_patients_improved_in_both,str(r._asdict().get('pass',False))] for r in gain.itertuples(index=False)]
 # pandas renames a keyword column such as pass in named tuples; render that column separately below.
 gainrows=[row[:-1]+[str(bool(gain.iloc[i]['pass']))] for i,row in enumerate(gainrows)]
 foldrows=[]
 for panel in PANELS:
  for sample in IDS:
   d=g[(g.phase=='lopo')&(g.model=='D1')&(g.panel==panel)&(g.sample_id==sample)];foldrows.append([panel,sample,f'{d.delta_encoder.median():+.3f}',f'{d.delta_mse.mean():+.3f}',f'{(d.mse_pathology.mean()/d.mse_resnet.mean()-1):+.1%}'])
 top=g[(g.phase=='primary')&(g.model=='D1')].sort_values('delta_encoder',ascending=False).head(15);tops=[[r.panel,r.gene,f'{r.delta_encoder:+.3f}',f'{r.delta_mse:+.3f}',str(r.strong_gene)] for r in top.itertuples()]
 controls=[[r['model'],r['control'],f"{r['mean_median_Pearson']:.3f}",f"{r['drop']:.3f}",str(r['pass'])] for r in gate['negative_controls']]
 patches=[[r['sample_id'],r['n_locations'],f"{r['estimated_FOV_um']:.2f}",f"{r['estimated_resized_mpp']:.4f}",f"{r['alignment_max_error_px']:.3f}",f"{r['seconds']:.1f}"] for r in extract]
 recommendation='STOP_TECHNICAL' if status=='FAIL_TECHNICAL' else ('PROCEED_TO_SPARSE_ST_RECONSTRUCTION' if status=='FAIL_PATHOLOGY_REPRESENTATION' else 'NEED_REPRESENTATION_REVIEW')
 answer='The predeclared representation gate passes on this development cohort; a minimal C4B spatial recheck is authorized.' if status=='PASS_PATHOLOGY_REPRESENTATION' else 'A clear representation pass was not established under the predeclared gate. C4B is not authorized; stop for review.'
 text=f'''C4A status: {status}
C4B status: NOT_RUN
Overall C4 status: {status}

# C4A — Frozen pathology representation benchmark

## 1. Does a pathology-specific encoder materially improve prediction?

**{answer}** This compares the specific frozen Phikon release with the exact C0 ResNet18 embeddings. It is not a conclusion about all pathology foundation models, and it is not external validation.

{table(['Panel','Model','Encoder','Median gene Pearson','Median Spearman','Mean MSE','Mean MAE'],rows)}

D0 is Ridge; D1 is a one-hidden-layer MLP. MLP values are mean ± sample SD across three training seeds on one primary test patient. Ridge is deterministic and is not padded with duplicate-seed outcomes. The complete tables also include mean/median Pearson, Spearman, MSE and MAE; fractions above 0/0.1/0.2/0.3 Pearson; per-gene target/prediction variance; and secondary coexpression diagnostics. The 280-gene setting is the **full shared assay panel**, not whole-transcriptome prediction. Errors use the unchanged C0 normalized target units.

![Primary comparison](figures/primary_comparison.png)

### Matched downstream design

Both encoders use the same five Ridge alphas `[0.1, 1, 10, 100, 1000]` and the same four MLP candidates: hidden width 128/256 crossed with dropout 0/0.1. Seed-0 validation median Pearson chooses the configuration, with validation MSE and a canonical config tie-break. The unchanged C1 AdamW, learning-rate, weight-decay, batch-size, patience and maximum-epoch settings apply to both representations. Three MLP primary seeds follow selection. This C4 ResNet18 probe is a new matched-budget benchmark; it does not replace or change frozen C1 results.

Feature and target scalers are fitted on the C4 training patients only. Gene lists/order, C0 target normalization and patient split remain unchanged. Model families, hidden widths and tuning budgets match; **absolute parameter counts differ because inputs are 512 versus 768 dimensions**. Ridge is an additional check across model families. [C4A_RESOURCE_PROFILE.csv](C4A_RESOURCE_PROFILE.csv) records the actual counts and training resources.

The specification and gate were hashed before new test evaluation. Only training/validation images were extracted initially; selected downstream configs were locked before extracting the test image cache. The test target cache was then used for final evaluation only. C0–C3 test results were already known, so this is sequential development rather than fresh confirmation.

## 2. Which panels and genes benefit most?

The primary scientific quantities are gene-paired encoder Pearson gain and MSE reduction. Positive values mean Phikon improves. Median paired gene gains need not equal the difference of panel median correlations.

{table(['Panel (MLP)','Median gene gain','Mean gene gain','Pearson improved','Both endpoints improved','Strong genes','Folds improved in both','Panel gate pass'],gainrows)}

A strong descriptive gene requires mean primary Pearson gain ≥0.05, positive gain in every MLP seed and non-increasing mean MSE. This is a development threshold, not biological significance. The following top gains are shown irrespective of whether every rule passes; all genes remain in [C4A_ENCODER_GAIN.csv](C4A_ENCODER_GAIN.csv) and [C4A_GENE_SUPPORT.csv](C4A_GENE_SUPPORT.csv).

{table(['Panel','Gene','Mean Pearson gain','Mean MSE reduction','Strong descriptive gene'],tops)}

![Per-gene gain and error](figures/gene_gain_vs_error.png)

![All-gene gain distributions](figures/gain_distributions.png)

No pathways, mechanisms or clinical conclusions are assigned to these gene identities.

## 3. Are gains reproducible across conditional held-out patients?

All four patients and all four panels were evaluated with both models and both encoders. Configurations remain primary-selected. MLP epochs are fixed per encoder/panel to the rounded median primary best epoch; the held-out fold never early-stops a model. Fold feature/target scaling uses only the other three patients. These folds overlap in training data and inherit primary-selected panels/hyperparameters; they are **conditional sensitivity analyses, not unbiased cross-validation**.

{table(['Panel (MLP)','Held-out sample','Median paired Pearson gain','Mean MSE reduction','Relative MSE increase'],foldrows)}

![Patient sensitivity](figures/patient_sensitivity.png)

### Pairing negative controls

Only training image/expression pairings are permuted; validation/test stay intact. NC1 permutes training labels, NC2 permutes training feature locations. Each uses three seeds for both D0 and D1 on Phikon S2-50.

{table(['Model','Control','Mean median Pearson','Drop from paired','Pass'],controls)}

![Negative controls](figures/negative_controls.png)

Each control must have mean median Pearson ≤0.20 and lose at least 0.15 versus the corresponding unpermuted model. Failing this integrity gate overrides the scientific status with FAIL_TECHNICAL.

## Encoder, patch compatibility and resources

UNI2-h and UNI returned HTTP 401 gated-weight responses with no configured access. Under the plan's fallback order, one public CPU-practical pathology model was used: [Phikon](https://huggingface.co/owkin/phikon/tree/057cc0295895c2df3dd7681a89680da6015cbefe), revision `057cc0295895c2df3dd7681a89680da6015cbefe`. No second encoder or outcome-driven encoder sweep was run. The published feature interface is a 768-dimensional CLS representation from a frozen ViT-B/16; C4 follows the official extraction convention without the unused pooler. [Official extraction instructions](https://github.com/owkin/HistoSSLscaling#update--phikon-release-on-hugging-face-).

The exact 224×224 RGB C0 pixel arrays are reused with the official image processor: scaling by 1/255 and ImageNet channel means/standard deviations; no stain normalization, color augmentation or recropping. Approximate physical calibration supports the usual 20× interpretation, but scanner-level calibration is not independently certified. Source patch fields of view and alignment errors are:

{table(['Sample','Locations','Source FOV µm','Resized µm/pixel','Max alignment error px','Phikon extraction seconds'],patches)}

**Preprocessing is part of the compared encoder pipeline.** Phikon consumes the whole 224-pixel source tile. The frozen ResNet18 official transform had resized it to 256 pixels then center-cropped 224, retaining about 87.5% of the source width (approximately 98 µm rather than 112 µm). Thus this benchmark does not isolate pretraining from architecture, feature dimension or effective crop field of view. No C0 patch/cache was altered to conceal that difference. Public-pretraining overlap has not been independently excluded at individual-slide level.

Extraction used CPU float32, batch size 16, eight PyTorch threads, on Intel Core i7-14700; no functioning NVIDIA driver was available. Actual extraction time sums to {sum(r['seconds'] for r in extract):.1f} seconds. Downstream training uses four CPU threads. The isolated `C4/deps` environment adds pinned Transformers/Hugging Face libraries without changing the C0–C3 environment. Model weight, processor/config, patch, cache, barcode-order, scaler and checkpoint hashes are retained. [C4A_ENCODER_MANIFEST.csv](C4A_ENCODER_MANIFEST.csv) records cache dimensions, bytes and extraction resources.

E0 reference cache files are wrappers containing bit-identical C0 embeddings plus the new provenance fields; ResNet18 was not extracted again. E1 embeddings were extracted once per location and reused. Every location retains sample, patient, barcode, coordinates, source HDF5 row, source RGB hash and encoder/preprocessing hash. Both cache directories are under `C4/embeddings/` and referenced from C4A manifests. Large caches, weights, models and predictions remain local with hashes and are excluded from Git.

## Frozen gate and technical audit

A multi-gene MLP panel passes only if its panel median correlation improves ≥0.02, median paired gene gain improves ≥0.02, mean MSE worsens no more than 2%, at least 20% of genes improve both endpoints, and at least `max(3, ceil(5% of panel size))` genes meet the strong-gene rule. At least two conditional patients must improve median gene correlation with MSE worsening no more than 2%; mean fold-relative MSE worsening must be ≤5%. The same-panel Ridge median gene gain must be nonnegative. Negative controls and technical integrity must pass. MIXED/FAIL fallbacks are specified in [C4A_SPEC_LOCK.json](C4A_SPEC_LOCK.json); thresholds were not changed after test outcomes.

Technical status: **{audit['status']}**. {audit['saved_models_replayed']} stored downstream models were replayed, including {audit['selection_metric_replays']} validation candidates and {audit['evaluation_predictions_replayed']} held-out outcomes; maximum prediction difference {audit['maximum_prediction_error']:.3g}. {audit['patch_hashes_verified']} source RGB patch hashes, {audit['embedding_caches_checked']} caches and {audit['scaler_training_only_checks']} training-only scaler checks passed. A sampled replay of {audit['sampled_embedding_rows']} encoder rows had maximum difference {audit['sampled_embedding_max_error']:.3g}. This sampled encoder replay is not a claim of re-extracting the full cohort twice. All frozen C0–C3 hashes are unchanged.

All extraction/fitting/evaluation events were appended to the existing registry. Its frozen header/prefix remains unchanged; additional encoder fields are in `C4A_REGISTRY_EXTENSIONS.csv`. Source/config locks precede test access; the final artifact manifest is separate and does not rewrite that history.

## 4. Does true topology become more informative?

C4B is **NOT_RUN** in this C4A report. Only a clear PASS_PATHOLOGY_REPRESENTATION authorizes the minimal spatial recheck. No new radius, multiresolution model, adaptive scale fusion, explanation or sparse-ST input was introduced.

## 5. What should happen next?

Current recommendation: **{recommendation}**. {'A clear representation pass authorizes C4B before assigning the final C4 recommendation.' if status=='PASS_PATHOLOGY_REPRESENTATION' else 'C4A did not clearly pass, so stop here for review. Do not add encoders or proceed to later checkpoints automatically.'}

A faithful CMRCNet auxiliary benchmark is deferred for the concrete compatibility reasons in [CROSS_MODAL_COMPATIBILITY.md](../CROSS_MODAL_COMPATIBILITY.md). No approximate CMRCNet result is reported. This checkpoint supports no causal, clinical, external-validation or publication-level claims.
'''
 (A/'C4A_REPORT.md').write_text(text)
 roottext=text.replace('](C4A_','](C4A/C4A_').replace('](figures/','](C4A/figures/').replace('](../CROSS_MODAL_COMPATIBILITY.md)','](CROSS_MODAL_COMPATIBILITY.md)')
 (OUT/'C4_REPORT.md').write_text(roottext)
 dump(OUT/'C4_STATUS.json',{'C4A_status':status,'C4B_status':'NOT_RUN','overall_C4_status':status,'recommendation':recommendation,'C4B_allowed':status=='PASS_PATHOLOGY_REPRESENTATION','C4A_gate':gate,'technical_audit':audit});print(status)
if __name__=='__main__':main()
