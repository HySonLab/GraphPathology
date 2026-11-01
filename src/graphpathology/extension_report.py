"""Patient-level descriptive reports; no architecture selection or pooled-study inference."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .assets import sha
from .journal_statistics import METRICS,contrast,summarize_paired
from .extension_evaluate import model_inventory
from .extension_train import metadata

KEYS=['cohort','patient','design','condition','budget','model_id','origin','variant','training_seed','route']

def run(root,require_complete=True):
 base=root/'experiments/journal_revision/transformer_egnn';folder=base/'evaluation';out=base/'summary';out.mkdir(exist_ok=True)
 items=model_inventory(root);patients=['TENX99','1160920F','CID4290','CID44971']
 recovered=base/'GEO_SELECTIVE_ACQUISITION.json'
 if recovered.exists() and json.loads(recovered.read_text())['status']=='COMPLETE':patients+=sorted({r['patient_id'] for r in json.loads(recovered.read_text())['selected_files']})
 missing=[(p,i['id']) for p in patients for i in items if not (folder/p/i['id']/'COMPLETE.json').exists()]
 if require_complete and missing:raise ValueError(f'Incomplete evaluation: {len(missing)} model/patient cases remain')
 frames=[];checks=0
 for p in sorted(folder.glob('*/*/COMPLETE.json')):
  record=json.loads(p.read_text())
  for name,digest in record['output_hashes'].items():assert sha(p.parent/name)==digest
  checks+=record['exact_label_corruption_checks'];frames.append(pd.read_csv(p.parent/'PANEL.csv.gz'))
 for p in sorted(folder.glob('*/MOLECULAR.csv')):frames.append(pd.read_csv(p))
 runs=pd.concat(frames,ignore_index=True)
 counts=runs.groupby(KEYS,dropna=False).agg(min_defined_pearson=('n_defined_pearson','min'),mean_defined_pearson=('n_defined_pearson','mean'),max_defined_pearson=('n_defined_pearson','max'),panel_genes=('n_genes','first')).reset_index();counts.to_csv(out/'VALID_CORRELATION_COUNTS.csv',index=False)
 secondary=runs.groupby(KEYS+['geometry'],dropna=False)[['mean_spearman','median_spearman']].mean().groupby(KEYS,dropna=False).mean().groupby([k for k in KEYS if k not in ['training_seed','model_id']],dropna=False).mean().reset_index();secondary.to_csv(out/'SECONDARY_SPEARMAN_PATIENT_METRICS.csv',index=False)
 # Mask seeds averaged equally within geometry, then geometries equally within patient; training seeds remain explicit.
 geometry=runs.groupby(KEYS+['geometry'],dropna=False)[METRICS+['n_defined_pearson']].mean().reset_index()
 trained=geometry.groupby(KEYS,dropna=False)[METRICS+['n_defined_pearson']].mean().reset_index()
 keys=[k for k in KEYS if k not in ['training_seed','model_id']]
 patient=trained.groupby(keys,dropna=False)[METRICS+['n_defined_pearson']].mean().reset_index()
 training_sd=trained.groupby(keys,dropna=False)[METRICS].std(ddof=1).reset_index();training_sd.to_csv(out/'TRAINING_SEED_SD.csv',index=False)
 patient.to_csv(out/'PATIENT_METRICS.csv',index=False);trained.to_csv(out/'TRAINING_SEED_PATIENT_METRICS.csv',index=False)
 paired=[];pairkeys=['cohort','patient','design','condition','budget']
 for (origin,variant,route),candidate in patient[patient.origin!='molecular'].groupby(['origin','variant','route']):
  ref=patient[patient.variant=='idw_k16_p2'];p=contrast(candidate,ref,pairkeys,variant+'__'+route+'__vs_IDW','extension_'+origin);p['candidate']=variant;p['route']=route;p['reference']='idw_k16_p2';p['information_match']='same anchors + images vs anchors' if route!='image' else 'unequal: images vs measured anchors';paired.append(p)
  controls=['gps_local','gps_no_relative','gps_self','gps_shuffled','gps_degree_length'] if variant=='gps_full' else (['distance_fixed','distance_disabled','egnn_self','egnn_shuffled','egnn_degree_length'] if variant=='egnn_latent' else [])
  for baseline in ['matched_local','matched_mean','matched_gcn','matched_sage']+controls:
   ref=patient[(patient.origin=='new')&(patient.variant==baseline)&(patient.route==route)]
   p=contrast(candidate,ref,pairkeys,variant+'__'+route+'__vs_'+baseline,'extension_'+origin);p['candidate']=variant;p['route']=route;p['reference']=baseline;p['information_match']='same route and anchors';paired.append(p)
 pairs=pd.concat(paired,ignore_index=True);pairs.to_csv(out/'PAIRED_PATIENT_EFFECTS.csv',index=False)
 summary=summarize_paired(pairs,['origin','cohort','design','condition','budget','candidate','route','reference','information_match','contrast']);summary.to_csv(out/'PATIENT_UNCERTAINTY.csv',index=False)
 # Region-exclusion minus distributed compares exactly the same hidden H and anchor count.
 b2=[]
 for (origin,variant,route),g in patient[patient.design=='B2_identical_H'].groupby(['origin','variant','route']):
  p=contrast(g[g.condition=='region_exclusion'],g[g.condition=='distributed'],['cohort','patient','design','budget'],variant+'__'+route+'__region_exclusion_minus_distributed','extension_B2');p['candidate']=variant;p['route']=route;p['model_origin']=origin;b2.append(p)
 b2=pd.concat(b2,ignore_index=True);b2.to_csv(out/'B2_PAIRED_PATIENT_EFFECTS.csv',index=False)
 summarize_paired(b2,['origin','cohort','design','budget','candidate','route','model_origin','contrast']).to_csv(out/'B2_PATIENT_UNCERTAINTY.csv',index=False)
 # All image-only predictions are mask-independent; H, not input, varies in historical comparisons.
 image=runs[(runs.design=='B2_identical_H')&(runs.route=='image')]
 audit=image.groupby(['patient','model_id'])[METRICS].agg(lambda x:float(x.max()-x.min()))
 assert (audit.fillna(0).to_numpy()==0).all();audit.to_csv(out/'B2_IMAGE_MASK_INDEPENDENCE.csv')
 fits=[json.loads(p.read_text()) for p in sorted((base/'fits').glob('*.json')) if not p.name.endswith('_FAILED.json')]
 manifest=[]
 for r in fits:
  manifest.append({k:r.get(k) for k in ['variant','seed','status','parameters','selected_epoch','validation_mse','epochs_ran','elapsed_seconds','peak_rss_mb','best_checkpoint_relative','best_checkpoint_sha256','protocol_commit']})
 pd.DataFrame(manifest).to_csv(out/'FIT_INVENTORY.csv',index=False)
 matrix=[{'model_id':i['id'],'origin':i['origin'],'variant':i['variant'],'seed':i['seed'],'status':'COMPLETE' if all((folder/p/i['id']/'COMPLETE.json').exists() for p in patients) else 'BLOCKED_INPUTS','expected_patients':len(patients),'completed_patients':sum((folder/p/i['id']/'COMPLETE.json').exists() for p in patients)} for i in items]
 matrix += [{'model_id':name,'origin':origin,'variant':name,'seed':-1,'status':status,'expected_patients':0,'completed_patients':0} for name,origin,status in [('anchor_conditioned_learning','optional','SKIPPED_PROTOCOL'),('full_image_equivariance','optional','SKIPPED_PROTOCOL'),('original_selected_primary_heads','historical','BLOCKED_INPUTS')]]
 pd.DataFrame(matrix).to_csv(out/'EXPERIMENT_MATRIX.csv',index=False)
 report=metadata(root,status='COMPLETE' if not missing and len(items)==78 else 'PARTIAL',models=len(items),new_fits=len(fits),patients=len(patients),evaluation_cases=len(frames)-len(list(folder.glob('*/MOLECULAR.csv'))),route_rows=len(runs),exact_label_corruption_checks=checks,missing_cases=missing,
  biological_units='patient; mask/training seeds averaged before inference; cohorts never pooled',output_hashes={p.name:sha(p) for p in out.glob('*.csv')})
 scientific_report(root,patient,pairs,summary,pd.DataFrame(manifest))
 (out/'STATUS.json').write_text(json.dumps(report,indent=2)+'\n');return report


def scientific_report(root,patient,pairs,uncertainty,fits):
 from .extension_paper import NAMES
 from .extension_models import VARIANTS
 def md(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
 lines=['# Post-v1 GraphTransformer/EGNN scientific report','',
 'All16 configurations and seeds0/1/2 were fixed in local protocol commit2e03fae636d6f6527d25db58f4cbab72f7c1a286 before fitting. This is exploratory work in previously inspected cohorts, not formal preregistration or fresh blind validation. No external label selected an architecture, epoch, seed or calibration coefficient.','',
 'Exact GEO byte ranges recovered all88 Bassiouni input files with original SHA256 matches (291,569,893 bytes transferred; no37.7GB whole-series transfer). The matched A/B1/B2 study completed19,800 runs and39,600 exact all-unobserved-label corruption checks. At20% B2 coverage, all22 patients deteriorated under region exclusion for both fixed estimators: IDW median paired Δr=-0.154080 and relative MSE=+44.9864%; uniform Δr=-0.150195 and relative MSE=+39.9128%. H and anchor counts are identical between placements. See [full baseline source tables](bassiouni_anchors/PATIENT_UNCERTAINTY.csv).','',
 'The neural matrix has48 new fits plus30 frozen journal-v2 fits evaluated in26 patients (development1, Wu3, Bassiouni22), using image, fixed0.2/0.8 and exploratory residual routes. All routes use the same targets/masks; only the image-only versus molecular comparison is unequal-information. Mask means precede geometry/patient means, then training seeds are averaged within patient. Relative MSE is normalized after these averages. Wu patients are shown individually below; only Bassiouni uses10,000 patient bootstrap resamples (seed20261101). No cohort is pooled.','',
 'The primary neural results are negative. At20% Bassiouni V1 coverage, GPS full residual reconstruction changed median paired correlation by−0.050421 and MSE by+13.4264%; EGNN latent changed correlation by−0.047994 and MSE by+13.0083%. Neither primary model improved Pearson/MSE/MAE in any of the22 patients. Nonlocal GPS versus its local ablation and latent EGNN versus fixed distance also worsened the primary point estimates. This does not support a stronger neural novelty claim under the tested conditions.','','## Validation and resources','',md(['Configuration','Parameters','Selected val MSE mean ± SD','Selected epoch range','Fit wall s median','Peak process MiB max'],[[NAMES[v],int((q:=fits[fits.variant==v]).parameters.iloc[0]),f'{q.validation_mse.mean():.5f} ± {q.validation_mse.std():.5f}',f'{q.selected_epoch.min()}–{q.selected_epoch.max()}',f'{q.elapsed_seconds.median():.2f}',f'{q.peak_rss_mb.max():.1f}'] for v in VARIANTS]),'',
 'Validation uses TENX95 only. All new models use two full-section source updates per epoch; older30 reference fits retain their different original minibatch regime. Parameter counts are nominal, not equal effective capacity. Training wall time includes concurrent workload and is not a controlled speed ranking. [Standardized CPU inference](diagnostics/CPU_BENCHMARK.csv) reports all seeds, four threads, two warmups and five timed forwards on the same7,760-node validation section, excluding frozen CNN/graph setup. No GPU or paid service was used.','',
 '## All architecture test metrics at20% historical V1 coverage','',
 'Each reported value is a patient summary after all three mask/training seeds; Bassiouni rows are medians across22 such summaries. Wu remains individual. Full10/20/40%/central-V2/B2/per-gene exports are retained, including undefined-correlation counts and secondary near-tie-sensitive Spearman. No model is selected from these tables.']
 groups=[('Development','TENX99'),('Wu','1160920F'),('Wu','CID4290'),('Wu','CID44971'),('Bassiouni',None)]
 for cohort,who in groups:
  subset=patient[(patient.cohort==cohort)&(patient.design=='historical_V1')&(patient.budget==.2)]
  if who:subset=subset[subset.patient==who]
  for route in ['image','fixed_20pct_image','exploratory_residual_anchor']:
   rows=[]
   for v in list(VARIANTS)+['idw_k16_p2','uniform_k16']:
    q=subset[(subset.variant==v)&((subset.route==route)|(subset.origin=='molecular'))]
    rows.append([NAMES.get(v,v),f'{q.median_pearson.median():.5f}',f'{q.mean_mse.median():.6f}',f'{q.mean_mae.median():.6f}',f'{q.n_defined_pearson.min():.1f}–{q.n_defined_pearson.max():.1f}'])
   lines += ['',f'### {cohort} / {who or "22 patient medians"} / {route}','',md(['Configuration','Pearson','MSE','MAE','Mean defined genes range'],rows)]
 lines+=['','## Matched patient effects and component controls','',
 'The following Bassiouni residual-route contrasts use the same20% V1 anchors and scored support. Effects are median paired patient differences; brackets are descriptive95% percentile intervals, unadjusted for multiple contrasts. They are not p-values or simultaneous bands. Every model/control/route/budget and each Wu patient effect is in the [paired exports](summary/PAIRED_PATIENT_EFFECTS.csv) and [uncertainty table](summary/PATIENT_UNCERTAINTY.csv).']
 rows=[]
 for v in VARIANTS:
  references=['idw_k16_p2','matched_local','matched_gcn','matched_sage']+(['gps_local','gps_no_relative','gps_self','gps_shuffled','gps_degree_length'] if v=='gps_full' else ['distance_fixed','distance_disabled','egnn_self','egnn_shuffled','egnn_degree_length'] if v=='egnn_latent' else [])
  for ref in references:
   q=uncertainty[(uncertainty.cohort=='Bassiouni')&(uncertainty.origin=='extension_new')&(uncertainty.candidate==v)&(uncertainty.reference==ref)&(uncertainty.route=='exploratory_residual_anchor')&(uncertainty.design=='historical_V1')&(uncertainty.condition=='historical')&(uncertainty.budget==.2)].set_index('effect')
   values=[]
   for metric,scale in [('delta_pearson',1),('relative_mse',100),('delta_mae',1)]:
    r=q.loc[metric];values.append(f'{scale*r.median_effect:+.4f} [{scale*r.ci_low:+.4f}, {scale*r.ci_high:+.4f}]')
   rows.append([NAMES[v],NAMES.get(ref,ref)]+values)
 lines+=['',md(['Candidate','Reference','Δr [interval]','Relative MSE% [interval]','ΔMAE [interval]'],rows),'',
 '## Geometric/long-range controls and limitations','',
 'All384 checkpoint-level rotation/reflection/translation/permutation checks passed with explicit precision-dependent tolerances. They are conditional on fixed embeddings and documented graph/landmark handling; no end-to-end CNN equivariance is claimed. The latent EGNN really updates private coordinates; its final coordinate update does not affect scalar prediction, and both the fixed-distance and latent models are E(2)-compatible. The latent model has128 extra nominal coefficients (last64 inactive for scalar loss). Full/local GPS has the same121,018 nominal parameters but changes key count and connection structure as well as range. Self/no-distance controls change effective function class; no-distance retains implicit graph/landmark geometry.','',
 '[Attention/reach and distant-feature diagnostics](diagnostics/) report physical versus assigned connectivity and attention mass beyond250/600/1200um. A remote perturbation tests dependency, not learned molecular correlation. Degree/length controls retain actual accepted-swap counts; zero/low swaps are limited evidence. Low-resolution external images (nominal112um native crop plus published encoder center crop), generic frozen features, source training in only two patients and transported calibration constrain morphology claims.','',
 '## Main text, supplement and novelty','',
 'The22-patient identical-H/equal-count placement experiment is the clearest new main-text evidence: it separates anchor placement from estimator choice and changed scored support. Keep all20 estimators, all neural routes/seeds, individual Wu effects, distance supports, resources and geometric checks fully available in the supplement/source tables. Keep a compact negative neural paired-effect/control summary in the main text to avoid hiding the failed hypotheses; move the full matrices and diagnostics to the supplement; it must retain errors alongside correlation and cannot establish biological long-range learning from attention reach.','',
 'The implemented architectures follow existing GraphGPS/EGNN ideas. This study supports a rigorous information-controlled reconstruction benchmark and conditional software properties, not a new equivariance theorem, general state-of-the-art claim, clinical benefit or demonstrated total assay savings. Any stronger geometry-specific neural novelty claim must be supported by the full ablation/error/patient matrix and future independent native-resolution cohorts.','',
 '## Preservation, statuses and author bundles','',
 'All407 protected prior files and74,253 archived files remain unchanged. The original29 publication targets/71 checks and810 frozen-model development primary-metric rows pass; new raw Bassiouni references are separately audited. The original selected primary heads remain missing (BLOCKED_INPUTS), and the78 post-v1 states are not replacements. Optional anchor-conditioned learning and full-image equivariance have SKIPPED_PROTOCOL status. All completed analyses remain visible regardless of direction. The [experiment matrix](summary/EXPERIMENT_MATRIX.csv) separates these statuses from negative findings.','',
 'Raw data/model binaries are excluded from Git. Full per-gene exports and48 new best/resume/optimizer states are preserved in author-local ZIP bundles with a [release-ready manifest](RELEASE_READY_MANIFEST.json). No public release upload or new credential exposure was performed; the manifest provides exact author paths/checksums for an approved release destination. The existing30-state bundle remains separate. PR3 stays draft and unmerged for author review.']
 path=root/'experiments/journal_revision/transformer_egnn/SCIENTIFIC_REPORT.md';path.write_text('\n'.join(lines)+'\n')
