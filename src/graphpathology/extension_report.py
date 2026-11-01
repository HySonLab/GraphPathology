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
  for baseline in ['matched_local','matched_mean','matched_gcn','matched_sage']:
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
  manifest.append({k:r.get(k) for k in ['variant','seed','status','parameters','selected_epoch','epochs_completed','elapsed_seconds','peak_rss_mb','best_checkpoint_relative','best_checkpoint_sha256','protocol_commit']})
 pd.DataFrame(manifest).to_csv(out/'FIT_INVENTORY.csv',index=False)
 report=metadata(root,status='COMPLETE' if not missing and len(items)==78 else 'PARTIAL',models=len(items),new_fits=len(fits),patients=len(patients),evaluation_cases=len(frames)-len(list(folder.glob('*/MOLECULAR.csv'))),route_rows=len(runs),exact_label_corruption_checks=checks,missing_cases=missing,
  biological_units='patient; mask/training seeds averaged before inference; cohorts never pooled',output_hashes={p.name:sha(p) for p in out.glob('*.csv')})
 (out/'STATUS.json').write_text(json.dumps(report,indent=2)+'\n');return report
