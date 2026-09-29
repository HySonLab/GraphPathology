"""Whole-location leakage, provenance, training-only scalers and stratified replay audit."""
from common import *
import torch
from threadpoolctl import threadpool_limits
from engine import *


def main():
 start=time.perf_counter();lock=check_lock();selection=json.loads((OUT/'SELECTION_LOCK.json').read_text());access=json.loads((OUT/'logs/TEST_ACCESS.json').read_text());assert lock['timestamp']<selection['timestamp']<access['timestamp'];assert selection['spec_sha256']==sha(OUT/'C5_SPEC_LOCK.json');assert access['selection_sha256']==sha(OUT/'SELECTION_LOCK.json');assert not selection['primary_hidden_labels_used'];reg=ROOT/'results/experiment_registry.csv'
 with reg.open('rb') as f:assert hashlib.sha256(f.read(lock['registry_prefix']['bytes'])).hexdigest()==lock['registry_prefix']['sha256']
 index={r['mask_id']:r for r in json.loads((OUT/'configs/masks.json').read_text())}
 for sample in IDS:
  with np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as target:
   for encoder in ENCODERS:
    d=image_data(encoder,sample);assert np.array_equal(target['barcodes'],d['barcodes']);assert np.array_equal(target['xy_um'],d['xy'])

 for r in index.values():
  assert sha(ROOT/r['path'])==r['checksum']==sha(ROOT/r['C0_source_path']);obs=mask_array(r);assert r['num_observed']==int(obs.sum()) and r['realized_coverage']==float(obs.mean())
  with np.load(ROOT/r['path']) as f:assert np.array_equal(f['barcodes'],image_data(ENCODERS[0],r['sample_id'])['barcodes'])
 means=0
 for path in (OUT/'models').glob('B0_*.npz'):
  cohort,panel=path.stem[len('B0_'):].split('_',1);rs=[next(r for r in masks(s) if r['regime']=='random' and r['nominal_coverage']==.6 and r['replicate']==0) for s in training(cohort)];expected=np.concatenate([anchors(r['sample_id'],panel,r)['values'] for r in rs]).mean(0,dtype=np.float64).astype(np.float32)
  with np.load(path) as f:assert np.array_equal(f['mean'],expected),path
  means+=1
 assert means==11
 scalers=0
 for path in (OUT/'models').glob('PCA_*.npz'):
  name=path.stem[len('PCA_'):];cohort,encoder=name.split('_',1);x=np.concatenate([image_data(encoder,s)['x'] for s in training(cohort)]);sc=core.Scaler().fit(x)
  with np.load(path) as f:assert np.array_equal(sc.mean,f['mean']) and np.array_equal(sc.scale,f['scale']);assert f['components'].shape==(8,x.shape[1]);assert (f['variance']>0).all()
  scalers+=1
 assert scalers==10
 for p in (OUT/'models').glob('B4_*.json'):
  r=json.loads(p.read_text());assert sha(ROOT/r['checkpoint_path'])==r['checkpoint_hash'];assert sha(ROOT/r['scaler_path'])==r['scaler_hash'];assert r['train_slides']==training(r['cohort'])
  assert all(index[m['mask_id']]['sample_id'] in r['train_slides'] and index[m['mask_id']]['checksum']==m['checksum'] for m in r['training_masks']);assert r['no_expression_scaler']
 records=[json.loads(f.read_text()) for f in (OUT/'logs/runs').glob('*.json')];assert len(records)==3456;replays=0;maximum=0.;checks=0
 for r in records:
  mask=index[r['mask_id']];observed=mask_array(mask);hidden=np.flatnonzero(~observed);enc=ENCODERS[0] if r['encoder']=='none' else r['encoder'];d=image_data(enc,r['sample_id']);assert d['feature_cache_hash']==r['feature_cache_hash'];assert mask['checksum']==r['mask_checksum'];assert sha(ROOT/r['prediction_path'])==r['prediction_sha256'];assert r['invariance_errors']==[0,0]
  assert r['config']['interpolation']==selection['simple'][r['panel']][r['regime']]
  assert r['scaler_hash']==sha(OUT/'models'/f'PCA_{r["cohort"]}_{enc}.npz')
  assert r['frozen_image_prediction_hash']==sha(OUT/'cache'/('_'.join([r['cohort'],enc,r['panel'],r['sample_id']])+'_B1.npz'))
  checkpoint=OUT/'SELECTION_LOCK.json'
  if r['model']=='B0':checkpoint=OUT/'models'/f'B0_{r["cohort"]}_{r["panel"]}.npz'
  elif r['model']=='B4':checkpoint=OUT/'models'/f'B4_{r["cohort"]}_{enc}_{r["panel"]}_a{selection["B4_alpha"][r["panel"]][enc]}.npz'
  assert r['checkpoint_hash']==sha(checkpoint)
  with np.load(ROOT/r['prediction_path']) as f:
   assert np.array_equal(f['observed_barcodes'],d['barcodes'][observed]) and np.array_equal(f['hidden_barcodes'],d['barcodes'][hidden]);assert list(f['genes'])==PANELS[r['panel']];assert str(f['mask_checksum'])==r['mask_checksum'];assert len(f['distance_um'])==len(hidden);assert sum(d['n_hidden'] for d in r['distance_bins'])==len(hidden);assert np.isfinite(f['distance_um']).all() and (f['distance_um']>=0).all()
   if r['control'] is not None or (r['nominal_coverage']==.2 and r['replicate']==0):
    # Every model/encoder, all regimes, every evaluated patient/panel, plus all 54 controls.
    s=scene(r['cohort'],enc,r['panel'],r['sample_id']);a=anchors(r['sample_id'],r['panel'],mask);s,a=control_inputs(s,a,r['control'],r['seed']);spec=selection['simple'][r['panel']][r['regime']];lam=spec['lambda'].get(r['encoder']);w=None
    if r['model']=='B4':
     alpha=selection['B4_alpha'][r['panel']][enc]
     with np.load(OUT/'models'/f'B4_{r["cohort"]}_{enc}_{r["panel"]}_a{alpha}.npz') as z:w=z['weights']
    with np.load(OUT/'models'/f'B0_{r["cohort"]}_{r["panel"]}.npz') as z:mean=z['mean']
    p=predict(s,a,hidden,r['model'],spec.get(r['model']),spec,lam,w,mean);err=float(np.max(abs(p-f['prediction'])));assert np.array_equal(p,f['prediction']),r['run_id'];maximum=max(maximum,err);y=labels(r['sample_id'],r['panel'])[hidden];m=core.gene_metrics(y,p)
    for k,v in m.items():
     if k in f.files:assert np.array_equal(v,f[k],equal_nan=True),r['run_id']
    replays+=1
  checks+=1
  if checks%500==0:print('AUDITED',checks,'replays',replays,flush=True)
 events=[json.loads(f.read_text()) for f in (OUT/'logs/registry').glob('*.json')];rr=pd.read_csv(reg);rr=rr[rr.run_id.str.startswith('C5_')];assert len(rr)==len(events) and not rr.duplicated(['event','run_id']).any();assert replays==351,replays
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C4/freeze.py'),'--verify'],check=True)
 result={'status':'PASS','prediction_files_hashed_and_identity_checked':checks,'hidden_label_corruption_prediction_comparisons':checks*2,'corruption_maximum_error':0.,'stratified_prediction_and_metric_replays':replays,'replay_maximum_error':maximum,'mask_files_verified':len(index),'train_only_image_scalers_verified':scalers,'training_anchor_means_verified':means,'expression_scalers_fitted':0,'B4_checkpoints_verified':len(list((OUT/'models').glob('B4_*.json'))),'registry_events':len(events),'C0_C4_unchanged':True,'seconds':time.perf_counter()-start};dump(OUT/'logs/INTEGRITY_AUDIT.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):main()
