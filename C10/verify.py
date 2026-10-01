"""Exact prediction/metric replay and frozen-parent/input integrity audit."""
from contract10 import *
from threadpoolctl import threadpool_limits

def main():
 lock=check_lock();obj=scaler();assert cfg()['calibration_patients']==[] and not set(SAMPLES)&set(obj['calibration_patients'])
 for key in ['image_assets_sha256','model_parents_sha256','raw_asset_sha256']:
  for p,h in lock[key].items():assert sha(ROOT/p)==h,p
 rows=[json.loads(p.read_text()) for p in sorted((B/'logs/runs').glob('*.json'))];assert len(rows)==len(SAMPLES)*132+18;replayed=0
 for s in SAMPLES:
  y=targets(s)
  with np.load(CACHE/(s+'_metadata.npz')) as f:xy=f['xy_um'];bc=f['barcodes']
  with np.load(CACHE/(s+'_image.npz')) as f:im=f['prediction'];assert np.array_equal(bc,f['barcodes'])
  for r in [r for r in rows if r['patient_id']==s]:
   assert sha(ROOT/r['prediction_path'])==r['prediction_sha256'];assert sha(ROOT/r['mask_path'])==r['mask_sha256'];assert r['calibration_patients']==[]
   with np.load(ROOT/r['prediction_path']) as f:
    obs=f['observed_indices'];h=f['hidden_indices'];assert not set(obs)&set(h);assert np.array_equal(bc[obs],f['observed_barcodes']) and np.array_equal(bc[h],f['hidden_barcodes']);assert list(f['genes'])==genes();pred,d=predict(xy,im,y,obs,h,r['regime'],obj,r['control'],r['control_seed']);assert np.array_equal(d,f['distance_um'])
    for m in pred:
     assert np.array_equal(pred[m],f[m]);sm,met=score(y[h],pred[m])
     for k in ['pearson','spearman','mse','mae']:np.testing.assert_allclose(met[k],f[m+'_'+k],rtol=0,atol=0,equal_nan=True)
     for k in sm:np.testing.assert_allclose(sm[k],r['summary'][m][k],rtol=0,atol=0,equal_nan=True)
     replayed+=1
  print('REPLAYED',s,flush=True)
 for path in (B/'masks').glob('*.json'):
  r=json.loads(path.read_text())
  with np.load(ROOT/r['mask_path']) as f:
   xy=f['xy_um'];forbidden=frozen_policy.region(xy,r['regime'],r['geometry']);order=frozen_policy.choose({'xy':xy},np.flatnonzero(~forbidden),parent.count(len(xy),.4),r['strategy'],r['seed']);assert np.array_equal(forbidden,f['forbidden']) and np.array_equal(order,f['order'])
 reg=lock['registry_prefix_before_C10']
 with (ROOT/'results/experiment_registry.csv').open('rb') as f:assert hashlib.sha256(f.read(reg['bytes'])).hexdigest()==reg['sha256']
 df=pd.DataFrame(rows);clean=df[df.control==''];assert len(clean)==len(SAMPLES)*132
 for s in SAMPLES:
  for b in [.1,.2,.4]:
   for regime in ['V1','V2']:
    for g in ([0] if regime=='V1' else [0,1,2]):assert set(clean[(clean.patient_id==s)&(clean.budget==b)&(clean.regime==regime)&(clean.geometry==g)&(clean.strategy=='random')].seed)==set(range(10))
 result={'status':'PASS','patients':len(SAMPLES),'replayed_method_predictions':replayed,'replayed_metric_arrays':replayed*4,'hidden_corruption_checks':sum(r['hidden_corruption_checks'] for r in rows),'no_C10_calibration_or_fitting':True,'prior_registry_prefix_unchanged':True,'source_image_model_hash_checks':True,'clean_runs':len(clean),'control_runs':18,'label_free_mask_replays':len(SAMPLES)*44};dump(B/'logs/INTEGRITY_AUDIT.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':
 with threadpool_limits(4):main()
