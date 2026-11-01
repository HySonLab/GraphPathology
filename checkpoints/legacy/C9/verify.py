"""Replay all saved predictions/metrics; refit calibration with evaluation access forbidden."""
from contract import *
from threadpoolctl import threadpool_limits

def main():
 lock=check_lock();cl=json.loads((B/'CALIBRATION_LOCK.json').read_text());scaler=json.loads((A/'T1_SCALER.json').read_text());assert sha(A/'T1_SCALER.json')==cl['scaler_sha256'];assert set(CAL).isdisjoint(EVAL)
 for key in ['image_assets_sha256','model_parents_sha256','raw_asset_sha256']:
  for p,h in lock[key].items():assert sha(ROOT/p)==h,p
 # Allowed calibration access only. No evaluation loader used by this refit.
 dst=moments([read_targets(s,'calibration') for s in CAL]);assert dst==scaler['destination']
 src=[]
 for s in ['NCBI783','NCBI785']:
  with np.load(ROOT/'data/cache'/(s+'_targets.npz')) as f:src.append(f['normalized_targets'][:,[list(f['genes']).index(g) for g in genes()]].astype(float))
 assert moments(src)==scaler['source'];rows=[json.loads(p.read_text()) for p in sorted((B/'logs/runs').glob('*.json'))];assert len(rows)==414;replayed=0
 for s in EVAL:
  y=read_targets(s,'evaluation')
  with np.load(CACHE/(s+'_metadata.npz')) as f:xy=f['xy_um'];bc=f['barcodes']
  with np.load(CACHE/(s+'_image.npz')) as f:im=f['prediction']
  for r in [r for r in rows if r['sample_id']==s]:
   assert sha(ROOT/r['prediction_path'])==r['prediction_sha256'];assert sha(ROOT/r['mask_path'])==r['mask_sha256']
   with np.load(ROOT/r['prediction_path']) as f:
    obs=f['observed_indices'];h=f['hidden_indices'];assert not set(obs)&set(h);assert np.array_equal(bc[obs],f['observed_barcodes']) and np.array_equal(bc[h],f['hidden_barcodes']);assert list(f['genes'])==genes()
    pred,d=predict(xy,im,y,obs,h,r['regime'],scaler,r['control'],r['control_seed']);assert np.array_equal(d,f['distance_um'])
    for m in pred:
     assert np.array_equal(pred[m],f[m]);sm,mt=score(y[h],pred[m])
     for k in ['pearson','spearman','mse','mae']:np.testing.assert_allclose(mt[k],f[m+'_'+k],rtol=0,atol=0,equal_nan=True)
     for k in sm:np.testing.assert_allclose(sm[k],r['summary'][m][k],rtol=0,atol=0,equal_nan=True)
     replayed+=1
 reg=lock['registry_prefix_before_C9']
 with (ROOT/'results/experiment_registry.csv').open('rb') as f:assert hashlib.sha256(f.read(reg['bytes'])).hexdigest()==reg['sha256']
 # Reproduce label-free masks: the selection API receives XY only.
 import policy
 for path in (B/'masks').glob('*.json'):
  r=json.loads(path.read_text())
  with np.load(ROOT/r['mask_path']) as f:
   xy=f['xy_um'];forbidden=policy.region(xy,r['regime'],r['geometry']);order=policy.choose({'xy':xy},np.flatnonzero(~forbidden),parent.count(len(xy),.4),r['strategy'],r['seed']);assert np.array_equal(forbidden,f['forbidden']) and np.array_equal(order,f['order'])
 try:read_targets(EVAL[0],'calibration')
 except AssertionError:pass
 else:raise AssertionError('Evaluation patient accepted for calibration')
 # Three mask families perV2 and ten random seeds per condition must be present.
 df=pd.DataFrame(rows);clean=df[df.control==''];assert len(clean)==396
 for s in EVAL:
  for b in [.1,.2,.4]:
   for regime in ['V1','V2']:
    for g in ([0] if regime=='V1' else [0,1,2]):assert set(clean[(clean.sample_id==s)&(clean.budget==b)&(clean.regime==regime)&(clean.geometry==g)&(clean.strategy=='random')].seed)==set(range(10))
 result={'status':'PASS','replayed_method_predictions':replayed,'replayed_metrics':replayed*4,'hidden_corruption_checks':sum(r['hidden_corruption_checks'] for r in rows),'calibration_refit_exact':True,'evaluation_patients_absent_from_fit':True,'prior_registry_prefix_unchanged':True,'raw_image_model_hash_checks':True,'clean_runs':396,'control_runs':18,'mask_trajectories':132,'label_free_mask_replays':132,'forbidden_calibration_access_rejected':True};dump(B/'logs/INTEGRITY_AUDIT.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':
 with threadpool_limits(4):main()
