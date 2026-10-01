"""Independent artifact, mask, target-transform and frozen-inference replay audit."""
from support import *
from policy import choose,region
from threadpoolctl import threadpool_limits

def verify():
 start=time.perf_counter();lock=check_lock();ml=json.loads((B/'MASK_LOCK.json').read_text());access=json.loads((B/'logs/EXTERNAL_LABEL_ACCESS.json').read_text());assert lock['timestamp']<ml['timestamp']<access['timestamp'];assert ml['spec_hash']==sha(A/'C7A_SPEC_LOCK.json') and access['mask_lock_sha256']==sha(B/'MASK_LOCK.json')
 for collection in ['qualified_assets_sha256','frozen_model_and_preprocessing_sha256']:
  for p,h in lock[collection].items():assert sha(ROOT/p)==h,p
 prefix=lock['registry_prefix']
 with (ROOT/prefix['path']).open('rb') as f:assert hashlib.sha256(f.read(prefix['bytes'])).hexdigest()==prefix['sha256']
 maskrows=[]
 for p,h in ml['mask_files'].items():
  assert sha(ROOT/p)==h;r=json.loads((ROOT/p).with_suffix('.json').read_text());maskrows.append(r)
  with np.load(ROOT/p) as f:xy=f['xy_um'];bc=f['barcodes'];forbidden=region(xy,r['regime'],r['geometry']);assert np.array_equal(forbidden,f['forbidden']);eligible=np.flatnonzero(~forbidden);k=count(len(bc),.4)
  for val in [np.full((len(bc),280),np.nan),np.random.default_rng(999).normal(size=(len(bc),280))]:
   order=choose({'xy':xy,'expression':val},eligible,k,r['strategy'],r['seed'])
   with np.load(ROOT/p) as f:assert np.array_equal(order,f['order']);assert np.array_equal(bc[order],f['ordered_barcodes'])
 assert len(maskrows)==352
 norm=json.loads((ROOT/'reports/C0/C0_NORMALIZATION.json').read_text());target_checks=0
 for sample in SAMPLES:
  with np.load(CACHE/(sample+'_targets.npz')) as t,np.load(CACHE/(sample+'_metadata.npz')) as m:
   assert np.array_equal(t['barcodes'],m['barcodes']);idx=[norm['genes'].index(g) for g in t['genes']];mean=np.array(norm['mean'])[idx];scale=np.array(norm['scale'])[idx];expected=((np.log1p(t['raw_counts'])-mean)/scale).astype('float32');assert np.array_equal(expected,t['normalized_targets']);assert np.array_equal(t['raw_counts'],np.rint(t['raw_counts']));obs=np.arange(len(expected))%5==0
   for val in [np.nan,9999.]:
    raw=t['raw_counts'].astype(float);raw[~obs]=val;observed=((np.log1p(raw[obs])-mean)/scale).astype('float32');assert np.array_equal(observed,t['normalized_targets'][obs]);target_checks+=1
 rows=[json.loads(p.read_text()) for p in (B/'logs/runs').glob('*.json')];assert len(rows)==8556;assert len({r['run_id'] for r in rows})==len(rows);replays=0;corrupt_replays=0;max_error=0
 for r in rows:
  path=ROOT/r['prediction_path'];assert sha(path)==r['prediction_sha256'];assert r['hidden_corruption_checks']==2
  with np.load(ROOT/r['mask_path']) as m,np.load(path) as f:
   bc=m['barcodes'];order=m['order'][:count(len(bc),r['nominal_budget'])];obs=np.zeros(len(bc),bool);obs[order]=True;h=np.flatnonzero(~obs if r['regime']=='V1' else m['forbidden']);assert not m['forbidden'][obs].any();assert np.array_equal(order,f['observed_indices']) and np.array_equal(bc[order],f['observed_barcodes']);assert np.array_equal(h,f['hidden_indices']) and np.array_equal(bc[h],f['hidden_barcodes']);assert list(f['genes'])==panels()[r['panel']];assert r['n_hidden']==len(h) and r['n_observed']==len(order);saved=f['prediction'].copy()
   replay=r['control'] is not None or (r['nominal_budget']==.2 and r['seed']==0)
   if not replay:continue
   s=scene(r['sample_id'],r['panel']);y=truth(r['sample_id'],r['panel']);a=c5.observed_payload(y,obs,s['xy']);ss=s
   if r['control']:ss,a=c5.control_inputs(s,a,r['control'],r['control_seed'])
   pred=predict(ss,a,h,r['panel'],r['regime'],r['model']);assert np.array_equal(saved,pred);met=core.gene_metrics(y[h],pred)
   for k,v in met.items():assert np.array_equal(v,f[k],equal_nan=True)
   for k,v in core.summarize(met).items():assert (v==r[k]) or (not np.isfinite(v) and r[k] is None)
   for val in [np.nan,12345.]:
    alt=y.copy();alt[~obs]=val;aa=c5.observed_payload(alt,obs,s['xy']);sss=s
    if r['control']:sss,aa=c5.control_inputs(s,aa,r['control'],r['control_seed'])
    assert np.array_equal(saved,predict(sss,aa,h,r['panel'],r['regime'],r['model']));corrupt_replays+=1
   replays+=1
 assert replays==620 and corrupt_replays==1240
 events=[json.loads(p.read_text()) for p in (B/'logs/registry').glob('*.json')];assert len(events)==8908 and len({r['run_id'] for r in events})==8908
 with (ROOT/'results/experiment_registry.csv').open('rb') as f:f.read(prefix['bytes']);tail=f.read().decode();assert len(tail.splitlines())==8908
 for r in events:assert r['run_id'] in tail
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C6/freeze.py'),'--verify'],check=True)
 result={'status':'PASS','mask_trajectories_replayed':352,'mask_corruption_comparisons':704,'prediction_files_hashed_and_identity_checked':8556,'recorded_prediction_corruption_comparisons':17112,'independent_prediction_and_metric_replays':replays,'additional_replay_corruption_comparisons':corrupt_replays,'raw_observed_normalization_corruption_comparisons':target_checks,'maximum_replay_error':max_error,'registry_events':8908,'C0_C6_unchanged':True,'source_archives_and_frozen_models_verified':True,'selection_preceded_external_labels':True,'seconds':time.perf_counter()-start};dump(B/'logs/INTEGRITY_AUDIT.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):verify()
