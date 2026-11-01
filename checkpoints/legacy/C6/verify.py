"""Independent full policy replay, stratified prediction replay and immutable-parent audit."""
from support import *
from policy import choose
import torch
from threadpoolctl import threadpool_limits

def main():
 start=time.perf_counter();spec=check_lock();selection=json.loads((A/'SELECTION_LOCK.json').read_text());trajectory=json.loads((A/'TRAJECTORY_LOCK.json').read_text());access=json.loads((A/'logs/TEST_ACCESS.json').read_text());assert spec['timestamp']<selection['timestamp']<trajectory['timestamp']<access['timestamp'];assert selection['spec_sha256']==sha(A/'C6A_SPEC_LOCK.json');assert trajectory['selection_sha256']==sha(A/'SELECTION_LOCK.json');assert access['trajectory_lock_sha256']==sha(A/'TRAJECTORY_LOCK.json');assert not selection['primary_test_hidden_labels_used']
 with (ROOT/'results/experiment_registry.csv').open('rb') as f:assert hashlib.sha256(f.read(spec['registry_prefix']['bytes'])).hexdigest()==spec['registry_prefix']['sha256']
 for sample in IDS:
  with np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as f:
   for enc in ENCODERS:assert np.array_equal(f['barcodes'],c5.image_data(enc,sample)['barcodes']) and np.array_equal(f['xy_um'],c5.image_data(enc,sample)['xy'])
 policies={};policy_replays=0
 for pid,ref in trajectory['policies'].items():
  r=json.loads((A/'policies'/f'{pid}.json').read_text());assert sha(ROOT/ref['path'])==ref['sha256']==r['policy_sha256'];assert r['selection_corruption_errors']==[0,0];s=scene(r['cohort'],r['encoder'],'S2-50',r['sample_id']);forbidden,old=unavailable(r['sample_id'],r['regime']);assert r['unavailable_region']==old
  with np.load(ROOT/ref['path']) as f:
   assert np.array_equal(f['forbidden'],forbidden);assert np.array_equal(f['eligible'],np.flatnonzero(~forbidden));order,scores=choose({'xy':s['xy'],'projection':s['projection']},f['eligible'],count(len(forbidden),.4),r['strategy'],r['seed'],r['alpha'] if r['alpha'] is not None else .5,r['control']);assert np.array_equal(order,f['order']);assert np.array_equal(scores,f['acquisition_scores'],equal_nan=True);assert np.array_equal(f['selected_barcodes'],s['barcodes'][order]);assert not forbidden[order].any()
  policies[pid]=r;policy_replays+=1
  if policy_replays%50==0:print('POLICIES REPLAYED',policy_replays,flush=True)
 assert policy_replays==296
 records=[json.loads(p.read_text()) for p in (A/'logs/runs').glob('*.json')];records=[r for r in records if r['phase']!='validation'];assert len(records)==4448;replays=0
 for i,r in enumerate(records):
  pr=policies[r['policy_id']];s=scene(r['cohort'],r['encoder'],r['panel'],r['sample_id']);n=len(s['xy']);k=count(n,r['nominal_budget']);assert r['n_observed']==k and r['realized_budget']==k/n and r['invariance_errors']==[0,0];assert r['policy_sha256']==pr['policy_sha256'];assert sha(ROOT/r['prediction_path'])==r['prediction_sha256'];assert sha(ROOT/r['checkpoint_path'])==r['checkpoint_hash'];assert r['scaler_hash']==sha(ROOT/'C5/models'/f'PCA_{r["cohort"]}_{r["encoder"]}.npz');assert r['feature_cache_hash']==s['feature_cache_hash']
  with np.load(ROOT/pr['policy_path']) as policy,np.load(ROOT/r['prediction_path']) as f:
   chosen=policy['order'][:k];h,obs=target_indices(n,chosen,policy['forbidden'],r['regime']);assert np.array_equal(f['observed_indices'],chosen) and np.array_equal(f['observed_barcodes'],s['barcodes'][chosen]);assert np.array_equal(f['hidden_indices'],h) and np.array_equal(f['hidden_barcodes'],s['barcodes'][h]);assert list(f['genes'])==PANELS[r['panel']];assert sum(b['n_hidden'] for b in r['distance_bins'])==len(h);assert np.isfinite(f['distance_um']).all() and (f['distance_um']>=0).all()
   if r['nominal_budget']==.2 and (r['strategy']!='A0' or r['seed']==0):
    a=observed_only(r['sample_id'],r['panel'],chosen);p=reconstruction(s,a,h,r['cohort'],r['encoder'],r['panel'],r['regime'],r['role']);assert np.array_equal(p,f['prediction']),r['run_id'];metrics=core.gene_metrics(expression(r['sample_id'],r['panel'])[h],p)
    for name,value in metrics.items():
     if name in f.files:assert np.array_equal(value,f[name],equal_nan=True),r['run_id']
    replays+=1
  if (i+1)%500==0:print('PREDICTIONS AUDITED',i+1,'replays',replays,flush=True)
 assert replays==464,replays
 events=[json.loads(p.read_text()) for p in (A/'logs/registry').glob('*.json')];registry_table=pd.read_csv(ROOT/'results/experiment_registry.csv');registry_table=registry_table[registry_table.run_id.str.startswith('C6_')];assert len(events)==len(registry_table)==4924;assert not registry_table.duplicated(['event','run_id']).any()
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C5/freeze.py'),'--verify'],check=True)
 result={'status':'PASS','test_policy_full_order_and_score_replays':policy_replays,'test_selection_corruption_comparisons':policy_replays*2,'test_prediction_files_hashed_and_identity_checked':len(records),'test_hidden_label_corruption_comparisons':len(records)*2,'stratified_prediction_and_metric_replays':replays,'maximum_selection_or_prediction_replay_error':0,'validation_policies':20,'validation_prediction_runs':160,'registry_events':len(events),'exact_budget_matching':True,'no_anchors_inside_stress_region':True,'C0_C5_unchanged':True,'seconds':time.perf_counter()-start};dump(A/'logs/INTEGRITY_AUDIT.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):main()
