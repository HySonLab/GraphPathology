"""Independent C8 replay, fold exclusion, frozen-baseline and artifact audit."""
from contract import *
import run as pipeline
from threadpoolctl import threadpool_limits

def main():
 start=time.perf_counter();lock=check_lock();fitted=json.loads((A/'FIT_LOCK.json').read_text());access=json.loads((A/'logs/OUTER_ACCESS.json').read_text());assert lock['timestamp']<fitted['timestamp']<access['timestamp'];assert access['fit_lock_sha256']==sha(A/'FIT_LOCK.json')
 for k in ['input_sha256','mask_sha256']:
  for p,h in lock[k].items():assert sha(ROOT/p)==h,p
 prefix=lock['registry_prefix']
 with (ROOT/prefix['path']).open('rb') as f:assert hashlib.sha256(f.read(prefix['bytes'])).hexdigest()==prefix['sha256'];tail=f.read().decode()
 scaler_checks=0;fit_selections=0;gate_replays=0
 for fold,held in cfg()['folds'].items():
  rec=fitted['folds'][fold];assert set(rec['train_patients'])==set(SAMPLES)-set(held);set_forbidden(held)
  for p,h in rec['scaler_files'].items():
   assert sha(ROOT/p)==h;old=json.loads((ROOT/p).read_text());assert not set(old['train_patients'])&set(held);assert old==fit_scaler(old['train_patients'],old['transform']);scaler_checks+=1
  df=pd.read_csv(A/'logs'/f'{fold}_INNER_METRICS.csv');assert not set(df.sample_id)&set(held)
  for r in df.to_dict('records'):
   # CSV records retain lists as strings; authoritative scaler JSON identifies fit patients.
   tr=json.loads((ROOT/r['scaler_path']).read_text())['train_patients'];assert r['sample_id'] not in tr and not set(tr)&set(held)
  chosen,_=select(df[df.model=='E2'].to_dict('records'),'transform');assert chosen==rec['selected_transform'];fit_selections+=1
  if (B/'C8B_SPEC_LOCK.json').exists():
   bl=json.loads((B/'C8B_SPEC_LOCK.json').read_text());op=ROOT/bl['fold_objects'][fold]['path'];assert sha(op)==bl['fold_objects'][fold]['sha256'];obj=json.loads(op.read_text());assert set(obj['train_patients'])==set(rec['train_patients']);index=json.loads((A/'logs'/f'{fold}_INNER_CACHE.json').read_text());records=[]
   for rr in index:
    if rr['transform']!=chosen:continue
    assert rr['sample_id'] in rec['train_patients'];assert sha(ROOT/rr['path'])==rr['sha256']
    with np.load(ROOT/rr['path']) as f:v={k:f[k].copy() for k in ['image','st','features','distance','scale','center']};t=f['truth'].copy()
    records.append({**rr,'v':v,'truth':t})
   gate=fit_gate(records);assert gate==obj['family_models']['F3'];assert gate['optimization_success'],gate;gate_replays+=1
   trials=pd.read_csv(B/'logs'/f'{fold}_CALIBRATION_METRICS.csv');win,_=select(trials.to_dict('records'),'candidate');assert win.split('_')[0]==obj['selected_family'];fit_selections+=1
 set_forbidden([])
 c7rows=[json.loads(p.read_text()) for p in (ROOT/'C7/C7B/logs/runs').glob('*.json')];baseline={(r['mask_id'],r['nominal_budget'],r['model']):r for r in c7rows if r['panel']=='S2-50' and not r['control']};count=0;replays=0;corrupt=0;zero=0;records=[]
 for stage in ['C8A','C8B','C8C']:
  if not (OUT/stage/'logs/runs').exists():continue
  for path in sorted((OUT/stage/'logs/runs').glob('*.json')):
   r=json.loads(path.read_text());records.append(r);pp=ROOT/r['prediction_path'];assert sha(pp)==r['prediction_sha256'];assert r['hidden_corruption_checks']==2;assert r['sample_id'] in cfg()['folds'][r['fold']];obj=json.loads((A/'scalers'/f'{r["fold"]}_{r["transform"]}.json').read_text());assert r['sample_id'] not in obj['train_patients'];assert sha(ROOT/r['mask_path'])==r['mask_sha256']
   with np.load(pp) as f,np.load(ROOT/r['mask_path']) as mask:
    obs,h,order=arrays(r);assert np.array_equal(h,f['hidden_indices']) and np.array_equal(order,f['observed_indices']);assert np.array_equal(mask['barcodes'][h],f['hidden_barcodes']) and np.array_equal(mask['barcodes'][order],f['observed_barcodes']);assert list(f['genes'])==GENES;assert np.array_equal(inverse(f['prediction_native'],obj),f['prediction_common']);dist=np.round(cKDTree(mask['xy_um'][order]).query(mask['xy_um'][h])[0],6);assert np.array_equal(dist,f['distance_um']);count+=1
    if stage=='C8A' and r['transform']=='T0':
     with np.load(ROOT/baseline[(r['mask_id'],r['budget'],r['model'])]['prediction_path']) as old:assert np.array_equal(old['prediction'],f['prediction_common']);zero+=1
    if r['seed']!=0 and not r['control']:continue
    with allowed([r['sample_id']]):y=target(r['sample_id'])
    v=inputs(r,obj,y)
    if stage=='C8A':pred=predict_A(v,obj,r['model'],r['regime']);model=None
    else:
     cal=json.loads((ROOT/r['calibration_object_path']).read_text());model=cal['selected_model'] if r['model']=='SELECTED' else cal['family_models'][r['model']];vv,mo=pipeline.control_inputs(v,model,r['control'],r['control_seed']) if r['control'] else (v,model);pred=pipeline.predict_B(vv,mo)
    assert np.array_equal(pred,f['prediction_native']);met=core.gene_metrics(convert(y[h],obj),pred)
    for k,val in met.items():assert np.array_equal(val,f[k],equal_nan=True)
    for bad in [np.nan,808.]:
     alt=y.copy();alt[~obs]=bad;vv=inputs(r,obj,alt)
     if stage=='C8A':pred2=predict_A(vv,obj,r['model'],r['regime'])
     else:
      vv,mo=pipeline.control_inputs(vv,model,r['control'],r['control_seed']) if r['control'] else (vv,model);pred2=pipeline.predict_B(vv,mo)
     assert np.array_equal(pred,pred2);corrupt+=1
    replays+=1
 events=[json.loads(p.read_text()) for p in (OUT/'logs/registry').glob('*.json')];assert len(events)==len({r['run_id'] for r in events});assert len(tail.splitlines())==len(events)
 for r in events:assert r['run_id'] in tail
 assert scaler_checks==84 and zero==1440 and count==len(records)
 result={'status':'PASS','prediction_files_verified':count,'exact_C7_T0_prediction_replays':zero,'independent_prediction_metric_replays':replays,'independent_corruption_checks':corrupt,'recorded_outer_corruption_checks':count*2,'fold_scalers_rebuilt_with_heldout_access_forbidden':scaler_checks,'inner_selections_replayed':fit_selections,'gates_refit_exactly':gate_replays,'registry_events':len(events),'original_registry_prefix_unchanged':True,'scaler_and_prediction_identity_and_distances_verified':True,'fit_objects_sealed_before_outer_evaluation':True,'seconds':time.perf_counter()-start};dump(OUT/'logs/INTEGRITY_AUDIT.json',result);print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 with threadpool_limits(4):main()
