"""Audit C1 inputs, split/scaler provenance, model replay and output completeness."""
import json,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits
from run import ROOT,OUT,CFG,MANIFEST,SPLITS,IDS,GENES,PANELS,sha,load,load_scaler,load_model,combine,panel_indices,dump
from core import predict,Scaler

def main():
 start=time.perf_counter()
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'scripts/verify_c0.py'),'--include-data'],check=True)
 lock=json.loads((OUT/'SELECTION_LOCK.json').read_text());lock_sha=sha(OUT/'SELECTION_LOCK.json')
 assert sha(ROOT/'reports/C0/C0_SPEC_LOCK.json')==lock['c0_lock_sha256']
 assert sha(OUT/'configs/protocol.json')==lock['protocol_sha256']
 for path,digest in {**lock['source_sha256'],**lock['scaler_sha256']}.items():assert sha(ROOT/path)==digest,path
 runs=[json.loads(p.read_text()) for p in (OUT/'logs/runs').glob('*.json')]
 evaluations=[json.loads(p.read_text()) for p in (OUT/'logs/evaluation').glob('*.json')]
 assert len(runs)==156, len(runs)
 assert len(evaluations)==112,len(evaluations)
 data={s:load(s) for s in IDS};by_id={r['run_id']:r for r in runs}
 seen=set();checked_scalers=set();replayed=0;max_error=0.
 for r in runs:
  assert r['run_id'] not in seen;seen.add(r['run_id'])
  assert sha(ROOT/r['model_path'])==r['model_sha256']
  if r['phase']=='lopo':
   assert len(r['evaluation_slides'])==1
   assert set(r['train_slides']).isdisjoint(r['evaluation_slides'])
   assert set(r['train_slides'])|set(r['evaluation_slides'])==set(IDS)
   if r['model']=='B2':assert r['best_epoch']==lock['selected'][r['panel']]['lopo_fixed_epochs']
  else:
   assert r['train_slides']==SPLITS['splits']['train']
   assert r['evaluation_slides']==SPLITS['splits']['validation']
  if r['preprocessing'] not in checked_scalers:
   folder=ROOT/r['preprocessing'];ix=panel_indices(r['panel'])
   for key,filename,indices in [('x','feature_scaler.npz',None),('y','target_scaler.npz',ix)]:
    values=combine(data,r['train_slides'],key)
    if indices is not None:values=values[:,indices]
    expected=Scaler().fit(values);actual=load_scaler(folder/filename)
    np.testing.assert_allclose(expected.mean,actual.mean,rtol=0,atol=1e-12)
    np.testing.assert_allclose(expected.scale,actual.scale,rtol=0,atol=1e-12)
   checked_scalers.add(r['preprocessing'])
 for ev in evaluations:
  r=by_id[ev['run_id']];sample=ev['sample_id']
  assert ev['selection_lock_sha256']==lock_sha
  assert sample not in r['train_slides']
  assert sha(ROOT/ev['prediction_path'])==ev['prediction_sha256']
  with np.load(ROOT/ev['prediction_path']) as f:
   stored=f['prediction'];assert np.array_equal(f['barcodes'],data[sample]['barcodes']);assert list(f['genes'])==PANELS[r['panel']]
  assert stored.shape==(len(data[sample]['x']),len(PANELS[r['panel']])) and np.isfinite(stored).all()
  folder=ROOT/r['preprocessing'];xs=load_scaler(folder/'feature_scaler.npz');ys=load_scaler(folder/'target_scaler.npz')
  replay=ys.inverse(predict(load_model(r),xs.transform(data[sample]['x']),r['model'])).astype(np.float32)
  error=float(np.max(abs(stored-replay)));max_error=max(max_error,error)
  np.testing.assert_allclose(replay,stored,rtol=1e-6,atol=1e-6);replayed+=1
 registry=pd.read_csv(ROOT/'results/experiment_registry.csv')
 registry=registry[registry.run_id.isin(seen)]
 assert len(registry)==268,len(registry)
 assert not registry.duplicated(['event','run_id']).any()
 assert set(registry[registry.event=='fit'].run_id)==seen
 assert set(registry[registry.event=='evaluation'].run_id)=={e['run_id'] for e in evaluations}
 access=json.loads((OUT/'logs/TEST_ACCESS.json').read_text())
 assert access['timestamp']>lock['timestamp']
 report={'status':'PASS','fit_runs':len(runs),'evaluated_runs':len(evaluations),'registry_events':len(registry),
         'scaler_pairs_verified':len(checked_scalers),'model_predictions_replayed':replayed,
         'maximum_replay_absolute_error':max_error,'c0_integrity':'261 frozen hashes match',
         'primary_test_access_after_selection':True,'selection_lock_sha256':lock_sha,'seconds':time.perf_counter()-start}
 if not (OUT/'C1_SPEC_LOCK.json').exists():
  dump(OUT/'logs/INTEGRITY_AUDIT.json',report)
 else:
  previous=json.loads((OUT/'logs/INTEGRITY_AUDIT.json').read_text())
  for key in report:
   if key!='seconds':assert report[key]==previous[key],key
 print(json.dumps(report,indent=2))
if __name__=='__main__':
 torch.set_num_threads(CFG['threads']);torch.use_deterministic_algorithms(True)
 with threadpool_limits(CFG['threads']):main()
