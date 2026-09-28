"""C2 integrity and model replay audit; never runs C1 prediction or fitting."""
import json,subprocess,time
import numpy as np
import pandas as pd
from scipy import sparse
import torch
from threadpoolctl import threadpool_limits
from run import ROOT,OUT,CFG,C1LOCK,SPLITS,IDS,PANELS,core,sha,load,features,scaler,model_load,dump,INPUT_CACHE
from graph import build_graph

def main():
 start=time.perf_counter()
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C1/freeze.py'),'--verify'],check=True)
 selection=json.loads((OUT/'SELECTION_LOCK.json').read_text());digest=sha(OUT/'SELECTION_LOCK.json')
 for path,h in selection['source_sha256'].items():assert sha(ROOT/path)==h,path
 assert sha(OUT/'configs/protocol.json')==selection['protocol_sha256']
 assert sha(ROOT/'C1/C1_SPEC_LOCK.json')==selection['c1_lock_sha256']
 analysis=json.loads((OUT/'ANALYSIS_LOCK.json').read_text());access=json.loads((OUT/'logs/TEST_ACCESS.json').read_text())
 assert analysis['timestamp']<access['timestamp'] and selection['timestamp']<access['timestamp']
 assert sha(OUT/'configs/analysis.json')==analysis['analysis_sha256']
 assert sha(OUT/'C2_VALIDATION_GENE_SHORTLIST.csv')==selection['validation_shortlist_sha256']
 runs=[json.loads(p.read_text()) for p in (OUT/'logs/runs').glob('*.json')];evs=[json.loads(p.read_text()) for p in (OUT/'logs/evaluation').glob('*.json')]
 assert len(runs)==237,len(runs);assert len(evs)==189,len(evs)
 data={s:load(s) for s in IDS};by_id={r['run_id']:r for r in runs};assert len(by_id)==len(runs)
 graph_count=0
 for path in (OUT/'graphs').glob('*.npz'):
  record=json.loads(path.with_suffix('.json').read_text());assert sha(path)==record['sha256']
  sample=record['sample_id'];a=sparse.load_npz(path);b=build_graph(data[sample]['xy'],record['config'])
  assert (a-b).nnz==0;assert not (a-a.T).nnz and not a.diagonal().any()
  assert len(data[sample]['x'])==a.shape[0]==a.shape[1];graph_count+=1
 expected_graphs={(sample,g['name']) for sample in SPLITS['splits']['train']+SPLITS['splits']['validation'] for g in CFG['graphs']}
 expected_graphs |= {(sample,c['graph']['name']) for sample in IDS for models in selection['selected'].values() for c in models.values()}
 actual_graphs={(json.loads(p.read_text())['sample_id'],json.loads(p.read_text())['config']['name']) for p in (OUT/'graphs').glob('*.json')}
 assert actual_graphs==expected_graphs,(actual_graphs,expected_graphs)
 for r in runs:
  assert sha(ROOT/r['model_path'])==r['model_sha256']
  folder=ROOT/r['preprocessing'];assert sha(folder/'feature_scaler.npz')==r['feature_scaler_sha256'];assert sha(folder/'target_scaler.npz')==r['target_scaler_sha256']
  assert abs(r['parameters']-r['local_parameters'])/r['local_parameters']<.01
  if r['phase']=='lopo':
   assert r['heldout'] not in r['train_slides'] and not r['validation_slides']
   assert set(r['train_slides'])|{r['heldout']}==set(IDS)
   assert r['best_epoch']==selection['selected'][r['panel']][r['model']]['lopo_epochs']
  else:
   assert r['train_slides']==SPLITS['splits']['train'];assert r['validation_slides']==SPLITS['splits']['validation']
 for panel in CFG['panels']:
  for model in CFG['models']:
   candidates=[r for r in runs if r['phase']=='selection' and r['panel']==panel and r['model']==model];assert len(candidates)==4
   best=min(candidates,key=lambda r:(-r['validation_median_pearson'],r['validation_mse'],r['graph']['name']))
   assert best['graph']==selection['selected'][panel][model]['graph']
 max_error=0.
 for ev in evs:
  r=by_id[ev['run_id']];s=ev['sample_id'];assert s not in r['train_slides'];assert ev['selection_lock_sha256']==digest
  assert sha(ROOT/ev['prediction_path'])==ev['prediction_sha256']
  with np.load(ROOT/ev['prediction_path']) as f:
   stored=f['prediction'];assert list(f['genes'])==PANELS[r['panel']];assert np.array_equal(f['barcodes'],data[s]['barcodes'])
  assert np.isfinite(stored).all()
  x=features(s,data[s],r['panel'],r['model'],r['graph'],r['heldout'],r['seed'] if r['shuffled'] else None)
  replay=scaler(ROOT/r['preprocessing'],'target').inverse(core.predict(model_load(r),x,'B2')).astype(np.float32)
  error=float(np.max(abs(replay-stored)));max_error=max(max_error,error);np.testing.assert_allclose(replay,stored,rtol=1e-6,atol=1e-6)
  INPUT_CACHE.clear()
 registry=pd.read_csv(ROOT/'results/experiment_registry.csv');registry=registry[registry.run_id.str.startswith('C2_')]
 assert len(registry)==426 and not registry.duplicated(['event','run_id']).any()
 assert set(registry[registry.event=='fit'].run_id)==set(by_id)
 assert set(registry[registry.event=='evaluation'].run_id)=={e['run_id'] for e in evs}
 report={'status':'PASS','completed_fits':237,'evaluations_replayed':189,'max_replay_error':max_error,'graphs_reconstructed':graph_count,
         'registry_events':426,'C0_C1_integrity':'unchanged, verified by frozen hashes only','C1_models_run':0,
         'validation_only_selection_verified':True,'capacity_difference_below_one_percent':True,'seconds':time.perf_counter()-start}
 if not (OUT/'C2_SPEC_LOCK.json').exists():dump(OUT/'logs/INTEGRITY_AUDIT.json',report)
 print(json.dumps(report,indent=2))
if __name__=='__main__':
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):main()
