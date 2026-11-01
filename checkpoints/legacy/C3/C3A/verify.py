"""Audit C3A without training or replaying C1/C2 models."""
from common import *
import time
from scipy import sparse
from graph_controls import randomize,audit,canonical_hash
spec=importlib.util.spec_from_file_location('c3_runtime',OUT/'run.py');run=importlib.util.module_from_spec(spec);spec.loader.exec_module(run)
import torch
from threadpoolctl import threadpool_limits

def verify():
 start=time.perf_counter();check_lock()
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C2/freeze.py'),'--verify'],check=True)
 data={s:load(s) for s in IDS};scaler_checks=0
 for heldout in [None]+IDS:
  train=[s for s in IDS if s!=heldout] if heldout else SPLITS['splits']['train']
  for panel in PANELS:
   folder=preprocessing(panel,heldout)
   for name,values in [('feature',np.concatenate([data[s]['x'] for s in train])),('target',np.concatenate([data[s]['y'][:,ix(panel)] for s in train]))]:
    expected=core.Scaler().fit(values);frozen=scaler(folder,name)
    assert np.allclose(expected.mean,frozen.mean,rtol=1e-7,atol=1e-7),(heldout,panel,name,'mean')
    assert np.allclose(expected.scale,frozen.scale,rtol=1e-7,atol=1e-7),(heldout,panel,name,'scale')
    scaler_checks+=1
 graphs=0
 for p in sorted((OUT/'graphs').glob('*.json')):
  r=json.loads(p.read_text());a=sparse.load_npz(ROOT/r['path']);reference=sparse.load_npz(ROOT/r['reference_path']);xy=data[r['sample_id']]['xy']
  assert sha(ROOT/r['path'])==r['file_sha256'];assert sha(ROOT/r['reference_path'])==r['reference_sha256'];assert canonical_hash(a)==r['graph_hash']
  stats=audit(a,reference,xy,r['graph_control_type'])
  for key in ['num_nodes','num_edges','degree_preserved','distance_distribution_preserved','graph_hash','changed_edge_fraction']:assert stats[key]==r[key],(p,key)
  if r['graph_control_type'] in ['A4','A5','A6']:
   rebuilt,_=randomize(reference,xy,r['graph_control_type'],r['derived_seed']);assert canonical_hash(rebuilt)==r['graph_hash'],p
  graphs+=1
  if graphs%20==0:print('RECONSTRUCTED',graphs,flush=True)
 paths=sorted((OUT/'logs/runs').glob('*.json'));assert len(paths)==480,len(paths)
 assert len(list((OUT/'models').glob('*.pt')))==480
 expected=set()
 for phase,heldouts,controls in [('primary',[None],CFG['primary_new_controls']),('lopo',IDS,CFG['conditional_lopo_new_controls'])]:
  for heldout in heldouts:
   for panel in PANELS:
    for control in controls:
     for seed in CFG['training_seeds']:
      for gs in ([-1] if control=='A1' else CFG['graph_randomization_seeds']):expected.add((phase,heldout,panel,control,seed,gs))
 actual=set()
 for path in paths:
  rec=json.loads(path.read_text());actual.add((rec['phase'],rec['heldout'],rec['panel'],rec['control'],rec['training_seed'],rec['graph_randomization_seed']))
  assert rec['params']==SELECTED[rec['panel']]['SAGE']['params']
  if rec['heldout']:assert rec['best_epoch']==SELECTED[rec['panel']]['SAGE']['lopo_epochs'] and rec['validation_slides']==[]
 assert actual==expected
 maxerr=0.;n=0
 for path in paths:
  r=json.loads(path.read_text());sample=r['heldout'] or SPLITS['splits']['test'][0];ev=json.loads((OUT/'logs/evaluation'/f'{r["run_id"]}_{sample}.json').read_text())
  assert sample not in r['train_slides'];assert sha(ROOT/r['preprocessing']/'feature_scaler.npz')==r['feature_scaler_sha256'];assert sha(ROOT/r['preprocessing']/'target_scaler.npz')==r['target_scaler_sha256']
  for s,g in r['graphs'].items():
   if r['control']!='A1':assert canonical_hash(sparse.load_npz(ROOT/g['path']))==g['graph_hash']
  x,g=run.features(sample,data[sample],r['panel'],r['control'],r['graph_randomization_seed'],r['heldout']);p=scaler(ROOT/r['preprocessing'],'target').inverse(core.predict(run.model_load(r),x,'B2'))
  assert sha(ROOT/ev['prediction_path'])==ev['prediction_sha256']
  with np.load(ROOT/ev['prediction_path']) as f:
   assert np.array_equal(f['barcodes'],data[sample]['barcodes']);assert list(f['genes'])==PANELS[r['panel']];error=float(np.max(abs(p-f['prediction'])))
  maxerr=max(maxerr,error);assert error<1e-7,(r['run_id'],error)
  expected=core.gene_metrics(data[sample]['y'][:,ix(r['panel'])],p);actual=pd.read_csv(ROOT/ev['gene_metrics_path'])
  for key,arr in expected.items():assert np.allclose(arr,actual[key],rtol=1e-10,atol=1e-10,equal_nan=True),(r['run_id'],key)
  n+=1
  if n%50==0:print('REPLAYED',n,'error',maxerr,flush=True)
 evals=[json.loads(p.read_text()) for p in (OUT/'logs/evaluation').glob('*.json')];assert len(evals)==660,len(evals)
 for r in evals:
  if 'reused_from' in r:
   frozen=pd.read_csv(ROOT/r['frozen_gene_metrics_path']);copied=pd.read_csv(ROOT/r['gene_metrics_path'])
   for key in ['pearson','spearman','mse','mae','prediction_variance','target_variance']:assert np.allclose(frozen[key],copied[key],rtol=1e-12,atol=1e-12,equal_nan=True)
   assert sha(ROOT/r['prediction_path'])==r['prediction_sha256']
 reg=pd.read_csv(ROOT/'results/experiment_registry.csv');new=reg[reg.run_id.str.startswith('C3A_')|reg.run_id.isin([p.stem.split('__',1)[1] for p in (OUT/'logs/registry').glob('*.json')])]
 assert not new.duplicated(['event','run_id']).any()
 events=list((OUT/'logs/registry').glob('*.json'));assert len(new)==len(events),(len(new),len(events))
 test=json.loads((OUT/'logs/TEST_ACCESS.json').read_text());lock=json.loads((OUT/'C3A_SPEC_LOCK.json').read_text());assert lock['timestamp']<test['timestamp'] and test['spec_sha256']==sha(OUT/'C3A_SPEC_LOCK.json')
 analysis=json.loads((OUT/'ANALYSIS_IMPLEMENTATION_LOCK.json').read_text());assert analysis['timestamp']<test['timestamp']
 for p,h in analysis['source_sha256'].items():assert sha(ROOT/p)==h
 result={'status':'PASS','new_models_replayed':n,'maximum_absolute_prediction_error':maxerr,'graphs_reconstructed':graphs,'frozen_scaler_training_only_checks':scaler_checks,'evaluations':len(evals),'reused_frozen_evaluations':180,'registry_events':len(events),'frozen_C0_C1_C2':'unchanged','C1_C2_models_rerun':0,'spec_predates_C3_test_evaluation':True,'seconds':time.perf_counter()-start}
 dump(OUT/'logs/INTEGRITY_AUDIT.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):verify()
