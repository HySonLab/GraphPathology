"""Replay C4B, verify exact graph reuse, schedules, references and frozen C4A artifacts."""
from common import *
import torch,time,io
from scipy import sparse
from threadpoolctl import threadpool_limits
spec=importlib.util.spec_from_file_location('c4_b_runtime',OUT/'b_spatial.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
B=b.B

def main():
 t=time.perf_counter();lock=b.check_lock()
 for p,h in lock['C4A_artifacts_sha256'].items():assert sha(ROOT/p)==h,p
 primary=json.loads((B/'PRIMARY_LOCK.json').read_text());access=json.loads((B/'logs/TEST_ACCESS.json').read_text());assert lock['timestamp']<primary['timestamp']<access['timestamp']
 data={s:b.local.load(b.ENC,s) for s in IDS};graphs=0
 for f in (B/'graphs').glob('*.json'):
  r=json.loads(f.read_text());assert sha(ROOT/r['file_path'])==r['file_sha256']==sha(ROOT/r['source_path']);assert sha(ROOT/r['reference_path'])==r['reference_sha256'];a=sparse.load_npz(ROOT/r['file_path']);ref=sparse.load_npz(ROOT/r['reference_path'])
  with np.load(ROOT/'data/cache'/f'{r["sample_id"]}_targets.npz') as z:xy=z['xy_um']
  stats=b.ga.audit(a,ref,xy,'A6' if r['control']=='MATCHED' else 'A3');assert stats['graph_hash']==r['graph_hash']
  if r['control']=='MATCHED':assert stats['changed_edge_fraction']>=.1 and stats['edge_length_ks_rounded_0_001_um']<=.2 and abs(stats['relative_mean_edge_length_change'])<=.1
  graphs+=1
 assert graphs==56,graphs
 records=[json.loads(f.read_text()) for f in (B/'logs/runs').glob('*.json')];assert len(records)==210;error=0.;vcount=0
 for i,r in enumerate(records):
  selected=b.SELECTED[r['panel']]['MEAN' if r['control']=='MEAN' else 'SAGE'];assert r['params']==selected['params'];folder=ROOT/r['preprocessing'];assert sha(folder/'feature.npz')==r['feature_scaler_sha256'] and sha(folder/'target.npz')==r['target_scaler_sha256']
  for sample,h in r['embedding_hashes'].items():assert data[sample]['cache_sha256']==h
  expected=[s for s in IDS if s!=r['heldout']] if r['heldout'] else SPLITS['splits']['train'];assert r['train_slides']==expected
  if r['heldout']:assert r['best_epoch']==selected['lopo_epochs'] and r['validation_slides']==[]
  else:
   assert r['validation_slides']==SPLITS['splits']['validation'];sample=r['validation_slides'][0];pred,_=b.predict(r,sample,data[sample]);m=core.gene_metrics(data[sample]['y'][:,ix(r['panel'])],pred);stored=pd.read_csv(B/'logs/validation'/f'{r["run_id"]}.csv')
   for k,v in m.items():assert np.allclose(v,stored[k],rtol=1e-9,atol=1e-9,equal_nan=True)
   vcount+=1
  sample=r['heldout'] or SPLITS['splits']['test'][0];assert sample not in r['train_slides'];pred,gr=b.predict(r,sample,data[sample]);ev=json.loads((B/'logs/evaluation'/f'{r["run_id"]}_{sample}.json').read_text());assert sha(ROOT/ev['prediction_path'])==ev['prediction_sha256'];assert ev['graph_hash']==gr['graph_hash'] and ev['embedding_sha256']==data[sample]['cache_sha256']
  with np.load(ROOT/ev['prediction_path']) as z:err=float(abs(pred-z['prediction']).max());assert err<1e-7;error=max(error,err);assert np.array_equal(z['barcodes'],data[sample]['barcodes']) and list(z['genes'])==PANELS[r['panel']]
  m=core.gene_metrics(data[sample]['y'][:,ix(r['panel'])],pred);stored=pd.read_csv(ROOT/ev['gene_metrics_path'])
  for k,v in m.items():assert np.allclose(v,stored[k],rtol=1e-9,atol=1e-9,equal_nan=True),(r['run_id'],k)
  if (i+1)%40==0:print('REPLAY',i+1,flush=True)
 evs=[json.loads(f.read_text()) for f in (B/'logs/evaluation').glob('*.json')];assert len(evs)==240;reused=0;csv_error=0.
 for ev in evs:
  if ev['control']!='LOCAL':continue
  source=json.loads((ROOT/ev['reused_from']).read_text());assert ev['prediction_sha256']==source['prediction_sha256']==sha(ROOT/ev['prediction_path']);one=pd.read_csv(ROOT/ev['gene_metrics_path']);two=pd.read_csv(ROOT/source['gene_metrics_path'])
  for col in ['pearson','spearman','mse','mae','target_variance','prediction_variance']:
   err=float(np.nanmax(abs(one[col]-two[col])));csv_error=max(csv_error,err);assert np.allclose(one[col],two[col],equal_nan=True,rtol=0,atol=1e-14)
  reused+=1
 assert reused==30;events=list((B/'logs/registry').glob('*.json'));reg=pd.read_csv(ROOT/'results/experiment_registry.csv');reg=reg[reg.run_id.str.startswith('C4B_')];assert len(reg)==len(events)==506 and not reg.duplicated(['event','run_id']).any()
 # Verify every C4A registry row still equals its frozen lossless event record.
 with (ROOT/'results/experiment_registry.csv').open(newline='') as f:
  reader=csv.DictReader(f);fields=reader.fieldnames;actual={(r['event'],r['run_id']):r for r in reader if r['run_id'].startswith('C4A_')}
 assert len(actual)==420
 for f in (A/'logs/registry').glob('*.json'):
  r=json.loads(f.read_text());row={**r,'model':r.get('model',r.get('encoder','')),'cpu':'Intel Core i7-14700','threads':r.get('threads',4)}
  for k in ['train_slides','evaluation_slides']:
   if isinstance(row.get(k),list):row[k]='|'.join(row[k])
  buffer=io.StringIO();writer=csv.DictWriter(buffer,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerow(row);buffer.seek(0);expected=next(csv.DictReader(buffer));assert expected==actual[(r['event'],r['run_id'])],r['run_id']
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C3/freeze.py'),'--verify'],check=True)
 result={'status':'PASS','saved_models_replayed':len(records),'heldout_predictions_replayed':210,'primary_validation_replays':vcount,'reused_local_evaluations_checked':reused,'reused_metric_CSV_roundtrip_max_error':csv_error,'reused_prediction_files_byte_identical':True,'maximum_prediction_error':error,'graphs_checked':graphs,'registry_events':len(events),'C4A_frozen_artifacts_checked':len(lock['C4A_artifacts_sha256']),'C0_C3':'unchanged','seconds':time.perf_counter()-t};dump(B/'logs/INTEGRITY_AUDIT.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):main()
