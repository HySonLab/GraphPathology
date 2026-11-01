"""C4 audit: saved-model replay, cache/order/patch checks and frozen-parent hashes."""
from common import *
import time,h5py,torch
from threadpoolctl import threadpool_limits
spec=importlib.util.spec_from_file_location('c4_runtime',OUT/'run.py');run=importlib.util.module_from_spec(spec);spec.loader.exec_module(run)

def verify():
 t=time.perf_counter();subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C3/freeze.py'),'--verify'],check=True)
 lock=json.loads((A/'C4A_SPEC_LOCK.json').read_text());assert sha(A/'C4A_SPEC_LOCK.json')==(A/'C4A_SPEC_LOCK.sha256').read_text().split()[0]
 for p,h in lock['source_sha256'].items():assert sha(ROOT/p)==h,p
 selected=json.loads((A/'SELECTION_LOCK.json').read_text());test=json.loads((A/'logs/TEST_ACCESS.json').read_text());assert lock['timestamp']<selected['timestamp']<test['timestamp']
 files=json.loads((A/'configs/encoder_files.json').read_text())
 for r in files['files']:assert sha(OUT/'weights/phikon'/r['name'])==r['sha256']
 cache_checks=0;patch_checks=0
 for sample in IDS:
  meta=json.loads((A/'logs/extraction'/f'{sample}.json').read_text());assert sha(ROOT/meta['cache_path'])==meta['cache_sha256'];assert sha(OUT/'embeddings/resnet18_reference'/f'{sample}.npz')==meta['reference_cache_sha256']
  with np.load(OUT/'embeddings/pathology_encoder_1'/f'{sample}.npz') as p,np.load(OUT/'embeddings/resnet18_reference'/f'{sample}.npz') as r,np.load(ROOT/'data/cache'/f'{sample}_features.npz') as old,np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as targets,h5py.File(ROOT/'data/raw/IDC/patches'/f'{sample}.h5') as patches:
   assert np.array_equal(r['embeddings'],old['embeddings']);assert np.array_equal(p['barcodes'],r['barcodes']) and np.array_equal(r['barcodes'],targets['barcodes']);assert np.array_equal(p['xy_um'],targets['xy_um']);assert p['embeddings'].shape==(len(p['barcodes']),768) and np.isfinite(p['embeddings']).all()
   assert np.array_equal(p['source_patch_sha256'],r['source_patch_sha256']);saved=p['source_patch_sha256']
   for begin in range(0,len(saved),256):
    actual=[hashlib.sha256(im.tobytes()).hexdigest() for im in patches['img'][begin:begin+256]];assert actual==list(saved[begin:begin+256]);patch_checks+=len(actual)
  cache_checks+=2
 scaler_checks=0;replays=0;error=0.;records=[json.loads(f.read_text()) for f in (A/'logs/runs').glob('*.json')]
 assert len(records)==244,len(records)
 for encoder in run.CFG['encoders']:
  data={s:run.load(encoder,s) for s in IDS}
  for heldout in [None]+IDS:
   train=[s for s in IDS if s!=heldout] if heldout else SPLITS['splits']['train']
   for panel in PANELS:
    folder=A/'preprocessing'/encoder/(heldout or 'primary')/panel
    for name,values in [('feature',np.concatenate([data[s]['x'] for s in train])),('target',np.concatenate([data[s]['y'][:,ix(panel)] for s in train]))]:
     expected=core.Scaler().fit(values);actual=scaler(folder/f'{name}.npz');assert np.allclose(expected.mean,actual.mean,rtol=1e-7,atol=1e-7) and np.allclose(expected.scale,actual.scale,rtol=1e-7,atol=1e-7);scaler_checks+=1
  for r in [r for r in records if r['encoder']==encoder]:
   for sample,digest in r['cache_hashes'].items():assert data[sample]['cache_sha256']==digest
   folder=ROOT/r['preprocessing'];assert sha(folder/'feature.npz')==r['feature_scaler_sha256'] and sha(folder/'target.npz')==r['target_scaler_sha256']
   if r['phase']=='selection':
    sample=SPLITS['splits']['validation'][0];p=run.predict(r,data[sample]);m=core.gene_metrics(data[sample]['y'][:,ix(r['panel'])],p);stored=pd.read_csv(A/'logs/validation'/f'{r["run_id"]}.csv')
   else:
    sample=r['heldout'] or SPLITS['splits']['test'][0];assert sample not in r['train_slides'];p=run.predict(r,data[sample]);ev=json.loads((A/'logs/evaluation'/f'{r["run_id"]}_{sample}.json').read_text());assert sha(ROOT/ev['prediction_path'])==ev['prediction_sha256']
    with np.load(ROOT/ev['prediction_path']) as f:delta=float(abs(f['prediction']-p).max());assert delta<1e-7;error=max(error,delta)
    m=core.gene_metrics(data[sample]['y'][:,ix(r['panel'])],p);stored=pd.read_csv(ROOT/ev['gene_metrics_path'])
   for k,v in m.items():assert np.allclose(v,stored[k],rtol=1e-9,atol=1e-9,equal_nan=True),(r['run_id'],k)
   if r['control']:
    path=A/'models'/f'{r["run_id"]}_permutation.npy';assert sha(path)==r['permutation_sha256'];perm=np.load(path);assert np.array_equal(np.sort(perm),np.arange(r['n_train']))
   replays+=1
   if replays%40==0:print('REPLAY',replays,'max error',error,flush=True)
 # Replay fixed embedding rows with the same batch partition (first batch of each slide).
 extraction_spec=importlib.util.spec_from_file_location('c4_extraction',OUT/'extract.py');extract=importlib.util.module_from_spec(extraction_spec);extraction_spec.loader.exec_module(extract);processor,net=extract.setup();embedding_error=0.
 for sample in IDS:
  with h5py.File(ROOT/'data/raw/IDC/patches'/f'{sample}.h5') as f:images=list(f['img'][:16])
  with torch.inference_mode():z=net(processor(images,return_tensors='pt').pixel_values).last_hidden_state[:,0].numpy()
  with np.load(OUT/'embeddings/pathology_encoder_1'/f'{sample}.npz') as f:err=float(abs(z-f['embeddings'][:16]).max());assert err<1e-6;embedding_error=max(embedding_error,err)
  print('EMBEDDING REPLAY',sample,err,flush=True)
 evals=list((A/'logs/evaluation').glob('*.json'));assert len(evals)==172,len(evals)
 events=list((A/'logs/registry').glob('*.json'));reg=pd.read_csv(ROOT/'results/experiment_registry.csv');reg=reg[reg.run_id.str.startswith('C4A_')];assert len(events)==len(reg) and not reg.duplicated(['event','run_id']).any()
 result={'status':'PASS','saved_models_replayed':replays,'evaluation_predictions_replayed':172,'selection_metric_replays':72,'maximum_prediction_error':error,'sampled_embedding_max_error':embedding_error,'sampled_embedding_rows':64,'patch_hashes_verified':patch_checks,'embedding_caches_checked':cache_checks,'scaler_training_only_checks':scaler_checks,'registry_events':len(events),'frozen_C0_C3':'unchanged','seconds':time.perf_counter()-t};dump(A/'logs/INTEGRITY_AUDIT.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):verify()
