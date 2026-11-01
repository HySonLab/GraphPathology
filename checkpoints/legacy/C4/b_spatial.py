"""Minimal C4B recheck using frozen C2 operators and exact C3 distance-matched graphs."""
from common import *
import argparse,time,shutil,torch
from scipy import sparse
from threadpoolctl import threadpool_limits
from graph import aggregate
spec=importlib.util.spec_from_file_location('c4_runtime',OUT/'run.py');local=importlib.util.module_from_spec(spec);spec.loader.exec_module(local)
spec=importlib.util.spec_from_file_location('c3_graph_audit',ROOT/'C3/C3A/graph_controls.py');ga=importlib.util.module_from_spec(spec);spec.loader.exec_module(ga)
B=OUT/'C4B';ENC='pathology_encoder_1';SELECTED=json.loads((ROOT/'C2/SELECTION_LOCK.json').read_text())['selected'];GCACHE={}
def config():return json.loads((B/'configs/protocol.json').read_text())
def registry(event,run_id,record):
 r={**record,'event':event,'run_id':run_id,'checkpoint':'C4B','timestamp':now(),'source_commit':commit()};p=B/'logs/registry'/f'{event}__{run_id}.json';assert not p.exists();dump(p,r)
 fields=next(csv.reader((ROOT/'results/experiment_registry.csv').open()));row={**r,'model':r.get('control',''),'seed':r.get('training_seed',r.get('graph_seed','')),'threads':4,'cpu':'Intel Core i7-14700'}
 for k in ['train_slides','evaluation_slides']:
  if isinstance(row.get(k),list):row[k]='|'.join(row[k])
 with (ROOT/'results/experiment_registry.csv').open('a',newline='') as f:csv.DictWriter(f,fieldnames=fields,extrasaction='ignore').writerow(row)
def check_lock():
 lock=json.loads((B/'C4B_SPEC_LOCK.json').read_text());assert sha(B/'C4B_SPEC_LOCK.json')==(B/'C4B_SPEC_LOCK.sha256').read_text().split()[0]
 for p,h in lock['source_sha256'].items():assert sha(ROOT/p)==h,p
 assert sha(A/'logs/INTEGRITY_AUDIT.json')==lock['C4A_integrity_sha256'];assert sha(A/'logs/GATE_RESULT.json')==lock['C4A_gate_sha256'];return lock

def graph(sample,panel,control,gseed):
 family='MEAN' if control=='MEAN' else 'SAGE';geometry=SELECTED[panel][family]['graph']['name'];tag=f'{sample}_{geometry}_{control}_g{gseed}'
 if tag in GCACHE:return GCACHE[tag]
 ref=ROOT/'C2/graphs'/f'{sample}_{geometry}.npz';source=ROOT/'C3/C3A/graphs'/f'{sample}_{geometry}_A6_g{gseed}.npz' if control=='MATCHED' else ref;path=B/'graphs'/f'{tag}.npz';meta=path.with_suffix('.json')
 if meta.exists():
  r=json.loads(meta.read_text());assert sha(path)==r['file_sha256'] and sha(source)==r['source_sha256'];a=sparse.load_npz(path)
 else:
  a=sparse.load_npz(source);reference=sparse.load_npz(ref)
  with np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as f:xy=f['xy_um']
  stats=ga.audit(a,reference,xy,'A6' if control=='MATCHED' else 'A3')
  if control=='MATCHED':assert stats['changed_edge_fraction']>=.1 and stats['edge_length_ks_rounded_0_001_um']<=.2 and abs(stats['relative_mean_edge_length_change'])<=.1
  shutil.copyfile(source,path);r={'sample_id':sample,'geometry':geometry,'control':control,'graph_seed':gseed,'source_path':str(source.relative_to(ROOT)),'source_sha256':sha(source),'file_path':str(path.relative_to(ROOT)),'file_sha256':sha(path),'reference_path':str(ref.relative_to(ROOT)),'reference_sha256':sha(ref),**stats};dump(meta,r);registry('graph_reference','C4B_graph_'+tag,r)
 GCACHE[tag]=(a,r);return a,r

def folder(panel,heldout=None):return A/'preprocessing'/ENC/(heldout or 'primary')/panel
def features(sample,d,panel,control,gseed,heldout=None):
 x=scaler(folder(panel,heldout)/'feature.npz').transform(d['x']);a,r=graph(sample,panel,control,gseed);return aggregate(a,x,'MEAN' if control=='MEAN' else 'SAGE'),r

def fit(phase,panel,control,seed,gseed,data,heldout=None):
 rid=f'C4B_{phase}_{heldout or "primary"}_{panel}_{control}_s{seed}_g{gseed}';path=B/'logs/runs'/f'{rid}.json'
 if path.exists():return json.loads(path.read_text())
 train=[s for s in IDS if s!=heldout] if heldout else SPLITS['splits']['train'];val=[] if heldout else SPLITS['splits']['validation'];p=SELECTED[panel]['MEAN' if control=='MEAN' else 'SAGE'];ys=scaler(folder(panel,heldout)/'target.npz');refs={};t=time.perf_counter()
 def make(ids):
  z=[]
  for s in ids:x,r=features(s,data[s],panel,control,gseed,heldout);z.append(x);refs[s]=r
  return np.concatenate(z)
 x=make(train);y=ys.transform(np.concatenate([data[s]['y'][:,ix(panel)] for s in train]));xv=make(val) if val else None;yv=ys.transform(np.concatenate([data[s]['y'][:,ix(panel)] for s in val])) if val else None;prep=time.perf_counter()-t
 with core.PeakRSS() as mem:net,info=core.train_mlp(x,y,xv,yv,p['params'],base.C1CFG['mlp'],seed,p['lopo_epochs'] if heldout else None)
 history=info.pop('history');dump(B/'logs/epochs'/f'{rid}.json',history);mp=B/'models'/f'{rid}.pt';torch.save(net.state_dict(),mp)
 if val:
  pred=ys.inverse(core.predict(net,xv,'B2'));metrics=core.gene_metrics(np.concatenate([data[s]['y'][:,ix(panel)] for s in val]),pred);pd.DataFrame({'gene':PANELS[panel],**metrics}).to_csv(B/'logs/validation'/f'{rid}.csv',index=False)
 r={'run_id':rid,'phase':phase,'panel':panel,'control':control,'training_seed':seed,'graph_seed':gseed,'heldout':heldout,'train_slides':train,'validation_slides':val,'params':p['params'],'fixed_fold_epochs':p['lopo_epochs'] if heldout else None,'input_dim':x.shape[1],'output_dim':y.shape[1],'model_path':str(mp.relative_to(ROOT)),'model_sha256':sha(mp),'preprocessing':str(folder(panel,heldout).relative_to(ROOT)),'feature_scaler_sha256':sha(folder(panel,heldout)/'feature.npz'),'target_scaler_sha256':sha(folder(panel,heldout)/'target.npz'),'embedding_hashes':{s:data[s]['cache_sha256'] for s in train+val},'graphs':refs,'prepare_seconds':prep,'peak_rss_mb':mem.peak/1024**2,'source_commit':commit(),**info};dump(path,r);registry('fit',rid,r);print('FIT',rid,'epoch',r['best_epoch'],flush=True);return r

def model_load(r):
 assert sha(ROOT/r['model_path'])==r['model_sha256'];net=core.mlp_model(r['input_dim'],r['output_dim'],r['params']);net.load_state_dict(torch.load(ROOT/r['model_path'],map_location='cpu',weights_only=True));return net.eval()
def predict(r,sample,d):
 x,ref=features(sample,d,r['panel'],r['control'],r['graph_seed'],r['heldout']);return scaler(ROOT/r['preprocessing']/'target.npz').inverse(core.predict(model_load(r),x,'B2')),ref

def evaluate(r,sample,d):
 path=B/'logs/evaluation'/f'{r["run_id"]}_{sample}.json'
 if path.exists():return
 t=time.perf_counter();pred,ref=predict(r,sample,d);y=d['y'][:,ix(r['panel'])];m=core.gene_metrics(y,pred);co,_,_=core.coexpression(y,pred);pp=B/'predictions'/f'{r["run_id"]}_{sample}.npz';np.savez_compressed(pp,prediction=pred,barcodes=d['barcodes'],genes=np.array(PANELS[r['panel']]))
 cols={k:r[k] for k in ['run_id','phase','panel','control','training_seed','graph_seed']};cols['sample_id']=sample;gp=B/'logs/gene_metrics'/f'{r["run_id"]}_{sample}.csv';pd.DataFrame({'gene':PANELS[r['panel']],**m}).assign(**cols).to_csv(gp,index=False)
 out={**cols,**core.summarize(m),**co,'prediction_path':str(pp.relative_to(ROOT)),'prediction_sha256':sha(pp),'gene_metrics_path':str(gp.relative_to(ROOT)),'graph_hash':ref['graph_hash'],'embedding_sha256':d['cache_sha256'],'inference_seconds':time.perf_counter()-t};dump(path,out);registry('evaluation',r['run_id'],{**out,'evaluation_slides':[sample]});print('EVAL',r['run_id'],out['median_pearson'],flush=True)

def reuse_local():
 for f in (A/'logs/evaluation').glob('*.json'):
  r=json.loads(f.read_text())
  if r['encoder']!=ENC or r['model']!='D1' or r['panel'] not in config()['panels'] or r['phase'] not in ['primary','lopo']:continue
  rid='C4B_reuse_'+r['run_id'];path=B/'logs/evaluation'/f'{rid}_{r["sample_id"]}.json'
  if path.exists():continue
  assert sha(ROOT/r['prediction_path'])==r['prediction_sha256'];cols={'run_id':rid,'phase':r['phase'],'panel':r['panel'],'control':'LOCAL','training_seed':r['seed'],'graph_seed':-1,'sample_id':r['sample_id']};g=pd.read_csv(ROOT/r['gene_metrics_path'])[['gene','pearson','spearman','mse','mae','prediction_variance','target_variance']].assign(**cols);gp=B/'logs/gene_metrics'/f'{rid}_{r["sample_id"]}.csv';g.to_csv(gp,index=False);out={**r,**cols,'gene_metrics_path':str(gp.relative_to(ROOT)),'reused_from':str(f.relative_to(ROOT))};dump(path,out);registry('reuse_evaluation',rid,out)

def primary():
 data={s:local.load(ENC,s) for s in SPLITS['splits']['train']+SPLITS['splits']['validation']}
 for panel in config()['panels']:
  for control in ['MEAN','TRUE','MATCHED']:
   for gs in range(5) if control=='MATCHED' else [-1]:
    for seed in range(3):fit('primary',panel,control,seed,gs,data)
 paths=sorted((B/'logs/runs').glob('*.json'));assert len(paths)==42
 dump(B/'PRIMARY_LOCK.json',{'timestamp':now(),'spec_sha256':sha(B/'C4B_SPEC_LOCK.json'),'primary_model_records':{str(p.relative_to(ROOT)):sha(p) for p in paths},'no_B_test_predictions':True})

def final():
 lock=json.loads((B/'PRIMARY_LOCK.json').read_text())
 for p,h in lock['primary_model_records'].items():assert sha(ROOT/p)==h
 if not (B/'logs/TEST_ACCESS.json').exists():dump(B/'logs/TEST_ACCESS.json',{'timestamp':now(),'primary_lock_sha256':sha(B/'PRIMARY_LOCK.json'),'C4A_test_outcomes_already_known':True})
 data={s:local.load(ENC,s) for s in IDS};reuse_local()
 for f in sorted((B/'logs/runs').glob('*.json')):
  r=json.loads(f.read_text())
  if r['phase']=='primary':evaluate(r,SPLITS['splits']['test'][0],data[SPLITS['splits']['test'][0]])
 for sample in IDS:
  for panel in config()['panels']:
   for control in ['MEAN','TRUE','MATCHED']:
    for gs in range(5) if control=='MATCHED' else [-1]:
     for seed in range(3):r=fit('lopo',panel,control,seed,gs,data,sample);evaluate(r,sample,data[sample])
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['primary','final']);args=parser.parse_args();check_lock();torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):globals()[args.stage]()
