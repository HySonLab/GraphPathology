"""C2 finite validation graph selection, frozen-test evaluation and conditional LOPO."""
import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','4');os.environ.setdefault('MKL_NUM_THREADS','4')
import argparse,csv,hashlib,importlib.util,json,subprocess,time,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
import torch
from threadpoolctl import threadpool_limits
from graph import build_graph,aggregate,shuffle_graph,hidden_width,parameter_count
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C2'
spec=importlib.util.spec_from_file_location('frozen_c1_core',ROOT/'C1/core.py');core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
CFG=json.loads((OUT/'configs/protocol.json').read_text());C1CFG=json.loads((ROOT/'C1/configs/protocol.json').read_text());C1LOCK=json.loads((ROOT/'C1/SELECTION_LOCK.json').read_text())
PANELS=C1LOCK['panels'];SPLITS=json.loads((ROOT/'reports/C0/C0_SPLITS.json').read_text());IDS=sum(SPLITS['splits'].values(),[]);GENES=PANELS['S2-280']
FIELDS=next(csv.reader((ROOT/'results/experiment_registry.csv').open()))
GRAPH_CACHE={};INPUT_CACHE={}

def sha(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def safe(value):
 if isinstance(value,dict):return {k:safe(v) for k,v in value.items()}
 if isinstance(value,(list,tuple)):return [safe(v) for v in value]
 if isinstance(value,(np.integer,np.floating)):value=value.item()
 if isinstance(value,float) and not np.isfinite(value):return None
 return value
def dump(path,value):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(safe(value),indent=2,allow_nan=False)+'\n')
def commit():return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
def registry(row):
 path=ROOT/'results/experiment_registry.csv'
 with path.open(newline='') as f:
  reader=csv.DictReader(f)
  assert reader.fieldnames==FIELDS
  if any(r['event']==row['event'] and r['run_id']==row['run_id'] for r in reader):raise ValueError('Duplicate registry event')
 row={**row,'cpu':'Intel Core i7-14700','threads':4,'source_commit':commit(),'timestamp':datetime.now(timezone.utc).isoformat()}
 with path.open('a',newline='') as f:csv.DictWriter(f,fieldnames=FIELDS,extrasaction='ignore').writerow(row)

def load(sample):
 with np.load(ROOT/'data/cache'/f'{sample}_features.npz') as f:x=f['embeddings'].copy();barcodes=f['barcodes'].copy()
 with np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as f:
  assert list(f['genes'])==GENES and np.array_equal(barcodes,f['barcodes'])
  y=f['normalized_targets'].copy();xy=f['xy_um'].copy()
 assert np.isfinite(x).all() and np.isfinite(y).all()
 return {'x':x,'y':y,'xy':xy,'barcodes':barcodes}
def ix(panel):return [GENES.index(g) for g in PANELS[panel]]
def scaler(folder,name):
 s=core.Scaler()
 with np.load(folder/f'{name}_scaler.npz') as f:s.mean=f['mean'];s.scale=f['scale']
 return s
def preprocessing(panel,heldout=None):return ROOT/'C1'/('preprocessing/primary' if heldout is None else f'lopo/{heldout}')/panel

def graph(sample,data,config):
 key=(sample,config['name'])
 if key in GRAPH_CACHE:return GRAPH_CACHE[key]
 path=OUT/'graphs'/f'{sample}_{config["name"]}.npz';meta=path.with_suffix('.json')
 if path.exists():
  record=json.loads(meta.read_text());assert sha(path)==record['sha256'];a=sparse.load_npz(path)
 else:
  t=time.perf_counter();a=build_graph(data['xy'],config);seconds=time.perf_counter()-t;sparse.save_npz(path,a)
  i,j=a.nonzero();length=np.linalg.norm(data['xy'][i]-data['xy'][j],axis=1);degree=np.asarray(a.sum(1)).ravel()
  record={'sample_id':sample,'config':config,'n_nodes':len(data['x']),'undirected_edges':a.nnz//2,
          'degree_min':int(degree.min()),'degree_mean':float(degree.mean()),'degree_max':int(degree.max()),
          'isolates':int((degree==0).sum()),'components':int(connected_components(a,directed=False,return_labels=False)),
          'edge_um_median':float(np.median(length)) if len(length) else None,'edge_um_max':float(length.max()) if len(length) else None,
          'construction_seconds':seconds,'sha256':sha(path),'barcodes_sha256':hashlib.sha256('\n'.join(data['barcodes']).encode()).hexdigest(),
          'coordinates_sha256':hashlib.sha256(data['xy'].tobytes()).hexdigest(),'labels_used':False}
  dump(meta,record)
 GRAPH_CACHE[key]=a;return a

def features(sample,data,panel,model,config,heldout=None,shuffle_seed=None):
 key=(sample,panel,model,config['name'],heldout,shuffle_seed)
 if key in INPUT_CACHE:return INPUT_CACHE[key]
 a=graph(sample,data,config)
 if shuffle_seed is not None:a,_=shuffle_graph(a,sample,shuffle_seed)
 x=scaler(preprocessing(panel,heldout),'feature').transform(data['x'])
 z=aggregate(a,x,model)
 # Cache only one panel/run cohort to bound memory during graph search.
 INPUT_CACHE[key]=z;return z

def fit(run_id,phase,panel,model,gconfig,seed,data,train_ids,val_ids,heldout=None,shuffle=False,fixed_epochs=None):
 path=OUT/'logs/runs'/f'{run_id}.json'
 if path.exists():raise ValueError('Run already completed')
 INPUT_CACHE.clear()
 folder=preprocessing(panel,heldout);ys=scaler(folder,'target');indices=ix(panel)
 start=time.perf_counter()
 x=np.concatenate([features(s,data[s],panel,model,gconfig,heldout,seed if shuffle else None) for s in train_ids])
 y=ys.transform(np.concatenate([data[s]['y'][:,indices] for s in train_ids]))
 if val_ids:
  xv=np.concatenate([features(s,data[s],panel,model,gconfig,heldout,seed if shuffle else None) for s in val_ids]);yv=ys.transform(np.concatenate([data[s]['y'][:,indices] for s in val_ids]))
 else:xv=yv=None
 prepare_seconds=time.perf_counter()-start
 local=C1LOCK['selected'][panel]['mlp'];params={'hidden':hidden_width(512,len(indices),local['hidden'],model),'dropout':local['dropout']}
 with core.PeakRSS() as mem:
  net,info=core.train_mlp(x,y,xv,yv,params,C1CFG['mlp'],seed,fixed_epochs)
  if val_ids:
   start=time.perf_counter();vp=ys.inverse(core.predict(net,xv,'B2'));info['validation_seconds']+=time.perf_counter()-start
   truth=np.concatenate([data[s]['y'][:,indices] for s in val_ids]);metrics=core.gene_metrics(truth,vp);summary=core.summarize(metrics)
   rows=pd.DataFrame({'gene':PANELS[panel],**metrics});vpath=OUT/'logs/validation'/f'{run_id}.csv';vpath.parent.mkdir(exist_ok=True);rows.to_csv(vpath,index=False)
  else:summary={}
  model_path=OUT/'models'/f'{run_id}.pt';torch.save(net.state_dict(),model_path)
 history=info.pop('history');dump(OUT/'logs/epochs'/f'{run_id}.json',history)
 record={'run_id':run_id,'phase':phase,'panel':panel,'model':model,'graph':gconfig,'seed':seed,'params':params,
         'train_slides':train_ids,'validation_slides':val_ids,'heldout':heldout,'shuffled':shuffle,
         'preprocessing':str(folder.relative_to(ROOT)),'feature_scaler_sha256':sha(folder/'feature_scaler.npz'),
         'target_scaler_sha256':sha(folder/'target_scaler.npz'),'model_path':str(model_path.relative_to(ROOT)),
         'model_sha256':sha(model_path),'input_dim':x.shape[1],'output_dim':y.shape[1],'n_train':len(x),
         'local_parameters':parameter_count(512,local['hidden'],len(indices)),
         'source_commit':commit(),'prepare_seconds':prepare_seconds,'peak_rss_mb':mem.peak/1024**2,**info,
         'validation_median_pearson':summary.get('median_pearson'),'validation_mse':summary.get('mean_mse')}
 dump(path,record)
 registry({**record,'event':'fit','train_slides':'|'.join(train_ids),'evaluation_slides':'|'.join(val_ids or [heldout]),
           'median_pearson':record['validation_median_pearson'],'mean_mse':record['validation_mse']})
 print(run_id,'val PCC',record['validation_median_pearson'],'epoch',info['best_epoch'],flush=True)
 INPUT_CACHE.clear();return record

def candidate_genes(main):
 records=[]
 for panel in CFG['panels']:
  if panel=='S1':continue
  for model in ['GCN','SAGE']:
   frames=[]
   for seed in CFG['seeds']:
    a=pd.read_csv(OUT/'logs/validation'/f'C2_primary_{panel}_{model}_s{seed}.csv').set_index('gene')
    b=pd.read_csv(OUT/'logs/validation'/f'C2_primary_{panel}_MEAN_s{seed}.csv').set_index('gene')
    assert a.index.equals(b.index)
    frames.append(pd.DataFrame({'gene':a.index,'seed':seed,'gain_over_mean':(a.pearson-b.pearson).values,'mse_change':(a.mse-b.mse).values}))
   table=pd.concat(frames).groupby('gene').agg(mean_gain=('gain_over_mean','mean'),minimum_seed_gain=('gain_over_mean','min'),mse_change=('mse_change','mean')).reset_index()
   table['validation_candidate']=(table.mean_gain>=.05)&(table.minimum_seed_gain>0)&(table.mse_change<=0)
   table['panel']=panel;table['model']=model;records.extend(table.to_dict('records'))
 pd.DataFrame(records).to_csv(OUT/'C2_VALIDATION_GENE_SHORTLIST.csv',index=False)
 return [r for r in records if r['validation_candidate']]

def select():
 if (OUT/'SELECTION_LOCK.json').exists():raise ValueError('Selection already frozen')
 train=SPLITS['splits']['train'];val=SPLITS['splits']['validation'];data={s:load(s) for s in train+val}
 selected={};main=[];tune=[]
 for panel in CFG['panels']:
  selected[panel]={}
  for model in CFG['models']:
   records=[]
   for g in CFG['graphs']:
    r=fit(f'C2_tune_{panel}_{model}_{g["name"]}','selection',panel,model,g,0,data,train,val);records.append(r);tune.append(r['run_id'])
   best=min(records,key=lambda r:(-r['validation_median_pearson'],r['validation_mse'],r['graph']['name']))
   selected[panel][model]={'graph':best['graph'],'params':best['params']}
   epochs=[]
   for seed in CFG['seeds']:
    r=fit(f'C2_primary_{panel}_{model}_s{seed}','primary',panel,model,best['graph'],seed,data,train,val);main.append(r['run_id']);epochs.append(r['best_epoch'])
   selected[panel][model]['lopo_epochs']=int(np.floor(np.median(epochs)+.5))
   if panel=='S2-50':
    for seed in CFG['seeds']:
     r=fit(f'C2_shuffle_{panel}_{model}_s{seed}','shuffled_graph',panel,model,best['graph'],seed,data,train,val,shuffle=True);main.append(r['run_id'])
 shortlist=candidate_genes(main)
 lock={'timestamp':datetime.now(timezone.utc).isoformat(),'source_commit':commit(),'test_loaded':False,'selected':selected,
       'primary_runs':main,'selection_runs':tune,'validation_shortlist':shortlist,
       'validation_shortlist_sha256':sha(OUT/'C2_VALIDATION_GENE_SHORTLIST.csv'),
       'c1_lock_sha256':sha(ROOT/'C1/C1_SPEC_LOCK.json'),'protocol_sha256':sha(OUT/'configs/protocol.json'),
       'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'run.py',OUT/'graph.py',ROOT/'C1/core.py']}}
 dump(OUT/'SELECTION_LOCK.json',lock);print('C2 SELECTION LOCKED; no test cache loaded',flush=True)

def model_load(r):
 path=ROOT/r['model_path'];assert sha(path)==r['model_sha256']
 net=core.mlp_model(r['input_dim'],r['output_dim'],r['params']);net.load_state_dict(torch.load(path,map_location='cpu',weights_only=True));net.eval();return net

def evaluate_one(r,sample,d,lock_sha):
 path=OUT/'logs/evaluation'/f'{r["run_id"]}_{sample}.json'
 if path.exists():raise ValueError('Evaluation already exists')
 with core.PeakRSS() as mem:
  start=time.perf_counter();x=features(sample,d,r['panel'],r['model'],r['graph'],r['heldout'],r['seed'] if r['shuffled'] else None)
  ys=scaler(ROOT/r['preprocessing'],'target');p=ys.inverse(core.predict(model_load(r),x,'B2'));seconds=time.perf_counter()-start
  y=d['y'][:,ix(r['panel'])];metrics=core.gene_metrics(y,p);summary=core.summarize(metrics)
 ppath=OUT/'predictions'/f'{r["run_id"]}_{sample}.npz';np.savez_compressed(ppath,prediction=p.astype(np.float32),genes=np.array(PANELS[r['panel']]),barcodes=d['barcodes'])
 rows=pd.DataFrame({'gene':PANELS[r['panel']],**metrics})
 for k,v in {'run_id':r['run_id'],'phase':r['phase'],'panel':r['panel'],'model':r['model'],'seed':r['seed'],'sample_id':sample}.items():rows[k]=v
 gpath=OUT/'logs/gene_metrics'/f'{r["run_id"]}_{sample}.csv';rows.to_csv(gpath,index=False)
 record={**{k:r[k] for k in ['run_id','phase','panel','model','seed']},'graph_name':r['graph']['name'],'sample_id':sample,'patient_id':SPLITS['patients'][sample],
         **summary,'n_locations':len(p),'inference_seconds':seconds,'inference_peak_rss_mb':mem.peak/1024**2,
         'prediction_path':str(ppath.relative_to(ROOT)),'prediction_sha256':sha(ppath),'gene_metrics_path':str(gpath.relative_to(ROOT)),
         'selection_lock_sha256':lock_sha}
 dump(path,record);registry({**r,**summary,'event':'evaluation','train_slides':'|'.join(r['train_slides']),'evaluation_slides':sample,'inference_seconds':seconds,'peak_rss_mb':mem.peak/1024**2})
 print('EVAL',r['run_id'],sample,'PCC',summary['median_pearson'],'MSE',summary['mean_mse'],flush=True)
 INPUT_CACHE.clear();return record

def check_selection():
 lock=json.loads((OUT/'SELECTION_LOCK.json').read_text())
 for p,digest in lock['source_sha256'].items():assert sha(ROOT/p)==digest,p
 assert sha(ROOT/'C1/C1_SPEC_LOCK.json')==lock['c1_lock_sha256']
 assert sha(OUT/'configs/protocol.json')==lock['protocol_sha256']
 return lock,sha(OUT/'SELECTION_LOCK.json')
def evaluate():
 lock,digest=check_selection();test=SPLITS['splits']['test'];data={s:load(s) for s in test}
 dump(OUT/'logs/TEST_ACCESS.json',{'timestamp':datetime.now(timezone.utc).isoformat(),'selection_lock_sha256':digest})
 for run_id in lock['primary_runs']:
  r=json.loads((OUT/'logs/runs'/f'{run_id}.json').read_text())
  for sample in test:evaluate_one(r,sample,data[sample],digest)

def lopo():
 lock,digest=check_selection();data={s:load(s) for s in IDS}
 for heldout in IDS:
  train=[s for s in IDS if s!=heldout]
  for panel in CFG['panels']:
   for model in CFG['models']:
    chosen=lock['selected'][panel][model]
    for seed in CFG['seeds']:
     r=fit(f'C2_lopo_{heldout}_{panel}_{model}_s{seed}','lopo',panel,model,chosen['graph'],seed,data,train,[],heldout=heldout,fixed_epochs=chosen['lopo_epochs'])
     evaluate_one(r,heldout,data[heldout],digest)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['select','evaluate','lopo']);args=parser.parse_args()
 for folder in ['graphs','models','predictions','logs/runs','logs/evaluation','logs/gene_metrics']:(OUT/folder).mkdir(parents=True,exist_ok=True)
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):globals()[args.stage]()
if __name__=='__main__':main()
