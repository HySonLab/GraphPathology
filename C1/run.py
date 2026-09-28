"""C1 staged execution: select/freeze on train+validation, then evaluate, then LOPO."""
import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
os.environ.setdefault('MKL_NUM_THREADS','4')
import argparse,csv,hashlib,json,platform,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits
from core import Scaler,PeakRSS,train_mlp,predict,permute_training,column_corr,finite_stat,gene_metrics,summarize,coexpression,mlp_model
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C1'
CFG=json.loads((OUT/'configs/protocol.json').read_text())
MANIFEST=json.loads((ROOT/'reports/C0/C0_DATA_MANIFEST.json').read_text())
SPLITS=json.loads((ROOT/'reports/C0/C0_SPLITS.json').read_text())
PANELS0=json.loads((ROOT/'reports/C0/C0_GENE_PANELS.json').read_text())['panels']
GENES=json.loads((ROOT/'reports/C0/C0_NORMALIZATION.json').read_text())['genes']
PANELS={'S1':PANELS0['1']['genes'],'S2-50':PANELS0['50']['genes'],'S2-250':PANELS0['250']['genes'],'S2-280':GENES}
PATIENTS=SPLITS['patients'];IDS=sum(SPLITS['splits'].values(),[])
CPU=next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),platform.processor())
FIELDS=['event','run_id','phase','panel','model','seed','train_slides','evaluation_slides','train_seconds','validation_seconds','inference_seconds','peak_rss_mb','parameters','input_dim','output_dim','n_train','cpu','threads','best_epoch','median_pearson','mean_mse','source_commit','timestamp']

def sha(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def safe(v):
 if isinstance(v,dict):return {k:safe(x) for k,x in v.items()}
 if isinstance(v,list):return [safe(x) for x in v]
 if isinstance(v,(np.integer,np.floating)):v=v.item()
 if isinstance(v,float) and not np.isfinite(v):return None
 return v
def dump(path,obj):
 Path(path).parent.mkdir(parents=True,exist_ok=True);Path(path).write_text(json.dumps(safe(obj),indent=2,allow_nan=False)+'\n')
def commit():return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
def registry(row):
 path=ROOT/'results/experiment_registry.csv';path.parent.mkdir(exist_ok=True)
 row={**row,'cpu':CPU,'threads':CFG['threads'],'source_commit':commit(),'timestamp':datetime.now(timezone.utc).isoformat()}
 if path.exists():
  with path.open(newline='') as f:
   reader=csv.DictReader(f)
   if reader.fieldnames!=FIELDS:raise ValueError('Registry schema differs; refusing to overwrite')
   if any(r['event']==row['event'] and r['run_id']==row['run_id'] for r in reader):raise ValueError('Duplicate registry event')
 with path.open('a',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=FIELDS,extrasaction='ignore')
  if f.tell()==0:writer.writeheader()
  writer.writerow(row)

def load(sample):
 record=next(r for r in MANIFEST['slides'] if r['sample_id']==sample)
 with np.load(ROOT/record['feature_cache'],allow_pickle=False) as f:
  x=f['embeddings'].copy();barcodes=f['barcodes'].copy()
 with np.load(ROOT/record['target_cache'],allow_pickle=False) as f:
  if not np.array_equal(barcodes,f['barcodes']) or list(f['genes'])!=GENES:raise ValueError('C0 cache alignment mismatch')
  y=f['normalized_targets'].copy();raw=f['raw_counts'].copy();xy=f['xy_um'].copy()
 if not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('Nonfinite C0 cache')
 return {'x':x,'y':y,'raw':raw,'xy':xy,'barcodes':barcodes}
def combine(data,ids,key):return np.concatenate([data[s][key] for s in ids])
def panel_indices(panel):return [GENES.index(g) for g in PANELS[panel]]
def scaler_pair(x,y,folder):
 folder.mkdir(parents=True,exist_ok=True)
 xs,ys=Scaler().fit(x),Scaler().fit(y)
 xs.save(folder/'feature_scaler.npz');ys.save(folder/'target_scaler.npz')
 return xs,ys

def fit_run(run_id,phase,panel,kind,seed,x,y,xv,yv,train_ids,eval_ids,params,preprocess,fixed_epochs=None,negative=None):
 record_path=OUT/'logs/runs'/f'{run_id}.json'
 if record_path.exists():raise RuntimeError(f'Run already exists: {run_id}')
 permutation=None
 if negative:x,y,permutation=permute_training(x,y,negative,seed)
 with PeakRSS() as memory:
  start=time.perf_counter()
  if kind=='B0':
   model=y.mean(0);info={'train_seconds':time.perf_counter()-start,'validation_seconds':0.,'parameters':y.shape[1],'best_epoch':0}
  elif kind=='B1':
   model=Ridge(alpha=params['alpha'],solver='cholesky').fit(x.astype(np.float64),y.astype(np.float64))
   info={'train_seconds':time.perf_counter()-start,'validation_seconds':0.,'parameters':int(model.coef_.size+model.intercept_.size),'best_epoch':0}
  else:model,info=train_mlp(x,y,xv,yv,params,CFG['mlp'],seed,fixed_epochs)
  if xv is not None:
   start=time.perf_counter();vp=predict(model,xv,kind);info['validation_seconds']+=time.perf_counter()-start
   info['validation_median_pearson']=finite_stat(column_corr(yv,vp),'median')
   info['validation_mse']=float(np.mean((yv-vp)**2))
  else:info.update(validation_median_pearson=None,validation_mse=None)
  model_path=OUT/'models'/f'{run_id}.npz'
  if kind=='B2':
   model_path=model_path.with_suffix('.pt');torch.save(model.state_dict(),model_path)
  elif kind=='B1':np.savez_compressed(model_path,coef=model.coef_,intercept=model.intercept_)
  else:np.savez_compressed(model_path,mean=model)
 info['peak_rss_mb']=memory.peak/1024**2
 if 'history' in info:dump(OUT/'logs/epochs'/f'{run_id}.json',info.pop('history'))
 record={'run_id':run_id,'phase':phase,'panel':panel,'model':kind,'seed':seed,'train_slides':train_ids,
         'evaluation_slides':eval_ids,'params':params,'preprocessing':str(preprocess.relative_to(ROOT)),
         'model_path':str(model_path.relative_to(ROOT)),'model_sha256':sha(model_path),
         'input_dim':x.shape[1],'output_dim':y.shape[1],'n_train':len(x),'negative_control':negative,
         'source_commit':commit(),**info}
 if permutation is not None:
  perm_path=OUT/'logs/permutations'/f'{run_id}.npy';perm_path.parent.mkdir(exist_ok=True);np.save(perm_path,permutation)
  record['permutation_path']=str(perm_path.relative_to(ROOT));record['permutation_sha256']=sha(perm_path)
 dump(record_path,record)
 registry({**record,'event':'fit','train_slides':'|'.join(train_ids),'evaluation_slides':'|'.join(eval_ids),
           'median_pearson':info['validation_median_pearson'],'mean_mse':info['validation_mse']})
 print(f'{run_id}: {info["train_seconds"]:.2f}s, epoch {info["best_epoch"]}, val PCC {info["validation_median_pearson"]}',flush=True)
 return record

def select():
 if (OUT/'SELECTION_LOCK.json').exists():raise RuntimeError('Selection already frozen')
 train_ids=SPLITS['splits']['train'];val_ids=SPLITS['splits']['validation']
 # Deliberately do not load the test cache in this stage.
 data={s:load(s) for s in train_ids+val_ids};x=combine(data,train_ids,'x');xv=combine(data,val_ids,'x')
 selected={};main_records=[];selection_records=[];map_genes={}
 for panel in CFG['panels']:
  ix=panel_indices(panel);y=combine(data,train_ids,'y')[:,ix];yv=combine(data,val_ids,'y')[:,ix]
  folder=OUT/'preprocessing/primary'/panel;xs,ys=scaler_pair(x,y,folder)
  xt,yt,xvt,yvt=xs.transform(x),ys.transform(y),xs.transform(xv),ys.transform(yv)
  variance=np.log1p(combine(data,train_ids,'raw')[:,ix].astype(float)).var(0)
  order=np.lexsort((np.array(PANELS[panel]),-variance));positive=order[variance[order]>0]
  chosen=[int(order[0]),int(order[len(order)//2]),int(positive[-1])]
  map_genes[panel]=list(dict.fromkeys(PANELS[panel][i] for i in chosen))
  candidates=[]
  for a in CFG['ridge_alphas']:
   r=fit_run(f'tune_{panel}_B1_a{a:g}','selection',panel,'B1',0,xt,yt,xvt,yvt,train_ids,val_ids,{'alpha':a},folder)
   candidates.append(r);selection_records.append(r['run_id'])
  best=min(candidates,key=lambda r:(-(r['validation_median_pearson'] if np.isfinite(r['validation_median_pearson']) else -np.inf),r['validation_mse'],r['params']['alpha']))
  ridge={'alpha':best['params']['alpha']}
  candidates=[]
  for i,params in enumerate(CFG['mlp_candidates']):
   r=fit_run(f'tune_{panel}_B2_c{i}','selection',panel,'B2',0,xt,yt,xvt,yvt,train_ids,val_ids,params,folder)
   candidates.append(r);selection_records.append(r['run_id'])
  best=min(candidates,key=lambda r:r['validation_mse']);mlp=best['params']
  selected[panel]={'ridge':ridge,'mlp':mlp,'selected_mlp_candidate':best['run_id']}
  for kind in ['B0','B1','B2']:
   seeds=CFG['seeds'] if kind=='B2' else [0]
   for seed in seeds:
    params=mlp if kind=='B2' else ridge if kind=='B1' else {}
    r=fit_run(f'primary_{panel}_{kind}_s{seed}','primary',panel,kind,seed,xt,yt,xvt,yvt,train_ids,val_ids,params,folder)
    main_records.append(r['run_id'])
  epochs=[json.loads((OUT/'logs/runs'/f'primary_{panel}_B2_s{s}.json').read_text())['best_epoch'] for s in CFG['seeds']]
  selected[panel]['lopo_fixed_epochs']=int(np.floor(np.median(epochs)+.5))
  if panel=='S2-50':
   for nc in ['NC1','NC2']:
    for kind in ['B1','B2']:
     for seed in CFG['seeds']:
      params=ridge if kind=='B1' else mlp
      r=fit_run(f'{nc}_{panel}_{kind}_s{seed}',nc,panel,kind,seed,xt,yt,xvt,yvt,train_ids,val_ids,params,folder,negative=nc)
      main_records.append(r['run_id'])
 lock={'timestamp':datetime.now(timezone.utc).isoformat(),'source_commit':commit(),'test_loaded':False,
       'selected':selected,'panels':PANELS,'map_genes':map_genes,'primary_and_control_runs':main_records,
       'selection_runs':selection_records,'training_slides':train_ids,'validation_slides':val_ids,
       'c0_lock_sha256':sha(ROOT/'reports/C0/C0_SPEC_LOCK.json'),'protocol_sha256':sha(OUT/'configs/protocol.json'),
       'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in OUT.glob('*.py')},
       'scaler_sha256':{str(p.relative_to(ROOT)):sha(p) for p in (OUT/'preprocessing').rglob('*.npz')}}
 dump(OUT/'SELECTION_LOCK.json',lock)
 print('SELECTION FROZEN; primary test data have not been loaded.',flush=True)

def load_scaler(path):
 s=Scaler()
 with np.load(path) as f:s.mean=f['mean'];s.scale=f['scale']
 return s

def load_model(record):
 path=ROOT/record['model_path']
 if sha(path)!=record['model_sha256']:raise ValueError('Model checksum mismatch')
 if record['model']=='B2':
  model=mlp_model(record['input_dim'],record['output_dim'],record['params']);model.load_state_dict(torch.load(path,map_location='cpu',weights_only=True));model.eval();return model
 with np.load(path) as f:
  if record['model']=='B0':return f['mean'].copy()
  model=Ridge();model.coef_=f['coef'].copy();model.intercept_=f['intercept'].copy();return model

def evaluate_record(record,sample,data,selection_sha):
 result_path=OUT/'logs/evaluation'/f'{record["run_id"]}_{sample}.json'
 if result_path.exists():raise RuntimeError('Evaluation already exists')
 folder=ROOT/record['preprocessing'];xs=load_scaler(folder/'feature_scaler.npz');ys=load_scaler(folder/'target_scaler.npz')
 ix=panel_indices(record['panel']);y=data['y'][:,ix]
 with PeakRSS() as memory:
  start=time.perf_counter();model=load_model(record)
  pred=ys.inverse(predict(model,xs.transform(data['x']),record['model']))
  seconds=time.perf_counter()-start
  m=gene_metrics(y,pred);summary=summarize(m)
  co,truth_corr,pred_corr=coexpression(y,pred)
  pred_path=OUT/'predictions'/f'{record["run_id"]}_{sample}.npz'
  np.savez_compressed(pred_path,prediction=pred.astype(np.float32),barcodes=data['barcodes'],genes=np.array(PANELS[record['panel']]))
 rows=pd.DataFrame({'gene':PANELS[record['panel']],**m})
 for k,v in {'run_id':record['run_id'],'phase':record['phase'],'panel':record['panel'],'model':record['model'],'seed':record['seed'],'sample_id':sample,'patient_id':PATIENTS[sample]}.items():rows[k]=v
 gene_path=OUT/'logs/gene_metrics'/f'{record["run_id"]}_{sample}.csv';gene_path.parent.mkdir(exist_ok=True);rows.to_csv(gene_path,index=False)
 result={**{k:record[k] for k in ['run_id','phase','panel','model','seed']},'sample_id':sample,'patient_id':PATIENTS[sample],
         'n_locations':len(y),**summary,**co,'inference_seconds':seconds,'inference_peak_rss_mb':memory.peak/1024**2,
         'selection_lock_sha256':selection_sha,'prediction_path':str(pred_path.relative_to(ROOT)),
         'prediction_sha256':sha(pred_path),'gene_metrics_path':str(gene_path.relative_to(ROOT))}
 dump(result_path,result)
 registry({**record,**summary,'event':'evaluation','train_slides':'|'.join(record['train_slides']),
           'evaluation_slides':sample,'inference_seconds':seconds,'peak_rss_mb':memory.peak/1024**2})
 print(f'EVAL {record["run_id"]} {sample}: median PCC={summary["median_pearson"]:.4f}, MSE={summary["mean_mse"]:.4f}',flush=True)
 return result

def evaluate():
 lock=json.loads((OUT/'SELECTION_LOCK.json').read_text());lock_sha=sha(OUT/'SELECTION_LOCK.json')
 for p,digest in lock['source_sha256'].items():
  if sha(ROOT/p)!=digest:raise ValueError('Model source changed after selection')
 test_ids=SPLITS['splits']['test'];data={s:load(s) for s in test_ids}
 dump(OUT/'logs/TEST_ACCESS.json',{'timestamp':datetime.now(timezone.utc).isoformat(),'selection_lock_sha256':lock_sha,'test_ids':test_ids})
 for run_id in lock['primary_and_control_runs']:
  record=json.loads((OUT/'logs/runs'/f'{run_id}.json').read_text())
  for sample in test_ids:evaluate_record(record,sample,data[sample],lock_sha)

def lopo():
 lock=json.loads((OUT/'SELECTION_LOCK.json').read_text());lock_sha=sha(OUT/'SELECTION_LOCK.json')
 data={s:load(s) for s in IDS}
 for heldout in IDS:
  train_ids=[s for s in IDS if s!=heldout]
  x=combine(data,train_ids,'x')
  for panel in CFG['panels']:
   ix=panel_indices(panel);y=combine(data,train_ids,'y')[:,ix]
   folder=OUT/'lopo'/heldout/panel;xs,ys=scaler_pair(x,y,folder);xt,yt=xs.transform(x),ys.transform(y)
   params=lock['selected'][panel]
   dump(folder/'config.json',{'label':'Conditional four-patient robustness sensitivity analysis — not unbiased validation',
        'train_slides':train_ids,'evaluation_slide':heldout,'gene_panel':PANELS[panel],
        'fixed_epochs':params['lopo_fixed_epochs'],'primary_parameters':params,
        'feature_scaler_sha256':sha(folder/'feature_scaler.npz'),'target_scaler_sha256':sha(folder/'target_scaler.npz')})
   results=[]
   for kind in ['B0','B1','B2']:
    seeds=CFG['seeds'] if kind=='B2' else [0]
    for seed in seeds:
     p=params['ridge'] if kind=='B1' else params['mlp'] if kind=='B2' else {}
     record=fit_run(f'lopo_{heldout}_{panel}_{kind}_s{seed}','lopo',panel,kind,seed,xt,yt,None,None,train_ids,[heldout],p,folder,fixed_epochs=params['lopo_fixed_epochs'] if kind=='B2' else None)
     results.append(evaluate_record(record,heldout,data[heldout],lock_sha))
   pd.DataFrame(results).to_csv(folder/'metrics.csv',index=False)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['select','evaluate','lopo']);args=parser.parse_args()
 for name in ['models','predictions','logs/runs','logs/evaluation','logs/gene_metrics']:(OUT/name).mkdir(parents=True,exist_ok=True)
 torch.set_num_threads(CFG['threads']);torch.use_deterministic_algorithms(True)
 with threadpool_limits(CFG['threads']):globals()[args.stage]()
if __name__=='__main__':main()
