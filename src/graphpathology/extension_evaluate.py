"""Matched-mask extension evaluation with per-model/patient checkpoints and exact leakage checks."""
import gzip,json,time
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy import sparse
from .assets import sha
from .core import Scaler,mlp_model,predict,gene_metrics,summarize
from .graph import build_graph,aggregate,shuffle_graph
from .journal_data import original_masks,WU_EVAL
from .journal_methods import identical_target_design
from .journal_learned import hybrids
from .policy import count
from .extension_models import ExtensionModel,graph_spec,model_inputs
from .extension_train import dataset,metadata

BINS=[0,150,300,600,1200,2400,np.inf]

def model_inventory(root):
 folder=root/'experiments/journal_revision/transformer_egnn/fits';items=[]
 for p in sorted(folder.glob('*.json')):
  r=json.loads(p.read_text())
  if r['status']=='COMPLETE':items.append({'id':f"new_{r['variant']}_s{r['seed']}",'origin':'new','variant':r['variant'],'seed':r['seed'],'checkpoint':r['best_checkpoint_relative'],'checkpoint_sha256':r['best_checkpoint_sha256'],'fit_record':str(p.relative_to(root))})
 old=json.loads((root/'experiments/journal_revision/results/hest_learned/FIT_MANIFEST.json').read_text())
 for r in old:items.append({'id':'frozen_'+r['run_id'],'origin':'frozen','variant':r['model']+'_'+r['control'],'seed':r['seed'],'checkpoint':r['checkpoint_relative_path'],'checkpoint_sha256':r['checkpoint_sha256'],'old_model':r['model'],'old_control':r['control'],'old_config':r['config']})
 return items

@lru_cache(maxsize=4)
def load_data(root,patient):
 if patient=='TENX99':return dataset(root,[patient])[patient]
 p=root/'data/journal_revision/transformer_egnn/wu'/(patient+'.npz');record=json.loads(p.with_suffix('.json').read_text());assert sha(p)==record['cache_sha256']
 with np.load(p) as f:d={k:f[k].copy() for k in f.files}
 with np.load(root/'checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz') as f:d['x']=((d['embeddings']-f['mean'])/f['scale']).astype(np.float32)
 return d

def prediction(root,d,patient,item):
 cp=root/item['checkpoint'];assert sha(cp)==item['checkpoint_sha256'];state=torch.load(cp,map_location='cpu',weights_only=True)
 if item['origin']=='new':
  spec=graph_spec(d['xy_um'],patient,item['seed'],item['variant']);model=ExtensionModel(item['variant']);model.load_state_dict(state['model']);model.eval()
  with torch.inference_mode():y=model(torch.tensor(model_inputs(d['x'],spec,item['variant'])),spec).numpy().astype(float)
 else:
  a=build_graph(d['xy_um'],{'kind':'delaunay','max_edge_um':250});kind=item['old_model'];control=item['old_control'];x=d['x']
  if kind!='local':
   if control=='self':a=sparse.eye(len(x),format='csr') if kind=='SAGE' else sparse.csr_matrix(a.shape)
   elif control=='shuffled':a,_=shuffle_graph(a,patient,item['seed'])
   x=aggregate(a,x,kind)
  model=mlp_model(x.shape[1],50,item['old_config']);model.load_state_dict(state['state_dict']);model.eval();y=predict(model,x,'B2').astype(float)
 if patient!='TENX99':
  cal=json.loads((root/'checkpoints/legacy/C9/C9A/T1_SCALER.json').read_text());np.testing.assert_array_equal(d['genes'],cal['genes'])
  y=np.array(cal['destination']['center'])+np.array(cal['destination']['scale'])*(y-np.array(cal['source']['center']))/np.array(cal['source']['scale'])
 return y

def masks(root,d,patient):
 n=len(d['targets']);records=[];cache=root/'data/journal_revision/transformer_egnn/masks'/patient;cache.mkdir(parents=True,exist_ok=True)
 if patient=='TENX99':
  for budget in [.1,.2,.4]:
   for seed in range(3):
    p=root/'checkpoints/legacy/reports/C0/C0_MASKS'/f'TENX99_random_q{budget:.2f}_r{seed}.npz'
    with np.load(p) as f:np.testing.assert_array_equal(f['barcodes'],d['barcodes']);obs=np.flatnonzero(f['observed']);hidden=np.flatnonzero(~f['observed'])
    records.append({'design':'historical_V1','condition':'distributed','budget':budget,'mask_seed':seed,'geometry':0,'observed':obs,'hidden':hidden,'mask_sha256':sha(p)})
 else:
  for r in original_masks(root,'wu',patient):
   if r['seed']>2 or (r['regime']=='V2' and r['geometry']!=1):continue
   with np.load(r['absolute_mask_path']) as f:
    np.testing.assert_array_equal(f['barcodes'],d['barcodes']);np.testing.assert_array_equal(f['xy_um'],d['xy_um']);order=f['order'].copy();forbidden=f['forbidden'].copy()
   for budget in [.1,.2,.4]:
    obs=order[:count(n,budget)];hidden=np.setdiff1d(np.arange(n),obs) if r['regime']=='V1' else np.flatnonzero(forbidden)
    records.append({'design':'historical_'+r['regime'],'condition':'historical','budget':budget,'mask_seed':r['seed'],'geometry':r['geometry'],'observed':obs,'hidden':hidden,'mask_sha256':r['mask_sha256']})
 for seed in range(3):
  region,hidden,orders,holdout_seed=identical_target_design(d['xy_um'],patient,seed)
  if patient in WU_EVAL:
   old=root/'data/journal_revision/processed/journal_v2_masks/wu'/patient/f'central_seed{seed}.npz'
   with np.load(old) as f:
    np.testing.assert_array_equal(hidden,f['hidden']);np.testing.assert_array_equal(region,f['region'])
    for name,order in orders.items():np.testing.assert_array_equal(order,f[name+'_order'])
  p=cache/f'B2_seed{seed}.npz'
  if not p.exists():np.savez_compressed(p,region=region,hidden=hidden,distributed_order=orders['distributed'],region_exclusion_order=orders['region_exclusion'],barcodes=d['barcodes'],xy_um=d['xy_um'])
  for budget in [.1,.2,.4]:
   for condition,order in orders.items():records.append({'design':'B2_identical_H','condition':condition,'budget':budget,'mask_seed':seed,'geometry':1,'observed':order[:count(n,budget)],'hidden':hidden,'mask_sha256':sha(p)})
 manifest=cache/'MANIFEST.json'
 clean=[{k:v for k,v in r.items() if k not in ['observed','hidden']}|{'n_observed':len(r['observed']),'n_scored':len(r['hidden'])} for r in records]
 if not manifest.exists():manifest.write_text(json.dumps(clean,indent=2)+'\n')
 return records

def score_set(d,image,mask,item,patient):
 y=d['targets'].astype(float);obs=mask['observed'];hidden=mask['hidden'];base={k:v for k,v in mask.items() if k not in ['observed','hidden']}
 base.update(patient=patient,cohort='Development' if patient=='TENX99' else 'Wu',model_id=item['id'],origin=item['origin'],variant=item['variant'],training_seed=item['seed'])
 outputs=hybrids(d['xy_um'],y,image,obs,hidden);unknown=np.setdiff1d(np.arange(len(y)),obs);checks=0
 for mode in ['nan','random']:
  corrupt=y.copy();corrupt[unknown]=np.nan if mode=='nan' else np.random.default_rng(31).normal(size=(len(unknown),50))
  alt=hybrids(d['xy_um'],corrupt,image,obs,hidden)
  for route in ['fixed_20pct_image','exploratory_residual_anchor']:np.testing.assert_array_equal(outputs[route],alt[route]);checks+=1
 from scipy.spatial import cKDTree
 distances=np.round(cKDTree(d['xy_um'][obs]).query(d['xy_um'][hidden],k=1)[0],6)
 panels=[];genes=[];bins=[]
 for route in ['image','fixed_20pct_image','exploratory_residual_anchor']:
  pred=outputs[route];m=gene_metrics(y[hidden],pred);attr={**base,'route':route}
  panels.append({**attr,'n_observed':len(obs),'n_scored':len(hidden),'mean_nearest_anchor_um':float(distances.mean()),**summarize(m)})
  for j,g in enumerate(d['genes']):genes.append({**attr,'gene':str(g),**{k:float(m[k][j]) for k in ['pearson','spearman','mse','mae']}})
  for low,high in zip(BINS[:-1],BINS[1:]):
   selected=(distances>=low)&(distances<high)
   if not selected.any():continue
   err=y[hidden][selected]-pred[selected];bins.append({**attr,'low_um':low,'high_um':high,'n_locations':int(selected.sum()),'mean_mse':float(np.mean(err**2)),'mean_mae':float(np.mean(abs(err)))})
 return panels,genes,bins,checks

def _evaluate_one(root,patient,item):
 out=root/'experiments/journal_revision/transformer_egnn/evaluation'/patient/item['id'];out.mkdir(parents=True,exist_ok=True);done=out/'COMPLETE.json'
 if done.exists():
  r=json.loads(done.read_text())
  for name,digest in r['output_hashes'].items():assert sha(out/name)==digest
  return r
 d=load_data(root,patient);ms=masks(root,d,patient);cache=root/'data/journal_revision/transformer_egnn/predictions'/patient;cache.mkdir(parents=True,exist_ok=True)
 p=cache/(item['id']+'.npy');meta=p.with_suffix('.json')
 input_path=root/'data/journal_revision/hest_training/processed/journal_v2_hest/TENX99.npz' if patient=='TENX99' else root/'data/journal_revision/transformer_egnn/wu'/(patient+'.npz')
 input_digest=sha(input_path)
 if p.exists() and meta.exists():
  r=json.loads(meta.read_text());assert r['checkpoint_sha256']==item['checkpoint_sha256'] and r['input_cache_sha256']==input_digest and sha(p)==r['prediction_sha256'];image=np.load(p)
 else:
  image=prediction(root,d,patient,item);np.save(p,image);meta.write_text(json.dumps({'checkpoint_sha256':item['checkpoint_sha256'],'input_cache_sha256':input_digest,'prediction_sha256':sha(p),'calibration_sha256':sha(root/'checkpoints/legacy/C9/C9A/T1_SCALER.json') if patient!='TENX99' else None},indent=2)+'\n')
 panels=[];genes=[];bins=[];checks=0;start=time.monotonic()
 for m in ms:
  a,b,c,k=score_set(d,image,m,item,patient);panels.extend(a);genes.extend(b);bins.extend(c);checks+=k
 for name,rows in [('PANEL.csv.gz',panels),('GENES.csv.gz',genes),('BINS.csv.gz',bins)]:pd.DataFrame(rows).to_csv(out/name,index=False)
 result=metadata(root,status='COMPLETE',patient=patient,model=item,input_cache_sha256=input_digest,prediction_sha256=sha(p),n_mask_sets=len(ms),n_route_rows=len(panels),exact_label_corruption_checks=checks,
  elapsed_seconds=time.monotonic()-start,output_hashes={name:sha(out/name) for name in ['PANEL.csv.gz','GENES.csv.gz','BINS.csv.gz']},calibration='none' if patient=='TENX99' else 'frozen original C9 T1; no refit')
 done.write_text(json.dumps(result,indent=2)+'\n');print('EVALUATED',patient,item['id'],len(panels),flush=True);return result

def baselines(root,patient):
 from .journal_methods import predict_anchors
 out=root/'experiments/journal_revision/transformer_egnn/evaluation'/patient;out.mkdir(parents=True,exist_ok=True)
 if (out/'MOLECULAR.csv').exists():return
 d=load_data(root,patient);y=d['targets'].astype(float);rows=[]
 for m in masks(root,d,patient):
  for method,cfg in [('idw_k16_p2',{'kind':'idw','k':16,'p':2}),('uniform_k16',{'kind':'uniform','k':16})]:
   pred,dist=predict_anchors(d['xy_um'],y,m['observed'],m['hidden'],cfg);metrics=gene_metrics(y[m['hidden']],pred)
   rows.append({**{k:v for k,v in m.items() if k not in ['observed','hidden']},'patient':patient,'cohort':'Development' if patient=='TENX99' else 'Wu','model_id':method,'origin':'molecular','variant':method,'training_seed':-1,'route':'anchor_only','n_observed':len(m['observed']),'n_scored':len(m['hidden']),**summarize(metrics)})
 pd.DataFrame(rows).to_csv(out/'MOLECULAR.csv',index=False)

def run(root):
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True);items=model_inventory(root)
 for patient in ['TENX99']+WU_EVAL:
  baselines(root,patient)
  for item in items:evaluate_one(root,patient,item)
 return {'status':'COMPLETE' if sum(i['origin']=='new' for i in items)==48 else 'PARTIAL_PENDING_FITS','models':len(items),'patients':4}


def evaluate_one(root,patient,item):
 import fcntl
 folder=root/'data/journal_revision/transformer_egnn/locks';folder.mkdir(parents=True,exist_ok=True)
 with (folder/(patient+'_'+item['id']+'.lock')).open('w') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  return _evaluate_one(root,patient,item)
