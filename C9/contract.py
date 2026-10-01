"""Prospective C9 contract: frozen parent models, patient-scoped labels, common units."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(k,'4')
import sys,json,csv,hashlib,time,subprocess
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C9';A=OUT/'C9A';B=OUT/'C9B';CACHE=OUT/'cache';RAW=OUT/'raw/extracted'
sys.path.insert(0,str(ROOT/'C7'));import support as parent
core=parent.core;sha=parent.sha;dump=parent.dump;ENC=parent.ENC
SAMPLES=['1142243F','1160920F','CID4290','CID4465','CID44971','CID4535']
CAL=['1142243F','CID4465','CID4535'];EVAL=['1160920F','CID4290','CID44971']
def now():return datetime.now(timezone.utc).isoformat()
def commit():return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
def cfg():return json.loads((OUT/'configs/protocol.json').read_text())
def genes():return cfg()['genes']
def check_lock():
 p=A/'C9A_SPEC_LOCK.json';assert sha(p)==p.with_suffix('.sha256').read_text().split()[0];r=json.loads(p.read_text())
 for name,h in r['source_sha256'].items():assert sha(ROOT/name)==h,name
 return r

def registry(event,ident,rec):
 p=B/'logs/registry'/(ident+'.json');assert not p.exists();r={**rec,'run_id':ident,'event':event,'checkpoint':'C9','source_commit':commit(),'timestamp':now()};dump(p,r)
 fields=next(csv.reader((ROOT/'results/experiment_registry.csv').open()));row={**r,'threads':4,'evaluation_slides':r.get('sample_id',''),'train_slides':'|'.join(CAL)}
 with (ROOT/'results/experiment_registry.csv').open('a',newline='') as f:csv.DictWriter(f,fieldnames=fields,extrasaction='ignore').writerow(row)
def read_targets(sample,purpose):
 check_lock();assert sample in (CAL if purpose=='calibration' else EVAL)
 if purpose=='evaluation':assert (B/'CALIBRATION_LOCK.json').exists() and (B/'MASK_LOCK.json').exists()
 from scipy.io import mmread
 import gzip
 d=RAW/'filtered_count_matrices'/sample
 # Source file layout is fixed by qualification.
 files=cfg()['count_files'][sample]
 with (ROOT/files['features']).open() as f:features=[l.strip() for l in f]
 with (ROOT/files['barcodes']).open() as f:bc=np.array([l.strip() for l in f])
 with (ROOT/files['matrix']).open('rb') as f:mat=mmread(f).tocsr()
 with np.load(CACHE/(sample+'_metadata.npz')) as f:kept=f['barcodes'].copy()
 idxbc=pd.Index(bc).get_indexer(kept);assert (idxbc>=0).all()
 raw=mat[[features.index(g) for g in genes()]].toarray().T[idxbc].astype(float)
 assert np.isfinite(raw).all() and (raw>=0).all() and np.array_equal(raw,np.rint(raw))
 norm=json.loads((ROOT/'reports/C0/C0_NORMALIZATION.json').read_text());idx=[norm['genes'].index(g) for g in genes()]
 y=(np.log1p(raw)-np.array(norm['mean'])[idx])/np.array(norm['scale'])[idx]
 return y

def moments(arrays):
 mean=np.mean([a.mean(0) for a in arrays],axis=0);var=np.mean([np.mean((a-mean)**2,axis=0) for a in arrays],axis=0);sd=np.sqrt(var)
 return {'center':mean.tolist(),'scale':np.where(sd<1e-8,1,sd).tolist(),'constant_fallback':(sd<1e-8).tolist()}
def predict(xy,image,y,observed,hidden,regime,scaler,control='',control_seed=0):
 # Slice immediately: no hidden target can enter any input feature.
 values=np.asarray(y)[observed].copy();positions=xy[observed].copy();im=image[hidden].copy();rng=np.random.default_rng(control_seed)
 if control=='NC1':values=values[rng.permutation(len(values))]
 if control=='NC2':positions=positions[rng.permutation(len(positions))]
 if control=='NC3':im=image[rng.permutation(len(image))][hidden]
 dist,idx=cKDTree(positions).query(xy[hidden],k=min(16,len(observed)))
 if idx.ndim==1:idx=idx[:,None];dist=dist[:,None]
 dist=np.round(dist,6);weights=1/np.maximum(dist,1e-6)**2 if regime=='V1' else np.ones_like(dist);weights/=weights.sum(1,keepdims=True)
 at=np.array(scaler['destination']['center']);bt=np.array(scaler['destination']['scale']);asrc=np.array(scaler['source']['center']);bs=np.array(scaler['source']['scale'])
 st=np.einsum('nk,nkg->ng',weights,values[idx],optimize=True)
 z=np.einsum('nk,nkg->ng',weights,((values-at)/bt)[idx],optimize=True);harm=at+bt*z
 calibrated_image=at+bt*((im-asrc)/bs)
 return {'M0':im,'M1':st,'M2':harm,'M3':.2*calibrated_image+.8*harm},dist[:,0]
def score(y,p):
 met=core.gene_metrics(y,p);r={}
 for key in ['pearson','spearman','mse','mae']:
  for stat in ['mean','median']:r[stat+'_'+key]=core.finite_stat(met[key],stat)
 for t in [0,.1,.2,.3]:r['fraction_pearson_gt_'+str(t)]=float(np.mean(np.isfinite(met['pearson'])&(met['pearson']>t)))
 r['n_valid_pearson']=int(np.isfinite(met['pearson']).sum());return r,met
