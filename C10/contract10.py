"""C10 frozen inference contract; no destination labels used for fitting."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(k,'4')
import sys,json,csv,hashlib,time,subprocess,importlib.util
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,h5py
from scipy.sparse import csc_matrix
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C10';A=OUT/'C10A';B=OUT/'C10B';CACHE=OUT/'cache';RAW=OUT/'raw'
spec=importlib.util.spec_from_file_location('frozen_c9_contract',ROOT/'C9/contract.py');c9=importlib.util.module_from_spec(spec);spec.loader.exec_module(c9)
parent=c9.parent;core=c9.core;sha=c9.sha;dump=c9.dump;ENC=c9.ENC
import policy as frozen_policy
SAMPLES=[r['patient_id'] for r in json.loads((A/'logs/SOURCE_SAMPLE_SELECTION.json').read_text())['selected']]
def now():return datetime.now(timezone.utc).isoformat()
def commit():return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
def cfg():return json.loads((OUT/'configs/protocol.json').read_text())
def genes():return cfg()['genes']
def check_lock():
 p=A/'SPEC_LOCK.json';assert sha(p)==p.with_suffix('.sha256').read_text().split()[0];r=json.loads(p.read_text())
 for name,h in r['source_sha256'].items():assert sha(ROOT/name)==h,name
 return r

def registry(event,ident,rec):
 p=B/'logs/registry'/(ident+'.json');assert not p.exists();r={**rec,'run_id':ident,'event':event,'checkpoint':'C10','source_commit':commit(),'timestamp':now()};dump(p,r)
 fields=next(csv.reader((ROOT/'results/experiment_registry.csv').open()));row={**r,'threads':4,'evaluation_slides':r.get('patient_id',''),'train_slides':'none; frozen C9 image calibration only'}
 with (ROOT/'results/experiment_registry.csv').open('a',newline='') as f:csv.DictWriter(f,fieldnames=fields,extrasaction='ignore').writerow(row)
def targets(patient):
 check_lock();assert patient in SAMPLES and (B/'MASK_LOCK.json').exists()
 assert (B/'logs/EVALUATION_LABEL_ACCESS.json').exists()
 p=ROOT/cfg()['inputs'][patient]['counts']
 with h5py.File(p) as f:
  group=f['matrix'];gg=group['features/name'][:].astype(str);bc=group['barcodes'][:].astype(str);mat=csc_matrix((group['data'][:],group['indices'][:],group['indptr'][:]),shape=tuple(group['shape'][:]))
 with np.load(CACHE/(patient+'_metadata.npz')) as f:kept=f['barcodes'].copy()
 rows=pd.Index(bc).get_indexer(kept);assert (rows>=0).all();idx=[list(gg).index(g) for g in genes()];raw=mat[idx][:,rows].toarray().T.astype(float)
 assert np.isfinite(raw).all() and (raw>=0).all() and np.array_equal(raw,np.rint(raw))
 norm=json.loads((ROOT/'reports/C0/C0_NORMALIZATION.json').read_text());jj=[norm['genes'].index(g) for g in genes()]
 return (np.log1p(raw)-np.array(norm['mean'])[jj])/np.array(norm['scale'])[jj]
def scaler():
 obj=json.loads((ROOT/'C9/C9A/T1_SCALER.json').read_text());idx=[obj['genes'].index(g) for g in genes()]
 for domain in ['source','destination']:
  for key in ['center','scale','constant_fallback']:obj[domain][key]=[obj[domain][key][i] for i in idx]
 obj['genes']=genes();return obj

def predict(xy,image,y,observed,hidden,regime,obj,control='',control_seed=0):
 preds,dist=c9.predict(xy,image,y,observed,hidden,regime,obj,control,control_seed)
 assert np.allclose(preds['M1'],preds['M2'],rtol=1e-12,atol=1e-12)
 return {m:preds[m] for m in ['M0','M1','M3']},dist
score=c9.score
