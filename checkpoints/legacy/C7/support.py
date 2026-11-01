"""Read-only frozen C5 method contract and isolated C7 provenance."""
import os,sys,json,csv,hashlib,time,subprocess
from pathlib import Path
from datetime import datetime,timezone
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','4');os.environ.setdefault('MKL_NUM_THREADS','4')
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C7';A=OUT/'C7A';B=OUT/'C7B';CACHE=OUT/'cache/her2'
sys.path.insert(0,str(ROOT/'C5'));import common as parent;import engine as c5
import numpy as np,pandas as pd,torch
from scipy.spatial import cKDTree
from graph import build_graph,operator
core=parent.core;sha=parent.sha;dump=parent.dump;ENC='resnet18_reference';SEL=json.loads((ROOT/'C5/SELECTION_LOCK.json').read_text());ORIGINAL=parent.PANELS
SAMPLES=[s+'1' for s in 'ABCDEFGH'];SCENES={};TRUTH={};WEIGHTS={}
def now():return datetime.now(timezone.utc).isoformat()
def cfg():return json.loads((OUT/'configs/protocol.json').read_text())
def panels():return cfg()['panels']
def count(n,b):return int(np.floor(n*b+.5))
def commit():return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
def rid(prefix,obj):return 'C7_'+prefix+'_'+hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()[:18]
def registry(event,run_id,record):
 p=B/'logs/registry'/f'{event}__{run_id}.json';assert not p.exists();r={**record,'event':event,'run_id':run_id,'timestamp':now(),'checkpoint':'C7','source_commit':commit()};dump(p,r)
 fields=next(csv.reader((ROOT/'results/experiment_registry.csv').open()));row={**r,'threads':4,'evaluation_slides':r.get('sample_id','')}
 with (ROOT/'results/experiment_registry.csv').open('a',newline='') as f:csv.DictWriter(f,fieldnames=fields,extrasaction='ignore').writerow(row)
def check_lock():
 p=A/'C7A_SPEC_LOCK.json';r=json.loads(p.read_text());assert sha(p)==(A/'C7A_SPEC_LOCK.sha256').read_text().split()[0];assert r['external_evaluation_authorized']
 for p,h in r['source_sha256'].items():assert sha(ROOT/p)==h,p
 return r

def scene(sample,panel):
 key=(sample,panel)
 if key not in SCENES:
  with np.load(CACHE/(sample+'_features.npz')) as f:x=f['embeddings'].copy();bc=f['barcodes'].copy();xy=f['xy_um'].copy()
  with np.load(CACHE/(sample+'_'+panel+'_image.npz')) as f:im=f['prediction'].copy();assert np.array_equal(f['barcodes'],bc)
  with np.load(ROOT/'C5/models'/f'PCA_primary_{ENC}.npz') as pc:z=((x-pc['mean'])/pc['scale']).astype(np.float32);z=((z-pc['pca_mean'])@pc['components'].T/np.sqrt(pc['variance'])).astype(np.float32)
  if sample not in c5.GRAPHS:c5.GRAPHS[sample]=operator(build_graph(xy,{'name':'knn6_cap250','kind':'knn','k':6,'max_edge_um':250}),'MEAN').tocsr()
  SCENES[key]={'xy':xy,'barcodes':bc,'image':im,'projection':z,'sample':sample,'feature_cache_hash':sha(CACHE/(sample+'_features.npz'))}
 return SCENES[key]
def predict(s,a,h,panel,regime,model):
 sp=SEL['simple'][panel]['random' if regime=='V1' else 'block'];lam=sp['lambda'][ENC]
 if model=='E0':return c5.predict(s,a,h,'B1',None)
 if model=='E1':return c5.predict(s,a,h,sp['best_family'],sp['best_param'])
 if model=='E2':return c5.predict(s,a,h,'B3',None,sp,lam)
 if panel not in WEIGHTS:
  p=ROOT/'C5/models'/f'B4_primary_{ENC}_{panel}_a{SEL["B4_alpha"][panel][ENC]}.npz';assert sha(p)==json.loads(p.with_suffix('.json').read_text())['checkpoint_hash']
  with np.load(p) as f:WEIGHTS[panel]=f['weights'][[ORIGINAL[panel].index(g) for g in panels()[panel]]].copy()
 return c5.predict(s,a,h,'B4',None,sp,lam,WEIGHTS[panel])
def truth(sample,panel):
 check_lock();assert (B/'MASK_LOCK.json').exists()
 if sample not in TRUTH:
  with np.load(CACHE/(sample+'_targets.npz')) as f:TRUTH[sample]={k:f[k].copy() for k in f.files}
 t=TRUTH[sample];idx=[list(t['genes']).index(g) for g in panels()[panel]];return t['normalized_targets'][:,idx].copy()
