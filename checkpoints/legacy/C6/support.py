"""C6 read-only C5 access, observation boundary, and append-only provenance."""
import sys,os,json,csv,hashlib,time,subprocess
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C6';A=OUT/'C6A'
sys.path.insert(0,str(ROOT/'C5'))
import common as parent
import engine as c5
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
core=parent.core;sha=parent.sha;dump=parent.dump;IDS=parent.IDS;ENCODERS=parent.ENCODERS;PANELS={p:parent.PANELS[p] for p in ['S2-50','S2-280']}
BUDGETS=[.05,.1,.2,.4];REGIMES=['open','contiguous'];SEL=json.loads((ROOT/'C5/SELECTION_LOCK.json').read_text());SCENES={};WEIGHTS={};TRUTH={}
def now():return datetime.now(timezone.utc).isoformat()
def commit():return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
def rid(prefix,obj):return 'C6_'+prefix+'_'+hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()[:18]
def config():return json.loads((A/'configs/protocol.json').read_text())
def count(n,b):return int(np.floor(n*b+.5))
def registry(event,run_id,record):
 path=A/'logs/registry'/f'{event}__{run_id}.json';assert not path.exists()
 r={**record,'event':event,'run_id':run_id,'timestamp':now(),'checkpoint':'C6A','source_commit':commit()};dump(path,r)
 fields=next(csv.reader((ROOT/'results/experiment_registry.csv').open()));row={**r,'threads':4,'cpu':'Intel Core i7-14700'}
 for key in ['train_slides','evaluation_slides']:
  if isinstance(row.get(key),list):row[key]='|'.join(row[key])
 with (ROOT/'results/experiment_registry.csv').open('a',newline='') as f:csv.DictWriter(f,fieldnames=fields,extrasaction='ignore').writerow(row)
def check_lock():
 p=A/'C6A_SPEC_LOCK.json';r=json.loads(p.read_text());assert sha(p)==(A/'C6A_SPEC_LOCK.sha256').read_text().split()[0]
 for name,h in r['source_sha256'].items():assert sha(ROOT/name)==h,name
 return r

def scene(cohort,enc,panel,sample):
 key=(cohort,enc,panel,sample)
 if key not in SCENES:
  # All C5 caches must exist: the parent code is used only in read-only mode.
  for p in [ROOT/'C5/models'/f'PCA_{cohort}_{enc}.npz',ROOT/'C5/cache'/('_'.join(key)+'_B1.npz')]:assert p.exists(),p
  SCENES[key]=c5.scene(*key)
 return SCENES[key]
def primary_model(enc):return config()['primary_model'][enc]
def model_spec(panel,regime,enc):
 sp=SEL['simple'][panel]['random' if regime=='open' else 'block'];return sp,sp['lambda'][enc]
def reconstruction(s,anchor,h,cohort,enc,panel,regime,role):
 sp,lam=model_spec(panel,regime,enc);model=primary_model(enc) if role=='primary' else sp['best_family'];w=None
 if model=='B4':
  key=(cohort,enc,panel)
  if key not in WEIGHTS:
   p=ROOT/'C5/models'/f'B4_{cohort}_{enc}_{panel}_a{SEL["B4_alpha"][panel][enc]}.npz';meta=json.loads(p.with_suffix('.json').read_text());assert sha(p)==meta['checkpoint_hash']
   with np.load(p) as f:WEIGHTS[key]=f['weights'].copy()
  w=WEIGHTS[key]
 return c5.predict(s,anchor,h,model,sp.get(model),sp,lam,w)
def expression(sample,panel,validation=False):
 if validation:assert sample==parent.SPLITS['splits']['validation'][0]
 else:assert (A/'TRAJECTORY_LOCK.json').exists(),'No held-out label access before every non-adaptive trajectory is frozen'
 key=(sample,panel)
 if key not in TRUTH:TRUTH[key]=parent.labels(sample,panel)
 return TRUTH[key]
def observed_only(sample,panel,selected,validation=False):
 y=expression(sample,panel,validation);obs=np.zeros(len(y),bool);obs[selected]=True
 return c5.observed_payload(y,obs,c5.image_data(ENCODERS[0],sample)['xy'])
def unavailable(sample,regime):
 n=len(c5.image_data(ENCODERS[0],sample)['xy'])
 if regime=='open':return np.zeros(n,bool),None
 r=next(r for r in parent.masks(sample) if r['regime']=='block' and r['nominal_coverage']==.6 and r['replicate']==0)
 return ~parent.mask_array(r),{'mask_id':r['mask_id'],'path':r['path'],'sha256':r['checksum']}
def target_indices(n,selected,forbidden,regime):
 obs=np.zeros(n,bool);obs[selected]=True
 return np.flatnonzero(forbidden if regime=='contiguous' else ~obs),obs

def provenance(cohort,enc,panel,regime,role,sample):
 sp,lam=model_spec(panel,regime,enc);model=primary_model(enc) if role=='primary' else sp['best_family'];path=ROOT/'C5/models'/f'B4_{cohort}_{enc}_{panel}_a{SEL["B4_alpha"][panel][enc]}.npz' if model=='B4' else ROOT/'C5/SELECTION_LOCK.json'
 return {'model':model,'config':{'interpolation':sp,'lambda':lam if model in ['B3','B4'] else None},'checkpoint_path':str(path.relative_to(ROOT)),'checkpoint_hash':sha(path),'image_prediction_hash':sha(ROOT/'C5/cache'/f'{cohort}_{enc}_{panel}_{sample}_B1.npz'),'scaler_hash':sha(ROOT/'C5/models'/f'PCA_{cohort}_{enc}.npz'),'feature_cache_hash':c5.image_data(enc,sample)['feature_cache_hash']}
