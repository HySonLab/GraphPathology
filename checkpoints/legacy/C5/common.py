"""C5 read-only parent imports and append-only run provenance."""
import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','4');os.environ.setdefault('MKL_NUM_THREADS','4')
import sys,json,csv,hashlib,subprocess,importlib.util,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C5'
sys.path.insert(0,str(ROOT/'C2'))
spec=importlib.util.spec_from_file_location('c5_parent',ROOT/'C2/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
core=base.core;PANELS={k:v for k,v in base.PANELS.items() if k!='S1'};IDS=base.IDS;SPLITS=base.SPLITS;ix=base.ix;sha=base.sha;dump=base.dump;commit=base.commit
ENCODERS=['resnet18_reference','pathology_encoder_1'];REGIMES=['random','block','capture_window']
def now():return datetime.now(timezone.utc).isoformat()
def rid(prefix,obj):return 'C5_'+prefix+'_'+hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()[:18]
def registry(event,run_id,record):
 p=OUT/'logs/registry'/f'{event}__{run_id}.json';assert not p.exists()
 r={**record,'event':event,'run_id':run_id,'timestamp':now(),'checkpoint':'C5','source_commit':commit()};dump(p,r)
 fields=next(csv.reader((ROOT/'results/experiment_registry.csv').open()));row={**r,'threads':4,'cpu':'Intel Core i7-14700'}
 for k in ['train_slides','evaluation_slides']:
  if isinstance(row.get(k),list):row[k]='|'.join(row[k])
 with (ROOT/'results/experiment_registry.csv').open('a',newline='') as f:csv.DictWriter(f,fieldnames=fields,extrasaction='ignore').writerow(row)
def cfg():return json.loads((OUT/'configs/protocol.json').read_text())
def masks(sample,fold=False):
 rows=json.loads((OUT/'configs/masks.json').read_text());return [r for r in rows if r['sample_id']==sample and (not fold or r['nominal_coverage'] in [.1,.2,.4])]
def mask_array(r):
 with np.load(ROOT/r['path']) as f:return f['observed'].copy()
def labels(sample,panel):
 with np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as f:return f['normalized_targets'][:,ix(panel)].copy()
def check_lock():
 p=OUT/'C5_SPEC_LOCK.json';r=json.loads(p.read_text());assert sha(p)==(OUT/'C5_SPEC_LOCK.sha256').read_text().split()[0]
 for name,digest in r['source_sha256'].items():assert sha(ROOT/name)==digest,name
 return r
