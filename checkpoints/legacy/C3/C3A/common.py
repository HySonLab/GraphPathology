import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','4');os.environ.setdefault('MKL_NUM_THREADS','4')
import sys,json,csv,hashlib,importlib.util,subprocess
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'C2'))
spec=importlib.util.spec_from_file_location('c2_readonly',ROOT/'C2/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
core=base.core;PANELS=base.PANELS;IDS=base.IDS;SPLITS=base.SPLITS;load=base.load;ix=base.ix;scaler=base.scaler;preprocessing=base.preprocessing;sha=base.sha;dump=base.dump;commit=base.commit
CFG=json.loads((OUT/'configs/protocol.json').read_text());SELECTED=json.loads((ROOT/'C2/SELECTION_LOCK.json').read_text())['selected']

def now():return datetime.now(timezone.utc).isoformat()
def registry(event,run_id,record):
 record={**record,'event':event,'run_id':run_id,'checkpoint':'C3A','source_commit':commit(),'timestamp':now()}
 fields=next(csv.reader((ROOT/'results/experiment_registry.csv').open()))
 row={**record,'model':record.get('control',record.get('model','')),'seed':record.get('training_seed',record.get('graph_randomization_seed','')),'threads':4,'cpu':'Intel Core i7-14700'}
 for k in ['train_slides','evaluation_slides']:
  if isinstance(row.get(k),list):row[k]='|'.join(row[k])
 # Event-level JSON is the lossless schema extension; global CSV remains append-only.
 path=OUT/'logs/registry'/f'{event}__{run_id}.json'
 if path.exists():raise ValueError(f'Duplicate event {event} {run_id}')
 dump(path,record)
 with (ROOT/'results/experiment_registry.csv').open('a',newline='') as f:csv.DictWriter(f,fieldnames=fields,extrasaction='ignore').writerow(row)

def check_lock():
 path=OUT/'C3A_SPEC_LOCK.json';lock=json.loads(path.read_text())
 assert sha(path)==(OUT/'C3A_SPEC_LOCK.sha256').read_text().split()[0]
 for p,h in lock['source_sha256'].items():assert sha(ROOT/p)==h,p
 assert sha(ROOT/'C2/C2_SPEC_LOCK.json')==lock['parent_C2_lock_sha256']
 return lock
