"""Shared C4 contracts; previous checkpoint code is imported read-only."""
import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','4');os.environ.setdefault('MKL_NUM_THREADS','4')
import sys,json,csv,hashlib,subprocess,importlib.util
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C4';A=OUT/'C4A'
sys.path.insert(0,str(ROOT/'C2'))
spec=importlib.util.spec_from_file_location('c4_base',ROOT/'C2/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
core=base.core;PANELS=base.PANELS;IDS=base.IDS;SPLITS=base.SPLITS;ix=base.ix;sha=base.sha;dump=base.dump;commit=base.commit

def now():return datetime.now(timezone.utc).isoformat()
def registry(event,run_id,record):
 r={**record,'event':event,'run_id':run_id,'checkpoint':'C4A','timestamp':now(),'source_commit':commit()};path=A/'logs/registry'/f'{event}__{run_id}.json'
 assert not path.exists(),path;dump(path,r)
 fields=next(csv.reader((ROOT/'results/experiment_registry.csv').open()))
 row={**r,'model':r.get('model',r.get('encoder','')),'cpu':'Intel Core i7-14700','threads':r.get('threads',4)}
 for k in ['train_slides','evaluation_slides']:
  if isinstance(row.get(k),list):row[k]='|'.join(row[k])
 with (ROOT/'results/experiment_registry.csv').open('a',newline='') as f:csv.DictWriter(f,fieldnames=fields,extrasaction='ignore').writerow(row)
def scaler(path):
 s=core.Scaler()
 with np.load(path) as f:s.mean=f['mean'];s.scale=f['scale']
 return s
