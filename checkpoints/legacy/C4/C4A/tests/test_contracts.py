import sys,importlib.util
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from common import core,A
import json

def test_equal_grids_and_seed_budget():
 c=json.loads((A/'configs/protocol.json').read_text());assert len(c['encoders'])==2 and len(c['mlp_candidates'])==4 and len(c['ridge_alphas'])==5 and c['training_seeds']==[0,1,2]

def test_scaling_uses_fit_sample_only():
 x=np.array([[1.,3.],[3.,5.]]);s=core.Scaler().fit(x);before=s.mean.copy();s.transform(np.array([[1e6,-1e6]]));assert np.array_equal(s.mean,before)

def test_permutation_preserves_marginals_breaks_pairing():
 x=np.arange(100).reshape(20,5).astype(np.float32);y=x[:,:2].copy()
 for mode in ['NC1','NC2']:
  a,b,p=core.permute_training(x,y,mode,2);assert sorted(p)==list(range(20)) and not np.array_equal(p,np.arange(20));assert np.array_equal(np.sort(a,axis=0),np.sort(x,axis=0));assert np.array_equal(np.sort(b,axis=0),np.sort(y,axis=0))

def test_gain_sign_conventions():
 # Perfect rank/linear association does not imply calibrated expression values.
 y=np.arange(1,6,dtype=float)[:,None];m=core.gene_metrics(y,2*y)
 assert np.isclose(m['pearson'][0],1.) and m['mse'][0]>0
