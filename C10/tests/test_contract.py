import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from contract10 import *
def test_frozen_c9_inference_and_hidden_invariance():
 rng=np.random.default_rng(43);xy=rng.normal(size=(80,2));y=rng.normal(size=(80,50));im=rng.normal(size=(80,50));obs=np.arange(0,80,4);hidden=np.setdiff1d(np.arange(80),obs);obj=scaler()
 for regime in ['V1','V2']:
  for control in ['','NC1','NC2','NC3']:
   p,_=predict(xy,im,y,obs,hidden,regime,obj,control,2);frozen,_=c9.predict(xy,im,y,obs,hidden,regime,obj,control,2)
   for mode in ['nan','random']:
    corrupt=y.copy();corrupt[hidden]=np.nan if mode=='nan' else rng.normal(size=corrupt[hidden].shape);q,_=predict(xy,im,corrupt,obs,hidden,regime,obj,control,2)
    for m in p:assert np.array_equal(p[m],q[m]) and np.array_equal(p[m],frozen[m])
def test_no_c10_calibration():
 assert cfg()['calibration_patients']==[]
 assert not set(SAMPLES)&set(scaler()['calibration_patients'])
