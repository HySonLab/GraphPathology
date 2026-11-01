from graphpathology.external import predict,moments,np
import sys
from pathlib import Path

def test_hidden_invariance_and_affine_identity():
 rng=np.random.default_rng(91);xy=rng.normal(size=(40,2));im=rng.normal(size=(40,7));y=rng.normal(size=(40,7));obs=np.arange(0,40,3);h=np.setdiff1d(np.arange(40),obs)
 obj={'source':moments([y[:10]]),'destination':moments([3+y[10:20]*np.arange(1,8)])}
 for regime in ['V1','V2']:
  for control in ['','NC1','NC2','NC3']:
   p,_=predict(xy,im,y,obs,h,regime,obj,control,2);z=y.copy();z[h]=np.nan;q,_=predict(xy,im,z,obs,h,regime,obj,control,2)
   assert all(np.array_equal(p[m],q[m]) for m in p)
   np.testing.assert_allclose(p['M1'],p['M2'],rtol=1e-12,atol=1e-12)

def test_patient_weighting():
 a=np.zeros((100,2));b=np.full((2,2),2.)
 obj=moments([a,b]);np.testing.assert_allclose(obj['center'],1);np.testing.assert_allclose(obj['scale'],1)
