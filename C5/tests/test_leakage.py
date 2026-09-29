import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from engine import *

def toy():
 rng=np.random.default_rng(17);n=12;xy=np.column_stack([np.arange(n)*100,np.zeros(n)]);obs=np.zeros(n,bool);obs[[0,5,11]]=True;y=rng.normal(size=(n,3)).astype(np.float32);s={'xy':xy,'sample':'toy','image':rng.normal(size=(n,3)).astype(np.float32),'projection':rng.normal(size=(n,8)).astype(np.float32)};GRAPHS['toy']=sparse.eye(n,format='csr');return y,obs,s

def test_whole_hidden_vectors_do_not_enter_any_model():
 y,obs,s=toy();h=np.flatnonzero(~obs);a=observed_payload(y,obs,s['xy']);spec={'best_family':'B2a','best_param':1};w=np.random.default_rng(1).normal(size=(3,35))
 for corruption in [np.nan,923.75]:
  changed=y.copy();changed[~obs]=corruption;b=observed_payload(changed,obs,s['xy']);assert np.array_equal(a['values'],b['values']);assert len(b['values'])==obs.sum()
  for model,param in [('B0',None),('B1',None),('B2a',1),('B2b',2),('B3',None),('B4',None)]:
   p=predict(s,a,h,model,param,spec,.5,w,np.zeros(3));q=predict(s,b,h,model,param,spec,.5,w,np.zeros(3));np.testing.assert_array_equal(p,q)

def test_indicator_and_expression_zero_fill():
 y,obs,s=toy();o,z,first,second,_,_=molecular(s,observed_payload(y,obs,s['xy']));np.testing.assert_array_equal(o,obs.astype(np.float32));assert np.all(z[~obs]==0);np.testing.assert_array_equal(z[obs],y[obs]);assert np.all(first[~obs]==0) and np.all(second[~obs]==0)

def test_controls_keep_hidden_locations_out_of_anchor_bank():
 y,obs,s=toy();a=observed_payload(y,obs,s['xy'])
 for name in ['NC1','NC2','NC3']:
  _,b=control_inputs(s,a,name,3);assert set(b['indices'])==set(np.flatnonzero(obs));assert np.isfinite(b['values']).all()
