import copy
import numpy as np
import pytest
import torch
from graphpathology.extension_models import VARIANTS,ExtensionModel,graph_spec,model_inputs,SCALE_UM

def fixture(variant):
 rng=np.random.default_rng(404);xy=rng.normal(size=(37,2))*80;x=rng.normal(size=(37,512)).astype(np.float32)
 s=graph_spec(xy,'synthetic',0,variant);return xy,x,s

@pytest.mark.parametrize('variant',list(VARIANTS))
def test_extension_forward_backward_serialization_offline(variant,tmp_path):
 torch.set_num_threads(1);torch.manual_seed(0);xy,x,s=fixture(variant);saved=xy.copy();p=s['p'].clone()
 model=ExtensionModel(variant);values=torch.tensor(model_inputs(x,s,variant));result=model(values,s)
 assert result.shape==(37,50) and torch.isfinite(result).all()
 result.square().mean().backward();assert any(q.grad is not None for q in model.parameters())
 assert all(q.grad is None or torch.isfinite(q.grad).all() for q in model.parameters())
 model.eval();path=tmp_path/'state.pt';torch.save(model.state_dict(),path);other=ExtensionModel(variant);other.load_state_dict(torch.load(path,weights_only=True));other.eval()
 with torch.inference_mode():torch.testing.assert_close(model(values,s),other(values,s),rtol=0,atol=0)
 np.testing.assert_array_equal(xy,saved);torch.testing.assert_close(p,s['p'],rtol=0,atol=0)
 if model.kind=='gps':
  with torch.inference_mode():_,_,weights=model(values,s,True)
  assert max(w.shape[-1] for w in weights)<=32

@pytest.mark.parametrize('variant',['gps_full','gps_local','gps_no_relative','distance_fixed','egnn_latent','distance_disabled'])
def test_conditional_e2_and_node_permutation(variant):
 torch.set_num_threads(1);torch.manual_seed(7);xy,x,s=fixture(variant);model=ExtensionModel(variant).double().eval();z=torch.tensor(x,dtype=torch.float64)
 with torch.inference_mode():base,p,_=model(z,s,True)
 angle=.731;rotation=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
 for matrix in [rotation,rotation@np.diag([-1.,1.]),np.eye(2)]:
  offset=np.array([391.2,-172.7]);q=xy@matrix.T+offset;t=graph_spec(q,'synthetic',0,variant)
  with torch.inference_mode():actual,latent,_=model(z,t,True)
  torch.testing.assert_close(base,actual,rtol=1e-8,atol=1e-8)
  old=p.numpy()*SCALE_UM+s['center_um'].numpy();new=latent.numpy()*SCALE_UM+t['center_um'].numpy()
  np.testing.assert_allclose(new,old@matrix.T+offset,rtol=1e-8,atol=1e-8)
 perm=np.random.default_rng(4).permutation(len(x));t=graph_spec(xy[perm],'synthetic',0,variant)
 with torch.inference_mode():actual=model(z[perm],t)
 torch.testing.assert_close(actual,base[perm],rtol=1e-8,atol=1e-8)


def test_degree_and_length_control_is_verified():
 _,_,s=fixture('gps_degree_length');assert s['info']['degree_preserved'] and s['info']['global_50um_bins_preserved']
 assert s['info']['successful_swaps']>0
