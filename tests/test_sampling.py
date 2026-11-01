from graphpathology.sampling import choose,policy_view,target_indices
import sys
from pathlib import Path
import numpy as np

def test_full_unrevealed_expression_never_changes_selection():
 rng=np.random.default_rng(0);raw={'xy':rng.normal(size=(90,2)),'projection':rng.normal(size=(90,8))};eligible=np.arange(20,90)
 for strategy in ['A0','A1','A2','A3']:
  order,score=choose(raw,eligible,36,strategy)
  for y in [np.full((90,280),np.nan),rng.normal(size=(90,280))]:
   alt,values=choose({**raw,'expression':y},eligible,36,strategy);assert np.array_equal(order,alt);assert np.array_equal(score,values,equal_nan=True)
  assert len(np.unique(order))==36 and np.isin(order,eligible).all()

def test_stress_region_always_fixed_and_unobserved():
 unavailable=np.zeros(100,bool);unavailable[:40]=True
 for k in [5,10,20,40]:
  selected=np.arange(40,40+k);h,obs=target_indices(100,selected,unavailable,'contiguous');assert np.array_equal(h,np.arange(40));assert not obs[h].any();assert obs.sum()==k

def test_permutation_controls_only_change_sampler_view():
 rng=np.random.default_rng(1);raw={'xy':rng.normal(size=(80,2)),'projection':rng.normal(size=(80,8))};before={k:v.copy() for k,v in raw.items()}
 for strategy,control in [('A1','NC3'),('A2','NC2'),('A3','NC2')]:choose(raw,np.arange(80),30,strategy,control=control)
 assert all(np.array_equal(raw[k],v) for k,v in before.items());assert 'expression' not in policy_view({**raw,'expression':np.ones((80,280))})
