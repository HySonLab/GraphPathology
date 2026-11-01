from graphpathology.graph import build_graph,operator,aggregate,shuffle_graph,hidden_width,parameter_count
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy import sparse

@pytest.mark.parametrize('cfg',[{'kind':'knn','k':6,'max_edge_um':250},{'kind':'radius','radius_um':150},{'kind':'radius','radius_um':250},{'kind':'delaunay','max_edge_um':250}])
def test_graph_physical_bounds_and_determinism(cfg):
 xy=np.stack(np.meshgrid(np.arange(6)*100,np.arange(5)*100),axis=-1).reshape(-1,2)
 a=build_graph(xy,cfg);assert not (a-a.T).nnz and not a.diagonal().any()
 assert (a-build_graph(xy,cfg)).nnz==0
 i,j=a.nonzero();assert np.max(np.linalg.norm(xy[i]-xy[j],axis=1))<=cfg.get('max_edge_um',cfg.get('radius_um'))+1e-8

def test_mean_and_symmetric_normalization():
 a=sparse.csr_matrix([[0,1,0],[1,0,0],[0,0,0]])
 np.testing.assert_allclose(operator(a,'MEAN').sum(1),np.ones((3,1)))
 assert (operator(a,'GCN')-operator(a,'GCN').T).nnz==0
 x=np.array([[1.,2.],[3.,4.],[7.,8.]],dtype=np.float32)
 np.testing.assert_allclose(aggregate(a,x,'MEAN'),[[2,3],[2,3],[7,8]])
 np.testing.assert_allclose(aggregate(a,x,'SAGE')[2],[7,8,0,0])

def test_gcn_precomputation_matches_message_passing():
 a=sparse.csr_matrix([[0,1,1],[1,0,0],[1,0,0]])
 rng=np.random.default_rng(2);x=rng.normal(size=(3,4));w=rng.normal(size=(4,2))
 np.testing.assert_allclose((operator(a,'GCN')@x)@w,operator(a,'GCN')@(x@w),atol=1e-12)

def test_shuffling_preserves_degrees_but_changes_edges():
 xy=np.stack(np.meshgrid(np.arange(5)*100,np.arange(5)*100),axis=-1).reshape(-1,2);a=build_graph(xy,{'kind':'radius','radius_um':150})
 b,p=shuffle_graph(a,'sample',0)
 assert (a-b).nnz>0;np.testing.assert_array_equal(np.sort(a.sum(1).A.ravel()),np.sort(b.sum(1).A.ravel()))
 np.testing.assert_array_equal(p,shuffle_graph(a,'sample',0)[1])

def test_capacity_match():
 for g,h in [(1,128),(50,128),(250,256),(280,256)]:
  actual=parameter_count(1024,hidden_width(512,g,h,'SAGE'),g);expected=parameter_count(512,h,g)
  assert abs(actual-expected)/expected<.01

def test_duplicate_and_nonfinite_coordinates_rejected():
 for xy in [np.zeros((4,2)),np.full((4,2),np.nan)]:
  with pytest.raises(ValueError):build_graph(xy,{'kind':'radius','radius_um':100})


def test_aggregation_is_local_and_cannot_cross_disconnected_regions():
 xy=np.array([[0.,0.],[100.,0.],[1000.,0.],[1100.,0.]])
 a=build_graph(xy,{'kind':'radius','radius_um':150})
 x=np.zeros((4,2),dtype=np.float32);changed=x.copy();changed[0]=1
 for model in ['MEAN','GCN','SAGE']:
  before=aggregate(a,x,model);after=aggregate(a,changed,model)
  np.testing.assert_array_equal(before[2:],after[2:])
  assert np.any(before[:2]!=after[:2])
