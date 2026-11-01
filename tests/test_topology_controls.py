from graphpathology.graph_controls import from_edges,randomize,audit,canonical_hash
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

def fixture():
 xy=np.array([(x*100,y*100) for x in range(12) for y in range(12)],float)
 edges=cKDTree(xy).query_pairs(150,output_type='ndarray')
 # Keep all cardinal edges and one diagonal orientation, making sparse local alternatives.
 delta=xy[edges[:,0]]-xy[edges[:,1]]
 edges=edges[(np.linalg.norm(delta,axis=1)<125)|(delta[:,0]*delta[:,1]>0)]
 return from_edges(edges,len(xy)),xy

def test_vertex_permutation():
 a,xy=fixture();b,_=randomize(a,xy,'A4',4);r=audit(b,a,xy,'A4');assert r['changed_edge_fraction']>.8

def test_degree_swaps():
 a,xy=fixture();b,_=randomize(a,xy,'A5',2);r=audit(b,a,xy,'A5');assert r['degree_preserved'] and r['changed_edge_fraction']>.8

def test_constrained_swaps():
 a,xy=fixture();b,_=randomize(a,xy,'A6',3);r=audit(b,a,xy,'A6');assert r['distance_distribution_preserved'] and r['degree_preserved']
 assert r['changed_edge_fraction']>0

def test_reproducible():
 a,xy=fixture()
 for control in ['A4','A5','A6']:
  b,_=randomize(a,xy,control,9);c,_=randomize(a,xy,control,9);assert canonical_hash(b)==canonical_hash(c)

def test_labels_not_input():
 import inspect
 assert list(inspect.signature(randomize).parameters)==['a','xy','control','seed']
