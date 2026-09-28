"""Label-independent, slide-local spatial graph operators for C2."""
import hashlib
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree,Delaunay


def build_graph(xy,config):
 xy=np.asarray(xy,dtype=np.float64);n=len(xy)
 if xy.ndim!=2 or xy.shape[1]!=2 or n<2 or not np.isfinite(xy).all():raise ValueError('Expected finite Nx2 coordinates')
 if len(np.unique(xy,axis=0))!=n:raise ValueError('Duplicate spatial coordinates')
 tree=cKDTree(xy);kind=config['kind']
 if kind=='radius':
  pairs=tree.query_pairs(config['radius_um'],output_type='ndarray')
 elif kind=='knn':
  k=min(config['k'],n-1);dist,_=tree.query(xy,k=k+1);edges=[]
  for i in range(n):
   candidates=np.array(tree.query_ball_point(xy[i],min(float(dist[i,-1])+1e-7,config['max_edge_um'])),dtype=int)
   candidates=candidates[candidates!=i]
   lengths=np.linalg.norm(xy[candidates]-xy[i],axis=1)
   order=np.lexsort((candidates,lengths))[:k]
   edges.extend((i,int(j)) for j in candidates[order])
  pairs=np.asarray(edges,dtype=int).reshape(-1,2)
 elif kind=='delaunay':
  triangles=Delaunay(xy).simplices
  pairs=np.concatenate([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[0,2]]])
  pairs=pairs[np.linalg.norm(xy[pairs[:,0]]-xy[pairs[:,1]],axis=1)<=config['max_edge_um']]
 else:raise ValueError('Unknown graph kind')
 if len(pairs):
  pairs=np.unique(np.sort(pairs,axis=1),axis=0)
  rows=np.r_[pairs[:,0],pairs[:,1]];cols=np.r_[pairs[:,1],pairs[:,0]]
 else:rows=cols=np.array([],dtype=int)
 a=sparse.csr_matrix((np.ones(len(rows),dtype=np.float64),(rows,cols)),shape=(n,n));a.sort_indices()
 if a.diagonal().any() or (a-a.T).nnz:raise AssertionError('Expected undirected graph without self loops')
 return a


def operator(a,model):
 if model in ['MEAN','GCN']:
  b=a+sparse.eye(a.shape[0],format='csr');degree=np.asarray(b.sum(1)).ravel()
  if model=='MEAN':return sparse.diags(1/degree)@b
  inv=sparse.diags(1/np.sqrt(degree));return inv@b@inv
 if model=='SAGE':
  degree=np.asarray(a.sum(1)).ravel();return sparse.diags(1/np.maximum(degree,1))@a
 raise ValueError('Unknown model')


def aggregate(a,x,model):
 if len(x)!=a.shape[0]:raise ValueError('Graph/feature row mismatch')
 neighbor=np.asarray(operator(a,model)@x,dtype=np.float32)
 return np.concatenate([x,neighbor],axis=1) if model=='SAGE' else neighbor


def shuffle_graph(a,sample,seed):
 digest=hashlib.sha256(f'C2|{sample}|{seed}'.encode()).digest()
 p=np.random.default_rng(int.from_bytes(digest[:8],'little')).permutation(a.shape[0])
 return a[p,:][:,p].tocsr(),p


def parameter_count(d,h,g):return (d+1)*h+(h+1)*g

def hidden_width(d,g,local_h,model):
 if model!='SAGE':return local_h
 target=parameter_count(d,local_h,g)
 h=max(1,int(round((target-g)/(2*d+1+g))))
 return h
