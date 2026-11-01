"""Reproducible label-free graph controls with exact audited constraints."""
import hashlib
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.stats import ks_2samp
from numba import njit

BINS=np.array([0.,175.,250.,np.inf])

def seed_for(sample,geometry,control,seed):
 return int.from_bytes(hashlib.sha256(f'C3A|{sample}|{geometry}|{control}|{seed}'.encode()).digest()[:4],'little')

def edge_array(a):return np.column_stack(sparse.triu(a,k=1).nonzero()).astype(np.int64)
def from_edges(edges,n):
 i,j=edges.T
 a=sparse.csr_matrix((np.ones(len(i)*2),(np.r_[i,j],np.r_[j,i])),shape=(n,n));a.sort_indices();return a

def lengths(edges,xy):return np.linalg.norm(xy[edges[:,0]]-xy[edges[:,1]],axis=1)

def canonical_hash(a):
 a=a.tocsr();a.sort_indices()
 return hashlib.sha256(np.asarray(a.shape,dtype='<i8').tobytes()+a.indptr.astype('<i8').tobytes()+a.indices.astype('<i8').tobytes()).hexdigest()

@njit(cache=True)
def _bin(u,v,xy,bins):
 d=np.sqrt(((xy[u]-xy[v])**2).sum())
 for b in range(len(bins)-1):
  if d<bins[b+1]:return b
 return len(bins)-2

@njit(cache=True)
def _swaps(edges,degree,xy,bins,cptr,cidx,seed,distance,attempt_factor,swap_factor):
 np.random.seed(seed);n=len(degree);m=len(edges)
 adj=np.full((n,int(degree.max())),-1,np.int64);eid=adj.copy();fill=np.zeros(n,np.int64)
 for e in range(m):
  u,v=edges[e]
  adj[u,fill[u]]=v;eid[u,fill[u]]=e;fill[u]+=1
  adj[v,fill[v]]=u;eid[v,fill[v]]=e;fill[v]+=1
 accepted=0;attempts=0
 for step in range(attempt_factor*m):
  attempts+=1;e1=np.random.randint(m);u,v=edges[e1]
  if np.random.random()<.5:u,v=v,u
  if distance:
   count=cptr[u+1]-cptr[u]
   if count==0:continue
   c=cidx[cptr[u]+np.random.randint(count)]
   if degree[c]==0:continue
   slot=np.random.randint(degree[c]);d=adj[c,slot];e2=eid[c,slot]
  else:
   e2=np.random.randint(m);c,d=edges[e2]
   if np.random.random()<.5:c,d=d,c
  if u==c or u==d or v==c or v==d:continue
  if distance:
   b1=_bin(u,v,xy,bins);b2=_bin(c,d,xy,bins);b3=_bin(u,c,xy,bins);b4=_bin(v,d,xy,bins)
   if min(b1,b2)!=min(b3,b4) or max(b1,b2)!=max(b3,b4):continue
  duplicate=False
  for s in range(degree[u]):
   if adj[u,s]==c:duplicate=True
  for s in range(degree[v]):
   if adj[v,s]==d:duplicate=True
  if duplicate:continue
  for s in range(degree[u]):
   if adj[u,s]==v:adj[u,s]=c;eid[u,s]=e1;break
  for s in range(degree[v]):
   if adj[v,s]==u:adj[v,s]=d;eid[v,s]=e2;break
  for s in range(degree[c]):
   if adj[c,s]==d:adj[c,s]=u;eid[c,s]=e1;break
  for s in range(degree[d]):
   if adj[d,s]==c:adj[d,s]=v;eid[d,s]=e2;break
  edges[e1,0]=u;edges[e1,1]=c;edges[e2,0]=v;edges[e2,1]=d
  accepted+=1
  if accepted>=swap_factor*m:break
 return edges,accepted,attempts

def randomize(a,xy,control,seed):
 if control=='A4':
  p=np.random.default_rng(seed).permutation(a.shape[0]);b=a[p,:][:,p].tocsr()
  return b,{'algorithm':'vertex permutation','accepted_swaps':0,'attempts':0}
 edges=edge_array(a);degree=np.diff(a.indptr).astype(np.int64)
 if control=='A6':
  pairs=cKDTree(xy).query_pairs(275.,output_type='ndarray');c=from_edges(pairs,len(xy));ptr=c.indptr.astype(np.int64);idx=c.indices.astype(np.int64)
 else:ptr=np.zeros(len(xy)+1,np.int64);idx=np.zeros(0,np.int64)
 e,accepted,attempts=_swaps(edges.copy(),degree,np.asarray(xy,dtype=np.float64),BINS,ptr,idx,seed,control=='A6',1000 if control=='A6' else 100,5)
 return from_edges(e,len(xy)),{'algorithm':'degree-preserving double-edge swaps'+(' constrained by exact length-bin counts' if control=='A6' else ''),'accepted_swaps':int(accepted),'attempts':int(attempts)}

def audit(a,reference,xy,control):
 assert a.shape==reference.shape and not a.diagonal().any() and (a-a.T).nnz==0 and np.all(a.data==1)
 edges=edge_array(a);le=lengths(edges,xy);deg=np.diff(a.indptr);refdeg=np.diff(reference.indptr)
 counts=np.histogram(le,BINS)[0];refcounts=np.histogram(lengths(edge_array(reference),xy),BINS)[0]
 refle=lengths(edge_array(reference),xy)
 ks=float(ks_2samp(np.round(le,3),np.round(refle,3)).statistic)
 overlap=a.multiply(reference).nnz/reference.nnz
 degree_exact=np.array_equal(deg,refdeg);degree_distribution=np.array_equal(np.sort(deg),np.sort(refdeg));bins_equal=np.array_equal(counts,refcounts)
 assert degree_distribution
 if control in ['A5','A6']:assert degree_exact
 if control=='A6':assert bins_equal
 return {'num_nodes':a.shape[0],'num_edges':len(edges),'mean_degree':float(deg.mean()),'degree_min':int(deg.min()),'degree_q25':float(np.quantile(deg,.25)),'degree_median':float(np.median(deg)),'degree_q75':float(np.quantile(deg,.75)),'degree_max':int(deg.max()),'num_components':int(connected_components(a,directed=False,return_labels=False)),
 'mean_edge_length_um':float(le.mean()),'median_edge_length_um':float(np.median(le)),'edge_length_q05_um':float(np.quantile(le,.05)),'edge_length_q95_um':float(np.quantile(le,.95)),'max_edge_length_um':float(le.max()),'edge_length_bin_counts':counts.tolist(),'degree_preserved':bool(degree_exact),'degree_distribution_preserved':bool(degree_distribution),'distance_distribution_preserved':bool(bins_equal),'edge_length_ks_rounded_0_001_um':ks,'relative_mean_edge_length_change':float(le.mean()/refle.mean()-1),'changed_edge_fraction':float(1-overlap),'graph_hash':canonical_hash(a),'constraints_pass':True}
