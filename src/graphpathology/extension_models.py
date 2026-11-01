"""post_v1_exploratory: bounded GraphGPS-inspired and 2D EGNN models.

Original observed coordinates are never updated. EGNN updates a private latent channel.
No model accepts expression labels; measured labels enter separate residual interpolation.
"""
import hashlib
import numpy as np
import torch
from torch import nn
from scipy import sparse
from .graph import build_graph, aggregate

SCALE_UM=250.
VARIANTS={
 'gps_full':('gps','spatial','global',True,False),
 'gps_local':('gps','spatial','local',True,False),
 'gps_no_relative':('gps','spatial','global',False,False),
 'gps_self':('gps','self','self',True,False),
 'gps_shuffled':('gps','shuffled','global',True,False),
 'gps_degree_length':('gps','degree_length','global',True,False),
 'distance_fixed':('egnn','spatial','none',True,False),
 'egnn_latent':('egnn','spatial','none',True,True),
 'distance_disabled':('egnn','spatial','none',False,False),
 'egnn_self':('egnn','self','none',True,True),
 'egnn_shuffled':('egnn','shuffled','none',True,True),
 'egnn_degree_length':('egnn','degree_length','none',True,True),
 'matched_local':('local','spatial','none',False,False),
 'matched_mean':('MEAN','spatial','none',False,False),
 'matched_gcn':('GCN','spatial','none',False,False),
 'matched_sage':('SAGE','spatial','none',False,False),
}


def stable_seed(sample,seed):
 return int.from_bytes(hashlib.sha256(f'extension|{sample}|{seed}'.encode()).digest()[:8],'little')


def landmarks(xy,k=32):
 """Geometry-only farthest-point landmarks; synthetic equivariance tests avoid ties."""
 xy=np.asarray(xy,float);n=len(xy);minimum=np.full(n,np.inf)
 first=int(np.argmax(np.sum((xy-xy.mean(0))**2,1)));order=[];j=first
 for _ in range(min(k,n)):
  order.append(j);minimum=np.minimum(minimum,np.sum((xy-xy[j])**2,1));minimum[order]=-np.inf;j=int(np.argmax(minimum))
 return np.array(order,dtype=np.int64)


def length_rewire(a,xy,seed):
 """Double-edge swaps preserve every degree and the global 50-um length-bin counts."""
 upper=sparse.triu(a,k=1).tocoo();edges=[tuple(e) for e in zip(upper.row.tolist(),upper.col.tolist())];existing=set(edges);rng=np.random.default_rng(seed)
 def bin_of(e):return min(4,int(np.linalg.norm(xy[e[0]]-xy[e[1]])/50))
 target=max(1,len(edges)//4);swaps=0;attempts=0
 for _ in range(max(1,20*len(edges))):
  if swaps>=target or len(edges)<2:break
  attempts+=1;i,j=rng.choice(len(edges),2,replace=False);u,v=edges[i];w,z=edges[j]
  if len({u,v,w,z})<4:continue
  proposed=[tuple(sorted((u,z))),tuple(sorted((w,v)))]
  if proposed[0] in existing or proposed[1] in existing:continue
  if any(np.linalg.norm(xy[e[0]]-xy[e[1]])>250 for e in proposed):continue
  if sorted([bin_of(edges[i]),bin_of(edges[j])])!=sorted(map(bin_of,proposed)):continue
  existing.remove(edges[i]);existing.remove(edges[j]);existing.update(proposed);edges[i],edges[j]=proposed;swaps+=1
 pairs=np.array(edges,dtype=int).reshape(-1,2);row=np.r_[pairs[:,0],pairs[:,1]];col=np.r_[pairs[:,1],pairs[:,0]]
 b=sparse.csr_matrix((np.ones(len(row)),(row,col)),shape=a.shape);b.sort_indices()
 np.testing.assert_array_equal(np.asarray(a.sum(1)),np.asarray(b.sum(1)))
 return b,{'successful_swaps':swaps,'target_swaps':target,'attempts':attempts,'degree_preserved':True,'global_50um_bins_preserved':True}


def graph_spec(xy,sample,seed,variant):
 kind,control,attention,geometry,latent=VARIANTS[variant];physical=np.asarray(xy,dtype=np.float64).copy();p=physical.copy()
 a=build_graph(p,{'kind':'delaunay','max_edge_um':250});info={'control':control}
 if control=='shuffled':
  perm=np.random.default_rng(stable_seed(sample,seed)).permutation(len(p));a=a[perm,:][:,perm].tocsr();p=p[perm]
 elif control=='degree_length':a,extra=length_rewire(a,p,stable_seed(sample,seed));info.update(extra)
 elif control=='self':a=sparse.csr_matrix(a.shape)
 a.sort_indices();row=np.repeat(np.arange(len(p)),np.diff(a.indptr));col=a.indices.copy();degree=np.maximum(np.diff(a.indptr),1)
 # Padding repeats self and is explicit; no N-by-N attention mask is allocated.
 local_links=np.tile(np.arange(len(p))[:,None],(1,8))
 for i in range(len(p)):
  js=a.indices[a.indptr[i]:a.indptr[i+1]];js=js[np.argsort(np.linalg.norm(p[js]-p[i],axis=1),kind='stable')][:7]
  local_links[i,1:len(js)+1]=js
 center=p.mean(0);spec={'row':torch.tensor(row,dtype=torch.long),'col':torch.tensor(col,dtype=torch.long),'degree':torch.tensor(degree,dtype=torch.float32),
  'p':torch.tensor((p-center)/SCALE_UM,dtype=torch.float64),'center_um':torch.tensor(center,dtype=torch.float64),
  'landmarks':torch.tensor(landmarks(p),dtype=torch.long),'local_links':torch.tensor(local_links,dtype=torch.long),
  'physical_xy':physical,'adjacency':a,'info':info}
 return spec


def mean_messages(messages,row,degree,n):
 out=messages.new_zeros((n,messages.shape[-1]));out.index_add_(0,row,messages);return out/degree.to(messages.dtype)[:,None]


def radial(value,maximum,width):
 centers=torch.linspace(0,maximum,8,dtype=value.dtype,device=value.device)
 return torch.exp(-.5*((value[...,None]-centers)/width)**2)


class GPSBlock(nn.Module):
 def __init__(self,d=64,heads=4,dropout=.1):
  super().__init__();self.d=d;self.heads=heads;self.norm=nn.LayerNorm(d);self.local=nn.Linear(d,d);self.gate=nn.Linear(8,d)
  self.local_out=nn.Linear(d,d);self.q=nn.Linear(d,d);self.k=nn.Linear(d,d);self.v=nn.Linear(d,d);self.out=nn.Linear(d,d);self.bias=nn.Linear(8,heads)
  self.ff=nn.Sequential(nn.LayerNorm(d),nn.Linear(d,2*d),nn.GELU(),nn.Dropout(dropout),nn.Linear(2*d,d));self.drop=nn.Dropout(dropout)
 def forward(self,h,spec,mode,geometry):
  u=self.norm(h);row=spec['row'];col=spec['col'];p=spec['p'].to(h.dtype);n=len(h);heads=self.heads;dh=self.d//heads
  distance=torch.linalg.vector_norm(p[row]-p[col],dim=1);rbf=radial(torch.log1p(distance),1.,.15)
  gate=torch.sigmoid(self.gate(rbf)) if geometry else torch.ones((len(row),self.d),dtype=h.dtype,device=h.device)
  local=self.local_out(mean_messages(self.local(u)[col]*gate,row,spec['degree'],n)+self.local(u))
  q=self.q(u).reshape(n,heads,dh)
  if mode=='global':
   ids=spec['landmarks'];k=self.k(u[ids]).reshape(len(ids),heads,dh);v=self.v(u[ids]).reshape(len(ids),heads,dh)
   score=torch.einsum('nhd,lhd->nhl',q,k)/np.sqrt(dh);dist=torch.linalg.vector_norm(p[:,None]-p[ids][None],dim=-1)
   bias=self.bias(radial(torch.log1p(dist),4.,.5)).permute(0,2,1)
   if geometry:score=score+bias
   weights=score.softmax(-1);att=torch.einsum('nhl,lhd->nhd',weights,v).reshape(n,self.d)
  else:
   ids=spec['local_links'] if mode=='local' else torch.arange(n,device=h.device)[:,None]
   k=self.k(u)[ids].reshape(n,ids.shape[1],heads,dh);v=self.v(u)[ids].reshape(n,ids.shape[1],heads,dh)
   score=torch.einsum('nhd,nlhd->nhl',q,k)/np.sqrt(dh);dist=torch.linalg.vector_norm(p[:,None]-p[ids],dim=-1)
   if geometry:score=score+self.bias(radial(torch.log1p(dist),4.,.5)).permute(0,2,1)
   weights=score.softmax(-1);att=torch.einsum('nhl,nlhd->nhd',weights,v).reshape(n,self.d)
  h=h+self.drop(local+self.out(att));return h+self.drop(self.ff(h)),weights


class EGBlock(nn.Module):
 def __init__(self,d=64,latent=False,dropout=.1):
  super().__init__();self.latent=latent;self.norm=nn.LayerNorm(d)
  self.edge=nn.Sequential(nn.Linear(2*d+8,d),nn.SiLU(),nn.Linear(d,d),nn.SiLU())
  self.node=nn.Sequential(nn.Linear(2*d,d),nn.SiLU(),nn.Dropout(dropout),nn.Linear(d,d))
  if latent:
   self.coordinate=nn.Linear(d,1,bias=False);nn.init.normal_(self.coordinate.weight,std=.001)
 def forward(self,h,p,spec,geometry):
  u=self.norm(h);r=spec['row'];c=spec['col'];delta=p[r]-p[c];d2=(delta**2).sum(-1)
  features=radial(d2,4.,.5) if geometry else h.new_zeros((len(r),8))
  m=self.edge(torch.cat([u[r],u[c],features],-1));agg=mean_messages(m,r,spec['degree'],len(h));out=h+self.node(torch.cat([u,agg],-1))
  if self.latent:
   shift=h.new_zeros((len(h),2));shift.index_add_(0,r,.1*torch.tanh(self.coordinate(m))*delta)
   p=p+shift/spec['degree'].to(h.dtype)[:,None]
  return out,p


class ExtensionModel(nn.Module):
 def __init__(self,variant,input_dim=512,genes=50):
  super().__init__();self.variant=variant;self.kind,self.control,self.attention,self.geometry,self.latent=VARIANTS[variant]
  if self.kind in ['local','MEAN','GCN','SAGE']:
   d=2*input_dim if self.kind=='SAGE' else input_dim;hidden=67 if self.kind=='SAGE' else 128
   self.head=nn.Sequential(nn.Linear(d,hidden),nn.GELU(),nn.Dropout(.1),nn.Linear(hidden,genes))
  else:
   self.input=nn.Linear(input_dim,64)
   self.blocks=nn.ModuleList([GPSBlock() for _ in range(2)]) if self.kind=='gps' else nn.ModuleList([EGBlock(latent=self.latent) for _ in range(2)])
   self.readout=nn.Sequential(nn.LayerNorm(64),nn.Linear(64,genes))
 def forward(self,x,spec,diagnostics=False):
  p=spec['p'].to(x.dtype).clone();weights=[]
  if self.kind in ['local','MEAN','GCN','SAGE']:
   result=self.head(x)
  else:
   h=self.input(x)
   for block in self.blocks:
    if self.kind=='gps':h,w=block(h,spec,self.attention,self.geometry);weights.append(w)
    else:h,p=block(h,p,spec,self.geometry)
   result=self.readout(h)
  return (result,p,weights) if diagnostics else result


def model_inputs(x,spec,variant):
 kind=VARIANTS[variant][0]
 return np.ascontiguousarray(aggregate(spec['adjacency'],x,kind)) if kind in ['MEAN','GCN','SAGE'] else np.ascontiguousarray(x)
