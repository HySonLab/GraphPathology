"""Observed-only C5 residual covariates, with explicit graph injection."""
import numpy as np
from scipy.spatial import cKDTree
from .graph import operator
GRAPHS={}
def graph(sample): return GRAPHS[sample]
def observed_payload(y,observed,xy):
 assert observed.dtype==bool and len(y)==len(observed)==len(xy)
 ids=np.flatnonzero(observed);values=np.array(y[ids],dtype=np.float32,copy=True);assert np.isfinite(values).all()
 return {'indices':ids,'xy':np.array(xy[ids],copy=True),'values':values}


def nearest(s,a):
 k=min(16,len(a['indices']));dist,index=cKDTree(a['xy']).query(s['xy'],k=k)
 if k==1:dist=dist[:,None];index=index[:,None]
 return np.round(dist,6),index


def interpolate(a,dist,index,family,param):
 if family=='B2a':return a['values'][index[:,:int(param)]].mean(1)
 w=1/np.maximum(dist,1e-6)**float(param);w/=w.sum(1,keepdims=True);return np.einsum('nk,nkg->ng',w,a['values'][index],optimize=True).astype(np.float32)


def molecular(s,a):
 n=len(s['xy']);o=np.zeros(n,np.float32);o[a['indices']]=1;z=np.zeros((n,a['values'].shape[1]),np.float32);z[a['indices']]=a['values'];op=graph(s['sample']);m1=np.asarray(op@o);m2=np.asarray(op@m1);v1=np.asarray(op@z);v2=np.asarray(op@v1);v1/=np.maximum(m1[:,None],1e-6);v2/=np.maximum(m2[:,None],1e-6)
 return o,z,v1.astype(np.float32),v2.astype(np.float32),m1,m2


def feature_parts(s,a,selected,lam):
 dist,index=nearest(s,a);interp=interpolate(a,dist,index,selected['best_family'],selected['best_param']);fuse=(lam*s['image']+(1-lam)*interp).astype(np.float32)
 o,z,v1,v2,m1,m2=molecular(s,a);op=graph(s['sample']);x=s['projection'];x1=np.asarray(op@x,dtype=np.float32);x2=np.asarray(op@x1,dtype=np.float32)
 shared=np.column_stack([x,x1,x2,o,m1,m2,np.log1p(dist[:,0]/100),np.sum(dist<=300,axis=1)/16]).astype(np.float32)
 gene=np.stack([s['image'],fuse,z,v1,v2],axis=2);return shared,gene,fuse,dist[:,0]


def design(shared,gene,g):return np.column_stack([shared,gene[:,g,:],np.ones(len(shared))]).astype(np.float64)


def b4_prediction(parts,weights):
 shared,gene,fuse,_=parts;y=np.empty_like(fuse)
 for g in range(fuse.shape[1]):y[:,g]=fuse[:,g]+design(shared,gene,g)@weights[g]
 return y


def predict(s,a,hidden,model,param,selected=None,lam=None,weights=None,mean=None):
 # No full expression matrix or hidden target is accepted by this function.
 if model=='B0':return np.broadcast_to(mean,(len(hidden),len(mean))).astype(np.float32).copy()
 if model=='B1':return s['image'][hidden].copy()
 if model in ['B2a','B2b']:
  d,i=nearest(s,a);return interpolate(a,d,i,model,param)[hidden]
 if model=='B3':
  d,i=nearest(s,a);p=interpolate(a,d,i,selected['best_family'],selected['best_param']);return (lam*s['image'][hidden]+(1-lam)*p[hidden]).astype(np.float32)
 if model=='B4':return b4_prediction(feature_parts(s,a,selected,lam),weights)[hidden]
 raise ValueError(model)


def control_inputs(s,a,name,seed):
 rng=np.random.default_rng(seed);s={**s};a={k:v.copy() for k,v in a.items()}
 if name=='NC1':a['values']=a['values'][rng.permutation(len(a['values']))]
 elif name=='NC2':p=rng.permutation(len(a['indices']));a['indices']=a['indices'][p];a['xy']=a['xy'][p]
 elif name=='NC3':p=rng.permutation(len(s['xy']));s['image']=s['image'][p];s['projection']=s['projection'][p]
 return s,a


def fit_residual(designs, residuals, alpha=10000.):
    """Summed-error ridge; the last (intercept) coefficient is unpenalized."""
    weights=[]
    for x,y in zip(designs,residuals):
        penalty=np.eye(x.shape[1])*alpha
        penalty[-1,-1]=0
        weights.append(np.linalg.solve(x.T@x+penalty,x.T@y))
    return np.asarray(weights)
