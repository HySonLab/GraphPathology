"""Frozen C9 numerical prediction contract, with no archive I/O or fitting."""
import numpy as np
from scipy.spatial import cKDTree
from . import core
def moments(arrays):
 mean=np.mean([a.mean(0) for a in arrays],axis=0);var=np.mean([np.mean((a-mean)**2,axis=0) for a in arrays],axis=0);sd=np.sqrt(var)
 return {'center':mean.tolist(),'scale':np.where(sd<1e-8,1,sd).tolist(),'constant_fallback':(sd<1e-8).tolist()}


def _predict(xy,image,y,observed,hidden,regime,scaler,control='',control_seed=0):
 # Slice immediately: no hidden target can enter any input feature.
 values=np.asarray(y)[observed].copy();positions=xy[observed].copy();im=image[hidden].copy();rng=np.random.default_rng(control_seed)
 if control=='NC1':values=values[rng.permutation(len(values))]
 if control=='NC2':positions=positions[rng.permutation(len(positions))]
 if control=='NC3':im=image[rng.permutation(len(image))][hidden]
 dist,idx=cKDTree(positions).query(xy[hidden],k=min(16,len(observed)))
 if idx.ndim==1:idx=idx[:,None];dist=dist[:,None]
 dist=np.round(dist,6);weights=1/np.maximum(dist,1e-6)**2 if regime=='V1' else np.ones_like(dist);weights/=weights.sum(1,keepdims=True)
 at=np.array(scaler['destination']['center']);bt=np.array(scaler['destination']['scale']);asrc=np.array(scaler['source']['center']);bs=np.array(scaler['source']['scale'])
 st=np.einsum('nk,nkg->ng',weights,values[idx],optimize=True)
 z=np.einsum('nk,nkg->ng',weights,((values-at)/bt)[idx],optimize=True);harm=at+bt*z
 calibrated_image=at+bt*((im-asrc)/bs)
 return {'M0':im,'M1':st,'M2':harm,'M3':.2*calibrated_image+.8*harm},dist[:,0]


def score(y,p):
 met=core.gene_metrics(y,p);r={}
 for key in ['pearson','spearman','mse','mae']:
  for stat in ['mean','median']:r[stat+'_'+key]=core.finite_stat(met[key],stat)
 for t in [0,.1,.2,.3]:r['fraction_pearson_gt_'+str(t)]=float(np.mean(np.isfinite(met['pearson'])&(met['pearson']>t)))
 r['n_valid_pearson']=int(np.isfinite(met['pearson']).sum());return r,met


def predict(xy,image,y,observed,hidden,regime,scaler,control='',control_seed=0):
    xy=np.asarray(xy);image=np.asarray(image);y=np.asarray(y)
    observed=np.asarray(observed);hidden=np.asarray(hidden)
    if xy.ndim!=2 or xy.shape[1]!=2 or not np.isfinite(xy).all():raise ValueError('Finite Nx2 coordinates required')
    if image.ndim!=2 or image.shape!=y.shape or len(image)!=len(xy) or not np.isfinite(image).all():raise ValueError('Aligned finite image predictions and target shape required')
    for name,ix in [('observed',observed),('hidden',hidden)]:
        if ix.ndim!=1 or ix.dtype.kind not in 'iu' or not len(ix) or len(np.unique(ix))!=len(ix) or ix.min()<0 or ix.max()>=len(xy):raise ValueError('Nonempty unique integer '+name+' indices required')
    if np.intersect1d(observed,hidden).size:raise ValueError('Observed/scored overlap')
    if regime not in ['V1','V2'] or control not in ['', 'NC1','NC2','NC3']:raise ValueError('Unknown regime or control')
    if not np.isfinite(y[observed]).all():raise ValueError('Observed vectors must be finite')
    for domain in ['source','destination']:
        for name in ['center','scale']:
            a=np.asarray(scaler[domain][name])
            if a.shape!=(y.shape[1],) or not np.isfinite(a).all() or name=='scale' and (a<=0).any():raise ValueError('Invalid calibration '+domain+' '+name)
    return _predict(xy,image,y,observed,hidden,regime,scaler,control,control_seed)
