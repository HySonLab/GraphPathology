"""HER2ST gate features use destination-transformed anchors only."""
import numpy as np
from scipy.spatial import cKDTree
from scipy.special import expit

def features(xy, observed, hidden, y, destination):
    values=(np.asarray(y)[observed]-destination['center'])/destination['scale']
    dist,idx=cKDTree(xy[observed]).query(xy[hidden],k=min(16,len(observed)))
    if idx.ndim==1: idx=idx[:,None]; dist=dist[:,None]
    dist=np.round(dist,6)
    count=cKDTree(xy[observed]).query_ball_point(xy[hidden],600,return_length=True)
    variance=values[idx].var(1).mean(1)
    return np.column_stack([np.log1p(dist[:,0]/200),np.log1p(dist.mean(1)/200),np.log1p(count),np.log1p(variance)])
def weights(v,model):
 typ=model['family']
 if typ=='F0':return np.zeros(len(v['st']))
 if typ=='F1':return np.full(len(v['st']),model['lambda'])
 if typ=='F2':return np.array(model['weights'])[np.digitize(v['distance'],[300,600])]
 z=(v['features']-model['feature_mean'])/model['feature_scale'];return expit(model['coef'][0]+z@np.array(model['coef'][1:]))


def predict_B(v,model):
 w=weights(v,model);return v['st']+w[:,None]*(v['image']-v['st'])

