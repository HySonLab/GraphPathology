"""External label-free nested trajectories, promoted from C7."""
import numpy as np
def count(n,b): return int(np.floor(n*b+.5))
def choose(raw,eligible,k,strategy,seed):
 if strategy not in ['random','spatial_FPS']:raise ValueError('Unknown sampling strategy')
 if k>len(eligible) or k<1:raise ValueError('Insufficient eligible locations')
 xy=np.asarray(raw['xy'],float);eligible=np.asarray(eligible,int)
 if strategy=='random':return np.random.default_rng(seed).permutation(eligible)[:k]
 z=xy[eligible];d=np.sum((z-z.mean(0))**2,axis=1);j=int(np.argmax(d));minimum=np.full(len(z),np.inf);order=[]
 for _ in range(k):
  order.append(int(eligible[j]));minimum=np.minimum(minimum,np.sum((z-z[j])**2,axis=1));minimum[j]=-np.inf;j=int(np.argmax(minimum))
 return np.array(order)


def region(xy,regime,geometry):
 if regime not in ['V1','V2'] or geometry not in [0,1,2]:raise ValueError('Invalid regional protocol')
 n=len(xy);forbidden=np.zeros(n,bool)
 if regime=='V1':return forbidden
 z=xy-xy.mean(0);_,_,v=np.linalg.svd(z,full_matrices=False);axis=v[0];axis=axis if axis[0]>=0 else -axis;p=z@axis;order=np.lexsort((np.arange(n),p));center=order[int(np.floor((n-1)*[.25,.5,.75][geometry]))];d=np.sum((xy-xy[center])**2,axis=1);forbidden[np.lexsort((np.arange(n),d))[:count(n,.4)]]=True;return forbidden

