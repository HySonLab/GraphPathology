import numpy as np
def policy_view(raw):
 # Deliberately whitelist only allowed image/coordinate inputs; expression is never forwarded.
 return {'xy':raw['xy'],'projection':raw['projection']}


def features(view,eligible,strategy,alpha=.5,control=None,seed=0):
 xy=np.asarray(view['xy'],np.float64).copy();m=np.asarray(view['projection'],np.float64).copy()/np.sqrt(8)
 rng=np.random.default_rng(seed)
 if control=='NC2':m=m[rng.permutation(len(m))]
 if control=='NC3':xy=xy[rng.permutation(len(xy))]
 xy-=xy.mean(0);xy/=max(np.sqrt(np.mean(np.sum(xy**2,axis=1))),1e-12)
 if strategy=='A1':return xy
 if strategy=='A2':return m
 return np.column_stack([np.sqrt(alpha)*xy,np.sqrt(1-alpha)*m])


def choose(raw,eligible,n,strategy,seed=0,alpha=.5,control=None):
 view=policy_view(raw);eligible=np.asarray(eligible,int);assert n<=len(eligible)
 if strategy=='A0':
  order=np.random.default_rng(seed).permutation(eligible)[:n];return order,np.full(n,np.nan)
 x=features(view,eligible,strategy,alpha,control,seed);z=x[eligible];center=z.mean(0);d=np.einsum('ij,ij->i',z-center,z-center);j=int(np.argmax(d));minimum=np.full(len(z),np.inf);order=np.empty(n,int);scores=np.empty(n,float)
 for k in range(n):
  order[k]=eligible[j];scores[k]=d[j] if k==0 else minimum[j];diff=z-z[j];minimum=np.minimum(minimum,np.einsum('ij,ij->i',diff,diff));minimum[j]=-np.inf;j=int(np.argmax(minimum))
 assert len(np.unique(order))==n;return order,scores

def target_indices(n,selected,forbidden,regime):
 obs=np.zeros(n,bool);obs[selected]=True
 return np.flatnonzero(forbidden if regime=='contiguous' else ~obs),obs
