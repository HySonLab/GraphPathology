"""Deterministic label-free nested coverage policies; exact ordered selection replay."""
from support import *
from scipy.sparse.csgraph import connected_components

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

GEOMETRY_CACHE={}
def geometry(s,selected,h,obs):
 key=(s["sample"],s["feature_cache_hash"],hashlib.sha256(s["projection"].tobytes()).hexdigest(),hashlib.sha256(selected.tobytes()+h.tobytes()).hexdigest())
 if key in GEOMETRY_CACHE:return GEOMETRY_CACHE[key]
 d=cKDTree(s['xy'][selected]).query(s['xy'][h])[0];m=s['projection']/np.sqrt(8);md=cKDTree(m[selected]).query(m[h])[0];unmeasured=np.flatnonzero(~obs);graph=c5.graph(s['sample'])[unmeasured][:,unmeasured];_,lab=connected_components(graph,directed=False);sizes=np.bincount(lab)
 result={'mean_anchor_distance_um':float(d.mean()),'median_anchor_distance_um':float(np.median(d)),'p90_anchor_distance_um':float(np.quantile(d,.9)),'max_anchor_distance_um':float(d.max()),'spatial_coverage_radius_um':float(d.max()),'unmeasured_components_at_least_50':int((sizes>=50).sum()),'morphology_mean_nearest_distance':float(md.mean()),'morphology_p90_nearest_distance':float(np.quantile(md,.9)),'observed_morphology_mean_pairwise_squared_distance':float(2*len(selected)/(len(selected)-1)*np.var(m[selected],axis=0).sum())},d
 GEOMETRY_CACHE[key]=result;return result

def make(cohort,sample,enc,regime,strategy,seed=0,alpha=.5,control=None,validation=False):
 identity=[cohort,sample,enc,regime,strategy,seed,alpha,control];pid=rid('policy',identity);meta=A/'policies'/f'{pid}.json'
 if meta.exists():return json.loads(meta.read_text())
 s=scene(cohort,enc,'S2-50',sample);forbidden,ref=unavailable(sample,regime);eligible=np.flatnonzero(~forbidden);n=count(len(forbidden),.4);start=time.perf_counter();raw={'xy':s['xy'],'projection':s['projection']};order,scores=choose(raw,eligible,n,strategy,seed,alpha,control)
 # At C6A selection time EVERY expression vector is unrevealed. Audit all rows, all genes.
 errors=[]
 for corruption in ['nan','random']:
  fake=np.full((len(forbidden),280),np.nan,np.float32) if corruption=='nan' else np.random.default_rng(817).normal(size=(len(forbidden),280)).astype(np.float32)
  alt,score=choose({**raw,'expression':fake},eligible,n,strategy,seed,alpha,control);assert np.array_equal(order,alt) and np.array_equal(scores,score,equal_nan=True);errors.append(0)
 path=meta.with_suffix('.npz');np.savez_compressed(path,order=order,selected_barcodes=s['barcodes'][order],eligible=eligible,forbidden=forbidden,acquisition_scores=scores)
 r={'policy_id':pid,'cohort':cohort,'sample_id':sample,'patient_id':parent.SPLITS['patients'][sample],'encoder':enc,'regime':regime,'strategy':strategy,'seed':seed,'alpha':alpha if strategy=='A3' else None,'control':control,'validation':validation,'n_locations':len(forbidden),'n_eligible':len(eligible),'n_maximum_selected':n,'unavailable_region':ref,'policy_path':str(path.relative_to(ROOT)),'policy_sha256':sha(path),'source_commit':commit(),'timestamp':now(),'feature_cache_hash':s['feature_cache_hash'],'scaler_hash':sha(ROOT/'C5/models'/f'PCA_{cohort}_{enc}.npz'),'selection_corruption_errors':errors,'selection_seconds':time.perf_counter()-start,'reveal_protocol':'No ST revealed for selection. Each budget prefix is revealed as one batch after the complete trajectory lock.'};dump(meta,r);registry('policy',pid,r);print('POLICY',cohort,enc,regime,strategy,seed,control,flush=True);return r
