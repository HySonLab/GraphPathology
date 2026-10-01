"""Label-free external observation masks with frozen geometric stress regions."""
from support import *
def choose(raw,eligible,k,strategy,seed):
 xy=np.asarray(raw['xy'],float);eligible=np.asarray(eligible,int)
 if strategy=='random':return np.random.default_rng(seed).permutation(eligible)[:k]
 z=xy[eligible];d=np.sum((z-z.mean(0))**2,axis=1);j=int(np.argmax(d));minimum=np.full(len(z),np.inf);order=[]
 for _ in range(k):
  order.append(int(eligible[j]));minimum=np.minimum(minimum,np.sum((z-z[j])**2,axis=1));minimum[j]=-np.inf;j=int(np.argmax(minimum))
 return np.array(order)
def region(xy,regime,geometry):
 n=len(xy);forbidden=np.zeros(n,bool)
 if regime=='V1':return forbidden
 z=xy-xy.mean(0);_,_,v=np.linalg.svd(z,full_matrices=False);axis=v[0];axis=axis if axis[0]>=0 else -axis;p=z@axis;order=np.lexsort((np.arange(n),p));center=order[int(np.floor((n-1)*[.25,.5,.75][geometry]))];d=np.sum((xy-xy[center])**2,axis=1);forbidden[np.lexsort((np.arange(n),d))[:count(n,.4)]]=True;return forbidden

def masks():
 check_lock();assert not (B/'MASK_LOCK.json').exists();records=[]
 for sample in SAMPLES:
  with np.load(CACHE/(sample+'_metadata.npz')) as f:xy=f['xy_um'];bc=f['barcodes']
  n=len(bc)
  for regime in ['V1','V2']:
   for geometry in ([0] if regime=='V1' else [0,1,2]):
    forbidden=region(xy,regime,geometry);eligible=np.flatnonzero(~forbidden)
    for strategy in ['random','spatial_FPS']:
     for seed in (range(10) if strategy=='random' else [0]):
      k=count(n,.4);order=choose({'xy':xy},eligible,k,strategy,seed)
      for mode in ['nan','random']:
       hidden=np.full((n,280),np.nan) if mode=='nan' else np.random.default_rng(919).normal(size=(n,280));assert np.array_equal(order,choose({'xy':xy,'expression':hidden},eligible,k,strategy,seed))
      pid=rid('mask',[sample,regime,geometry,strategy,seed]);p=B/'masks'/(pid+'.npz');p.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(p,order=order,ordered_barcodes=bc[order],forbidden=forbidden,barcodes=bc,xy_um=xy)
      r={'mask_id':pid,'sample_id':sample,'patient_id':'HER2ST_'+sample[0],'regime':regime,'geometry':geometry,'strategy':strategy,'seed':seed,'mask_path':str(p.relative_to(ROOT)),'mask_sha256':sha(p),'n_locations':n,'n_forbidden':int(forbidden.sum()),'selection_corruption_checks':2};dump(p.with_suffix('.json'),r);records.append(r);registry('mask',pid,r)
 assert len(records)==352;dump(B/'MASK_LOCK.json',{'timestamp':now(),'source_commit':commit(),'spec_hash':sha(A/'C7A_SPEC_LOCK.json'),'mask_files':{r['mask_path']:r['mask_sha256'] for r in records},'all_masks_frozen_before_expression_access':True});print('FROZEN',len(records),'mask trajectories',flush=True)
if __name__=='__main__':masks()
