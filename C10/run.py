"""Frozen C9 masks, calibration, then held-out evaluation with per-run leakage tests."""
from contract10 import *
import argparse
from threadpoolctl import threadpool_limits
import policy as frozen_policy

def masks():
 check_lock();records=[];assert not (B/'MASK_LOCK.json').exists()
 for s in SAMPLES:
  with np.load(CACHE/(s+'_metadata.npz')) as f:xy=f['xy_um'];bc=f['barcodes']
  n=len(bc)
  for regime in ['V1','V2']:
   for geometry in ([0] if regime=='V1' else [0,1,2]):
    forbidden=frozen_policy.region(xy,regime,geometry);eligible=np.flatnonzero(~forbidden)
    for strategy in ['random','spatial_FPS']:
     for seed in (range(10) if strategy=='random' else [0]):
      order=frozen_policy.choose({'xy':xy},eligible,parent.count(n,.4),strategy,seed)
      ident='C10mask_'+hashlib.sha256(json.dumps([s,regime,geometry,strategy,seed]).encode()).hexdigest()[:18];p=B/'masks'/(ident+'.npz');p.parent.mkdir(exist_ok=True)
      np.savez_compressed(p,order=order,forbidden=forbidden,barcodes=bc,xy_um=xy)
      rec={'mask_id':ident,'sample_id':s,'patient_id':s,'cohort_role':'evaluation','regime':regime,'geometry':geometry,'strategy':strategy,'seed':seed,'mask_path':str(p.relative_to(ROOT)),'mask_sha256':sha(p),'n_locations':n,'label_free':True}
      dump(p.with_suffix('.json'),rec);registry('mask',ident,rec);records.append(rec)
 dump(B/'MASK_LOCK.json',{'timestamp':now(),'source_commit':commit(),'spec_sha256':sha(A/'SPEC_LOCK.json'),'files':{r['mask_path']:r['mask_sha256'] for r in records},'count':len(records)})
 print('MASKS',len(records),flush=True)

def evaluate():
 check_lock();obj=scaler();assert not cfg()['calibration_patients']
 rows=[json.loads(p.read_text()) for p in sorted((B/'masks').glob('*.json'))]
 stamp=B/'logs/EVALUATION_LABEL_ACCESS.json'
 if not stamp.exists():dump(stamp,{'timestamp':now(),'source_commit':commit(),'spec_sha256':sha(A/'SPEC_LOCK.json'),'mask_sha256':sha(B/'MASK_LOCK.json'),'frozen_C9_scaler_sha256':sha(ROOT/'C9/C9A/T1_SCALER.json'),'purpose':'Locked final evaluation, no tuning'})
 for s in SAMPLES:
  y=targets(s);np.savez_compressed(CACHE/(s+'_targets.npz'),normalized_targets=y,genes=genes())
  with np.load(CACHE/(s+'_metadata.npz')) as f:xy=f['xy_um'];bc=f['barcodes']
  with np.load(CACHE/(s+'_image.npz')) as f:im=f['prediction'];assert np.array_equal(bc,f['barcodes'])
  for r in [r for r in rows if r['sample_id']==s]:
   with np.load(ROOT/r['mask_path']) as f:order=f['order'];forbidden=f['forbidden'];assert np.array_equal(bc,f['barcodes'])
   assert sha(ROOT/r['mask_path'])==r['mask_sha256']
   for budget in [.1,.2,.4]:
    observed=order[:parent.count(len(y),budget)];obs=np.zeros(len(y),bool);obs[observed]=True;hidden=np.flatnonzero(~obs if r['regime']=='V1' else forbidden);assert not forbidden[observed].any()
    controls=[('',0)]
    if s==SAMPLES[0] and budget==.2 and r['strategy']=='random' and r['seed']==0 and r['geometry']==(0 if r['regime']=='V1' else 1):controls += [(c,j) for c in ['NC1','NC2','NC3'] for j in range(3)]
    for control,cs in controls:
     ident='C10eval_'+hashlib.sha256(json.dumps([r['mask_id'],budget,control,cs]).encode()).hexdigest()[:18];log=B/'logs/runs'/(ident+'.json')
     if log.exists():continue
     start=time.perf_counter()
     with core.PeakRSS() as mem:
      preds,dist=predict(xy,im,y,observed,hidden,r['regime'],obj,control,cs)
      for mode in ['nan','random']:
       altered=y.copy();altered[~obs]=np.nan if mode=='nan' else np.random.default_rng(921).normal(size=altered[~obs].shape)
       alt,_=predict(xy,im,altered,observed,hidden,r['regime'],obj,control,cs)
       for m in preds:assert np.array_equal(preds[m],alt[m]),(ident,m,'hidden leakage')
      metrics={};summary={};bins=[]
      for m,p in preds.items():
       assert np.isfinite(p).all();summary[m],metrics[m]=score(y[hidden],p)
       for lo,hi in zip([0,150,300,600,1200,2400],[150,300,600,1200,2400,np.inf]):
        take=(dist>=lo)&(dist<hi)
        if not take.any():continue
        err=y[hidden][take]-p[take];bins.append({'method':m,'lower_um':lo,'upper_um':hi if np.isfinite(hi) else None,'n_hidden':int(take.sum()),'mse':float(np.mean(err**2)),'mae':float(np.mean(abs(err))),'median_pearson':core.finite_stat(core.column_corr(y[hidden][take],p[take]),'median') if take.sum()>=20 else None})
      path=B/'predictions'/(ident+'.npz');path.parent.mkdir(exist_ok=True)
      arrays={**preds,'observed_indices':observed,'hidden_indices':hidden,'observed_barcodes':bc[observed],'hidden_barcodes':bc[hidden],'genes':np.array(genes()),'distance_um':dist}
      for m,met in metrics.items():arrays.update({m+'_'+k:v for k,v in met.items()})
      np.savez_compressed(path,**arrays)
     rec={**r,'run_id':ident,'budget':budget,'panel':'S2-50','control':control,'control_seed':cs,'n_observed':len(observed),'n_hidden':len(hidden),'realized_coverage':len(observed)/len(y),'summary':summary,'distance_bins':bins,'hidden_corruption_checks':6,'affine_equivalence_checked':True,'prediction_path':str(path.relative_to(ROOT)),'prediction_sha256':sha(path),'image_prediction_sha256':sha(CACHE/(s+'_image.npz')),'feature_cache_sha256':sha(CACHE/(s+'_features.npz')),'scaler_sha256':sha(ROOT/'C9/C9A/T1_SCALER.json'),'spec_sha256':sha(A/'SPEC_LOCK.json'),'runtime_seconds':time.perf_counter()-start,'peak_rss_mb':mem.peak/1024**2,'calibration_patients':[],'frozen_Wu_scaler_patients':obj['calibration_patients'],'image_weight':.2,'source_commit':commit()};dump(log,rec);registry('evaluation',ident,rec)
  print('EVALUATED',s,flush=True)
 print('COMPLETE',len(list((B/'logs/runs').glob('*.json'))),flush=True)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['masks','evaluate']);args=parser.parse_args()
 with threadpool_limits(4):globals()[args.stage]()
