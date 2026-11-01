"""Locked external target reveal, C5 inference and all-run leakage checks."""
from support import *
import argparse,gzip
from threadpoolctl import threadpool_limits
sys.path.insert(0,str(OUT/'deps'));import pyzipper

def targets():
 check_lock();assert (B/'MASK_LOCK.json').exists();stamp=B/'logs/EXTERNAL_LABEL_ACCESS.json'
 if not stamp.exists():dump(stamp,{'timestamp':now(),'spec_sha256':sha(A/'C7A_SPEC_LOCK.json'),'mask_lock_sha256':sha(B/'MASK_LOCK.json'),'source_commit':commit(),'purpose':'Final evaluation only; no model/preprocessing/sampling changes permitted after this timestamp.'})
 norm=json.loads((ROOT/'reports/C0/C0_NORMALIZATION.json').read_text());genes=list(dict.fromkeys(sum(panels().values(),[])));idx=[norm['genes'].index(g) for g in genes]
 with pyzipper.AESZipFile(OUT/'raw/count-matrices.zip') as z:
  for sample in SAMPLES:
   p=CACHE/(sample+'_targets.npz')
   if p.exists():continue
   with np.load(CACHE/(sample+'_metadata.npz')) as f:bc=f['barcodes'].copy()
   frame=pd.read_csv(gzip.GzipFile(fileobj=z.open('count-matrices/'+sample+'.tsv.gz',pwd=b'zNLXkYk3Q9znUseS')),sep='\t',index_col=0);assert np.array_equal(frame.index.astype(str),bc);raw=frame[genes].to_numpy();assert np.isfinite(raw).all() and (raw>=0).all() and np.array_equal(raw,np.rint(raw)),'External target contract incompatible: expected integer counts'
   y=((np.log1p(raw)-np.array(norm['mean'])[idx])/np.array(norm['scale'])[idx]).astype('float32');np.savez_compressed(p,raw_counts=raw,normalized_targets=y,genes=np.array(genes),barcodes=bc);dump(p.with_suffix('.json'),{'sha256':sha(p),'raw_archive_sha256':sha(OUT/'raw/count-matrices.zip'),'normalization_sha256':sha(ROOT/'reports/C0/C0_NORMALIZATION.json'),'labels_used_only_after_locks':True})

def run(r,panel,budget,model,control=None,control_seed=0):
 run_id=rid('eval',[r['mask_id'],panel,budget,model,control,control_seed]);path=B/'logs/runs'/(run_id+'.json')
 if path.exists():return json.loads(path.read_text())
 start=time.perf_counter();sample=r['sample_id'];s=scene(sample,panel);y=truth(sample,panel);n=len(y)
 with np.load(ROOT/r['mask_path']) as f:order=f['order'][:count(n,budget)];forbidden=f['forbidden'];assert np.array_equal(f['barcodes'],s['barcodes'])
 obs=np.zeros(n,bool);obs[order]=True;assert not forbidden[obs].any();h=np.flatnonzero(~obs if r['regime']=='V1' else forbidden)
 a=c5.observed_payload(y,obs,s['xy']);orig=s
 if control:s,a=c5.control_inputs(s,a,control,control_seed)
 with core.PeakRSS() as mem:
  pred=predict(s,a,h,panel,r['regime'],model);assert np.isfinite(pred).all()
  for mode in ['nan','random']:
   corrupt=y.copy();corrupt[~obs]=np.nan if mode=='nan' else np.random.default_rng(716).normal(size=corrupt[~obs].shape);alt=c5.observed_payload(corrupt,obs,orig['xy']);ss=orig
   if control:ss,alt=c5.control_inputs(ss,alt,control,control_seed)
   assert np.array_equal(pred,predict(ss,alt,h,panel,r['regime'],model)),'Hidden-label leakage'
  metrics=core.gene_metrics(y[h],pred);summary=core.summarize(metrics);d=cKDTree(orig['xy'][order]).query(orig['xy'][h])[0];bins=[]
  for lo,hi in zip([0,150,300,600,1200,2400],[150,300,600,1200,2400,np.inf]):
   take=(d>=lo)&(d<hi);num=int(take.sum())
   if not num:continue
   err=y[h][take]-pred[take];corr=core.column_corr(y[h][take],pred[take]) if num>=20 else np.full(len(panels()[panel]),np.nan);bins.append({'lower_um':lo,'upper_um':hi if np.isfinite(hi) else None,'n_hidden':num,'mse':float(np.mean(err**2)),'mae':float(np.mean(abs(err))),'median_pearson':core.finite_stat(corr,'median')})
  pp=B/'predictions'/(run_id+'.npz');pp.parent.mkdir(exist_ok=True);np.savez(pp,prediction=pred,hidden_indices=h,hidden_barcodes=orig['barcodes'][h],observed_indices=order,observed_barcodes=orig['barcodes'][order],genes=np.array(panels()[panel]),distance_um=d,**metrics)
 rec={**r,'run_id':run_id,'panel':panel,'model':model,'encoder':ENC,'nominal_budget':budget,'realized_budget':len(order)/n,'n_observed':len(order),'n_hidden':len(h),'n_unmeasured':int((~obs).sum()),'control':control,'control_seed':control_seed,**summary,'distance_bins':bins,'hidden_corruption_checks':2,'feature_cache_hash':orig['feature_cache_hash'],'image_prediction_hash':sha(CACHE/(sample+'_'+panel+'_image.npz')),'model_contract_hash':sha(A/'C7A_SPEC_LOCK.json'),'target_cache_hash':sha(CACHE/(sample+'_targets.npz')),'prediction_path':str(pp.relative_to(ROOT)),'prediction_sha256':sha(pp),'runtime_seconds':time.perf_counter()-start,'peak_rss_mb':mem.peak/1024**2,'source_commit':commit()};dump(path,rec);registry('evaluation',run_id,rec);return rec

def evaluate():
 check_lock();targets();rows=[json.loads(p.read_text()) for p in sorted((B/'masks').glob('*.json'))];done=0
 for r in rows:
  for panel in panels():
   for budget in [.1,.2,.4]:
    for model in ['E0','E1','E2','E3']:run(r,panel,budget,model);done+=1
  if done%240==0:print('EVALUATED',done,'/8448',flush=True)
 nc=0
 for r in rows:
  if r['sample_id']!='A1' or r['strategy']!='random' or r['seed']!=0 or r['geometry']!=(0 if r['regime']=='V1' else 1):continue
  for panel in panels():
   for control in ['NC1','NC2','NC3']:
    for seed in range(3):
     for model in ['E1','E2','E3']:run(r,panel,.2,model,control,seed);nc+=1
 assert done==8448 and nc==108;print('FINAL COMPLETE',done,nc,flush=True)
if __name__=='__main__':
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):evaluate()
