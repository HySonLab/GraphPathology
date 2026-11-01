"""Freeze all label-free trajectories before test reveals, then score immutable C5 models."""
from support import *
from policy import make,geometry
import argparse,torch
from threadpoolctl import threadpool_limits

def evaluate(r,panel,budget,role,validation=False):
 run_id=rid('eval',[r['policy_id'],panel,budget,role]);path=A/'logs/runs'/f'{run_id}.json'
 if path.exists():return json.loads(path.read_text())
 start=time.perf_counter();cohort,enc,sample,regime=[r[k] for k in ['cohort','encoder','sample_id','regime']];s=scene(cohort,enc,panel,sample);n=len(s['xy']);k=count(n,budget)
 assert sha(ROOT/r['policy_path'])==r['policy_sha256']
 with np.load(ROOT/r['policy_path']) as f:selected=f['order'][:k].copy();forbidden=f['forbidden'].copy()
 h,obs=target_indices(n,selected,forbidden,regime);assert not forbidden[selected].any();a=observed_only(sample,panel,selected,validation)
 with core.PeakRSS() as memory:
  pred=reconstruction(s,a,h,cohort,enc,panel,regime,role);assert np.isfinite(pred).all();truth=expression(sample,panel,validation);errors=[]
  for name in ['nan','random']:
   changed=truth.copy();changed[~obs]=np.nan if name=='nan' else np.random.default_rng(916).normal(size=changed[~obs].shape)
   alt=c5.observed_payload(changed,obs,s['xy']);q=reconstruction(s,alt,h,cohort,enc,panel,regime,role);assert np.array_equal(pred,q);errors.append(0)
  metrics=core.gene_metrics(truth[h],pred);summary=core.summarize(metrics);diag,distance=geometry(s,selected,h,obs);bins=[]
  for lo,hi in zip([0,150,300,600,1200,2400],[150,300,600,1200,2400,np.inf]):
   take=(distance>=lo)&(distance<hi);num=int(take.sum())
   if not num:continue
   err=truth[h][take]-pred[take];corr=core.column_corr(truth[h][take],pred[take]) if num>=20 else np.full(len(PANELS[panel]),np.nan)
   bins.append({'lower_um':lo,'upper_um':hi if np.isfinite(hi) else None,'n_hidden':num,'mse':float(np.mean(err**2)),'mae':float(np.mean(abs(err))),'median_pearson':core.finite_stat(corr,'median')})
  pp=A/'predictions'/f'{run_id}.npz';np.savez(pp,prediction=pred,hidden_indices=h,hidden_barcodes=s['barcodes'][h],observed_indices=selected,observed_barcodes=s['barcodes'][selected],genes=np.array(PANELS[panel]),distance_um=distance,**metrics)
  rec={**{k:r[k] for k in ['policy_id','cohort','sample_id','patient_id','encoder','regime','strategy','seed','alpha','control','policy_path','policy_sha256']},'run_id':run_id,'phase':'validation' if validation else ('primary' if cohort=='primary' else 'sensitivity'),'panel':panel,'nominal_budget':budget,'realized_budget':k/n,'n_observed':k,'n_hidden':len(h),'n_unmeasured':int((~obs).sum()),'role':role,**provenance(cohort,enc,panel,regime,role,sample),**summary,**diag,'distance_bins':bins,'invariance_errors':errors,'source_commit':commit(),'prediction_path':str(pp.relative_to(ROOT)),'prediction_sha256':sha(pp),'runtime_seconds':time.perf_counter()-start,'peak_rss_mb':memory.peak/1024**2,'observed_label_reveal_sequence':'Exact prefix in observed_indices / observed_barcodes; no earlier molecular policy access'}
 dump(path,rec);registry('evaluation',run_id,{**rec,'evaluation_slides':[sample]});return rec

def select():
 assert not (A/'SELECTION_LOCK.json').exists();sample=parent.SPLITS['splits']['validation'][0];choices={};scores=[]
 for enc in ENCODERS:
  for strategy in ['A1','A2','A3']:
   for alpha in ([.25,.5,.75] if strategy=='A3' else [.5]):
    vals=[]
    for regime in REGIMES:
     r=make('primary',sample,enc,regime,strategy,alpha=alpha,validation=True)
     for panel in PANELS:
      for budget in BUDGETS:vals.append(evaluate(r,panel,budget,'primary',True)['mean_mse'])
    scores.append({'encoder':enc,'strategy':strategy,'alpha':alpha if strategy=='A3' else None,'validation_mean_MSE':float(np.mean(vals)),'n_settings':len(vals)})
  a=min([r for r in scores if r['encoder']==enc and r['strategy']=='A3'],key=lambda r:(r['validation_mean_MSE'],r['alpha']));best=min([r for r in scores if r['encoder']==enc and (r['strategy']!='A3' or r['alpha']==a['alpha'])],key=lambda r:(r['validation_mean_MSE'],r['strategy']));choices[enc]={'alpha':a['alpha'],'best_nonrandom_validation_strategy':best['strategy']}
 dump(A/'SELECTION_LOCK.json',{'timestamp':now(),'source_commit':commit(),'spec_sha256':sha(A/'C6A_SPEC_LOCK.json'),'choices':choices,'validation_scores':scores,'validation_patient':sample,'primary_test_hidden_labels_used':False});print('SELECTION LOCKED',choices,flush=True)

def trajectories():
 assert not (A/'TRAJECTORY_LOCK.json').exists();selection=json.loads((A/'SELECTION_LOCK.json').read_text());records=[]
 for cohort in ['primary']+IDS:
  sample=parent.SPLITS['splits']['test'][0] if cohort=='primary' else cohort
  for enc in ENCODERS:
   alpha=selection['choices'][enc]['alpha']
   for regime in REGIMES:
    for strategy in ['A0','A1','A2','A3']:
     for seed in (range(10) if strategy=='A0' else [0]):records.append(make(cohort,sample,enc,regime,strategy,seed,alpha if strategy=='A3' else .5))
 # Controls change only the sampling inputs; reconstruction and scoring use true images/coordinates.
 for enc in ENCODERS:
  for regime in REGIMES:
   for strategy,control in [('A1','NC3'),('A2','NC2'),('A3','NC2')]:
    for seed in range(3):records.append(make('primary',parent.SPLITS['splits']['test'][0],enc,regime,strategy,seed,selection['choices'][enc]['alpha'] if strategy=='A3' else .5,control))
 assert len(records)==296
 dump(A/'TRAJECTORY_LOCK.json',{'timestamp':now(),'source_commit':commit(),'selection_sha256':sha(A/'SELECTION_LOCK.json'),'all_primary_trajectories_frozen_before_hidden_label_access':True,'policies':{r['policy_id']:{'path':r['policy_path'],'sha256':r['policy_sha256']} for r in records}});print('ALL 296 TEST TRAJECTORIES FROZEN',flush=True)

def final():
 lock=json.loads((A/'TRAJECTORY_LOCK.json').read_text());access=A/'logs/TEST_ACCESS.json'
 if not access.exists():dump(access,{'timestamp':now(),'trajectory_lock_sha256':sha(A/'TRAJECTORY_LOCK.json'),'selection_sha256':sha(A/'SELECTION_LOCK.json')})
 records=[json.loads((A/'policies'/f'{key}.json').read_text()) for key in lock['policies']];done=0
 for r in records:
  for panel in PANELS:
   for budget in ([.1,.2] if r['control'] else BUDGETS):
    for role in ['primary','baseline']:evaluate(r,panel,budget,role);done+=1
  if done%64==0:print('EVALUATED',done,'of 4448',r['cohort'],r['strategy'],flush=True)
 assert done==4448;print('FINAL COMPLETE',done,flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['select','trajectories','final']);args=parser.parse_args();check_lock();torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):globals()[args.stage]()
