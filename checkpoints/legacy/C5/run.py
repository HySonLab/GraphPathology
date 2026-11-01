"""Validation-only selection, supervised training on training patients, then held-out reconstruction."""
from common import *
import argparse,torch
from threadpoolctl import threadpool_limits
from engine import *

def training_mean(cohort,panel):
 path=OUT/'models'/f'B0_{cohort}_{panel}.npz'
 if not path.exists():
  rs=[next(r for r in masks(s) if r['regime']=='random' and r['nominal_coverage']==.6 and r['replicate']==0) for s in training(cohort)];y=np.concatenate([anchors(r['sample_id'],panel,r)['values'] for r in rs]);np.savez_compressed(path,mean=y.mean(0,dtype=np.float64).astype(np.float32));registry('fit',rid('B0',[cohort,panel]),{'model':'B0','panel':panel,'train_slides':training(cohort),'observed_training_masks':[{r['mask_id']:r['checksum']} for r in rs],'checkpoint_path':str(path.relative_to(ROOT)),'checkpoint_hash':sha(path)})
 with np.load(path) as f:return f['mean'].copy()
def selection_record(model,panel,regime,encoder,param,values):
 r={'model':model,'panel':panel,'regime':regime,'encoder':encoder,'parameter':param,'validation_mean_MSE':float(np.mean(values)),'n_validation_masks':len(values),'evaluation_slides':SPLITS['splits']['validation']};registry('selection',rid('select',r),r);return r['validation_mean_MSE']
def select_simple():
 sample=SPLITS['splits']['validation'][0];selected={}
 for panel in PANELS:
  selected[panel]={};truth=labels(sample,panel)
  for regime in REGIMES:
   rows=[r for r in masks(sample) if r['regime']==regime];cache=[]
   for r in rows:
    s=scene('primary',ENCODERS[0],panel,sample);a=anchors(sample,panel,r);d,i=nearest(s,a);h=np.flatnonzero(~mask_array(r));cache.append((r,a,d,i,h))
   candidates={}
   for family,grid in [('B2a',cfg()['B2a']['k_grid']),('B2b',cfg()['B2b']['power_grid'])]:
    scored=[]
    for param in grid:
     losses=[float(np.mean((interpolate(a,d,i,family,param)[h]-truth[h])**2)) for r,a,d,i,h in cache];score=selection_record(family,panel,regime,'none',param,losses);scored.append((score,param))
    best=min(scored);candidates[family]={'param':best[1],'mse':best[0]}
   family=min(candidates,key=lambda k:(candidates[k]['mse'],k));param=candidates[family]['param'];record={'B2a':candidates['B2a']['param'],'B2b':candidates['B2b']['param'],'best_family':family,'best_param':param,'lambda':{}}
   for enc in ENCODERS:
    image=scene('primary',enc,panel,sample)['image'];scores=[]
    for lam in cfg()['B3']['lambda_grid']:
     losses=[float(np.mean((lam*image[h]+(1-lam)*interpolate(a,d,i,family,param)[h]-truth[h])**2)) for r,a,d,i,h in cache];scores.append((selection_record('B3',panel,regime,enc,lam,losses),lam))
    record['lambda'][enc]=min(scores)[1]
   selected[panel][regime]=record;print('SIMPLE',panel,regime,record,flush=True)
 dump(OUT/'SIMPLE_SELECTION.json',selected);return selected

def fit_b4(cohort,encoder,panel,selected,alphas):
 t=time.perf_counter();shared=[];gene=[];response=[];used=[]
 for sample in training(cohort):
  s=scene(cohort,encoder,panel,sample);truth=labels(sample,panel)
  for r in masks(sample):
   a=anchors(sample,panel,r);spec,lam=fusion_spec(selected,panel,r['regime'],encoder);parts=feature_parts(s,a,spec,lam);h=np.flatnonzero(~mask_array(r));rng=np.random.default_rng(int(r['seed']));h=rng.choice(h,min(256,len(h)),replace=False);shared.append(parts[0][h]);gene.append(parts[1][h]);response.append(truth[h]-parts[2][h]);used.append({'mask_id':r['mask_id'],'checksum':r['checksum'],'supervised_rows_sha256':hashlib.sha256(h.tobytes()).hexdigest(),'n_supervised':len(h)})
 sh=np.concatenate(shared);ge=np.concatenate(gene);target=np.concatenate(response);dim=sh.shape[1]+ge.shape[2]+1;weights={alpha:np.empty((len(PANELS[panel]),dim),np.float64) for alpha in alphas}
 with core.PeakRSS() as memory:
  penalty=np.eye(dim);penalty[-1,-1]=0
  for g in range(target.shape[1]):
   x=design(sh,ge,g);gram=x.T@x;rhs=x.T@target[:,g]
   for alpha in alphas:weights[alpha][g]=np.linalg.solve(gram+float(alpha)*penalty,rhs)
  for alpha in alphas:
   path=OUT/'models'/f'B4_{cohort}_{encoder}_{panel}_a{alpha}.npz';np.savez_compressed(path,weights=weights[alpha]);record={'model':'B4','panel':panel,'encoder':encoder,'cohort':cohort,'alpha':alpha,'train_slides':training(cohort),'n_train':len(target),'training_masks':used,'source_commit':commit(),'checkpoint_path':str(path.relative_to(ROOT)),'checkpoint_hash':sha(path),'scaler_path':str((OUT/'models'/f'PCA_{cohort}_{encoder}.npz').relative_to(ROOT)),'scaler_hash':sha(OUT/'models'/f'PCA_{cohort}_{encoder}.npz'),'parameters':int(weights[alpha].size),'train_seconds':time.perf_counter()-t,'peak_rss_mb':memory.peak/1024**2,'seed':0,'no_expression_scaler':True};dump(path.with_suffix('.json'),record);registry('fit',rid('B4',[cohort,encoder,panel,alpha]),record)
 return weights

def select():
 selected=select_simple();b4={};sample=SPLITS['splits']['validation'][0]
 for panel in PANELS:
  training_mean('primary',panel);b4[panel]={};truth=labels(sample,panel)
  for encoder in ENCODERS:
   weights=fit_b4('primary',encoder,panel,selected,cfg()['B4']['alpha_grid']);losses={a:[] for a in weights};s=scene('primary',encoder,panel,sample)
   for r in masks(sample):
    a=anchors(sample,panel,r);sp,lam=fusion_spec(selected,panel,r['regime'],encoder);parts=feature_parts(s,a,sp,lam);h=np.flatnonzero(~mask_array(r))
    for alpha,w in weights.items():p=b4_prediction(parts,w)[h];losses[alpha].append(float(np.mean((p-truth[h])**2)))
   scores=[(selection_record('B4',panel,'all',encoder,alpha,vals),alpha) for alpha,vals in losses.items()];b4[panel][encoder]=min(scores)[1];print('B4 SELECT',panel,encoder,b4[panel][encoder],flush=True)
 lock={'timestamp':now(),'source_commit':commit(),'simple':selected,'B4_alpha':b4,'spec_sha256':sha(OUT/'C5_SPEC_LOCK.json'),'primary_hidden_labels_used':False};dump(OUT/'SELECTION_LOCK.json',lock)

def weights_for(cohort,encoder,panel,selection):
 alpha=selection['B4_alpha'][panel][encoder];path=OUT/'models'/f'B4_{cohort}_{encoder}_{panel}_a{alpha}.npz'
 with np.load(path) as f:w=f['weights'].copy()
 r=json.loads(path.with_suffix('.json').read_text());assert sha(path)==r['checkpoint_hash'];return w,r

def evaluate(cohort,sample,panel,r,model,encoder,selection,control=None,seed=0):
 identity=[cohort,sample,panel,r['mask_id'],model,encoder,control,seed];run_id=rid('eval',identity);path=OUT/'logs/runs'/f'{run_id}.json'
 if path.exists():return
 t=time.perf_counter();observed=mask_array(r);h=np.flatnonzero(~observed);enc=ENCODERS[0] if encoder=='none' else encoder;s=scene(cohort,enc,panel,sample);a=anchors(sample,panel,r);spec=selection['simple'][panel][r['regime']];lam=spec['lambda'].get(encoder);mean=training_mean(cohort,panel);w,mr=weights_for(cohort,encoder,panel,selection) if model=='B4' else (None,None);param=spec.get(model);base_s=s
 s,a=control_inputs(s,a,control,seed)
 with core.PeakRSS() as memory:
  p=predict(s,a,h,model,param,spec,lam,w,mean);assert p.shape==(len(h),len(PANELS[panel])) and np.isfinite(p).all()
  # Access full held-out expression only inside this audit boundary and metric computation.
  y=labels(sample,panel);errors=[]
  for corruption in ['nan','random']:
   changed=y.copy();changed[h]=np.nan if corruption=='nan' else np.random.default_rng(3901).normal(size=(len(h),y.shape[1]))
   ca=observed_payload(changed,observed,base_s['xy']);cs,ca=control_inputs(base_s,ca,control,seed);q=predict(cs,ca,h,model,param,spec,lam,w,mean);error=float(np.max(abs(p-q)));assert np.array_equal(p,q),('CORRUPTION INVARIANCE FAILED',run_id,corruption,error);errors.append(error)
  metrics=core.gene_metrics(y[h],p);summary=core.summarize(metrics);dist=nearest(base_s,anchors(sample,panel,r))[0][h,0];edges=[0,150,300,600,1200,2400,np.inf];distance=[]
  for lo,hi in zip(edges[:-1],edges[1:]):
   take=(dist>=lo)&(dist<hi);n=int(take.sum())
   if not n:continue
   truth=y[h][take];pred=p[take];pc=core.column_corr(truth,pred) if n>=20 else np.full(y.shape[1],np.nan)
   distance.append({'lower_um':lo,'upper_um':hi if np.isfinite(hi) else None,'n_hidden':n,'mse':float(np.mean((truth-pred)**2)),'mae':float(np.mean(abs(truth-pred))),'median_pearson':core.finite_stat(pc,'median')})
  pp=OUT/'predictions'/f'{run_id}.npz';np.savez(pp,prediction=p,hidden_barcodes=base_s['barcodes'][h],observed_barcodes=base_s['barcodes'][observed],genes=np.array(PANELS[panel]),mask_checksum=np.array(r['checksum']),distance_um=dist,**metrics)
  d={'run_id':run_id,'cohort':cohort,'phase':'primary' if cohort=='primary' else 'sensitivity','sample_id':sample,'patient_id':SPLITS['patients'][sample],'panel':panel,'mask_id':r['mask_id'],'mask_checksum':r['checksum'],'regime':r['regime'],'nominal_coverage':r['nominal_coverage'],'realized_coverage':r['realized_coverage'],'replicate':r['replicate'],'model':model,'encoder':encoder,'control':control,'seed':seed,'config':{'interpolation':spec,'lambda':lam,'alpha':selection['B4_alpha'][panel].get(encoder) if model=='B4' else None},'feature_cache_hash':base_s['feature_cache_hash'],'scaler_hash':sha(OUT/'models'/f'PCA_{cohort}_{enc}.npz'),'checkpoint_hash':mr['checkpoint_hash'] if mr else sha(OUT/'models'/f'B0_{cohort}_{panel}.npz') if model=='B0' else sha(OUT/'SELECTION_LOCK.json'),'frozen_image_prediction_hash':sha(OUT/'cache'/('_'.join([cohort,enc,panel,sample])+'_B1.npz')),'source_commit':commit(),'prediction_path':str(pp.relative_to(ROOT)),'prediction_sha256':sha(pp),'n_hidden':len(h),'n_observed':int(observed.sum()),'invariance_errors':errors,'distance_bins':distance,**summary,'inference_audit_metric_seconds':time.perf_counter()-t,'peak_rss_mb':memory.peak/1024**2}
 dump(path,d);registry('evaluation',run_id,{**d,'evaluation_slides':[sample]})
 if len(list((OUT/'logs/runs').glob('*.json')))%30==0:print('EVALUATED',len(list((OUT/'logs/runs').glob('*.json'))),cohort,sample,panel,model,flush=True)

def final():
 selection=json.loads((OUT/'SELECTION_LOCK.json').read_text());dump(OUT/'logs/TEST_ACCESS.json',{'timestamp':now(),'selection_sha256':sha(OUT/'SELECTION_LOCK.json'),'prior_C0_C4_results_known':True})
 for cohort in ['primary']+IDS:
  panels=list(PANELS) if cohort=='primary' else ['S2-50','S2-280'];sample=SPLITS['splits']['test'][0] if cohort=='primary' else cohort
  for panel in panels:
   training_mean(cohort,panel)
   if cohort!='primary':
    for enc in ENCODERS:fit_b4(cohort,enc,panel,selection['simple'],[selection['B4_alpha'][panel][enc]])
   for r in masks(sample,cohort!='primary'):
    for model in ['B0','B2a','B2b']:evaluate(cohort,sample,panel,r,model,'none',selection)
    for encoder in ENCODERS:
     for model in ['B1','B3','B4']:evaluate(cohort,sample,panel,r,model,encoder,selection)
  SCENES.clear();PRED.clear()
 sample=SPLITS['splits']['test'][0];r=next(r for r in masks(sample) if r['regime']=='random' and r['nominal_coverage']==.2 and r['replicate']==0)
 for control in ['NC1','NC2','NC3']:
  for seed in range(3):
   if control!='NC3':
    for model in ['B2a','B2b']:evaluate('primary',sample,'S2-50',r,model,'none',selection,control,seed)
   for encoder in ENCODERS:
    for model in (['B1','B3','B4'] if control=='NC3' else ['B3','B4']):evaluate('primary',sample,'S2-50',r,model,encoder,selection,control,seed)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['select','final']);args=parser.parse_args();check_lock();torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):globals()[args.stage]()
