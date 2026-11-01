"""Nested patient-only fitting, sealed objects and held-out C8 evaluations."""
from contract import *
import argparse
from threadpoolctl import threadpool_limits

def fit_A():
 check_lock();assert not (A/'FIT_LOCK.json').exists();allfolds={}
 for fold,held in cfg()['folds'].items():
  set_forbidden(held);train=[s for s in SAMPLES if s not in held];assert not set(train)&set(held);selection=[];cache=[];fitted=[]
  for kind in ['T0','T1','T2']:
   obj=fit_scaler(train,kind);path=A/'scalers'/f'{fold}_{kind}.json';dump(path,obj);fitted.append(str(path.relative_to(ROOT)));registry('target_fit',rid([fold,kind,'outer_scaler']),{'phase':'C8A_fit','fold':fold,'transform':kind,'train_patients':train,'scaler_path':str(path.relative_to(ROOT)),'sha256':sha(path)})
   for sample in train:
    inner=[s for s in train if s!=sample];obj=fit_scaler(inner,kind);sp=A/'scalers'/f'{fold}_{kind}_inner_{sample}.json';dump(sp,obj);fitted.append(str(sp.relative_to(ROOT)))
    with allowed([sample]):y=target(sample)
    for r in masks(sample,inner=True):
     v=inputs(r,obj,y);t=convert(y[v['hidden']],obj);ident=rid(['inner',fold,kind,r['mask_id'],r['budget']]);p=A/'inner_cache'/(ident+'.npz');p.parent.mkdir(exist_ok=True);np.savez_compressed(p,**{k:v[k] for k in ['image','st','features','distance','scale','center']},truth=t)
     cache.append({'fold':fold,'transform':kind,'sample_id':sample,'regime':r['regime'],'geometry':r['geometry'],'seed':r['seed'],'budget':r['budget'],'mask_path':r['mask_path'],'path':str(p.relative_to(ROOT)),'sha256':sha(p),'scaler_path':str(sp.relative_to(ROOT))})
     for model in ['E0','E1','E2']:
      pred=predict_A(v,obj,model,r['regime']);met=fast_score(t,pred,v['scale']);rec={'fold':fold,'transform':kind,'sample_id':sample,'regime':r['regime'],'geometry':r['geometry'],'seed':r['seed'],'budget':r['budget'],'model':model,**met,'phase':'C8A_inner','train_patients':inner,'scaler_path':str(sp.relative_to(ROOT)),'cache_path':str(p.relative_to(ROOT))};selection.append(rec);registry('inner_evaluation',rid([ident,model]),rec)
  chosen,summary=select([r for r in selection if r['model']=='E2'],'transform');pd.DataFrame(selection).to_csv(A/'logs'/f'{fold}_INNER_METRICS.csv',index=False);dump(A/'logs'/f'{fold}_INNER_CACHE.json',cache);allfolds[fold]={'held_out_patients':held,'train_patients':train,'selected_transform':chosen,'inner_selection_summary':summary,'scaler_files':{p:sha(ROOT/p) for p in fitted},'inner_cache_index_sha256':sha(A/'logs'/f'{fold}_INNER_CACHE.json')};print('FITTED A',fold,chosen,flush=True)
 dump(A/'FIT_LOCK.json',{'timestamp':now(),'spec_sha256':sha(A/'C8A_SPEC_LOCK.json'),'folds':allfolds,'outer_evaluation_started':False});print('ALL A FOLD OBJECTS SEALED',flush=True)

def save_outer(stage,fold,r,obj,model,pred,v,y,object_path,control=None,control_seed=0,model_obj=None):
 start=time.perf_counter();path=OUT/stage;ident=rid([stage,fold,r['mask_id'],r['budget'],obj['transform'],model,control,control_seed]);log=path/'logs/runs'/(ident+'.json')
 if log.exists():return json.loads(log.read_text())
 obs=np.zeros(len(y),bool);obs[v['observed']]=True
 for val in [None,123]:
  alt=y.copy();alt[~obs]=np.nan if val is None else np.random.default_rng(val).normal(size=alt[~obs].shape);vv=inputs(r,obj,alt)
  if stage=='C8A':pp=predict_A(vv,obj,model,r['regime'])
  else:
   if control:vv,mo=control_inputs(vv,model_obj,control,control_seed);pp=predict_B(vv,mo)
   else:pp=predict_B(vv,model_obj)
  assert np.array_equal(pred,pp),'Hidden-label leakage'
 h=v['hidden'];t=convert(y[h],obj);metrics,gm=score(t,pred,v['scale']);common=inverse(pred,obj);cm=core.gene_metrics(y[h],common);ppath=path/'predictions'/(ident+'.npz');ppath.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(ppath,prediction_native=pred,prediction_common=common,observed_indices=v['observed'],hidden_indices=h,observed_barcodes=v['barcodes'][v['observed']],hidden_barcodes=v['barcodes'][h],genes=np.array(GENES),distance_um=v['distance'],**gm,common_gene_mse=cm['mse'],common_gene_mae=cm['mae'])
 bins=[]
 for lo,hi in zip([0,150,300,600,1200,2400],[150,300,600,1200,2400,np.inf]):
  take=(v['distance']>=lo)&(v['distance']<hi)
  if take.any():bins.append({'lower_um':lo,'upper_um':hi if np.isfinite(hi) else None,'n':int(take.sum()),'common_mse':float(np.mean((y[h][take]-common[take])**2))})
 rec={'phase':stage,'fold':fold,'sample_id':r['sample_id'],'patient_id':'HER2ST_'+r['sample_id'][0],'panel':'S2-50_40','transform':obj['transform'],'model':model,'regime':r['regime'],'geometry':r['geometry'],'seed':r['seed'],'budget':r['budget'],'mask_id':r['mask_id'],'mask_path':r['mask_path'],'mask_sha256':sha(ROOT/r['mask_path']),'scaler_sha256':sha(A/'scalers'/f'{fold}_{obj["transform"]}.json'),'calibration_object_path':str(object_path.relative_to(ROOT)),'calibration_object_sha256':sha(object_path),'prediction_path':str(ppath.relative_to(ROOT)),'prediction_sha256':sha(ppath),'n_hidden':len(h),'n_observed':int(obs.sum()),'control':control,'control_seed':control_seed,'hidden_corruption_checks':2,**metrics,'target_mean':float(t.mean()),'mean_target_variance':float(gm['target_variance'].mean()),'distance_bins':bins,'runtime_seconds':time.perf_counter()-start,'run_id':ident};dump(log,rec);registry('outer_evaluation' if not control else 'negative_control',ident,rec);return rec

def evaluate_A():
 check_lock();set_forbidden([]);lock=json.loads((A/'FIT_LOCK.json').read_text());stamp=A/'logs/OUTER_ACCESS.json'
 if not stamp.exists():dump(stamp,{'timestamp':now(),'fit_lock_sha256':sha(A/'FIT_LOCK.json'),'purpose':'Only evaluate held-out predictions; transforms already selected using inner train-only metrics'})
 n=0
 for fold,held in cfg()['folds'].items():
  for kind in ['T0','T1','T2']:
   path=A/'scalers'/f'{fold}_{kind}.json';assert sha(path)==lock['folds'][fold]['scaler_files'][str(path.relative_to(ROOT))];obj=json.loads(path.read_text());assert not set(obj['train_patients'])&set(held)
   for sample in held:
    with allowed([sample]):y=target(sample)
    for r in masks(sample):
     v=inputs(r,obj,y)
     for model in ['E0','E1','E2']:pred=predict_A(v,obj,model,r['regime']);save_outer('C8A',fold,r,obj,model,pred,v,y,path);n+=1
  print('EVALUATED A',fold,n,flush=True)
 assert n==4320

def fit_B():
 check_lock();assert json.loads((A/'C8A_STATUS.json').read_text())['status']=='PASS_TARGET_HARMONIZATION';assert not (B/'C8B_SPEC_LOCK.json').exists();folds={}
 for fold,held in cfg()['folds'].items():
  a=json.loads((A/'FIT_LOCK.json').read_text())['folds'][fold];kind=a['selected_transform'];index=json.loads((A/'logs'/f'{fold}_INNER_CACHE.json').read_text());records=[]
  for r in index:
   if r['transform']!=kind:continue
   assert r['sample_id'] in a['train_patients'] and r['sample_id'] not in held;assert sha(ROOT/r['path'])==r['sha256']
   with np.load(ROOT/r['path']) as f:v={k:f[k].copy() for k in ['image','st','features','distance','scale','center']};t=f['truth'].copy()
   records.append({**r,'v':v,'truth':t})
  cands=candidates(records);trials=[]
  for i,model in enumerate(cands):
   ident=f'{model["family"]}_{i:02}'
   for r in records:trials.append({**{k:r[k] for k in ['sample_id','regime','geometry','seed','budget']},'candidate':ident,'family':model['family'],**fast_score(r['truth'],predict_B(r['v'],model),r['v']['scale'])})
  best,summary=select(trials,'candidate');byfam={}
  for fam in ['F0','F1','F2','F3']:
   name,_=select([r for r in trials if r['family']==fam],'candidate');idx=int(name.rsplit('_',1)[1]);byfam[fam]=cands[idx]
  win=cands[int(best.rsplit('_',1)[1])];obj={'fold':fold,'train_patients':a['train_patients'],'held_out_patients':held,'transform':kind,'selected_family':win['family'],'selected_model':win,'family_models':byfam,'selection_summary':summary,'fit_timestamp':now(),'inner_cache_index_sha256':sha(A/'logs'/f'{fold}_INNER_CACHE.json')};p=B/'calibration'/f'{fold}.json';dump(p,obj);pd.DataFrame(trials).to_csv(B/'logs'/f'{fold}_CALIBRATION_METRICS.csv',index=False)
  for i,model in enumerate(cands):registry('fusion_candidate_fit',rid(['fit_B',fold,i]),{'phase':'C8B_fit','fold':fold,'model':model['family'],'candidate':i,'train_patients':a['train_patients'],'parameters':model,'calibration_object_sha256':sha(p)})
  folds[fold]={'path':str(p.relative_to(ROOT)),'sha256':sha(p)};print('FITTED B',fold,kind,win,flush=True)
 dump(B/'C8B_SPEC_LOCK.json',{'timestamp':now(),'source_spec_sha256':sha(A/'C8A_SPEC_LOCK.json'),'A_gate_sha256':sha(A/'C8A_STATUS.json'),'fold_objects':folds,'no_outer_patient_labels_used_for_fitting':True});(B/'C8B_SPEC_LOCK.sha256').write_text(sha(B/'C8B_SPEC_LOCK.json')+'  C8B_SPEC_LOCK.json\n')

def control_inputs(v,model,control,seed):
 vv={k:(a.copy() if isinstance(a,np.ndarray) else a) for k,a in v.items()};rng=np.random.default_rng(941+seed);idx=rng.permutation(len(v['st']));mo=dict(model)
 if control=='NC3':vv['image']=v['image'][idx]
 elif control=='NC2':vv['distance']=v['distance'][idx];vv['features'][:,:3]=v['features'][idx,:3]
 elif control=='NC1':
  # A lookup gate exactly reproduces a shuffled empirical weight distribution.
  w=weights(v,model)[idx];vv['random_weight_override']=w;mo={'family':'explicit_control','weights':w.tolist()}
 return vv,mo
# Explicit control weights are input-only, not a fitted candidate.
_original_predict_B=predict_B
def predict_B(v,model):
 if model['family']=='explicit_control':return v['st']+np.array(model['weights'])[:,None]*(v['image']-v['st'])
 return _original_predict_B(v,model)

def evaluate_B(stress=False):
 check_lock();set_forbidden([]);lock=json.loads((B/'C8B_SPEC_LOCK.json').read_text());stage='C8C' if stress else 'C8B';path=OUT/stage
 if stress:assert json.loads((B/'C8B_STATUS.json').read_text())['status']=='PASS_FUSION_CALIBRATION'
 dump(path/'logs/OUTER_ACCESS.json',{'timestamp':now(),'fit_lock_sha256':sha(B/'C8B_SPEC_LOCK.json'),'stress':stress})
 n=0
 for fold,held in cfg()['folds'].items():
  op=ROOT/lock['fold_objects'][fold]['path'];assert sha(op)==lock['fold_objects'][fold]['sha256'];cal=json.loads(op.read_text());obj=json.loads((A/'scalers'/f'{fold}_{cal["transform"]}.json').read_text())
  for sample in held:
   with allowed([sample]):y=target(sample)
   for r in masks(sample,stress=stress):
    if stress and (r['regime']!='V2' or r['budget']==.2):continue
    v=inputs(r,obj,y)
    models={'F0':cal['family_models']['F0'],'SELECTED':cal['selected_model']} if stress else {**cal['family_models'],'SELECTED':cal['selected_model']}
    for name,model in models.items():
     pred=predict_B(v,model);save_outer(stage,fold,r,obj,name,pred,v,y,op,model_obj=model);n+=1
     if not stress and r['budget']==.2 and r['seed']==0 and r['geometry']==(0 if r['regime']=='V1' else 1) and name in ['F1','F2','F3','SELECTED']:
      for control in ['NC1','NC2','NC3']:
       for seed in range(3):vv,mo=control_inputs(v,model,control,seed);pred=predict_B(vv,mo);save_outer(stage,fold,r,obj,name,pred,v,y,op,control,seed,model)
  print('EVALUATED',stage,fold,n,flush=True)
 assert n==(960 if stress else 2400)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['fit_A','evaluate_A','fit_B','evaluate_B','evaluate_C']);args=parser.parse_args()
 with threadpool_limits(4):globals()['evaluate_B' if args.stage=='evaluate_C' else args.stage](**({'stress':True} if args.stage=='evaluate_C' else {}))
