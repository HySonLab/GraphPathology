"""Matched validation-only Ridge/MLP benchmark on two immutable embedding caches."""
from common import *
import argparse,time
import torch
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits
CFG=json.loads((A/'configs/protocol.json').read_text())

def load(encoder,sample):
 path=OUT/'embeddings'/encoder/f'{sample}.npz'
 with np.load(path) as f:x=f['embeddings'].copy();barcodes=f['barcodes'].copy()
 with np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as f:
  assert np.array_equal(barcodes,f['barcodes']) and list(f['genes'])==PANELS['S2-280'];y=f['normalized_targets'].copy()
 assert np.isfinite(x).all();return {'x':x,'y':y,'barcodes':barcodes,'cache_sha256':sha(path)}

def prep(encoder,panel,data,heldout=None):
 folder=A/'preprocessing'/encoder/(heldout or 'primary')/panel
 if not folder.exists():
  folder.mkdir(parents=True);train=[s for s in IDS if s!=heldout] if heldout else SPLITS['splits']['train']
  core.Scaler().fit(np.concatenate([data[s]['x'] for s in train])).save(folder/'feature.npz');core.Scaler().fit(np.concatenate([data[s]['y'][:,ix(panel)] for s in train])).save(folder/'target.npz')
 return folder

def fit(encoder,panel,model,params,seed,data,phase,heldout=None,fixed_epochs=None,control=None):
 cfgname=('alpha'+str(params['alpha'])) if model=='D0' else f'h{params["hidden"]}d{params["dropout"]}'
 run_id=f'C4A_{phase}_{encoder}_{heldout or "primary"}_{panel}_{model}_{cfgname}_s{seed}'+('_'+control if control else '')
 path=A/'logs/runs'/f'{run_id}.json'
 if path.exists():return json.loads(path.read_text())
 train=[s for s in IDS if s!=heldout] if heldout else SPLITS['splits']['train'];val=[] if heldout else SPLITS['splits']['validation'];folder=prep(encoder,panel,data,heldout);xs=scaler(folder/'feature.npz');ys=scaler(folder/'target.npz')
 x=xs.transform(np.concatenate([data[s]['x'] for s in train]));y=ys.transform(np.concatenate([data[s]['y'][:,ix(panel)] for s in train]));xv=xs.transform(np.concatenate([data[s]['x'] for s in val])) if val else None;yv=ys.transform(np.concatenate([data[s]['y'][:,ix(panel)] for s in val])) if val else None
 permhash=None
 if control:
  x,y,perm=core.permute_training(x,y,control,seed);pp=A/'models'/f'{run_id}_permutation.npy';np.save(pp,perm);permhash=sha(pp)
 with core.PeakRSS() as mem:
  if model=='D0':
   t=time.perf_counter();net=Ridge(alpha=params['alpha'],solver='cholesky').fit(x,y);info={'train_seconds':time.perf_counter()-t,'best_epoch':0,'parameters':int(net.coef_.size+net.intercept_.size),'validation_seconds':0};mp=A/'models'/f'{run_id}.npz';np.savez_compressed(mp,coef=net.coef_,intercept=net.intercept_)
  else:
   net,info=core.train_mlp(x,y,xv,yv,params,base.C1CFG['mlp'],seed,fixed_epochs);history=info.pop('history');dump(A/'logs/epochs'/f'{run_id}.json',history);mp=A/'models'/f'{run_id}.pt';torch.save(net.state_dict(),mp)
  if val:
   pred=ys.inverse(core.predict(net,xv,'B1' if model=='D0' else 'B2'));m=core.gene_metrics(np.concatenate([data[s]['y'][:,ix(panel)] for s in val]),pred);summary=core.summarize(m);pd.DataFrame({'gene':PANELS[panel],**m}).to_csv(A/'logs/validation'/f'{run_id}.csv',index=False)
  else:summary={}
 r={'run_id':run_id,'encoder':encoder,'panel':panel,'model':model,'params':params,'seed':seed,'phase':phase,'heldout':heldout,'control':control,'permutation_sha256':permhash,'train_slides':train,'validation_slides':val,'input_dim':x.shape[1],'output_dim':y.shape[1],'n_train':len(x),'preprocessing':str(folder.relative_to(ROOT)),'feature_scaler_sha256':sha(folder/'feature.npz'),'target_scaler_sha256':sha(folder/'target.npz'),'model_path':str(mp.relative_to(ROOT)),'model_sha256':sha(mp),'cache_hashes':{s:data[s]['cache_sha256'] for s in train+val},'source_commit':commit(),'peak_rss_mb':mem.peak/1024**2,'validation_median_pearson':summary.get('median_pearson'),'validation_mse':summary.get('mean_mse'),**info}
 dump(path,r);registry('fit',run_id,{**r,'evaluation_slides':val or [heldout]});print('FIT',run_id,'val',r['validation_median_pearson'],flush=True);return r

def model_load(r):
 assert sha(ROOT/r['model_path'])==r['model_sha256']
 if r['model']=='D0':
  net=Ridge(alpha=r['params']['alpha'])
  with np.load(ROOT/r['model_path']) as f:net.coef_=f['coef'];net.intercept_=f['intercept']
  net.n_features_in_=r['input_dim']
 else:
  net=core.mlp_model(r['input_dim'],r['output_dim'],r['params']);net.load_state_dict(torch.load(ROOT/r['model_path'],map_location='cpu',weights_only=True));net.eval()
 return net

def predict(r,d):
 folder=ROOT/r['preprocessing'];x=scaler(folder/'feature.npz').transform(d['x']);return scaler(folder/'target.npz').inverse(core.predict(model_load(r),x,'B1' if r['model']=='D0' else 'B2'))

def evaluate(r,sample,d):
 path=A/'logs/evaluation'/f'{r["run_id"]}_{sample}.json'
 if path.exists():return
 t=time.perf_counter();p=predict(r,d);y=d['y'][:,ix(r['panel'])];m=core.gene_metrics(y,p);summary=core.summarize(m);co,_,_=core.coexpression(y,p);pp=A/'predictions'/f'{r["run_id"]}_{sample}.npz';np.savez_compressed(pp,prediction=p,genes=np.array(PANELS[r['panel']]),barcodes=d['barcodes'])
 cols={k:r[k] for k in ['run_id','encoder','panel','model','seed','phase','control']};cols['sample_id']=sample;gp=A/'logs/gene_metrics'/f'{r["run_id"]}_{sample}.csv';pd.DataFrame({'gene':PANELS[r['panel']],**m}).assign(**cols).to_csv(gp,index=False)
 out={**cols,**summary,**co,'prediction_path':str(pp.relative_to(ROOT)),'prediction_sha256':sha(pp),'gene_metrics_path':str(gp.relative_to(ROOT)),'inference_seconds':time.perf_counter()-t,'selection_lock_sha256':sha(A/'SELECTION_LOCK.json')};dump(path,out);registry('evaluation',r['run_id'],{**out,'evaluation_slides':[sample]});print('EVAL',r['run_id'],summary['median_pearson'],flush=True)

def select():
 assert not (A/'SELECTION_LOCK.json').exists();selected={};primary=[]
 for encoder in CFG['encoders']:
  data={s:load(encoder,s) for s in SPLITS['splits']['train']+SPLITS['splits']['validation']};selected[encoder]={}
  for panel in PANELS:
   selected[encoder][panel]={}
   for model in CFG['models']:
    grid=[{'alpha':a} for a in CFG['ridge_alphas']] if model=='D0' else CFG['mlp_candidates'];rs=[fit(encoder,panel,model,p,0,data,'selection') for p in grid];best=min(rs,key=lambda r:(-r['validation_median_pearson'],r['validation_mse'],json.dumps(r['params'],sort_keys=True)));epochs=[]
    for seed in ([0] if model=='D0' else CFG['training_seeds']):
     r=fit(encoder,panel,model,best['params'],seed,data,'primary');primary.append(r['run_id']);epochs.append(r['best_epoch'])
    selected[encoder][panel][model]={'params':best['params'],'lopo_epochs':int(np.floor(np.median(epochs)+.5))}
   if encoder=='pathology_encoder_1' and panel=='S2-50':
    for model in CFG['models']:
     for control in ['NC1','NC2']:
      for seed in CFG['training_seeds']:
       r=fit(encoder,panel,model,selected[encoder][panel][model]['params'],seed,data,'negative',control=control);primary.append(r['run_id'])
 lock={'timestamp':now(),'source_commit':commit(),'selected':selected,'primary_runs':primary,'spec_sha256':sha(A/'C4A_SPEC_LOCK.json'),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'run.py',OUT/'common.py',OUT/'extract.py',A/'configs/protocol.json']},'test_loaded':False};dump(A/'SELECTION_LOCK.json',lock);print('SELECTION LOCKED',flush=True)

def final():
 lock=json.loads((A/'SELECTION_LOCK.json').read_text())
 for p,h in lock['source_sha256'].items():assert sha(ROOT/p)==h
 if not (A/'logs/TEST_ACCESS.json').exists():dump(A/'logs/TEST_ACCESS.json',{'timestamp':now(),'selection_lock_sha256':sha(A/'SELECTION_LOCK.json')})
 for encoder in CFG['encoders']:
  data={s:load(encoder,s) for s in IDS}
  for run_id in lock['primary_runs']:
   r=json.loads((A/'logs/runs'/f'{run_id}.json').read_text())
   if r['encoder']==encoder:evaluate(r,SPLITS['splits']['test'][0],data[SPLITS['splits']['test'][0]])
  for heldout in IDS:
   for panel in PANELS:
    for model in CFG['models']:
     p=lock['selected'][encoder][panel][model]
     for seed in ([0] if model=='D0' else CFG['training_seeds']):
      r=fit(encoder,panel,model,p['params'],seed,data,'lopo',heldout,p['lopo_epochs']);evaluate(r,heldout,data[heldout])

def main():
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['select','final']);args=parser.parse_args()
 for folder in ['logs/runs','logs/epochs','logs/validation','logs/evaluation','logs/gene_metrics']:(A/folder).mkdir(exist_ok=True)
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):globals()[args.stage]()
if __name__=='__main__':main()
