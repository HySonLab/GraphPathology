"""Versioned full-section extension fitting; existing journal-v2 fits are read-only."""
import copy,datetime,json,time,subprocess
from pathlib import Path
import numpy as np
import torch
from .assets import sha
from .core import Scaler,PeakRSS
from .journal_learned import load_prepared
from .extension_models import VARIANTS,ExtensionModel,graph_spec,model_inputs

PROTOCOL_COMMIT='2e03fae636d6f6527d25db58f4cbab72f7c1a286'

def metadata(root,**extra):
 return {'label':'post_v1_exploratory','protocol_commit':PROTOCOL_COMMIT,'code_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
 'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'code_hashes':{str(p.relative_to(root)):sha(p) for p in (root/'src/graphpathology').glob('extension*.py')},
 'environment':json.loads((root/'experiments/journal_revision/ENVIRONMENT.json').read_text()),'device':'cpu','threads':4,**extra}

def dataset(root,samples):
 source=root/'data/journal_revision/hest_training';data=load_prepared(source,samples)
 with np.load(root/'checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz') as f:mean=f['mean'];scale=f['scale']
 for sample,d in data.items():d['x']=((d['embeddings']-mean)/scale).astype(np.float32)
 return data

def prepare_variant(data,variant,seed):
 specs={s:graph_spec(d['xy_um'],s,seed,variant) for s,d in data.items()}
 xs={s:torch.from_numpy(model_inputs(d['x'],specs[s],variant)) for s,d in data.items()}
 return specs,xs

def preflight(root):
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 data=dataset(root,['NCBI783','NCBI785','TENX95']);records=[]
 for variant in VARIANTS:
  specs,xs=prepare_variant(data,variant,0);torch.manual_seed(0);model=ExtensionModel(variant)
  # Reduced real-data source-only backward; never evaluate TENX95 expression in preflight.
  small={k:(v[:256] if isinstance(v,np.ndarray) and len(v)==len(data['NCBI783']['x']) else v) for k,v in data['NCBI783'].items()}
  ss=graph_spec(small['xy_um'],'source_smoke',0,variant);xx=torch.tensor(model_inputs(small['x'],ss,variant))
  loss=(model(xx,ss)-torch.tensor(small['targets'])) .square().mean();loss.backward();model.zero_grad(set_to_none=True)
  start=time.monotonic()
  with PeakRSS() as memory:
   for sample in ['NCBI783','NCBI785']:
    loss=(model(xs[sample],specs[sample])-torch.tensor(data[sample]['targets'])).square().mean();loss.backward();model.zero_grad(set_to_none=True)
   model.eval()
   with torch.inference_mode():model(xs['TENX95'],specs['TENX95'])
  elapsed=time.monotonic()-start
  records.append({'variant':variant,'parameters':sum(p.numel() for p in model.parameters()),'source_backward_and_unscored_validation_forward_seconds':elapsed,'peak_rss_mb':memory.peak/2**20,
                  'per_trial80_epoch_ceiling_seconds':elapsed*80,'geometry_controls':{s:specs[s]['info'] for s in specs}})
  print('PREFLIGHT',variant,round(elapsed,3),'seconds/epoch',flush=True)
 estimate=sum(x['per_trial80_epoch_ceiling_seconds'] for x in records)*3
 epochs=80 if estimate<=6*3600 else 40
 status='COMPLETE' if estimate*epochs/80<=6*3600 and max(x['peak_rss_mb'] for x in records)<8192 else 'BLOCKED_RESOURCES'
 result=metadata(root,status=status,records=records,estimated_total80_epoch_hours=estimate/3600,common_max_epochs=epochs,resource_rule='all16 fresh configurations use the same selected ceiling; no validation losses inspected')
 out=root/'experiments/journal_revision/transformer_egnn/PREFLIGHT.json'
 if out.exists():raise FileExistsError(out)
 out.write_text(json.dumps(result,indent=2)+'\n');return result

def train(root,variant,seed):
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 cfg=json.loads((root/'configs/journal_revision/transformer_egnn/protocol.json').read_text());pf=json.loads((root/'experiments/journal_revision/transformer_egnn/PREFLIGHT.json').read_text())
 if pf['status']!='COMPLETE':raise RuntimeError('Preflight resource gate failed')
 if pf['common_max_epochs']!=80 and not (root/'experiments/journal_revision/transformer_egnn/RESOURCE_ADJUSTMENT.json').exists():raise RuntimeError('Commit common resource adjustment before fitting')
 ident=f'{variant}_seed{seed}';out=root/'experiments/journal_revision/transformer_egnn/fits';out.mkdir(exist_ok=True)
 cp=root/'checkpoints/journal_revision/transformer_egnn'/ident;cp.mkdir(parents=True,exist_ok=True);record_path=out/(ident+'.json')
 if record_path.exists():
  r=json.loads(record_path.read_text());assert sha(cp/'best.pt')==r['best_checkpoint_sha256'];return r
 data=dataset(root,['NCBI783','NCBI785','TENX95']);specs,xs=prepare_variant(data,variant,seed);ys={s:torch.tensor(d['targets']) for s,d in data.items()}
 torch.manual_seed(seed);model=ExtensionModel(variant);optimizer=torch.optim.AdamW(model.parameters(),lr=cfg['lr'],weight_decay=cfg['weight_decay'])
 begin=time.monotonic();best=float('inf');best_epoch=0;stale=0;history=[]
 resume=cp/'resume.pt';start_epoch=1
 if resume.exists():
  saved=torch.load(resume,map_location='cpu',weights_only=True);model.load_state_dict(saved['model']);optimizer.load_state_dict(saved['optimizer']);torch.set_rng_state(saved['rng'])
  start_epoch=saved['epoch']+1;best=saved['best'];best_epoch=saved['best_epoch'];stale=saved['stale'];history=saved['history']
 with PeakRSS() as memory:
  for epoch in range(start_epoch,pf['common_max_epochs']+1):
   model.train();losses=[]
   for sample in np.random.default_rng(seed*10000+epoch).permutation(['NCBI783','NCBI785']):
    optimizer.zero_grad(set_to_none=True);prediction=model(xs[sample],specs[sample]);loss=(prediction-ys[sample]).square().mean()
    if not torch.isfinite(loss):raise FloatingPointError('Nonfinite extension fit loss')
    loss.backward();optimizer.step();losses.append(float(loss.detach()))
   model.eval()
   with torch.inference_mode():val=float((model(xs['TENX95'],specs['TENX95'])-ys['TENX95']).square().mean())
   if not np.isfinite(val):raise FloatingPointError('Nonfinite validation loss')
   history.append({'epoch':epoch,'training_mse':float(np.mean(losses)),'validation_mse':val})
   if val<best-cfg['min_delta']:
    best=val;best_epoch=epoch;stale=0;torch.save({'model':model.state_dict(),'optimizer':optimizer.state_dict(),'epoch':epoch,'variant':variant,'seed':seed,'config':cfg},cp/'best.pt')
   else:stale+=1
   torch.save({'model':model.state_dict(),'optimizer':optimizer.state_dict(),'epoch':epoch,'rng':torch.get_rng_state(),'best':best,'best_epoch':best_epoch,'stale':stale,'history':history},cp/'resume.tmp');(cp/'resume.tmp').replace(resume)
   if epoch%10==0:print('FIT',ident,epoch,'best_epoch',best_epoch,flush=True)
   if stale>=cfg['patience']:break
 result=metadata(root,status='COMPLETE',variant=variant,seed=seed,parameters=sum(p.numel() for p in model.parameters()),selected_epoch=best_epoch,epochs_ran=len(history),validation_mse=best,
  elapsed_seconds=time.monotonic()-begin,peak_rss_mb=memory.peak/2**20,best_checkpoint_relative=str((cp/'best.pt').relative_to(root)),best_checkpoint_sha256=sha(cp/'best.pt'),
  input_cache_hashes={s:sha(root/'data/journal_revision/hest_training/processed/journal_v2_hest'/(s+'.npz')) for s in data},feature_scaler_sha256=sha(root/'checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz'),
  graph_controls={s:specs[s]['info'] for s in specs},training=cfg,history=history,optimizer_state_included=True,selection_patient='TENX95',fitting_patients=['NCBI783','NCBI785'])
 record_path.write_text(json.dumps(result,indent=2)+'\n');print('FIT COMPLETE',ident,best_epoch,flush=True);return result
