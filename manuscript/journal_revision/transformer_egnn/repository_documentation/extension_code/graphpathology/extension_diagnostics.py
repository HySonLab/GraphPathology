"""Conditional geometric, attention-reach and distant-feature diagnostics; no biological claims."""
import copy,json,time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from .extension_models import ExtensionModel,VARIANTS,graph_spec,model_inputs,SCALE_UM
from .extension_evaluate import model_inventory,load_data
from .extension_train import metadata
from .assets import sha

def permuted(spec,perm):
 out=copy.deepcopy(spec);inv=np.empty(len(perm),int);inv[perm]=np.arange(len(perm))
 for key in ['p','degree']:out[key]=spec[key][perm]
 out['row']=torch.tensor(inv[spec['row'].numpy()],dtype=torch.long);out['col']=torch.tensor(inv[spec['col'].numpy()],dtype=torch.long)
 out['landmarks']=torch.tensor(inv[spec['landmarks'].numpy()],dtype=torch.long);out['local_links']=torch.tensor(inv[spec['local_links'].numpy()[perm]],dtype=torch.long)
 out['adjacency']=spec['adjacency'][perm,:][:,perm].tocsr();out['physical_xy']=spec['physical_xy'][perm]
 return out

def geometric_test(item,root):
 seed=item['seed'];variant=item['variant'];torch.manual_seed(19);rng=np.random.default_rng(808)
 xy=rng.normal(size=(41,2))*80;x=rng.normal(size=(41,512)).astype(np.float32)
 spec=graph_spec(xy,'diagnostic',seed,variant);saved=spec['p'].clone();physical=xy.copy();model=ExtensionModel(variant)
 model.load_state_dict(torch.load(root/item['checkpoint'],weights_only=True)['model']);model.eval()
 records=[];angle=.731;rotation=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
 for precision in ['float64','float32']:
  net=copy.deepcopy(model).double() if precision=='float64' else copy.deepcopy(model).float()
  # Fixed-coordinate graph baselines aggregate in float32; their scalar diagnostics use float32 bounds.
  tolerance=1e-8 if precision=='float64' and net.kind in ['gps','egnn'] else 1e-4
  z=torch.tensor(model_inputs(x,spec,variant),dtype=torch.float64 if precision=='float64' else torch.float32)
  with torch.inference_mode():base,latent,_=net(z,spec,True)
  for name,R in [('rotation',rotation),('reflection',rotation@np.diag([-1.,1.])),('translation',np.eye(2))]:
   t=copy.deepcopy(spec);offset=np.array([391.2,-172.7]);t['p']=spec['p']@torch.tensor(R.T);t['center_um']=spec['center_um']@torch.tensor(R.T)+torch.tensor(offset)
   t['physical_xy']=xy@R.T+offset
   with torch.inference_mode():yp,lp,_=net(z,t,True)
   scalar=float((base-yp).abs().max());old=latent.numpy()*SCALE_UM+spec['center_um'].numpy();new=lp.numpy()*SCALE_UM+t['center_um'].numpy();coordinate=float(np.max(abs(new-(old@R.T+offset))))
   np.testing.assert_allclose(yp.numpy(),base.numpy(),rtol=tolerance,atol=tolerance)
   # Physical latent coordinates inherit the scale factor; explicit tolerance is reported in um.
   coordinate_atol=tolerance*SCALE_UM*4
   np.testing.assert_allclose(new,old@R.T+offset,rtol=tolerance,atol=coordinate_atol)
   records.append({'model_id':item['id'],'variant':variant,'seed':seed,'precision':precision,'transform':name,'max_scalar_abs_delta':scalar,'max_latent_um_abs_delta':coordinate,'rtol':tolerance,'scalar_atol':tolerance,'latent_um_atol':coordinate_atol,'status':'COMPLETE'})
  perm=rng.permutation(len(x));t=permuted(spec,perm);v=torch.tensor(model_inputs(x[perm],t,variant),dtype=z.dtype)
  with torch.inference_mode():yp=net(v,t)
  delta=float((yp-base[perm]).abs().max());torch.testing.assert_close(yp,base[perm],rtol=tolerance,atol=tolerance)
  records.append({'model_id':item['id'],'variant':variant,'seed':seed,'precision':precision,'transform':'node_permutation','max_scalar_abs_delta':delta,'rtol':tolerance,'scalar_atol':tolerance,'status':'COMPLETE'})
 torch.testing.assert_close(spec['p'],saved,rtol=0,atol=0);np.testing.assert_array_equal(xy,physical)
 return records

def attention(root,patient,item):
 d=load_data(root,patient);spec=graph_spec(d['xy_um'],patient,item['seed'],item['variant']);net=ExtensionModel(item['variant']);net.load_state_dict(torch.load(root/item['checkpoint'],weights_only=True)['model']);net.eval()
 with torch.inference_mode():_,_,weights=net(torch.tensor(d['x']),spec,True)
 if net.attention=='global':ids=spec['landmarks'].numpy();physical=np.linalg.norm(d['xy_um'][:,None]-d['xy_um'][ids][None],axis=-1);assigned=np.linalg.norm(spec['p'].numpy()[:,None]-spec['p'].numpy()[ids][None],axis=-1)*SCALE_UM
 else:
  ids=spec['local_links'].numpy() if net.attention=='local' else np.arange(len(d['x']))[:,None];physical=np.linalg.norm(d['xy_um'][:,None]-d['xy_um'][ids],axis=-1);assigned=np.linalg.norm(spec['p'].numpy()[:,None]-spec['p'].numpy()[ids],axis=-1)*SCALE_UM
 records=[]
 for layer,w in enumerate(weights):
  a=w.numpy()
  for head in range(a.shape[1]):
   for threshold in [250,600,1200]:
    far=physical>threshold;far_a=assigned>threshold
    records.append({'patient':patient,'model_id':item['id'],'variant':item['variant'],'seed':item['seed'],'layer':layer+1,'head':head,'threshold_um':threshold,
     'n_nodes':len(physical),'key_slots':physical.shape[1],'n_connections':physical.size,'physical_fraction_beyond':float(far.mean()),'assigned_fraction_beyond':float(far_a.mean()),'attention_mass_beyond':float((a[:,head,:]*far).sum(1).mean()),
     'note':'Connectivity/attention is not evidence of molecular correlation','graph_control':json.dumps(spec['info'],sort_keys=True)})
 return records

def remote_test(root,items):
 patient='TENX99';d=load_data(root,patient);full=next(i for i in items if i['variant']=='gps_full' and i['seed']==0);s=graph_spec(d['xy_um'],patient,0,'gps_full');landmarks=s['landmarks'].numpy()
 distances=np.linalg.norm(d['xy_um'][:,None]-d['xy_um'][landmarks][None],axis=-1);query,j=np.unravel_index(np.argmax(distances),distances.shape);remote=int(landmarks[j]);distance=float(distances[query,j]);rows=[]
 for variant in ['gps_full','gps_local','egnn_latent']:
  item=next(i for i in items if i['variant']==variant and i['seed']==0);spec=graph_spec(d['xy_um'],patient,0,variant);model=ExtensionModel(variant);model.load_state_dict(torch.load(root/item['checkpoint'],weights_only=True)['model']);model.eval()
  x=torch.tensor(d['x']);altered=x.clone();altered[remote]=0
  with torch.inference_mode():a=model(x,spec)[query];b=model(altered,spec)[query]
  rms=float(torch.mean((a-b)**2).sqrt());maximum=float((a-b).abs().max())
  if variant in ['gps_local','egnn_latent'] and distance>1200:torch.testing.assert_close(a,b,rtol=0,atol=0)
  rows.append({'patient':patient,'variant':variant,'query_index':int(query),'perturbed_landmark_index':remote,'distance_um':distance,'rms_prediction_change':rms,'max_abs_prediction_change':maximum,'status':'COMPLETE' if distance>1200 else 'BLOCKED_INPUTS','interpretation':'distant feature sensitivity, not biological correlation'})
 return rows

def run(root):
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True);out=root/'experiments/journal_revision/transformer_egnn/diagnostics';out.mkdir(exist_ok=True)
 items=[i for i in model_inventory(root) if i['origin']=='new'];assert len(items)==48
 rows=[]
 for item in items:rows.extend(geometric_test(item,root))
 pd.DataFrame(rows).to_csv(out/'GEOMETRIC_INVARIANCE.csv',index=False)
 reaches=[]
 for patient in ['TENX99','1160920F','CID4290','CID44971']:
  for item in items:
   if item['variant'].startswith('gps_'):reaches.extend(attention(root,patient,item))
  print('ATTENTION DIAGNOSTICS',patient,flush=True)
 pd.DataFrame(reaches).to_csv(out/'ATTENTION_REACH.csv',index=False);pd.DataFrame(remote_test(root,items)).to_csv(out/'DISTANT_FEATURE_PERTURBATION.csv',index=False)
 r=metadata(root,status='COMPLETE',geometric_checks=len(rows),attention_rows=len(reaches),image_pipeline_equivariance='SKIPPED_PROTOCOL; no end-to-end claim',output_hashes={p.name:sha(p) for p in out.glob('*.csv')})
 (out/'STATUS.json').write_text(json.dumps(r,indent=2)+'\n');return r


def recovered_patient(root,patient):
 """Additional recovered-cohort diagnostics, preserving the first completed bundle."""
 out=root/'experiments/journal_revision/transformer_egnn/diagnostics/recovered'/patient;out.mkdir(parents=True,exist_ok=True);done=out/'COMPLETE.json'
 if done.exists():
  r=json.loads(done.read_text())
  for name,digest in r['output_hashes'].items():assert sha(out/name)==digest
  return r
 torch.set_num_threads(1);items=[i for i in model_inventory(root) if i['origin']=='new'];rows=[];latent=[]
 for item in items:
  if item['variant'].startswith('gps_'):rows.extend(attention(root,patient,item))
  if VARIANTS[item['variant']][-1]:
   d=load_data(root,patient);spec=graph_spec(d['xy_um'],patient,item['seed'],item['variant']);torch.manual_seed(item['seed']);initial=ExtensionModel(item['variant']);net=ExtensionModel(item['variant']);net.load_state_dict(torch.load(root/item['checkpoint'],weights_only=True)['model']);net.eval()
   with torch.inference_mode():_,p,_=net(torch.tensor(d['x']),spec,True)
   shift=np.linalg.norm((p.numpy()-spec['p'].numpy())*SCALE_UM,axis=1)
   latent.append({'patient':patient,'variant':item['variant'],'seed':item['seed'],'mean_latent_displacement_um':float(shift.mean()),'max_latent_displacement_um':float(shift.max()),'first_coordinate_weight_change_norm':float(torch.linalg.vector_norm(net.blocks[0].coordinate.weight-initial.blocks[0].coordinate.weight).detach()),'last_coordinate_weight_change_norm':float(torch.linalg.vector_norm(net.blocks[1].coordinate.weight-initial.blocks[1].coordinate.weight).detach()),'last_update_affects_scalar_readout':False})
 pd.DataFrame(rows).to_csv(out/'ATTENTION_REACH.csv',index=False);pd.DataFrame(latent).to_csv(out/'LATENT_UPDATES.csv',index=False)
 r=metadata(root,status='COMPLETE',threads=1,patient=patient,output_hashes={p.name:sha(p) for p in out.glob('*.csv')});done.write_text(json.dumps(r,indent=2)+'\n');print('RECOVERED DIAGNOSTICS',patient,flush=True);return r
