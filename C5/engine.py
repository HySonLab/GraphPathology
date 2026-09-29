"""Sparse-ST prediction contract: scene has no labels; anchors contain observed rows only."""
from common import *
import torch
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from graph import operator
SCENES={};IMAGES={};PRED={};GRAPHS={};PCA_MODELS={}
def training(cohort):return SPLITS['splits']['train'] if cohort=='primary' else [s for s in IDS if s!=cohort]
def image_data(encoder,sample):
 key=(encoder,sample)
 if key not in IMAGES:
  path=ROOT/'C4/embeddings'/encoder/f'{sample}.npz'
  with np.load(path) as f:IMAGES[key]={'x':f['embeddings'].copy(),'barcodes':f['barcodes'].copy(),'xy':f['xy_um'].copy(),'feature_cache_hash':sha(path)}
 return IMAGES[key]
def pca_model(cohort,encoder):
 key=(cohort,encoder)
 if key in PCA_MODELS:return PCA_MODELS[key]
 path=OUT/'models'/f'PCA_{cohort}_{encoder}.npz'
 if not path.exists():
  x=np.concatenate([image_data(encoder,s)['x'] for s in training(cohort)]);sc=core.Scaler().fit(x);z=sc.transform(x);p=PCA(n_components=8,svd_solver='randomized',random_state=0,whiten=True).fit(z)
  np.savez_compressed(path,mean=sc.mean,scale=sc.scale,components=p.components_,pca_mean=p.mean_,variance=p.explained_variance_)
  registry('fit',rid('PCA',key),{'model':'image_PCA','encoder':encoder,'train_slides':training(cohort),'checkpoint_path':str(path.relative_to(ROOT)),'checkpoint_hash':sha(path),'expression_access':False})
 with np.load(path) as f:PCA_MODELS[key]={k:f[k].copy() for k in f.files}
 return PCA_MODELS[key]
def image_prediction(cohort,encoder,panel,sample):
 key=(cohort,encoder,panel,sample)
 if key in PRED:return PRED[key]
 path=OUT/'cache'/('_'.join(key)+'_B1.npz');meta=path.with_suffix('.json')
 if path.exists():
  r=json.loads(meta.read_text());assert sha(path)==r['sha256']
  with np.load(path) as f:y=f['prediction'].copy()
 else:
  rs=[]
  for p in (ROOT/'C4/C4A/logs/runs').glob('*.json'):
   r=json.loads(p.read_text())
   if r['encoder']==encoder and r['panel']==panel and r['model']=='D1' and ((cohort=='primary' and r['phase']=='primary') or (cohort!='primary' and r['phase']=='lopo' and r['heldout']==cohort)):rs.append(r)
  assert len(rs)==3;pred=[];refs=[]
  for r in sorted(rs,key=lambda r:r['seed']):
   ev=ROOT/'C4/C4A/logs/evaluation'/f'{r["run_id"]}_{sample}.json'
   if ev.exists():
    e=json.loads(ev.read_text());assert sha(ROOT/e['prediction_path'])==e['prediction_sha256']
    with np.load(ROOT/e['prediction_path']) as f:z=f['prediction'].copy();assert np.array_equal(f['barcodes'],image_data(encoder,sample)['barcodes'])
    mode='reused saved prediction'
   else:
    assert sha(ROOT/r['model_path'])==r['model_sha256'];net=core.mlp_model(r['input_dim'],r['output_dim'],r['params']);net.load_state_dict(torch.load(ROOT/r['model_path'],map_location='cpu',weights_only=True));net.eval()
    with np.load(ROOT/r['preprocessing']/'feature.npz') as f:x=((image_data(encoder,sample)['x']-f['mean'])/f['scale']).astype(np.float32)
    with np.load(ROOT/r['preprocessing']/'target.npz') as f:z=core.predict(net,x,'B2').astype(float)*f['scale']+f['mean']
    mode='read-only frozen-head inference for C5 training/validation'
   pred.append(z);refs.append({'run_id':r['run_id'],'model_path':r['model_path'],'model_hash':r['model_sha256'],'mode':mode,'preprocessing':r['preprocessing']})
  y=np.mean(pred,axis=0).astype(np.float32);np.savez_compressed(path,prediction=y);dump(meta,{'sha256':sha(path),'parents':refs,'source_commit':commit(),'target_labels_accessed':False})
 PRED[key]=y;return y

def scene(cohort,encoder,panel,sample):
 key=(cohort,encoder,panel,sample)
 if key not in SCENES:
  d=image_data(encoder,sample);pc=pca_model(cohort,encoder);z=((d['x']-pc['mean'])/pc['scale']).astype(np.float32);z=((z-pc['pca_mean'])@pc['components'].T/np.sqrt(pc['variance'])).astype(np.float32)
  SCENES[key]={'xy':d['xy'],'barcodes':d['barcodes'],'image':image_prediction(cohort,encoder,panel,sample),'projection':z,'feature_cache_hash':d['feature_cache_hash'],'sample':sample}
 return SCENES[key]
def graph(sample):
 if sample not in GRAPHS:GRAPHS[sample]=operator(sparse.load_npz(ROOT/'C2/graphs'/f'{sample}_knn6_cap250.npz'),'MEAN').tocsr()
 return GRAPHS[sample]
def observed_payload(y,observed,xy):
 assert observed.dtype==bool and len(y)==len(observed)==len(xy)
 ids=np.flatnonzero(observed);values=np.array(y[ids],dtype=np.float32,copy=True);assert np.isfinite(values).all()
 return {'indices':ids,'xy':np.array(xy[ids],copy=True),'values':values}
def anchors(sample,panel,r):
 # Full labels are immediately reduced to observed rows at this I/O boundary.
 return observed_payload(labels(sample,panel),mask_array(r),image_data(ENCODERS[0],sample)['xy'])
def nearest(s,a):
 k=min(16,len(a['indices']));dist,index=cKDTree(a['xy']).query(s['xy'],k=k)
 if k==1:dist=dist[:,None];index=index[:,None]
 return np.round(dist,6),index

def interpolate(a,dist,index,family,param):
 if family=='B2a':return a['values'][index[:,:int(param)]].mean(1)
 w=1/np.maximum(dist,1e-6)**float(param);w/=w.sum(1,keepdims=True);return np.einsum('nk,nkg->ng',w,a['values'][index],optimize=True).astype(np.float32)
def molecular(s,a):
 n=len(s['xy']);o=np.zeros(n,np.float32);o[a['indices']]=1;z=np.zeros((n,a['values'].shape[1]),np.float32);z[a['indices']]=a['values'];op=graph(s['sample']);m1=np.asarray(op@o);m2=np.asarray(op@m1);v1=np.asarray(op@z);v2=np.asarray(op@v1);v1/=np.maximum(m1[:,None],1e-6);v2/=np.maximum(m2[:,None],1e-6)
 return o,z,v1.astype(np.float32),v2.astype(np.float32),m1,m2

def fusion_spec(selection,panel,regime,encoder):return selection[panel][regime],selection[panel][regime]['lambda'][encoder]
def feature_parts(s,a,selected,lam):
 dist,index=nearest(s,a);interp=interpolate(a,dist,index,selected['best_family'],selected['best_param']);fuse=(lam*s['image']+(1-lam)*interp).astype(np.float32)
 o,z,v1,v2,m1,m2=molecular(s,a);op=graph(s['sample']);x=s['projection'];x1=np.asarray(op@x,dtype=np.float32);x2=np.asarray(op@x1,dtype=np.float32)
 shared=np.column_stack([x,x1,x2,o,m1,m2,np.log1p(dist[:,0]/100),np.sum(dist<=300,axis=1)/16]).astype(np.float32)
 gene=np.stack([s['image'],fuse,z,v1,v2],axis=2);return shared,gene,fuse,dist[:,0]
def design(shared,gene,g):return np.column_stack([shared,gene[:,g,:],np.ones(len(shared))]).astype(np.float64)
def b4_prediction(parts,weights):
 shared,gene,fuse,_=parts;y=np.empty_like(fuse)
 for g in range(fuse.shape[1]):y[:,g]=fuse[:,g]+design(shared,gene,g)@weights[g]
 return y

def predict(s,a,hidden,model,param,selected=None,lam=None,weights=None,mean=None):
 # No full expression matrix or hidden target is accepted by this function.
 if model=='B0':return np.broadcast_to(mean,(len(hidden),len(mean))).astype(np.float32).copy()
 if model=='B1':return s['image'][hidden].copy()
 if model in ['B2a','B2b']:
  d,i=nearest(s,a);return interpolate(a,d,i,model,param)[hidden]
 if model=='B3':
  d,i=nearest(s,a);p=interpolate(a,d,i,selected['best_family'],selected['best_param']);return (lam*s['image'][hidden]+(1-lam)*p[hidden]).astype(np.float32)
 if model=='B4':return b4_prediction(feature_parts(s,a,selected,lam),weights)[hidden]
 raise ValueError(model)
def control_inputs(s,a,name,seed):
 rng=np.random.default_rng(seed);s={**s};a={k:v.copy() for k,v in a.items()}
 if name=='NC1':a['values']=a['values'][rng.permutation(len(a['values']))]
 elif name=='NC2':p=rng.permutation(len(a['indices']));a['indices']=a['indices'][p];a['xy']=a['xy'][p]
 elif name=='NC3':p=rng.permutation(len(s['xy']));s['image']=s['image'][p];s['projection']=s['projection'][p]
 return s,a
