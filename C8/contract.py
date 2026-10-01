"""Patient-scoped affine target harmonization and frozen image/ST inference."""
import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','4');os.environ.setdefault('MKL_NUM_THREADS','4')
import sys,json,csv,hashlib,time,subprocess,itertools
from pathlib import Path
from contextlib import contextmanager
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from scipy.special import expit
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C8';A=OUT/'C8A';B=OUT/'C8B'
sys.path.insert(0,str(ROOT/'C7'));import support as c7
core=c7.core;sha=c7.sha;dump=c7.dump;GENES=c7.panels()['S2-50'];SAMPLES=c7.SAMPLES;ALLOW=None;FORBIDDEN=set();SOURCE=['NCBI783','NCBI785'];CACHE={}
def now():return datetime.now(timezone.utc).isoformat()
def cfg():return json.loads((OUT/'configs/protocol.json').read_text())
def commit():return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
def rid(obj):return 'C8_'+hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()[:22]
def registry(event,ident,rec):
 p=OUT/'logs/registry'/(ident+'.json');assert not p.exists();r={**rec,'event':event,'run_id':ident,'checkpoint':'C8','source_commit':commit(),'timestamp':now()};dump(p,r)
 fields=next(csv.reader((ROOT/'results/experiment_registry.csv').open()));row={**r,'threads':4,'evaluation_slides':r.get('sample_id',''),'train_slides':'|'.join(r.get('train_patients',[]))}
 with (ROOT/'results/experiment_registry.csv').open('a',newline='') as f:csv.DictWriter(f,fieldnames=fields,extrasaction='ignore').writerow(row)
@contextmanager
def allowed(patients):
 global ALLOW
 old=ALLOW;ALLOW=frozenset(patients)
 try:yield
 finally:ALLOW=old
def target(sample):
 assert ALLOW is not None and sample in ALLOW and sample not in FORBIDDEN,('Forbidden label access',sample,ALLOW,FORBIDDEN)
 if sample in SOURCE:
  with np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as f:return f['normalized_targets'][:,[list(f['genes']).index(g) for g in GENES]].astype(float)
 return c7.truth(sample,'S2-50').astype(float)
def scene(sample):return c7.scene(sample,'S2-50')
def masks(sample,inner=False,stress=False):
 if 'mask_rows' not in CACHE:CACHE['mask_rows']=[json.loads(p.read_text()) for p in (ROOT/'C7/C7B/masks').glob('*.json')]
 rows=CACHE['mask_rows']
 for r in sorted(rows,key=lambda r:(r['regime'],r['geometry'],r['seed'])):
  if r['sample_id']!=sample or r['strategy']!='random' or (inner and r['seed']>2):continue
  for budget in ([.1,.2,.4] if r['regime']=='V1' or stress else [.2]):yield {**r,'budget':budget}
def arrays(r):
 with np.load(ROOT/r['mask_path']) as f:
  n=len(f['barcodes']);obs=np.zeros(n,bool);order=f['order'][:c7.count(n,r['budget'])];obs[order]=True;h=np.flatnonzero(~obs if r['regime']=='V1' else f['forbidden']);return obs,h,order

def moments(patients,robust=False):
 ys=[target(s) for s in patients];x=np.concatenate(ys);w=np.concatenate([np.full(len(y),1/len(y)/len(ys)) for y in ys]);mean=np.sum(x*w[:,None],axis=0);sd=np.sqrt(np.sum((x-mean)**2*w[:,None],axis=0));fallback=np.zeros(len(GENES),bool)
 if robust:
  qs=[]
  for q in [.25,.5,.75]:
   val=[]
   for j in range(x.shape[1]):
    ind=np.argsort(x[:,j],kind='stable');cw=np.cumsum(w[ind]);val.append(x[ind[min(np.searchsorted(cw,q),len(ind)-1)],j])
   qs.append(np.array(val))
  center=qs[1];scale=qs[2]-qs[0];fallback=scale<1e-8;scale=np.where(fallback,sd,scale)
 else:center=mean;scale=sd
 scale=np.where(scale<1e-8,1.,scale)
 return {'center':center.tolist(),'scale':scale.tolist(),'robust_SD_fallback':fallback.tolist(),'patient_weighting':'equal patients'}
def fit_scaler(train,kind):
 assert not set(train)&set(SOURCE);assert set(train)<=set(SAMPLES)
 if kind=='T0':src=dst={'center':[0.]*len(GENES),'scale':[1.]*len(GENES),'robust_SD_fallback':[False]*len(GENES)}
 else:
  with allowed(SOURCE):src=moments(SOURCE,kind=='T2')
  with allowed(train):dst=moments(train,kind=='T2')
 return {'transform':kind,'genes':GENES,'train_patients':list(train),'source_train_patients':SOURCE,'source':src,'destination':dst}
def convert(y,obj):return (np.asarray(y)-obj['destination']['center'])/obj['destination']['scale']
def inverse(y,obj):return np.asarray(y)*obj['destination']['scale']+obj['destination']['center']
def image_native(image,obj):return (np.asarray(image)-obj['source']['center'])/obj['source']['scale']
def inputs(r,obj,y):
 s=scene(r['sample_id']);obs,h,order=arrays(r);xy=s['xy'];values=convert(y[order],obj);k=min(16,len(order));dist,idx=cKDTree(xy[order]).query(xy,k=k)
 if k==1:dist=dist[:,None];idx=idx[:,None]
 dist=np.round(dist,6)
 if r['regime']=='V1':w=1/np.maximum(dist,1e-6)**2;w/=w.sum(1,keepdims=True);st=np.einsum('nk,nkg->ng',w,values[idx],optimize=True)
 else:st=values[idx].mean(1)
 # Use the original float32 C7 arithmetic at T0 for exact baseline replay.
 if obj['transform']=='T0':
  a=c7.c5.observed_payload(y.astype('float32'),obs,xy);st=c7.predict(s,a,h,'S2-50',r['regime'],'E1').astype(float)
 else:st=st[h]
 im=image_native(s['image'][h],obj);count=np.array(cKDTree(xy[order]).query_ball_point(xy[h],600,return_length=True));var=values[idx[h]].var(1).mean(1)
 feat=np.column_stack([np.log1p(dist[h,0]/200),np.log1p(dist[h].mean(1)/200),np.log1p(count),np.log1p(var)])
 return {'image':im,'st':st,'features':feat,'distance':dist[h,0],'hidden':h,'observed':order,'barcodes':s['barcodes'],'scale':np.array(obj['destination']['scale']),'center':np.array(obj['destination']['center'])}
def predict_A(v,obj,model,regime):
 if model=='E0':return v['image']
 if model=='E1':return v['st']
 lam=.25 if regime=='V1' else .5
 if obj['transform']=='T0':return (lam*v['image'].astype('float32')+(1-lam)*v['st'].astype('float32')).astype('float32').astype(float)
 return lam*v['image']+(1-lam)*v['st']
def weights(v,model):
 typ=model['family']
 if typ=='F0':return np.zeros(len(v['st']))
 if typ=='F1':return np.full(len(v['st']),model['lambda'])
 if typ=='F2':return np.array(model['weights'])[np.digitize(v['distance'],[300,600])]
 z=(v['features']-model['feature_mean'])/model['feature_scale'];return expit(model['coef'][0]+z@np.array(model['coef'][1:]))
def predict_B(v,model):
 w=weights(v,model);return v['st']+w[:,None]*(v['image']-v['st'])
def score(y,p,scale=None):
 met=core.gene_metrics(y,p);row=core.summarize(met)
 if scale is not None:
  row.update({'common_mse':float(np.mean(met['mse']*scale**2)),'common_mae':float(np.mean(met['mae']*scale))})
 return row,met

def select(rows,key):
 df=pd.DataFrame(rows);cols=['median_pearson','common_mse'];patient=df.groupby([key,'sample_id','regime','budget'],as_index=False)[cols].mean().groupby([key,'sample_id'],as_index=False)[cols].mean();summary=patient.groupby(key)[cols].mean();best=summary.median_pearson.max();eligible=summary[summary.median_pearson>=best-.005];chosen=eligible.sort_values(['common_mse'],kind='stable').index[0];return chosen,summary.reset_index().to_dict('records')
def fit_gate(records):
 # Equal patient/condition/run/location contribution, normalized by each F0 error.
 groups={}
 for r in records:groups.setdefault((r['sample_id'],r['regime'],r['budget']),[]).append(r)
 features=[];delta=[];error=[];mass=[]
 for rr in groups.values():
  for r in rr:
   v=r['v'];s=v['scale'];d=(v['image']-v['st'])*s;err=(v['st']-r['truth'])*s;den=max(float(np.mean(err**2)),1e-8);features.append(v['features']);delta.append(d);error.append(err);mass.append(np.full(len(d),1/(len(groups)*len(rr)*len(d)*den)))
 x=np.concatenate(features);d=np.concatenate(delta);e=np.concatenate(error);w=np.concatenate(mass);mean=np.average(x,axis=0,weights=w);sd=np.sqrt(np.average((x-mean)**2,axis=0,weights=w));sd=np.where(sd<1e-8,1.,sd);z=np.column_stack([np.ones(len(x)),(x-mean)/sd]);aa=np.mean(d*d,1);bb=np.mean(e*d,1);cc=np.mean(e*e,1)
 def loss(beta):
  q=expit(z@beta);value=float(np.sum(w*(cc+2*q*bb+q*q*aa))+.01*np.sum(beta[1:]**2)+1e-6*beta[0]**2);grad=z.T@(w*2*(bb+q*aa)*q*(1-q));grad[1:]+=.02*beta[1:];grad[0]+=2e-6*beta[0];return value,grad
 result=minimize(loss,np.array([-2.,0,0,0,0]),jac=True,method='L-BFGS-B',bounds=[(-12,12)]*5,options={'maxiter':300,'ftol':1e-12,'gtol':1e-7});assert np.isfinite(result.fun)
 return {'family':'F3','coef':result.x.tolist(),'feature_mean':mean.tolist(),'feature_scale':sd.tolist(),'optimization_success':bool(result.success),'optimization_message':str(result.message),'iterations':int(result.nit),'objective':float(result.fun),'parameters':5}
def candidates(records):
 return [{'family':'F0'}]+[{'family':'F1','lambda':v/10} for v in range(11)]+[{'family':'F2','weights':list(t)} for t in itertools.combinations_with_replacement([0.,.1,.25,.5,1.],3)]+[fit_gate(records)]
def check_lock():
 p=A/'C8A_SPEC_LOCK.json';r=json.loads(p.read_text());assert sha(p)==p.with_suffix('.sha256').read_text().split()[0]
 for name,h in r['source_sha256'].items():assert sha(ROOT/name)==h,name
 return r

def fast_score(y,p,scale):
 return {'median_pearson':core.finite_stat(core.column_corr(y,p),'median'),'common_mse':float(np.mean(((y-p)*scale)**2))}

def set_forbidden(patients):
 global FORBIDDEN
 FORBIDDEN=set(patients)
