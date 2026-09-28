"""Execute locked C3A controls; never train or modify a C0-C2 artifact."""
from common import *
import argparse,time
from scipy import sparse
from threadpoolctl import threadpool_limits
import torch
from graph_controls import randomize,audit,seed_for,canonical_hash
from graph import aggregate
GCACHE={}

def graph(sample,geometry,control,gseed):
 key=(sample,geometry,control,gseed)
 if key in GCACHE:return GCACHE[key]
 reference=sparse.load_npz(ROOT/'C2/graphs'/f'{sample}_{geometry}.npz')
 with np.load(ROOT/'data/cache'/f'{sample}_targets.npz') as f:xy=f['xy_um']
 tag=f'{sample}_{geometry}_{control}_g{gseed}';path=OUT/'graphs'/f'{tag}.npz';mpath=path.with_suffix('.json')
 if mpath.exists():
  r=json.loads(mpath.read_text());assert sha(path)==r['file_sha256'];a=sparse.load_npz(path)
 else:
  start=time.perf_counter();derived=seed_for(sample,geometry,control,gseed)
  if control in ['A4','A5','A6']:a,info=randomize(reference,xy,control,derived)
  else:a=reference;info={'algorithm':'frozen C2 graph','accepted_swaps':0,'attempts':0}
  stats=audit(a,reference,xy,control)
  sparse.save_npz(path,a)
  r={'sample_id':sample,'geometry':geometry,'graph_control_type':control,'graph_randomization_seed':gseed,'derived_seed':derived,**stats,**info,'construction_seconds':time.perf_counter()-start,'path':str(path.relative_to(ROOT)),'file_sha256':sha(path),'reference_path':str((ROOT/'C2/graphs'/f'{sample}_{geometry}.npz').relative_to(ROOT)),'reference_sha256':sha(ROOT/'C2/graphs'/f'{sample}_{geometry}.npz'),'coordinates_sha256':hashlib.sha256(xy.tobytes()).hexdigest(),'source_commit':commit()}
  dump(mpath,r);registry('randomization' if control in ['A4','A5','A6'] else 'graph_reference',tag,r)
  print('GRAPH',tag,'changed',round(r['changed_edge_fraction'],3),flush=True)
 GCACHE[key]=(a,r);return a,r

def features(sample,data,panel,control,gseed,heldout=None):
 x=scaler(preprocessing(panel,heldout),'feature').transform(data['x'])
 if control=='A1':return np.concatenate([x,x],axis=1),{'graph_hash':'SELF_ONLY_DUPLICATED_LOCAL','graph_control_type':control,'num_nodes':len(x),'num_edges':0}
 geometry=SELECTED[panel]['SAGE']['graph']['name'];a,r=graph(sample,geometry,control,gseed)
 return aggregate(a,x,'SAGE'),r

def fit(phase,panel,control,seed,gseed,data,heldout=None):
 run_id=f'C3A_{phase}_{heldout or "primary"}_{panel}_{control}_s{seed}_g{gseed}'
 path=OUT/'logs/runs'/f'{run_id}.json'
 if path.exists():return json.loads(path.read_text())
 train=[s for s in IDS if s!=heldout] if heldout else SPLITS['splits']['train'];val=[] if heldout else SPLITS['splits']['validation']
 folder=preprocessing(panel,heldout);ys=scaler(folder,'target');selected=SELECTED[panel]['SAGE'];start=time.perf_counter();refs={}
 def make(ids):
  arrays=[]
  for s in ids:
   z,ref=features(s,data[s],panel,control,gseed,heldout);arrays.append(z);refs[s]=ref
  return np.concatenate(arrays)
 x=make(train);y=ys.transform(np.concatenate([data[s]['y'][:,ix(panel)] for s in train]));xv=make(val) if val else None;yv=ys.transform(np.concatenate([data[s]['y'][:,ix(panel)] for s in val])) if val else None
 prepare=time.perf_counter()-start
 with core.PeakRSS() as mem:net,info=core.train_mlp(x,y,xv,yv,selected['params'],base.C1CFG['mlp'],seed,selected['lopo_epochs'] if heldout else None)
 modelpath=OUT/'models'/f'{run_id}.pt';torch.save(net.state_dict(),modelpath)
 history=info.pop('history');dump(OUT/'logs/epochs'/f'{run_id}.json',history)
 if val:
  p=ys.inverse(core.predict(net,xv,'B2'));m=core.gene_metrics(np.concatenate([data[s]['y'][:,ix(panel)] for s in val]),p)
  pd.DataFrame({'gene':PANELS[panel],**m}).to_csv(OUT/'logs/validation'/f'{run_id}.csv',index=False)
 r={'run_id':run_id,'phase':phase,'panel':panel,'control':control,'training_seed':seed,'graph_randomization_seed':gseed,'heldout':heldout,'params':selected['params'],'input_dim':x.shape[1],'output_dim':y.shape[1],'train_slides':train,'validation_slides':val,'model_path':str(modelpath.relative_to(ROOT)),'model_sha256':sha(modelpath),'preprocessing':str(folder.relative_to(ROOT)),'feature_scaler_sha256':sha(folder/'feature_scaler.npz'),'target_scaler_sha256':sha(folder/'target_scaler.npz'),'graphs':refs,'prepare_seconds':prepare,'peak_rss_mb':mem.peak/1024**2,'source_commit':commit(),'spec_sha256':sha(OUT/'C3A_SPEC_LOCK.json'),**info}
 dump(path,r);registry('fit',run_id,{**r,'evaluation_slides':val or [heldout],'graph_hash':'|'.join(refs[s]['graph_hash'] for s in train+val),'graph_control_type':control})
 print('FIT',run_id,'epoch',r['best_epoch'],flush=True);return r

def model_load(r):
 assert sha(ROOT/r['model_path'])==r['model_sha256']
 net=core.mlp_model(r['input_dim'],r['output_dim'],r['params']);net.load_state_dict(torch.load(ROOT/r['model_path'],map_location='cpu',weights_only=True));net.eval();return net

def evaluate(r,sample,d):
 path=OUT/'logs/evaluation'/f'{r["run_id"]}_{sample}.json'
 if path.exists():return
 start=time.perf_counter();x,g=features(sample,d,r['panel'],r['control'],r['graph_randomization_seed'],r['heldout']);ys=scaler(ROOT/r['preprocessing'],'target');p=ys.inverse(core.predict(model_load(r),x,'B2'));y=d['y'][:,ix(r['panel'])]
 metrics=core.gene_metrics(y,p);summary=core.summarize(metrics);co,_,_=core.coexpression(y,p)
 ppath=OUT/'predictions'/f'{r["run_id"]}_{sample}.npz';np.savez_compressed(ppath,prediction=p,genes=np.array(PANELS[r['panel']]),barcodes=d['barcodes'])
 columns={k:r[k] for k in ['run_id','phase','panel','control','training_seed','graph_randomization_seed']};columns['sample_id']=sample
 frame=pd.DataFrame({'gene':PANELS[r['panel']],**metrics}).assign(**columns);gp=OUT/'logs/gene_metrics'/f'{r["run_id"]}_{sample}.csv';frame.to_csv(gp,index=False)
 out={**columns,**summary,**co,'prediction_path':str(ppath.relative_to(ROOT)),'prediction_sha256':sha(ppath),'gene_metrics_path':str(gp.relative_to(ROOT)),'graph_hash':g['graph_hash'],'graph_control_type':r['control'],'inference_seconds':time.perf_counter()-start,'spec_sha256':sha(OUT/'C3A_SPEC_LOCK.json')}
 dump(path,out);registry('evaluation',r['run_id'],{**out,'evaluation_slides':[sample],'graph':g})
 print('EVAL',r['run_id'],round(summary['median_pearson'],4),flush=True)

def reuse(data):
 # Import exact frozen metrics. Coexpression is evaluated from saved arrays only; no old model runs.
 for checkpoint,controls in [('C1',{'B2':'A0'}),('C2',{'MEAN':'A2','SAGE':'A3'})]:
  frames=pd.read_csv(ROOT/checkpoint/f'{checkpoint}_GENE_METRICS.csv')
  chosen=frames[frames.model.isin(controls)&frames.phase.isin(['primary','lopo'])]
  for (run_id,sample),df in chosen.groupby(['run_id','sample_id']):
   control=controls[df.model.iloc[0]];cid=f'C3A_reuse_{run_id}';path=OUT/'logs/evaluation'/f'{cid}_{sample}.json'
   if path.exists():continue
   r=json.loads((ROOT/checkpoint/'logs/evaluation'/f'{run_id}_{sample}.json').read_text());panel=r['panel']
   assert sha(ROOT/r['prediction_path'])==r['prediction_sha256']
   with np.load(ROOT/r['prediction_path']) as f:p=f['prediction'].copy();assert np.array_equal(f['barcodes'],data[sample]['barcodes']);assert list(f['genes'])==PANELS[panel]
   if 'coexpression_pearson' not in r:co,_,_=core.coexpression(data[sample]['y'][:,ix(panel)],p);r.update(co)
   columns={'run_id':cid,'phase':r['phase'],'panel':panel,'control':control,'training_seed':r['seed'],'graph_randomization_seed':-1,'sample_id':sample}
   original=ROOT/r['gene_metrics_path'];frame=pd.read_csv(original);frame=frame[['gene','pearson','spearman','mse','mae','prediction_variance','target_variance']].assign(**columns)
   gp=OUT/'logs/gene_metrics'/f'{cid}_{sample}.csv';frame.to_csv(gp,index=False)
   if control in ['A2','A3']:
    geom=SELECTED[panel]['MEAN' if control=='A2' else 'SAGE']['graph']['name'];_,g=graph(sample,geom,control,-1)
   else:g={'graph_hash':'NO_GRAPH','graph_control_type':control}
   out={**r,**columns,'reused_from':str((ROOT/checkpoint/'logs/evaluation'/f'{run_id}_{sample}.json').relative_to(ROOT)),'frozen_gene_metrics_path':str(original.relative_to(ROOT)),'gene_metrics_path':str(gp.relative_to(ROOT)),'graph_hash':g['graph_hash'],'graph_control_type':control,'spec_sha256':sha(OUT/'C3A_SPEC_LOCK.json')}
   dump(path,out);registry('reuse_evaluation',cid,{**out,'evaluation_slides':[sample],'graph':g})

def prepare():
 for sample in IDS:
  for geometry in sorted({SELECTED[p]['SAGE']['graph']['name'] for p in PANELS}):
   for control in ['A4','A5','A6']:
    for gseed in CFG['graph_randomization_seeds']:graph(sample,geometry,control,gseed)

def train_primary():
 data={s:load(s) for s in SPLITS['splits']['train']+SPLITS['splits']['validation']}
 for panel in PANELS:
  for control in CFG['primary_new_controls']:
   for gseed in ([-1] if control=='A1' else CFG['graph_randomization_seeds']):
    for seed in CFG['training_seeds']:fit('primary',panel,control,seed,gseed,data)

def final_evaluation():
 if not (OUT/'logs/TEST_ACCESS.json').exists():dump(OUT/'logs/TEST_ACCESS.json',{'timestamp':now(),'spec_sha256':sha(OUT/'C3A_SPEC_LOCK.json'),'prior_C2_test_results_known':True})
 data={s:load(s) for s in IDS};reuse(data)
 for path in sorted((OUT/'logs/runs').glob('*.json')):
  r=json.loads(path.read_text())
  if r['phase']=='primary':evaluate(r,SPLITS['splits']['test'][0],data[SPLITS['splits']['test'][0]])
 for heldout in IDS:
  for panel in PANELS:
   for control in CFG['conditional_lopo_new_controls']:
    for gseed in ([-1] if control=='A1' else CFG['graph_randomization_seeds']):
     for seed in CFG['training_seeds']:
      r=fit('lopo',panel,control,seed,gseed,data,heldout);evaluate(r,heldout,data[heldout])

def main():
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','train_primary','final_evaluation']);args=parser.parse_args();check_lock()
 for folder in ['logs/runs','logs/epochs','logs/validation','logs/gene_metrics','logs/evaluation']:(OUT/folder).mkdir(exist_ok=True)
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
 with threadpool_limits(4):globals()[args.stage]()
if __name__=='__main__':main()
