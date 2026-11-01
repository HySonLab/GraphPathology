"""Score the validation-selected rate, preserving every candidate and old result."""
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
import torch

from .assets import sha
from .backbone_data import STUDY,config,record,source_path,write
from .backbone_train import load,design
from .backbone_evaluate import support_manifest,array_sha
from .core import mlp_model,predict
from .extension_evaluate import score_set


def score(root,patient,item):
    base=root/f'experiments/journal_revision/{STUDY}';out=base/'primary_evaluation'/patient/item['run_id'];done=out/'COMPLETE.json';started=time.monotonic()
    if done.exists():
        r=json.loads(done.read_text())
        for name,digest in r['output_hashes'].items():assert sha(out/name)==digest
        return r
    cfg=config(root);original,ms=support_manifest(root,patient);data=load(root,item['backbone'],patient)
    for name in ['barcodes','genes','xy_um','targets']:np.testing.assert_array_equal(original[name],data[name])
    scaler_path=root/item['feature_scaler'];assert sha(scaler_path)==item['feature_scaler_sha256']
    with np.load(scaler_path) as f:x=((data['embeddings']-f['mean'])/f['scale']).astype(np.float32)
    xs,graphs=design({patient:data},{patient:x},item['model'],item['control'],item['seed']);features=xs[patient]
    checkpoint=root/item['checkpoint'];assert sha(checkpoint)==item['checkpoint_sha256'];state=torch.load(checkpoint,map_location='cpu',weights_only=True)
    net=mlp_model(features.shape[1],50,state['config']);net.load_state_dict(state['state_dict']);net.eval();raw=predict(net,features,'B2').astype(float)
    cache=root/f'data/journal_revision/{STUDY}/selected_predictions'/patient;cache.mkdir(parents=True,exist_ok=True)
    np.save(cache/(item['run_id']+'_raw.npy'),raw);image=raw.copy();calibration=None
    if patient!='TENX99':
        cp=root/'checkpoints/legacy/C9/C9A/T1_SCALER.json';cal=json.loads(cp.read_text());np.testing.assert_array_equal(data['genes'],cal['genes'])
        image=np.array(cal['destination']['center'])+np.array(cal['destination']['scale'])*(raw-np.array(cal['source']['center']))/np.array(cal['source']['scale']);calibration=sha(cp)
    np.save(cache/(item['run_id']+'_calibrated.npy'),image)
    meta={'id':item['run_id'],'origin':'matched_'+item['backbone'],'variant':item['model']+'_'+item['control'],'seed':item['seed']}
    panels=[];genes=[];bins=[];checks=0;metric_cache={};fixed_h={};independent=0
    for mask in ms:
        a,b,c,k=score_set(data,image,mask,meta,patient,metric_cache);panels.extend(a);genes.extend(b);bins.extend(c);checks+=k
        if mask['design']=='B2_identical_H':
            key=array_sha(mask['hidden']);v=image[mask['hidden']]
            if key in fixed_h:np.testing.assert_array_equal(v,fixed_h[key]);independent+=1
            else:fixed_h[key]=v.copy()
    out.mkdir(parents=True,exist_ok=True)
    for name,rows in [('PANEL.csv.gz',panels),('GENES.csv.gz',genes),('BINS.csv.gz',bins)]:pd.DataFrame(rows).to_csv(out/name,index=False)
    feature_path=(root/f'data/journal_revision/{STUDY}/phikon'/(patient+'.npz')) if item['backbone']=='phikon' else source_path(root,patient)
    r=record(root,status='COMPLETE',patient=patient,backbone=item['backbone'],model=item['model'],control=item['control'],training_seed=item['seed'],run_id=item['run_id'],learning_rate=item['learning_rate'],selected_fit_manifest_sha256=sha(base/'SELECTED_FIT_MANIFEST.json'),checkpoint_sha256=item['checkpoint_sha256'],feature_scaler_sha256=sha(scaler_path),input_cache_sha256=sha(feature_path),raw_prediction_sha256=sha(cache/(item['run_id']+'_raw.npy')),calibrated_prediction_sha256=sha(cache/(item['run_id']+'_calibrated.npy')),calibration_sha256=calibration,calibration_applications=0 if calibration is None else 1,exact_label_corruption_checks=checks,image_anchor_independence_checks=independent,mask_sets=len(ms),graphs=graphs,prediction_and_scoring_seconds=time.monotonic()-started,panel=str((out/'PANEL.csv.gz').relative_to(root)),panel_sha256=sha(out/'PANEL.csv.gz'),genes=str((out/'GENES.csv.gz').relative_to(root)),genes_sha256=sha(out/'GENES.csv.gz'),output_hashes={name:sha(out/name) for name in ['PANEL.csv.gz','GENES.csv.gz','BINS.csv.gz']})
    write(done,r);print('SELECTED MODEL EVALUATED',patient,item['run_id'],flush=True);return r


def evaluate(root,available_only=False):
    cfg=config(root);torch.set_num_threads(cfg['head_threads']);torch.use_deterministic_algorithms(True)
    base=root/f'experiments/journal_revision/{STUDY}';items=json.loads((base/'SELECTED_FIT_MANIFEST.json').read_text());assert len(items)==36
    prior={(r['patient'],r['run_id']):r for r in json.loads((base/'REUSED_EVALUATION_MANIFEST.json').read_text())}
    records=[]
    for patient in [cfg['development']]+cfg['external']:
        for item in items:
            if item['backbone']=='phikon' and available_only and not (root/f'data/journal_revision/{STUDY}/phikon'/(patient+'.json')).exists():continue
            if item['artifact_origin']=='historical_resnet_lr1e3':
                r=prior[(patient,item['run_id'])].copy()
                for name in ['panel','genes']:assert sha(root/r[name])==r[name+'_sha256']
                r['learning_rate']=.001;r['selected_fit_manifest_sha256']=sha(base/'SELECTED_FIT_MANIFEST.json')
            else:r=score(root,patient,item)
            records.append(r)
            write(base/'PRIMARY_EVALUATION_MANIFEST.json',records)
        if sum(r.get('prediction_and_scoring_seconds',0) for r in records)>cfg['evaluation_max_wall_hours']*3600:raise RuntimeError('Prespecified cumulative selected-model evaluation ceiling exceeded')
    return records

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--available-only',action='store_true');a=p.parse_args();print(json.dumps(evaluate(a.root.resolve(),a.available_only),indent=2))
