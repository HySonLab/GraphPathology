"""Exact-support reuse and new matched-head scoring, without external selection."""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from .assets import sha
from .backbone_data import STUDY,config,record,source,write,source_path
from .backbone_train import load,design
from .core import mlp_model,predict
from .extension_evaluate import masks,score_set,cohort


def array_sha(value):
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def support_manifest(root,patient):
    data,meta=source(root,patient);rows=masks(root,data,patient)
    records=[]
    for row in rows:
        assert not np.intersect1d(row['observed'],row['hidden']).size
        records.append({k:v for k,v in row.items() if k not in ['observed','hidden']}|{'n_observed':len(row['observed']),'n_scored':len(row['hidden']),'observed_index_sha256':array_sha(row['observed']),'hidden_index_sha256':array_sha(row['hidden'])})
    report=record(root,status='COMPLETE',patient=patient,source_cache_sha256=meta['cache_sha256'],barcodes_sha256=array_sha(data['barcodes']),genes_sha256=array_sha(data['genes']),coordinates_sha256=array_sha(data['xy_um']),targets_sha256=array_sha(data['targets']),masks=records)
    write(root/f'experiments/journal_revision/{STUDY}/supports'/(patient+'.json'),report)
    return data,rows


def reuse(root):
    cfg=config(root);items=json.loads((root/f'experiments/journal_revision/{STUDY}/REUSED_FIT_MANIFEST.json').read_text());records=[]
    for patient in [cfg['development']]+cfg['external']:
        data,ms=support_manifest(root,patient);lookup={(m['design'],m['condition'],m['budget'],m['mask_seed'],m['geometry']):m for m in ms}
        for item in items:
            folder=root/'experiments/journal_revision/transformer_egnn/evaluation'/patient/('frozen_'+item['run_id'])
            old=json.loads((folder/'COMPLETE.json').read_text());assert old['model']['checkpoint_sha256']==item['checkpoint_sha256']
            for name,digest in old['output_hashes'].items():assert sha(folder/name)==digest
            panel=pd.read_csv(folder/'PANEL.csv.gz')
            assert len(panel)==len(ms)*3
            for row in panel.itertuples():
                m=lookup[(row.design,row.condition,row.budget,row.mask_seed,row.geometry)]
                assert row.mask_sha256==m['mask_sha256'] and row.n_observed==len(m['observed']) and row.n_scored==len(m['hidden'])
            # Recompute control effectiveness using only existing coordinates and the exact prior shuffle namespace.
            graph_check={}
            if item['control']=='shuffled':
                _,g=design({patient:data},{patient:data['embeddings']},item['model'],item['control'],item['seed']);graph_check=g[patient]
            prediction=root/'data/journal_revision/transformer_egnn/predictions'/patient/('frozen_'+item['run_id']+'.npy')
            if prediction.exists():assert sha(prediction)==old['prediction_sha256']
            r=record(root,status='REUSED_VERIFIED',patient=patient,backbone='resnet18',model=item['model'],control=item['control'],training_seed=item['seed'],run_id=item['run_id'],source_record=str((folder/'COMPLETE.json').relative_to(root)),source_record_sha256=sha(folder/'COMPLETE.json'),panel=str((folder/'PANEL.csv.gz').relative_to(root)),panel_sha256=sha(folder/'PANEL.csv.gz'),genes=str((folder/'GENES.csv.gz').relative_to(root)),genes_sha256=sha(folder/'GENES.csv.gz'),prediction_available=prediction.exists(),prediction_sha256=old['prediction_sha256'],prior_exact_label_corruption_checks=old['exact_label_corruption_checks'],mask_sets=len(ms),graph_check=graph_check)
            records.append(r)
        print('REUSED MATCHED RESNET',patient,18,flush=True)
    write(root/f'experiments/journal_revision/{STUDY}/REUSED_EVALUATION_MANIFEST.json',records)
    return records


def evaluate(root,available_only=False):
    cfg=config(root);torch.set_num_threads(cfg['head_threads']);torch.use_deterministic_algorithms(True)
    items=json.loads((root/f'experiments/journal_revision/{STUDY}/FIT_MANIFEST.json').read_text());assert len(items)==18
    scaler_path=root/f'checkpoints/journal_revision/{STUDY}/phikon_feature_scaler.npz'
    with np.load(scaler_path) as f:mean=f['mean'];scale=f['scale']
    start=time.monotonic();records=[]
    for patient in [cfg['development']]+cfg['external']:
        feature=root/f'data/journal_revision/{STUDY}/phikon'/(patient+'.npz')
        if available_only and not feature.with_suffix('.json').exists():continue
        base,ms=support_manifest(root,patient);data=load(root,'phikon',patient)
        for name in ['barcodes','genes','xy_um','targets']:
            np.testing.assert_array_equal(base[name],data[name])
        scaled=((data['embeddings']-mean)/scale).astype(np.float32)
        folder=root/f'experiments/journal_revision/{STUDY}/evaluation'/patient
        cache=root/f'data/journal_revision/{STUDY}/predictions'/patient;cache.mkdir(parents=True,exist_ok=True)
        for item in items:
            ident=item['run_id'];out=folder/ident;done=out/'COMPLETE.json';bundle_start=time.monotonic()
            if done.exists():
                r=json.loads(done.read_text())
                for name,digest in r['output_hashes'].items():assert sha(out/name)==digest
                records.append(r);continue
            state=torch.load(root/item['checkpoint'],map_location='cpu',weights_only=True)
            assert sha(root/item['checkpoint'])==item['checkpoint_sha256'] and sha(scaler_path)==item['feature_scaler_sha256']
            xs,graphs=design({patient:data},{patient:scaled},item['model'],item['control'],item['seed']);x=xs[patient]
            net=mlp_model(x.shape[1],50,state['config']);net.load_state_dict(state['state_dict']);net.eval()
            raw=predict(net,x,'B2').astype(float);np.save(cache/(ident+'_raw.npy'),raw)
            image=raw.copy();calibration=None
            if patient!=cfg['development']:
                cal_path=root/'checkpoints/legacy/C9/C9A/T1_SCALER.json';cal=json.loads(cal_path.read_text());np.testing.assert_array_equal(data['genes'],cal['genes'])
                image=np.array(cal['destination']['center'])+np.array(cal['destination']['scale'])*(raw-np.array(cal['source']['center']))/np.array(cal['source']['scale']);calibration=sha(cal_path)
            np.save(cache/(ident+'_calibrated.npy'),image)
            meta={'id':ident,'origin':'matched_phikon','variant':item['model']+'_'+item['control'],'seed':item['seed']}
            panels=[];genes=[];bins=[];checks=0;metric_cache={};fixed_h_images={};anchor_independence_checks=0
            for mask in ms:
                a,b,c,k=score_set(data,image,mask,meta,patient,metric_cache);panels.extend(a);genes.extend(b);bins.extend(c);checks+=k
                if mask['design']=='B2_identical_H':
                    key=array_sha(mask['hidden']);value=image[mask['hidden']]
                    if key in fixed_h_images:np.testing.assert_array_equal(value,fixed_h_images[key]);anchor_independence_checks+=1
                    else:fixed_h_images[key]=value.copy()
            out.mkdir(parents=True,exist_ok=True)
            for name,rows in [('PANEL.csv.gz',panels),('GENES.csv.gz',genes),('BINS.csv.gz',bins)]:pd.DataFrame(rows).to_csv(out/name,index=False)
            r=record(root,status='COMPLETE',patient=patient,backbone='phikon',model=item['model'],control=item['control'],training_seed=item['seed'],run_id=ident,checkpoint_sha256=item['checkpoint_sha256'],feature_scaler_sha256=sha(scaler_path),input_cache_sha256=sha(root/f'data/journal_revision/{STUDY}/phikon'/(patient+'.npz')),prediction_and_scoring_seconds=time.monotonic()-bundle_start,raw_prediction_sha256=sha(cache/(ident+'_raw.npy')),calibrated_prediction_sha256=sha(cache/(ident+'_calibrated.npy')),calibration_sha256=calibration,calibration_applications=0 if calibration is None else 1,exact_label_corruption_checks=checks,image_anchor_independence_checks=anchor_independence_checks,graphs=graphs,mask_sets=len(ms),output_hashes={name:sha(out/name) for name in ['PANEL.csv.gz','GENES.csv.gz','BINS.csv.gz']})
            write(done,r);records.append(r);print('MATCHED EVALUATED',patient,ident,flush=True)
            if (time.monotonic()-start)/3600>cfg['evaluation_max_wall_hours']:raise RuntimeError('Prespecified evaluation resource ceiling exceeded; completed bundles preserved')
        write(root/f'experiments/journal_revision/{STUDY}/EVALUATION_MANIFEST.json',records)
    return records


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['reuse','evaluate']);parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);parser.add_argument('--available-only',action='store_true');a=parser.parse_args()
    print(json.dumps(reuse(a.root.resolve()) if a.action=='reuse' else evaluate(a.root.resolve(),a.available_only),indent=2))
