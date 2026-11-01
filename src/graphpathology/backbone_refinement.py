"""Versioned two-rate validation refinement; never reads evaluation outcomes."""
import json
import math
import time
from pathlib import Path

import numpy as np
import torch

from .assets import sha
from .backbone_data import STUDY,config,record,write
from .backbone_train import load,design
from .core import Scaler,train_mlp,mlp_model,predict,PeakRSS


def initial_trials(root):
    base=root/f'experiments/journal_revision/{STUDY}';cfg=config(root)
    old={r['run_id']:r for r in json.loads((root/'experiments/journal_revision/results/hest_learned/FIT_MANIFEST.json').read_text())}
    rows=[]
    for r in json.loads((base/'REUSED_FIT_MANIFEST.json').read_text()):
        original=old[r['run_id']]
        rows.append({**r,'parameters':r['parameters'],'config':original['config'],'learning_rate':.001,'validation_mse':original['validation_mse'],'feature_scaler':'checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz','artifact_origin':'historical_resnet_lr1e3','source_fit_record':f"experiments/journal_revision/results/hest_learned/{r['run_id']}.json",'optimizer_updates':14*r['epochs_ran']})
    for r in json.loads((base/'FIT_MANIFEST.json').read_text()):
        rows.append({**r,'learning_rate':.001,'validation_mse':r['history'][r['selected_epoch']-1]['validation_mse'],'feature_scaler':f'checkpoints/journal_revision/{STUDY}/phikon_feature_scaler.npz','artifact_origin':'initial_phikon_lr1e3','source_fit_record':str((base/'fits'/(r['run_id']+'.json')).relative_to(root))})
    assert len(rows)==36
    return rows


def train(root):
    cfg=config(root);assert cfg['training_refinement']['version']==2
    torch.set_num_threads(cfg['head_threads']);torch.use_deterministic_algorithms(True)
    out=root/f'experiments/journal_revision/{STUDY}/trials/lr1e4';out.mkdir(parents=True,exist_ok=True)
    cp=root/f'checkpoints/journal_revision/{STUDY}/trials/lr1e4';cp.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();records=[]
    for backbone in cfg['backbones']:
        data={s:load(root,backbone,s) for s in cfg['source_fit']+[cfg['validation']]}
        scaler_path=root/('checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz' if backbone=='resnet18' else f'checkpoints/journal_revision/{STUDY}/phikon_feature_scaler.npz')
        calculated=Scaler().fit(np.concatenate([data[s]['embeddings'] for s in cfg['source_fit']]))
        with np.load(scaler_path) as f:np.testing.assert_array_equal(f['mean'],calculated.mean);np.testing.assert_array_equal(f['scale'],calculated.scale)
        scaled={s:calculated.transform(d['embeddings']) for s,d in data.items()};n=sum(len(data[s]['targets']) for s in cfg['source_fit'])
        for seed in cfg['seeds']:
            for c in cfg['configurations']:
                model,control=c['model'],c['control'];ident=f'{backbone}_{model.lower()}_{control}_lr1e4_seed{seed}';path=out/(ident+'.json');checkpoint=cp/(ident+'.pt')
                if path.exists():
                    r=json.loads(path.read_text());assert sha(root/r['checkpoint'])==r['checkpoint_sha256'];records.append(r);continue
                x,graphs=design(data,scaled,model,control,seed);xt=np.concatenate([x[s] for s in cfg['source_fit']]);yt=np.concatenate([data[s]['targets'] for s in cfg['source_fit']])
                head={'hidden':cfg['hidden'][backbone][model],'dropout':cfg['dropout']};training=cfg['training']|{'learning_rate':.0001};begin=time.monotonic()
                with PeakRSS() as memory:net,fit=train_mlp(xt,yt,x[cfg['validation']],data[cfg['validation']]['targets'],head,training,seed)
                torch.save({'state_dict':net.state_dict(),'config':head,'model':model,'control':control,'seed':seed,'training':training,'backbone':backbone,'feature_scaler_sha256':sha(scaler_path),'genes':data[cfg['validation']]['genes'].tolist()},checkpoint)
                restored=mlp_model(xt.shape[1],50,head);restored.load_state_dict(torch.load(checkpoint,map_location='cpu',weights_only=True)['state_dict']);restored.eval();np.testing.assert_array_equal(predict(net,xt[:32],'B2'),predict(restored,xt[:32],'B2'))
                r=record(root,status='COMPLETE',backbone=backbone,run_id=ident,model=model,control=control,seed=seed,learning_rate=.0001,config=head,training=training,parameters=fit['parameters'],selected_epoch=fit['best_epoch'],epochs_ran=fit['epochs_run'],validation_mse=fit['history'][fit['best_epoch']-1]['validation_mse'],optimizer_updates=math.ceil(n/512)*fit['epochs_run'],selected_optimizer_updates=math.ceil(n/512)*fit['best_epoch'],checkpoint=str(checkpoint.relative_to(root)),checkpoint_sha256=sha(checkpoint),feature_scaler=str(scaler_path.relative_to(root)),feature_scaler_sha256=sha(scaler_path),source_fit=cfg['source_fit'],validation=cfg['validation'],history=fit['history'],graphs=graphs,fit_seconds=time.monotonic()-begin,peak_rss_mb=memory.peak/2**20,serialization_exact=True,artifact_origin='refinement_lr1e4',protocol_refinement_sha256=sha(root/f'experiments/journal_revision/{STUDY}/PROTOCOL_V2_TRAINING_REFINEMENT.md'),source_fit_record=str(path.relative_to(root)),input_cache_hashes={s:sha(root/(f'data/journal_revision/{STUDY}/phikon' if backbone=='phikon' else 'data/journal_revision/hest_training/processed/journal_v2_hest')/(s+'.npz')) for s in data})
                write(path,r);records.append(r);print('REFINEMENT FIT COMPLETE',ident,fit['best_epoch'],flush=True)
                if (time.monotonic()-started)/3600>cfg['head_max_wall_hours']:raise RuntimeError('Common head-fitting resource ceiling exceeded')
    write(root/f'experiments/journal_revision/{STUDY}/REFINEMENT_FIT_MANIFEST.json',records);return records


def validation_selection(rows):
    """One rate per backbone/configuration, retaining ALL seeds; no evaluation data."""
    groups={}
    for r in rows:groups.setdefault((r['backbone'],r['model'],r['control']),{}).setdefault(r['learning_rate'],[]).append(r)
    selected=[];audit=[]
    for key,rates in sorted(groups.items()):
        assert set(rates)=={.001,.0001}
        for values in rates.values():assert {r['seed'] for r in values}=={0,1,2} and len(values)==3
        means={lr:float(np.mean([r['validation_mse'] for r in values])) for lr,values in rates.items()}
        chosen=min(means,key=lambda lr:(means[lr],lr))
        selected.extend([{**r,'selected_learning_rate':chosen,'selection_mean_validation_mse':means[chosen]} for r in rates[chosen]])
        audit.append({'backbone':key[0],'model':key[1],'control':key[2],'selected_learning_rate':chosen,'mean_validation_mse_lr1e3':means[.001],'mean_validation_mse_lr1e4':means[.0001],'selection_data':'TENX95 only; mean across all three seed-specific selected losses'})
    return selected,audit


def select(root):
    base=root/f'experiments/journal_revision/{STUDY}';rows=initial_trials(root)+json.loads((base/'REFINEMENT_FIT_MANIFEST.json').read_text());assert len(rows)==72
    selected,audit=validation_selection(rows);assert len(selected)==36
    write(base/'ALL_CANDIDATE_FIT_MANIFEST.json',rows);write(base/'SELECTED_FIT_MANIFEST.json',selected);write(base/'VALIDATION_SELECTION.json',record(root,status='COMPLETE',selection_uses_evaluation_expression=False,records=audit))
    return audit


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('action',choices=['train','select']);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);a=p.parse_args();print(json.dumps((train if a.action=='train' else select)(a.root.resolve()),indent=2))
