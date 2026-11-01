"""Matched cached-aggregation heads; old ResNet18 states remain immutable."""
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from .assets import sha
from .backbone_data import STUDY,config,record,source,write
from .core import Scaler,PeakRSS,train_mlp,predict,mlp_model
from .graph import build_graph,aggregate,shuffle_graph,parameter_count


def load(root, backbone, patient):
    if backbone=='resnet18':return source(root,patient)[0]
    path=root/f'data/journal_revision/{STUDY}/phikon'/(patient+'.npz')
    meta=json.loads(path.with_suffix('.json').read_text());assert sha(path)==meta['cache_sha256']
    with np.load(path) as f:data={k:f[k].copy() for k in f.files}
    return data


def design(data,scaled,model,control,seed):
    xs={};graphs={}
    for patient,d in data.items():
        x=scaled[patient]
        if model=='local':xs[patient]=x;continue
        a=build_graph(d['xy_um'],{'kind':'delaunay','max_edge_um':250});original=a
        info={'edges':a.nnz//2}
        if control=='shuffled':
            a,permutation=shuffle_graph(a,patient,seed)
            changed=1-a.multiply(original).nnz/original.nnz
            if changed<=.5:raise ValueError('Ineffective shuffled-graph control')
            info.update(changed_edge_fraction=changed,permutation_sha256=__import__('hashlib').sha256(permutation.tobytes()).hexdigest())
        info['adjacency_sha256']=__import__('hashlib').sha256(a.indices.tobytes()+a.indptr.tobytes()).hexdigest()
        graphs[patient]=info;xs[patient]=aggregate(a,x,model)
    return xs,graphs


def audit(root):
    """Audit contracts and hashes without opening new head evaluation outcomes."""
    cfg=config(root);old=json.loads((root/'experiments/journal_revision/results/hest_learned/FIT_MANIFEST.json').read_text())
    rows=[]
    for r in old:
        if (r['model'],r['control']) not in [(c['model'],c['control']) for c in cfg['configurations']]:continue
        cp=root/r['checkpoint_relative_path'];assert sha(cp)==r['checkpoint_sha256']
        rows.append({'backbone':'resnet18','run_id':r['run_id'],'status':'REUSED_VERIFIED','model':r['model'],'control':r['control'],'seed':r['seed'],'embedding_dimension':512,'hidden':r['config']['hidden'],'parameters':r['n_parameters'],'checkpoint':r['checkpoint_relative_path'],'checkpoint_sha256':r['checkpoint_sha256'],'selection_patient':r['selection_patient'],'fitting_patients':r['fit_patient_ids'],'training':r['training'],'selected_epoch':r['selected_epoch'],'epochs_ran':r['epochs_ran'],'feature_scaler_sha256':r['feature_scaler_sha256'],'fit_seconds':r['fit_seconds'],'peak_rss_mb':r['peak_rss_mb']})
    assert len(rows)==18
    out=root/f'experiments/journal_revision/{STUDY}'
    write(out/'REUSED_FIT_MANIFEST.json',rows)
    text='''# Existing backbone comparison audit

Source snapshot: ab0811134bb1d62719206ab11a43bd7c4dd5d813. Contract audit precedes new matched outcome inspection. Exact raw records and checksums are in REUSED_FIT_MANIFEST.json; source protocols and original outputs remain immutable.

## Usable matched ResNet18 configurations

Every row below uses frozen torchvision ResNet18 IMAGENET1K_V1 (official f37072fd weights,512-D), official resize256/center-crop224 normalization, identical registered112um-field RGB patches (effective ResNet field87.5% of that width, approximately98um nominal), source-only scaler fit on NCBI783/NCBI785, C0-normalized50-gene targets without another head target scaler, AdamW0.001/decay0.0001, pooled512-row minibatches/max80/patience10/min-delta0.0001. TENX95 alone selects minimum MSE. All external heads receive the frozen C9 T1 map once. The exact source/development/external masks0/1/2 and all three routes are available in extension evaluation bundles. Those are distinct from original ten-mask headline analyses. Checkpoint, prediction, support and export hashes are checked before reuse; no arrays will be invented.

| Configuration | Seed | Parameters | Selected/run epochs | Checkpoint | Reuse |
|---|---:|---:|---:|---|---|
'''
    for r in rows:text+=f"| {r['model']} / {r['control']} | {r['seed']} | {r['parameters']} | {r['selected_epoch']}/{r['epochs_ran']} | `{r['checkpoint']}` | REUSED_VERIFIED |\n"
    text+='''
## Other existing configurations and limits

- Original three selected ResNet18 width256 image heads: historical primary ensemble,144178 parameters each; binaries remain missing. Exported original results remain preserved, but these are not substitutes for the matched width128 minibatch comparison.
- Later full-section ResNet18 local/mean/GCN/SAGE fits: all48 fits and advanced-model controls remain retained. One update per training section per epoch differs from pooled minibatches; these are not optimization-matched to the18 reused states.
- Historical C4 Phikon local MLP/ridge and C6 spatial sensitivity: final-LayerNorm CLS768-D, pinned owkin/phikon057cc0295895c2df3dd7681a89680da6015cbefe, official224px ViTImageProcessor/ImageNet normalization. The local selected width256 and an additional head-target scaler used a validation-Pearson selection grid; spatial widths128/67 were inherited. Development/conditional holdout exports, preprocessing and extraction registries exist, but exact embedding arrays and original primary model binaries are absent. These are valid historical development evidence, not the new matched contract or external predictions. Their registered FOV was approximately112um; no native-resolution recovery is claimed.
- New Phikon matched configurations: local,MEAN,GCN,SAGE,GCN-shuffled,SAGE-shuffled × seeds0/1/2. Source-only scaler is separately fit on the new768-D embeddings, direct common C0 targets, same pooled optimization/selection opportunities as reused ResNet. Local/MEAN/GCN128 (104882 parameters), SAGE66 (104792). Inputs and weights are hash-verified; these are NEW_FITS_REQUIRED until completed, never recovered historical heads.

## Common comparison contract

Fit NCBI783/NCBI785, select TENX95, report TENX99 separately, Wu1160920F/CID4290/CID44971 and the locked22 Bassiouni patients separately. No external expression determines tuning or calibration. Predictions/routes: image-only; fixed fusion0.2g+0.8IDW separately; observed-residual g+IDW(y_O-g_O). Use identical 10/20/40% historical neural masks and B2 supports across all comparators; primary20% distributed. The molecular IDW baseline always remains in anchor-assisted comparisons. All per-run/per-gene/prediction artifacts are identified in the study manifests. Encoder changes include pretraining, dimension and transforms; only within-backbone contrasts isolate added context. Fixed aggregation+MLP heads are not end-to-end encoder training or new message-passing architectures.
'''
    (out/'EXISTING_BACKBONE_COMPARISON_AUDIT.md').write_text(text)
    return rows


def train(root):
    cfg=config(root);torch.set_num_threads(cfg['head_threads']);torch.use_deterministic_algorithms(True)
    data={s:load(root,'phikon',s) for s in cfg['source_fit']+[cfg['validation']]}
    for d in data.values():np.testing.assert_array_equal(d['genes'],data[cfg['validation']]['genes'])
    scaler=Scaler().fit(np.concatenate([data[s]['embeddings'] for s in cfg['source_fit']]))
    cp=root/f'checkpoints/journal_revision/{STUDY}';cp.mkdir(parents=True,exist_ok=True)
    scaler_path=cp/'phikon_feature_scaler.npz'
    if scaler_path.exists():
        with np.load(scaler_path) as f:
            np.testing.assert_array_equal(f['mean'],scaler.mean);np.testing.assert_array_equal(f['scale'],scaler.scale)
    else:scaler.save(scaler_path)
    scaled={s:scaler.transform(d['embeddings']) for s,d in data.items()}
    out=root/f'experiments/journal_revision/{STUDY}/fits';out.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();records=[];n=sum(len(data[s]['targets']) for s in cfg['source_fit'])
    for seed in cfg['seeds']:
        for item in cfg['configurations']:
            model,control=item['model'],item['control'];run_id=f'phikon_{model.lower()}_{control}_seed{seed}'
            path=out/(run_id+'.json');state_path=cp/(run_id+'.pt')
            if path.exists():
                r=json.loads(path.read_text());assert sha(root/r['checkpoint'])==r['checkpoint_sha256'];records.append(r);continue
            x,graphs=design(data,scaled,model,control,seed)
            xt=np.concatenate([x[s] for s in cfg['source_fit']]);yt=np.concatenate([data[s]['targets'] for s in cfg['source_fit']])
            head={'hidden':cfg['hidden']['phikon'][model],'dropout':cfg['dropout']}
            # Source-only implementation check, before validation-based fitting.
            torch.manual_seed(seed);smoke=mlp_model(xt.shape[1],50,head)
            smoke(torch.tensor(xt[:16])).square().mean().backward()
            begin=time.monotonic()
            with PeakRSS() as memory:net,fit=train_mlp(xt,yt,x[cfg['validation']],data[cfg['validation']]['targets'],head,cfg['training'],seed)
            torch.save({'state_dict':net.state_dict(),'config':head,'model':model,'control':control,'seed':seed,'training':cfg['training'],'feature_scaler_sha256':sha(scaler_path),'genes':data[cfg['validation']]['genes'].tolist()},state_path)
            restored=mlp_model(xt.shape[1],50,head);restored.load_state_dict(torch.load(state_path,map_location='cpu',weights_only=True)['state_dict']);restored.eval()
            np.testing.assert_array_equal(predict(net,xt[:32],'B2'),predict(restored,xt[:32],'B2'))
            assert fit['parameters']==parameter_count(xt.shape[1],head['hidden'],50)
            r=record(root,status='COMPLETE',backbone='phikon',run_id=run_id,model=model,control=control,seed=seed,config=head,training=cfg['training'],parameters=fit['parameters'],selected_epoch=fit['best_epoch'],epochs_ran=fit['epochs_run'],optimizer_updates=math.ceil(n/cfg['training']['batch_size'])*fit['epochs_run'],selected_optimizer_updates=math.ceil(n/cfg['training']['batch_size'])*fit['best_epoch'],checkpoint=str(state_path.relative_to(root)),checkpoint_sha256=sha(state_path),feature_scaler_sha256=sha(scaler_path),input_cache_hashes={s:sha(root/f'data/journal_revision/{STUDY}/phikon'/(s+'.npz')) for s in data},source_fit=cfg['source_fit'],validation=cfg['validation'],graphs=graphs,fit_seconds=time.monotonic()-begin,peak_rss_mb=memory.peak/2**20,history=fit['history'],serialization_exact=True)
            r['protocol_json_sha256']=sha(root/f'configs/journal_revision/{STUDY}/protocol.json')
            r['implementation_hashes']={name:sha(root/'src/graphpathology'/name) for name in ['core.py','graph.py','backbone_train.py','backbone_data.py']}
            write(path,r);records.append(r);print('MATCHED FIT COMPLETE',run_id,fit['best_epoch'],flush=True)
            if (time.monotonic()-start)/3600>cfg['head_max_wall_hours'] or memory.peak/2**30>cfg['max_rss_gib']:raise RuntimeError('Declared head resource ceiling exceeded')
    write(root/f'experiments/journal_revision/{STUDY}/FIT_MANIFEST.json',records);return records


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['audit','train']);parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);a=parser.parse_args()
    print(json.dumps((audit if a.action=='audit' else train)(a.root.resolve()),indent=2))
