"""Small explicit inference and fixture routes; never writes scientific archives."""
import json
from pathlib import Path
import numpy as np
from .assets import new_output,sha

def smoke(root,output):
    import torch
    from .core import mlp_model,Scaler,gene_metrics,summarize
    from .graph import build_graph,aggregate
    from .external import predict,moments
    torch.set_num_threads(1);torch.manual_seed(7)
    out=new_output(root,output)
    rng=np.random.default_rng(7);xy=np.array([(x*100,y*100) for y in range(4) for x in range(5)],float)
    x=rng.normal(size=(20,512)).astype('float32');y=rng.normal(size=(20,50))
    a=build_graph(xy,{'kind':'knn','k':6,'max_edge_um':250})
    u=aggregate(a,x,'SAGE');head=mlp_model(1024,50,{'hidden':67,'dropout':.1})
    optimizer=torch.optim.AdamW(head.parameters(),lr=.001)
    for _ in range(2):
        optimizer.zero_grad();loss=(head(torch.from_numpy(u))-torch.from_numpy(y).float()).square().mean();loss.backward();optimizer.step()
    head.eval();torch.save(head.state_dict(),out/'synthetic_graph_head.pt')
    reload=mlp_model(1024,50,{'hidden':67,'dropout':.1});reload.load_state_dict(torch.load(out/'synthetic_graph_head.pt',map_location='cpu',weights_only=True));reload.eval()
    with torch.inference_mode():
        graph=head(torch.from_numpy(u)).numpy();assert np.array_equal(graph,reload(torch.from_numpy(u)).numpy())
        heads=[mlp_model(512,50,{'hidden':256,'dropout':.1}).eval() for _ in range(3)]
        image=np.mean([h(torch.from_numpy(x)).numpy() for h in heads],axis=0)
    observed=np.array([0,4,15,19]);hidden=np.setdiff1d(np.arange(20),observed)
    scaler={'source':moments([y[:5],y[5:10]]),'destination':moments([y[10:15],y[15:]])}
    predictions,_=predict(xy,image,y,observed,hidden,'V1',scaler)
    rows=[]
    for method,p in [('graph_fixture',graph[hidden]),('local_image_fixture',predictions['M0']),('molecular_fixture',predictions['M1']),('fusion_fixture',predictions['M3'])]:rows.append({'method':method,**summarize(gene_metrics(y[hidden],p))})
    import pandas as pd
    pd.DataFrame(rows).to_csv(out/'fixture_metrics.csv',index=False)
    result={'fixture_only':True,'original_model_replay':False,'training_steps':2,'locations':20,'genes':50,'output':str(out)}
    (out/'fixture.json').write_text(json.dumps(result,indent=2)+'\n');return result

def infer(root,config_path,dry_run=False):
    """One aligned section; bounded to configured maximum locations before prediction."""
    cfg=json.loads(config_path.read_text());base=config_path.resolve().parent
    def path(s):return (base/s).resolve()
    if not cfg.get('max_locations') or not 1<=cfg['max_locations']<=10000:raise ValueError('max_locations must be 1..10000; no cohort replay in this entrypoint')
    from .assets import inventory
    inv={r['original_path']:r for r in inventory(root)}
    expected=[f'C4/C4A/models/C4A_primary_resnet18_reference_primary_S2-50_D1_h256d0.1_s{i}.pt' for i in range(3)]
    if cfg['models']!=expected:raise ValueError('Primary model identities/order must be seeds 0, 1, 2 of the selected C4 ResNet18 head')
    models=[]
    for old in cfg['models']:
        if old not in inv:raise ValueError('Model must have an original manifest identity: '+old)
        record=inv[old];p=root/record['path']
        if not p.is_file():raise FileNotFoundError(f'Original state required: {p}; obtain the original from the author, SHA256 {record["sha256"]}')
        if sha(p)!=record['sha256']:raise ValueError('Original model checksum mismatch: '+str(p))
        models.append(p)
    if len(models)!=3:raise ValueError('Exactly three primary original states required')
    required=['input','feature_scaler','target_scaler','calibration']
    for key in required:
        if not path(cfg[key]).is_file():raise FileNotFoundError(f'{key}: {path(cfg[key])}; install the aligned original asset described in docs/DATA_AND_MODELS.md')
    with np.load(path(cfg['input']),allow_pickle=False) as f:d={k:f[k] for k in f.files}
    needed=['xy_um','embeddings','targets','barcodes','genes','observed','scored']
    if any(k not in d for k in needed):raise ValueError('Input NPZ requires '+', '.join(needed))
    n=len(d['xy_um']);g=len(d['genes'])
    if n>cfg['max_locations'] or d['embeddings'].shape!=(n,512) or d['targets'].shape!=(n,g):raise ValueError('Bounded input shape mismatch')
    if g!=50 or len(np.unique(d['barcodes']))!=n:raise ValueError('Expected 50 genes and unique barcodes')
    for key in ['observed','scored']:
        if d[key].dtype!=bool or d[key].shape!=(n,):raise ValueError(key+' must be a boolean N-vector')
    if not d['observed'].any() or not d['scored'].any() or (d['observed']&d['scored']).any():raise ValueError('Observed/scored indices must be nonempty and disjoint')
    obj=json.loads(path(cfg['calibration']).read_text())
    if not np.array_equal(d['genes'],obj['genes']):raise ValueError('Gene order differs from calibration')
    for key in ['feature_scaler','target_scaler','calibration']:
        if sha(path(cfg[key]))!=cfg[key+'_sha256']:raise ValueError(key+' checksum mismatch')
    if cfg['regime'] not in ['V1','V2']:raise ValueError('Unknown regime')
    if dry_run:return {'validated':True,'locations':n,'prediction_started':False}
    import torch,pandas as pd
    from .core import mlp_model,Scaler,gene_metrics,summarize
    from .vision import ensemble
    from .external import predict
    def scaler(p):
        with np.load(p) as f:s=Scaler();s.mean=f['mean'];s.scale=f['scale']
        return s
    heads=[]
    for p in models:
        head=mlp_model(512,50,{'hidden':256,'dropout':.1});head.load_state_dict(torch.load(p,map_location='cpu',weights_only=True));heads.append(head.eval())
    image=ensemble(heads,d['embeddings'],scaler(path(cfg['feature_scaler'])),scaler(path(cfg['target_scaler'])))
    hidden=np.flatnonzero(d['scored']);preds,_=predict(d['xy_um'],image,d['targets'],np.flatnonzero(d['observed']),hidden,cfg['regime'],obj)
    out=new_output(root,path(cfg['output']))
    np.savez_compressed(out/'predictions.npz',**{m:preds[m] for m in ['M0','M1','M3']},barcodes=d['barcodes'][hidden],genes=d['genes'])
    pd.DataFrame([{'method':m,**summarize(gene_metrics(d['targets'][hidden],preds[m]))} for m in ['M0','M1','M3']]).to_csv(out/'metrics.csv',index=False)
    return {'output':str(out),'locations':n,'original_states':True}
