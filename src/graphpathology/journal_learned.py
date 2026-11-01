"""New, explicitly labelled CPU development fits; never substitutes historical checkpoints."""
import copy
import json
from pathlib import Path
import time
import h5py
import numpy as np
import pandas as pd
import torch
from torchvision.models import resnet18, ResNet18_Weights
from PIL import Image
from scipy import sparse
from .assets import sha
from .data import read_slide, validate_patient_splits
from .core import Scaler, train_mlp, gene_metrics, summarize, PeakRSS
from .graph import build_graph, aggregate, shuffle_graph, parameter_count
from .journal_data import original_genes, normalized
from .journal_methods import predict_anchors

SAMPLES = ['NCBI783', 'NCBI785', 'TENX95', 'TENX99']
SPLITS = {'train': SAMPLES[:2], 'validation': ['TENX95'], 'test': ['TENX99']}
CONFIGURATIONS = [('local', 'local')] + [(model, control) for model in ['MEAN', 'GCN', 'SAGE'] for control in ['spatial', 'self', 'shuffled']]
TRAINING = {'learning_rate': 1e-3, 'weight_decay': 1e-4, 'batch_size': 512, 'max_epochs': 80, 'patience': 10, 'min_delta': 1e-4}


def prepare(root, data_root):
    from .journal import provenance
    cache = data_root / 'processed/journal_v2_hest'; cache.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((data_root / 'DATA_ACQUISITION_MANIFEST.json').read_text())
    for r in manifest['files']:
        if r['cohort'] in ['hest', 'pretrained_encoder']:
            if r['status'] != 'AVAILABLE' or sha(data_root / r['local']) != r['sha256']: raise ValueError('Unverified HEST/encoder input')
    raw = data_root / 'raw/hest'; metadata = pd.read_csv(raw / 'HEST_v1_1_0.csv').set_index('id')
    validate_patient_splits(SPLITS, {s: str(metadata.loc[s, 'patient']) for s in SAMPLES})
    torch.set_num_threads(4); torch.use_deterministic_algorithms(True); torch.manual_seed(0)
    weight = data_root / 'models/resnet18-f37072fd.pth'; weight_sha = sha(weight)
    if not weight_sha.startswith('f37072fd'): raise ValueError('Official ResNet18 weight prefix mismatch')
    net = resnet18(weights=None); net.load_state_dict(torch.load(weight, map_location='cpu', weights_only=True)); net.fc = torch.nn.Identity(); net.eval()
    for p in net.parameters(): p.requires_grad_(False)
    transform = ResNet18_Weights.IMAGENET1K_V1.transforms(); genes = original_genes(root)
    records = []
    for sample in SAMPLES:
        dest = cache / (sample + '.npz'); meta = dest.with_suffix('.json')
        if dest.exists():
            record = json.loads(meta.read_text())
            if sha(dest) != record['cache_sha256']: raise ValueError('New feature cache modified')
            records.append(record); print('FEATURES RESUMED', sample, flush=True); continue
        start = time.monotonic(); slide = read_slide(raw, sample, float(metadata.loc[sample, 'pixel_size_um_estimated']))
        idx = [list(slide.genes).index(g) for g in genes]; raw_counts = slide.counts[:, idx].astype(float)
        y = normalized(root, raw_counts).astype(np.float32); features = []
        with h5py.File(slide.patch_path) as f:
            for begin in range(0, len(slide.barcodes), 64):
                images = f['img'][begin:begin + 64]
                batch = torch.stack([transform(Image.fromarray(im)) for im in images])
                with torch.inference_mode(): features.append(net(batch).numpy())
                if begin % 1280 == 0: print('FEATURES', sample, begin, '/', len(slide.barcodes), flush=True)
        x = np.concatenate(features)
        if x.shape != (len(y), 512) or not np.isfinite(x).all(): raise ValueError('Invalid frozen features')
        np.savez_compressed(dest, embeddings=x, targets=y, xy_um=slide.xy_um, barcodes=slide.barcodes, genes=np.array(genes))
        record = provenance(root, sample=sample, n_locations=len(y), cache_sha256=sha(dest), encoder_sha256=weight_sha,
                            max_alignment_error_px=slide.coordinate_error_px, excluded_expression_locations=slide.excluded_expression_locations,
                            elapsed_seconds=time.monotonic()-start, encoder='frozen ResNet18 IMAGENET1K_V1; newly computed cache',
                            source_expression_sha256=sha(raw / 'IDC/adata' / (sample + '.h5ad')), source_patches_sha256=sha(slide.patch_path),
                            normalization_sha256=sha(root/'checkpoints/legacy/reports/C0/C0_NORMALIZATION.json'))
        meta.write_text(json.dumps(record, indent=2)+'\n'); records.append(record); print('FEATURES COMPLETE', sample, flush=True)
    (cache/'MANIFEST.json').write_text(json.dumps(records, indent=2)+'\n'); return records


def load_prepared(data_root):
    cache = data_root/'processed/journal_v2_hest'; data={}
    for sample in SAMPLES:
        p=cache/(sample+'.npz'); meta=json.loads(p.with_suffix('.json').read_text())
        if sha(p)!=meta['cache_sha256']: raise ValueError('Prepared HEST cache checksum mismatch')
        with np.load(p) as f:data[sample]={k:f[k].copy() for k in f.files}
    return data


def representation(data, scaled, graphs, model, control, seed):
    result={}
    for sample in SAMPLES:
        x=scaled[sample]
        if model=='local': result[sample]=x;continue
        a=graphs[sample]
        if control=='self': a=sparse.eye(len(x),format='csr') if model=='SAGE' else sparse.csr_matrix(a.shape)
        elif control=='shuffled': a,_=shuffle_graph(a,sample,seed)
        result[sample]=aggregate(a,x,model)
    return result


def train(root,data_root,output,checkpoint_root):
    from .journal import provenance
    output.mkdir(parents=True,exist_ok=True);checkpoint_root.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
    data=load_prepared(data_root); feature_scaler=Scaler().fit(np.concatenate([data[s]['embeddings'] for s in SAMPLES[:2]]))
    scaler_path=checkpoint_root/'journal_v2_feature_scaler.npz'
    if not scaler_path.exists():feature_scaler.save(scaler_path)
    else:
        with np.load(scaler_path) as f:np.testing.assert_array_equal(f['mean'],feature_scaler.mean);np.testing.assert_array_equal(f['scale'],feature_scaler.scale)
    scaled={s:feature_scaler.transform(data[s]['embeddings']) for s in SAMPLES}
    graphs={s:build_graph(data[s]['xy_um'],{'kind':'delaunay','max_edge_um':250}) for s in SAMPLES}
    graph_records={s:{'n_vertices':len(data[s]['xy_um']),'undirected_edges':graphs[s].nnz//2,'adjacency_sha256':__import__('hashlib').sha256(graphs[s].indices.tobytes()+graphs[s].indptr.tobytes()).hexdigest()} for s in SAMPLES}
    histories=[];run_records=[]
    # Execute seed0 local first as specified, then every controlled configuration/seed.
    for seed in [0,1,2]:
        for model,control in CONFIGURATIONS:
            run_id=f'journal_v2_{model.lower()}_{control}_seed{seed}'; record_path=output/(run_id+'.json');checkpoint=checkpoint_root/(run_id+'.pt')
            if record_path.exists():
                record=json.loads(record_path.read_text())
                if sha(checkpoint)!=record['checkpoint_sha256']:raise ValueError('Learned checkpoint modified')
                run_records.append(record);print('FIT RESUMED',run_id,flush=True);continue
            x=representation(data,scaled,graphs,model,control,seed);xt=np.concatenate([x[s] for s in SAMPLES[:2]]);yt=np.concatenate([data[s]['targets'] for s in SAMPLES[:2]])
            config={'hidden':67 if model=='SAGE' else 128,'dropout':.1};start=time.monotonic()
            with PeakRSS() as memory:net,fit=train_mlp(xt,yt,x['TENX95'],data['TENX95']['targets'],config,TRAINING,seed)
            net.eval()
            # Save new learned state before opening evaluation outcomes.
            torch.save({'state_dict':net.state_dict(),'family':'new_journal_v2','model':model,'control':control,'seed':seed,'config':config,'training':TRAINING,'genes':data['TENX99']['genes'].tolist(),'feature_scaler_sha256':sha(scaler_path)},checkpoint)
            pd.DataFrame(fit['history']).to_csv(output/(run_id+'_history.csv'),index=False)
            record=provenance(root,run_id=run_id,status='TRAINED',model=model,control=control,seed=seed,config=config,training=TRAINING,
                              checkpoint_relative_path=str(checkpoint.relative_to(root)),checkpoint_sha256=sha(checkpoint),feature_scaler_sha256=sha(scaler_path),
                              n_parameters=sum(p.numel() for p in net.parameters()),selected_epoch=fit['best_epoch'],epochs_ran=len(fit['history']),
                              validation_mse=fit['history'][fit['best_epoch']-1]['validation_mse'],fit_patient_ids=SAMPLES[:2],selection_patient='TENX95',evaluation_patient='TENX99',
                              fit_seconds=time.monotonic()-start,peak_rss_mb=memory.peak/2**20,graph_manifest=graph_records)
            record_path.write_text(json.dumps(record,indent=2)+'\n');run_records.append(record);print('FIT COMPLETE',run_id,'epoch',fit['best_epoch'],flush=True)
    (output/'FIT_MANIFEST.json').write_text(json.dumps(run_records,indent=2)+'\n');return run_records


def hybrids(xy,y,image,observed,hidden):
    config={'kind':'idw','k':16,'p':2}
    anchor,_=predict_anchors(xy,y,observed,hidden,config)
    residual_values=np.full_like(y,np.nan,dtype=float);residual_values[observed]=y[observed]-image[observed]
    residual,_=predict_anchors(xy,residual_values,observed,hidden,config)
    return {'image':image[hidden], 'anchor_idw_k16_p2':anchor,
            'fixed_20pct_image':.2*image[hidden]+.8*anchor, 'exploratory_residual_anchor':image[hidden]+residual}


def evaluate(root,data_root,output,checkpoint_root):
    from .journal import provenance
    from .core import mlp_model,predict
    data=load_prepared(data_root); records=json.loads((output/'FIT_MANIFEST.json').read_text())
    scaler=Scaler()
    with np.load(checkpoint_root/'journal_v2_feature_scaler.npz') as f:scaler.mean=f['mean'];scaler.scale=f['scale']
    scaled={s:scaler.transform(data[s]['embeddings']) for s in SAMPLES}
    graphs={s:build_graph(data[s]['xy_um'],{'kind':'delaunay','max_edge_um':250}) for s in SAMPLES}
    torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
    masks=[];sample=data['TENX99'];n=len(sample['targets'])
    for budget in [.1,.2,.4]:
        for replicate in range(3):
            p=root/'checkpoints/legacy/reports/C0/C0_MASKS'/f'TENX99_random_q{budget:.2f}_r{replicate}.npz'
            with np.load(p) as f:
                np.testing.assert_array_equal(f['barcodes'],sample['barcodes']);observed=np.flatnonzero(f['observed']);hidden=np.flatnonzero(~f['observed'])
            masks.append({'budget':budget,'replicate':replicate,'observed':observed,'hidden':hidden,'mask_sha256':sha(p)})
    panel=[];genes=[];predictions={};checks=0;start=time.monotonic()
    for record in records:
        run_id=record['run_id'];checkpoint=checkpoint_root/(run_id+'.pt')
        if sha(checkpoint)!=record['checkpoint_sha256']:raise ValueError('Checkpoint mismatch at evaluation')
        x=representation(data,scaled,graphs,record['model'],record['control'],record['seed'])['TENX99']
        net=mlp_model(x.shape[1],len(sample['genes']),record['config'])
        state=torch.load(checkpoint,map_location='cpu',weights_only=True);net.load_state_dict(state['state_dict']);net.eval()
        image=predict(net,x,'B2').astype(float);predictions[(record['model'],record['control'],record['seed'])]=image
        for mask in masks:
            obs=mask['observed'];hid=mask['hidden'];base={k:v for k,v in mask.items() if k not in ['observed','hidden']}
            y=sample['targets'].astype(float);pred=hybrids(sample['xy_um'],y,image,obs,hid)
            for mode in ['nan','random']:
                corrupt=y.copy();corrupt[hid]=np.nan if mode=='nan' else np.random.default_rng(19).normal(size=(len(hid),y.shape[1]))
                alt=hybrids(sample['xy_um'],corrupt,image,obs,hid)
                for route in ['fixed_20pct_image','exploratory_residual_anchor']:
                    np.testing.assert_array_equal(pred[route],alt[route]);checks+=1
            for route,yp in pred.items():
                m=gene_metrics(y[hid],yp);attributes={**base,'run_id':run_id,'model':record['model'],'control':record['control'],'training_seed':record['seed'],'route':route}
                panel.append({**attributes,'n_observed':len(obs),'n_hidden':len(hid),**summarize(m)})
                for j,gene in enumerate(sample['genes']):genes.append({**attributes,'gene':gene,**{k:float(m[k][j]) for k in ['pearson','spearman','mse','mae']}})
        print('EVALUATED',run_id,flush=True)
    # Identity controls must equal the matching trained local model, not count as independent evidence.
    identity_checks=[]
    for seed in [0,1,2]:
        for model in ['MEAN','GCN']:
            np.testing.assert_array_equal(predictions[('local','local',seed)],predictions[(model,'self',seed)])
            identity_checks.append({'seed':seed,'control':model+'_self','bitwise_equal_local':True})
    table=pd.DataFrame(panel);table.to_csv(output/'EVALUATION_RUN_METRICS.csv',index=False)
    pd.DataFrame(genes).to_csv(output/'EVALUATION_GENE_METRICS.csv.gz',index=False)
    means=table.groupby(['model','control','training_seed','budget','route'])[['median_pearson','mean_mse','mean_mae']].mean().reset_index()
    means.to_csv(output/'SEED_METRICS.csv',index=False)
    seeds=means.groupby(['model','control','budget','route'])[['median_pearson','mean_mse','mean_mae']].agg(['mean','std','min','max'])
    seeds.columns=['_'.join(c) for c in seeds.columns];seeds.reset_index().to_csv(output/'TRAINING_SEED_SUMMARY.csv',index=False)
    result=provenance(root,status='COMPLETE_DEVELOPMENT_ONLY',protocol_E_commit='7e3cc167f',fit_models=len(records),n_evaluation_patients=1,
                      n_evaluation_rows=len(panel),hidden_corruption_checks=checks,identity_checks=identity_checks,elapsed_seconds=time.monotonic()-start,
                      warning='TENX99 was already used in v1 development; masks/training seeds/genes are not biological replicates. No model or hybrid selected from outcomes.')
    (output/'EVALUATION_STATUS.json').write_text(json.dumps(result,indent=2)+'\n');return result
