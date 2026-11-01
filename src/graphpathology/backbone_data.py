"""Frozen Phikon embeddings on exactly the existing registered RGB patches."""
import datetime
import hashlib
import json
import subprocess
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import requests
import torch
from PIL import Image

from .assets import sha
from .core import PeakRSS
from .journal_learned import load_prepared

STUDY='backbone_graph_comparison'
PROTOCOL_COMMIT='60842a6b5'


def config(root):
    return json.loads((root/f'configs/journal_revision/{STUDY}/protocol.json').read_text())


def record(root, **fields):
    return {'classification':'exploratory_follow_up','protocol_commit':PROTOCOL_COMMIT,
            'code_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
            'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'code_hashes':{str(p.relative_to(root)):sha(p) for p in (root/'src/graphpathology').glob('backbone*.py')}, **fields}


def write(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n')


def source_path(root, patient):
    from .extension_evaluate import feature_path
    if patient in ['NCBI783','NCBI785','TENX95']:
        return root/'data/journal_revision/hest_training/processed/journal_v2_hest'/(patient+'.npz')
    return feature_path(root,patient)


def source(root, patient):
    path=source_path(root,patient)
    meta=json.loads(path.with_suffix('.json').read_text())
    assert sha(path)==meta['cache_sha256']
    with np.load(path) as f:
        data={k:f[k].copy() for k in f.files}
    assert len(set(data['barcodes']))==len(data['barcodes'])
    return data,meta


def acquire(root):
    cfg=config(root);dest=root/f'data/journal_revision/{STUDY}/encoder'
    expected=json.loads((root/'checkpoints/legacy/C4/C4A/configs/encoder_files.json').read_text())
    rows=[]
    for r in expected['files']:
        path=dest/r['name'];url=f"https://huggingface.co/owkin/phikon/resolve/{cfg['encoder_revision']}/{r['name']}"
        if not path.exists():
            path.parent.mkdir(parents=True,exist_ok=True)
            with requests.get(url,stream=True,timeout=(30,180)) as response:
                response.raise_for_status()
                with path.with_suffix(path.suffix+'.partial').open('wb') as f:
                    for chunk in response.iter_content(1024*1024):f.write(chunk)
            path.with_suffix(path.suffix+'.partial').replace(path)
        assert sha(path)==r['sha256'] and path.stat().st_size==r['bytes'],r['name']
        rows.append({'path':str(path.relative_to(root)),'url':url,'bytes':path.stat().st_size,'sha256':sha(path)})
    license_url='https://raw.githubusercontent.com/owkin/HistoSSLscaling/main/LICENSE.txt'
    license_path=dest/'LICENSE.txt'
    if not license_path.exists():
        response=requests.get(license_url,timeout=30);response.raise_for_status();license_path.write_bytes(response.content)
    report=record(root,status='COMPLETE',files=rows,license_url=license_url,license_sha256=sha(license_path),
                  use='Author-authorized noncommercial research; ungated public files; no access terms accepted or encoder redistributed.')
    write(root/f'experiments/journal_revision/{STUDY}/ENCODER_ACQUISITION.json',report)
    return report


def patches(root, patient, data, meta, begin, end):
    """Return the exact pre-encoder RGB arrays already used by the ResNet study."""
    if patient in ['NCBI783','NCBI785','TENX95','TENX99']:
        from .data import read_slide
        raw=root/'data/journal_revision/hest_training/raw/hest'
        md=pd.read_csv(raw/'HEST_v1_1_0.csv').set_index('id')
        slide=read_slide(raw,patient,float(md.loc[patient,'pixel_size_um_estimated']))
        np.testing.assert_array_equal(data['barcodes'],slide.barcodes)
        np.testing.assert_array_equal(data['xy_um'],slide.xy_um)
        assert sha(slide.patch_path)==meta['source_patches_sha256']
        # This context is entered once by the streaming extractor below.
        return slide.patch_path
    from .journal_data import WU_EVAL
    if patient in WU_EVAL:
        base=root/'data/journal_revision/raw/wu/extracted/spatial'/(patient+'_spatial')
        image_path=base/'tissue_hires_image.png';position_path=base/'tissue_positions_list.csv';sf=base/'scalefactors_json.json'
    else:
        acquisition=json.loads((root/'experiments/journal_revision/transformer_egnn/GEO_SELECTIVE_ACQUISITION.json').read_text())
        rows=[r for r in acquisition['selected_files'] if r['patient_id']==patient]
        raw=root/'data/journal_revision/transformer_egnn/geo_ranges/raw/bassiouni'
        def locate(term):return raw/Path(next(r['extracted_path'] for r in rows if term in r['path'])).name
        image_path=locate('tissue_hires_image');position_path=locate('tissue_positions');sf=locate('scalefactors')
    assert sha(image_path)==meta['source_image_sha256'] and sha(position_path)==meta['source_positions_sha256']
    positions=pd.read_csv(position_path,header=None,names=['barcode','in_tissue','array_row','array_col','y','x']).set_index('barcode').loc[data['barcodes']]
    factor=json.loads(sf.read_text())['tissue_hires_scalef']
    width=meta['native_patch_pixels'];image=Image.open(image_path).convert('RGB')
    return image,positions[['x','y']].to_numpy(float)*factor,width


def prepare(root):
    from transformers import ViTModel,ViTImageProcessor
    cfg=config(root);acquire(root)
    torch.set_num_threads(cfg['encoder_threads']);torch.use_deterministic_algorithms(True)
    path=root/f'data/journal_revision/{STUDY}/encoder'
    processor=ViTImageProcessor.from_pretrained(path,local_files_only=True)
    net,loading=ViTModel.from_pretrained(path,add_pooling_layer=False,attn_implementation='sdpa',local_files_only=True,output_loading_info=True)
    assert not loading['missing_keys'] and not loading.get('mismatched_keys') and not loading.get('error_msgs'),loading
    assert set(loading['unexpected_keys'])<= {'pooler.dense.weight','pooler.dense.bias','embeddings.mask_token'},loading
    net.eval()
    write(root/f'experiments/journal_revision/{STUDY}/ENCODER_LOADING.json',record(root,status='COMPLETE',loading=loading))
    for parameter in net.parameters():parameter.requires_grad_(False)
    out=root/f'data/journal_revision/{STUDY}/phikon';out.mkdir(parents=True,exist_ok=True)
    samples=['NCBI783','NCBI785','TENX95','TENX99']+cfg['external']
    manifest=[];started=time.monotonic();batch_size=cfg['encoder_batch_size']
    total=sum(len(source(root,s)[0]['barcodes']) for s in samples)
    for patient in samples:
        target=out/(patient+'.npz');done=target.with_suffix('.json')
        if target.exists() and done.exists():
            r=json.loads(done.read_text());assert sha(target)==r['cache_sha256'];manifest.append(r);print('PHIKON RESUMED',patient,flush=True);continue
        data,meta=source(root,patient);n=len(data['barcodes']);npy=out/(patient+'.partial.npy')
        state=out/(patient+'.progress.json');offset=0
        if npy.exists() and state.exists():
            r=json.loads(state.read_text());assert r['source_cache_sha256']==meta['cache_sha256'];offset=r['completed_rows'];x=np.lib.format.open_memmap(npy,mode='r+')
        else:x=np.lib.format.open_memmap(npy,mode='w+',dtype=np.float32,shape=(n,768))
        inputs=patches(root,patient,data,meta,0,n);is_hest=isinstance(inputs,Path)
        handle=h5py.File(inputs,'r') if is_hest else None
        begin_time=time.monotonic();timed_rows=0;preflight_path=root/f'experiments/journal_revision/{STUDY}/ENCODER_PREFLIGHT.json'
        preflight_done=preflight_path.exists()
        if preflight_done:assert json.loads(preflight_path.read_text())['status']=='COMPLETE'
        with PeakRSS() as memory:
            try:
                for begin in range(offset,n,batch_size):
                    end=min(n,begin+batch_size)
                    if is_hest:arrays=handle['img'][begin:end]
                    else:
                        image,xy,width=inputs;arrays=[]
                        for xx,yy in xy[begin:end]:
                            left=int(np.floor(xx-width/2+.5));top=int(np.floor(yy-width/2+.5))
                            arrays.append(np.asarray(image.crop((left,top,left+width,top+width)).resize((224,224),Image.Resampling.LANCZOS)))
                    if 'patch_sha256' in data:
                        hashes=np.array([hashlib.sha256(np.asarray(a).tobytes()).hexdigest() for a in arrays])
                        np.testing.assert_array_equal(hashes,data['patch_sha256'][begin:end])
                    with torch.inference_mode():
                        values=net(**processor(images=[Image.fromarray(a) for a in arrays],return_tensors='pt')).last_hidden_state[:,0,:].numpy()
                    assert values.shape==(end-begin,768) and np.isfinite(values).all()
                    x[begin:end]=values;x.flush();timed_rows+=end-begin
                    write(state,{'completed_rows':end,'source_cache_sha256':meta['cache_sha256']})
                    if not preflight_done and timed_rows>=64:
                        estimate=(time.monotonic()-begin_time)/timed_rows*total/3600
                        report=record(root,status='COMPLETE' if estimate<=cfg['encoder_max_wall_hours'] and memory.peak/2**30<cfg['max_rss_gib'] else 'BLOCKED_RESOURCES',source_patient=patient,timed_rows=timed_rows,total_rows=total,estimated_hours=estimate,peak_rss_gib=memory.peak/2**30,threads=cfg['encoder_threads'])
                        write(root/f'experiments/journal_revision/{STUDY}/ENCODER_PREFLIGHT.json',report)
                        if report['status']!='COMPLETE':raise RuntimeError('Predeclared encoder resource ceiling exceeded')
                        preflight_done=True
                    if end%256==0 or end==n:print('PHIKON',patient,end,'/',n,'seconds',round(time.monotonic()-begin_time,1),flush=True)
                    if (time.monotonic()-started)/3600>cfg['encoder_max_wall_hours']:raise RuntimeError('Encoder wall-hour ceiling exceeded; resumable partial cache retained')
                    if memory.peak/2**30>cfg['max_rss_gib']:raise MemoryError('Predeclared RSS limit exceeded')
            finally:
                if handle:handle.close()
        keep={k:v for k,v in data.items() if k!='embeddings'}
        np.savez_compressed(target,embeddings=np.asarray(x),**keep)
        r=record(root,status='COMPLETE',patient=patient,encoder_revision=cfg['encoder_revision'],encoder_sha256=cfg['encoder_sha256'],source_cache_sha256=meta['cache_sha256'],cache_sha256=sha(target),rows=n,embedding_dimension=768,elapsed_seconds_this_resume=time.monotonic()-begin_time,resumed_rows=offset,peak_rss_gib=memory.peak/2**30,input_contract='Exact registered RGB arrays; official processor; frozen final-LayerNorm CLS without pooler',physical_field_of_view=meta.get('actual_FOV_um','Archived HEST C0 patch geometry; nominal112um'))
        write(done,r);manifest.append(r);npy.unlink();state.unlink();print('PHIKON COMPLETE',patient,flush=True)
        write(root/f'experiments/journal_revision/{STUDY}/FEATURE_MANIFEST.json',manifest)
    return manifest


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['acquire','prepare']);parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);a=parser.parse_args()
    print(json.dumps((acquire if a.action=='acquire' else prepare)(a.root.resolve()),indent=2))
