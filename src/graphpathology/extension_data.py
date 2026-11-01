"""Verified external image caches and immutable extension masks; no refitting on Wu."""
import json,time
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import torch
from PIL import Image
from torchvision.models import resnet18,ResNet18_Weights
from .assets import sha
from .journal_data import load_patient,WU_EVAL
from .extension_train import metadata


def prepare_wu(root):
 out=root/'data/journal_revision/transformer_egnn/wu';out.mkdir(parents=True,exist_ok=True)
 weight=root/'data/journal_revision/hest_training/models/resnet18-f37072fd.pth';digest=sha(weight)
 expected=next(f['sha256'] for f in json.loads((root/'experiments/journal_revision/ACQUIRED_INPUTS.json').read_text())['files'] if 'resnet18-f37072fd.pth' in f['url'])
 if digest!=expected:raise ValueError('Unverified encoder')
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True);net=resnet18(weights=None);net.load_state_dict(torch.load(weight,weights_only=True));net.fc=torch.nn.Identity();net.eval();transform=ResNet18_Weights.IMAGENET1K_V1.transforms()
 records=[]
 for patient in WU_EVAL:
  dest=out/(patient+'.npz');meta=dest.with_suffix('.json')
  if dest.exists() and meta.exists():
   r=json.loads(meta.read_text());assert sha(dest)==r['cache_sha256'];records.append(r);continue
  d=load_patient(root,root/'data/journal_revision','wu',patient);n=len(d['barcodes']);base=root/'data/journal_revision/raw/wu/extracted/spatial'/(patient+'_spatial')
  positions=pd.read_csv(base/'tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','y','x']).set_index('barcode').loc[d['barcodes']]
  xy=positions[['x','y']].to_numpy(float);sf=json.loads((base/'scalefactors_json.json').read_text());hires=sf['tissue_hires_scalef']
  lattice=np.column_stack([positions.array_col.to_numpy()*50,positions.array_row.to_numpy()*np.sqrt(3)*50]);design=np.column_stack([lattice,np.ones(n)])
  coef=np.linalg.lstsq(design,xy,rcond=None)[0];mpp=1/float(np.mean(np.linalg.svd(coef[:2],compute_uv=False)))
  np.testing.assert_allclose(xy*mpp,d['xy'],rtol=0,atol=1e-8)
  width=int(np.floor(112/(mpp/hires)+.5))
  legacy=root/'checkpoints/legacy/C9/C9A/logs'/(patient+'_EXTRACTION.json')
  if legacy.exists():
   original=json.loads(legacy.read_text());assert width==original['native_patch_pixels'] and sha(base/'tissue_hires_image.png')==original['source_image_sha256']
  image=Image.open(base/'tissue_hires_image.png').convert('RGB');vectors=[];start=time.monotonic();patch_hashes=[]
  for begin in range(0,n,64):
   patches=[]
   for x,y in xy[begin:begin+64]*hires:
    left=int(np.floor(x-width/2+.5));top=int(np.floor(y-width/2+.5));patch=image.crop((left,top,left+width,top+width)).resize((224,224),Image.Resampling.LANCZOS)
    patch_hashes.append(__import__('hashlib').sha256(np.asarray(patch).tobytes()).hexdigest());patches.append(transform(patch))
   with torch.inference_mode():vectors.append(net(torch.stack(patches)).numpy())
   if begin%1280==0:print('WU FEATURES',patient,begin,'/',n,flush=True)
  embedding=np.concatenate(vectors)
  np.savez_compressed(dest,embeddings=embedding,targets=d['targets'],xy_um=d['xy'],barcodes=d['barcodes'],genes=d['genes'],patch_sha256=np.array(patch_hashes))
  r=metadata(root,status='COMPLETE',patient=patient,cache_sha256=sha(dest),encoder_sha256=digest,source_hashes=d['source_hashes'],source_image_sha256=sha(base/'tissue_hires_image.png'),source_positions_sha256=sha(base/'tissue_positions_list.csv'),
    native_patch_pixels=width,hires_mpp=mpp/hires,nominal_FOV_um=112,actual_FOV_um=width*mpp/hires,elapsed_seconds=time.monotonic()-start,
    limitations='15-25 native pixels per112um; previously inspected patients; no full-image rotation-equivariance claim')
  meta.write_text(json.dumps(r,indent=2)+'\n');records.append(r);print('WU FEATURES COMPLETE',patient,flush=True)
 report=root/'experiments/journal_revision/transformer_egnn/WU_FEATURE_MANIFEST.json';report.write_text(json.dumps(records,indent=2)+'\n');return records
