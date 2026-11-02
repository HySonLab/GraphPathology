"""Auditable driver for a fresh checkout; imports only the tagged checkout's code.

The probe uses the existing pinned dependency environment and two registered real
inputs read locally. It does not claim independent hardware or fresh data acquisition.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--checkout',type=Path,required=True)
p.add_argument('--source-inputs',type=Path,required=True)
p.add_argument('--downloads',type=Path,required=True)
p.add_argument('--encoder',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();root=a.checkout.resolve();source=a.source_inputs.resolve()
if a.output.exists():raise FileExistsError(a.output)
sys.path.insert(0,str(root/'src'))
import graphpathology
assert Path(graphpathology.__file__).resolve().is_relative_to(root/'src')
import h5py
import numpy as np
from PIL import Image
import torch
from torchvision.models import resnet18,ResNet18_Weights
from graphpathology.assets import sha
from graphpathology.reproduction_primary import image_predictions
from graphpathology.reproduction_release import restore,preservation
from graphpathology.extension_evaluate import model_inventory,load_data,prediction
from graphpathology.reproduction_artifacts import matched_prediction

started=time.monotonic();torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
manifest=json.loads((a.downloads/'graphpathology-model-and-results-manifest-v1.json').read_text())
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
assert commit==manifest['code_commit']=='44918ef8f9e9876234ceb02a8b4a815d8b4072d4'
restored=restore(root,a.downloads,manifest)
expected_env=json.loads((root/'reproduction/ENVIRONMENT.json').read_text())
for package,version in expected_env['packages'].items():assert importlib.metadata.version(package)==version,(package,version)
# Every restored public PT file, including optimizer/resume snapshots, is loaded safely.
pt={name for asset in manifest['assets'] for name in asset.get('files',{}) if name.endswith('.pt')}
assert len(pt)==147
for name in sorted(pt):torch.load(root/name,map_location='cpu',weights_only=True)

# This real554-location cache is copied locally, not downloaded or redistributed.
patient='BASS_P19';relative=Path('data/journal_revision/transformer_egnn/bassiouni')/(patient+'.npz')
for path in [relative,relative.with_suffix('.json')]:
 target=root/path;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/path,target)
data=load_data(root,patient);assert len(data['barcodes'])==554
seen=set();forwards=[]
for item in model_inventory(root):
 image=prediction(root,data,patient,item)
 assert image.shape==(554,50) and np.isfinite(image).all()
 seen.add(item['checkpoint']);forwards.append({'checkpoint':item['checkpoint'],'rows':554,'finite':True})
items=json.loads((root/'experiments/journal_revision/backbone_graph_comparison/ALL_CANDIDATE_FIT_MANIFEST.json').read_text())
for item in items:
 if item['backbone']!='resnet18' or item['checkpoint'] in seen:continue
 _,image=matched_prediction(root,patient,item)
 assert image.shape==(554,50) and np.isfinite(image).all()
 seen.add(item['checkpoint']);forwards.append({'checkpoint':item['checkpoint'],'rows':554,'finite':True})
assert len(seen)==96

# Fresh official encoder bytes were downloaded into the clean-room workspace.
encoder_sha=sha(a.encoder)
assert encoder_sha=='f37072fd47e89c5e827621c5baffa7500819f7896bbacec160b1a16c560e07ec'
net=resnet18(weights=None);net.load_state_dict(torch.load(a.encoder,map_location='cpu',weights_only=True));net.fc=torch.nn.Identity();net.eval();net.requires_grad_(False)
base=source/'data/journal_revision/hest_training'
cache=base/'processed/journal_v2_hest/NCBI783.npz'
meta=json.loads(cache.with_suffix('.json').read_text())
patch=base/'raw/hest/IDC/patches/NCBI783.h5';assert sha(patch)==meta['source_patches_sha256']
with h5py.File(patch) as f:
 rgb=f['img'][:8];barcodes=f['barcode'][:8].reshape(-1).astype(str)
transform=ResNet18_Weights.IMAGENET1K_V1.transforms()
batch=torch.stack([transform(Image.fromarray(im)) for im in rgb])
with torch.inference_mode():features=net(batch).numpy()
with np.load(cache,allow_pickle=False) as f:
 reference=f['embeddings'][:8];np.testing.assert_array_equal(barcodes,f['barcodes'][:8])
# Locked before this probe: tolerance for float32 batch-shape convolution reductions.
np.testing.assert_allclose(features,reference,atol=1e-5,rtol=1e-5)
folder=root/'checkpoints/reproduced/primary_resnet18/primary-20261101-v1'
image,_=image_predictions(folder,features)
assert image.shape==(8,50) and image.dtype==np.float32 and np.isfinite(image).all()
env={**os.environ,'PYTHONPATH':str(root/'src')}
doctor=subprocess.run([sys.executable,'-m','graphpathology','--root',str(root),'doctor','--require-models'],cwd=root,env=env,capture_output=True,text=True)
assert doctor.returncode==1
state=json.loads(doctor.stdout);assert len(state['primary_heads'])==3 and not any(r['available'] for r in state['primary_heads'])
protected=preservation(root);assert protected['status']=='PASS'
subprocess.run(['git','diff','--exit-code','HEAD','--'],cwd=root,check=True,capture_output=True)
result={'status':'PASS','code_commit':commit,'remote':'https://github.com/HySonLab/GraphPathology.git',
        'code_imported_from_clean_checkout':True,'dependency_environment':'Existing pinned venv reused on the same machine; not an independent hardware/environment replication.',
        'restore':restored,'safe_public_pt_files':len(pt),'recovered_unique_heads_forwarded':len(seen),'recovered_forward_input':'Actual BASS_P19,554 registered locations; cache copied locally and hash-verified.',
        'encoder_url':'https://download.pytorch.org/models/resnet18-f37072fd.pth','encoder_bytes':a.encoder.stat().st_size,'encoder_sha256':encoder_sha,
        'primary_real_image_rows':8,'primary_output_shape':list(image.shape),'primary_output_finite':True,
        'primary_output_array_sha256':hashlib.sha256(image.tobytes()).hexdigest(),
        'fresh_encoder_feature_max_abs_delta':float(np.max(abs(features-reference))),
        'fresh_encoder_feature_tolerance':{'atol':1e-5,'rtol':1e-5},'real_input_patch_sha256':sha(patch),
        'original_doctor_expected_missing_exit_code':doctor.returncode,'preservation':protected,
        'tracked_checkout_files_unchanged':True,'probe_script_sha256':sha(Path(__file__)),
        'elapsed_seconds':time.monotonic()-started,'forwards':forwards}
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['forwards','preservation']},indent=2))
