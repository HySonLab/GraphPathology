"""Fetch the pinned, public Phikon release; no gated-model workaround."""
import json,hashlib,time
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C4';rev='057cc0295895c2df3dd7681a89680da6015cbefe';folder=OUT/'weights/phikon';folder.mkdir(exist_ok=True)
metadata=OUT/'C4A/configs/encoder_files.json';existing=json.loads(metadata.read_text()) if metadata.exists() else None
expected={r['name']:r for r in existing['files']} if existing else {}
records=[]
for name in ['config.json','preprocessor_config.json','README.md','model.safetensors']:
 p=folder/name;t=time.perf_counter()
 if not p.exists():
  with requests.get(f'https://huggingface.co/owkin/phikon/resolve/{rev}/{name}',stream=True,timeout=120) as r:
   r.raise_for_status()
   with p.open('wb') as f:
    for block in r.iter_content(1024*1024):f.write(block)
 with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 if existing:assert digest==expected[name]['sha256'],f'Pinned weight/config hash mismatch: {name}'
 records.append({'name':name,'sha256':digest,'bytes':p.stat().st_size,'seconds':time.perf_counter()-t});print(name,p.stat().st_size,flush=True)
if not existing:metadata.write_text(json.dumps({'model':'owkin/phikon','revision':rev,'source':'https://huggingface.co/owkin/phikon','files':records},indent=2)+'\n')
