from pathlib import Path
import requests,hashlib,json
p=Path(__file__).resolve().parent/'raw/images.pdf';url='https://zenodo.org/records/4739739/files/images.pdf?download=1'
if not p.exists():
 with requests.get(url,stream=True,timeout=120) as r:
  r.raise_for_status()
  with p.open('wb') as f:
   for b in r.iter_content(1024*1024):f.write(b)
with p.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()=='7a44782b0c376e5b4cb8bb1127821e6f'
with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
(p.parents[1]/'C9A/logs/IMAGE_DOWNLOAD.json').write_text(json.dumps({'url':url,'source_md5':'7a44782b0c376e5b4cb8bb1127821e6f','sha256':h,'bytes':p.stat().st_size},indent=2)+'\n')
