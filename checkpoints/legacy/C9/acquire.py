"""Acquire public source archives; qualification does not read expression values."""
from pathlib import Path
import requests,hashlib,json,tarfile
R=Path(__file__).resolve().parent
FILES={'metadata.tar.gz':'1d0b34ead70635c094ac74ac58a88d68','spatial.tar.gz':'29b02ae433d459d19dc52b75a529e8e9','filtered_count_matrices.tar.gz':'ea1220f0606c4d4e9307468fb2d5427a'}
rows=[]
for name,expected in FILES.items():
 p=R/'raw'/name;url='https://zenodo.org/records/4739739/files/'+name+'?download=1'
 if not p.exists():
  with requests.get(url,stream=True,timeout=120) as r:
   r.raise_for_status()
   with p.open('wb') as f:
    for b in r.iter_content(1024*1024):f.write(b)
 with p.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==expected
 with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
 rows.append({'url':url,'file':name,'source_md5':expected,'sha256':h,'bytes':p.stat().st_size})
 with tarfile.open(p) as tf:
  tf.extractall(R/'raw/extracted',filter='data');print(name,tf.getnames()[:15],flush=True)
(R/'C9A/logs/DOWNLOADS.json').write_text(json.dumps(rows,indent=2)+'\n')
