"""Recover the three published HER2ST archives without reading expression values.
The URLs/checksums are the source metadata captured during qualification.
"""
from pathlib import Path
import json,hashlib,requests
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C7'
def main():
 metadata=json.loads((OUT/'raw/HER2ST_metadata.json').read_text())
 for item in metadata['files']:
  if item['key'] not in ['images.zip','count-matrices.zip','spot-selections.zip']:continue
  path=OUT/'raw'/item['key']
  if not path.exists():
   url=item['links'].get('self') or item['links']['download']
   with requests.get(url,stream=True,timeout=(30,180)) as response:
    response.raise_for_status()
    with path.with_suffix('.zip.part').open('wb') as f:
     for chunk in response.iter_content(2**20):f.write(chunk)
   path.with_suffix('.zip.part').replace(path)
  with path.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==item['checksum'].split(':')[-1]
  print('Verified source archive:',path.name,flush=True)
if __name__=='__main__':main()
