"""Public GEO acquisition; never loads expression values during qualification."""
from pathlib import Path
import requests,re,json,hashlib,concurrent.futures,argparse,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C7';RAW=OUT/'raw'
IDS=['GSM9509172','GSM9509160','GSM9509170','GSM9509165']
def fetch(url):
 p=RAW/url.rsplit('/',1)[-1];t=time.time()
 if not p.exists():
  with requests.get(url,stream=True,timeout=(30,180)) as r:
   r.raise_for_status()
   with p.with_suffix(p.suffix+'.part').open('wb') as f:
    for b in r.iter_content(2**20):f.write(b)
  p.with_suffix(p.suffix+'.part').replace(p)
 with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
 row={'url':url,'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':h,'seconds':time.time()-t}
 q=OUT/'C7A/logs/downloads';q.mkdir(exist_ok=True);(q/(p.name+'.json')).write_text(json.dumps(row,indent=2)+'\n');print(p.name,p.stat().st_size,flush=True)
 return row
def main():
 a=argparse.ArgumentParser();a.add_argument('stage',choices=['qualification','targets']);args=a.parse_args();urls=[]
 for s in IDS:
  text=(RAW/f'{s}.txt').read_text()
  for u in re.findall(r'!Sample_supplementary_file = (\S+)',text):
   suffix=['cell_feature_matrix.h5','morphology_focus.ome.tif.gz'] if args.stage=='qualification' else ['transcripts.parquet.gz']
   if any(u.endswith(x) for x in suffix):urls.append(u.replace('ftp://','https://'))
 if args.stage=='qualification':
  text=(RAW/'GSE311609.txt').read_text()
  urls.extend(u.replace('ftp://','https://') for u in re.findall(r'!Series_supplementary_file = (\S+)',text) if '_hne_breast_' in u and u.endswith('.ndpi'))
 else:raise RuntimeError('GSE311609 was not qualified; target acquisition is disabled for C7')
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(fetch,urls))
if __name__=='__main__':main()
