"""Read public ZIP members with HTTP ranges; qualification never reads counts."""
from pathlib import Path
import requests,io,zipfile,json,hashlib,argparse,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C7';RAW=OUT/'raw'
BUNDLE='https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/3.0.0/Xenium_Prime_Breast_Cancer_FFPE/Xenium_Prime_Breast_Cancer_FFPE_outs.zip'
BASE='https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_Breast_Cancer_FFPE/Xenium_Prime_Breast_Cancer_FFPE'
class Remote(io.RawIOBase):
 def __init__(self,url):
  self.url=url;self.pos=0;self.s=requests.Session();r=self.s.head(url,timeout=60);r.raise_for_status();self.size=int(r.headers['Content-Length'])
 def seek(self,offset,whence=0):
  self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset;return self.pos
 def tell(self):return self.pos
 def read(self,n=-1):
  n=min(self.size-self.pos,n if n>=0 else self.size-self.pos)
  if n<=0:return b''
  r=self.s.get(self.url,headers={'Range':f'bytes={self.pos}-{self.pos+n-1}'},timeout=(30,180));r.raise_for_status();assert r.status_code==206,(r.status_code,self.pos,n);assert len(r.content)==n;self.pos+=n;return r.content
 def seekable(self):return True
 def readable(self):return True

def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['metadata','image','targets']);args=ap.parse_args()
 if args.stage=='image':
  from download import fetch
  fetch(BASE+'_he_image.ome.tif');fetch(BASE+'_he_imagealignment.csv');return
 if args.stage=='targets':
  lock=json.loads((OUT/'C7A/C7A_SPEC_LOCK.json').read_text());assert lock['external_evaluation_authorized']
 with zipfile.ZipFile(Remote(BUNDLE)) as z:
  entries=[{'name':i.filename,'bytes':i.file_size,'compressed_bytes':i.compress_size,'crc':i.CRC} for i in z.infolist()];(RAW/'prime_bundle_inventory.json').write_text(json.dumps(entries,indent=2)+'\n')
  print('members',len(entries));print([e for e in entries if any(e['name'].endswith(s) for s in ['gene_panel.json','experiment.xenium','transcripts.parquet'])],flush=True)
  names=[i.filename for i in z.infolist() if i.filename.endswith(('gene_panel.json','experiment.xenium'))] if args.stage=='metadata' else [i.filename for i in z.infolist() if i.filename.endswith('transcripts.parquet')]
  for name in names:
   dest=RAW/('prime_'+Path(name).name)
   if not dest.exists():
    t=time.time()
    with z.open(name) as src,dest.with_suffix(dest.suffix+'.part').open('wb') as dst:
     while b:=src.read(2**22):dst.write(b)
    dest.with_suffix(dest.suffix+'.part').replace(dest);print(name,dest.stat().st_size,time.time()-t,flush=True)
   with dest.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
   (OUT/'C7A/logs/downloads'/('prime_'+Path(name).name+'.json')).write_text(json.dumps({'url':BUNDLE,'member':name,'path':str(dest.relative_to(ROOT)),'sha256':h,'bytes':dest.stat().st_size},indent=2)+'\n')
if __name__=='__main__':main()
