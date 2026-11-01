"""Metadata-first acquisition of one first-numbered section per fresh TNBC patient."""
from pathlib import Path
import requests,re,json,hashlib,gzip
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C10';LOG=OUT/'C10A/logs';RAW=OUT/'raw'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 p=LOG/'GSE210616_SAMPLES.soft'
 if not p.exists():
  r=requests.get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi',params={'acc':'GSE210616','targ':'gsm','form':'text','view':'full'},timeout=120);r.raise_for_status();p.write_text(r.text)
 samples=[];allrows=[]
 for block in p.read_text().split('^SAMPLE = ')[1:]:
  gsm=block.splitlines()[0].strip();title=re.search(r'!Sample_title = (.*)',block).group(1).strip();m=re.search(r'patient (\d+), section (\d+)',title);assert m,(gsm,title)
  patient,section=map(int,m.groups());urls=re.findall(r'!Sample_supplementary_file_\d+ = (.*)',block);urls=[u.strip().replace('ftp://','https://') for u in urls if '.cloupe' not in u];row={'gsm':gsm,'patient_id':f'BASS_P{patient:02}','section':section,'title':title,'urls':urls};allrows.append(row)
  if section==1:samples.append(row)
 assert len(samples)==22 and len({r['patient_id'] for r in samples})==22
 samples=sorted(samples,key=lambda r:r['patient_id']);(LOG/'SOURCE_SAMPLE_SELECTION.json').write_text(json.dumps({'rule':'All22patients, first numbered section only, fixed from GEO titles before outcome/feature matrix access','selected':samples,'all_sections':allrows},indent=2)+'\n')
 def fetch(task):
  row,url=task;name=url.split('/')[-1];p=RAW/name
  if not p.exists():
   tmp=p.with_suffix(p.suffix+'.part')
   with requests.get(url,stream=True,timeout=120) as r:
    r.raise_for_status()
    with tmp.open('wb') as f:
     for b in r.iter_content(1024*1024):f.write(b)
   tmp.rename(p)
  rec={'patient_id':row['patient_id'],'gsm':row['gsm'],'url':url,'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)}
  if name.endswith('.gz'):
   dst=RAW/name[:-3]
   if not dst.exists():
    with gzip.open(p,'rb') as f:dst.write_bytes(f.read())
   rec['extracted_path']=str(dst.relative_to(ROOT));rec['extracted_sha256']=sha(dst)
  print('ACQUIRED',row['patient_id'],name,flush=True);return rec
 with ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fetch,[(r,u) for r in samples for u in r['urls']]))
 (LOG/'DOWNLOADS.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
