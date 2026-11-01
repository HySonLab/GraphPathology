"""Selective public GEO TAR ranges; never download the 37.7GB whole-series bundle."""
import datetime,hashlib,json,math,shutil,time
from pathlib import Path
import requests
from .assets import sha

URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210616/suppl/GSE210616_RAW.tar'
TOTAL=37690060800
CAP=1024**3

def acquire(root):
 report=root/'experiments/journal_revision/transformer_egnn/GEO_SELECTIVE_ACQUISITION.json'
 if report.exists():
  previous=json.loads(report.read_text())
  if previous['status']=='COMPLETE':
   complete=True
   for record in previous['selected_files']:
    path=root/record['cache_relative']
    if path.exists() and sha(path)!=record['sha256']:raise ValueError('Completed acquisition cache changed')
    if not path.exists():complete=False
   if complete:return previous
 out=root/'data/journal_revision/transformer_egnn/geo_ranges';(out/'raw/bassiouni').mkdir(parents=True,exist_ok=True)
 original=json.loads((root/'checkpoints/legacy/C10/C10A/logs/DOWNLOADS.json').read_text());wanted={Path(x['path']).name:x for x in original};records=[];headers=[];total_read=0;session=requests.Session()
 def read_range(start,size):
  nonlocal total_read
  end=start+size-1
  if total_read+size>CAP:raise RuntimeError('Selective download exceeds explicit1GiB cap; stop, do not fall back to full bundle')
  for attempt in range(3):
   try:
    with session.get(URL,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'},stream=True,timeout=(15,60)) as r:
     expected=f'bytes {start}-{end}/{TOTAL}'
     if r.status_code!=206 or r.headers.get('Content-Range')!=expected:raise ValueError('Range not honored; whole-bundle fallback prohibited')
     content=bytearray()
     for chunk in r.iter_content(2**20):
      content.extend(chunk)
      if len(content)>size:raise ValueError('Server exceeded requested range')
     if len(content)!=size:raise ValueError('Incomplete byte range')
    total_read+=size;return bytes(content)
   except requests.RequestException:
    if attempt==2:raise
  raise RuntimeError('Range retry exhausted')
 offset=0
 while offset+512<=TOTAL:
  block=read_range(offset,512)
  if block==b'\0'*512:break
  stated=int(block[148:156].strip(b'\0 ') or b'0',8);actual=sum(block[:148])+8*32+sum(block[156:])
  if stated!=actual:raise ValueError('TAR header checksum failed')
  name=block[:100].split(b'\0')[0].decode();prefix=block[345:500].split(b'\0')[0].decode();name=(prefix+'/' if prefix else '')+name
  size=int(block[124:136].strip(b'\0 ') or b'0',8);kind=block[156:157];base=Path(name).name;headers.append({'offset':offset,'name':name,'size':size,'type':kind.decode(errors='replace')})
  if base in wanted and kind in [b'0',b'\0']:
   expected=wanted[base];dest=out/'raw/bassiouni'/base
   if size!=expected['bytes']:raise ValueError('Member size differs from original manifest: '+base)
   if dest.exists():
    if sha(dest)!=expected['sha256']:raise ValueError('Existing selected file checksum mismatch')
   else:
    existing=root/'data/journal_revision/transformer_egnn/acquisition_recheck/raw/bassiouni'/base
    if existing.exists() and sha(existing)==expected['sha256']:shutil.copyfile(existing,dest)
    else:
     payload=read_range(offset+512,size)
     if hashlib.sha256(payload).hexdigest()!=expected['sha256']:raise ValueError('Selected member SHA256 mismatch: '+base)
     part=dest.with_suffix(dest.suffix+'.part');part.write_bytes(payload);part.replace(dest)
   records.append({**expected,'selected_tar_header_offset':offset,'cache_relative':str(dest.relative_to(root)),'verified_sha256':sha(dest)})
   print('GEO SELECTED VERIFIED',base,size,flush=True)
  offset+=512+((size+511)//512)*512
 missing=sorted(set(wanted)-{Path(x['path']).name for x in records})
 result={'status':'COMPLETE' if not missing else 'BLOCKED_INPUTS','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_url':URL,'advertised_whole_bundle_bytes':TOTAL,
  'actual_range_bytes_downloaded':total_read,'selective_cap_bytes':CAP,'whole_bundle_downloaded':False,'selected_files':records,'missing':missing,'headers_scanned':len(headers),'license_source':'https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html'}
 target=report if not report.exists() else out/'REACQUISITION.json';target.write_text(json.dumps(result,indent=2)+'\n')
 (out/'TAR_HEADER_INDEX.json').write_text(json.dumps(headers,indent=2)+'\n');return result
