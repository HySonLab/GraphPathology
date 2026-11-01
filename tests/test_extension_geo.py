import hashlib,json,tarfile
from pathlib import Path
import pytest
from graphpathology import extension_geo


def setup(tmp_path,monkeypatch,ignore_range=False):
 payload=b'exact original input';info=tarfile.TarInfo('GSM_test_counts.h5');info.size=len(payload);archive=info.tobuf()+payload+b'\0'*(512-len(payload))+b'\0'*1024;calls=[]
 folder=tmp_path/'checkpoints/legacy/C10/C10A/logs';folder.mkdir(parents=True);(folder/'DOWNLOADS.json').write_text(json.dumps([{'path':'C10/raw/GSM_test_counts.h5','bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest()}]));(tmp_path/'experiments/journal_revision/transformer_egnn').mkdir(parents=True)
 class Response:
  def __init__(self,start,end):self.start=start;self.end=end;self.status_code=200 if ignore_range else 206;self.headers={} if ignore_range else {'Content-Range':f'bytes {start}-{end}/{len(archive)}'}
  def __enter__(self):return self
  def __exit__(self,*args):pass
  def iter_content(self,*args):
   assert not ignore_range,'Must refuse range-ignoring server before reading any body'
   yield archive[self.start:self.end+1]
 class Session:
  def get(self,url,headers,**kwargs):
   start,end=map(int,headers['Range'][6:].split('-'));calls.append((start,end));return Response(start,end)
 monkeypatch.setattr(extension_geo,'TOTAL',len(archive));monkeypatch.setattr(extension_geo.requests,'Session',Session)
 return payload,calls


def test_only_required_member_ranges_and_exact_checksum(tmp_path,monkeypatch):
 payload,calls=setup(tmp_path,monkeypatch);r=extension_geo.acquire(tmp_path)
 assert r['status']=='COMPLETE' and not r['whole_bundle_downloaded']
 assert r['actual_range_bytes_downloaded']==1024+len(payload)
 assert calls==[(0,511),(512,512+len(payload)-1),(1024,1535)]
 count=len(calls);assert extension_geo.acquire(tmp_path)==r and len(calls)==count


def test_ignored_range_never_reads_large_body(tmp_path,monkeypatch):
 setup(tmp_path,monkeypatch,ignore_range=True)
 with pytest.raises(ValueError,match='whole-bundle fallback prohibited'):extension_geo.acquire(tmp_path)
