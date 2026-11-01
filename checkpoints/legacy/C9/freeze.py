"""Seal C9; read-only verification includes the full frozen parent chain."""
import json,hashlib,argparse,subprocess
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C9'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();p=OUT/'C9_MANIFEST.json';dp=OUT/'C9_MANIFEST.sha256';sp=OUT/'C9A/C9A_SPEC_LOCK.json';assert sha(sp)==sp.with_suffix('.sha256').read_text().split()[0];spec=json.loads(sp.read_text())
 for name,h in spec['source_sha256'].items():assert sha(ROOT/name)==h,name
 if args.verify:
  subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C8/freeze.py'),'--verify'],check=True)
  assert sha(p)==dp.read_text().split()[0];obj=json.loads(p.read_text())
  for name,h in obj['artifact_sha256'].items():assert sha(ROOT/name)==h,name
  reg=obj['registry_prefix']
  with (ROOT/reg['path']).open('rb') as f:assert hashlib.sha256(f.read(reg['bytes'])).hexdigest()==reg['sha256']
  print('PASS:',len(obj['artifact_sha256']),'C9 artifacts, scientific lock, registry prefix, and unchanged C0-C8');return
 assert not p.exists(),'C9 already frozen';status=json.loads((OUT/'C9_STATUS.json').read_text());assert status['technical_audit']['status']=='PASS' and status['C0_C8_unchanged'];assert 'PASS: 27448 C8 artifacts' in (OUT/'C9B/logs/PARENT_VERIFICATION.txt').read_text()
 skip={'__pycache__','.pytest_cache','.git'};files={str(f.relative_to(ROOT)):sha(f) for f in OUT.rglob('*') if f.is_file() and not skip.intersection(f.relative_to(OUT).parts) and f.name not in ['C9_MANIFEST.json','C9_MANIFEST.sha256']};reg=ROOT/'results/experiment_registry.csv'
 obj={'timestamp':datetime.now(timezone.utc).isoformat(),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'parent_C8_manifest_sha256':sha(ROOT/'C8/C8_MANIFEST.json'),'spec_sha256':sha(sp),'status':status,'artifact_sha256':files,'registry_prefix':{'path':str(reg.relative_to(ROOT)),'bytes':reg.stat().st_size,'sha256':sha(reg)},'local_assets':'Raw archives/extracted matrices, image feature caches and predictions are retained locally under C9 and ignored by Git; all are covered by these SHA-256 hashes. Source/metrics/configs/scalers/masks/logs/reports/figures are committed.'};p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');dp.write_text(sha(p)+'  C9_MANIFEST.json\n');print('FROZEN',len(files),'artifacts',status['overall'])
if __name__=='__main__':main()
