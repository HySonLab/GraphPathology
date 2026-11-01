"""Hash final C6 artifacts and verify the complete frozen checkpoint chain."""
import argparse,json,hashlib,subprocess
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C6'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();path=OUT/'C6_MANIFEST.json';digest=OUT/'C6_MANIFEST.sha256'
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C5/freeze.py'),'--verify'],check=True)
 if args.verify:
  assert sha(path)==digest.read_text().split()[0];lock=json.loads(path.read_text())
  for p,h in lock['artifact_sha256'].items():assert sha(ROOT/p)==h,p
  prefix=lock['registry_prefix']
  with (ROOT/prefix['path']).open('rb') as f:assert hashlib.sha256(f.read(prefix['bytes'])).hexdigest()==prefix['sha256']
  print('PASS:',len(lock['artifact_sha256']),'C6 artifacts, registry prefix, and all frozen parents');return
 assert not path.exists(),'C6 already frozen'
 status=json.loads((OUT/'C6_STATUS.json').read_text());assert status['technical_audit']['status']=='PASS'
 assert status['C0_C5_unchanged'] and not status['external_validation_executed']
 if status['C6A_status']=='PASS_SAMPLING_EFFECT':assert status['C6B_status'] in ['PASS_ADAPTIVE_SAMPLING','MIXED_ADAPTIVE_SAMPLING','NO_ADAPTIVE_ADVANTAGE']
 else:assert not (OUT/'C6B').exists()
 skip={'__pycache__','.pytest_cache','deps','.git','CMRCNet'}
 files={str(p.relative_to(ROOT)):sha(p) for p in OUT.rglob('*') if p.is_file() and not skip.intersection(p.relative_to(OUT).parts) and p.name not in ['C6_MANIFEST.json','C6_MANIFEST.sha256']}
 reg=ROOT/'results/experiment_registry.csv'
 lock={'timestamp':datetime.now(timezone.utc).isoformat(),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'parent_C5_manifest_sha256':sha(ROOT/'C5/C5_MANIFEST.json'),'status':status,'artifact_sha256':files,'registry_prefix':{'path':str(reg.relative_to(ROOT)),'bytes':reg.stat().st_size,'sha256':sha(reg)},'local_binary_note':'Large prediction arrays and uncompressed gene/selection tables remain local and are hashed. Gzip tables, sampling traces, code, figures, reports and provenance are committed. Frozen reconstruction artifacts remain under C5.'}
 path.write_text(json.dumps(lock,indent=2,allow_nan=False)+'\n');digest.write_text(sha(path)+'  C6_MANIFEST.json\n');print('Frozen',status['C6_status'],len(files),'artifacts')
if __name__=='__main__':main()
