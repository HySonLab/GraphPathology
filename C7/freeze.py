"""Freeze or read-only verify final C7 artifacts and the entire parent chain."""
import argparse,json,hashlib,subprocess
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C7'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();path=OUT/'C7_MANIFEST.json';digest=OUT/'C7_MANIFEST.sha256'
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C6/freeze.py'),'--verify'],check=True)
 spec=OUT/'C7A/C7A_SPEC_LOCK.json';assert sha(spec)==spec.with_suffix('.sha256').read_text().split()[0];s=json.loads(spec.read_text())
 for p,h in s['source_sha256'].items():assert sha(ROOT/p)==h,p
 if args.verify:
  assert sha(path)==digest.read_text().split()[0];lock=json.loads(path.read_text())
  for p,h in lock['artifact_sha256'].items():assert sha(ROOT/p)==h,p
  prefix=lock['registry_prefix']
  with (ROOT/prefix['path']).open('rb') as f:assert hashlib.sha256(f.read(prefix['bytes'])).hexdigest()==prefix['sha256']
  print('PASS:',len(lock['artifact_sha256']),'C7 artifacts, registry prefix, locked scientific sources, and frozen C0–C6');return
 assert not path.exists(),'C7 already frozen'
 status=json.loads((OUT/'C7_STATUS.json').read_text());assert status['technical_audit']['status']=='PASS' and status['distance_audit']['status']=='PASS';assert status['C0_C6_unchanged'] and status['no_post_outcome_tuning'];assert status['main_evaluations']==8448 and status['negative_control_evaluations']==108;assert status['C7_status']=='MIXED_EXTERNAL_VALIDATION';assert not (OUT/'C7C').exists()
 skip={'__pycache__','.pytest_cache','deps','.git'}
 files={str(p.relative_to(ROOT)):sha(p) for p in OUT.rglob('*') if p.is_file() and not skip.intersection(p.relative_to(OUT).parts) and p.name not in ['C7_MANIFEST.json','C7_MANIFEST.sha256']}
 reg=ROOT/'results/experiment_registry.csv'
 lock={'timestamp':datetime.now(timezone.utc).isoformat(),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'parent_C6_manifest_sha256':sha(ROOT/'C6/C6_MANIFEST.json'),'scientific_spec_sha256':sha(spec),'status':status,'artifact_sha256':files,'registry_prefix':{'path':str(reg.relative_to(ROOT)),'bytes':reg.stat().st_size,'sha256':sha(reg)},'local_binary_note':'Raw binary data, external caches, predictions, and the uncompressed gene CSV remain local and are hashed. Compressed metrics, masks, metadata, code, figures and reports are committed.'}
 path.write_text(json.dumps(lock,indent=2,allow_nan=False)+'\n');digest.write_text(sha(path)+'  C7_MANIFEST.json\n');print('Frozen',status['C7_status'],len(files),'artifacts')
if __name__=='__main__':main()
