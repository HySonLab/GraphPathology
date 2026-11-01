"""Freeze/read-only verify C8 and all earlier checkpoint artifacts."""
import argparse,json,hashlib,subprocess
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C8'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();p=OUT/'C8_MANIFEST.json';digest=OUT/'C8_MANIFEST.sha256'
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C7/freeze.py'),'--verify'],check=True)
 spec=OUT/'C8A/C8A_SPEC_LOCK.json';assert sha(spec)==spec.with_suffix('.sha256').read_text().split()[0];s=json.loads(spec.read_text())
 for name,h in s['source_sha256'].items():assert sha(ROOT/name)==h,name
 if args.verify:
  assert sha(p)==digest.read_text().split()[0];lock=json.loads(p.read_text())
  for name,h in lock['artifact_sha256'].items():assert sha(ROOT/name)==h,name
  pr=lock['registry_prefix']
  with (ROOT/pr['path']).open('rb') as f:assert hashlib.sha256(f.read(pr['bytes'])).hexdigest()==pr['sha256']
  print('PASS:',len(lock['artifact_sha256']),'C8 artifacts, immutable scientific sources, registry prefix, and complete frozen C0–C7');return
 assert not p.exists(),'C8 already frozen';status=json.loads((OUT/'C8_STATUS.json').read_text());assert status['technical_audit']['status']=='PASS';assert status['C0_C7_unchanged'] and status['C7_zero_shot_preserved_exactly'];assert status['C8B_status']=='MIXED_FUSION_CALIBRATION' and not (OUT/'C8C').exists()
 skip={'__pycache__','.pytest_cache','.git'};files={str(f.relative_to(ROOT)):sha(f) for f in OUT.rglob('*') if f.is_file() and not skip.intersection(f.relative_to(OUT).parts) and f.name not in ['C8_MANIFEST.json','C8_MANIFEST.sha256']};reg=ROOT/'results/experiment_registry.csv'
 lock={'timestamp':datetime.now(timezone.utc).isoformat(),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'parent_C7_manifest_sha256':sha(ROOT/'C7/C7_MANIFEST.json'),'scientific_spec_sha256':sha(spec),'status':status,'artifact_sha256':files,'registry_prefix':{'path':str(reg.relative_to(ROOT)),'bytes':reg.stat().st_size,'sha256':sha(reg)},'local_artifacts':'Prediction/inner-cache arrays and uncompressed gene CSV are local, ignored by Git, and included in these hashes. Compressed metrics, source, objects, logs, reports and figures are committed.'};p.write_text(json.dumps(lock,indent=2,allow_nan=False)+'\n');digest.write_text(sha(p)+'  C8_MANIFEST.json\n');print('Frozen',status['C8_status'],len(files),'artifacts')
if __name__=='__main__':main()
