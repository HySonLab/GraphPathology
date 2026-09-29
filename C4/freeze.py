"""Hash final C4 artifacts and verify the complete frozen checkpoint chain."""
import argparse,json,hashlib,subprocess
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C4'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();path=OUT/'C4_MANIFEST.json';digest=OUT/'C4_MANIFEST.sha256'
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C3/freeze.py'),'--verify'],check=True)
 if args.verify:
  assert sha(path)==digest.read_text().split()[0];lock=json.loads(path.read_text())
  for p,h in lock['artifact_sha256'].items():assert sha(ROOT/p)==h,p
  prefix=lock['registry_prefix']
  with (ROOT/prefix['path']).open('rb') as f:assert hashlib.sha256(f.read(prefix['bytes'])).hexdigest()==prefix['sha256']
  print('PASS:',len(lock['artifact_sha256']),'C4 artifacts, registry prefix, and all frozen parents');return
 assert not path.exists(),'C4 already frozen'
 status=json.loads((OUT/'C4_STATUS.json').read_text());assert status['technical_audit']['status']=='PASS'
 skip={'__pycache__','.pytest_cache','deps','.git','CMRCNet'}
 files={str(p.relative_to(ROOT)):sha(p) for p in OUT.rglob('*') if p.is_file() and not skip.intersection(p.relative_to(OUT).parts) and p.name not in ['C4_MANIFEST.json','C4_MANIFEST.sha256']}
 reg=ROOT/'results/experiment_registry.csv'
 lock={'timestamp':datetime.now(timezone.utc).isoformat(),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'parent_C3_manifest_sha256':sha(ROOT/'C3/C3_MANIFEST.json'),'status':status,'artifact_sha256':files,'registry_prefix':{'path':str(reg.relative_to(ROOT)),'bytes':reg.stat().st_size,'sha256':sha(reg)},'local_binary_note':'Weights, embedding caches, model checkpoints and predictions remain local and are hashed here; isolated dependencies are versioned in execution_environment.json and excluded from the artifact hashes.'}
 path.write_text(json.dumps(lock,indent=2,allow_nan=False)+'\n');digest.write_text(sha(path)+'  C4_MANIFEST.json\n');print('Frozen',status['overall_C4_status'],len(files),'artifacts')
if __name__=='__main__':main()
