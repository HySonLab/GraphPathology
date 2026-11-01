"""Seal C10 or verify it read-only, including every frozen parent checkpoint."""
from contract10 import *
import argparse

def main():
 p=OUT/'C10_MANIFEST.json';dp=p.with_suffix('.sha256');args=argparse.ArgumentParser();args.add_argument('--verify',action='store_true');verify=args.parse_args().verify;check_lock()
 if verify:
  subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C9/freeze.py'),'--verify'],check=True)
  assert sha(p)==dp.read_text().split()[0];obj=json.loads(p.read_text())
  for name,h in obj['artifact_sha256'].items():assert sha(ROOT/name)==h,name
  reg=obj['registry_prefix']
  with (ROOT/reg['path']).open('rb') as f:assert hashlib.sha256(f.read(reg['bytes'])).hexdigest()==reg['sha256']
  print('PASS:',len(obj['artifact_sha256']),'C10 artifacts, immutable scientific sources, registry prefix, and all C0–C9 frozen parents');return
 assert not p.exists(),'C10 already frozen';status=json.loads((OUT/'C10_STATUS.json').read_text());assert status['technical_audit']['status']=='PASS' and status['C0_C9_unchanged'];assert json.loads((OUT/'MANUSCRIPT/PACKAGE_AUDIT.json').read_text())['status']=='PASS';assert 'PASS: 1829 C9 artifacts' in (B/'logs/PARENT_VERIFICATION_AFTER.log').read_text()
 skip={'__pycache__','.pytest_cache','.git'};files={str(f.relative_to(ROOT)):sha(f) for f in OUT.rglob('*') if f.is_file() and not skip.intersection(f.relative_to(OUT).parts) and f.name not in ['C10_MANIFEST.json','C10_MANIFEST.sha256']};reg=ROOT/'results/experiment_registry.csv'
 obj={'timestamp':now(),'source_commit':commit(),'parent_C9_manifest_sha256':sha(ROOT/'C9/C9_MANIFEST.json'),'spec_sha256':sha(A/'SPEC_LOCK.json'),'status':status,'artifact_sha256':files,'registry_prefix':{'path':str(reg.relative_to(ROOT)),'bytes':reg.stat().st_size,'sha256':sha(reg)},'local_assets':'Raw downloads, extracted counts/images, caches, predictions and uncompressed gene metrics are retained locally under C10 and ignored by Git. All are hashed. Source,compressed gene metrics,patient metrics,masks,logs,figures,tables,manuscript and ledger are committed.'};dump(p,obj);dp.write_text(sha(p)+'  C10_MANIFEST.json\n');print('FROZEN',len(files),'artifacts',status['overall'])
if __name__=='__main__':main()
