"""Freeze C2 outputs and verify them without changing C0 or C1."""
import argparse,hashlib,json,platform,subprocess,sys
from importlib.metadata import version
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C2'

def sha(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
 path=OUT/'C2_SPEC_LOCK.json';digest_path=OUT/'C2_SPEC_LOCK.sha256'
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'C1/freeze.py'),'--verify'],check=True)
 if args.verify:
  assert sha(path)==digest_path.read_text().split()[0]
  lock=json.loads(path.read_text());assert sha(ROOT/'C1/C1_SPEC_LOCK.json')==lock['parent_C1_lock_sha256']
  for rel,digest in lock['artifact_sha256'].items():assert sha(ROOT/rel)==digest,rel
  prefix=lock['registry_prefix']
  with (ROOT/prefix['path']).open('rb') as f:assert hashlib.sha256(f.read(prefix['bytes'])).hexdigest()==prefix['sha256']
  print('PASS: C2 lock,',len(lock['artifact_sha256']),'artifacts, registry prefix and frozen parent inputs match');return
 if path.exists():raise SystemExit('C2 is already frozen')
 audit=json.loads((OUT/'logs/INTEGRITY_AUDIT.json').read_text());assert audit['status']=='PASS'
 decision=json.loads((OUT/'logs/C2_DECISION.json').read_text());selection=json.loads((OUT/'SELECTION_LOCK.json').read_text())
 artifacts={str(p.relative_to(ROOT)):sha(p) for p in OUT.rglob('*') if p.is_file() and '__pycache__' not in str(p)
            and '.pytest_cache' not in str(p) and p.name not in ['C2_SPEC_LOCK.json','C2_SPEC_LOCK.sha256']}
 registry=ROOT/'results/experiment_registry.csv'
 runs=[json.loads(p.read_text()) for p in (OUT/'logs/runs').glob('*.json')]
 lock={'checkpoint':'C2','status':decision['status'],'timestamp':datetime.now(timezone.utc).isoformat(),
       'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
       'fit_source_commits':sorted(set(r['source_commit'] for r in runs)),
       'parent_C1_lock_sha256':sha(ROOT/'C1/C1_SPEC_LOCK.json'),'parent_C0_lock_sha256':sha(ROOT/'reports/C0/C0_SPEC_LOCK.json'),
       'selection_lock_sha256':sha(OUT/'SELECTION_LOCK.json'),'analysis_lock_sha256':sha(OUT/'ANALYSIS_LOCK.json'),
       'selected':selection['selected'],'gene_panels':json.loads((ROOT/'C1/SELECTION_LOCK.json').read_text())['panels'],
       'primary_splits':json.loads((ROOT/'reports/C0/C0_SPLITS.json').read_text()),
       'protocol':json.loads((OUT/'configs/protocol.json').read_text()),'gene_analysis_rules':json.loads((OUT/'configs/analysis.json').read_text()),
       'frozen_scalers':{r['preprocessing']:{'feature_sha256':r['feature_scaler_sha256'],'target_sha256':r['target_scaler_sha256']} for r in runs},
       'software':{'python':sys.version,**{p:version(p) for p in ['numpy','scipy','pandas','torch','scikit-learn','matplotlib','pytest']}},
       'hardware':{'platform':platform.platform(),'CPU':'Intel Core i7-14700','threads':4,'device':'CPU'},
       'integrity_audit':audit,'decision':decision,'artifact_sha256':artifacts,
       'registry_prefix':{'path':'results/experiment_registry.csv','bytes':registry.stat().st_size,'sha256':sha(registry)},
       'local_binary_note':'C2/models and C2/predictions are local and hashed but excluded from Git. Frozen C1 is never retrained or replayed by C2.'}
 path.write_text(json.dumps(lock,indent=2,allow_nan=False)+'\n');digest_path.write_text(sha(path)+'  C2_SPEC_LOCK.json\n')
 print('Frozen',decision['status'],'with',len(artifacts),'C2 file hashes; lock',sha(path))
if __name__=='__main__':main()
