"""Freeze or read-only verify the complete C1 reproducibility contract."""
import argparse,hashlib,json,platform,subprocess,sys
from importlib.metadata import version
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C1'

def sha(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
 path=OUT/'C1_SPEC_LOCK.json';digest_path=OUT/'C1_SPEC_LOCK.sha256'
 if args.verify:
  assert sha(path)==digest_path.read_text().split()[0],'Lock hash changed'
  lock=json.loads(path.read_text())
  for rel,digest in {**lock['artifact_sha256'],**lock['c0_input_hashes']}.items():
   assert sha(ROOT/rel)==digest,rel
  prefix=lock['registry_prefix']
  with (ROOT/prefix['path']).open('rb') as stream:
   assert hashlib.sha256(stream.read(prefix['bytes'])).hexdigest()==prefix['sha256'],'C1 registry prefix changed'
  print(f"PASS: C1 lock and {len(lock['artifact_sha256'])} C1 artifacts plus {len(lock['c0_input_hashes'])} C0 inputs match")
  return
 if path.exists():raise SystemExit('C1 already frozen; refusing to overwrite')
 selection=json.loads((OUT/'SELECTION_LOCK.json').read_text())
 audit=json.loads((OUT/'logs/INTEGRITY_AUDIT.json').read_text())
 decision=json.loads((OUT/'logs/C1_DECISION.json').read_text())
 assert audit['status']=='PASS' and audit['evaluated_runs']==112
 c0_lock=json.loads((ROOT/'reports/C0/C0_SPEC_LOCK.json').read_text())
 c0_manifest=json.loads((ROOT/'reports/C0/C0_DATA_MANIFEST.json').read_text())
 c0_hashes={**c0_lock['source_sha256'],**c0_lock['artifact_sha256'],
            'reports/C0/C0_SPEC_LOCK.json':selection['c0_lock_sha256']}
 for slide in c0_manifest['slides']:
  c0_hashes[slide['feature_cache']]=slide['feature_sha256'];c0_hashes[slide['target_cache']]=slide['target_sha256']
 for record in c0_manifest['downloads']['files']:c0_hashes[record['path']]=record['sha256']
 for rel,digest in c0_hashes.items():assert sha(ROOT/rel)==digest,rel
 scalers={}
 for folder in [OUT/'preprocessing',OUT/'lopo']:
  for p in folder.rglob('*_scaler.npz'):
   with np.load(p) as f:state={k:f[k].tolist() for k in f.files}
   scalers[str(p.relative_to(ROOT))]={'sha256':sha(p),'parameters':state}
 artifacts={str(p.relative_to(ROOT)):sha(p) for p in OUT.rglob('*') if p.is_file()
            and '__pycache__' not in str(p) and '.pytest_cache' not in str(p)
            and p.name not in ['C1_SPEC_LOCK.json','C1_SPEC_LOCK.sha256']}
 registry_path=ROOT/'results/experiment_registry.csv'
 registry_prefix={'path':'results/experiment_registry.csv','bytes':registry_path.stat().st_size,'sha256':sha(registry_path)}
 source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
 lock={'checkpoint':'C1','status':decision['status'],'timestamp':datetime.now(timezone.utc).isoformat(),
       'source_commit':source_commit,'fitting_source_commit':selection['source_commit'],
       'c0_source_commit':c0_lock['source_commit'],'c0_input_hashes':c0_hashes,
       'primary_splits':json.loads((ROOT/'reports/C0/C0_SPLITS.json').read_text()),
       'selection_lock_sha256':sha(OUT/'SELECTION_LOCK.json'),'gene_panels':selection['panels'],
       'single_gene':selection['panels']['S1'][0],'selected_configurations':selection['selected'],
       'protocol':json.loads((OUT/'configs/protocol.json').read_text()),'preprocessing':scalers,
       'software':{'python':sys.version,**{p:version(p) for p in ['numpy','scipy','pandas','scikit-learn','torch','torchvision','matplotlib','h5py','pytest']}},
       'cpu_environment':{'platform':platform.platform(),'model':'Intel Core i7-14700','threads':4,'device':'cpu'},
       'integrity_audit':audit,'decision':decision,'artifact_sha256':artifacts,'registry_prefix':registry_prefix,
       'local_binary_note':'C1/models and C1/predictions are locally retained and hashed but excluded from Git; regenerate from frozen inputs/code if absent.'}
 path.write_text(json.dumps(lock,indent=2,allow_nan=False)+'\n')
 digest_path.write_text(sha(path)+'  C1_SPEC_LOCK.json\n')
 print(f"Frozen {decision['status']}: {len(artifacts)} artifact hashes; lock {sha(path)}")
if __name__=='__main__':main()
