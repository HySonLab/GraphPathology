"""Freeze C4B only after C4A scientific and technical gates pass."""
from common import *
B=OUT/'C4B'
if __name__=='__main__':
 assert json.loads((A/'logs/GATE_RESULT.json').read_text())['status']=='PASS_PATHOLOGY_REPRESENTATION'
 assert json.loads((A/'logs/INTEGRITY_AUDIT.json').read_text())['status']=='PASS'
 path=B/'C4B_SPEC_LOCK.json';assert not path.exists()
 files=[OUT/'b_spatial.py',OUT/'b_analyze.py',OUT/'common.py',OUT/'run.py',B/'configs/protocol.json']
 artifacts={str(p.relative_to(ROOT)):sha(p) for p in A.rglob('*') if p.is_file() and not {'__pycache__','.pytest_cache'}.intersection(p.parts)}
 dump(path,{'timestamp':now(),'source_commit':commit(),'C4A_gate_sha256':sha(A/'logs/GATE_RESULT.json'),'C4A_integrity_sha256':sha(A/'logs/INTEGRITY_AUDIT.json'),'parent_C3_manifest_sha256':sha(ROOT/'C3/C3_MANIFEST.json'),'C4A_artifacts_sha256':artifacts,'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'protocol':json.loads((B/'configs/protocol.json').read_text()),'prior_C4A_test_results_known':True,'C4B_test_predictions_started':False})
 (B/'C4B_SPEC_LOCK.sha256').write_text(sha(path)+'  C4B_SPEC_LOCK.json\n');print('C4B locked',sha(path))
