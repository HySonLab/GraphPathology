"""Create the immutable pre-test specification, separate from the final artifact manifest."""
from common import *
if __name__=='__main__':
 path=OUT/'C3A_SPEC_LOCK.json'
 assert not path.exists(),'Specification already locked'
 files=list(OUT.glob('*.py'))+[OUT/'configs/protocol.json',ROOT/'C3/C3_EXECUTION_PLAN_UPDATED.md',ROOT/'C1/core.py',ROOT/'C2/run.py',ROOT/'C2/graph.py']
 dump(path,{'timestamp':now(),'source_commit':commit(),'parent_C2_lock_sha256':sha(ROOT/'C2/C2_SPEC_LOCK.json'),'C3_new_test_evaluation_started':False,'prior_C1_C2_test_results_known':True,'protocol':CFG,'selected_C2':SELECTED,'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files}})
 (OUT/'C3A_SPEC_LOCK.sha256').write_text(sha(path)+'  C3A_SPEC_LOCK.json\n')
 print('C3A specification frozen:',sha(path))
