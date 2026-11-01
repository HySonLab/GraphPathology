from common import *
if __name__=='__main__':
 path=A/'C4A_SPEC_LOCK.json';assert not path.exists()
 files=[OUT/f for f in ['common.py','run.py','extract.py','analyze.py']]+[A/'configs/protocol.json',A/'configs/extraction.json',A/'configs/encoder_files.json']
 dump(path,{'timestamp':now(),'source_commit':commit(),'parent_C3_manifest_sha256':sha(ROOT/'C3/C3_MANIFEST.json'),'protocol':json.loads((A/'configs/protocol.json').read_text()),'encoder':json.loads((A/'configs/extraction.json').read_text()),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'test_evaluation_started':False,'prior_checkpoint_test_results_known':True})
 (A/'C4A_SPEC_LOCK.sha256').write_text(sha(path)+'  C4A_SPEC_LOCK.json\n');print('Locked',sha(path))
