from common import *
if __name__=='__main__':
 p=OUT/'C5_SPEC_LOCK.json';assert not p.exists();files=[OUT/f for f in ['common.py','engine.py','run.py','analyze.py','configs/protocol.json','configs/masks.json','C5_MASK_INDEX.csv']];reg=ROOT/'results/experiment_registry.csv'
 dump(p,{'timestamp':now(),'source_commit':commit(),'parent_C4_manifest_sha256':sha(ROOT/'C4/C4_MANIFEST.json'),'source_sha256':{str(f.relative_to(ROOT)):sha(f) for f in files},'protocol':cfg(),'primary_hidden_labels_used_for_C5_selection':False,'prior_C0_C4_test_results_known':True,'registry_prefix':{'bytes':reg.stat().st_size,'sha256':sha(reg)}});(OUT/'C5_SPEC_LOCK.sha256').write_text(sha(p)+'  C5_SPEC_LOCK.json\n');print('C5 locked',sha(p))
