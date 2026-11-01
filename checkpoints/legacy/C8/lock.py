"""Seal all C8 scientific choices before fitting/evaluating C8 outcomes."""
from contract import *
def main():
 p=A/'C8A_SPEC_LOCK.json';assert not p.exists();files=['C8/contract.py','C8/run.py','C8/analyze.py','C8/configs/protocol.json'];reg=ROOT/'results/experiment_registry.csv';assets={}
 for sample in SAMPLES:
  for suffix in ['_targets.npz','_S2-50_image.npz','_metadata.npz']:
   q=c7.CACHE/(sample+suffix);assets[str(q.relative_to(ROOT))]=sha(q)
 for sample in SOURCE:assets[f'data/cache/{sample}_targets.npz']=sha(ROOT/'data/cache'/f'{sample}_targets.npz')
 masks_index={r['mask_path']:sha(ROOT/r['mask_path']) for sample in SAMPLES for r in masks(sample)}
 lock={'timestamp':now(),'source_commit':commit(),'protocol':cfg(),'genes':GENES,'source_sha256':{p:sha(ROOT/p) for p in files},'input_sha256':assets,'mask_sha256':masks_index,'parent_C7_manifest_sha256':sha(ROOT/'C7/C7_MANIFEST.json'),'registry_prefix':{'path':str(reg.relative_to(ROOT)),'bytes':reg.stat().st_size,'sha256':sha(reg)},'HER2ST_is_now_development':True,'previous_C7_outcomes_known':True,'C8_outcomes_used_for_spec':False};dump(p,lock);p.with_suffix('.sha256').write_text(sha(p)+'  C8A_SPEC_LOCK.json\n');print('C8 scientific lock',sha(p))
if __name__=='__main__':main()
