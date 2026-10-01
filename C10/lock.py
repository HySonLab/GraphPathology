"""Freeze second-cohort patient scope, scientific code and H&E inputs before outcomes."""
from contract10 import *
def main():
 assert not (A/'SPEC_LOCK.json').exists();assert len(SAMPLES)==22 and len(genes())==50
 files=[OUT/n for n in ['contract10.py','protocol.py','run.py','analyze.py','qualify.py','prepare_images.py','lock.py','configs/protocol.json','tests/test_contract.py']]+[ROOT/'C9/contract.py',ROOT/'C9/C9A/T1_SCALER.json',ROOT/'C7/policy.py',ROOT/'C7/support.py',ROOT/'reports/C0/C0_NORMALIZATION.json',ROOT/'data/weights/hub/checkpoints/resnet18-f37072fd.pth']
 assets={};models={};raw={}
 for s in SAMPLES:
  for suffix in ['_metadata.npz','_features.npz','_image.npz']:p=CACHE/(s+suffix);assets[str(p.relative_to(ROOT))]=sha(p)
  r=json.loads((A/'logs'/(s+'_EXTRACTION.json')).read_text())
  for m in r['parents']:
   for key in ['checkpoint','feature_scaler','target_scaler']:models[m[key]]=sha(ROOT/m[key])
  for name in cfg()['inputs'][s].values():raw[name]=sha(ROOT/name)
 images=[raw[cfg()['inputs'][s]['image']] for s in SAMPLES];assert len(set(images))==22
 previous=set()
 for checkpoint in ['C7','C8','C9']:previous.update(json.loads((ROOT/checkpoint/(checkpoint+'_MANIFEST.json')).read_text())['artifact_sha256'].values())
 assert not set(images)&previous
 prior=subprocess.run(['git','grep','-l','-F','-e','GSE210616','-e','GSM6433585','--',':!C10'],cwd=ROOT,text=True,capture_output=True);assert prior.returncode==1 and not prior.stdout
 dump(A/'logs/FRESHNESS_AUDIT.json',{'status':'PASS','distinct_source_patient_ids':22,'distinct_first_section_images':22,'parent_image_sha256_matches':0,'previous_project_identifier_matches':0,'identity_scope':'Source-study patient IDs plus available image hashes, not patient reidentification across unrelated anonymized studies','C10_target_calibration':'none','C9_scaler_frozen_and_external_to_C10':True})
 lock={'timestamp':now(),'source_commit':commit(),'parent_C9_manifest_sha256':sha(ROOT/'C9/C9_MANIFEST.json'),'protocol':cfg(),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'image_assets_sha256':assets,'model_parents_sha256':models,'raw_asset_sha256':raw,'evaluation_labels_not_accessed':True,'registry_prefix_before_C10':{'bytes':(ROOT/'results/experiment_registry.csv').stat().st_size,'sha256':sha(ROOT/'results/experiment_registry.csv')}}
 dump(A/'SPEC_LOCK.json',lock);(A/'SPEC_LOCK.sha256').write_text(sha(A/'SPEC_LOCK.json')+'  SPEC_LOCK.json\n');print('LOCKED',sha(A/'SPEC_LOCK.json'))
if __name__=='__main__':main()
