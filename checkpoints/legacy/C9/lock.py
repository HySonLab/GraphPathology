"""Seal C9 specifications and image inputs before calibration/evaluation access."""
from contract import *
def main():
 assert not (A/'C9A_SPEC_LOCK.json').exists();protocol=cfg();assert len(genes())==50 and set(CAL).isdisjoint(EVAL)
 files=[OUT/n for n in ['contract.py','run.py','analyze.py','qualify.py','prepare_images.py','prepare_protocol.py','lock.py','configs/protocol.json','tests/test_contract.py']]
 files += [ROOT/'C7/policy.py',ROOT/'C7/support.py',ROOT/'reports/C0/C0_NORMALIZATION.json',ROOT/'data/weights/hub/checkpoints/resnet18-f37072fd.pth']
 assets={};modelparents={}
 for s in SAMPLES:
  for suffix in ['_metadata.npz','_features.npz','_image.npz']:
   p=CACHE/(s+suffix);assets[str(p.relative_to(ROOT))]=sha(p)
  rec=json.loads((A/'logs'/(s+'_EXTRACTION.json')).read_text())
  for r in rec['parents']:
   for key in ['checkpoint','feature_scaler','target_scaler']:p=ROOT/r[key];modelparents[str(p.relative_to(ROOT))]=sha(p)
 parentmanifest=json.loads((ROOT/'C8/C8_MANIFEST.json').read_text());oldhash=set(parentmanifest['artifact_sha256'].values());sourcefiles=list((OUT/'raw').glob('*.tar.gz'))+[OUT/'raw/images.pdf']+list((RAW/'spatial').glob('*/tissue_hires_image.png'))
 rawhash={str(p.relative_to(ROOT)):sha(p) for p in sourcefiles};assert not oldhash.intersection(rawhash.values()),'Duplicate frozen asset'
 prior=subprocess.run(['git','grep','-l','-F','-e','CID44971','-e','CID4535','-e','zenodo.4739739','--',':!C9'],cwd=ROOT,text=True,capture_output=True);assert prior.returncode==1 and not prior.stdout
 dump(A/'logs/FRESHNESS_AUDIT.json',{'status':'PASS','prior_tracked_artifact_identity_search':'No study/selected patient identifiers found before C9','raw_asset_duplicate_hash_matches':0,'prior_manifest_sha256':sha(ROOT/'C8/C8_MANIFEST.json'),'excluded':['development NCBI783/785 TENX95/99 and related patients','all HER2ST A-H sections','GSE311609 diagnostic cohort','10x Prime breast diagnostic release'],'identity_scope':'Source study and reported patient identity, source archive checksums and available parent asset hashes; no cross-study patient reidentification claim'})
 lock={'timestamp':now(),'source_commit':commit(),'parent_C8_manifest_sha256':sha(ROOT/'C8/C8_MANIFEST.json'),'protocol':protocol,'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'image_assets_sha256':assets,'model_parents_sha256':modelparents,'raw_asset_sha256':rawhash,'calibration_and_evaluation_not_accessed':True,'registry_prefix_before_C9':{'bytes':(ROOT/'results/experiment_registry.csv').stat().st_size,'sha256':sha(ROOT/'results/experiment_registry.csv')}}
 dump(A/'C9A_SPEC_LOCK.json',lock);(A/'C9A_SPEC_LOCK.sha256').write_text(sha(A/'C9A_SPEC_LOCK.json')+'  C9A_SPEC_LOCK.json\n');(OUT/'C9A_SPEC_LOCK.json').write_bytes((A/'C9A_SPEC_LOCK.json').read_bytes());(OUT/'C9A_SPEC_LOCK.sha256').write_bytes((A/'C9A_SPEC_LOCK.sha256').read_bytes());print('LOCKED',sha(A/'C9A_SPEC_LOCK.json'),flush=True)
if __name__=='__main__':main()
