"""Freeze all scientific settings before the first external expression evaluation."""
from support import *
def main():
 p=A/'C7A_SPEC_LOCK.json';assert not p.exists();assert json.loads((A/'C7A_STATUS.json').read_text())['status']=='PASS_EXTERNAL_DATASET_QUALIFICATION';assert not (B/'logs/EXTERNAL_LABEL_ACCESS.json').exists()
 md5s={}
 for rec in json.loads((OUT/'raw/HER2ST_metadata.json').read_text())['files']:
  if rec['key'] not in ['images.zip','count-matrices.zip','spot-selections.zip']:continue
  path=OUT/'raw'/rec['key']
  with path.open('rb') as f:digest=hashlib.file_digest(f,'md5').hexdigest()
  assert digest==rec['checksum'].split(':')[-1];md5s[rec['key']]=digest
 dump(A/'logs/SOURCE_ARCHIVE_MD5_AUDIT.json',{'status':'PASS','source_archives':md5s})
 code=['C7/support.py','C7/policy.py','C7/execute.py','C7/extract.py','C7/analyze.py','C7/configs/protocol.json','C7/her2_metadata.py','C7/prepare.py'];assets={}
 for s in SAMPLES:
  for suffix in ['.jpg','_metadata.npz','_patches.h5','_features.npz','_S2-50_image.npz','_S2-280_image.npz']:path=CACHE/(s+suffix);assets[str(path.relative_to(ROOT))]=sha(path)
 for name in ['images.zip','count-matrices.zip','spot-selections.zip']:assets['C7/raw/'+name]=sha(OUT/'raw'/name)
 parents={}
 for path in [ROOT/'C5/SELECTION_LOCK.json',ROOT/'reports/C0/C0_NORMALIZATION.json',ROOT/'data/weights/hub/checkpoints/resnet18-f37072fd.pth',ROOT/'C5/models'/f'PCA_primary_{ENC}.npz']+[ROOT/'C5/models'/f'B4_primary_{ENC}_{panel}_a{SEL["B4_alpha"][panel][ENC]}.npz' for panel in panels()]:parents[str(path.relative_to(ROOT))]=sha(path)
 for s in SAMPLES:
  for panel in panels():
   for r in json.loads((CACHE/(s+'_'+panel+'_image.json')).read_text())['parents']:
    parents[r['checkpoint']]=r['sha256'];parents[r['preprocessing']+'/feature.npz']=r['feature_scaler_sha256'];parents[r['preprocessing']+'/target.npz']=r['target_scaler_sha256']
 reg=ROOT/'results/experiment_registry.csv';lock={'timestamp':now(),'source_commit':commit(),'external_evaluation_authorized':True,'qualification_status':'PASS_EXTERNAL_DATASET_QUALIFICATION','protocol':cfg(),'source_sha256':{s:sha(ROOT/s) for s in code},'qualified_assets_sha256':assets,'frozen_model_and_preprocessing_sha256':parents,'parent_C6_manifest_sha256':sha(ROOT/'C6/C6_MANIFEST.json'),'registry_prefix':{'path':'results/experiment_registry.csv','bytes':reg.stat().st_size,'sha256':sha(reg)},'external_expression_values_used_for_selection':False,'qualification_limitations':['Independent breast study but different sequencing ST platform','Common 40/50 and 226/280 panel intersection','Source tissue/gene filtering inherited','One section per patient; 8 biological replicates','Only frozen ResNet18 branch']};dump(p,lock);(A/'C7A_SPEC_LOCK.sha256').write_text(sha(p)+'  C7A_SPEC_LOCK.json\n');print('C7A SPEC LOCKED',sha(p))
if __name__=='__main__':main()
