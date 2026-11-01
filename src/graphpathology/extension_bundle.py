"""Release-ready author bundles; no public upload, credentials, raw data or Git binaries."""
import json,zipfile,datetime
from pathlib import Path
from .assets import sha
from .extension_evaluate import model_inventory


def run(root,destination):
 base=root/'experiments/journal_revision/transformer_egnn';destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
 # The unchanged relative paths let the author extract the full exports directly into a checkout.
 files=sorted(list((base/'evaluation').glob('*/*/GENES.csv.gz'))+list((base/'evaluation').glob('*/MOLECULAR_GENES.csv.gz'))+list((base/'bassiouni_anchors').glob('*/GENE_METRICS.csv.gz')))
 assert len(files)==2028+26+22
 metric=destination/'GraphPathology_Transformer_EGNN_Full_Gene_Exports.zip'
 with zipfile.ZipFile(metric,'w',compression=zipfile.ZIP_STORED) as z:
  for p in files:z.write(p,p.relative_to(root))
 weights=destination/'GraphPathology_Transformer_EGNN_New_Heads.zip';inventory=model_inventory(root)
 with zipfile.ZipFile(weights,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for item in inventory:
   if item['origin']!='new':continue
   p=root/item['checkpoint'];assert sha(p)==item['checkpoint_sha256']
   for name in ['best.pt','resume.pt']:
    f=p.parent/name;z.write(f,f.relative_to(root))
   f=root/item['fit_record'];z.write(f,f.relative_to(root))
  f=root/'checkpoints/journal_revision/hest_learned/journal_v2_feature_scaler.npz';z.write(f,f.relative_to(root))
 r={'status':'RELEASE_READY_AUTHOR_LOCAL','label':'post_v1_exploratory','public_upload_performed':False,'raw_datasets_included':False,'original_missing_weights_recovered':False,
  'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'author_directory':str(destination.resolve()),'bundles':[{'filename':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in [metric,weights]],
  'full_gene_files':[{'path':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files],
  'models':inventory,'instructions':'Author may review local bundles and approve a licensed GitHub release/DOI destination. No public release is claimed. Extract full gene exports into repository root; new-head bundle contains48 best/optimizer states plus resumable final states and the unchanged source feature scaler. The30 older heads remain a distinct previously prepared journal-v2 bundle.'}
 (base/'RELEASE_READY_MANIFEST.json').write_text(json.dumps(r,indent=2)+'\n');return r
