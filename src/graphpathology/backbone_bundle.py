"""Author-local review bundles; never upload or redistribute the frozen encoder."""
import json
import zipfile
from pathlib import Path
from .assets import sha
from .backbone_data import STUDY,record,write


def weights(root,output):
    base=root/f'experiments/journal_revision/{STUDY}';rows=json.loads((base/'ALL_CANDIDATE_FIT_MANIFEST.json').read_text());selected=json.loads((base/'SELECTED_FIT_MANIFEST.json').read_text())
    assert len(rows)==72 and len(selected)==36
    files={}
    for r in rows:
        for name,digest in [('checkpoint',r['checkpoint_sha256']),('feature_scaler',r['feature_scaler_sha256'])]:
            path=root/r[name];assert sha(path)==digest;files[str(path.relative_to(root))]=path
    for name in ['PROTOCOL.md','PROTOCOL_V2_TRAINING_REFINEMENT.md','RESOURCE_ADJUSTMENT.md','ALL_CANDIDATE_FIT_MANIFEST.json','SELECTED_FIT_MANIFEST.json','VALIDATION_SELECTION.json','ENVIRONMENT.json','ENCODER_ACQUISITION.json']:
        path=base/name;files[str(path.relative_to(root))]=path
    config=root/f'configs/journal_revision/{STUDY}/protocol.json';files[str(config.relative_to(root))]=config
    license_path=root/f'data/journal_revision/{STUDY}/encoder/LICENSE.txt';files['licenses/PHIKON_LICENSE.txt']=license_path
    manifest=[{'path':name,'bytes':path.stat().st_size,'sha256':sha(path)} for name,path in sorted(files.items())]
    output=Path(output).resolve();output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,path in sorted(files.items()):z.write(path,name)
        z.writestr('BUNDLE_MANIFEST.json',json.dumps({'candidate_heads':72,'validation_selected_heads':36,'files':manifest,'pretrained_encoder_included':False,'original_missing_primary_heads_included':False,'availability':'Author-local review bundle; no public deposition.'},indent=2)+'\n')
        z.writestr('README.txt','Author-review only. This bundle contains72 candidate head states and36 validation-selected head identifiers, source-only feature scalers and provenance. The pretrained encoders and raw research inputs are excluded. The Phikon license notice is included; distribution permissions require author review. Original three primary image-head binaries remain missing. Extract at the repository root to restore the documented relative paths. No release or DOI has been created.\n')
    with zipfile.ZipFile(output) as z:assert z.testzip() is None
    result=record(root,status='COMPLETE',bundle='author_local_weights',path=str(output),bytes=output.stat().st_size,sha256=sha(output),candidate_heads=72,selected_heads=36,pretrained_encoder_included=False,files=manifest)
    write(base/'WEIGHT_BUNDLE_MANIFEST.json',result);return result


def numerical(root,output):
    base=root/f'experiments/journal_revision/{STUDY}';analysis=json.loads((base/'ANALYSIS_STATUS.json').read_text());assert analysis['status']=='COMPLETE'
    rows=json.loads((base/'PRIMARY_EVALUATION_MANIFEST.json').read_text());assert len(rows)==936
    files={}
    for r in rows:
        for key in ['panel','genes']:
            p=root/r[key];assert sha(p)==r[key+'_sha256'];files[str(p.relative_to(root))]=p
        if r['status']=='COMPLETE':
            p=(root/r['panel']).parent/'COMPLETE.json';files[str(p.relative_to(root))]=p
    for patient in sorted({r['patient'] for r in rows}):
        for name in ['MOLECULAR.csv','MOLECULAR_GENES.csv.gz']:
            p=root/'experiments/journal_revision/transformer_egnn/evaluation'/patient/name;assert p.is_file();files[str(p.relative_to(root))]=p
    for directory in ['summary','fixed_lr1e3_resnet_summary','supports','fits','trials']:
        for path in (base/directory).rglob('*'):
            if path.is_file():files[str(path.relative_to(root))]=path
    for path in base.glob('*.json'):
        if path.name!='NUMERICAL_BUNDLE_MANIFEST.json':files[str(path.relative_to(root))]=path
    for path in base.glob('*.md'):files[str(path.relative_to(root))]=path
    manifest=[{'path':name,'bytes':path.stat().st_size,'sha256':sha(path)} for name,path in sorted(files.items())]
    output=Path(output).resolve();output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,path in sorted(files.items()):z.write(path,name)
        z.writestr('NUMERICAL_BUNDLE_MANIFEST.json',json.dumps({'files':manifest,'patient_model_bundles':936,'publication_status':'Author-local; not a public release.'},indent=2)+'\n')
    with zipfile.ZipFile(output) as z:assert z.testzip() is None
    result=record(root,status='COMPLETE',bundle='author_local_numerical',path=str(output),bytes=output.stat().st_size,sha256=sha(output),patient_model_bundles=936,files=manifest)
    write(base/'NUMERICAL_BUNDLE_MANIFEST.json',result);return result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('kind',choices=['weights','numerical']);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=(weights if a.kind=='weights' else numerical)(a.root.resolve(),a.output);print(json.dumps({k:v for k,v in r.items() if k!='files'},indent=2))
