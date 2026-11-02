"""Versioned artifact packaging and checksum-first, non-overwriting restoration."""
import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import shutil
import stat
import subprocess
import zipfile

from .assets import sha
from .reproduction_primary import write_json

RESTORE_PREFIXES = ('checkpoints/reproduced/', 'checkpoints/journal_revision/',
                    'experiments/journal_revision/', 'reproduction/', 'configs/reproduction/')


def member_path(root, name):
    pure = PurePosixPath(name)
    if pure.is_absolute() or '..' in pure.parts or '\\' in name or not name.startswith(RESTORE_PREFIXES):
        raise ValueError('Disallowed artifact member: ' + name)
    target = root.joinpath(*pure.parts)
    if root.resolve() not in target.resolve().parents:
        raise ValueError('Artifact member escapes checkout through a symlink')
    return target


def stream_sha(stream):
    h = hashlib.sha256()
    for block in iter(lambda: stream.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()


def verify_assets(directory, manifest):
    for asset in manifest['assets']:
        name = asset['name']
        if Path(name).name != name: raise ValueError('Asset filename must be a basename')
        path = directory / name
        if path.stat().st_size != asset['bytes'] or sha(path) != asset['sha256']:
            raise ValueError('Release asset size/checksum mismatch: ' + name)
    return {'status':'PASS','assets':len(manifest['assets']),
            'bytes':sum(a['bytes'] for a in manifest['assets'])}


def restore(root, directory, manifest):
    """Preflight every member and existing destination before writing any file."""
    verified = verify_assets(directory, manifest)
    pending = []; destinations = {}
    for asset in manifest['assets']:
        if 'files' not in asset: continue
        with zipfile.ZipFile(directory / asset['name']) as archive:
            names = [r.filename for r in archive.infolist() if not r.is_dir()]
            if len(names) != len(set(names)) or set(names) != set(asset['files']):
                raise ValueError('Archive member set differs from release file index')
            for info in archive.infolist():
                if info.is_dir(): continue
                if stat.S_ISLNK(info.external_attr >> 16): raise ValueError('Symlink members are forbidden')
                target = member_path(root, info.filename); expected = asset['files'][info.filename]
                with archive.open(info) as stream:
                    if stream_sha(stream) != expected: raise ValueError('Member checksum mismatch')
                if target in destinations and destinations[target] != expected:
                    raise ValueError('Conflicting duplicate member across assets')
                destinations[target] = expected
                if target.exists():
                    if not target.is_file() or sha(target) != expected:
                        raise FileExistsError('Refusing to overwrite existing content: ' + str(target))
                else: pending.append((asset['name'], info.filename, target))
    # All preconditions have passed. Existing identical files are left byte-identical.
    written = 0
    for name, member, target in pending:
        if target.exists(): continue  # Identical shared member already restored by another archive.
        target.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(directory / name) as archive, archive.open(member) as source, target.open('xb') as dest:
            shutil.copyfileobj(source, dest)
        written += 1
    return {**verified,'restored_files':written,'existing_identical_or_shared_files':len(destinations)-written}


def preservation(root):
    baseline=json.loads((root/'reproduction/FROZEN_FILE_SHA256.json').read_text())
    failures=[name for name,digest in baseline['files'].items() if not (root/name).is_file() or sha(root/name)!=digest]
    base=baseline['base_commit'];trees={}
    for name in ['manuscript','source_data','checkpoints','experiments','archive','archives','results','figures','tables']:
        before=subprocess.check_output(['git','ls-tree',base,'--',name],cwd=root,text=True).strip()
        after=subprocess.check_output(['git','ls-tree','HEAD','--',name],cwd=root,text=True).strip()
        trees[name]={'before':before,'after':after,'identical':before==after}
    original_changes=subprocess.check_output(['git','diff','--name-only','--diff-filter=MD',base,'--'],cwd=root,text=True).splitlines()
    forbidden_original_changes=[name for name in original_changes if name!='.gitignore']
    return {'status':'PASS' if not failures and not forbidden_original_changes and all(t['identical'] for t in trees.values()) else 'FAIL',
            'baseline_tracked_modified_paths':original_changes,'allowed_original_changes':['.gitignore'],
            'base_commit':base,'checked_files':len(baseline['files']),'differences':failures,'git_trees':trees,
            'manuscript_rebuilt_in_place':False,'historical_results_overwritten':False}


def add_archive(assets, path, entries):
    """Entries map restored path to bytes or an existing file; preserve member bytes."""
    files={}
    with zipfile.ZipFile(path,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for name, value in sorted(entries.items()):
            member_path(Path('/tmp/reproduction-checkout'),name)
            payload=value if isinstance(value,bytes) else value.read_bytes()
            info=zipfile.ZipInfo(name,date_time=(2026,11,1,0,0,0));info.external_attr=0o100644<<16
            info.compress_type=zipfile.ZIP_STORED if name.endswith(('.gz','.zip')) else zipfile.ZIP_DEFLATED
            z.writestr(info,payload);files[name]=hashlib.sha256(payload).hexdigest()
    assets.append({'name':path.name,'bytes':path.stat().st_size,'sha256':sha(path),'files':files})


def build(root, destination, version='v1'):
    if destination.exists(): raise FileExistsError('Fresh release staging directory required')
    if root in destination.parents and destination.relative_to(root).parts[0] not in ['work','outputs']:
        raise ValueError('Stage release assets outside tracked scientific directories')
    public=destination/'public';held=destination/'held-phikon';public.mkdir(parents=True);held.mkdir()
    index=json.loads((root/'reproduction/MODEL_INDEX.json').read_text())
    training=root/'checkpoints/reproduced/primary_resnet18/primary-20261101-v1'
    prefix='reproduction/release-'+version+'/'
    common={prefix+'LICENSE_ARTIFACTS.txt':root/'reproduction/licenses/CC-BY-NC-SA-4.0.txt',
            prefix+'UPSTREAM_NOTICES.md':root/'reproduction/LICENSE_AUDIT.md'}
    assets=[]
    entries={str(p.relative_to(root)):p for p in training.iterdir() if p.is_file()};entries.update(common)
    add_archive(assets,public/f'graphpathology-retrained-primary-resnet18-heads-{version}.zip',entries)
    for backbone,folder in [('resnet18',public),('phikon',held)]:
        models=[m for m in index['models'] if m['backbone']==backbone]
        entries={}
        for model in models:
            path=root/model['checkpoint']
            if sha(path)!=model['sha256']:raise ValueError('Model changed after verification')
            entries[model['checkpoint']]=path
            scaler=root/model['feature_scaler'];entries[model['feature_scaler']]=scaler
            if path.name=='best.pt':
                for sibling in path.parent.iterdir():
                    if sibling.is_file():entries[str(sibling.relative_to(root))]=sibling
        subset={**index,'models':models,'unique_best_states':len(models),
                'logical_references':sum(len(m['references']) for m in models)}
        entries[prefix+'MODEL_INDEX_'+backbone+'.json']=(json.dumps(subset,indent=2)+'\n').encode()
        entries.update(common if backbone=='resnet18' else {prefix+'OWKIN-Phikon-LICENSE.txt':root/'reproduction/licenses/OWKIN-Phikon-LICENSE.txt',prefix+'RIGHTS_REVIEW.md':root/'reproduction/LICENSE_AUDIT.md'})
        target=folder/(f'graphpathology-exploratory-neural-heads-{version}.zip' if backbone=='resnet18' else f'graphpathology-phikon-heads-RIGHTS-BLOCKED-{version}.zip')
        target_assets=assets if backbone=='resnet18' else []
        add_archive(target_assets,target,entries)
        if backbone=='phikon':write_json(held/'HELD_MODEL_MANIFEST.json',{'status':'NOT_UPLOADED_RIGHTS_UNCLEAR','assets':target_assets})
    # Source exports are copied byte-for-byte; no raw images/counts/encoders are selected.
    entries={}
    sources=[(root.parent/'GraphPathology_Transformer_EGNN_Full_Gene_Exports.zip',False),
             (root.parent/'GraphPathology_Matched_Backbone_Numerical_Results.zip',True)]
    for source,filter_matched in sources:
        with zipfile.ZipFile(source) as z:
            for name in z.namelist():
                if name.endswith('/'):continue
                if filter_matched and not any('journal_v2_' in part or 'resnet18_' in part for part in PurePosixPath(name).parts):continue
                if 'phikon' in name.lower() or not name.endswith(('.csv','.csv.gz','.json')):
                    raise ValueError('Unexpected source-data release member: '+name)
                payload=z.read(name)
                if name in entries and entries[name]!=payload:raise ValueError('Conflicting recovered export')
                entries[name]=payload
    # New rerun metrics, histories and primary provenance are separate from historical paths.
    for folder in ['tier-a-v1','tier-b-v2']:
        for path in (root/'reproduction/runs'/folder).rglob('*'):
            if path.is_file() and path.suffix in ['.json','.csv','.gz']:
                entries[str(path.relative_to(root))]=path
    # Phikon-derived binary/result archives remain local. Export only the ResNet matched rows here.
    import pandas as pd
    p=pd.read_csv(root/'reproduction/runs/matched-v2/PANEL.csv.gz')
    entries[prefix+'MATCHED_RESNET18_REPLAY.csv']=p[p.backbone.eq('resnet18')].to_csv(index=False).encode()
    entries.update(common)
    add_archive(assets,public/f'graphpathology-source-data-results-{version}.zip',entries)
    # Original mixed matched numerical archive remains intact in local rights-review staging.
    mixed=root.parent/'GraphPathology_Matched_Backbone_Numerical_Results.zip'
    shutil.copyfile(mixed,held/mixed.name)
    for source in ['REPRODUCTION_REPORT.md','RESULTS_COMPARISON.csv','README_RELEASE.md','LICENSE_AUDIT.md','requirements-v1.txt','ENVIRONMENT.json','REPLAY_PROTOCOL.md','PRESERVATION_REPORT.json','VALIDATION.json']:
        shutil.copyfile(root/'reproduction'/source,public/source)
        path=public/source;assets.append({'name':source,'bytes':path.stat().st_size,'sha256':sha(path)})
    manifest={'schema_version':1,'version':version,'tag':'research-artifacts-'+version,
              'code_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
              'original_main':'67738be6f0aaace07b55afac44db3d3d23c41c9b','code_is_on_unmerged_pr_branch':True,
              'assets':assets,'artifact_license':'CC-BY-NC-SA-4.0; see scoped LICENSE_AUDIT.md',
              'original_primary_recovered':False,'replacement_primary_heads':3,'recovered_public_best_heads':96,
              'additional_resume_files':48,'phikon_best_heads_held_locally':36,
              'pretrained_encoders_included':False,'raw_patient_data_included':False,
              'dependencies':json.loads((root/'reproduction/ENVIRONMENT.json').read_text()),
              'training_environment':json.loads((training/'MANIFEST.json').read_text())['environment'],
              'excluded_from_self_hash':'This manifest, checksum index and later verification receipt cannot hash themselves; upload receipts record their GitHub digests.'}
    manifest_path=public/f'graphpathology-model-and-results-manifest-{version}.json';write_json(manifest_path,manifest)
    checks=[f"{a['sha256']}  {a['name']}" for a in assets]+[f'{sha(manifest_path)}  {manifest_path.name}']
    (public/'SHA256SUMS.txt').write_text('\n'.join(checks)+'\n')
    verify_assets(public,manifest)
    return {'public_assets':len(assets)+2,'public_bytes':sum(p.stat().st_size for p in public.iterdir()),
            'manifest':str(manifest_path),'held_directory':str(held)}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('build');q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('preserve');q.add_argument('--output',type=Path,required=True)
    for action in ['verify','restore']:
        q=sub.add_parser(action);q.add_argument('--assets',type=Path,required=True);q.add_argument('--manifest',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve()
    if a.action=='build':result=build(root,a.output.resolve())
    elif a.action=='preserve':
        target=a.output.resolve()
        if root not in target.parents or target.relative_to(root).parts[0] not in ['work','outputs','reproduction'] or target.exists():
            raise ValueError('Preservation output must be a fresh work/outputs/reproduction report')
        result=preservation(root);write_json(target,result)
        if result['status']!='PASS':raise SystemExit('Scientific preservation failed')
    else:
        manifest=json.loads(a.manifest.read_text())
        result=restore(root,a.assets.resolve(),manifest) if a.action=='restore' else verify_assets(a.assets.resolve(),manifest)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
