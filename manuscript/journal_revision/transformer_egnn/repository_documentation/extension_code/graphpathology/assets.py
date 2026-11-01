"""Read-only migration resolution and integrity checks."""
from pathlib import Path
import csv, hashlib, json

def sha(path, length=None):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        if length is not None: h.update(f.read(length))
        else:
            for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def inventory(root):
    return list(csv.DictReader((Path(root)/'checkpoints/manifest.csv').open()))

def resolve(root, old_path):
    root=Path(root).resolve()
    p=Path(old_path)
    if p.is_absolute() or '..' in p.parts: raise ValueError('Expected repository-relative historical path')
    # Missing model installation paths are also in the inventory.
    for row in inventory(root):
        if row['original_path']==old_path: return root/row['path']
    raise FileNotFoundError(f'Historical path not inventoried: {old_path}')

def verify(root, historical=False):
    root=Path(root).resolve(); counts={'checked':0,'changed':0,'missing':0};details=[]
    paths=list(csv.DictReader((root/'checkpoints/path_map.csv').open()))
    for row in paths:
        p=root/row['new_path']; counts['checked']+=1
        if not p.is_file(): counts['missing']+=1;details.append({'path':row['old_path'],'status':'migration_missing'})
        elif sha(p)!=row['sha256']:counts['changed']+=1;details.append({'path':row['old_path'],'status':'migration_changed'})
    out={'migration':counts,'migration_issues':details}
    if historical:
        mapping={x['old_path']:root/x['new_path'] for x in paths}
        results=[]
        for row in csv.DictReader((root/'checkpoints/historical_expectations.csv').open()):
            p=mapping.get(row['path']); status='missing'
            if p is not None and p.is_file():
                actual=sha(p,int(row['prefix_bytes']) if row['scope']=='prefix' else None)
                status='match' if actual==row['expected_sha256'] else 'pre_existing_mismatch'
            results.append({**row,'current_status':status,'changed_from_start':status!=row['starting_status']})
        from collections import Counter
        out['historical_counts']=dict(Counter(x['current_status'] for x in results))
        out['historical_new_issues']=[x for x in results if x['changed_from_start']]
        out['historical_details']=results
    return out

def new_output(root, output):
    root=Path(root).resolve(); out=Path(output).resolve()
    if out==root or root in out.parents and out.relative_to(root).parts[0] not in ['outputs','work']:
        raise ValueError('Use a new outputs/ directory; archive, manuscript and source writes are forbidden')
    if out.exists(): raise FileExistsError(f'Refusing to overwrite output: {out}')
    out.mkdir(parents=True)
    return out
