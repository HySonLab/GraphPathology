"""Read-only integrity check against the frozen source, artifact and data hashes."""
import argparse
import hashlib
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--include-data',action='store_true',help='Also read/hash about 6 GB of local assets')
    args=parser.parse_args()
    lock=json.loads((ROOT/'reports/C0/C0_SPEC_LOCK.json').read_text())
    expected={**lock['source_sha256'],**lock['artifact_sha256']}
    if args.include_data:
        manifest=json.loads((ROOT/'reports/C0/C0_DATA_MANIFEST.json').read_text())
        expected.update({r['path']:r['sha256'] for r in manifest['downloads']['files']})
        for slide in manifest['slides']:
            expected[slide['feature_cache']]=slide['feature_sha256']
            expected[slide['target_cache']]=slide['target_sha256']
    bad=[]
    for path,digest in expected.items():
        p=ROOT/path
        if not p.exists():
            bad.append(path+' missing'); continue
        with p.open('rb') as stream:
            actual=hashlib.file_digest(stream,'sha256').hexdigest()
        if actual!=digest: bad.append(path+' changed')
    if bad: raise SystemExit('\n'.join(bad))
    print(f'PASS: {len(expected)} frozen file hashes match')
if __name__=='__main__': main()
