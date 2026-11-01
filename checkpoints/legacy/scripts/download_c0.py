"""Fetch revision-pinned HEST IDC assets; never requires or prints credentials."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[1]
REV = 'b039b5965ed7092f280d4f336da7972ef9a10e1a'
IDS = ['NCBI783', 'NCBI785', 'TENX95', 'TENX99']

def fetch(item):
    rel, url = item
    dest = ROOT / 'data/raw' / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        tmp = dest.with_suffix(dest.suffix + '.part')
        with requests.get(url, stream=True, timeout=(30, 120)) as response:
            response.raise_for_status()
            with tmp.open('wb') as out:
                for block in response.iter_content(1024 * 1024):
                    out.write(block)
        tmp.replace(dest)
    with dest.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    print(f'{rel}: {dest.stat().st_size} bytes', flush=True)
    return {'path': str(dest.relative_to(ROOT)), 'url': url,
            'sha256': digest, 'bytes': dest.stat().st_size}

def main():
    rev = '3ddb5eaf5bd2a8133e0c0e8015816489a3d99dc3'
    items = [('HEST_v1_1_0.csv', f'https://raw.githubusercontent.com/mahmoodlab/HEST/{rev}/assets/HEST_v1_1_0.csv')]
    for sample in IDS:
        for rel in [f'IDC/adata/{sample}.h5ad', f'IDC/patches/{sample}.h5']:
            items.append((rel, f'https://huggingface.co/datasets/MahmoodLab/hest-bench/resolve/{REV}/{rel}'))
    previous_path = ROOT / 'reports/C0/C0_DOWNLOADS.json'
    previous = json.loads(previous_path.read_text()) if previous_path.exists() else None
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(fetch, items))
    out = ROOT / 'reports/C0/C0_DOWNLOADS.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    if previous is not None and records != previous['files']:
        raise RuntimeError('Existing data differ from frozen download checksums')
    out.write_text(json.dumps({'hest_bench_revision': REV, 'hest_metadata_revision': rev, 'files': records}, indent=2) + '\n')

if __name__ == '__main__':
    main()
