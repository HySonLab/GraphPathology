"""Explicit, resumable, checksum-verified journal data acquisition; no historical writes."""
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import ftplib
from urllib.parse import urlparse
import time
import requests
from .assets import sha

PROTOCOL_COMMIT = '2e70b30722b56317e82005999ed90d8faeb2537a'


def fetch_file(task, data_root, quota_gib=20):
    dest = data_root / task['local']; dest.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic(); errors = []
    for attempt in range(3):
        try:
            if not dest.exists():
                part = dest.with_name(dest.name + '.part'); offset = part.stat().st_size if part.exists() else 0
                url = task.get('download_url', task['url'])
                if url.startswith('ftp://'):
                    parsed = urlparse(url)
                    with ftplib.FTP(parsed.hostname, timeout=45) as ftp:
                        ftp.login()
                        with part.open('ab' if offset else 'wb') as f:
                            ftp.retrbinary('RETR ' + parsed.path, f.write, blocksize=2**20, rest=offset or None)
                else:
                    headers = {'Range': f'bytes={offset}-'} if offset else {}
                    with requests.get(url, headers=headers, stream=True, timeout=(20, 90)) as response:
                        response.raise_for_status()
                        append = offset and response.status_code == 206
                        with part.open('ab' if append else 'wb') as f:
                            for chunk in response.iter_content(2**20):
                                f.write(chunk)
                                if f.tell() > quota_gib * 2**30: raise ValueError('Single-file download quota exceeded')
                part.rename(dest)
            digest = sha(dest)
            if task.get('sha256_prefix') and not digest.startswith(task['sha256_prefix']): raise ValueError('Source SHA256 prefix mismatch')
            if task.get('sha256') and digest != task['sha256']: raise ValueError('SHA256 mismatch: ' + str(dest))
            if task.get('md5'):
                with dest.open('rb') as f: md5 = hashlib.file_digest(f, 'md5').hexdigest()
                if md5 != task['md5']: raise ValueError('Source MD5 mismatch: ' + str(dest))
            result = {**task, 'status': 'AVAILABLE', 'sha256': digest, 'bytes': dest.stat().st_size,
                      'seconds': time.monotonic() - started, 'errors': errors}
            print('ACQUIRED', task['local'], result['bytes'], flush=True)
            return result
        except (requests.RequestException, OSError, ValueError, ftplib.Error) as exc:
            errors.append(str(exc))
            # A complete bad digest is not resumable; retain it visibly rather than silently deleting/replacing it.
            if dest.exists(): break
    print('BLOCKED', task['local'], errors[-1], flush=True)
    return {**task, 'status': 'BLOCKED', 'errors': errors, 'seconds': time.monotonic() - started}


def plan(root, cohorts):
    old = root / 'checkpoints/legacy'; tasks = []
    if 'bassiouni' in cohorts:
        selected = json.loads((old / 'C10/C10A/logs/SOURCE_SAMPLE_SELECTION.json').read_text())['selected']
        original = json.loads((old / 'C10/C10A/logs/DOWNLOADS.json').read_text())
        hashes = {r['url']: r['sha256'] for r in original}
        for patient in selected:
            for url in patient['urls']:
                if any(s in url for s in ['filtered_feature_bc_matrix', 'tissue_positions_list', 'scalefactors_json']):
                    tasks.append({'cohort': 'bassiouni', 'patient': patient['patient_id'], 'url': url,
                                  'local': 'raw/bassiouni/' + url.split('/')[-1], 'sha256': hashes[url],
                                  'download_url': url.replace('https://', 'ftp://'),
                                  'license': 'Public GEO research data; NCBI imposes no restrictions, subject to submitter rights; no raw redistribution',
                                  'license_source': 'https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html'})
    if 'wu' in cohorts:
        metadata = requests.get('https://zenodo.org/api/records/4739739', timeout=45); metadata.raise_for_status(); record = metadata.json()
        if record['metadata'].get('license', {}).get('id') != 'cc-by-4.0': raise ValueError('Wu license changed; inspect before download')
        for f in record['files']:
            if f['key'] in ['filtered_count_matrices.tar.gz', 'spatial.tar.gz', 'metadata.tar.gz']:
                tasks.append({'cohort': 'wu', 'url': f['links']['self'], 'local': 'raw/wu/' + f['key'],
                              'md5': f['checksum'].split(':')[1], 'license': 'CC-BY-4.0', 'license_source': 'https://zenodo.org/records/4739739'})
    if 'hest' in cohorts:
        original = json.loads((old / 'reports/C0/C0_DOWNLOADS.json').read_text())
        for f in original['files']:
            if '/IDC/' in f['path'] or f['path'].endswith('HEST_v1_1_0.csv'):
                suffix = f['path'].split('data/raw/', 1)[-1]
                tasks.append({'cohort': 'hest', 'url': f['url'], 'local': 'raw/hest/' + suffix,
                              'sha256': f['sha256'], 'license': 'CC-BY-NC-SA-4.0; academic research only; raw data not redistributed',
                              'license_source': 'https://github.com/mahmoodlab/HEST/blob/main/LICENSE.md'})
    if 'her2st' in cohorts:
        response = requests.get('https://zenodo.org/api/records/4751624', timeout=45); response.raise_for_status(); record = response.json()
        if record['metadata'].get('license', {}).get('id') != 'cc-by-4.0': raise ValueError('HER2ST license changed')
        for f in record['files']:
            if f['key'] in ['count-matrices.zip', 'spot-selections.zip']:
                tasks.append({'cohort': 'her2st', 'url': f['links']['self'], 'local': 'raw/her2st/' + f['key'],
                              'md5': f['checksum'].split(':')[1], 'license': 'CC-BY-4.0', 'license_source': 'https://zenodo.org/records/4751624'})
    if 'resnet18' in cohorts:
        tasks.append({'cohort': 'pretrained_encoder', 'url': 'https://download.pytorch.org/models/resnet18-f37072fd.pth',
                      'local': 'models/resnet18-f37072fd.pth', 'sha256_prefix': 'f37072fd', 'license': 'Official torchvision ImageNet1K_V1 weights; research use under upstream/ImageNet terms, not redistributed',
                      'license_source': 'https://docs.pytorch.org/vision/stable/models/generated/torchvision.models.resnet18.html'})
    return tasks


def acquire(root, data_root, cohorts, dry_run=False):
    data_root = Path(data_root).resolve(); tasks = plan(root, cohorts)
    if dry_run: return {'tasks': tasks, 'downloads_started': False}
    data_root.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(data_root).free < 20 * 2**30: raise ValueError('Need at least 20 GiB disk free before acquisition')
    path = data_root / 'DATA_ACQUISITION_MANIFEST.json'
    previous = json.loads(path.read_text()).get('files', []) if path.exists() else []
    with ThreadPoolExecutor(max_workers=2) as pool: results = list(pool.map(lambda t: fetch_file(t, data_root), tasks))
    combined = {r['local']: r for r in previous}
    combined.update({r['local']: r for r in results})
    result = {'protocol_commit': PROTOCOL_COMMIT, 'files': list(combined.values()), 'data_root': str(data_root),
              'raw_data_committed': False, 'original_states_recovered': False}
    path.write_text(json.dumps(result, indent=2) + '\n')
    return result
