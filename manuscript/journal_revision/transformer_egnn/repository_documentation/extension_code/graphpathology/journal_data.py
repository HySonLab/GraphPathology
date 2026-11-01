"""Journal-only public-data readers with strict original identity/coordinate checks."""
import gzip
import io
import json
from pathlib import Path
import tarfile
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
from .assets import sha

WU_EVAL = ['1160920F', 'CID4290', 'CID44971']


def original_genes(root):
    return json.loads((root / 'checkpoints/legacy/C9/configs/protocol.json').read_text())['genes']


def normalized(root, raw):
    norm = json.loads((root / 'checkpoints/legacy/reports/C0/C0_NORMALIZATION.json').read_text())
    idx = [norm['genes'].index(g) for g in original_genes(root)]
    return (np.log1p(raw) - np.asarray(norm['mean'])[idx]) / np.asarray(norm['scale'])[idx]


def original_masks(root, cohort, patient=None):
    cp = {'wu': 'C9', 'bassiouni': 'C10'}[cohort]; old = root / 'checkpoints/legacy'
    records = []
    for p in sorted((old / cp / (cp + 'B') / 'masks').glob('*.json')):
        r = json.loads(p.read_text())
        if r['strategy'] != 'random' or r['seed'] not in range(10): continue
        if cohort == 'wu' and r['sample_id'] not in WU_EVAL: continue
        if patient is not None and r['sample_id'] != patient: continue
        r['absolute_mask_path'] = str(old / r['mask_path']); r['mask_sha256'] = sha(old / r['mask_path'])
        records.append(r)
    return records


def lattice_coordinates(pos):
    xy = pos[['pixel_col', 'pixel_row']].to_numpy(float)
    lattice = np.column_stack([pos.array_col.to_numpy() * 50, pos.array_row.to_numpy() * np.sqrt(3) * 50])
    design = np.column_stack([lattice, np.ones(len(pos))])
    coef = np.linalg.lstsq(design, xy, rcond=None)[0]
    mpp = 1 / float(np.mean(np.linalg.svd(coef[:2], compute_uv=False)))
    xy_um = xy * mpp
    calibrated = np.linalg.lstsq(design, xy_um, rcond=None)[0]
    residual = np.linalg.norm(design @ calibrated - xy_um, axis=1)
    singular = np.linalg.svd(calibrated[:2], compute_uv=False)
    if residual.max() >= 2 or max(abs(singular - 1)) >= .03: raise ValueError('Original lattice calibration check failed')
    return xy_um


def _plain_or_gzip(path):
    data = path.read_bytes()
    return gzip.decompress(data) if data[:2] == b'\x1f\x8b' else data


def _unpack_wu(data_root):
    target = data_root / 'raw/wu/extracted'; marker = target / '.verified_unpack.json'
    manifest = json.loads((data_root / 'DATA_ACQUISITION_MANIFEST.json').read_text())
    records = {r['local']: r for r in manifest['files']}
    needed = ['filtered_count_matrices.tar.gz', 'spatial.tar.gz']
    hashes = {}
    for name in needed:
        p = data_root / 'raw/wu' / name; rec = records.get('raw/wu/' + name, {})
        if rec.get('status') != 'AVAILABLE' or not p.exists() or sha(p) != rec['sha256']: raise FileNotFoundError('Verified Wu source missing: ' + name)
        hashes[name] = sha(p)
    if not marker.exists():
        target.mkdir(parents=True, exist_ok=True)
        for name in needed:
            with tarfile.open(data_root / 'raw/wu' / name) as archive: archive.extractall(target, filter='data')
        marker.write_text(json.dumps(hashes, sort_keys=True) + '\n')
    elif json.loads(marker.read_text()) != hashes: raise ValueError('Extraction source hashes changed')
    return target, hashes


def load_patient(root, data_root, cohort, patient):
    genes = original_genes(root); source_files = {}
    if cohort == 'wu':
        extracted, source_files = _unpack_wu(data_root)
        folders = list(extracted.rglob(patient + '_filtered_count_matrix'))
        if len(folders) != 1: raise ValueError('Ambiguous Wu count directory')
        folder = folders[0]; features = _plain_or_gzip(folder / 'features.tsv.gz').decode().splitlines()
        barcodes = np.asarray(_plain_or_gzip(folder / 'barcodes.tsv.gz').decode().splitlines())
        with threadpool_limits(limits=4): mat = mmread(io.BytesIO(_plain_or_gzip(folder / 'matrix.mtx.gz'))).tocsr()
        raw = mat[[features.index(g) for g in genes]].toarray().T.astype(float)
        positions = next(extracted.rglob(patient + '_spatial')) / 'tissue_positions_list.csv'
        pos = pd.read_csv(positions, header=None, names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pixel_row', 'pixel_col']).set_index('barcode').loc[barcodes]
        keep = (pos.in_tissue == 1).to_numpy(); pos = pos.loc[keep]; barcodes = barcodes[keep]; raw = raw[keep]
    elif cohort == 'bassiouni':
        selected = json.loads((root / 'checkpoints/legacy/C10/C10A/logs/SOURCE_SAMPLE_SELECTION.json').read_text())['selected']
        sample = next(p for p in selected if p['patient_id'] == patient)
        files = {url.split('/')[-1]: data_root / 'raw/bassiouni' / url.split('/')[-1] for url in sample['urls']}
        old_hashes = {r['url'].split('/')[-1]: r['sha256'] for r in json.loads((root / 'checkpoints/legacy/C10/C10A/logs/DOWNLOADS.json').read_text())}
        needed = [name for name in files if 'filtered_feature_bc_matrix' in name or 'tissue_positions_list' in name]
        for name in needed:
            if not files[name].exists(): raise FileNotFoundError(str(files[name]))
            if sha(files[name]) != old_hashes[name]: raise ValueError('Bassiouni source checksum mismatch')
            source_files[name] = sha(files[name])
        h5 = next(files[n] for n in needed if 'filtered_feature' in n)
        with h5py.File(h5) as f:
            m = f['matrix']; mat = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:]))
            features = [x.decode() for x in m['features/name'][:]]; barcodes = np.asarray([x.decode() for x in m['barcodes'][:]])
            raw = mat[[features.index(g) for g in genes]].toarray().T.astype(float)
        positions = next(files[n] for n in needed if 'tissue_positions_list' in n)
        pos = pd.read_csv(positions, header=None, names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pixel_row', 'pixel_col']).set_index('barcode').loc[barcodes]
        keep = (pos.in_tissue == 1).to_numpy(); pos = pos.loc[keep]; barcodes = barcodes[keep]; raw = raw[keep]
    else: raise ValueError('Unknown cohort')
    if not np.isfinite(raw).all() or (raw < 0).any() or not np.array_equal(raw, np.rint(raw)): raise ValueError('Expected unnormalized nonnegative integer counts')
    xy = lattice_coordinates(pos); y = normalized(root, raw)
    masks = original_masks(root, cohort, patient)
    if not masks: raise ValueError('No preserved random masks for patient')
    with np.load(masks[0]['absolute_mask_path']) as f:
        if not np.array_equal(barcodes, f['barcodes']): raise ValueError('Original barcode ordering mismatch')
        np.testing.assert_allclose(xy, f['xy_um'], atol=1e-8, rtol=0)
        derived_coordinate_delta = float(np.max(abs(xy - f['xy_um'])))
        # The byte-preserved, label-free original coordinates govern exact neighbor ties.
        xy = f['xy_um'].copy()
    return {'patient': patient, 'cohort': cohort, 'genes': np.asarray(genes), 'barcodes': barcodes,
            'xy': xy, 'targets': y, 'raw_counts': raw, 'source_hashes': source_files,
            'derived_coordinate_max_abs_delta': derived_coordinate_delta,
            'normalization_sha256': sha(root / 'checkpoints/legacy/reports/C0/C0_NORMALIZATION.json')}
