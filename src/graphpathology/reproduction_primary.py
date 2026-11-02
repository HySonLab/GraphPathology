"""Explicit replacement-primary workflow; never resolves missing historical identities.

Run with ``python -m graphpathology.reproduction_primary --help``. Scientific inputs
and the frozen historical inference commands are read-only.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess
import time

import h5py
import numpy as np
import torch
from threadpoolctl import threadpool_limits

from .assets import sha
from .core import Scaler, PeakRSS, mlp_model, predict, train_mlp, gene_metrics, summarize
from .journal_data import normalized, original_genes
from .journal_learned import load_prepared

ORIGIN = 'RETRAINED_NOT_ORIGINAL'
SCALERS = {
    'feature_scaler': '6335381832baf11e87366152a9942a8167ba988891faca99e240cebd83851820',
    'target_scaler': '53c082753b4756db62f28fb2ef91f010d84359979245ed1a7eb89edbf4e001a4',
}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def load_scaler(path):
    scaler = Scaler()
    with np.load(path, allow_pickle=False) as f:
        scaler.mean, scaler.scale = f['mean'].copy(), f['scale'].copy()
    if not np.isfinite(scaler.mean).all() or not np.isfinite(scaler.scale).all() or (scaler.scale <= 0).any():
        raise ValueError('Invalid source scaler')
    return scaler


def checked_namespace(root, path, namespace):
    root, path = root.resolve(), path.resolve()
    allowed = (root / namespace).resolve()
    if path == allowed or allowed not in path.parents:
        raise ValueError('Output must be a fresh child of ' + str(allowed))
    return path


def config(root, path):
    cfg = json.loads(path.read_text())
    if cfg['origin'] != ORIGIN or cfg['train'] != ['NCBI783', 'NCBI785'] or cfg['validation'] != ['TENX95']:
        raise ValueError('Historical source split/origin cannot be changed in this route')
    if cfg['seeds'] != [0, 1, 2] or cfg['model'] != {'hidden': 256, 'dropout': .1}:
        raise ValueError('This route implements the locked three-seed primary, not a model search')
    expected = {'learning_rate': .001, 'weight_decay': .0001, 'batch_size': 512,
                'max_epochs': 80, 'patience': 10, 'min_delta': .0001}
    if cfg['training'] != expected or cfg['threads'] != 4:
        raise ValueError('Training parameters differ from the locked primary protocol')
    for key, digest in SCALERS.items():
        if sha(root / cfg[key]) != digest:
            raise ValueError('Historical scaler identity mismatch: ' + key)
    return cfg


def audit_inputs(root, cfg):
    """Verify raw count identities and reconstruct the *float32* C0 transform.

    The archived CountTransform.transform casts to float32 before the separate
    head-target scaler. Replacing this with a float64 intermediate is different.
    """
    base = root / 'data/journal_revision/hest_training'
    names = cfg['train'] + cfg['validation']
    data = load_prepared(base, names)
    records = []
    genes = original_genes(root)
    for sample in names:
        cache = base / 'processed/journal_v2_hest' / (sample + '.npz')
        meta = json.loads(cache.with_suffix('.json').read_text())
        expression = base / 'raw/hest/IDC/adata' / (sample + '.h5ad')
        patches = base / 'raw/hest/IDC/patches' / (sample + '.h5')
        for path, key in [(expression, 'source_expression_sha256'), (patches, 'source_patches_sha256')]:
            if sha(path) != meta[key]:
                raise ValueError('Raw source identity mismatch: ' + str(path))
        with h5py.File(expression) as f, h5py.File(patches) as p:
            ids = f['obs/_index'][:].astype(str)
            patch_ids = p['barcode'][:].reshape(-1).astype(str)
            np.testing.assert_array_equal(patch_ids, data[sample]['barcodes'])
            lookup = {v: i for i, v in enumerate(ids)}
            rows = np.array([lookup[v] for v in patch_ids])
            columns = [list(f['var/_index'][:].astype(str)).index(g) for g in genes]
            raw = f['X'][:][rows][:, columns].astype(np.float64)
            expected = normalized(root, raw).astype(np.float32)
            np.testing.assert_array_equal(expected, data[sample]['targets'])
            np.testing.assert_array_equal(data[sample]['genes'], genes)
            if p['img'].shape != (len(expected), 224, 224, 3):
                raise ValueError('Unexpected source RGB patch geometry')
            error = float(np.max(abs(p['coords'][:][:, ::-1] - f['obsm/spatial'][:][rows])))
            if error > 1.01:
                raise ValueError('Patch/count coordinate mismatch')
        if meta['encoder_sha256'] != 'f37072fd47e89c5e827621c5baffa7500819f7896bbacec160b1a16c560e07ec':
            raise ValueError('Unexpected pretrained encoder')
        records.append({'sample': sample, 'locations': len(expected), 'cache_sha256': sha(cache),
                        'expression_sha256': sha(expression), 'patches_sha256': sha(patches),
                        'barcodes_sha256': __import__('hashlib').sha256(patch_ids.tobytes()).hexdigest(),
                        'patch_coordinate_error_px': error, 'targets_exact_float32_C0': True})
    sx, sy = (load_scaler(root / cfg[k]) for k in ['feature_scaler', 'target_scaler'])
    deltas = {}
    for key, values, old in [('feature', np.concatenate([data[s]['embeddings'] for s in cfg['train']]), sx),
                             ('target', np.concatenate([data[s]['targets'] for s in cfg['train']]), sy)]:
        rebuilt = Scaler().fit(values)
        deltas[key] = {}
        for field in ['mean', 'scale']:
            a, b = getattr(rebuilt, field), getattr(old, field)
            np.testing.assert_allclose(a, b, atol=1e-6, rtol=1e-6)
            deltas[key][field + '_max_abs_delta'] = float(np.max(abs(a - b)))
    return data, {'status': 'PASS', 'inputs': records, 'source_scaler_reconstruction_deltas': deltas,
                  'scalers_used': SCALERS, 'target_precision': 'C0 float32, then archived head-target scaler to float32',
                  'feature_provenance': 'Recovered, independently recomputed frozen official ResNet18 caches; historical cache containers are unavailable.'}


def load_heads(folder):
    manifest = json.loads((folder / 'MANIFEST.json').read_text())
    if manifest['origin'] != ORIGIN or [r['seed'] for r in manifest['heads']] != [0,1,2]:
        raise ValueError('Replacement identity/ensemble size mismatch')
    heads = []
    for item in manifest['heads']:
        path = folder / item['file']
        if sha(path) != item['sha256']:
            raise ValueError('Replacement checkpoint checksum mismatch')
        payload = torch.load(path, map_location='cpu', weights_only=True)
        if payload['origin'] != ORIGIN or payload['seed'] != item['seed'] or payload['genes'] != manifest['genes']:
            raise ValueError('Replacement metadata mismatch')
        net = mlp_model(512, 50, manifest['config']['model'])
        net.load_state_dict(payload['state_dict'], strict=True)
        heads.append(net.eval())
    for name, digest in manifest['scalers'].items():
        if sha(folder / name) != digest:
            raise ValueError('Replacement scaler checksum mismatch')
    return manifest, heads, load_scaler(folder / 'feature.npz'), load_scaler(folder / 'target.npz')


def image_predictions(folder, features):
    features=np.asarray(features)
    if features.ndim!=2 or features.shape[1]!=512 or len(features)==0 or not np.isfinite(features).all():
        raise ValueError('Finite nonempty N by512 frozen feature matrix required')
    _, heads, sx, sy = load_heads(folder)
    x = sx.transform(features)
    individual = [sy.inverse(predict(h, x, 'B2')) for h in heads]
    return np.mean(individual, axis=0).astype(np.float32), individual


def train(root, cfg_path, folder, output):
    folder = checked_namespace(root, folder, 'checkpoints/reproduced/primary_resnet18')
    output = checked_namespace(root, output, 'reproduction/runs')
    if folder.exists() or output.exists():
        raise FileExistsError('Fresh run identities required; refusing to replace any existing run')
    cfg = config(root, cfg_path)
    folder.mkdir(parents=True); output.mkdir(parents=True)
    started = time.monotonic()
    data, audit = audit_inputs(root, cfg)
    write_json(output / 'INPUT_AUDIT.json', audit)
    for key, name in [('feature_scaler', 'feature.npz'), ('target_scaler', 'target.npz')]:
        shutil.copyfile(root / cfg[key], folder / name)
    sx, sy = load_scaler(folder / 'feature.npz'), load_scaler(folder / 'target.npz')
    x = sx.transform(np.concatenate([data[s]['embeddings'] for s in cfg['train']]))
    y = sy.transform(np.concatenate([data[s]['targets'] for s in cfg['train']]))
    xv = sx.transform(data['TENX95']['embeddings']); yv = sy.transform(data['TENX95']['targets'])
    torch.set_num_threads(cfg['threads']); torch.use_deterministic_algorithms(True)
    manifest = {'origin': ORIGIN, 'created_utc': datetime.now(timezone.utc).isoformat(),
                'code_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
                'code_hashes': {str(p.relative_to(root)): sha(p) for p in [Path(__file__).resolve(), root / 'src/graphpathology/core.py']},
                'protocol_commit': '1b19cf9d9', 'config': cfg, 'config_sha256': sha(cfg_path),
                'genes': original_genes(root), 'heads': [], 'input_audit': audit,
                'scalers': {n: sha(folder / n) for n in ['feature.npz', 'target.npz']},
                'environment': {'python': platform.python_version(), 'torch': str(torch.__version__),
                                'numpy': np.__version__, 'platform': platform.platform(), 'logical_cpu': __import__('os').cpu_count(),
                                'device': 'cpu', 'dtype': 'float32', 'threads': cfg['threads'], 'deterministic_algorithms': True,
                                'cuda_available': torch.cuda.is_available()},
                'ensemble': 'Inverse head target scaler per seed; float64 arithmetic mean then float32 cast; C9 transport only once at external fusion.'}
    with threadpool_limits(cfg['threads']), PeakRSS() as memory:
        for seed in cfg['seeds']:
            net, fit = train_mlp(x, y, xv, yv, cfg['model'], cfg['training'], seed)
            name = 'replacement_seed' + str(seed) + '.pt'
            torch.save({'origin': ORIGIN, 'seed': seed, 'genes': manifest['genes'],
                        'config': cfg['model'], 'state_dict': net.state_dict()}, folder / name)
            validation = summarize(gene_metrics(data['TENX95']['targets'], sy.inverse(predict(net, xv, 'B2'))))
            record = {'seed': seed, 'file': name, 'sha256': sha(folder / name), 'fit': fit,
                      'validation_source_units': validation}
            manifest['heads'].append(record)
            write_json(output / ('FIT_SEED' + str(seed) + '.json'), record)
            print('FITTED', seed, fit['best_epoch'], validation['median_pearson'], flush=True)
        manifest['peak_rss_bytes'] = memory.peak
    manifest['elapsed_seconds'] = time.monotonic() - started
    write_json(folder / 'MANIFEST.json', manifest)
    write_json(output / 'TRAINING_MANIFEST.json', manifest)
    load_heads(folder)  # Safe round-trip of actual selected states.
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('train'); p.add_argument('--config', type=Path, required=True)
    p.add_argument('--checkpoints', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('predict'); p.add_argument('--checkpoints', type=Path, required=True)
    p.add_argument('--features', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); root = args.root.resolve()
    if args.action == 'train':
        train(root, args.config.resolve(), args.checkpoints.resolve(), args.output.resolve())
    else:
        output = checked_namespace(root, args.output, 'reproduction/runs')
        if output.exists(): raise FileExistsError(output)
        with np.load(args.features, allow_pickle=False) as f: features = f['embeddings']
        torch.set_num_threads(4); torch.use_deterministic_algorithms(True)
        with threadpool_limits(4): image, _ = image_predictions(args.checkpoints, features)
        output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(output, prediction=image, origin=np.array(ORIGIN))

if __name__ == '__main__': main()
