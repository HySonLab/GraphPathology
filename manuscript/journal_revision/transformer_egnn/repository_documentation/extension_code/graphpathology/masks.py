"""Coordinate-only observation masks, where True means measured transcriptomics."""
import hashlib
import numpy as np
from scipy.spatial import cKDTree

REGIMES = ('random', 'block', 'capture_window')

def stable_seed(sample_id, regime, fraction, replicate, seed=20260927):
    key = f'{seed}|{sample_id}|{regime}|{fraction:.8f}|{replicate}'.encode()
    return int.from_bytes(hashlib.sha256(key).digest()[:8], 'little')


def make_mask(xy_um, regime, fraction, seed):
    xy = np.asarray(xy_um, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or len(xy) < 2 or not np.isfinite(xy).all():
        raise ValueError('Finite Nx2 coordinates with at least two locations required')
    if regime not in REGIMES or not 0 < fraction < 1:
        raise ValueError('Unknown regime or invalid coverage')
    rng = np.random.default_rng(seed)
    n = len(xy)
    k = min(n - 1, max(1, int(round(fraction * n))))
    meta = {'requested_spot_coverage': fraction, 'regime': regime, 'seed': seed}
    mask = np.zeros(n, dtype=bool)
    if regime == 'random':
        mask[rng.choice(n, k, replace=False)] = True
        meta['geometry'] = 'uniform random locations without replacement'
    else:
        center = xy[rng.integers(n)]
        distance = np.max(np.abs(xy - center), axis=1) if regime == 'block' else np.linalg.norm(xy-center, axis=1)
        count = n-k if regime == 'block' else k
        radius = float(np.partition(distance, count-1)[count-1])
        inside = distance <= radius  # Keep all boundary ties; report actual coverage.
        mask = ~inside if regime == 'block' else inside
        meta.update(center_um=center.tolist(), radius_um=radius,
                    geometry='hidden axis-aligned square' if regime == 'block' else 'observed circular window')
        meta['nominal_window_area_um2'] = (4 if regime == 'block' else np.pi) * radius**2
        meta['tissue_area_coverage'] = None
        meta['area_note'] = 'No tissue polygon supplied; spot fraction is not tissue-area fraction. Simulated, not platform-specific.'
    if not mask.any() or mask.all():
        raise ValueError('Geometry produced an empty observed or hidden set')
    meta['actual_spot_coverage'] = float(mask.mean())
    meta['observed_count'] = int(mask.sum())
    meta['hidden_count'] = int((~mask).sum())
    return mask, meta


def anchor_distances(xy_um, mask):
    mask = np.asarray(mask)
    if mask.dtype != bool or mask.shape != (len(xy_um),):
        raise ValueError('Invalid mask')
    if not mask.any():
        return np.full(len(mask), np.inf)
    return cKDTree(np.asarray(xy_um)[mask]).query(xy_um)[0]
