"""Preregistered anchor-only estimators and paired patient uncertainty, journal V1."""
import hashlib
import numpy as np
from scipy.spatial import cKDTree
from .policy import region, count


def candidates():
    rows = [{'method': 'nearest', 'kind': 'uniform', 'k': 1}]
    rows += [{'method': f'uniform_k{k}', 'kind': 'uniform', 'k': k} for k in [4, 8, 16, 32]]
    rows += [{'method': f'idw_k{k}_p{p}', 'kind': 'idw', 'k': k, 'p': p} for k in [4, 8, 16, 32] for p in [1, 2, 3]]
    rows += [{'method': f'gaussian_k{k}', 'kind': 'gaussian', 'k': k} for k in [8, 16, 32]]
    return rows


def predict_anchors(xy, y, observed, queries, config):
    """Only observed y is read; queries may contain ANY unobserved locations."""
    xy = np.asarray(xy, dtype=float); observed = np.asarray(observed, dtype=int); queries = np.asarray(queries, dtype=int)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all(): raise ValueError('Finite Nx2 coordinates required')
    if len(observed) == 0 or len(set(observed)) != len(observed): raise ValueError('Nonempty unique anchors required')
    if np.intersect1d(observed, queries).size: raise ValueError('Scored queries overlap anchors')
    if min(observed.min(), queries.min(initial=0)) < 0 or max(observed.max(), queries.max(initial=0)) >= len(xy): raise ValueError('Index out of range')
    values = np.asarray(y, dtype=float)[observed].copy()
    if values.ndim != 2 or not np.isfinite(values).all() or len(y) != len(xy): raise ValueError('Finite observed targets and aligned shapes required')
    tree = cKDTree(xy[observed]); k = min(config['k'], len(observed))
    dist, idx = tree.query(xy[queries], k=k)
    if k == 1: dist, idx = dist[:, None], idx[:, None]
    dist = np.round(dist, 6)
    if config['kind'] == 'idw': weights = 1 / np.maximum(dist, 1e-6) ** config['p']
    elif config['kind'] == 'uniform': weights = np.ones_like(dist)
    elif config['kind'] == 'gaussian':
        if len(observed) == 1: bandwidth = 1e-6
        else:
            other_k = min(config['k'], len(observed) - 1)
            anchor_dist, _ = tree.query(xy[observed], k=other_k + 1)
            bandwidth = max(float(np.median(np.round(anchor_dist[:, -1], 6))), 1e-6)
        squared = dist ** 2
        weights = np.exp(-(squared - squared.min(1, keepdims=True)) / (2 * bandwidth ** 2))
    else: raise ValueError('Unknown anchor operator')
    weights /= weights.sum(1, keepdims=True)
    pred = np.einsum('nk,nkg->ng', weights, values[idx], optimize=True)
    return pred, dist[:, 0]


def identical_target_design(xy, patient, seed):
    """Freeze central R, fixed H, and two nested anchor trajectories using coordinates only."""
    # Existing region API returns a boolean forbidden vector.
    r = region(xy, 'V2', 1)
    r_indices = np.flatnonzero(r)
    holdout_seed = int.from_bytes(hashlib.sha256(('H|' + patient).encode()).digest()[:8], 'little')
    h = np.sort(np.random.default_rng(holdout_seed).permutation(r_indices)[:max(1, count(len(r_indices), .2))])
    non_h = np.setdiff1d(np.arange(len(xy)), h)
    outside = np.flatnonzero(~r)
    orders = {'distributed': np.random.default_rng(seed).permutation(non_h),
              'region_exclusion': np.random.default_rng(seed).permutation(outside)}
    return r, h, orders, holdout_seed


def paired_effects(candidate, reference):
    a = np.asarray(candidate, dtype=float); b = np.asarray(reference, dtype=float)
    return a - b


def patient_bootstrap(values, n_resamples=10000, seed=20261101, direction='positive'):
    """Median paired effect; patients are the resampling units, never genes/spots/seeds."""
    values = np.asarray(values, dtype=float); values = values[np.isfinite(values)]
    if not len(values): return {'n_patients': 0, 'median_effect': float('nan'), 'ci_low': float('nan'), 'ci_high': float('nan'), 'fraction_improved': float('nan')}
    rng = np.random.default_rng(seed)
    resampled = np.median(values[rng.integers(0, len(values), size=(n_resamples, len(values)))], axis=1)
    low, high = np.percentile(resampled, [2.5, 97.5])
    return {'n_patients': len(values), 'median_effect': float(np.median(values)), 'ci_low': float(low), 'ci_high': float(high),
            'fraction_improved': float(np.mean(values > 0 if direction == 'positive' else values < 0)),
            'bootstrap_resamples': n_resamples, 'bootstrap_seed': seed}
