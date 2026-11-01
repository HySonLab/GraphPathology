"""A/B1/B2 measured-data runner with original-mask pairing and all-unobserved corruption checks."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .assets import sha
from .core import gene_metrics, summarize, PeakRSS
from .journal_data import load_patient, original_masks, WU_EVAL
from .journal_methods import candidates, predict_anchors, identical_target_design
from .journal_statistics import contrast, summarize_paired, METRICS
from .policy import count

REFERENCE = 'idw_k16_p2'
OPERATORS = [r for r in candidates() if r['method'] in [REFERENCE, 'uniform_k16']]
BINS = [0, 150, 300, 600, 1200, 2400, np.inf]


def tensor_hash(x):
    x = np.asarray(x)
    return hashlib.sha256(str((x.shape, str(x.dtype))).encode() + x.tobytes(order='C')).hexdigest()


def run_patient(root, data_root, cohort, patient, output):
    from .journal import provenance, replay
    output.mkdir(parents=True, exist_ok=True)
    replay_path = output / 'REPLAY.json'
    if not replay_path.exists(): replay(root, data_root, cohort, patient, replay_path)
    replay_record = json.loads(replay_path.read_text())
    if not replay_record['status'].startswith('PASS_PRIMARY'): raise ValueError('Primary real-data replay must pass first')
    done = output / 'COMPLETE.json'
    if done.exists():
        r = json.loads(done.read_text())
        for name, digest in r['output_hashes'].items():
            if sha(output / name) != digest: raise ValueError('Completed output modified')
        return r
    started = time.monotonic(); data = load_patient(root, data_root, cohort, patient); y = data['targets']; xy = data['xy']; n = len(xy)
    panel_rows = []; gene_rows = []; bin_rows = []; design_records = []; corruption_checks = 0
    mask_cache = data_root / 'processed/journal_v2_masks' / cohort / patient; mask_cache.mkdir(parents=True, exist_ok=True)

    def evaluate(observed, hidden, config, attributes):
        nonlocal corruption_checks
        start = time.monotonic(); p, distances = predict_anchors(xy, y, observed, hidden, config)
        unknown = np.setdiff1d(np.arange(n), observed)
        altered = y.copy(); altered[unknown] = np.nan
        q, _ = predict_anchors(xy, altered, observed, hidden, config); np.testing.assert_array_equal(p, q)
        altered[unknown] = np.random.default_rng(17).normal(size=(len(unknown), y.shape[1]))
        q, _ = predict_anchors(xy, altered, observed, hidden, config); np.testing.assert_array_equal(p, q); corruption_checks += 2
        m = gene_metrics(y[hidden], p); row = {**attributes, 'method': config['method'], 'patient': patient,
             'cohort': 'Wu' if cohort == 'wu' else 'Bassiouni', 'n_observed': len(observed), 'n_scored': len(hidden),
             'n_unobserved_including_unscored': len(unknown), **summarize(m), 'mean_nearest_anchor_um': float(distances.mean()),
             'median_nearest_anchor_um': float(np.median(distances)), 'seconds_including_corruption': time.monotonic() - start}
        panel_rows.append(row)
        for j, gene in enumerate(data['genes']):
            gene_rows.append({**attributes, 'patient': patient, 'method': config['method'], 'gene': gene,
                              **{k: float(m[k][j]) for k in ['pearson', 'spearman', 'mse', 'mae']}})
        for low, high in zip(BINS[:-1], BINS[1:]):
            selected = (distances >= low) & (distances < high)
            if not selected.any(): continue
            err = y[hidden][selected] - p[selected]
            bin_rows.append({**attributes, 'patient': patient, 'method': config['method'], 'lower_um': low, 'upper_um': high,
                             'n_locations': int(selected.sum()), 'mean_mse': float(np.mean(err ** 2)), 'mean_mae': float(np.mean(abs(err))),
                             'mean_nearest_anchor_um': float(distances[selected].mean())})

    with PeakRSS() as memory:
        for r in original_masks(root, cohort, patient):
            with np.load(r['absolute_mask_path']) as f:
                order = f['order'].copy(); forbidden = f['forbidden'].copy()
                np.testing.assert_array_equal(data['barcodes'], f['barcodes']); np.testing.assert_array_equal(xy, f['xy_um'])
            design_records.append({k: r[k] for k in ['mask_id', 'mask_sha256', 'regime', 'geometry', 'seed']})
            for budget in [.1, .2, .4]:
                observed = order[:count(n, budget)]
                hidden = np.setdiff1d(np.arange(n), observed) if r['regime'] == 'V1' else np.flatnonzero(forbidden)
                base = {'mask_id': r['mask_id'], 'mask_sha256': r['mask_sha256'], 'budget': budget, 'regime': r['regime'], 'geometry': r['geometry'], 'seed': r['seed'], 'condition': 'historical'}
                for config in (candidates() if r['regime'] == 'V1' else OPERATORS):
                    evaluate(observed, hidden, config, {**base, 'experiment': 'A' if r['regime'] == 'V1' else 'B1'})
        for seed in range(10):
            region, hidden, orders, holdout_seed = identical_target_design(xy, patient, seed)
            central = next(r for r in original_masks(root,cohort,patient) if r['regime']=='V2' and r['geometry']==1 and r['seed']==0)
            with np.load(central['absolute_mask_path']) as f: np.testing.assert_array_equal(region,f['forbidden'])
            path = mask_cache / f'central_seed{seed}.npz'
            if not path.exists(): np.savez_compressed(path, region=region, hidden=hidden, distributed_order=orders['distributed'], region_exclusion_order=orders['region_exclusion'], barcodes=data['barcodes'], xy_um=xy)
            else:
                with np.load(path) as f:
                    for key, value in [('region', region), ('hidden', hidden), ('distributed_order', orders['distributed']), ('region_exclusion_order', orders['region_exclusion']), ('barcodes', data['barcodes']), ('xy_um', xy)]: np.testing.assert_array_equal(f[key], value)
            digest = sha(path); design_records.append({'experiment': 'B2', 'seed': seed, 'holdout_seed': holdout_seed, 'mask_sha256': digest,
                                                       'relative_path': str(path.relative_to(data_root)), 'n_region': int(region.sum()), 'n_holdout': len(hidden)})
            for budget in [.1, .2, .4]:
                k = count(n, budget)
                if any(len(order) < k for order in orders.values()): raise ValueError('Paired B2 capacity exclusion required; no asymmetric count reduction')
                for condition, order in orders.items():
                    for config in OPERATORS:
                        evaluate(order[:k], hidden, config, {'experiment': 'B2', 'mask_id': f'{patient}_central_seed{seed}', 'mask_sha256': digest,
                                 'budget': budget, 'regime': 'identical_H', 'geometry': 1, 'seed': seed, 'condition': condition})
    panel = pd.DataFrame(panel_rows); panel.to_csv(output / 'PANEL_METRICS.csv.gz', index=False)
    pd.DataFrame(gene_rows).to_csv(output / 'GENE_METRICS.csv.gz', index=False)
    pd.DataFrame(bin_rows).to_csv(output / 'DISTANCE_BINS.csv.gz', index=False)
    (output / 'MASK_MANIFEST.json').write_text(json.dumps(design_records, indent=2) + '\n')
    hashes = {name: sha(output / name) for name in ['PANEL_METRICS.csv.gz', 'GENE_METRICS.csv.gz', 'DISTANCE_BINS.csv.gz', 'MASK_MANIFEST.json', 'REPLAY.json']}
    result = provenance(root, patient=patient, cohort=cohort, status='COMPLETE', n_run_rows=len(panel), corruption_checks=corruption_checks,
                        n_locations=n, n_genes=y.shape[1], input_source_hashes=data['source_hashes'], normalization_sha256=data['normalization_sha256'],
                        ordered_barcodes_sha256=tensor_hash(data['barcodes']), ordered_genes_sha256=tensor_hash(data['genes']),
                        coordinate_array_sha256=tensor_hash(xy), target_array_sha256=tensor_hash(y), output_hashes=hashes,
                        elapsed_seconds=time.monotonic() - started, peak_rss_mb=memory.peak / 2**20,
                        secondary_spearman='Near-tied rank sensitivity; not used for historical exactness or primary claims')
    done.write_text(json.dumps(result, indent=2) + '\n'); print('PATIENT COMPLETE', patient, len(panel), flush=True)
    return result


def summarize_all(output):
    panels = [pd.read_csv(p) for p in sorted(output.glob('*/PANEL_METRICS.csv.gz')) if (p.parent / 'COMPLETE.json').exists()]
    if not panels: raise ValueError('No completed patients')
    runs = pd.concat(panels, ignore_index=True)
    keys = ['cohort', 'patient', 'experiment', 'budget', 'regime', 'condition', 'geometry', 'method']
    # Equal seed mean within geometry, then equal geometry mean within patient.
    by_geometry = runs.groupby(keys, dropna=False)[METRICS + ['mean_nearest_anchor_um', 'median_nearest_anchor_um']].mean().reset_index()
    patient = by_geometry.groupby([k for k in keys if k != 'geometry'], dropna=False)[METRICS + ['mean_nearest_anchor_um', 'median_nearest_anchor_um']].mean().reset_index()
    pairs = []
    for method in [r['method'] for r in candidates() if r['method'] != REFERENCE]:
        a = patient[(patient.experiment == 'A') & (patient.method == method)]
        ref = patient[(patient.experiment == 'A') & (patient.method == REFERENCE)]
        p = contrast(a, ref, ['cohort', 'patient', 'budget', 'regime'], method + '_vs_' + REFERENCE, 'new_A_matched_information'); pairs.append(p)
    # V1 operator rows already appear in A; add them to B1 without rerunning or double counting A summaries.
    b1 = pd.concat([patient[patient.experiment == 'B1'], patient[(patient.experiment == 'A') & patient.method.isin([REFERENCE, 'uniform_k16'])]], ignore_index=True)
    p = contrast(b1[b1.method == 'uniform_k16'], b1[b1.method == REFERENCE], ['cohort', 'patient', 'budget', 'regime'], 'uniform_k16_vs_idw_k16_p2', 'new_B1_operator_factorial'); pairs.append(p)
    for method in [REFERENCE, 'uniform_k16']:
        b2 = patient[(patient.experiment == 'B2') & (patient.method == method)]
        p = contrast(b2[b2.condition == 'region_exclusion'], b2[b2.condition == 'distributed'], ['cohort', 'patient', 'budget', 'regime'], method + '_region_exclusion_vs_distributed', 'new_B2_identical_targets'); pairs.append(p)
    paired = pd.concat(pairs, ignore_index=True)
    summary = summarize_paired(paired, ['origin', 'cohort', 'contrast', 'budget', 'regime'])
    patient.to_csv(output / 'PATIENT_METRICS.csv', index=False); paired.to_csv(output / 'PAIRED_PATIENT_EFFECTS.csv', index=False); summary.to_csv(output / 'PATIENT_UNCERTAINTY.csv', index=False)
    by_geometry.to_csv(output / 'GEOMETRY_METRICS.csv', index=False)
    distance = patient[patient.experiment == 'B2'][['cohort','patient','budget','method','condition','mean_nearest_anchor_um','median_nearest_anchor_um']]
    distance.to_csv(output / 'B2_PATIENT_ANCHOR_DISTANCE.csv', index=False)
    valid = runs.groupby(['cohort','patient','experiment','budget','regime','condition','method']).agg(min_defined_pearson=('n_defined_pearson','min'),mean_defined_pearson=('n_defined_pearson','mean'),max_defined_pearson=('n_defined_pearson','max'),n_panel_genes=('n_genes','first')).reset_index()
    valid.to_csv(output / 'VALID_CORRELATION_COUNTS.csv',index=False)
    return {'completed_patients': patient.groupby('cohort').patient.nunique().to_dict(), 'run_rows': len(runs), 'paired_rows': len(paired), 'summary_rows': len(summary)}


def run(root, data_root, cohort, output):
    output.mkdir(parents=True, exist_ok=True)
    if cohort == 'wu': patients = WU_EVAL
    else: patients = [p['patient_id'] for p in json.loads((root / 'checkpoints/legacy/C10/C10A/logs/SOURCE_SAMPLE_SELECTION.json').read_text())['selected']]
    for patient in patients: run_patient(root, data_root, cohort, patient, output / patient)
    result = summarize_all(output); (output / 'SUMMARY.json').write_text(json.dumps(result, indent=2) + '\n'); return result
