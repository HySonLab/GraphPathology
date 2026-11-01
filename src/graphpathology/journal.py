"""Explicit post-v1 journal research commands, independent of immutable publication tools."""
import argparse
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from .assets import sha
from .core import gene_metrics, summarize
from .journal_acquire import PROTOCOL_COMMIT
from .journal_data import load_patient, original_masks
from .journal_methods import predict_anchors


def provenance(root, **extra):
    return {'protocol_commit': PROTOCOL_COMMIT, 'code_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
            'code_hashes': {str(p.relative_to(root)): sha(p) for p in sorted((root / 'src/graphpathology').glob('journal*.py'))},
            'device': 'cpu', 'environment': json.loads((root / 'experiments/journal_revision/ENVIRONMENT.json').read_text()), **extra}


def replay(root, data_root, cohort, patient, output):
    started = time.monotonic(); data = load_patient(root, data_root, cohort, patient)
    records = [r for r in original_masks(root, cohort, patient) if r['regime'] == 'V1' and r['seed'] == 0]
    if len(records) != 1: raise ValueError('Expected one seed0 V1 trajectory')
    record = records[0]
    with np.load(record['absolute_mask_path']) as f: observed = f['order'][:int(np.floor(.2 * len(data['xy']) + .5))]
    hidden = np.setdiff1d(np.arange(len(data['xy'])), observed)
    prediction, _ = predict_anchors(data['xy'], data['targets'], observed, hidden, {'kind': 'idw', 'k': 16, 'p': 2})
    metrics = gene_metrics(data['targets'][hidden], prediction); summary = summarize(metrics)
    cp = {'wu': 'C9', 'bassiouni': 'C10'}[cohort]; old = root / 'checkpoints/legacy' / cp / (cp + 'B')
    logs = []
    for p in sorted((old / 'logs/runs').glob('*.json')):
        r = json.loads(p.read_text())
        if r['mask_id'] == record['mask_id'] and r['budget'] == .2 and not r['control']: logs.append(r)
    if len(logs) != 1: raise ValueError('Ambiguous original run')
    original = logs[0]; deltas = {}
    for key, expected in original['summary']['M1'].items():
        actual = summary.get(key, summary.get('n_defined_pearson') if key == 'n_valid_pearson' else None)
        if actual is None: raise ValueError('Unrecognized original summary key ' + key)
        deltas[key] = float(abs(actual - expected))
        if 'spearman' not in key: np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10, equal_nan=True)
    csv = old / (cp + 'B_GENE_METRICS.csv' + ('.gz' if cohort == 'bassiouni' else ''))
    gene_rows = pd.read_csv(csv); gene_rows = gene_rows[(gene_rows.run_id == original['run_id']) & (gene_rows.method == 'M1')].set_index('gene').loc[data['genes']]
    max_gene_delta = {}
    for key in ['pearson', 'spearman', 'mse', 'mae']:
        expected = gene_rows[key].to_numpy(); actual = metrics[key]
        if key != 'spearman': np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10, equal_nan=True)
        max_gene_delta[key] = float(np.nanmax(abs(actual - expected)))
    corrupt = data['targets'].copy(); corrupt[hidden] = np.nan
    check, _ = predict_anchors(data['xy'], corrupt, observed, hidden, {'kind': 'idw', 'k': 16, 'p': 2})
    np.testing.assert_array_equal(check, prediction)
    from scipy.stats import rankdata
    from .core import column_corr
    rounded_spearman = column_corr(rankdata(np.round(data['targets'][hidden], 12), axis=0), rankdata(np.round(prediction, 12), axis=0))
    tie_sensitivity = float(np.nanmax(abs(metrics['spearman'] - rounded_spearman)))
    result = provenance(root, status='PASS_PRIMARY_WITH_SECONDARY_RANK_SENSITIVITY', cohort=cohort, patient=patient, original_run_id=original['run_id'], mask_sha256=record['mask_sha256'],
                        source_hashes=data['source_hashes'], derived_coordinate_max_abs_delta=data['derived_coordinate_max_abs_delta'],
                        spearman_status='NOT_NUMERICALLY_REPLAYED: near-tied floating-point ranks; no historical Spearman exactness claim',
                        spearman_round12_sensitivity_max= tie_sensitivity, primary_metrics='Pearson, MSE and MAE pass rtol=atol=1e-10 per gene', normalization_sha256=data['normalization_sha256'],
                        max_summary_abs_deltas=deltas, max_gene_abs_deltas=max_gene_delta, n_locations=len(data['xy']), genes=data['genes'].tolist(),
                        hidden_corruption_invariant=True, elapsed_seconds=time.monotonic() - started)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists(): raise FileExistsError(output)
    output.write_text(json.dumps(result, indent=2) + '\n'); print(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--data', type=Path, default=Path('data/journal_revision'))
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('acquire'); p.add_argument('--cohorts', nargs='+', choices=['wu', 'bassiouni', 'her2st', 'hest', 'resnet18'], required=True); p.add_argument('--dry-run', action='store_true')
    p = sub.add_parser('original-uncertainty'); p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('replay'); p.add_argument('--cohort', choices=['wu', 'bassiouni'], required=True); p.add_argument('--patient', required=True); p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('anchors'); p.add_argument('--cohort', choices=['wu', 'bassiouni'], required=True); p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); root = args.root.resolve(); data = args.data.resolve()
    if args.action == 'acquire':
        from .journal_acquire import acquire
        print(json.dumps(acquire(root, data, args.cohorts, args.dry_run), indent=2))
    elif args.action == 'original-uncertainty':
        from .journal_statistics import original_uncertainty
        result = original_uncertainty(root, args.output)
        (args.output / 'PROVENANCE.json').write_text(json.dumps(provenance(root, **result), indent=2) + '\n'); print(json.dumps(result, indent=2))
    elif args.action == 'anchors':
        from .journal_anchors import run
        print(json.dumps(run(root, data, args.cohort, args.output), indent=2))
    elif args.action == 'replay': replay(root, data, args.cohort, args.patient, args.output)

if __name__ == '__main__': main()
