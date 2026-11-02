"""New, isolated replay outputs; all historical research artifacts remain read-only."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits

from .assets import sha
from .core import gene_metrics, summarize
from .journal_data import load_patient, original_masks, WU_EVAL
from .policy import count
from .reproduction_primary import checked_namespace, image_predictions, write_json

PRIMARY = ['median_pearson', 'mean_mse', 'mean_mae']
PANEL_METRICS = ['mean_pearson', 'median_pearson', 'mean_mse', 'median_mse', 'mean_mae', 'median_mae']


def patients(root):
    selected = json.loads((root / 'checkpoints/legacy/C10/C10A/logs/SOURCE_SAMPLE_SELECTION.json').read_text())['selected']
    return [('wu', p) for p in WU_EVAL] + [('bassiouni', p['patient_id']) for p in selected]


def data_root(root, cohort):
    return root / ('data/journal_revision' if cohort == 'wu' else 'data/journal_revision/transformer_egnn/geo_ranges')


def old_anchor_path(root, cohort, patient):
    return root / ('experiments/journal_revision/results/wu_anchors' if cohort == 'wu' else
                   'experiments/journal_revision/transformer_egnn/bassiouni_anchors') / patient


def compare_tables(actual, expected, keys, metrics, tolerance=1e-10):
    joined = actual[keys + metrics].merge(expected[keys + metrics], on=keys, how='outer',
                                         suffixes=('_rerun', '_reference'), validate='one_to_one', indicator=True)
    if not joined['_merge'].eq('both').all():
        raise ValueError('Replay reference support differs')
    deltas = {}; failures = []
    for metric in metrics:
        a, b = joined[metric + '_rerun'].to_numpy(), joined[metric + '_reference'].to_numpy()
        if not np.allclose(a, b, rtol=tolerance, atol=tolerance, equal_nan=True): failures.append(metric)
        difference = abs(a - b); finite = difference[np.isfinite(difference)]
        deltas[metric] = float(finite.max()) if len(finite) else None
    return {'status': 'NOT_REPRODUCED' if failures else 'EXACT_REPLAY', 'rows': len(joined),
            'failed_metrics': failures, 'max_abs_deltas': deltas, 'atol': tolerance, 'rtol': tolerance}


def tier_a(root, output):
    from .journal_anchors import run_patient, summarize_all
    output = checked_namespace(root, output, 'reproduction/runs')
    output.mkdir(parents=True, exist_ok=True)
    start = time.monotonic(); reports = []
    for cohort, patient in patients(root):
        folder = output / patient
        result = run_patient(root, data_root(root, cohort), cohort, patient, folder)
        expected = old_anchor_path(root, cohort, patient)
        record = {'patient': patient, 'cohort': cohort, 'elapsed_seconds': result['elapsed_seconds'],
                  'corruption_checks': result['corruption_checks'], 'source_hashes': result['input_source_hashes']}
        keys = ['experiment', 'mask_id', 'budget', 'regime', 'geometry', 'seed', 'condition', 'method']
        for scope, name, metrics in [('panel', 'PANEL_METRICS.csv.gz', PANEL_METRICS),
                                      ('gene', 'GENE_METRICS.csv.gz', ['pearson', 'mse', 'mae'])]:
            a, b = pd.read_csv(folder / name), pd.read_csv(expected / name)
            record[scope] = compare_tables(a, b, keys + (['gene'] if scope == 'gene' else []), metrics)
            record[scope]['reference_sha256'] = sha(expected / name)
        reports.append(record)
        write_json(output / 'REGRESSION.json', reports)
        print('STRICT REPLAY', patient, record['panel']['status'], record['gene']['status'], flush=True)
    summary = summarize_all(output)
    summary.update(elapsed_seconds=time.monotonic() - start,
                   status='EXACT_REPLAY' if all(r[s]['status']=='EXACT_REPLAY' for r in reports for s in ['panel','gene']) else 'PARTIAL',
                   primary_gene_values_checked=sum(r['gene']['rows'] * 3 for r in reports),
                   panel_rows_checked=sum(r['panel']['rows'] for r in reports),
                   corruption_checks=sum(r['corruption_checks'] for r in reports),
                   secondary_spearman='Near-tie sensitivity; no strict secondary-rank exactness claim.')
    write_json(output / 'SUMMARY.json', summary)
    return summary



def validate_external_alignment(cached, real):
    """External targets are float64; source HEST C0 targets have a different contract."""
    if cached['targets'].dtype != np.float64 or real['targets'].dtype != np.float64:
        raise ValueError('External targets must retain the historical float64 precision')
    for a,b in [('barcodes','barcodes'),('genes','genes'),('targets','targets')]:
        np.testing.assert_array_equal(cached[a],real[b])
    np.testing.assert_allclose(cached['xy_um'],real['xy'],rtol=0,atol=1e-8)


def tier_b(root, checkpoint, output):
    from .backbone_data import source
    from .external import predict as external_predict
    output = checked_namespace(root, output, 'reproduction/runs')
    if output.exists(): raise FileExistsError(output)
    output.mkdir(parents=True); start = time.monotonic()
    cal = json.loads((root / 'checkpoints/legacy/C9/C9A/T1_SCALER.json').read_text())
    panels = []; provenance = []; checks = 0
    for cohort, patient in patients(root):
        real = load_patient(root, data_root(root, cohort), cohort, patient)
        cached, meta = source(root, patient)
        validate_external_alignment(cached,real)
        ensemble, individual = image_predictions(checkpoint, cached['embeddings'])
        images = [('ensemble', ensemble)] + [('seed' + str(i), image) for i, image in enumerate(individual)]
        masks = [r for r in original_masks(root, cohort, patient) if r['regime'] == 'V1']
        if len(masks) != 10: raise ValueError('Expected exactly ten original random V1 masks')
        for mask in masks:
            with np.load(mask['absolute_mask_path']) as f: order = f['order']
            n = len(real['xy']); obs = order[:count(n, .2)]; hidden = np.setdiff1d(np.arange(n), obs)
            for training_seed, image in images:
                predictions, _ = external_predict(real['xy'], image, real['targets'], obs, hidden, 'V1', cal)
                for mode in ['nan', 'random']:
                    corrupt = real['targets'].copy()
                    corrupt[hidden] = np.nan if mode == 'nan' else np.random.default_rng(31).normal(size=(len(hidden),50))
                    alternate, _ = external_predict(real['xy'], image, corrupt, obs, hidden, 'V1', cal)
                    for method in ['M0', 'M3']:
                        np.testing.assert_array_equal(predictions[method], alternate[method]); checks += 1
                for method in ['M0', 'M3']:
                    panels.append({'cohort': cohort, 'patient': patient, 'method': method, 'mask': mask['mask_id'],
                                   'mask_sha256': mask['mask_sha256'], 'mask_seed': mask['seed'], 'training_seed': training_seed,
                                   'budget': .2, 'n_scored': len(hidden), **summarize(gene_metrics(real['targets'][hidden], predictions[method]))})
        provenance.append({'patient': patient, 'cohort': cohort, 'cache_sha256': meta['cache_sha256'],
                           'source_hashes': real['source_hashes'], 'normalization_sha256': real['normalization_sha256'],
                           'calibration_sha256': sha(root / 'checkpoints/legacy/C9/C9A/T1_SCALER.json'),
                           'calibration_applications_M0': 0, 'calibration_applications_M3': 1})
        print('PRIMARY EVALUATED', patient, flush=True)
    panel = pd.DataFrame(panels); panel.to_csv(output / 'PANEL.csv.gz', index=False)
    old = []
    for cohort, cp in [('wu', 'C9'), ('bassiouni', 'C10')]:
        d = pd.read_csv(root / f'checkpoints/legacy/{cp}/{cp}B/{cp}B_PANEL_METRICS.csv')
        d = d[(d.regime=='V1') & d.strategy.eq('random') & d.budget.eq(.2) & d.control.fillna('').eq('') & d.method.isin(['M0','M3'])].copy()
        d['patient'] = d['sample_id'] if 'sample_id' in d else d['patient_id']
        d['cohort'] = cohort
        d = d[d.patient.isin([p for _,p in patients(root)])]
        old.append(d)
    reference = pd.concat(old, ignore_index=True)
    keys = ['cohort','patient','method']
    ref_patient = reference.groupby(keys)[PRIMARY].mean().reset_index()
    actual_patient = panel.groupby(keys+['training_seed'])[PRIMARY].mean().reset_index()
    joined = actual_patient.merge(ref_patient, on=keys, suffixes=('_rerun','_paper'), validate='many_to_one')
    joined.to_csv(output / 'PATIENT_COMPARISON.csv', index=False)
    cohort_rows = []
    for (cohort, method, seed), group in joined.groupby(['cohort','method','training_seed']):
        for metric in PRIMARY:
            paper, rerun = float(group[metric+'_paper'].median()), float(group[metric+'_rerun'].median())
            delta = rerun-paper; relative = delta/abs(paper) if paper else None
            passed = abs(delta)<=.05 if metric=='median_pearson' else relative is not None and abs(relative)<=.15
            cohort_rows.append({'cohort':cohort,'method':method,'training_seed':seed,'metric':metric,'paper_value':paper,
                                'rerun_value':rerun,'delta':delta,'absolute_delta':abs(delta),'relative_delta':relative,
                                'n_patients':len(group),'aggregation':'mean over ten masks within patient; median over patients',
                                'status':'APPROXIMATE_RETRAIN' if passed else 'NOT_REPRODUCED','criterion_pass':passed})
    pd.DataFrame(cohort_rows).to_csv(output / 'COHORT_COMPARISON.csv', index=False)
    # Whole-section development replay never feeds selection, thresholds or calibration.
    development, _ = source(root, 'TENX99')
    ensemble, individual = image_predictions(checkpoint, development['embeddings'])
    dev = [{'training_seed': seed, **summarize(gene_metrics(development['targets'], image))}
           for seed, image in [('ensemble',ensemble)]+[('seed'+str(i),v) for i,v in enumerate(individual)]]
    write_json(output / 'DEVELOPMENT.json', dev)
    write_json(output / 'INPUT_PROVENANCE.json', provenance)
    summary = {'origin':'RETRAINED_NOT_ORIGINAL','patients':25,'training_seeds':3,'mask_seeds':10,
               'panel_rows':len(panel),'hidden_corruption_checks':checks,'elapsed_seconds':time.monotonic()-start,
               'ensemble_criteria_pass':all(r['criterion_pass'] for r in cohort_rows if r['training_seed']=='ensemble'),
               'scope':'Historical primary M0 and M3,20% distributed V1 only; other original neural branches are not replayed.'}
    write_json(output / 'SUMMARY.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    sub = parser.add_subparsers(dest='action',required=True)
    for action in ['tier-a','tier-b']:
        p=sub.add_parser(action);p.add_argument('--output',type=Path,required=True)
        if action=='tier-b':p.add_argument('--checkpoints',type=Path,required=True)
    args=parser.parse_args();root=args.root.resolve()
    torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
    with threadpool_limits(4):
        result=tier_a(root,args.output.resolve()) if args.action=='tier-a' else tier_b(root,args.checkpoints.resolve(),args.output.resolve())
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
