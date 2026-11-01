"""Patient-level paired effect tables; original-export and new-prediction analyses stay distinct."""
import numpy as np
import pandas as pd
from .journal_methods import patient_bootstrap

METRICS = ['median_pearson', 'mean_mse', 'mean_mae']


def contrast(candidate, reference, keys, name, origin):
    if candidate.duplicated(keys).any() or reference.duplicated(keys).any(): raise ValueError('Patient pairing is ambiguous')
    paired = candidate[keys + METRICS].merge(reference[keys + METRICS], on=keys, suffixes=('_candidate', '_reference'), validate='one_to_one', how='outer', indicator=True)
    if not (paired['_merge'] == 'both').all(): raise ValueError('Unmatched patient effects; no silent exclusion')
    paired = paired.drop(columns='_merge'); paired['contrast'] = name; paired['origin'] = origin
    paired['delta_pearson'] = paired.median_pearson_candidate - paired.median_pearson_reference
    paired['delta_mae'] = paired.mean_mae_candidate - paired.mean_mae_reference
    paired['relative_mse'] = (paired.mean_mse_candidate - paired.mean_mse_reference) / paired.mean_mse_reference
    paired.loc[paired.mean_mse_reference <= 0, 'relative_mse'] = np.nan
    return paired


def summarize_paired(paired, grouping):
    output = []
    for key, group in paired.groupby(grouping, dropna=False, sort=True):
        attributes = dict(zip(grouping, key if isinstance(key, tuple) else [key]))
        if group.patient.duplicated().any(): raise ValueError('Resampling unit must be one row per patient')
        for metric in ['delta_pearson', 'delta_mae', 'relative_mse']:
            values = group[metric].to_numpy(); finite = values[np.isfinite(values)]
            if attributes['cohort'] == 'Bassiouni':
                result = patient_bootstrap(values, direction='positive' if metric == 'delta_pearson' else 'negative')
            else:
                result = {'n_patients': len(finite), 'median_effect': float(np.median(finite)) if len(finite) else float('nan'),
                          'ci_low': float('nan'), 'ci_high': float('nan'), 'bootstrap_resamples': 0,
                          'fraction_improved': float(np.mean(finite > 0 if metric == 'delta_pearson' else finite < 0)) if len(finite) else float('nan')}
            output.append({**attributes, 'effect': metric, **result, 'uncertainty_unit': 'patient',
                           'interpretation': 'post-v1 descriptive; Wu individual patient effects are primary'})
    return pd.DataFrame(output)


def original_uncertainty(root, output):
    source = root / 'manuscript/revised/source_data/SuppTable2_per_patient.csv'
    data = pd.read_csv(source)
    data = data[data.cohort.isin(['Wu', 'Bassiouni'])].copy()
    paired = []
    for candidate, reference, name in [('M1', 'M0', 'anchors_vs_image'), ('M3', 'M1', 'fusion_vs_anchors')]:
        paired.append(contrast(data[data.method == candidate], data[data.method == reference], ['cohort', 'patient', 'coverage', 'regime'], name, 'original_v1_export_reanalysis'))
    low = data[(data.method == 'M1') & (data.coverage == .1) & (data.regime == 'V1')].copy()
    high = data[(data.method == 'M1') & (data.coverage == .4) & (data.regime == 'V1')].copy()
    c = contrast(high, low, ['cohort', 'patient'], 'coverage_40_vs_10_V1', 'original_v1_export_reanalysis'); c['coverage'] = .4; c['regime'] = 'V1'; paired.append(c)
    v1 = data[(data.method == 'M1') & (data.coverage == .2) & (data.regime == 'V1')]
    v2 = data[(data.method == 'M1') & (data.coverage == .2) & (data.regime == 'V2')]
    c = contrast(v2, v1, ['cohort', 'patient'], 'historical_V2_vs_V1_compound_protocol', 'original_v1_export_reanalysis'); c['coverage'] = .2; c['regime'] = 'V2_vs_V1'; paired.append(c)
    pairs = pd.concat(paired, ignore_index=True)
    summary = summarize_paired(pairs, ['origin', 'cohort', 'contrast', 'coverage', 'regime'])
    output.mkdir(parents=True, exist_ok=False)
    pairs.to_csv(output / 'ORIGINAL_V1_PAIRED_PATIENT_EFFECTS.csv', index=False)
    summary.to_csv(output / 'ORIGINAL_V1_PATIENT_BOOTSTRAP.csv', index=False)
    return {'source': str(source), 'paired_rows': len(pairs), 'summary_rows': len(summary), 'patient_counts': data.groupby('cohort').patient.nunique().to_dict()}
