"""Publication regression against the frozen author package and independent numerical contracts.

This verifies exported results, not original-model predictions. No training/replay occurs.
"""
import csv
import json
from pathlib import Path
import re
import shutil
import numpy as np
import pandas as pd
from .assets import sha

RTOL = ATOL = 1e-12


def compare_trace(actual, expected, path='trace'):
    """Keep ordering, shapes, nulls and masks; compare floats without PDF metadata."""
    if isinstance(expected, dict):
        if actual.keys() != expected.keys(): raise ValueError('Trace keys differ: ' + path)
        for key in expected: compare_trace(actual[key], expected[key], path + '/' + key)
    elif isinstance(expected, list):
        if len(actual) != len(expected): raise ValueError('Trace length differs: ' + path)
        for i, (a, b) in enumerate(zip(actual, expected)): compare_trace(a, b, path + '/' + str(i))
    elif isinstance(expected, (float, int)) and not isinstance(expected, bool):
        if not np.isclose(actual, expected, rtol=RTOL, atol=ATOL, equal_nan=True):
            raise ValueError(f'Trace value differs: {path}: {actual} != {expected}')
    elif actual != expected:
        raise ValueError('Trace label/null/mask differs: ' + path)


def graph_checks(package):
    from .graph import build_graph
    nodes = pd.read_csv(package / 'source_data/Figure2_synthetic_locations.csv')
    edges = pd.read_csv(package / 'source_data/Figure2_synthetic_edges.csv')
    xy = nodes[['x_um', 'y_um']].to_numpy()
    if len(xy) != 36 or not np.array_equal(nodes.location_id, np.arange(36)):
        raise ValueError('Synthetic location order/count changed')
    settings = [({'kind': 'knn', 'k': 6, 'max_edge_um': 250}, 94),
                ({'kind': 'radius', 'radius_um': 150}, 22),
                ({'kind': 'radius', 'radius_um': 250}, 98),
                ({'kind': 'delaunay', 'max_edge_um': 250}, 85)]
    for (_, records), (config, count) in zip(edges.groupby('graph', sort=False), settings):
        a = build_graph(xy, config)
        i, j = a.nonzero(); actual = {(int(u), int(v)) for u, v in zip(i, j) if u < v}
        expected = set(zip(records.source_id, records.target_id))
        if len(expected) != count or actual != expected: raise ValueError('Synthetic topology differs')
        distance = np.linalg.norm(xy[records.source_id] - xy[records.target_id], axis=1)
        # Exported illustration coordinates and distances have six decimals; allow rounding only.
        np.testing.assert_allclose(distance, records.distance_um, rtol=0, atol=2e-6)
        for seed in [0, 1, 23]:
            order = np.random.default_rng(seed).permutation(36)
            b = build_graph(xy[order], config)
            u, v = b.nonzero()
            remapped = {tuple(sorted((int(order[x]), int(order[y])))) for x, y in zip(u, v) if x < y}
            if remapped != expected: raise ValueError('Graph depends on input node order')
    return {'nodes': 36, 'edges': [94, 22, 98, 85], 'permutations_checked': 3, 'distance_atol_um': 2e-6}


def numeric_rows(tex):
    return [[s.strip().replace('\\_', '_') for s in line.split('&')]
            for line in tex.splitlines() if '&' in line and re.search(r'\d', line)
            and '\\caption' not in line and '\\begin' not in line]


def number(s):
    return float(s.replace(',', '').replace('\\%', '').replace('\\', '').strip())


def table_checks(package):
    def read(name): return pd.read_csv(package / name, low_memory=False)
    checks = []
    t2 = read('tables/Table2_primary_20pct.csv')
    cells = [row for row in numeric_rows((package / 'tables/Table2.tex').read_text())
             if len(row) == 7 and row[3] in ['Image only', 'Spatial interpolation', 'Fixed fusion']]
    if len(cells) != len(t2): raise ValueError('Table2 row count changed')
    for row, (_, original) in zip(cells, t2.iterrows()):
        expected = [f'{original.median_pearson:.3f}', f'{original.mean_mse:.4f}', f'{original.mean_mae:.4f}']
        if [row[4], row[5], row[6].replace('\\', '').strip()] != expected: raise ValueError('Table2 printed value differs')
    checks.append('Table2: all 36 printed metric cells and 12-row order')
    # Independently reconstruct cohort medians from the patient summaries.
    patient = read('source_data/SuppTable2_per_patient.csv')
    for _, r in t2.iterrows():
        values = patient[(patient.cohort == r.cohort) & (patient.coverage == r.coverage)
                         & (patient.regime == r.regime) & (patient.method == r.method)]
        np.testing.assert_allclose(values[['median_pearson', 'mean_mse', 'mean_mae']].median(),
                                   r[['median_pearson', 'mean_mse', 'mean_mae']].to_numpy(float), rtol=RTOL, atol=ATOL)
    checks.append('Table2: patient-first cohort medians recomputed from exports')
    supp = (package / 'supplementary_content.tex').read_text()
    genes = read('methods_data/target_and_calibration_constants.csv').gene.tolist()
    ranked = {}
    for row in numeric_rows(supp):
        if len(row) == 4 and row[0].isdigit() and row[2].isdigit() and re.fullmatch(r'[A-Z0-9]+', row[1]):
            ranked[int(row[0])] = row[1]; ranked[int(row[2])] = row[3].replace('\\', '').strip()
    if [ranked[i] for i in range(1, 51)] != genes: raise ValueError('S1 gene order differs')
    checks.append('TableS1: all 50 gene identifiers/ranks in original order')
    fresh = patient[(patient.coverage == .2) & (patient.regime == 'V1')
                    & patient.cohort.isin(['Wu fresh validation', 'Bassiouni second replication'])]
    seen = []
    for row in numeric_rows(supp):
        if len(row) == 7 and re.fullmatch(r'(1160920F|CID\d+|BASS_P\d+)', row[0]):
            seen.append(row[0]); p = fresh[fresh.patient == row[0]].set_index('method').loc[['M0', 'M1', 'M3']]
            expected = [f'{x:.3f}' for x in [*p.median_pearson, *p.mean_mse]]
            if [s.replace('\\', '').strip() for s in row[1:]] != expected: raise ValueError('TableS2 patient metric differs: ' + row[0])
    expected_patients = ['1160920F', 'CID4290', 'CID44971'] + sorted(x for x in fresh.patient.unique() if x.startswith('BASS_'))
    if seen != expected_patients: raise ValueError('S2 patient ordering differs')
    checks.append('TableS2: 25 patients, 150 printed metric cells and exact ordering')
    coverage = read('tables/Table3_coverage_geometry.csv')
    cells = [row for row in numeric_rows((package / 'tables/Table3.tex').read_text()) if len(row) == 6 and '\\%' in row[1]]
    for row, (_, r) in zip(cells, coverage.iterrows()):
        expected = [f'{r.V1_r:.3f}', f'{r.V1_mse:.4f}', f'{r.V2_r:.3f}', f'{r.V2_mse:.4f}']
        if [s.replace('\\', '').strip() for s in row[2:]] != expected: raise ValueError('S9 printed metric differs')
        for regime in ['V1', 'V2']:
            p = patient[(patient.cohort == r.cohort) & (patient.coverage == r.coverage)
                        & (patient.regime == regime) & (patient.method == 'M1')]
            np.testing.assert_allclose(p[['median_pearson', 'mean_mse']].median(), [r[regime + '_r'], r[regime + '_mse']], rtol=RTOL, atol=ATOL)
    if len(cells) != 6: raise ValueError('S9 coverage row count differs')
    inline = [row for row in numeric_rows(supp) if len(row) == 6 and '\\%' in row[1]]
    if inline != cells: raise ValueError('Printed inline S9 differs from its editable Table3 source')
    checks.append('TableS9: all coverage cells and within-cohort aggregation')
    cohorts = read('tables/Table1_cohorts.csv')
    cells = [row for row in numeric_rows((package / 'tables/Table1.tex').read_text()) if len(row) == 6 and '/' in row[2]]
    if len(cells) != 4: raise ValueError('Table1 cohort row count differs')
    for row, (_, r) in zip(cells, cohorts.iterrows()):
        for value, expected in [(row[2], [r.patient_count, r.evaluation_patient_count]),
                                (row[3], [r.location_count, r.evaluation_locations])]:
            if [number(x) for x in value.split('/')] != expected: raise ValueError('Table1 counts differ')
    checks.append('Table1: all patient/location count cells')
    audit = read('source_data/SuppTable5_integrity.csv')
    fields = [('Mask trajectories', 'mask_trajectories', 'label_free_mask_replays'), ('Clean run bundles', 'clean_runs', 'clean_runs'),
              ('Control run bundles', 'control_runs', 'control_runs'), ('Method prediction replays', 'replayed_method_predictions', 'replayed_method_predictions'),
              ('Metric array replays', 'replayed_metrics', 'replayed_metric_arrays'), ('Hidden-label corruption checks', 'hidden_corruption_checks', 'hidden_corruption_checks')]
    for label, wu, bass in fields:
        row = next(row for row in numeric_rows(supp) if row[0] == label)
        if [number(x) for x in row[1:]] != [audit.iloc[0][wu], audit.iloc[1][bass]]: raise ValueError('S5 audit cells differ')
    checks.append('TableS5: all 12 recorded audit cells; historical replay not rerun')
    trends = read('source_data/Bassiouni_distance_trends.csv')
    for label, field in [('Higher MSE in farthest than nearest occupied bin', 'farthest_minus_nearest_mse'),
                         ('Strictly increasing occupied-bin MSE', 'strictly_increasing_mse')]:
        row = next(row for row in numeric_rows(supp) if row[0] == label)
        expected = [int((trends[trends.regime == reg][field] > 0).sum()) for reg in ['V1', 'V2']]
        actual = [int(x.split('/')[0]) for x in row[1:]]
        if actual != expected: raise ValueError('S7 support counts differ')
    checks.append('TableS7: distance support counts from per-patient exported trends')
    from .graph import parameter_count
    block = supp[supp.index('Selected supporting predictor architectures'):supp.index('\\label{tab:coverage}')]
    for d, h, g in re.findall(r'(\d+)\\rightarrow(\d+)\\rightarrow(\d+)', block):
        count = parameter_count(int(d), int(h), int(g))
        if count <= 0: raise ValueError('Invalid S8 dimensions')
    for row in numeric_rows(block):
        dimensions = re.findall(r'(\d+)\\rightarrow(\d+)\\rightarrow(\d+)', ' '.join(row))
        printed = re.search(r'([\d,]+) head parameters', ' '.join(row))
        if len(dimensions) == 1 and printed:
            d, h, g = map(int, dimensions[0])
            if parameter_count(d, h, g) != int(printed.group(1).replace(',', '')): raise ValueError('S8 printed parameter count differs from production architecture')
    for count in [25650, 72114, 72075, 144178, 209714, 1750]:
        if f'{count:,}' not in block: raise ValueError('S8 architecture constant missing')
    checks.append('TableS8: architecture dimensions and printed parameter-count schema')
    for title, tokens in [('Per-gene results and metric definitions', ['Constant targets', 'average', 'Panel error']),
                           ('Negative-control families', ['Anchor-expression', 'Anchor-position', 'Image-location']),
                           ('Observation masks and sampling protocols', ['0.10,0.20,0.40', '0.5', '40\\%'])]:
        index = supp.index(title); block = supp[index:supp.index('\\end{table}', index)]
        if not all(t in block for t in tokens): raise ValueError('Definition table schema changed: ' + title)
    checks.append('TablesS3/S4/S6: definitions/masking schemas; complete digital exports hash checked')
    return checks


def constants_checks(root, package):
    old = root / 'checkpoints/legacy'
    constants = pd.read_csv(package / 'methods_data/target_and_calibration_constants.csv')
    normalization = json.loads((old / 'reports/C0/C0_NORMALIZATION.json').read_text())
    genes = constants.gene.tolist(); index = [normalization['genes'].index(g) for g in genes]
    np.testing.assert_allclose(constants.log_count_center, np.array(normalization['mean'])[index], rtol=RTOL, atol=ATOL)
    np.testing.assert_allclose(constants.log_count_scale, np.array(normalization['scale'])[index], rtol=RTOL, atol=ATOL)
    calibration = json.loads((old / 'C9/C9A/T1_SCALER.json').read_text())
    if genes != calibration['genes']: raise ValueError('Calibration gene order changed')
    for domain, prefix in [('source', 'source'), ('destination', 'Wu')]:
        for key, suffix in [('center', 'center'), ('scale', 'scale')]:
            np.testing.assert_allclose(constants[prefix + '_calibration_' + suffix], calibration[domain][key], rtol=RTOL, atol=ATOL)
    folder = old / 'C4/C4A/preprocessing/resnet18_reference/primary/S2-50'
    with np.load(folder / 'target.npz') as f:
        for key in ['center', 'scale']:
            np.testing.assert_allclose(constants['head_target_' + key], f['mean' if key == 'center' else key], rtol=RTOL, atol=ATOL)
    features = pd.read_csv(package / 'methods_data/image_feature_scaler.csv')
    with np.load(folder / 'feature.npz') as f:
        np.testing.assert_allclose(features.iloc[:, 1], f['mean'], rtol=RTOL, atol=ATOL)
        np.testing.assert_allclose(features.iloc[:, 2], f['scale'], rtol=RTOL, atol=ATOL)
    from .assets import inventory
    identities = {r['original_path']: r['sha256'] for r in inventory(root)}
    for row in csv.DictReader((package / 'repository_documentation/MODEL_FILES_TO_RELEASE.csv').open()):
        if identities.get(row['expected_repository_path']) != row['sha256']:
            raise ValueError('Primary fitted-head identity differs from original registry')
    return ['50-gene normalization/head-target/source-Wu calibration, 512-feature scaler and three original head identities agree with archived originals']


def independent_figures(package, traces):
    checks = []
    for figure, name in [('Figure4', 'Fig5_source.csv'), ('Figure5', 'Fig6_source.csv')]:
        d = pd.read_csv(package / 'source_data' / name)
        t = traces[figure]
        for ax, metric in [(1, 'median_pearson'), (2, 'mean_mse')]:
            for line, method in zip(t[ax]['lines'], ['M0', 'M1', 'M3']):
                p = d[(d.regime == 'V1') & (d.method == method)].groupby('coverage')[metric].median()
                np.testing.assert_allclose(line['y']['values'], p.to_numpy(), rtol=RTOL, atol=ATOL)
        p = d[(d.regime == 'V1') & (d.coverage == .2)].pivot(index='patient', columns='method', values='median_pearson')
        for col, first, second in zip(t[3]['collections'], ['M1', 'M3'], ['M0', 'M1']):
            np.testing.assert_allclose(np.array(col['offsets']['values'])[:, 1], p[first] - p[second], rtol=RTOL, atol=ATOL)
        checks.append(figure + ': independent patient median curves and within-patient paired effects')
    maps = pd.read_csv(package / 'source_data/Fig3_map_source.csv')
    cap = maps.absolute_error.quantile(.98)
    for ax, regime in zip(traces['Figure6'][2:4], ['random', 'block']):
        selected = maps[(maps.regime == regime) & ~maps.observed]
        col = ax['collections'][0]
        np.testing.assert_allclose(col['clim'], [0, cap], rtol=RTOL, atol=ATOL)
        np.testing.assert_allclose(col['offsets']['values'], selected[['x_um', 'y_um']], rtol=RTOL, atol=ATOL)
        np.testing.assert_allclose(col['array']['values'], selected.absolute_error, rtol=RTOL, atol=ATOL)
    if not (maps.regime == 'capture_window').any(): raise ValueError('Full capture-window map values missing')
    checks.append('Figure6: full-export 98th percentile cap, exact point ordering/values including preserved capture-window inputs')
    sampling = pd.read_csv(package / 'source_data/SFig6_source.csv')
    if len(sampling) != 4448: raise ValueError('Incomplete sampling export')
    checks.append('4448-row sampling export retained; all source rows/nulls hash checked')
    return checks, float(cap)


def validate(root, generated, reference, render_reference):
    from .publication import copy_package, render
    checks = []
    # The reference is the unmodified author package whose complete bytes are pinned by ZIP provenance.
    provenance = json.loads((root / 'manuscript/PACKAGE_PROVENANCE.json').read_text())
    for name, expected in provenance['files'].items():
        if sha(reference / name) != expected: raise ValueError('Reference package bytes differ: ' + name)
    checks.append(f'Unmodified author baseline: {len(provenance["files"])} file hashes')
    for line in (reference / 'SCIENTIFIC_EXPORT_SHA256SUMS.txt').read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        if sha(generated / name.strip()) != expected: raise ValueError('Scientific source bytes differ: ' + name)
        checks.append('Scientific export bytes/schema/row order/null pattern: ' + name.strip())
    baseline = copy_package(reference, render_reference)
    render(baseline)
    actual = json.loads((generated / 'NUMERICAL_TRACES.json').read_text())
    expected = json.loads((baseline / 'NUMERICAL_TRACES.json').read_text())
    # Every main/SI figure is traced, including Figure1 atomic in-memory exports.
    names = {'Figure' + str(n) for n in range(1, 7)} | {'FigureS' + str(n) for n in range(1, 13)}
    if not names <= actual.keys() or actual.keys() != expected.keys(): raise ValueError('Incomplete plotted-value traces')
    for figure in sorted(names):
        compare_trace(actual[figure], expected[figure], figure)
        checks.append(figure + ': all ordered plotted arrays/masks, scalar colors, color limits and labels agree with author renderer')
    if sha(generated / 'scripts/draw_architecture.py') != sha(reference / 'scripts/draw_architecture.py'):
        raise ValueError('Figure1 schematic renderer changed')
    checks.append('Figure1: exact author schematic renderer and synthetic coordinates preserved; no patient measurements')
    checks += table_checks(generated)
    checks += constants_checks(root, generated)
    graph = graph_checks(generated); checks.append('Figure2: production graph equivalence, caps/symmetry, rounded distances and permutation invariance')
    independent, cap = independent_figures(generated, actual); checks += independent
    return {'status': 'PASS', 'targets': 29, 'checks': checks, 'rtol': RTOL, 'atol': ATOL,
            'graph': graph, 'figure6_color_cap': cap, 'pdf_byte_equality_required': False,
            'baseline': 'Unmodified author-supplied PREFLIGHT_FIXED package, pinned by ZIP/file SHA256',
            'original_prediction_replay': False,
            'scope': 'All author-rendered plotted arrays and all table sources checked; independent checks for printed numeric cells, definitions, ordering, graphs and constants. No new scientific results.'}
