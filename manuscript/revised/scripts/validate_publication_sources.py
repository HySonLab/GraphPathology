#!/usr/bin/env python3
"""Bounded publication preflight (stdlib only; no research experiments or redrawing)."""
from collections import Counter
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
errors = []

def check(condition, msg):
    if not condition:
        errors.append(msg)

manifest = ROOT / 'SCIENTIFIC_EXPORT_SHA256SUMS.txt'
checks = 0
for line in manifest.read_text().splitlines():
    digest, name = line.split(maxsplit=1)
    file = ROOT / name.strip()
    check(file.is_file(), f'Missing export: {name}')
    if file.is_file():
        check(hashlib.sha256(file.read_bytes()).hexdigest() == digest, f'Export changed: {name}')
    checks += 1

for i in range(1, 7):
    for ext in ('tex', 'pdf', 'svg'):
        check((ROOT / 'figures' / f'Figure{i}.{ext}').is_file(), f'Missing Figure{i}.{ext}')
for i in range(1, 13):
    check((ROOT / 'supplementary_figures' / f'FigureS{i}.pdf').is_file(), f'Missing FigureS{i}.pdf')

with (ROOT/'repository_documentation'/'REPRODUCTION_TARGETS.csv').open(newline='') as f:
    target_rows = list(csv.DictReader(f))
check(len(target_rows) == 29, f'Expected 29 publication targets, got {len(target_rows)}')

with (ROOT/'source_data'/'Figure2_synthetic_locations.csv').open(newline='') as f:
    locations = list(csv.DictReader(f))
with (ROOT/'source_data'/'Figure2_synthetic_edges.csv').open(newline='') as f:
    edges = list(csv.DictReader(f))
check(len(locations) == 36, 'Expected 36 synthetic locations')
counts = Counter(r['graph'] for r in edges)
expected = {
    'kNN (k=6, capped at 250 µm)': 94,
    'Radius (150 µm)': 22,
    'Radius (250 µm)': 98,
    'Delaunay (capped at 250 µm)': 85,
}
check(dict(counts) == expected, f'Unexpected synthetic edge counts: {dict(counts)}')
for row in edges:
    a, b = int(row['source_id']), int(row['target_id'])
    check(0 <= a < b < 36, f'Invalid undirected edge: {a}, {b}')
    cap = 150 if row['graph'] == 'Radius (150 µm)' else 250
    check(float(row['distance_um']) <= cap + 1e-6,
          f'Exceeded spatial cap for {row["graph"]}: {a}, {b}')

# Ensure none of the legacy plotters will overwrite the synthetic graph artwork.
files = {
    'redraw_development_figures.py': ['Figure3', 'Figure6'],
    'redraw_external_figures.py': ['Figure4', 'Figure5'],
}
for script, outputs in files.items():
    s = (ROOT/'scripts'/script).read_text()
    check(all(repr(x) in s for x in outputs), f'Wrong figure output mapping in {script}')
    check("'Figure2'" not in s and '"Figure2"' not in s,
          f'{script} may overwrite synthetic Figure 2')

if errors:
    raise SystemExit('Preflight FAILED:\n' + '\n'.join('- '+x for x in errors))
print(f'Preflight PASS: {checks} export hashes, 29 publication targets, '
      f'6 main figures, 12 supplementary figures; synthetic graph edges '
      f'{[counts[x] for x in expected]}; redraw mapping safe.')
