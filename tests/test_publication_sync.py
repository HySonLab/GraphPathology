"""Offline publication contracts; full renderer/TeX integration is a separate bounded command."""
from pathlib import Path
import copy
import json
import subprocess
import sys
import numpy as np
import pytest
from graphpathology.publication import check_sources, arxiv_dependencies
from graphpathology.publication_regression import compare_trace, graph_checks, table_checks, constants_checks

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'manuscript/revised'


def test_author_package_and_all_target_paths():
    package, rows = check_sources(ROOT)
    assert package == PACKAGE
    assert len(rows) == 29
    assert {r['target'] for r in rows} == {'Figure' + str(i) for i in range(1, 7)} | {'FigureS' + str(i) for i in range(1, 13)} | {'Table1', 'Table2'} | {'TableS' + str(i) for i in range(1, 10)}


def test_synthetic_graph_equivalence_and_node_ordering():
    assert graph_checks(PACKAGE)['edges'] == [94, 22, 98, 85]


def test_printed_cells_ordering_and_definition_tables():
    assert len(table_checks(PACKAGE)) == 10


def test_method_constants_match_archived_preprocessing():
    assert len(constants_checks(ROOT, PACKAGE)) == 1


def test_arxiv_dependency_closure_has_one_root_and_no_numeric_dependencies():
    names = arxiv_dependencies(PACKAGE)
    assert 'main.tex' in names and 'main.bbl' in names
    assert 'article_only.tex' not in names and 'supplementary.tex' not in names
    assert sum(n.startswith('figures/') and n.endswith('.pdf') for n in names) == 6
    assert sum(n.startswith('supplementary_figures/') and n.endswith('.pdf') for n in names) == 12
    assert not any(n.startswith(('source_data/', 'methods_data/', 'anc/')) for n in names)
    assert not any(n in ['main.pdf', 'article_only.pdf', 'supplementary.pdf'] for n in names)


@pytest.mark.parametrize('mutation', ['value', 'mask', 'ordering', 'null', 'label'])
def test_regression_rejects_material_changes(mutation):
    expected = {'values': [1., float('nan'), 2.], 'mask': [False, True, False], 'label': 'FOXA1', 'scale': None}
    actual = copy.deepcopy(expected)
    if mutation == 'value': actual['values'][0] = 1.01
    elif mutation == 'mask': actual['mask'][1] = False
    elif mutation == 'ordering': actual['values'] = [2., float('nan'), 1.]
    elif mutation == 'null': actual['scale'] = 0
    elif mutation == 'label': actual['label'] = 'ERBB2'
    with pytest.raises(ValueError): compare_trace(actual, expected)


def test_regression_accepts_roundoff_and_nan_identity():
    compare_trace({'a': [1 + 1e-14, float('nan')]}, {'a': [1., float('nan')]})


@pytest.mark.parametrize('command', ['publication', 'arxiv'])
def test_publication_dry_run_outside_checkout_with_spaces(command, tmp_path):
    wd = tmp_path / 'directory with spaces'; wd.mkdir()
    result = subprocess.run([sys.executable, '-m', 'graphpathology', '--root', str(ROOT), command, '--dry-run'],
                            cwd=wd, text=True, capture_output=True, timeout=15)
    assert result.returncode == 0, result.stderr


def test_missing_current_source_is_specific_before_rendering(tmp_path):
    (tmp_path / 'manuscript').mkdir()
    (tmp_path / 'manuscript/PACKAGE_PROVENANCE.json').write_text(json.dumps({'files': {'main.tex': '0' * 64}}))
    with pytest.raises(ValueError, match='main.tex'): check_sources(tmp_path)
