"""Isolated publication regeneration, document builds and arXiv source packaging."""
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile
from .assets import new_output, sha

RENDERERS = ['draw_architecture.py', 'draw_graph_architecture.py',
             'redraw_development_figures.py', 'redraw_external_figures.py',
             'redraw_supplementary_figures.py']


def targets(root):
    return list(csv.DictReader((root / 'docs/REPRODUCTION_TARGETS.csv').open()))


def check_sources(root):
    package = root / 'manuscript/revised'
    provenance = json.loads((root / 'manuscript/PACKAGE_PROVENANCE.json').read_text())
    for name, digest in provenance['files'].items():
        if not (package / name).is_file() or sha(package / name) != digest:
            raise ValueError('Author package changed or missing: ' + name)
    rows = targets(root)
    if len(rows) != 29 or len({r['target'] for r in rows}) != 29:
        raise ValueError('Expected 29 distinct publication targets')
    for row in rows:
        for key in ['latex_or_graphic_source', 'numerical_or_synthetic_source']:
            for name in row[key].split(';'):
                if not (package / name.strip()).exists():
                    raise FileNotFoundError(row['target'] + ': ' + name)
    return package, rows


def copy_package(package, destination):
    shutil.copytree(package, destination)
    return destination


def render(package):
    env = dict(os.environ, MPLCONFIGDIR=str(package / '.mpl-cache'),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for script in RENDERERS:
        result = subprocess.run([sys.executable, '-m', 'graphpathology._publication_worker',
                                 '--package', str(package), '--script', script],
                                text=True, capture_output=True, timeout=120, env=env)
        (package / (script + '.log')).write_text(result.stdout + result.stderr)
        if result.returncode:
            raise ValueError(f'Renderer failed: {script}; see {package / (script + ".log")}')
        print('Rendered ' + script, flush=True)


def run(root, args):
    package, rows = check_sources(root)
    if args.target != 'all':
        raise ValueError('Publication regeneration is atomic across all 29 targets; use --target all')
    if args.dry_run:
        print(json.dumps({'validated': True, 'targets': 29, 'model_inputs_required': False})); return 0
    if not args.output:
        raise ValueError('--output is required')
    out = new_output(root, args.output)
    generated = copy_package(package, out / 'publication')
    render(generated)
    # Only record success after numerical regression; no partial views masquerade as baseline reproduction.
    from .publication_regression import validate
    report = validate(root, generated, reference=package, render_reference=out / 'baseline')
    (out / 'REGRESSION.json').write_text(json.dumps(report, indent=2) + '\n')
    statuses = [{'target': r['target'], 'status': 'validated', 'source': r['numerical_or_synthetic_source'],
                 'output': r['latex_or_graphic_source']} for r in rows]
    (out / 'TARGET_STATUS.json').write_text(json.dumps(statuses, indent=2) + '\n')
    print(json.dumps({'targets_validated': 29, 'output': str(out), 'checks': len(report['checks'])}))
    return 0


def compile_documents(package, documents=('main', 'article_only', 'supplementary'), bibliography=True):
    if not shutil.which('pdflatex'):
        raise FileNotFoundError('Install TeX Live with pdflatex, latex-extra, fonts-recommended and bibtex8/bibtex')
    bib = shutil.which('bibtex8') or shutil.which('bibtex')
    if bibliography and not bib:
        raise FileNotFoundError('Install bibtex8 or bibtex')
    reports = {}
    for name in documents:
        commands = [['pdflatex', '-no-shell-escape', '-recorder', '-interaction=nonstopmode', '-halt-on-error', name + '.tex']]
        if bibliography and name != 'supplementary':
            commands.append([bib, name])
        commands += [commands[0], commands[0]]
        for index, cmd in enumerate(commands):
            result = subprocess.run(cmd, cwd=package, text=True, capture_output=True, timeout=90)
            (package / f'{name}.pass{index}.txt').write_text(result.stdout + result.stderr)
            if result.returncode:
                raise ValueError(f'Build failed for {name}; see {package / (name + ".log")}')
        log = (package / (name + '.log')).read_text(errors='replace')
        diagnostics = [line for line in log.splitlines() if re.search(r'Warning|undefined|Overfull|Underfull|destination.*duplicate', line)]
        if diagnostics:
            raise ValueError(name + ' final TeX diagnostics: ' + '; '.join(diagnostics))
        info = subprocess.check_output(['pdfinfo', str(package / (name + '.pdf'))], text=True)
        reports[name] = {'pages': int(re.search(r'Pages:\s+(\d+)', info).group(1)), 'diagnostics': []}
    return reports


def check_documents(package):
    text = subprocess.check_output(['pdftotext', '-layout', str(package / 'main.pdf'), '-'], text=True)
    figures = set(re.findall(r'\bFigure\s+(S?\d+)\s*:', text))
    tables = set(re.findall(r'\bTable\s+(S?\d+)\s*:', text))
    if figures != {str(n) for n in range(1, 7)} | {'S' + str(n) for n in range(1, 13)}:
        raise ValueError('Incorrect combined figure numbering')
    if tables != {'1', '2'} | {'S' + str(n) for n in range(1, 10)}:
        raise ValueError('Incorrect combined table numbering')
    if len(re.findall(r'\\bibitem', (package / 'main.bbl').read_text())) != 29:
        raise ValueError('Expected 29 bibliography entries')
    supplement = next(i for i, page in enumerate(text.split('\f')) if 'Supplementary Information' in page and 'Contents and source-data guide' in page)
    before = '\f'.join(text.split('\f')[:supplement])
    if '[29]' not in before:
        raise ValueError('Supplement must follow the completed References')
    for index, page in enumerate(text.split('\f')):
        if not page.strip():
            continue
        if page.strip().splitlines()[-1].strip() != str(index + 1):
            raise ValueError('Combined page numbering is not continuous')
    return {'figures': 18, 'tables': 11, 'references': 29, 'supplement_start_page': supplement + 1,
            'continuous_pages': True, 'duplicate_destinations': False}


def build(root, args):
    package, _ = check_sources(root)
    if args.dry_run:
        if not all(shutil.which(x) for x in ['pdflatex', 'pdfinfo', 'pdftotext']):
            raise FileNotFoundError('Install pdflatex and Poppler (pdfinfo/pdftotext)')
        print(json.dumps({'validated': True, 'documents': ['main', 'article_only', 'supplementary']})); return 0
    if not args.output:
        raise ValueError('--output is required')
    out = new_output(root, args.output)
    project = copy_package(package, out / 'project')
    if args.publication:
        generated = Path(args.publication).resolve() / 'publication'
        if not (Path(args.publication) / 'REGRESSION.json').is_file():
            raise ValueError('--publication must identify a numerically validated regeneration output')
        for folder in ['figures', 'supplementary_figures']:
            for src in (generated / folder).glob('*.pdf'):
                shutil.copy2(src, project / folder / src.name)
    report = {'documents': compile_documents(project), 'structure': check_documents(project)}
    for name in ['main', 'article_only', 'supplementary']:
        shutil.copy2(project / (name + '.pdf'), out / (name + '.pdf'))
    (out / 'BUILD_REPORT.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report)); return 0


def arxiv_dependencies(package):
    pending = ['main.tex']; selected = set()
    while pending:
        name = pending.pop()
        if name in selected:
            continue
        path = package / name
        if not path.is_file():
            raise FileNotFoundError('TeX dependency missing: ' + name)
        selected.add(name)
        if path.suffix == '.tex':
            source = path.read_text()
            for command, dep in re.findall(r'\\(input|include|includegraphics|bibliography|bibliographystyle)(?:\[[^]]*\])?\{([^}]+)\}', source):
                suffix = '.pdf' if command == 'includegraphics' else '.bib' if command == 'bibliography' else '.bst' if command == 'bibliographystyle' else '.tex'
                pending.append(dep if Path(dep).suffix else dep + suffix)
    selected.add('main.bbl')
    return selected

def arxiv(root, args):
    package, _ = check_sources(root)
    if args.dry_run:
        print('Validated: one main.tex root; frozen figures; main.bbl; numerical ancillary data under anc/'); return 0
    if not args.output:
        raise ValueError('--output is required')
    out = new_output(root, args.output)
    staging = out / 'source'; staging.mkdir()
    # Compute the TeX input/graphics dependency closure. \\path strings are prose, not dependencies.
    selected = arxiv_dependencies(package)
    for name in sorted(selected):
        dst = staging / name; dst.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(package / name, dst)
    for folder in ['source_data', 'methods_data']:
        shutil.copytree(package / folder, staging / 'anc' / folder)
    if any('anc/' in (staging / n).read_text() for n in selected if n.endswith('.tex')):
        raise ValueError('Ancillary numerical files must not be TeX dependencies')
    zip_path = out / 'GraphPathology_arxiv_source.zip'
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
        for src in sorted(staging.rglob('*')):
            if src.is_file(): z.write(src, src.relative_to(staging))
    verification = out / 'clean-extracted'
    with zipfile.ZipFile(zip_path) as z: z.extractall(verification)
    report = {'single_tex_root': 'main.tex', 'tex_dependencies': len(selected),
              'clean_extracted_build': compile_documents(verification, ('main',), bibliography=False),
              'structure': check_documents(verification)}
    (out / 'ARXIV_REPORT.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'zip': str(zip_path), **report})); return 0
