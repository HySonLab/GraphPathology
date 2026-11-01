#!/usr/bin/env python3
"""Export the editable result-figure SVG files to vector PDF.

Requirements: Python 3 and Inkscape (available on PATH).
Run from any directory:
    python scripts/prepare_existing_figures.py

This exports the supplied artwork without recalculating results or modifying
data marks. Regenerate the native SVGs with their corresponding drawing scripts first,
if needed. Figure 2 is the new synthetic-geometry diagram, not a development plot.
"""

from pathlib import Path
import shutil
import subprocess


def main():
    package = Path(__file__).resolve().parents[1]
    renderer = shutil.which('inkscape')
    if renderer is None:
        raise SystemExit('Install Inkscape to export SVG artwork. Supplied PDFs are ready to compile.')
    figures = [package / 'figures' / f'Figure{i}.svg' for i in range(1, 7)]
    figures += [package / 'supplementary_figures' / f'FigureS{i}.svg' for i in range(1, 13)]
    for source in figures:
        if not source.is_file():
            raise FileNotFoundError(source)
        destination = source.with_suffix('.pdf')
        subprocess.run([renderer, str(source), '--export-type=pdf',
                        '--export-filename=' + str(destination)], check=True)
        print(destination.relative_to(package))


if __name__ == '__main__':
    main()
