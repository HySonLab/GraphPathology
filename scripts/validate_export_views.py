"""Validate current publication outputs (the initial 26-view format is superseded)."""
import argparse
import contextlib
import sys
import json
from pathlib import Path
from graphpathology.publication import check_sources
from graphpathology.publication_regression import validate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--output', type=Path, required=True, help='Existing all-29-target publication output')
    args = p.parse_args()
    out = args.output.resolve()
    if not (out / 'REGRESSION.json').is_file():
        p.error('Expected current publication output; regenerate the superseded partial export views first')
    reference, _ = check_sources(args.root.resolve())
    baseline = out / 'revalidation-baseline'
    if baseline.exists():
        p.error('Refusing to overwrite an existing revalidation-baseline')
    with contextlib.redirect_stdout(sys.stderr):
        report = validate(args.root.resolve(), out / 'publication', reference, baseline)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
