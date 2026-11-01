"""Run one reviewed author renderer and capture the actual numerical artists."""
import argparse
import json
import runpy
import sys
from pathlib import Path
import numpy as np


def array(value):
    if value is None:
        return None
    a = np.ma.asarray(value, dtype=float)
    return {'values': a.filled(np.nan).tolist(), 'mask': np.ma.getmaskarray(a).tolist()}


def snapshot(fig):
    axes = []
    for ax in fig.axes:
        axes.append({
            'title': ax.get_title(), 'xlabel': ax.get_xlabel(), 'ylabel': ax.get_ylabel(),
            'xticklabels': [t.get_text() for t in ax.get_xticklabels()],
            'yticklabels': [t.get_text() for t in ax.get_yticklabels()],
            'xscale': ax.get_xscale(), 'yscale': ax.get_yscale(),
            'lines': [{'x': array(l.get_xdata()), 'y': array(l.get_ydata())} for l in ax.lines],
            'collections': [{'offsets': array(c.get_offsets()), 'array': array(c.get_array()),
                             'clim': list(c.get_clim())} for c in ax.collections],
            'images': [{'array': array(i.get_array()), 'clim': list(i.get_clim())} for i in ax.images],
            'texts': [t.get_text() for t in ax.texts],
            'patches': [{'vertices': array(p.get_path().vertices), 'transform': array(p.get_transform().get_matrix())} for p in ax.patches],
        })
    return axes


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--package', type=Path, required=True)
    p.add_argument('--script', required=True)
    args = p.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.figure import Figure
    original = Figure.savefig
    traces = {}

    def save(fig, fname, *a, **kw):
        # Capture both file exports and the schematic's atomic in-memory PDF export.
        if isinstance(fname, (str, Path)) and Path(fname).suffix == '.pdf':
            traces[Path(fname).stem] = snapshot(fig)
        elif kw.get('format') == 'pdf' and args.script == 'draw_architecture.py':
            traces['Figure1'] = snapshot(fig)
        return original(fig, fname, *a, **kw)

    Figure.savefig = save
    sys.argv = [str(args.package / 'scripts' / args.script)]
    runpy.run_path(sys.argv[0], run_name='__main__')
    path = args.package / 'NUMERICAL_TRACES.json'
    previous = json.loads(path.read_text()) if path.exists() else {}
    previous.update(traces)
    path.write_text(json.dumps(previous, allow_nan=True) + '\n')


if __name__ == '__main__':
    main()
