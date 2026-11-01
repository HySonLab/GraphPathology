#!/usr/bin/env python3
"""Redraw Figures 4 and 5 from the supplied full-precision result tables.

The original plotted values and aggregation are retained. Only the layout,
lettering, and text size are changed to make the six panels readable at a
160-mm manuscript width. No model fitting or metric recomputation occurs.

Requirements: Python 3, NumPy, pandas, and Matplotlib.
Default invocation, from any directory:
    python scripts/redraw_external_figures.py
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


PACKAGE = Path(__file__).resolve().parents[1]
COLORS = {'M0': '#7b7b7b', 'M1': '#087e8b', 'M3': '#b55628'}
LABELS = {'M0': 'H&E only', 'M1': 'Sparse ST', 'M3': 'H&E + ST'}


def coverage_values(data, metric):
    """Original figure operation: median of patient summaries at each budget."""
    result = {}
    for method, group in data.groupby('method'):
        if method in COLORS:
            result[method] = group.groupby('coverage')[metric].median()
    return result


def distance_values(distance_data, cohort):
    """Original seed -> geometry -> patient aggregation, then cohort median."""
    selected = distance_data[
        (distance_data.cohort == cohort)
        & (distance_data.coverage == 0.2)
        & (distance_data.method == 'M1')
        & distance_data.regime.isin(['V1', 'V2'])
    ]
    result = {}
    for regime, group in selected.groupby('regime'):
        patient = (
            group.groupby(['patient', 'geometry', 'lower_um']).mse.mean()
            .groupby(['patient', 'lower_um']).mean().reset_index()
        )
        result[regime] = patient.groupby('lower_um').mse.median()
    return result


def draw_coverage(ax, values, ylabel, title):
    for method, series in values.items():
        ax.plot(series.index * 100, series.values, 'o-', color=COLORS[method],
                label=LABELS[method], ms=4, lw=1.3)
    ax.set(xlabel='Measured locations (%)', ylabel=ylabel, title=title)
    ax.legend(frameon=False, fontsize=8.8, handlelength=1.6)


def draw_external(data, distance_data, cohort, destination):
    fig, grid = plt.subplots(3, 2, figsize=(7.2, 7.8))
    axs = grid.ravel()
    primary = data[(data.regime == 'V1') & (data.coverage == 0.2)]

    # a. Patient observations and cohort medians at 20% distributed coverage.
    ax = axs[0]
    for column, method in enumerate(['M0', 'M1', 'M3']):
        values = primary[primary.method == method].median_pearson
        ax.scatter(np.repeat(column, len(values)), values, s=20,
                   color=COLORS[method], alpha=0.65)
        ax.plot([column - .2, column + .2], [values.median()] * 2,
                lw=2.5, color=COLORS[method])
    ax.set_xticks(range(3), ['H&E', 'ST', 'H&E + ST'])
    ax.set(ylabel='Median gene Pearson r',
           title=f'20% V1 · {primary.patient.nunique()} patients · 50 genes')

    # b, c. Identical cohort medians from the original coverage routines.
    distributed = data[data.regime == 'V1']
    draw_coverage(axs[1], coverage_values(distributed, 'median_pearson'),
                  'Median gene Pearson r', 'V1 coverage response')
    draw_coverage(axs[2], coverage_values(distributed, 'mean_mse'),
                  'Mean gene MSE', 'V1 error response')
    axs[2].set_yscale('log')

    # d. Effects are paired within patient before their cohort median.
    paired = primary.pivot(index='patient', columns='method', values='median_pearson')
    ax = axs[3]
    for column, (first, second) in enumerate([('M1', 'M0'), ('M3', 'M1')]):
        values = paired[first] - paired[second]
        ax.scatter(np.repeat(column, len(values)), values, s=20,
                   alpha=.6, color=COLORS[first])
        ax.plot([column - .2, column + .2], [values.median()] * 2,
                color=COLORS[first], lw=2.5)
    ax.axhline(0, color='black', lw=.7)
    ax.set_xticks([0, 1], ['ST − H&E', 'Fusion − ST'])
    ax.set(ylabel='Paired Δ median Pearson', title='Patient effects · 20% V1')

    # e. Each line connects the same patient's two protocol summaries.
    paired_regimes = data[(data.method == 'M1') & (data.coverage == .2)].pivot(
        index='patient', columns='regime', values='median_pearson')
    ax = axs[4]
    for _, row in paired_regimes.iterrows():
        ax.plot([0, 1], row[['V1', 'V2']], 'o-', alpha=.4,
                color=COLORS['M1'], ms=3, lw=1.2)
    ax.set_xticks([0, 1], ['V1', 'V2'])
    ax.set(ylabel='Median gene Pearson r', title='Paired patient protocol contrast')

    # f. Occupied-bin medians; no interpolation over missing distance bins.
    ax = axs[5]
    for regime, series in distance_values(distance_data, cohort).items():
        ax.plot(series.index, series.values, 'o-', label=regime,
                color=COLORS['M1'] if regime == 'V1' else '#8d4ca6',
                ms=4, lw=1.3)
    ax.set(xlabel='Distance-bin lower bound (µm)', ylabel='Mean squared error',
           title='Nearest-anchor distance · ST only')
    ax.legend(frameon=False, fontsize=9, handlelength=1.6)

    for letter, ax in zip('abcdef', axs):
        ax.text(-.22, 1.13, letter, transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='top', ha='left')
        ax.tick_params(axis='both', labelsize=9)
    fig.tight_layout(pad=1.25, w_pad=1.8, h_pad=2.1)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Fixed canvas retains the intended physical font size on PDF/SVG export.
    fig.savefig(destination.with_suffix('.pdf'))
    fig.savefig(destination.with_suffix('.svg'))
    plt.close(fig)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=PACKAGE / 'source_data')
    parser.add_argument('--distance-data', type=Path,
                        help='Optional override for the distance-bin CSV.')
    parser.add_argument('--output-dir', type=Path, default=PACKAGE / 'figures')
    args = parser.parse_args()
    distance_path = args.distance_data or args.data_dir / 'SuppTable7_distance_bins.csv'
    distance_data = pd.read_csv(distance_path, low_memory=False)
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.titlesize': 10, 'axes.labelsize': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
        'svg.fonttype': 'none', 'pdf.fonttype': 42,
        'savefig.dpi': 300,
    })
    tasks = [('Fig5_source.csv', 'Wu fresh validation', 'Figure4'),
             ('Fig6_source.csv', 'Bassiouni second replication', 'Figure5')]
    for filename, cohort, figure in tasks:
        data = pd.read_csv(args.data_dir / filename)
        output = draw_external(data, distance_data, cohort, args.output_dir / figure)
        print(f'{output.name}: PDF and SVG exported from {filename}')


if __name__ == '__main__':
    main()
