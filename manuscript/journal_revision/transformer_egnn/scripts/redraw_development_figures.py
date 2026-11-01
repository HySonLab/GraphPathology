#!/usr/bin/env python3
"""Redraw publication Figures 3 and 6 from archived scientific CSV exports.

Requirements: Python 3, NumPy, pandas, and Matplotlib.
Run from any directory:
    python scripts/redraw_development_figures.py

All selections, patient-summary aggregation, plotted values, spatial
coordinates, and map color limits follow the original figure routines.
Only presentation is revised. The spatial maps retain their original
rasterized point representation within otherwise vector PDF/SVG figures.
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
    return {method: group.groupby('coverage')[metric].median()
            for method, group in data.groupby('method') if method in COLORS}


def distance_values(data):
    selected = data[(data.cohort == 'Development') & (data.coverage == .2)
                    & (data.method == 'M1') & data.regime.isin(['V1', 'V2'])]
    result = {}
    for regime, group in selected.groupby('regime'):
        patient = (group.groupby(['patient', 'geometry', 'lower_um']).mse.mean()
                   .groupby(['patient', 'lower_um']).mean().reset_index())
        result[regime] = patient.groupby('lower_um').mse.median()
    return result


def panel_letters(fig, axes, map_panels=False):
    """Place letters in fixed subplot margins, including equal-aspect maps."""
    for index, ax in enumerate(axes.ravel()):
        if map_panels and index >= 2:
            position = ax.get_subplotspec().get_position(fig)
            fig.text(position.x0 - .055, position.y1 + .012, chr(97 + index),
                     fontsize=12, fontweight='bold', va='top')
        else:
            ax.text(-.22, 1.13, chr(97 + index), transform=ax.transAxes,
                    fontsize=12, fontweight='bold', va='top', ha='left')


def export(fig, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination.with_suffix('.pdf'))
    fig.savefig(destination.with_suffix('.svg'))
    plt.close(fig)


def draw_development(data, destination):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.7))
    for column, regime in enumerate(['V1', 'V2']):
        for row, metric in enumerate(['median_pearson', 'mean_mse']):
            ax = axes[row, column]
            values = coverage_values(data[data.regime == regime], metric)
            for method, series in values.items():
                ax.plot(series.index * 100, series.values, 'o-',
                        color=COLORS[method], label=LABELS[method], ms=4, lw=1.3)
            ax.set(xlabel='Measured locations (%)',
                   ylabel='Median gene Pearson r' if row == 0 else 'Mean gene MSE',
                   title=('Random' if regime == 'V1' else 'Block') + ' masks · TENX99')
            ax.legend(frameon=False, fontsize=8.5, handlelength=1.6)
            ax.tick_params(labelsize=9)
    panel_letters(fig, axes)
    fig.tight_layout(pad=1.25, w_pad=1.8, h_pad=2.1)
    export(fig, destination)


def draw_geometry(distance, sampling, maps, destination):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0),
                             gridspec_kw={'height_ratios': [1, 1.1]})
    # a. The original patient/geometry/seed-aggregated distance values.
    ax = axes[0, 0]
    for regime, series in distance_values(distance).items():
        ax.plot(series.index, series.values, 'o-',
                label={'V1': 'Random', 'V2': 'Block'}[regime],
                color='#087e8b' if regime == 'V1' else '#8d4ca6', ms=4, lw=1.3)
    ax.set(xlabel='Distance-bin lower bound (µm)', ylabel='Mean squared error',
           title='Nearest-anchor distance · ST only')
    ax.legend(frameon=False, fontsize=8.8, handlelength=1.6)

    # b. Original baseline-policy selection and arithmetic means.
    chosen = sampling[
        (sampling.phase == 'primary') & (sampling.panel == 'S2-50')
        & (sampling.encoder == 'resnet18_reference') & (sampling.role == 'baseline')
        & sampling.control.isna() & (sampling.regime == 'open')
        & sampling.strategy.isin(['A0', 'A1'])
    ]
    ax = axes[0, 1]
    for policy, group in chosen.groupby('strategy'):
        series = group.groupby('nominal_budget').mean_mse.mean()
        ax.plot(series.index * 100, series.values, 'o-', ms=4, lw=1.3,
                color={'A0': '#1f77b4', 'A1': '#ff7f0e'}[policy],
                label={'A0': 'Random · 10-seed mean',
                       'A1': 'Spatial FPS · deterministic'}[policy])
    ax.set(xlabel='Measured locations (%)', ylabel='Mean gene MSE',
           title='Secondary sampling comparison')
    ax.legend(frameon=False, fontsize=8.2, handlelength=1.5)

    # c,d. Shared original 98th-percentile cap uses the entire exported map
    # table, including the additional capture-window source values.
    cap = maps.absolute_error.quantile(.98)
    for ax, regime in zip(axes[1], ['random', 'block']):
        selected = maps[maps.regime == regime]
        observed = selected[selected.observed]
        hidden = selected[~selected.observed]
        points = ax.scatter(hidden.x_um, hidden.y_um, c=hidden.absolute_error,
                            s=1, cmap='magma', vmin=0, vmax=cap, rasterized=True)
        ax.scatter(observed.x_um, observed.y_um, s=.4, color='#53d5d0',
                   alpha=.4, rasterized=True)
        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set(xlabel='x (µm)', ylabel='y (µm)',
               title=f'{regime.capitalize()} · {hidden.gene.iloc[0]} · fixed fusion')
        ax.title.set_position((.5, 1.06))
        colorbar = fig.colorbar(points, ax=ax, pad=.035, fraction=.052,
                               label='Absolute error (98th percentile cap)')
        colorbar.ax.tick_params(labelsize=8.5)
        colorbar.ax.yaxis.label.set_size(8.5)

    for ax in axes.ravel():
        ax.tick_params(labelsize=9)
    fig.tight_layout(pad=1.25, w_pad=1.8, h_pad=2.1)
    panel_letters(fig, axes, map_panels=True)
    export(fig, destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=PACKAGE / 'source_data')
    parser.add_argument('--output-dir', type=Path, default=PACKAGE / 'figures')
    args = parser.parse_args()
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.titlesize': 10, 'axes.labelsize': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'savefig.dpi': 300,
    })
    load = lambda name: pd.read_csv(args.data_dir / name, low_memory=False)
    draw_development(load('development_patient.csv'), args.output_dir / 'Figure3')
    draw_geometry(load('SuppTable7_distance_bins.csv'), load('SFig6_source.csv'),
                  load('Fig3_map_source.csv'), args.output_dir / 'Figure6')
    print('Figures 3 and 6: PDF and SVG exported from exact source tables.')


if __name__ == '__main__':
    main()
