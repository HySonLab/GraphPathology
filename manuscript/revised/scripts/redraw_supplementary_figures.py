#!/usr/bin/env python3
"""Regenerate supplementary figures from the accompanying exact CSV exports.

No model fitting, inference, resampling, or additional statistical tests are run.
The data filters, group means/medians, patient pairing, and plotted values follow
the original exported figure routines. Changes concern figure size, typography,
scientific labels, and panel arrangement only.

Usage (from any directory):
    python scripts/redraw_supplementary_figures.py
    python scripts/redraw_supplementary_figures.py --figure S2

Dependencies: numpy, pandas, matplotlib.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "source_data"
OUTPUT = ROOT / "supplementary_figures"
COLORS = {"M0": "#7b7b7b", "M1": "#087e8b", "M3": "#b55628"}
TEAL = COLORS["M1"]
ORANGE = COLORS["M3"]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.labelsize": 9.5,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.dpi": 300,
})

TRACES: dict[str, list[dict]] = {}


def load(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / (name + ".csv"), low_memory=False)


def plot_values(ax: plt.Axes) -> dict:
    """Record only numerical artist arrays, independently of presentation."""
    return {
        "lines": [
            {"x": np.asarray(line.get_xdata(), dtype=float).tolist(),
             "y": np.asarray(line.get_ydata(), dtype=float).tolist()}
            for line in ax.lines
        ],
        "points": [
            np.ma.asarray(collection.get_offsets(), dtype=float).filled(np.nan).tolist()
            for collection in ax.collections
        ],
    }


def label_panels(axes: list[plt.Axes], letters: str | None = None) -> None:
    letters = letters or "abcdefghijklmnopqrstuvwxyz"[:len(axes)]
    for ax, letter in zip(axes, letters):
        offset = -0.085 if ax.get_position().width > 0.6 else -0.19
        ax.text(offset, 1.065, letter, transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="bottom", ha="left")


def save(fig: plt.Figure, axes: list[plt.Axes], number: int) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    TRACES[f"S{number}"] = [plot_values(ax) for ax in axes]
    for extension in ("pdf", "svg"):
        fig.savefig(OUTPUT / f"FigureS{number}.{extension}", facecolor="white")
    plt.close(fig)


def category_label(label: str) -> str:
    """Wrap labels without changing their sorting or their associated values."""
    names = {
        "shuffled_graph / GCN": "Shuffled\nGCN",
        "shuffled_graph / MEAN": "Shuffled\nMEAN",
        "shuffled_graph / SAGE": "Shuffled\nSAGE",
        "Local MLP": "Local\nMLP",
        "Degree rewire": "Degree\nrewire",
        "Degree+length": "Degree +\nlength",
        "Mean pool": "Mean\npool",
        "Self-only": "Self-\nonly",
        "True SAGE": "True\nSAGE",
        "Vertex perm.": "Vertex\nperm.",
    }
    return names.get(label, label.replace(" / ", "\n"))


def grouped(ax: plt.Axes, data: pd.DataFrame, group: str,
            metric: str, ylabel: str) -> None:
    values, labels = [], []
    # This is the original grouping and arithmetic-mean operation.
    for label, frame in data.groupby(group, dropna=False):
        q = frame[metric].dropna()
        if len(q):
            values.append(q)
            labels.append(str(label))
    for i, q in enumerate(values):
        ax.scatter(np.repeat(i, len(q)), q, s=18, alpha=0.5, color=TEAL)
        ax.plot([i - 0.2, i + 0.2], [q.mean()] * 2, c=ORANGE, lw=2)
    ax.set_xticks(range(len(labels)), [category_label(x) for x in labels])
    ax.set_ylabel(ylabel)
    ax.margins(x=0.09)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))


def simple_layout(title: str, height: float = 3.65):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, height))
    fig.suptitle(title, y=0.975, fontsize=9.5)
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.19,
                        top=0.79, wspace=0.40)
    return fig, list(axes)


def her2_calibration() -> None:
    # Select the fixed-fusion calibration comparison by its scientific role.
    d = load("Fig4_source")
    d = d[(d.coverage == 0.2) & (d.cohort == "HER2ST calibrated")]
    fig, grid = plt.subplots(2, 2, figsize=(7.2, 5.8))
    fig.subplots_adjust(left=0.115, right=0.985, bottom=0.08,
                        top=0.91, hspace=0.63, wspace=0.43)
    for j, regime in enumerate(["V1", "V2"]):
        for i, metric in enumerate(["median_pearson", "mean_mse"]):
            ax = grid[i, j]
            for k, transform in enumerate(["T0", "T1"]):
                medians = []
                for z, model in enumerate(["E0", "E1", "E2"]):
                    q = d[(d.regime == regime) &
                          (d.method == transform + "_" + model)][metric]
                    medians.append(q.median())
                    ax.scatter(np.repeat(z + (k - 0.5) * 0.24, len(q)), q,
                               s=15, alpha=0.4, color=["#7b7b7b", ORANGE][k])
                ax.plot(np.arange(3) + (k - 0.5) * 0.24, medians, "o-",
                        label="Zero-shot" if k == 0 else "T1 calibrated",
                        color=["#7b7b7b", ORANGE][k], markersize=4)
            ax.set_xticks(range(3), ["H&E", "ST", "H&E + ST"])
            ax.set_title(regime + " · 8 patients · 40 genes", pad=9)
            ax.set_ylabel("Median gene Pearson r" if i == 0 else
                          "Mean gene MSE\nCommon target units")
            ax.legend(frameon=False, handlelength=1.4, labelspacing=0.35)
            if metric == "mean_mse":
                ax.set_yscale("log")
    axes = list(grid.ravel())
    label_panels(axes)
    save(fig, axes, 1)


def cross_cohort() -> None:
    # This cleaned table contains the exact cross-cohort source values.
    d = load("SuppTable2_per_patient")
    cohorts = ["Development", "HER2ST zero-shot", "HER2ST calibrated",
               "Wu fresh validation", "Bassiouni second replication"]
    short_titles = ["Development", "HER2ST zero-shot", "HER2ST calibrated",
                    "Wu fresh validation", "Bassiouni replication"]
    fig, grid = plt.subplots(5, 2, figsize=(7.2, 8.5))
    fig.subplots_adjust(left=0.115, right=0.985, bottom=0.065,
                        top=0.96, hspace=0.92, wspace=0.47)
    for row, cohort in enumerate(cohorts):
        p = d[(d.cohort == cohort) & (d.coverage == 0.2) &
              (d.regime == "V1")].pivot(index="patient", columns="method",
                                        values="median_pearson")
        ax = grid[row, 0]
        for k, (a, b) in enumerate([("M1", "M0"), ("M3", "M1")]):
            values = p[a] - p[b]
            ax.scatter(np.repeat(k, len(values)), values, s=16,
                       c=COLORS[a], alpha=0.65)
            ax.plot([k - 0.2, k + 0.2], [values.median()] * 2,
                    c=COLORS[a], lw=2)
        ax.axhline(0, lw=0.6, c="k")
        ax.set_xticks([0, 1], ["ST − H&E", "Fusion − ST"])
        genes = 40 if cohort.startswith("HER2") else 50
        ax.set_title(f"{short_titles[row]} · {genes} genes · n={len(p)}",
                     fontsize=9.2, pad=7)
        ax.set_ylabel("Paired Δ\nmedian Pearson", fontsize=9)
        p = d[(d.cohort == cohort) & (d.method == "M1") &
              (d.regime == "V1") & d.coverage.isin([0.1, 0.2, 0.4])]
        ax = grid[row, 1]
        for _, group in p.groupby("patient"):
            ax.plot(group.coverage * 100, group.median_pearson, "o-",
                    c=TEAL, alpha=0.35, ms=3)
        ax.set_ylabel("Median gene\nPearson r", fontsize=9)
        ax.set_xticks([10, 20, 40])
        if row == 4:
            ax.set_xlabel("Measured locations (%)", fontsize=9)
        for ax in grid[row]:
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    # Preserve the original a-e effects / f-j coverage panel identities.
    axes = list(grid[:, 0]) + list(grid[:, 1])
    label_panels(axes)
    save(fig, axes, 2)


def development_and_controls(number: int) -> None:
    source_number = number - 2
    d = load(f"SFig{source_number}_source")
    titles = {
        3: "Morphology baselines · TENX99 · 50 genes\nModel seeds and controls",
        4: "Graph context · TENX99 · 50 genes\nModel seeds and graph controls",
        5: "Topology controls · TENX99 · 50 genes\nTraining-seed means; five graph seeds",
        6: "Image representations · TENX99 · local MLP\n50 genes · three training seeds",
        7: "Sparse-reconstruction controls · 20% · 50 genes\nOne random mask · three permutation seeds",
        8: "Sampling policies · 20% coverage\nA0 random · A1 spatial · A2 morphology · A3 joint",
        9: "HER2ST zero-shot · 8 patients · 40 genes\nE0 H&E · E1 ST · E2 fusion",
        10: "HER2ST selected fusion · 8 development patients\n40 genes · 20% coverage",
        11: "Wu controls · one patient · 20% coverage\nThree seeds × two regimes",
    }
    if number == 5:
        fig, grid = plt.subplots(2, 1, figsize=(7.2, 5.5))
        axes = list(grid)
        fig.suptitle(titles[number], y=0.975, fontsize=9.5)
        fig.subplots_adjust(left=0.13, right=0.985, bottom=0.10,
                            top=0.85, hspace=0.57)
    else:
        fig, axes = simple_layout(titles[number])
    if source_number in [1, 2]:
        d = d[d.panel == "S2-50"]
        a = d.copy()
        a["label"] = a.model
        if source_number == 2:
            base = load("SFig1_source")
            base = base[(base.panel == "S2-50") & (base.model == "B2")].copy()
            base["label"] = "Local MLP"
            a = pd.concat([base, a])
        grouped(axes[0], a, "label", "median_pearson", "Median gene Pearson r")
        c = load(f"SFig{source_number}_controls_source")
        c = c[c.panel == "S2-50"].copy()
        c["label"] = c.phase.astype(str) + " / " + c.model.astype(str)
        grouped(axes[1], c, "label", "mean_mse", "Mean gene MSE")
    elif source_number == 3:
        d = load("SFig3_panel_source")
        d = d[(d.phase == "primary") & (d.panel == "S2-50")]
        p = d.groupby(["control", "graph_randomization_seed"])[
            ["median_pearson", "mean_mse"]].mean().reset_index()
        p["label"] = p.control.map({
            "A0": "Local", "A1": "Self-only", "A2": "Mean pool",
            "A3": "True SAGE", "A4": "Vertex perm.",
            "A5": "Degree rewire", "A6": "Degree+length"})
        grouped(axes[0], p, "label", "median_pearson", "Median gene Pearson r")
        grouped(axes[1], p, "label", "mean_mse", "Mean gene MSE")
        for ax in axes:
            ax.margins(x=0.025)
    elif source_number == 4:
        d = d[(d.panel == "S2-50") & (d.phase == "primary") &
              (d.model == "D1")].copy()
        d["label"] = d.encoder.map({"pathology_encoder_1": "Phikon",
                                    "resnet18_reference": "ResNet18"})
        grouped(axes[0], d, "label", "median_pearson", "Median gene Pearson r")
        grouped(axes[1], d, "label", "mean_mse", "Mean gene MSE")
    elif source_number == 5:
        d = d[(d.panel == "S2-50") & (d.encoder == "resnet18_reference") &
              (d.nominal_coverage == 0.2)].copy()
        d["label"] = d.control.astype(str) + " / " + d.model.astype(str)
        grouped(axes[0], d, "label", "median_pearson", "Median gene Pearson r")
        grouped(axes[1], d, "label", "mean_mse", "Mean gene MSE")
    elif source_number == 6:
        d = d[(d.panel == "S2-50") & (d.phase == "primary") &
              (d.encoder == "resnet18_reference") & (d.nominal_budget == 0.2) &
              (d.role == "baseline") & d.control.isna()]
        for ax, regime in zip(axes, ["open", "contiguous"]):
            grouped(ax, d[d.regime == regime], "strategy", "median_pearson",
                    "Median gene Pearson r")
            ax.set_title(regime, pad=8)
    elif source_number == 7:
        d = d[(d.panel == "S2-50") & (d.nominal_budget == 0.2) &
              (d.strategy == "random") & d.model.isin(["E0", "E1", "E2"])]
        for ax, regime in zip(axes, ["V1", "V2"]):
            grouped(ax, d[d.regime == regime], "model", "median_pearson",
                    "Median gene Pearson r")
            ax.set_title(regime, pad=8)
    elif source_number == 8:
        d = d[d.budget == 0.2]
        grouped(axes[0], d, "regime", "gain_pearson",
                "Selected fusion − ST\nΔ median Pearson")
        grouped(axes[1], d, "regime", "relative_mse_reduction",
                "Relative MSE reduction\nSelected fusion vs ST")
    elif source_number == 9:
        d = d[d.method.isin(["M1", "M3"])].copy()
        d["label"] = d.control + " / " + d.method
        grouped(axes[0], d, "label", "delta_pearson",
                "Control − matched clean\nΔ Pearson")
        grouped(axes[1], d, "label", "mse_reduction",
                "Matched clean − control\nMean MSE")
    label_panels(axes)
    save(fig, axes, number)


def replication_controls() -> None:
    fig = plt.figure(figsize=(7.2, 5.6))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1],
                           left=0.13, right=0.985, bottom=0.13,
                           top=0.84, hspace=0.75, wspace=0.42)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]),
            fig.add_subplot(grid[1, :])]
    fig.suptitle("Bassiouni controls · one patient · 20% coverage\n"
                 "Three seeds × two regimes; image quality: 22 patients",
                 y=0.975, fontsize=9.5)
    d = load("SFig10_source")
    d = d[d.method.isin(["M1", "M3"])].copy()
    d["label"] = d.control + " / " + d.method
    grouped(axes[0], d, "label", "delta_pearson",
            "Control − matched clean\nΔ Pearson")
    grouped(axes[1], d, "label", "mse_reduction",
            "Matched clean − control\nMean MSE")
    q = load("SFig10_QC_source")
    axes[2].scatter(range(len(q)), q.native_patch_pixels, s=22, color=TEAL)
    axes[2].set(xlabel="Patient (source order)",
                ylabel="Native pixels per\n112-µm patch",
                title="Public image resolution")
    axes[2].set_xticks(range(len(q)), q.patient_id.str.replace("BASS_", ""),
                      rotation=90)
    axes[2].yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    label_panels(axes[:2])
    axes[2].text(-0.085, 1.065, "c", transform=axes[2].transAxes,
                 fontsize=12, fontweight="bold", va="bottom", ha="left")
    save(fig, axes, 12)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", default="all", help="all or S1–S12")
    parser.add_argument("--audit-json", type=Path,
                        help="Optional exact plotted-array record")
    args = parser.parse_args()
    selected = args.figure.upper().removeprefix("S")
    for number in range(1, 13):
        if selected not in {"ALL", str(number)}:
            continue
        if number == 1:
            her2_calibration()
        elif number == 2:
            cross_cohort()
        elif number == 12:
            replication_controls()
        else:
            development_and_controls(number)
    if not TRACES:
        parser.error("--figure must be all or S1–S12")
    if args.audit_json:
        args.audit_json.parent.mkdir(parents=True, exist_ok=True)
        args.audit_json.write_text(json.dumps(TRACES, indent=2) + "\n")
    print("Regenerated " + ", ".join(TRACES) + " as vector PDF and SVG.")


if __name__ == "__main__":
    main()
