"""Plot the surface demo outputs.

Reads the CSV files written by ``surface_demo`` and produces three figures:

1. ``fig_surface.png``    calibrated total-variance surface + market quotes
2. ``fig_residuals.png``  fit residuals by expiry and model-vs-market scatter
3. ``fig_sweep.png``      out-of-sample error and arbitrage violations vs noise

Run:  python3 plot.py --data-dir outputs --out-dir outputs
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PALETTE = ["#2b6cb0", "#c05621", "#2f855a", "#6b46c1"]
SCHEME_COLORS = {"uniform": "#2b6cb0", "vega": "#c05621", "wing": "#2f855a"}


def style(ax: plt.Axes) -> None:
    ax.grid(True, color="#e2e8f0", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def plot_surface(grid: pd.DataFrame, quotes: pd.DataFrame, path: Path) -> None:
    pivot = grid.pivot_table(
        index="log_moneyness", columns="time", values="implied_vol"
    )
    moneyness = pivot.index.to_numpy()
    times = pivot.columns.to_numpy()
    log_m, time_grid = np.meshgrid(moneyness, times)

    fig = plt.figure(figsize=(7.8, 5.4))
    ax = fig.add_subplot(111, projection="3d")
    surface = ax.plot_surface(
        time_grid,
        log_m,
        pivot.to_numpy().T,
        cmap="viridis",
        linewidth=0,
        antialiased=True,
        shade=False,
        alpha=0.9,
    )
    ax.scatter(
        quotes["time"],
        quotes["log_moneyness"],
        quotes["iv_market"],
        s=9,
        c="#111111",
        depthshade=False,
        label="market quotes",
    )
    ax.set_xlabel("time to expiry (years)", labelpad=8)
    ax.set_ylabel("log-moneyness", labelpad=8)
    ax.set_zlabel("implied vol", labelpad=8)
    ax.view_init(elev=26, azim=-118)
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    fig.colorbar(surface, ax=ax, shrink=0.55, pad=0.08, label="implied vol")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_residuals(quotes: pd.DataFrame, path: Path) -> None:
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.4, 4.2))

    for index, (time, group) in enumerate(quotes.groupby("time")):
        group = group.sort_values("log_moneyness")
        left.plot(
            group["log_moneyness"],
            (group["iv_fit"] - group["iv_market"]) * 100.0,
            marker="o",
            markersize=3,
            linewidth=1.2,
            color=PALETTE[index % len(PALETTE)],
            label=f"T={time:g}",
        )
    left.axhline(0.0, color="#4a5568", linewidth=0.8, linestyle="--")
    left.set_xlabel("log-moneyness")
    left.set_ylabel("fit - market (vol points)")
    left.set_title("Residuals by expiry")
    left.legend(frameon=False, fontsize=8, ncol=2)
    style(left)

    for index, (time, group) in enumerate(quotes.groupby("time")):
        group = group.sort_values("log_moneyness")
        color = PALETTE[index % len(PALETTE)]
        right.plot(
            group["log_moneyness"],
            group["iv_fit"] * 100.0,
            linewidth=1.4,
            color=color,
            label=f"T={time:g}",
        )
        right.scatter(
            group["log_moneyness"],
            group["iv_market"] * 100.0,
            s=10,
            color=color,
            alpha=0.55,
        )
    right.set_xlabel("log-moneyness")
    right.set_ylabel("implied vol (%)")
    right.set_title("Fitted smiles (line) vs market (dots)")
    right.legend(frameon=False, fontsize=8, ncol=2)
    style(right)

    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_sweep(weights: pd.DataFrame, sparse: pd.DataFrame, path: Path) -> None:
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.4, 4.2))

    for scheme, group in weights.groupby("scheme"):
        group = group.sort_values("noise")
        left.plot(
            group["noise"] * 100.0,
            group["rmse_true_oos"] * 100.0,
            marker="o",
            markersize=4,
            linewidth=1.6,
            color=SCHEME_COLORS.get(scheme, "#4a5568"),
            label=scheme,
        )

    left.set_xlabel("quote noise (vol points)")
    left.set_ylabel("out-of-sample RMSE (vol points)")
    left.set_title("Surface recovery on held-out strikes")
    left.legend(frameon=False, fontsize=8, title="weights")
    style(left)

    for index, (noise, group) in enumerate(sparse.groupby("noise")):
        group = group.sort_values("n_quotes")
        right.plot(
            group["n_quotes"],
            group["rmse_true_oos"] * 100.0,
            marker="o",
            markersize=4,
            linewidth=1.6,
            color=PALETTE[index % len(PALETTE)],
            label=f"noise {noise * 100:g} vol pts",
        )

    right.set_xlabel("quotes per expiry")
    right.set_ylabel("median out-of-sample RMSE (vol points)")
    right.set_title("How many quotes does a slice need?")
    right.set_xticks(sorted(sparse["n_quotes"].unique()))
    right.legend(frameon=False, fontsize=7)
    style(right)

    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="outputs")
    parser.add_argument("--out-dir", default="outputs")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    grid = pd.read_csv(data_dir / "surface_grid.csv")
    quotes = pd.read_csv(data_dir / "market_vs_model.csv")
    weights = pd.read_csv(data_dir / "sweep_weights.csv")
    sparse = pd.read_csv(data_dir / "sweep_sparse.csv")

    plot_surface(grid, quotes, out_dir / "fig_surface.png")
    plot_residuals(quotes, out_dir / "fig_residuals.png")
    plot_sweep(weights, sparse, out_dir / "fig_sweep.png")

    print(f"figures written to {out_dir.resolve()}")
    print("\nOut-of-sample RMSE (vol points) by noise level:")
    table = weights.pivot(index="noise", columns="scheme", values="rmse_true_oos")
    print((table * 100.0).round(3).to_string())
    print("\nButterfly violation rate (%) by weight scheme:")
    butterfly = weights.pivot(index="noise", columns="scheme", values="butterfly_rate")
    print((butterfly * 100.0).round(3).to_string())
    print("\nButterfly violation rate (%) by quote count:")
    sparse_table = sparse.pivot(
        index="n_quotes", columns="noise", values="butterfly_rate"
    )
    print((sparse_table * 100.0).round(3).to_string())
    print("\nOut-of-sample RMSE (vol points) by quote count:")
    sparse_rmse = sparse.pivot(
        index="n_quotes", columns="noise", values="rmse_true_oos"
    )
    print((sparse_rmse * 100.0).round(3).to_string())


if __name__ == "__main__":
    main()