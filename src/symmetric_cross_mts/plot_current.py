from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from .model import SymmetricCrossParams, current_grid, ensure_output_dir


def plot_geometry(params: SymmetricCrossParams, output_dir):
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    ax.add_patch(
        Rectangle(
            (-params.px_mm / 2, -params.py_mm / 2),
            params.px_mm,
            params.py_mm,
            facecolor="#e9e9e9",
            edgecolor="#444444",
            linewidth=1.2,
        )
    )
    ax.add_patch(
        Rectangle(
            (-params.lx_mm / 2, -params.w_mm / 2),
            params.lx_mm,
            params.w_mm,
            facecolor="#f6a21a",
            edgecolor="#9a5a00",
        )
    )
    ax.add_patch(
        Rectangle(
            (-params.w_mm / 2, -params.ly_mm / 2),
            params.w_mm,
            params.ly_mm,
            facecolor="#f6a21a",
            edgecolor="#9a5a00",
        )
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-params.px_mm / 2 - 2, params.px_mm / 2 + 2)
    ax.set_ylim(-params.py_mm / 2 - 2, params.py_mm / 2 + 2)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Symmetric cross MTS unit cell")
    ax.grid(True, alpha=0.25)
    fig.savefig(output_dir / "geometry.png", dpi=220)
    plt.close(fig)


def plot_current(params: SymmetricCrossParams, output_dir):
    xx, yy, jx, jy = current_grid(params)
    mag = np.sqrt(jx**2 + jy**2)

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    im = ax.pcolormesh(xx, yy, mag, shading="auto", cmap="viridis")
    step = 28
    ax.quiver(
        xx[::step, ::step],
        yy[::step, ::step],
        jx[::step, ::step],
        jy[::step, ::step],
        color="#f28e2b",
        scale=15,
        width=0.004,
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Analytical surface current model")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("a.u.")
    fig.savefig(output_dir / "current_distribution.png", dpi=220)
    plt.close(fig)


def main():
    params = SymmetricCrossParams()
    output_dir = ensure_output_dir()
    plot_geometry(params, output_dir)
    plot_current(params, output_dir)
    print(f"Wrote plots to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
