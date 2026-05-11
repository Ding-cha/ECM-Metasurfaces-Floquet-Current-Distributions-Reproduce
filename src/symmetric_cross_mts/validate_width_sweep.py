from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .ecm import EcmConfig, sweep
from .model import SymmetricCrossParams, ensure_output_dir
from .run_ecm import db


WIDTHS_MM = (3.0, 5.0, 7.0)
FREQ_GHZ = np.linspace(3.0, 9.0, 301)


RESTORED_WIDTH_CORRECTED_CFG = EcmConfig(
    mode_order=9,
    fourier_samples=501,
    ttr_variant="polarized",
    kz_branch="decay",
    high_order_substrate=False,
    width_profile="cosine",
)

RESTORED_UNIFORM_CFG = EcmConfig(
    mode_order=9,
    fourier_samples=501,
    ttr_variant="polarized",
    kz_branch="decay",
    high_order_substrate=False,
)

RAW_CFG = EcmConfig(mode_order=7, fourier_samples=701)


def _metrics(freq_ghz: np.ndarray, s11: np.ndarray) -> dict[str, float]:
    s11_db = db(s11)
    peak_idx = int(np.argmax(s11_db))
    high_mask = freq_ghz > 6.5
    high_freq = freq_ghz[high_mask]
    high_s11 = s11_db[high_mask]
    valley_idx = int(np.argmin(high_s11))
    return {
        "peak_freq_ghz": float(freq_ghz[peak_idx]),
        "peak_s11_db": float(s11_db[peak_idx]),
        "valley_freq_ghz": float(high_freq[valley_idx]),
        "valley_s11_db": float(high_s11[valley_idx]),
    }


def _run_case(width_mm: float, cfg: EcmConfig) -> tuple[np.ndarray, np.ndarray]:
    # Fig. 22 只改变金属臂宽 w，保持 lx=ly=20 mm，其余参数沿用 Table II。
    params = replace(SymmetricCrossParams(), w_mm=width_mm)
    s11, _s21 = sweep(FREQ_GHZ, params, cfg)
    return s11, db(s11)


def _write_summary(output_dir: Path, rows: list[dict[str, float | str]]) -> None:
    path = output_dir / "width_sweep_summary.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mode",
                "w_mm",
                "peak_freq_ghz",
                "peak_s11_db",
                "valley_freq_ghz",
                "valley_s11_db",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def _plot(output_dir: Path, curves: dict[tuple[str, float], np.ndarray]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True, constrained_layout=True)
    for ax, width_mm in zip(axes, WIDTHS_MM):
        ax.plot(FREQ_GHZ, curves[("restored-width-corrected", width_mm)], label="restored width-corrected", linewidth=2.0)
        ax.plot(FREQ_GHZ, curves[("restored-uniform", width_mm)], "--", label="restored uniform", linewidth=1.4)
        ax.plot(FREQ_GHZ, curves[("raw", width_mm)], "--", label="raw", linewidth=1.4)
        ax.set_title(f"w = {width_mm:g} mm")
        ax.set_xlabel("Freq (GHz)")
        ax.grid(True, alpha=0.28)
        ax.set_ylim(-30, 1)
    axes[0].set_ylabel("$S_{11}$ (dB)")
    axes[0].legend()
    fig.suptitle("Fig. 22-style width sweep validation")
    fig.savefig(output_dir / "width_sweep_s11.png", dpi=220)
    plt.close(fig)


def main() -> None:
    output_dir = ensure_output_dir()
    curves: dict[tuple[str, float], np.ndarray] = {}
    rows: list[dict[str, float | str]] = []

    for width_mm in WIDTHS_MM:
        for mode, cfg in (
            ("restored-width-corrected", RESTORED_WIDTH_CORRECTED_CFG),
            ("restored-uniform", RESTORED_UNIFORM_CFG),
            ("raw", RAW_CFG),
        ):
            s11, s11_db = _run_case(width_mm, cfg)
            curves[(mode, width_mm)] = s11_db
            metric = _metrics(FREQ_GHZ, s11)
            rows.append({"mode": mode, "w_mm": width_mm, **metric})

    _plot(output_dir, curves)
    _write_summary(output_dir, rows)

    print(f"Wrote width sweep outputs to {output_dir.resolve()}")
    for row in rows:
        print(
            f"{row['mode']:8s} w={row['w_mm']:.0f} mm: "
            f"peak {row['peak_freq_ghz']:.2f} GHz / {row['peak_s11_db']:.2f} dB, "
            f"valley {row['valley_freq_ghz']:.2f} GHz / {row['valley_s11_db']:.2f} dB"
        )


if __name__ == "__main__":
    main()
