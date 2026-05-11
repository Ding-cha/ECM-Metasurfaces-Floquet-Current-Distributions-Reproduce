from __future__ import annotations

import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np

from .calibrated import calibrated_sweep, target_anchor_data
from .ecm import EcmConfig, sweep
from .model import SymmetricCrossParams, ensure_output_dir


def db(value: np.ndarray) -> np.ndarray:
    return 20.0 * np.log10(np.maximum(np.abs(value), 1e-12))


def phase_deg(value: np.ndarray) -> np.ndarray:
    return np.rad2deg(np.unwrap(np.angle(value)))


def write_csv(path, freq_ghz, s11, s21):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["freq_ghz", "s11_db", "s21_db", "s11_phase_deg", "s21_phase_deg"])
        for row in zip(freq_ghz, db(s11), db(s21), phase_deg(s11), phase_deg(s21)):
            writer.writerow([f"{item:.8g}" for item in row])


def plot_sparameters(output_dir, freq_ghz, s11, s21, *, show_anchors: bool = False, title_suffix: str = ""):
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)
    axes[0].plot(freq_ghz, db(s11), label="$S_{11}$ ECM")
    axes[0].plot(freq_ghz, db(s21), label="$S_{21}$ ECM")
    if show_anchors:
        anchors = target_anchor_data()
        axes[0].scatter(*anchors["s11_db"], marker="o", s=18, label="$S_{11}$ Fig. 15 anchors")
        axes[0].scatter(*anchors["s21_db"], marker="s", s=18, label="$S_{21}$ Fig. 15 anchors")
    if title_suffix:
        axes[0].set_title(title_suffix)
    axes[0].set_ylabel("Magnitude (dB)")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend()

    axes[1].plot(freq_ghz, phase_deg(s11), label="$S_{11}$ ECM")
    axes[1].plot(freq_ghz, phase_deg(s21), label="$S_{21}$ ECM")
    if show_anchors:
        anchors = target_anchor_data()
        axes[1].scatter(*anchors["s11_phase_deg"], marker="o", s=18, label="$S_{11}$ Fig. 15 anchors")
        axes[1].scatter(*anchors["s21_phase_deg"], marker="s", s=18, label="$S_{21}$ Fig. 15 anchors")
    axes[1].set_xlabel("Frequency (GHz)")
    axes[1].set_ylabel("Phase (deg)")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend()
    fig.savefig(output_dir / "s_parameters.png", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Run the symmetric-cross MTS sweep.")
    parser.add_argument(
        "--mode",
        choices=["restored", "calibrated", "raw"],
        default="restored",
        help="Select restored analytical candidate, Fig. 15 calibration, or original raw ECM.",
    )
    parser.add_argument("--raw", action="store_true", help="Compatibility alias for --mode raw.")
    parser.add_argument("--calibrated", action="store_true", help="Compatibility alias for --mode calibrated.")
    parser.add_argument("--anchors", action="store_true", help="Overlay the hand-digitized Fig. 15 anchor points.")
    args = parser.parse_args()

    mode = args.mode
    if args.raw:
        mode = "raw"
    if args.calibrated:
        mode = "calibrated"

    # 601 个频点可以解析出 S21 在 5.8 GHz 附近的窄陷波。
    freq_ghz = np.linspace(3.0, 9.0, 601)
    if mode == "raw":
        # raw：尽量直接按论文公式实现，用来展示最初差异来自哪里。
        params = SymmetricCrossParams()
        cfg = EcmConfig(mode_order=7, fourier_samples=1201)
        s11, s21 = sweep(freq_ghz, params, cfg)
        title = "Raw analytical ECM"
    elif mode == "calibrated":
        # calibrated：用作者 Fig. 15 手工锚点插值，只用于图形复现对照。
        s11, s21 = calibrated_sweep(freq_ghz)
        title = "Fig. 15-calibrated reproduction"
    else:
        # restored：当前最接近作者曲线的解析候选，不依赖锚点插值。
        params = SymmetricCrossParams()
        cfg = EcmConfig(
            mode_order=9,
            fourier_samples=501,
            ttr_variant="polarized",
            kz_branch="decay",
            high_order_substrate=False,
        )
        s11, s21 = sweep(freq_ghz, params, cfg)
        title = "Restored ECM candidate"

    output_dir = ensure_output_dir()
    write_csv(output_dir / "s_parameters.csv", freq_ghz, s11, s21)
    plot_sparameters(output_dir, freq_ghz, s11, s21, show_anchors=args.anchors, title_suffix=title)
    print(f"Wrote ECM outputs to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
