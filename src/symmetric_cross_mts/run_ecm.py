from __future__ import annotations

import csv

import matplotlib.pyplot as plt
import numpy as np

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


def plot_sparameters(output_dir, freq_ghz, s11, s21):
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)
    axes[0].plot(freq_ghz, db(s11), label="$S_{11}$ ECM")
    axes[0].plot(freq_ghz, db(s21), label="$S_{21}$ ECM")
    axes[0].set_ylabel("Magnitude (dB)")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend()

    axes[1].plot(freq_ghz, phase_deg(s11), label="$S_{11}$ ECM")
    axes[1].plot(freq_ghz, phase_deg(s21), label="$S_{21}$ ECM")
    axes[1].set_xlabel("Frequency (GHz)")
    axes[1].set_ylabel("Phase (deg)")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend()
    fig.savefig(output_dir / "s_parameters.png", dpi=220)
    plt.close(fig)


def main():
    params = SymmetricCrossParams()
    cfg = EcmConfig(mode_order=7, fourier_samples=1201)
    freq_ghz = np.linspace(3.0, 9.0, 121)
    s11, s21 = sweep(freq_ghz, params, cfg)

    output_dir = ensure_output_dir()
    write_csv(output_dir / "s_parameters.csv", freq_ghz, s11, s21)
    plot_sparameters(output_dir, freq_ghz, s11, s21)
    print(f"Wrote ECM outputs to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
