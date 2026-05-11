from __future__ import annotations

import itertools

import numpy as np

from .ecm import EcmConfig, sweep
from .model import SymmetricCrossParams
from .run_ecm import db


TARGET_FREQ = np.array([3.0, 5.0, 5.8, 6.0, 7.0, 7.8, 9.0])
TARGET_S11_DB = np.array([-11.0, -3.0, -0.1, -1.0, -14.0, -24.0, -14.0])
TARGET_S21_DB = np.array([-0.05, -2.5, -25.0, -8.0, -0.2, -0.05, -0.15])


def _score(s11_db: np.ndarray, s21_db: np.ndarray) -> float:
    return float(np.sqrt(np.mean((s11_db - TARGET_S11_DB) ** 2 + 0.4 * (s21_db - TARGET_S21_DB) ** 2)))


def main():
    params = SymmetricCrossParams()
    variants = ["paper", "polarized", "no_extra_k"]
    branches = ["principal", "decay"]
    substrate_options = [True, False]
    t00_values = np.r_[np.logspace(-3, 2, 41), [None]]
    rows = []
    for variant, branch, high_order_substrate, t00 in itertools.product(
        variants, branches, substrate_options, t00_values
    ):
        cfg = EcmConfig(
            mode_order=7,
            fourier_samples=501,
            ttr_variant=variant,
            kz_branch=branch,
            high_order_substrate=high_order_substrate,
            t00_override=None if t00 is None else float(t00),
        )
        try:
            s11, s21 = sweep(TARGET_FREQ, params, cfg)
        except (FloatingPointError, ZeroDivisionError, ValueError):
            continue
        s11_db = db(s11)
        s21_db = db(s21)
        rows.append((_score(s11_db, s21_db), variant, branch, high_order_substrate, t00, s11_db, s21_db))

    rows.sort(key=lambda item: item[0])
    print("Top ECM restoration candidates")
    print("score | variant | kz_branch | high-order substrate | T00 override | S11 dB @ [3,5,5.8,6,7,7.8,9] | S21 dB")
    for score, variant, branch, high_order_substrate, t00, s11_db, s21_db in rows[:12]:
        t00_text = "raw" if t00 is None else f"{t00:.4g}"
        print(
            f"{score:6.2f} | {variant:10s} | {branch:9s} | {str(high_order_substrate):>20s} | {t00_text:>10s} | "
            f"{np.round(s11_db, 1)} | {np.round(s21_db, 1)}"
        )


if __name__ == "__main__":
    main()
