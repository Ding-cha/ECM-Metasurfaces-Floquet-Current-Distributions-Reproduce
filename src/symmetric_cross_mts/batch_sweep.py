from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .ecm import EcmConfig, sweep
from .model import SymmetricCrossParams, ensure_output_dir
from .run_ecm import _restored_defaults, db, phase_deg


def _float_list(text: str) -> list[float]:
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise ValueError("Parameter list must contain at least one value.")
    return values


def _inclusive_range(start: float, stop: float, step: float) -> list[float]:
    if step <= 0.0:
        raise ValueError("Range step must be positive.")
    count = int(np.floor((stop - start) / step + 0.5)) + 1
    values = [start + idx * step for idx in range(count)]
    return [value for value in values if value <= stop + 1e-9]


def _case_name(params: SymmetricCrossParams) -> str:
    return (
        f"lx{params.lx_mm:g}_ly{params.ly_mm:g}_"
        f"w{params.w_mm:g}_h{params.h_mm:g}_eps{params.epsilon_r:g}"
    ).replace(".", "p")


def _write_case_csv(path: Path, freq_ghz: np.ndarray, s11: np.ndarray, s21: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["freq_ghz", "s11_db", "s21_db", "s11_phase_deg", "s21_phase_deg"])
        for row in zip(freq_ghz, db(s11), db(s21), phase_deg(s11), phase_deg(s21)):
            writer.writerow([f"{item:.8g}" for item in row])


def _features(freq_ghz: np.ndarray, s11: np.ndarray, s21: np.ndarray) -> dict[str, float]:
    s11_db = db(s11)
    s21_db = db(s21)
    s11_dip_idx = int(np.argmin(s11_db))
    s11_peak_idx = int(np.argmax(s11_db))
    s21_notch_idx = int(np.argmin(s21_db))
    return {
        "s21_notch_freq_ghz": float(freq_ghz[s21_notch_idx]),
        "s21_notch_db": float(s21_db[s21_notch_idx]),
        "s11_dip_freq_ghz": float(freq_ghz[s11_dip_idx]),
        "s11_dip_db": float(s11_db[s11_dip_idx]),
        "s11_peak_freq_ghz": float(freq_ghz[s11_peak_idx]),
        "s11_peak_db": float(s11_db[s11_peak_idx]),
    }


def _write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "case_id",
        "lx_mm",
        "ly_mm",
        "w_mm",
        "px_mm",
        "py_mm",
        "h_mm",
        "epsilon_r",
        "mode_order",
        "theta_deg",
        "phi_deg",
        "fourier_samples",
        "current_scale",
        "t00_override",
        "ttr_variant",
        "kz_branch",
        "width_profile",
        "high_order_substrate",
        "epsilon_eff_scale",
        "arm_length_correction_mm",
        "s21_notch_freq_ghz",
        "s21_notch_db",
        "s11_dip_freq_ghz",
        "s11_dip_db",
        "s11_peak_freq_ghz",
        "s11_peak_db",
        "case_csv",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_batch_config(path: Path, args: argparse.Namespace, lx_values: list[float], w_values: list[float], h_values: list[float]) -> None:
    payload = {
        "lx_values_mm": lx_values,
        "ly_equals_lx": args.ly_equals_lx,
        "w_values_mm": w_values,
        "h_values_mm": h_values,
        "px_mm": args.px_mm,
        "py_mm": args.py_mm,
        "epsilon_r": args.epsilon_r,
        "frequency_ghz": {
            "start": args.freq_start_ghz,
            "stop": args.freq_stop_ghz,
            "points": args.freq_points,
        },
        "restored_profile": args.restored_profile,
    }
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)


def run_batch(args: argparse.Namespace) -> Path:
    lx_values = _inclusive_range(args.lx_start_mm, args.lx_stop_mm, args.lx_step_mm)
    w_values = _float_list(args.w_mm_list)
    h_values = _float_list(args.h_mm_list)
    freq_ghz = np.linspace(args.freq_start_ghz, args.freq_stop_ghz, args.freq_points)

    output_dir = ensure_output_dir(args.output_dir)
    case_dir = ensure_output_dir(output_dir / "cases")
    summary_rows: list[dict[str, object]] = []

    for lx_mm in lx_values:
        ly_mm = lx_mm if args.ly_equals_lx else args.ly_mm
        for w_mm in w_values:
            for h_mm in h_values:
                params = SymmetricCrossParams(
                    lx_mm=lx_mm,
                    ly_mm=ly_mm,
                    w_mm=w_mm,
                    px_mm=args.px_mm,
                    py_mm=args.py_mm,
                    h_mm=h_mm,
                    epsilon_r=args.epsilon_r,
                )
                cfg = _restored_defaults(params, args.restored_profile)
                if args.mode_order is not None:
                    cfg = EcmConfig(**{**asdict(cfg), "mode_order": args.mode_order})
                if args.fourier_samples is not None:
                    cfg = EcmConfig(**{**asdict(cfg), "fourier_samples": args.fourier_samples})

                s11, s21 = sweep(freq_ghz, params, cfg)
                case_id = _case_name(params)
                case_path = case_dir / f"{case_id}.csv"
                _write_case_csv(case_path, freq_ghz, s11, s21)

                summary_rows.append(
                    {
                        "case_id": case_id,
                        **params.as_dict(),
                        **asdict(cfg),
                        **_features(freq_ghz, s11, s21),
                        "case_csv": str(case_path.relative_to(output_dir)),
                    }
                )

    _write_summary(output_dir / "summary.csv", summary_rows)
    _write_batch_config(output_dir / "batch_config.json", args, lx_values, w_values, h_values)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch sweep symmetric-cross MTS parameters.")
    parser.add_argument("--output-dir", default="outputs/batch_lx_w_h_test")
    parser.add_argument("--lx-start-mm", type=float, default=20.0)
    parser.add_argument("--lx-stop-mm", type=float, default=25.0)
    parser.add_argument("--lx-step-mm", type=float, default=1.0)
    parser.add_argument("--ly-equals-lx", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--ly-mm", type=float, default=20.0, help="Used only when --no-ly-equals-lx is set.")
    parser.add_argument("--w-mm-list", default="3,4,5")
    parser.add_argument("--h-mm-list", default="1,2,3")
    parser.add_argument("--px-mm", type=float, default=30.0)
    parser.add_argument("--py-mm", type=float, default=30.0)
    parser.add_argument("--epsilon-r", type=float, default=4.4)
    parser.add_argument("--freq-start-ghz", type=float, default=3.0)
    parser.add_argument("--freq-stop-ghz", type=float, default=9.0)
    parser.add_argument("--freq-points", type=int, default=301)
    parser.add_argument("--restored-profile", choices=["auto", "thin", "substrate"], default="auto")
    parser.add_argument("--mode-order", type=int, default=None)
    parser.add_argument("--fourier-samples", type=int, default=None)
    args = parser.parse_args()

    if args.freq_points < 2:
        raise ValueError("--freq-points must be >= 2")
    if args.freq_stop_ghz <= args.freq_start_ghz:
        raise ValueError("--freq-stop-ghz must be larger than --freq-start-ghz")

    output_dir = run_batch(args)
    print(f"Wrote batch sweep outputs to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
