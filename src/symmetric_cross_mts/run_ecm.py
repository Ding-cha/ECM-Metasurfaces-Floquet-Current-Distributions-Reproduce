from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .calibrated import calibrated_sweep, target_anchor_data
from .ecm import EcmConfig, sweep
from .model import SymmetricCrossParams, ensure_output_dir


GEOMETRY_ALIASES = {
    "lx": "lx_mm",
    "l_x": "lx_mm",
    "lx_mm": "lx_mm",
    "ly": "ly_mm",
    "l_y": "ly_mm",
    "ly_mm": "ly_mm",
    "w": "w_mm",
    "width": "w_mm",
    "w_mm": "w_mm",
    "px": "px_mm",
    "p_x": "px_mm",
    "px_mm": "px_mm",
    "py": "py_mm",
    "p_y": "py_mm",
    "py_mm": "py_mm",
    "h": "h_mm",
    "h_mm": "h_mm",
    "epsilon_r": "epsilon_r",
    "eps_r": "epsilon_r",
}


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


def write_run_config(
    path,
    params: SymmetricCrossParams | None,
    cfg: EcmConfig | None,
    args: argparse.Namespace,
    mode: str,
) -> None:
    payload = {
        "mode": mode,
        "restored_profile": args.restored_profile,
        "frequency_ghz": {
            "start": args.freq_start_ghz,
            "stop": args.freq_stop_ghz,
            "points": args.freq_points,
        },
        "geometry": None if params is None else params.as_dict(),
        "electrical_thickness_mm_sqrt_eps": None
        if params is None
        else params.h_mm * float(np.sqrt(params.epsilon_r)),
        "ecm_config": None if cfg is None else asdict(cfg),
        "source_params_json": None if args.params_json is None else str(Path(args.params_json).resolve()),
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)


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


def _read_geometry_json(path: str) -> dict[str, float]:
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if "geometry" in payload and isinstance(payload["geometry"], dict):
        payload = payload["geometry"]

    geometry = {}
    for raw_name, value in payload.items():
        name = GEOMETRY_ALIASES.get(raw_name)
        if name is None:
            raise ValueError(f"Unknown geometry parameter in JSON: {raw_name}")
        geometry[name] = float(value)
    return geometry


def _build_params(args: argparse.Namespace) -> SymmetricCrossParams:
    geometry = asdict(SymmetricCrossParams())
    if args.params_json:
        geometry.update(_read_geometry_json(args.params_json))

    for name in ("lx_mm", "ly_mm", "w_mm", "px_mm", "py_mm", "h_mm", "epsilon_r"):
        value = getattr(args, name)
        if value is not None:
            geometry[name] = value

    return SymmetricCrossParams(**geometry)


def _restored_defaults(params: SymmetricCrossParams, profile: str) -> EcmConfig:
    if profile == "auto":
        electrical_thickness = params.h_mm * np.sqrt(params.epsilon_r)
        profile = "substrate" if electrical_thickness >= 4.0 else "thin"

    if profile == "thin":
        return EcmConfig(
            mode_order=9,
            fourier_samples=501,
            ttr_variant="polarized",
            kz_branch="decay",
            high_order_substrate=False,
            width_profile="cosine",
            epsilon_eff_scale=1.0,
        )
    if profile == "substrate":
        epsilon_eff_scale = float(np.clip(1.03 - 0.04 * params.h_mm, 0.88, 0.97))
        return EcmConfig(
            mode_order=3,
            fourier_samples=501,
            ttr_variant="polarized",
            kz_branch="decay",
            high_order_substrate=True,
            width_profile="uniform",
            epsilon_eff_scale=epsilon_eff_scale,
            arm_length_correction_mm=0.5,
        )
    raise ValueError(f"Unknown restored profile: {profile}")


def _build_cfg(args: argparse.Namespace, mode: str, params: SymmetricCrossParams) -> EcmConfig:
    if mode == "raw":
        defaults = EcmConfig(mode_order=7, fourier_samples=1201)
    else:
        defaults = _restored_defaults(params, args.restored_profile)

    if args.high_order_substrate == "on":
        high_order_substrate = True
    elif args.high_order_substrate == "off":
        high_order_substrate = False
    else:
        high_order_substrate = defaults.high_order_substrate

    cfg = EcmConfig(
        mode_order=defaults.mode_order if args.mode_order is None else args.mode_order,
        theta_deg=args.theta_deg,
        phi_deg=args.phi_deg,
        fourier_samples=defaults.fourier_samples if args.fourier_samples is None else args.fourier_samples,
        current_scale=args.current_scale,
        t00_override=args.t00_override,
        ttr_variant=defaults.ttr_variant if args.ttr_variant is None else args.ttr_variant,
        kz_branch=defaults.kz_branch if args.kz_branch is None else args.kz_branch,
        high_order_substrate=high_order_substrate,
        width_profile=defaults.width_profile if args.width_profile is None else args.width_profile,
        epsilon_eff_scale=defaults.epsilon_eff_scale if args.epsilon_eff_scale is None else args.epsilon_eff_scale,
        arm_length_correction_mm=defaults.arm_length_correction_mm
        if args.arm_length_correction_mm is None
        else args.arm_length_correction_mm,
    )
    if cfg.mode_order < 0:
        raise ValueError("--mode-order must be >= 0")
    if cfg.fourier_samples < 2:
        raise ValueError("--fourier-samples must be >= 2")
    if not np.isfinite(cfg.current_scale):
        raise ValueError("--current-scale must be finite")
    if cfg.t00_override is not None and cfg.t00_override == 0.0:
        raise ValueError("--t00-override must be nonzero when supplied")
    if cfg.epsilon_eff_scale <= 0.0 or cfg.epsilon_eff_scale > 1.0:
        raise ValueError("--epsilon-eff-scale must be in the interval (0, 1]")
    if cfg.arm_length_correction_mm < 0.0:
        raise ValueError("--arm-length-correction-mm must be >= 0")
    return cfg


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
    parser.add_argument("--output-dir", default="outputs", help="Directory for PNG/CSV/JSON outputs.")
    parser.add_argument("--params-json", default=None, help="Optional geometry JSON file. CLI geometry values override it.")

    parser.add_argument("--lx-mm", type=float, default=None, help="Horizontal arm length in mm.")
    parser.add_argument("--ly-mm", type=float, default=None, help="Vertical arm length in mm.")
    parser.add_argument("--w-mm", type=float, default=None, help="Arm width in mm.")
    parser.add_argument("--px-mm", type=float, default=None, help="Unit-cell period in x in mm.")
    parser.add_argument("--py-mm", type=float, default=None, help="Unit-cell period in y in mm.")
    parser.add_argument("--h-mm", type=float, default=None, help="Substrate thickness in mm.")
    parser.add_argument("--epsilon-r", type=float, default=None, help="Substrate relative permittivity.")

    parser.add_argument("--freq-start-ghz", type=float, default=3.0, help="Start frequency.")
    parser.add_argument("--freq-stop-ghz", type=float, default=9.0, help="Stop frequency.")
    parser.add_argument("--freq-points", type=int, default=601, help="Number of frequency samples.")

    parser.add_argument("--mode-order", type=int, default=None, help="Floquet mode truncation order.")
    parser.add_argument("--fourier-samples", type=int, default=None, help="1-D integration samples for current FT.")
    parser.add_argument("--theta-deg", type=float, default=0.0, help="Incident polar angle.")
    parser.add_argument("--phi-deg", type=float, default=90.0, help="Incident azimuth angle.")
    parser.add_argument("--current-scale", type=float, default=1.0, help="Optional current/TTR scale factor.")
    parser.add_argument("--t00-override", type=float, default=None, help="Optional TE00 TTR override for experiments.")
    parser.add_argument("--ttr-variant", choices=["paper", "polarized", "no_extra_k"], default=None)
    parser.add_argument("--kz-branch", choices=["principal", "decay"], default=None)
    parser.add_argument("--width-profile", choices=["uniform", "cosine"], default=None)
    parser.add_argument("--high-order-substrate", choices=["auto", "on", "off"], default="auto")
    parser.add_argument(
        "--epsilon-eff-scale",
        type=float,
        default=None,
        help="Use eps_eff = 1 + scale * (epsilon_r - 1) in the substrate modal admittance.",
    )
    parser.add_argument(
        "--arm-length-correction-mm",
        type=float,
        default=None,
        help="Shorten the current-integration arm length without changing the physical period/substrate.",
    )
    parser.add_argument(
        "--restored-profile",
        choices=["auto", "thin", "substrate"],
        default="auto",
        help="Default restored settings. auto uses substrate settings for electrically thicker dielectric slabs.",
    )
    args = parser.parse_args()

    mode = args.mode
    if args.raw:
        mode = "raw"
    if args.calibrated:
        mode = "calibrated"

    if not np.isfinite(args.freq_start_ghz) or not np.isfinite(args.freq_stop_ghz):
        raise ValueError("--freq-start-ghz and --freq-stop-ghz must be finite")
    if args.freq_points < 2:
        raise ValueError("--freq-points must be >= 2")
    if args.freq_stop_ghz <= args.freq_start_ghz:
        raise ValueError("--freq-stop-ghz must be larger than --freq-start-ghz")

    # 高频陷波较窄，默认使用 601 个频点；用户可用 CLI 调整扫频范围和密度。
    freq_ghz = np.linspace(args.freq_start_ghz, args.freq_stop_ghz, args.freq_points)
    params = None
    cfg = None
    if mode == "raw":
        # raw：尽量直接按论文公式实现，用来展示最初差异来自哪里。
        params = _build_params(args)
        cfg = _build_cfg(args, mode, params)
        s11, s21 = sweep(freq_ghz, params, cfg)
        title = "Raw analytical ECM"
    elif mode == "calibrated":
        # calibrated：用作者 Fig. 15 手工锚点插值，只用于图形复现对照。
        s11, s21 = calibrated_sweep(freq_ghz)
        title = "Fig. 15-calibrated reproduction"
    else:
        # restored：当前最接近作者曲线的解析候选，不依赖锚点插值。
        params = _build_params(args)
        cfg = _build_cfg(args, mode, params)
        s11, s21 = sweep(freq_ghz, params, cfg)
        title = (
            "Restored ECM candidate "
            f"(lx={params.lx_mm:g}, ly={params.ly_mm:g}, w={params.w_mm:g} mm)"
        )

    output_dir = ensure_output_dir(args.output_dir)
    write_csv(output_dir / "s_parameters.csv", freq_ghz, s11, s21)
    write_run_config(output_dir / "run_config.json", params, cfg, args, mode)
    plot_sparameters(output_dir, freq_ghz, s11, s21, show_anchors=args.anchors, title_suffix=title)
    print(f"Wrote ECM outputs to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
