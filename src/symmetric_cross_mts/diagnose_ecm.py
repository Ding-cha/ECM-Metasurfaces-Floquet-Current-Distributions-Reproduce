from __future__ import annotations

import numpy as np

from .ecm import (
    EcmConfig,
    _loaded_admittance,
    _modal_quantities,
    _turn_ratios,
    solve_frequency,
)
from .fourier import current_fourier
from .model import SymmetricCrossParams
from .run_ecm import db, phase_deg


def _collect_terms(freq_ghz: float, params: SymmetricCrossParams, cfg: EcmConfig):
    freq_hz = freq_ghz * 1e9
    z_ab = 0.0j
    t00_te = None
    b00_te = None
    z00_te_plus = None

    for m in range(-cfg.mode_order, cfg.mode_order + 1):
        for n in range(-cfg.mode_order, cfg.mode_order + 1):
            kx, ky, _kz_air, kz_diel, y_te_p, y_tm_p, y_te_m, y_tm_m = _modal_quantities(
                freq_hz, params, m, n, cfg
            )
            b_te = _loaded_admittance(y_te_p, y_te_m, kz_diel, params.h_m)
            b_tm = _loaded_admittance(y_tm_p, y_tm_m, kz_diel, params.h_m)
            jx, jy = current_fourier(kx, ky, params, cfg.fourier_samples)
            t_te, t_tm = _turn_ratios(kx, ky, jx, jy, params, cfg)

            if m == 0 and n == 0:
                t00_te = t_te
                b00_te = b_te
                z00_te_plus = 1.0 / y_te_p
            else:
                z_ab += t_te / (y_te_p + b_te)
            z_ab += t_tm / (y_tm_p + b_tm)

    if t00_te is None or b00_te is None or z00_te_plus is None:
        raise RuntimeError("TE00 mode was not evaluated.")
    z_in_bare_slab = 1.0 / b00_te
    s11_bare_slab = (z_in_bare_slab - z00_te_plus) / (z_in_bare_slab + z00_te_plus)
    return z_ab, t00_te, b00_te, z00_te_plus, s11_bare_slab


def main():
    params = SymmetricCrossParams()
    cfg = EcmConfig(mode_order=7, fourier_samples=1201)
    freq_ghz = 5.8

    z_ab, t00_te, b00_te, z00_te_plus, s11_bare_slab = _collect_terms(freq_ghz, params, cfg)
    s11, s21 = solve_frequency(freq_ghz * 1e9, params, cfg)

    print(f"Diagnostic frequency: {freq_ghz:.3f} GHz")
    print(f"T00_TE = {t00_te:.6e}")
    print(f"Z_AB = {z_ab:.6e} ohm")
    print(f"Z_AB / T00_TE = {z_ab / t00_te:.6e} ohm")
    print(f"1 / B00_TE = {1.0 / b00_te:.6e} ohm")
    print(f"Z00_TE+ = {z00_te_plus:.6e} ohm")
    print()
    print(f"Current code S11 = {db(np.array([s11]))[0]:.2f} dB, {phase_deg(np.array([s11]))[0]:.2f} deg")
    print(f"Current code S21 = {db(np.array([s21]))[0]:.2f} dB, {phase_deg(np.array([s21]))[0]:.2f} deg")
    print(
        f"Bare dielectric S11 = {db(np.array([s11_bare_slab]))[0]:.2f} dB, "
        f"{phase_deg(np.array([s11_bare_slab]))[0]:.2f} deg"
    )
    print()
    print("Interpretation:")
    print("- T00_TE is extremely small because it is taken directly from the area-normalized Fourier integral.")
    print("- As a result, Z_AB / T00_TE becomes very large and the sheet branch contributes almost nothing.")
    print("- The computed S-parameters are therefore close to the bare dielectric slab response, not Fig. 15.")
    print("- Matching the paper requires the missing incident-current/TTR normalization and the exact TE00 limit.")


if __name__ == "__main__":
    main()
