from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fourier import current_fourier
from .model import SymmetricCrossParams

EPS0 = 8.8541878128e-12
MU0 = 4.0e-7 * np.pi
C0 = 1.0 / np.sqrt(EPS0 * MU0)


@dataclass(frozen=True)
class EcmConfig:
    mode_order: int = 7
    theta_deg: float = 0.0
    phi_deg: float = 90.0
    fourier_samples: int = 1201
    current_scale: float = 1.0


def _sqrt_branch(value: complex) -> complex:
    root = np.sqrt(value + 0j)
    if root.imag < 0:
        root = -root
    return root


def _modal_quantities(freq_hz: float, params: SymmetricCrossParams, m: int, n: int, cfg: EcmConfig):
    omega = 2.0 * np.pi * freq_hz
    k0 = omega / C0
    theta = np.deg2rad(cfg.theta_deg)
    phi = np.deg2rad(cfg.phi_deg)

    kx = k0 * np.sin(theta) * np.cos(phi) + 2.0 * np.pi * m / params.px_m
    ky = k0 * np.sin(theta) * np.sin(phi) + 2.0 * np.pi * n / params.py_m

    kz_air = _sqrt_branch(k0**2 - kx**2 - ky**2)
    kz_diel = _sqrt_branch(k0**2 * params.epsilon_r - kx**2 - ky**2)

    y_te_plus = kz_air / (omega * MU0)
    y_tm_plus = omega * EPS0 / kz_air
    y_te_minus = kz_diel / (omega * MU0)
    y_tm_minus = omega * EPS0 * params.epsilon_r / kz_diel
    return kx, ky, kz_air, kz_diel, y_te_plus, y_tm_plus, y_te_minus, y_tm_minus


def _loaded_admittance(y_plus: complex, y_minus: complex, kz_diel: complex, h_m: float) -> complex:
    zc = 1.0 / y_minus
    z_load = 1.0 / y_plus
    tan_term = np.tan(kz_diel * h_m)
    z_in = zc * (z_load + 1j * zc * tan_term) / (zc + 1j * z_load * tan_term)
    return 1.0 / z_in


def _turn_ratios(
    kx: complex,
    ky: complex,
    jx: complex,
    jy: complex,
    params: SymmetricCrossParams,
    cfg: EcmConfig,
) -> tuple[complex, complex]:
    krho = np.sqrt(kx**2 + ky**2 + 0j)
    factor = (4.0 * np.pi**2) ** 2 / (params.px_m * params.py_m)

    if abs(krho) < 1e-12:
        # Normal-incidence TE00 limit for x-polarized excitation.
        j_te = jx
        j_tm = 0.0j
        return cfg.current_scale**2 * factor * j_te**2, cfg.current_scale**2 * factor * j_tm**2
    else:
        ft_te = jx * ky - jy * kx
        ft_tm = jx * kx + jy * ky
        t_te = factor * ft_te**2 * ky**2 / (kx**2 + ky**2)
        t_tm = factor * ft_tm**2 * kx**2 / (kx**2 + ky**2)
        return cfg.current_scale**2 * t_te, cfg.current_scale**2 * t_tm


def solve_frequency(freq_hz: float, params: SymmetricCrossParams, cfg: EcmConfig) -> tuple[complex, complex]:
    """Return approximate TE00 S11 and S21 for the symmetric-cross MTS."""

    z_ab = 0.0j
    t00_te = None
    b00_te = None
    z00_te_plus = None
    kz00_diel = None
    y00_te_plus = None
    y00_te_minus = None

    for m in range(-cfg.mode_order, cfg.mode_order + 1):
        for n in range(-cfg.mode_order, cfg.mode_order + 1):
            kx, ky, kz_air, kz_diel, y_te_p, y_tm_p, y_te_m, y_tm_m = _modal_quantities(
                freq_hz, params, m, n, cfg
            )
            b_te = _loaded_admittance(y_te_p, y_te_m, kz_diel, params.h_m)
            b_tm = _loaded_admittance(y_tm_p, y_tm_m, kz_diel, params.h_m)
            y_te = y_te_p + b_te
            y_tm = y_tm_p + b_tm
            jx, jy = current_fourier(kx, ky, params, cfg.fourier_samples)
            t_te, t_tm = _turn_ratios(kx, ky, jx, jy, params, cfg)

            if m == 0 and n == 0:
                t00_te = t_te
                b00_te = b_te
                z00_te_plus = 1.0 / y_te_p
                kz00_diel = kz_diel
                y00_te_plus = y_te_p
                y00_te_minus = y_te_m
            else:
                z_ab += t_te / y_te

            z_ab += t_tm / y_tm

    if t00_te is None or b00_te is None or z00_te_plus is None:
        raise RuntimeError("TE00 mode was not evaluated.")

    z_in = 1.0 / (1.0 / (z_ab / t00_te) + b00_te)
    s11 = (z_in - z00_te_plus) / (z_in + z00_te_plus)

    zc00 = 1.0 / y00_te_minus
    zload00 = 1.0 / y00_te_plus
    gamma = (zload00 - zc00) / (zload00 + zc00)
    phase = kz00_diel * params.h_m
    s21 = (2.0 * z_in / (z_in + z00_te_plus)) * (1.0 + gamma) / (
        np.exp(1j * phase) + gamma * np.exp(-1j * phase)
    )
    return s11, s21


def sweep(freq_ghz: np.ndarray, params: SymmetricCrossParams, cfg: EcmConfig) -> tuple[np.ndarray, np.ndarray]:
    s11 = np.empty(freq_ghz.shape, dtype=complex)
    s21 = np.empty(freq_ghz.shape, dtype=complex)
    for idx, f_ghz in enumerate(freq_ghz):
        s11[idx], s21[idx] = solve_frequency(f_ghz * 1e9, params, cfg)
    return s11, s21
