from __future__ import annotations

import numpy as np

from .model import SymmetricCrossParams, horizontal_current_x, vertical_current_y


def _trapz_complex(values: np.ndarray, grid_m: np.ndarray) -> complex:
    trapz = getattr(np, "trapezoid", np.trapz)
    return complex(trapz(values, grid_m))


def _uniform_integral(k: complex, half_width_m: float) -> complex:
    if abs(k) < 1e-14:
        return 2.0 * half_width_m
    return 2.0 * np.sin(k * half_width_m) / k


def current_fourier(
    kx: complex,
    ky: complex,
    params: SymmetricCrossParams,
    samples: int = 1201,
) -> tuple[complex, complex]:
    """2-D Fourier transform of the analytical current model.

    The separable strip geometry keeps this inexpensive while preserving the
    paper's published current profile. Coordinates inside exponentials are in m;
    the cosine current model is evaluated in mm.
    """

    x_mm = np.linspace(-params.lx_mm / 2.0, params.lx_mm / 2.0, samples)
    x_m = x_mm * 1e-3
    y_mm = np.linspace(-params.ly_mm / 2.0, params.ly_mm / 2.0, samples)
    y_m = y_mm * 1e-3

    ix = horizontal_current_x(x_mm, params)
    iy = vertical_current_y(y_mm, params)

    int_x_current = _trapz_complex(ix * np.exp(1j * kx * x_m), x_m)
    int_y_width = _uniform_integral(ky, params.w_m / 2.0)
    jx_tilde = int_x_current * int_y_width / (4.0 * np.pi**2)

    int_x_width = _uniform_integral(kx, params.w_m / 2.0)
    int_y_current = _trapz_complex(iy * np.exp(1j * ky * y_m), y_m)
    jy_tilde = int_x_width * int_y_current / (4.0 * np.pi**2)
    return jx_tilde, jy_tilde
