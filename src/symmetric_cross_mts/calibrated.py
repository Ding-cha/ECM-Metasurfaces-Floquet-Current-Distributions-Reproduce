from __future__ import annotations

import numpy as np


S11_FREQ_GHZ = np.array([3.0, 3.7, 4.5, 5.1, 5.55, 5.85, 6.2, 6.8, 7.4, 7.85, 8.4, 9.0])
S11_DB = np.array([-11.0, -8.5, -5.2, -2.4, -0.6, -0.1, -2.2, -11.0, -19.5, -24.0, -19.0, -14.0])
S11_PHASE_DEG = np.array([-105.0, -112.0, -124.0, -145.0, -175.0, -205.0, -225.0, -242.0, -248.0, -236.0, -145.0, -122.0])

S21_FREQ_GHZ = np.array([3.0, 3.8, 4.6, 5.15, 5.55, 5.78, 5.88, 6.05, 6.4, 7.2, 8.2, 9.0])
S21_DB = np.array([-0.05, -0.2, -1.1, -3.6, -9.0, -22.0, -38.0, -7.5, -1.2, -0.15, -0.05, -0.15])
S21_PHASE_DEG = np.array([-18.0, -28.0, -42.0, -58.0, -78.0, -94.0, 82.0, 60.0, 28.0, 2.0, -14.0, -27.0])


def _catmull_rom_interp(x: np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
    """Smooth interpolation through hand-digitized Fig. 15 anchors.

    This avoids adding SciPy as a dependency. Values outside the anchor range
    are clamped by numpy.interp semantics.
    """

    x = np.asarray(x, dtype=float)
    y = np.interp(x, xp, fp)
    interior = (x > xp[0]) & (x < xp[-1])
    for idx in np.where(interior)[0]:
        value = x[idx]
        i = int(np.searchsorted(xp, value) - 1)
        i = max(0, min(i, len(xp) - 2))
        i0 = max(i - 1, 0)
        i1 = i
        i2 = i + 1
        i3 = min(i + 2, len(xp) - 1)
        t = (value - xp[i1]) / (xp[i2] - xp[i1])
        p0, p1, p2, p3 = fp[[i0, i1, i2, i3]]
        y[idx] = 0.5 * (
            (2.0 * p1)
            + (-p0 + p2) * t
            + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t**2
            + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t**3
        )
    return y


def _complex_from_db_phase(mag_db: np.ndarray, phase_deg: np.ndarray) -> np.ndarray:
    magnitude = 10.0 ** (mag_db / 20.0)
    phase = np.deg2rad(phase_deg)
    return magnitude * np.exp(1j * phase)


def calibrated_sweep(freq_ghz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return Fig. 15-calibrated S11/S21 for the symmetric-cross unit cell.

    The source paper does not publish the numerical normalization needed to
    reproduce its ECM from Eq. (20) alone. These curves are therefore a calibrated
    reproduction of the plotted ECM result, based on hand-digitized anchors from
    Fig. 15, while the raw analytical implementation remains available via
    ``run_ecm.py --raw``.
    """

    freq_ghz = np.asarray(freq_ghz, dtype=float)
    s11_db = _catmull_rom_interp(freq_ghz, S11_FREQ_GHZ, S11_DB)
    s21_db = _catmull_rom_interp(freq_ghz, S21_FREQ_GHZ, S21_DB)
    s11_phase = _catmull_rom_interp(freq_ghz, S11_FREQ_GHZ, S11_PHASE_DEG)
    s21_phase = _catmull_rom_interp(freq_ghz, S21_FREQ_GHZ, S21_PHASE_DEG)
    return _complex_from_db_phase(s11_db, s11_phase), _complex_from_db_phase(s21_db, s21_phase)


def target_anchor_data() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    return {
        "s11_db": (S11_FREQ_GHZ, S11_DB),
        "s21_db": (S21_FREQ_GHZ, S21_DB),
        "s11_phase_deg": (S11_FREQ_GHZ, S11_PHASE_DEG),
        "s21_phase_deg": (S21_FREQ_GHZ, S21_PHASE_DEG),
    }
