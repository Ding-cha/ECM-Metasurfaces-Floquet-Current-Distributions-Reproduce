from __future__ import annotations

from dataclasses import dataclass
from dataclasses import asdict
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class SymmetricCrossParams:
    """Geometry from Fig. 13 and Table II of the paper.

    Length dimensions are stored in millimetres because the published current
    model uses millimetres in its cosine coefficients.
    """

    lx_mm: float = 20.0
    ly_mm: float = 20.0
    w_mm: float = 3.0
    px_mm: float = 30.0
    py_mm: float = 30.0
    h_mm: float = 1.524
    epsilon_r: float = 3.02

    def __post_init__(self) -> None:
        values = self.as_dict()
        for name, value in values.items():
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be a positive finite value.")
        if self.lx_mm > self.px_mm:
            raise ValueError("lx_mm must not exceed px_mm.")
        if self.ly_mm > self.py_mm:
            raise ValueError("ly_mm must not exceed py_mm.")
        if self.w_mm > min(self.lx_mm, self.ly_mm):
            raise ValueError("w_mm must not exceed the shorter arm length.")

    @property
    def lx_m(self) -> float:
        return self.lx_mm * 1e-3

    @property
    def ly_m(self) -> float:
        return self.ly_mm * 1e-3

    @property
    def w_m(self) -> float:
        return self.w_mm * 1e-3

    @property
    def px_m(self) -> float:
        return self.px_mm * 1e-3

    @property
    def py_m(self) -> float:
        return self.py_mm * 1e-3

    @property
    def h_m(self) -> float:
        return self.h_mm * 1e-3

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


def horizontal_current_x(
    x_mm: np.ndarray | float,
    params: SymmetricCrossParams | None = None,
) -> np.ndarray | float:
    """Paper Eq. (20a), scaled to arbitrary horizontal arm length.

    The published coefficients are fitted for lx=20 mm. For generalized runs,
    x is mapped to the equivalent coordinate of a 20 mm arm so that the current
    maxima/minima stay at the same normalized locations along the arm.
    """

    scale = 1.0 if params is None else 20.0 / params.lx_mm
    x_scaled = x_mm * scale
    return 0.95 * np.cos(0.12 * x_scaled) - 0.2 * np.cos(0.64 * x_scaled)


def vertical_current_y(y_mm: np.ndarray | float, params: SymmetricCrossParams) -> np.ndarray | float:
    """Paper Eq. (20b), with y in mm and sign for current direction."""

    half = params.ly_mm / 2.0
    magnitude = 0.5 * horizontal_current_x(0.0, params) * np.clip(1.0 - np.abs(y_mm) / half, 0.0, 1.0)
    return np.sign(y_mm) * magnitude


def current_grid(params: SymmetricCrossParams, samples: int = 501) -> tuple[np.ndarray, ...]:
    """Return x, y, Jx, Jy grids over the unit cell in millimetres."""

    x = np.linspace(-params.px_mm / 2.0, params.px_mm / 2.0, samples)
    y = np.linspace(-params.py_mm / 2.0, params.py_mm / 2.0, samples)
    xx, yy = np.meshgrid(x, y)

    horizontal = (np.abs(xx) <= params.lx_mm / 2.0) & (np.abs(yy) <= params.w_mm / 2.0)
    vertical = (np.abs(xx) <= params.w_mm / 2.0) & (np.abs(yy) <= params.ly_mm / 2.0)

    jx = np.zeros_like(xx)
    jy = np.zeros_like(xx)
    jx[horizontal] = horizontal_current_x(xx[horizontal], params)
    jy[vertical] = vertical_current_y(yy[vertical], params)
    return xx, yy, jx, jy


def ensure_output_dir(path: str | Path = "outputs") -> Path:
    output = Path(path)
    output.mkdir(parents=True, exist_ok=True)
    return output
