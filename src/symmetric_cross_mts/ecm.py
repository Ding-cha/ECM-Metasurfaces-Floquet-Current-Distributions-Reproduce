from __future__ import annotations

from dataclasses import replace
from dataclasses import dataclass

import numpy as np

from .fourier import current_fourier
from .model import SymmetricCrossParams

EPS0 = 8.8541878128e-12
MU0 = 4.0e-7 * np.pi
C0 = 1.0 / np.sqrt(EPS0 * MU0)


@dataclass(frozen=True)
class EcmConfig:
    """ECM 求解配置。

    `restored` 模式会使用一组通过探索得到的配置；这里仍保留多个开关，
    用来复现论文公式、对比原始实现、以及排查高阶 Floquet 模态的影响。
    """

    mode_order: int = 7
    theta_deg: float = 0.0
    phi_deg: float = 90.0
    fourier_samples: int = 1201
    current_scale: float = 1.0
    t00_override: float | None = None
    ttr_variant: str = "paper"
    kz_branch: str = "principal"
    high_order_substrate: bool = True
    width_profile: str = "uniform"
    epsilon_eff_scale: float = 1.0
    arm_length_correction_mm: float = 0.0


def _sqrt_branch(value: complex, branch: str = "principal") -> complex:
    """选择纵向波数 kz 的复数平方根分支。"""

    root = np.sqrt(value + 0j)
    # 对 evanescent 模态，decay 分支把 +j alpha 改成 -j alpha。
    # 这个符号会直接改变高阶模态的等效电抗，是恢复作者曲线时最敏感的因素之一。
    if branch == "decay" and abs(root.real) < 1e-14 and root.imag > 0:
        root = -root
    return root


def _modal_quantities(freq_hz: float, params: SymmetricCrossParams, m: int, n: int, cfg: EcmConfig):
    omega = 2.0 * np.pi * freq_hz
    k0 = omega / C0
    theta = np.deg2rad(cfg.theta_deg)
    phi = np.deg2rad(cfg.phi_deg)

    kx = k0 * np.sin(theta) * np.cos(phi) + 2.0 * np.pi * m / params.px_m
    ky = k0 * np.sin(theta) * np.sin(phi) + 2.0 * np.pi * n / params.py_m

    kz_air = _sqrt_branch(k0**2 - kx**2 - ky**2, cfg.kz_branch)
    eps_eff = 1.0 + cfg.epsilon_eff_scale * (params.epsilon_r - 1.0)
    kz_diel = _sqrt_branch(k0**2 * eps_eff - kx**2 - ky**2, cfg.kz_branch)

    y_te_plus = kz_air / (omega * MU0)
    y_tm_plus = omega * EPS0 / kz_air
    y_te_minus = kz_diel / (omega * MU0)
    y_tm_minus = omega * EPS0 * eps_eff / kz_diel
    return kx, ky, kz_air, kz_diel, y_te_plus, y_tm_plus, y_te_minus, y_tm_minus


def _loaded_admittance(y_plus: complex, y_minus: complex, kz_diel: complex, h_m: float) -> complex:
    """介质层传输线输入导纳，对应论文 Eq. (8) 的 B_mn。"""

    zc = 1.0 / y_minus
    z_load = 1.0 / y_plus
    tan_term = np.tan(kz_diel * h_m)
    z_in = zc * (z_load + 1j * zc * tan_term) / (zc + 1j * z_load * tan_term)
    return 1.0 / z_in


def _width_integral(k: complex, half_width_m: float, profile: str) -> complex:
    """金属条宽度方向的积分。

    uniform 是旧 restored 的假设；cosine 用一个中心强、边缘弱的横向电流
    profile 来近似宽度方向电流变化。Fig. 22 宽度扫描显示，这个修正能把
    高频反射谷拉回到更接近作者结果的位置。
    """

    if profile == "uniform":
        if abs(k) < 1e-14:
            return 2.0 * half_width_m
        return 2.0 * np.sin(k * half_width_m) / k
    if profile != "cosine":
        raise ValueError(f"Unknown width profile: {profile}")

    u = np.linspace(-half_width_m, half_width_m, 241)
    weight = 0.5 + 0.5 * np.cos(np.pi * u / half_width_m)
    weight = weight * (2.0 * half_width_m / np.trapz(weight, u))
    return complex(np.trapz(weight * np.exp(1j * k * u), u))


def _current_fourier_for_cfg(
    kx: complex,
    ky: complex,
    params: SymmetricCrossParams,
    cfg: EcmConfig,
) -> tuple[complex, complex]:
    current_params = params
    if cfg.arm_length_correction_mm:
        lx_eff = max(params.w_mm, params.lx_mm - cfg.arm_length_correction_mm)
        ly_eff = max(params.w_mm, params.ly_mm - cfg.arm_length_correction_mm)
        current_params = replace(params, lx_mm=lx_eff, ly_mm=ly_eff)

    if cfg.width_profile == "uniform":
        return current_fourier(kx, ky, current_params, cfg.fourier_samples)

    x_mm = np.linspace(-current_params.lx_mm / 2.0, current_params.lx_mm / 2.0, cfg.fourier_samples)
    x_m = x_mm * 1e-3
    y_mm = np.linspace(-current_params.ly_mm / 2.0, current_params.ly_mm / 2.0, cfg.fourier_samples)
    y_m = y_mm * 1e-3

    from .model import horizontal_current_x, vertical_current_y

    ix = horizontal_current_x(x_mm, current_params)
    iy = vertical_current_y(y_mm, current_params)
    trapz = getattr(np, "trapezoid", np.trapz)
    jx = trapz(ix * np.exp(1j * kx * x_m), x_m) * _width_integral(ky, params.w_m / 2.0, cfg.width_profile)
    jy = _width_integral(kx, params.w_m / 2.0, cfg.width_profile) * trapz(
        iy * np.exp(1j * ky * y_m), y_m
    )
    return jx / (4.0 * np.pi**2), jy / (4.0 * np.pi**2)


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
        # 法向入射 TE00 是 kx=ky=0 的极限点。论文没有展开该极限的
        # 归一化细节，这里保留一个可运行的 x 极化近似。
        j_te = jx
        j_tm = 0.0j
        return cfg.current_scale**2 * factor * j_te**2, cfg.current_scale**2 * factor * j_tm**2

    # FT_TE / FT_TM 是矢量电流傅里叶变换在 TE/TM 极化方向上的投影。
    ft_te = jx * ky - jy * kx
    ft_tm = jx * kx + jy * ky
    if cfg.ttr_variant == "polarized":
        # restored 候选使用这个形式：先投影到单位极化方向，再构造 TTR。
        # 它比论文排版式 Eq. (10) 更接近 Fig. 15 的数值曲线。
        t_te = factor * ft_te**2 / (kx**2 + ky**2)
        t_tm = factor * ft_tm**2 / (kx**2 + ky**2)
        return cfg.current_scale**2 * t_te, cfg.current_scale**2 * t_tm
    if cfg.ttr_variant == "no_extra_k":
        # 对照变体：不额外除以横向波数，主要用于探索量纲/归一化误差。
        return cfg.current_scale**2 * factor * ft_te**2, cfg.current_scale**2 * factor * ft_tm**2
    if cfg.ttr_variant == "paper":
        # 论文 Eq. (10) 的直接实现。当前项目保留它作为 raw 模式对照。
        t_te = factor * ft_te**2 * ky**2 / (kx**2 + ky**2)
        t_tm = factor * ft_tm**2 * kx**2 / (kx**2 + ky**2)
        return cfg.current_scale**2 * t_te, cfg.current_scale**2 * t_tm
    raise ValueError(f"Unknown TTR variant: {cfg.ttr_variant}")


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
            if cfg.high_order_substrate or (m == 0 and n == 0):
                b_te = _loaded_admittance(y_te_p, y_te_m, kz_diel, params.h_m)
                b_tm = _loaded_admittance(y_tm_p, y_tm_m, kz_diel, params.h_m)
            else:
                # restored 候选只保留基模的介质层负载；高阶 evanescent 模态用
                # 空气侧等效导纳近似。这样主谐振会回到作者图中的 5.8 GHz 附近。
                b_te = y_te_p
                b_tm = y_tm_p
            y_te = y_te_p + b_te
            y_tm = y_tm_p + b_tm
            jx, jy = _current_fourier_for_cfg(kx, ky, params, cfg)
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

    if cfg.t00_override is not None:
        # 探索脚本用它测试 TE00 端口归一化对曲线的影响。
        t00_te = cfg.t00_override

    # 论文 Fig. 1 / Eq. (12)-(13)：先求 TE00 输入阻抗，再换算 S11。
    z_in = 1.0 / (1.0 / (z_ab / t00_te) + b00_te)
    s11 = (z_in - z00_te_plus) / (z_in + z00_te_plus)

    # 论文 Eq. (14) 的传输线电压分配形式。这里仍是近似实现，
    # 主要用于与恢复候选和校准曲线做趋势对比。
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
