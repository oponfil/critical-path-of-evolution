#!/usr/bin/env python3
"""Numerical verification of calculations in critical_path_ru.tex (GHZ ring model)."""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

import numpy as np
from scipy.integrate import dblquad, quad
from scipy.optimize import brentq
from scipy.stats import norm

# --- Paper constants (Section 4) -------------------------------------------

MIN_PER_YEAR = 365.25 * 24 * 60
T = 1.38e10  # years
TLOUD = 990  # years
P0 = 2e13
P_II = 1e26
P_I = 1e16

LY_PER_KPC = 1000 / 0.306601  # 1 kpc in light-years
R_MIN_KPC = 4.0
R_MAX_KPC = 14.0
R_SUN_KPC = 8.0
H_LY = 1e3

R_MIN_LY = R_MIN_KPC * LY_PER_KPC
R_MAX_LY = R_MAX_KPC * LY_PER_KPC
R_SUN_LY = R_SUN_KPC * LY_PER_KPC

# V_GHZ = pi (r_max^2 - r_min^2) h with radii in ly (silence integral, note to Table 2)
V_GHZ = math.pi * (R_MAX_LY**2 - R_MIN_LY**2) * H_LY

GAMMA43 = math.gamma(4 / 3)

# Expected table values from critical_path_ru.tex
TABLE1 = [
    (0.015, 417, 1964),
    (0.020, 314, 1477),
    (0.025, 252, 1184),
    (0.030, 210, 990),
    (0.035, 181, 850),
    (0.040, 158, 746),
]

TABLE2 = [
    (1e2, 6960, 2170),
    (1e3, 2690, 1010),
    (1e4, 1290, 470),
    (2.9e4, 990, 330),
    (1e5, 775, 220),
    (1e6, 550, 101),
    (1e7, 430, 50),
    (1e8, 350, 22),
]

TABLE3 = [
    (1e2, 6960, 2170, 0.42, 1840, 3.9e12),
    (1e3, 2690, 1010, 1.1, 276, 2.6e13),
    (1e4, 1290, 470, 2.3, 63, 1.1e14),
    (2.9e4, 990, 330, 3.0, 37, 1.9e14),
    (1e5, 775, 220, 3.8, 23, 3.2e14),
    (1e6, 550, 101, 5.3, 12, 6.3e14),
    (1e7, 430, 50, 6.8, 7.0, 1.0e15),
    (1e8, 350, 22, 8.4, 4.7, 1.6e15),
]

N_MC = 400_000
RD_KPC = 3.0


@dataclass
class Result:
    name: str
    ok: bool
    detail: str = ""


results: list[Result] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    results.append(Result(name, ok, detail))


def near(a: float, b: float, *, atol: float = 0.0, rtol: float = 0.0) -> bool:
    return abs(a - b) <= max(atol, rtol * abs(b))


def ell_mean_distance(n: float) -> float:
    """Table 2 caption: ell = Gamma(4/3) (3 (r_max^2-r_min^2) h / (4 n))^(1/3), r in ly, no pi."""
    area = R_MAX_LY**2 - R_MIN_LY**2
    return GAMMA43 * (3 * area * H_LY / (4 * n)) ** (1 / 3)


def mu_min_from_sigma(sigma_years: float) -> float:
    return (sigma_years**2 / T) * MIN_PER_YEAR


def N_from_sigma(sigma_years: float) -> float:
    return (T / sigma_years) ** 2


def growth_pct_from_sigma(sigma_years: float) -> float:
    """Growth rate (percent/year) such that t_loud = sigma."""
    ratio = P_II / P0
    return (math.exp(math.log(ratio) / sigma_years) - 1.0) * 100.0


def t_kardashev(p_target: float, rate: float) -> float:
    return math.log(p_target / P0) / math.log(1.0 + rate)


def K_moment(a: float, y: float) -> float:
    apy = a + y
    return 0.5 * (
        (1 + a * a - y * y) * norm.cdf(-apy)
        + (y - a) * norm.pdf(apy)
    )


def W_closed(a: float, beta: float) -> float:
    value, _ = quad(lambda y: K_moment(a, y), 0.0, beta, limit=200)
    return value


def E_V_closed(n: float, sigma_years: float) -> float:
    """Closed form from the note after Eq. (5): sigma_r = c sigma, c = 1 ly/year."""
    sigma_r = sigma_years
    a = TLOUD / sigma_years
    beta = H_LY / (2.0 * sigma_r)
    return (n / V_GHZ) * 4.0 * math.pi * sigma_r**3 * W_closed(a, beta)


def solve_sigma(n: float) -> float:
    return brentq(lambda s: E_V_closed(n, s) - 1.0, 50.0, 50_000.0)


def E_V_ring_numeric(n: float, sigma_years: float) -> float:
    """Full azimuthal integral; observer at (R_sun, 0, 0) in ly."""

    def integrand(r: float, z: float) -> float:
        def phi_integrand(phi: float) -> float:
            s = math.sqrt(
                r * r
                + R_SUN_LY * R_SUN_LY
                - 2.0 * r * R_SUN_LY * math.cos(phi)
                + z * z
            )
            return norm.cdf(-(TLOUD + s) / sigma_years)

        phi_val, _ = quad(phi_integrand, 0.0, 2.0 * math.pi, limit=200)
        return r * phi_val

    val, _ = dblquad(
        integrand,
        -H_LY / 2.0,
        H_LY / 2.0,
        lambda _z: R_MIN_LY,
        lambda _z: R_MAX_LY,
        epsabs=1e-3,
        epsrel=1e-3,
    )
    return (n / V_GHZ) * val


def _sample_r_uniform_volume(rng: np.random.Generator, count: int) -> np.ndarray:
    u = rng.random(count)
    return np.sqrt(R_MIN_LY**2 + u * (R_MAX_LY**2 - R_MIN_LY**2))


def _sample_r_exponential(rng: np.random.Generator, count: int) -> np.ndarray:
    rd = RD_KPC * LY_PER_KPC
    grid = np.linspace(R_MIN_LY, R_MAX_LY, 4000)
    weights = grid * np.exp(-grid / rd)
    weights /= np.trapezoid(weights, grid)
    cdf = np.cumsum(weights)
    cdf /= cdf[-1]
    u = rng.random(count)
    return np.interp(u, cdf, grid)


def _distances_from_r(r: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    phi = rng.uniform(0.0, 2.0 * math.pi, r.size)
    z = rng.uniform(-H_LY / 2.0, H_LY / 2.0, r.size)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return np.sqrt((x - R_SUN_LY) ** 2 + y**2 + z**2)


def mc_sigma_shift(reference_n: float = 2.9e4) -> tuple[float, float, float]:
    """Return (sigma uniform MC, sigma exp MC, sigma closed) at E[V]=1."""
    rng = np.random.default_rng(0)

    def mc_ev(n: float, sigma_years: float, distances: np.ndarray) -> float:
        return float(n * np.mean(norm.cdf(-(TLOUD + distances) / sigma_years)))

    sigma_closed = solve_sigma(reference_n)

    def root_sigma(distances: np.ndarray) -> float:
        objective = lambda s: mc_ev(reference_n, s, distances) - 1.0
        return brentq(objective, 0.85 * sigma_closed, 1.15 * sigma_closed)

    s_ring = root_sigma(_distances_from_r(_sample_r_uniform_volume(rng, N_MC), rng))
    s_exp = root_sigma(_distances_from_r(_sample_r_exponential(rng, N_MC), rng))
    return s_ring, s_exp, sigma_closed


def run_checks() -> int:
    # --- Table 1 ---
    for rate, t1_paper, t2_paper in TABLE1:
        t1 = t_kardashev(P_I, rate)
        t2 = t_kardashev(P_II, rate)
        check(
            f"Table 1: {100 * rate:.1f}% -> Type I",
            near(t1, t1_paper, atol=2),
            f"calc={t1:.0f}, paper={t1_paper}",
        )
        check(
            f"Table 1: {100 * rate:.1f}% -> Type II (t_loud)",
            near(t2, t2_paper, atol=2),
            f"calc={t2:.0f}, paper={t2_paper}",
        )

    check("t_loud at 3%/year", near(t_kardashev(P_II, 0.03), TLOUD, atol=2))

    # --- Table 2 (silence integral) ---
    for n, sigma_paper, ell_paper in TABLE2:
        sigma_calc = solve_sigma(n)
        check(
            f"Table 2: n={n:.0e} sigma",
            near(sigma_calc, sigma_paper, atol=5),
            f"calc={sigma_calc:.0f}, paper={sigma_paper}",
        )
        ell_calc = ell_mean_distance(n)
        check(
            f"Table 2: n={n:.0e} ell",
            near(ell_calc, ell_paper, atol=5 if n >= 1e5 else 12),
            f"calc={ell_calc:.0f}, paper={ell_paper}",
        )
        if near(sigma_paper, sigma_calc, atol=2):
            check(
                f"Table 2: n={n:.0e} E[V] at paper sigma",
                near(E_V_closed(n, sigma_paper), 1.0, rtol=0.08, atol=0.08),
                f"E[V]={E_V_closed(n, sigma_paper):.4f}",
            )

    # Reference row: n such that sigma = t_loud
    n_ref = brentq(lambda n: solve_sigma(n) - TLOUD, 1e3, 1e5)
    check(
        "Reference n (sigma = t_loud)",
        near(n_ref, 2.9e4, rtol=0.02),
        f"calc={n_ref:.3g}",
    )

    # --- Table 3 ---
    for n, sigma_p, ell_p, growth_p, mu_p, N_p in TABLE3:
        check(
            f"Table 3: n={n:.0e} mu",
            near(mu_min_from_sigma(sigma_p), mu_p, rtol=0.02, atol=1.5),
            f"calc={mu_min_from_sigma(sigma_p):.1f}, paper={mu_p}",
        )
        check(
            f"Table 3: n={n:.0e} N",
            near(N_from_sigma(sigma_p), N_p, rtol=0.05),
            f"calc={N_from_sigma(sigma_p):.2e}, paper={N_p:.2e}",
        )
        check(
            f"Table 3: n={n:.0e} growth %",
            near(growth_pct_from_sigma(sigma_p), growth_p, atol=0.4),
            f"calc={growth_pct_from_sigma(sigma_p):.2f}, paper={growth_p}",
        )

    # --- Core inversions ---
    sigma_40 = math.sqrt((40 / MIN_PER_YEAR) * T)
    check(
        "mu=40 min -> sigma",
        near(sigma_40, 1024, atol=5),
        f"calc={sigma_40:.1f}",
    )
    check(
        "sigma=990 -> mu (min)",
        near(mu_min_from_sigma(990), 37.4, atol=0.5),
        f"calc={mu_min_from_sigma(990):.2f}",
    )
    check(
        "sigma=990 -> N",
        near(N_from_sigma(990), 1.94e14, rtol=0.02),
        f"calc={N_from_sigma(990):.3e}",
    )

    mu_bio_years = 40 / MIN_PER_YEAR
    cv = 990 / math.sqrt(mu_bio_years * T)
    check("c_v from channels", near(cv, 0.97, atol=0.02), f"calc={cv:.3f}")

    # --- Table 4 endpoints ---
    sigma_lo = solve_sigma(1e5)
    sigma_hi = solve_sigma(1e4)
    check("Table 4: sigma in [775, 1290]", 775 <= sigma_lo <= 785 and 1285 <= sigma_hi <= 1295)
    check(
        "Table 4: mu range",
        20 <= mu_min_from_sigma(sigma_hi) <= 80
        and 20 <= mu_min_from_sigma(sigma_lo) <= 80,
    )
    check(
        "Table 4: ell range (ly)",
        200 <= ell_mean_distance(1e5) <= 480
        and 200 <= ell_mean_distance(1e4) <= 480,
    )

    # --- Weak sigma(n) scaling (Section 4.5) ---
    ratio_sigma = solve_sigma(1e4) / solve_sigma(1e5)
    check(
        "sigma(1e4)/sigma(1e5) ~ 1.7",
        near(ratio_sigma, 1290 / 775, atol=0.08),
        f"calc={ratio_sigma:.3f}",
    )

    # --- Closed form vs numeric ring integral (reference point) ---
    ev_closed = E_V_closed(2.9e4, 990)
    ev_numeric = E_V_ring_numeric(2.9e4, 990)
    check(
        "Closed form vs numeric ring integral",
        near(ev_closed, ev_numeric, rtol=0.05),
        f"closed={ev_closed:.4f}, numeric={ev_numeric:.4f}",
    )

    # --- Monte Carlo (geometry caveat, Section 4.2) ---
    s_ring, s_exp, s_closed = mc_sigma_shift()
    check(
        "MC uniform ring matches closed sigma",
        near(s_ring, s_closed, rtol=0.01),
        f"mc={s_ring:.1f}, closed={s_closed:.1f}",
    )
    check(
        "MC exp profile lowers sigma vs uniform",
        s_exp < s_ring,
        f"uniform={s_ring:.1f}, exp={s_exp:.1f}, shift={(s_exp - s_ring) / s_ring * 100:.1f}%",
    )

    # --- Kardashev K ---
    k = (math.log10(P0) - 6.0) / 10.0
    check("Kardashev K at P0", near(k, 0.73, atol=0.01), f"calc={k:.3f}")

    passed = sum(r.ok for r in results)
    total = len(results)
    print(f"\n{passed}/{total} checks passed\n")
    for r in results:
        status = "OK" if r.ok else "FAIL"
        line = f"[{status}] {r.name}"
        if r.detail:
            line += f"  ({r.detail})"
        print(line)
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(run_checks())
