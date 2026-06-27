"""FRB burst-statistics — the gap-filling tool (GATE 1: FRB burst-statistics).

The open-source ecosystem has no lightweight, pip-installable library for the standard
fast-radio-burst burst-statistics analyses; they live as one-off scripts inside individual
papers, and FRBSTATS is web-only (see ``survey/github-landscape.md``). This module provides the
three core analyses as pure-NumPy/SciPy functions, composable from a notebook, a CLI, or an
Airflow task:

1. **Wait-time clustering** — fit a Weibull distribution to the inter-burst waiting times of a
   repeater. Shape ``k < 1`` means clustered (the empirical finding for hyperactive repeaters),
   ``k = 1`` is a memoryless Poisson process. Bootstrap confidence interval on ``k``.
2. **Energy/fluence distribution** — the differential power-law index ``gamma`` of the burst
   fluence distribution (dN/dF ∝ F^-gamma) via the Clauset/Hill maximum-likelihood estimator.
3. **Repeater vs non-repeater** — two-sample KS tests on DM / fluence / width, the live question
   of whether the two classes are drawn from the same population.

All functions take plain NumPy arrays so they are trivially testable against the synthetic
catalogue (:func:`synthetic_catalog`) with known ground truth — no network, no GPU.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

__all__ = [
    "BurstStats",
    "PowerLawFit",
    "WeibullFit",
    "compare_populations",
    "fit_power_law",
    "fit_weibull_waits",
    "summarise",
    "synthetic_catalog",
    "wait_times",
]


@dataclass(frozen=True)
class WeibullFit:
    """Weibull fit to inter-burst waiting times. ``k`` is the shape parameter."""

    k: float
    scale: float
    k_ci_low: float
    k_ci_high: float
    n_waits: int

    @property
    def clustered(self) -> bool:
        """True if the 95% CI on the shape excludes 1 from above (k < 1: clustering)."""
        return self.k_ci_high < 1.0


@dataclass(frozen=True)
class PowerLawFit:
    """Power-law fit to the fluence distribution, dN/dF ∝ F**(-gamma)."""

    gamma: float
    gamma_err: float
    f_min: float
    n_tail: int


@dataclass(frozen=True)
class BurstStats:
    """Top-level result bundle the pipeline serialises to ``results/metrics.json``."""

    weibull: WeibullFit
    energy: PowerLawFit
    ks: dict[str, dict[str, float]]
    n_bursts: int
    n_repeater_bursts: int


def wait_times(mjds: np.ndarray) -> np.ndarray:
    """Return sorted inter-burst waiting times (days) from burst arrival MJDs."""
    mjds = np.sort(np.asarray(mjds, dtype=float))
    dt = np.diff(mjds)
    return dt[dt > 0]


def fit_weibull_waits(mjds: np.ndarray, *, n_boot: int = 1000, seed: int | None = 0) -> WeibullFit:
    """Fit a Weibull to the inter-burst waiting times, with a bootstrap CI on the shape.

    Uses ``scipy.stats.weibull_min`` with the location fixed at 0 (waiting times are positive).
    The shape ``k`` is the clustering diagnostic: ``k < 1`` clustered, ``k = 1`` Poisson.
    """
    waits = wait_times(mjds)
    if waits.size < 3:
        raise ValueError("need at least 3 positive waiting times to fit a Weibull")
    k, _loc, scale = stats.weibull_min.fit(waits, floc=0.0)

    rng = np.random.default_rng(seed)
    boot_k = np.empty(n_boot)
    for i in range(n_boot):
        sample = rng.choice(waits, size=waits.size, replace=True)
        boot_k[i], _, _ = stats.weibull_min.fit(sample, floc=0.0)
    lo, hi = np.percentile(boot_k, [2.5, 97.5])
    return WeibullFit(
        k=float(k), scale=float(scale), k_ci_low=float(lo), k_ci_high=float(hi), n_waits=waits.size
    )


def fit_power_law(fluences: np.ndarray, *, f_min: float | None = None) -> PowerLawFit:
    """Maximum-likelihood differential power-law index of the fluence distribution.

    For dN/dF ∝ F**(-gamma) above a lower bound ``f_min``, the Clauset et al. (2009) / Hill
    estimator is ``gamma = 1 + n / sum(ln(F_i / f_min))`` with standard error
    ``(gamma - 1) / sqrt(n)``. If ``f_min`` is None it defaults to the smallest positive fluence
    (the whole sample is treated as the power-law tail).
    """
    f = np.asarray(fluences, dtype=float)
    f = f[np.isfinite(f) & (f > 0)]
    if f_min is None:
        f_min = float(f.min())
    tail = f[f >= f_min]
    n = tail.size
    if n < 2:
        raise ValueError("need at least 2 fluences above f_min for a power-law fit")
    gamma = 1.0 + n / np.sum(np.log(tail / f_min))
    return PowerLawFit(
        gamma=float(gamma), gamma_err=float((gamma - 1.0) / np.sqrt(n)), f_min=f_min, n_tail=n
    )


def compare_populations(
    repeater: dict[str, np.ndarray], oneoff: dict[str, np.ndarray], keys=("dm", "fluence", "width")
) -> dict[str, dict[str, float]]:
    """Two-sample KS tests comparing repeater vs non-repeater burst properties.

    Returns ``{key: {"ks": D, "p": p, "n_rep": .., "n_one": ..}}`` for each available key.
    A small p-value is evidence the two classes differ in that property.
    """
    out: dict[str, dict[str, float]] = {}
    for key in keys:
        if key not in repeater or key not in oneoff:
            continue
        a = np.asarray(repeater[key], dtype=float)
        b = np.asarray(oneoff[key], dtype=float)
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        if a.size < 2 or b.size < 2:
            continue
        d, p = stats.ks_2samp(a, b)
        out[key] = {"ks": float(d), "p": float(p), "n_rep": int(a.size), "n_one": int(b.size)}
    return out


def summarise(catalog: dict[str, np.ndarray]) -> BurstStats:
    """Run all three analyses over a catalogue dict and bundle the result.

    ``catalog`` keys: ``mjd``, ``fluence``, ``dm``, ``width`` (arrays), and ``repeater`` (bool
    array). The Weibull wait-time fit uses the repeater bursts (the bursts with a recurrence
    structure); the energy and population comparisons use the whole catalogue.
    """
    repeater_mask = np.asarray(catalog["repeater"], dtype=bool)
    weibull = fit_weibull_waits(np.asarray(catalog["mjd"])[repeater_mask])
    energy = fit_power_law(np.asarray(catalog["fluence"]))
    rep = {
        k: np.asarray(catalog[k])[repeater_mask] for k in ("dm", "fluence", "width") if k in catalog
    }
    one = {
        k: np.asarray(catalog[k])[~repeater_mask]
        for k in ("dm", "fluence", "width")
        if k in catalog
    }
    ks = compare_populations(rep, one)
    return BurstStats(
        weibull=weibull,
        energy=energy,
        ks=ks,
        n_bursts=int(repeater_mask.size),
        n_repeater_bursts=int(repeater_mask.sum()),
    )


def synthetic_catalog(
    n_repeater: int = 400,
    n_oneoff: int = 200,
    k_true: float = 0.7,
    scale_days: float = 0.05,
    gamma_true: float = 2.0,
    fluence_min: float = 0.4,
    seed: int | None = 0,
) -> dict[str, np.ndarray]:
    """Generate a synthetic FRB catalogue with known ground truth (offline fixture).

    Repeater bursts have Weibull-distributed (shape ``k_true``) waiting times and power-law
    fluences (index ``gamma_true``); non-repeaters are a separate Poisson/lognormal draw with a
    shifted DM so the population comparison has signal. Used by the tests and as the offline
    fallback when the real CHIME catalogue cannot be downloaded.
    """
    rng = np.random.default_rng(seed)
    # Repeater arrival times: cumulative Weibull waits.
    waits = stats.weibull_min.rvs(k_true, scale=scale_days, size=n_repeater - 1, random_state=rng)
    rep_mjd = np.concatenate([[58000.0], 58000.0 + np.cumsum(waits)])
    # Power-law fluences via inverse-CDF: F = f_min * (1-u)^(-1/(gamma-1)).
    u = rng.random(n_repeater)
    rep_fluence = fluence_min * (1.0 - u) ** (-1.0 / (gamma_true - 1.0))
    rep_dm = rng.normal(500.0, 120.0, n_repeater)
    rep_width = 10 ** rng.normal(0.0, 0.4, n_repeater)  # ms, lognormal

    one_mjd = rng.uniform(58000.0, 58400.0, n_oneoff)
    uo = rng.random(n_oneoff)
    one_fluence = fluence_min * (1.0 - uo) ** (-1.0 / (gamma_true - 1.0))
    one_dm = rng.normal(650.0, 150.0, n_oneoff)  # shifted DM -> KS has signal
    one_width = 10 ** rng.normal(0.1, 0.4, n_oneoff)

    return {
        "mjd": np.concatenate([rep_mjd, one_mjd]),
        "fluence": np.concatenate([rep_fluence, one_fluence]),
        "dm": np.concatenate([rep_dm, one_dm]),
        "width": np.concatenate([rep_width, one_width]),
        "repeater": np.concatenate([np.ones(n_repeater, bool), np.zeros(n_oneoff, bool)]),
    }
