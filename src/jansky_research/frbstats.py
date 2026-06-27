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
    "grouped_wait_times",
    "select_xmin",
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
    n_repeater_sources: int


def wait_times(mjds: np.ndarray) -> np.ndarray:
    """Return sorted inter-burst waiting times (days) from burst arrival MJDs."""
    mjds = np.sort(np.asarray(mjds, dtype=float))
    dt = np.diff(mjds)
    return dt[dt > 0]


def grouped_wait_times(mjds: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Inter-burst waiting times computed *within* each source, then pooled.

    Pooling arrival times across different repeaters is meaningless — a "wait" is only defined
    between consecutive bursts of the *same* source. This computes ``wait_times`` per group label
    and concatenates them, so a catalogue with many repeaters yields one combined wait sample.
    """
    mjds = np.asarray(mjds, dtype=float)
    groups = np.asarray(groups)
    out = [wait_times(mjds[groups == g]) for g in np.unique(groups)]
    out = [w for w in out if w.size]
    return np.concatenate(out) if out else np.array([])


def fit_weibull_waits(
    mjds: np.ndarray, *, groups: np.ndarray | None = None, n_boot: int = 1000, seed: int | None = 0
) -> WeibullFit:
    """Fit a Weibull to the inter-burst waiting times, with a bootstrap CI on the shape.

    With ``groups`` (per-burst source labels) the waits are computed *within* each source and
    pooled (the correct treatment for a multi-repeater catalogue); without it the bursts are
    treated as a single source. Uses ``scipy.stats.weibull_min`` with location fixed at 0. The
    shape ``k`` is the clustering diagnostic: ``k < 1`` clustered, ``k = 1`` Poisson.
    """
    waits = grouped_wait_times(mjds, groups) if groups is not None else wait_times(mjds)
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


def _hill(tail: np.ndarray, f_min: float) -> float:
    """Hill/Clauset differential power-law index for ``tail`` above ``f_min``."""
    return 1.0 + tail.size / np.sum(np.log(tail / f_min))


def _ks_distance(tail: np.ndarray, f_min: float, gamma: float) -> float:
    """KS distance between the empirical CCDF of ``tail`` and the fitted power law."""
    x = np.sort(tail)
    cdf_emp = np.arange(1, x.size + 1) / x.size
    cdf_fit = 1.0 - (x / f_min) ** (-(gamma - 1.0))
    return float(np.max(np.abs(cdf_emp - cdf_fit)))


def select_xmin(fluences: np.ndarray, *, min_tail: int = 30) -> float:
    """Choose the power-law lower bound by the Clauset-Shalizi-Newman KS criterion.

    Scans candidate lower bounds and returns the one minimising the KS distance between the
    empirical tail and its fitted power law — i.e. where the data actually become power-law,
    above the survey's incompleteness. ``min_tail`` guards against tiny, unstable tails.
    """
    f = np.asarray(fluences, dtype=float)
    f = f[np.isfinite(f) & (f > 0)]
    candidates = np.unique(f)
    best_xmin, best_d = float(f.min()), np.inf
    for xm in candidates:
        tail = f[f >= xm]
        if tail.size < min_tail:
            break  # candidates are sorted ascending; tails only shrink from here
        d = _ks_distance(tail, xm, _hill(tail, xm))
        if d < best_d:
            best_d, best_xmin = d, float(xm)
    return best_xmin


def fit_power_law(
    fluences: np.ndarray, *, f_min: float | None = None, auto_xmin: bool = False
) -> PowerLawFit:
    """Maximum-likelihood differential power-law index of the fluence distribution.

    For dN/dF ∝ F**(-gamma) above a lower bound ``f_min``, the Clauset et al. (2009) / Hill
    estimator is ``gamma = 1 + n / sum(ln(F_i / f_min))`` with standard error
    ``(gamma - 1) / sqrt(n)``. With ``auto_xmin`` (and ``f_min`` unset) the lower bound is chosen
    by :func:`select_xmin` so the fit is over the genuine power-law tail, not the incomplete faint
    end — essential for a real survey catalogue. Otherwise ``f_min`` defaults to the smallest
    positive fluence.
    """
    f = np.asarray(fluences, dtype=float)
    f = f[np.isfinite(f) & (f > 0)]
    if f_min is None:
        f_min = select_xmin(f) if auto_xmin else float(f.min())
    tail = f[f >= f_min]
    n = tail.size
    if n < 2:
        raise ValueError("need at least 2 fluences above f_min for a power-law fit")
    gamma = _hill(tail, f_min)
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

    ``catalog`` keys: ``mjd``, ``fluence``, ``dm``, ``width`` (arrays), ``repeater`` (bool mask),
    and optionally ``repeater_name`` (per-burst source labels). The Weibull wait-time fit uses the
    repeater bursts, with waits computed *within* each source when labels are available; the energy
    and population comparisons use the whole catalogue.
    """
    repeater_mask = np.asarray(catalog["repeater"], dtype=bool)
    rep_mjd = np.asarray(catalog["mjd"])[repeater_mask]
    if "repeater_name" in catalog:
        labels = np.asarray(catalog["repeater_name"])[repeater_mask]
        weibull = fit_weibull_waits(rep_mjd, groups=labels)
        n_sources = int(np.unique(labels).size)
    else:
        weibull = fit_weibull_waits(rep_mjd)
        n_sources = 1
    energy = fit_power_law(np.asarray(catalog["fluence"]), auto_xmin=True)
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
        n_repeater_sources=n_sources,
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

    # The repeater bursts are one synthetic source; non-repeaters get the catalogue sentinel.
    names = np.array(["SYN-R1"] * n_repeater + ["-9999"] * n_oneoff)
    return {
        "mjd": np.concatenate([rep_mjd, one_mjd]),
        "fluence": np.concatenate([rep_fluence, one_fluence]),
        "dm": np.concatenate([rep_dm, one_dm]),
        "width": np.concatenate([rep_width, one_width]),
        "repeater": np.concatenate([np.ones(n_repeater, bool), np.zeros(n_oneoff, bool)]),
        "repeater_name": names,
    }
