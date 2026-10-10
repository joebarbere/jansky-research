"""Coverage audit of per-burst FRB redshift posteriors p(z|DM) (plan 99).

Most FRBs have no host, so their redshifts come from p(z|DM). This module asks whether
the *intervals* such posteriors quote are calibrated: on bursts whose host redshift is
known, does the central 68% interval contain z_true 68% of the time?

Pieces (all offline, NumPy/SciPy/astropy):

* **Sample assembly** — dedupe the FRBs/FRB JSONs, split production vs certification
  bursts against the fit sample of the zdm parameter state under audit (plan 99 step 0).
* **E2, a minimal posterior** — Macquart et al. 2020 (Nature 581, 391, Methods eq. 4)
  p(DM_cosmic|z) with sigma_DM = F z^-1/2 and alpha = beta = 3 (the form zdm's
  ``pcosmic.py`` implements), convolved with a log-normal host DM redshifted by 1/(1+z),
  under a prior uniform in comoving volume, with no survey selection. It exists as a
  comparison arm for zdm (E1), not as a replacement.
* **PIT and coverage** — PIT_b = P(z < z_true | DM_b); a central q interval covers z_true
  iff PIT lies in [(1-q)/2, (1+q)/2].
* **The calibration rule and its controls** — C0 measures the rule's false-fail rate on
  E2's own generative model (calibrated by construction) before the rule is frozen; C1
  measures its power against a planted misspecification (host mean x2, or F x2).

Discretisation: DM_EG is binned (``dm_step``); the sum of a cosmic and a host bin is put
in the bin of the summed lower edges, a <= one-bin (2 pc cm^-3) offset that is negligible
against the >= tens of pc cm^-3 scatter. z is a geometric grid; PIT interpolates its CDF.
C0 measures whether either approximation leaves the null visibly miscalibrated.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM, Planck18
from scipy import optimize, signal, special, stats

__all__ = [
    "E2Params",
    "HOFFMANN_EMIN25",
    "mean_dm_cosmic",
    "solve_c0",
    "cosmic_bin_probs",
    "host_bin_probs",
    "LikelihoodTable",
    "build_table",
    "posterior",
    "pit",
    "central_interval",
    "coverage",
    "CalibrationRule",
    "null_rule",
    "simulate_sample",
    "run_controls",
    "dedupe",
    "provenance_side",
]

# --- E2 model ------------------------------------------------------------------------


@dataclass(frozen=True)
class E2Params:
    """Parameters of the minimal posterior. Defaults = zdm state ``HoffmannEmin25``
    (Hoffmann et al. 2025, PASA 42, e017; zdm ``states.py``): host log10 mean/sigma,
    F fixed at 0.32 (their Sec. 4), DM_halo = 50, H0 = 70.23; Omega_b h^2 = 0.0224
    held fixed as in zdm ``cosmology.py`` (DEF_Omega_b_h2)."""

    lmean: float = 2.18
    lsigma: float = 0.42
    F: float = 0.32
    dm_halo: float = 50.0
    H0: float = 70.23
    ombh2: float = 0.0224
    f_d: float = 0.844  # fraction of baryons in the diffuse IGM (Macquart+2020, Methods)
    f_e: float = 0.875  # electrons per baryon, fully ionised H + He with Y = 0.25

    def cosmology(self) -> FlatLambdaCDM:
        h = self.H0 / 100.0
        return FlatLambdaCDM(
            H0=self.H0,
            Om0=Planck18.Om0,
            Ob0=self.ombh2 / h**2,
            Tcmb0=Planck18.Tcmb0,
            Neff=Planck18.Neff,
            m_nu=Planck18.m_nu,
        )


HOFFMANN_EMIN25 = E2Params()

_ALPHA = 3.0
_BETA = 3.0


def mean_dm_cosmic(z: np.ndarray, p: E2Params = HOFFMANN_EMIN25) -> np.ndarray:
    """<DM_cosmic>(z) in pc cm^-3: 3 c H0 Omega_b f_d f_e / (8 pi G m_p) int (1+z)/E dz.

    Integrated cumulatively on a fine grid and interpolated, so ``z`` may be any shape.
    """
    cosmo = p.cosmology()
    pref = (
        (3 * const.c * cosmo.H0 * cosmo.Ob0 * p.f_d * p.f_e / (8 * np.pi * const.G * const.m_p))
        .to(u.pc / u.cm**3)
        .value
    )
    z = np.asarray(z, dtype=float)
    zmax = max(float(np.max(z)), 1e-3) * 1.001
    zz = np.linspace(0.0, zmax, 20001)
    integrand = (1 + zz) / cosmo.efunc(zz)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(zz))])
    return pref * np.interp(z, zz, cum)


def _delta_grid(n: int = 6000) -> np.ndarray:
    return np.geomspace(1e-3, 200.0, n)


def _pdelta_unnorm(delta: np.ndarray, sigma: float, c0: float) -> np.ndarray:
    """Macquart eq.-4 shape, scaled so its maximum is 1 (log-space: no underflow)."""
    logp = -_BETA * np.log(delta) - (delta**-_ALPHA - c0) ** 2 / (2 * _ALPHA**2 * sigma**2)
    return np.exp(logp - logp.max())


def solve_c0(sigma: float, delta: np.ndarray | None = None) -> float:
    """C0 such that <Delta> = 1 under the Macquart eq.-4 shape (bracketed root).

    <Delta> falls monotonically with C0; the lower bracket is extended until it brackets
    (C0 ~ -6 at sigma = 3.2, ~ -1900 at sigma = 40). The z grid nevertheless starts at
    z = 0.01: below it the cosmic DM (<~ 9 pc/cc) is under the DM_ISM/halo uncertainty,
    the comoving-volume prior there holds < 1e-5 of its mass, and no certification burst
    lies below z = 0.019.
    """
    d = _delta_grid() if delta is None else delta

    def mean_minus_one(c0: float) -> float:
        pdf = _pdelta_unnorm(d, sigma, c0)
        return float(np.trapezoid(pdf * d, d) / np.trapezoid(pdf, d) - 1.0)

    lo = -1.0
    while mean_minus_one(lo) < 0:
        lo *= 2
        if lo < -1e4:
            raise ValueError(f"no C0 gives <Delta> = 1 at sigma = {sigma:.3g}")
    return float(optimize.brentq(mean_minus_one, lo, 50.0, xtol=1e-10))


def cosmic_bin_probs(z: float, edges: np.ndarray, p: E2Params = HOFFMANN_EMIN25) -> np.ndarray:
    """P(DM_cosmic in each [edges[k], edges[k+1])) at redshift z (sums to <= 1)."""
    sigma = p.F * z**-0.5
    d = _delta_grid()
    pdf = _pdelta_unnorm(d, sigma, solve_c0(sigma, d))
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(d))])
    cdf /= cdf[-1]
    mu = float(mean_dm_cosmic(np.array([z]), p)[0])
    c = np.interp(edges / mu, d, cdf, left=0.0, right=1.0)
    return np.diff(c)  # not renormalised: mass beyond the last edge is the Delta^-3 tail


def host_bin_probs(z: float, edges: np.ndarray, p: E2Params = HOFFMANN_EMIN25) -> np.ndarray:
    """P(DM_host/(1+z) in each bin) for log10 DM_host ~ N(lmean, lsigma) (exact CDF)."""
    with np.errstate(divide="ignore"):
        x = np.log10(edges * (1 + z))
    c = special.ndtr((x - p.lmean) / p.lsigma)
    probs = np.diff(c)
    return probs / probs.sum()


@dataclass
class LikelihoodTable:
    """P(DM_EG bin | z) on a z grid, plus the prior weights per z node."""

    z: np.ndarray  # (nz,)
    edges: np.ndarray  # (nbins+1,)
    prob: np.ndarray  # (nz, nbins), rows sum to <= 1 (beyond-grid tail excluded)
    prior_w: np.ndarray  # (nz,) prior mass per node (sums to 1)
    params: E2Params = field(default_factory=E2Params)

    @property
    def dm_step(self) -> float:
        return float(self.edges[1] - self.edges[0])

    def bin_of(self, dm_eg: np.ndarray) -> np.ndarray:
        k = np.floor(np.asarray(dm_eg, dtype=float) / self.dm_step).astype(int)
        return np.clip(k, 0, self.prob.shape[1] - 1)


def _node_widths(z: np.ndarray) -> np.ndarray:
    mid = 0.5 * (z[1:] + z[:-1])
    lo = np.concatenate([[z[0] - (mid[0] - z[0])], mid])
    hi = np.concatenate([mid, [z[-1] + (z[-1] - mid[-1])]])
    return np.clip(hi - lo, 0, None)


def comoving_volume_prior(z: np.ndarray, p: E2Params = HOFFMANN_EMIN25) -> np.ndarray:
    """Prior mass per z node, uniform in comoving volume: dV/dz ~ D_C^2 / E(z)."""
    cosmo = p.cosmology()
    dc = cosmo.comoving_distance(z).value
    w = dc**2 / cosmo.efunc(z) * _node_widths(z)
    return w / w.sum()


def build_table(
    p: E2Params = HOFFMANN_EMIN25,
    zmin: float = 0.01,
    zmax: float = 4.0,
    nz: int = 600,
    dm_step: float = 2.0,
    dm_max: float = 20000.0,
) -> LikelihoodTable:
    """Tabulate P(DM_EG | z) = cosmic (*) host on a geometric z grid."""
    z = np.geomspace(zmin, zmax, nz)
    edges = np.arange(0.0, dm_max + dm_step, dm_step)
    cos = np.stack([cosmic_bin_probs(zi, edges, p) for zi in z])
    host = np.stack([host_bin_probs(zi, edges, p) for zi in z])
    conv = signal.fftconvolve(cos, host, axes=1)[:, : edges.size - 1]
    conv = np.clip(conv, 0.0, None)  # FFT round-off; rows sum to <= 1 (tail kept out)
    return LikelihoodTable(
        z=z, edges=edges, prob=conv, prior_w=comoving_volume_prior(z, p), params=p
    )


# --- posterior, PIT, intervals -------------------------------------------------------


def posterior(table: LikelihoodTable, dm_eg: np.ndarray | float) -> np.ndarray:
    """Posterior mass per z node for each DM_EG; shape (nz,) or (nz, n)."""
    k = table.bin_of(np.atleast_1d(dm_eg))
    w = table.prob[:, k] * table.prior_w[:, None]
    tot = w.sum(axis=0, keepdims=True)
    w = np.divide(w, tot, out=np.zeros_like(w), where=tot > 0)
    return w[:, 0] if np.ndim(dm_eg) == 0 else w


def _cdf_at(table: LikelihoodTable, post: np.ndarray, zq: np.ndarray) -> np.ndarray:
    """Posterior CDF at zq, treating each node's mass as uniform over its cell."""
    w = _node_widths(table.z)
    lo = table.z - 0.5 * w
    post2 = post[:, None] if post.ndim == 1 else post
    zq = np.atleast_1d(zq)
    frac = np.clip((zq[None, :] - lo[:, None]) / np.where(w > 0, w, 1.0)[:, None], 0.0, 1.0)
    return np.sum(post2 * frac, axis=0)


def pit(table: LikelihoodTable, dm_eg: np.ndarray, z_true: np.ndarray) -> np.ndarray:
    """PIT_b = P(z < z_true,b | DM_EG,b) under E2."""
    post = posterior(table, np.atleast_1d(dm_eg))
    return _cdf_at(table, post, np.atleast_1d(z_true))


def central_interval(table: LikelihoodTable, dm_eg: float, q: float) -> tuple[float, float]:
    """Central q credible interval [z_lo, z_hi] for one burst."""
    post = posterior(table, float(dm_eg))
    w = _node_widths(table.z)
    edges = np.concatenate([[table.z[0] - 0.5 * w[0]], table.z + 0.5 * w])
    cdf = np.concatenate([[0.0], np.cumsum(post)])
    lo, hi = (1 - q) / 2, (1 + q) / 2
    return float(np.interp(lo, cdf, edges)), float(np.interp(hi, cdf, edges))


def coverage(pits: np.ndarray, q: float) -> float:
    """Fraction of PIT values inside the central q band."""
    pits = np.asarray(pits)
    lo, hi = (1 - q) / 2, (1 + q) / 2
    return float(np.mean((pits >= lo) & (pits <= hi)))


# --- the calibration rule and its controls --------------------------------------------


@dataclass
class CalibrationRule:
    """Pass iff coverage at each q lies in its band AND KS(PIT vs U(0,1)) p > ks_p."""

    bands: dict[float, tuple[float, float]]
    ks_p: float = 0.01

    def evaluate(self, pits: np.ndarray) -> dict:
        cov = {q: coverage(pits, q) for q in self.bands}
        in_band = {q: self.bands[q][0] <= cov[q] <= self.bands[q][1] for q in self.bands}
        ks = float(stats.kstest(np.asarray(pits), "uniform").pvalue)
        med = float(np.median(pits))
        passed = all(in_band.values()) and ks > self.ks_p
        if passed:
            verdict = "CALIBRATED"
        elif any(cov[q] < self.bands[q][0] for q in self.bands):
            verdict = "OVERCONFIDENT"
        elif any(cov[q] > self.bands[q][1] for q in self.bands):
            verdict = "UNDERCONFIDENT"
        elif abs(med - 0.5) > 0.1:
            verdict = "BIASED"
        else:
            verdict = "KS_FAIL"  # in-band coverage, KS fails, |median-0.5| <= 0.1: unnamed
        return {
            "coverage": {str(q): cov[q] for q in cov},
            "ks_p": ks,
            "median_pit": med,
            "passed": passed,
            "verdict": verdict,
        }


def simulate_sample(
    gen: LikelihoodTable, n: int, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    """Draw n (z_true, DM_EG) pairs from a table's generative model: z from the prior
    (uniform within the node's cell), DM_EG from that node's bin probabilities
    (uniform within the bin)."""
    i = rng.choice(gen.z.size, size=n, p=gen.prior_w)
    w = _node_widths(gen.z)
    z = gen.z[i] + (rng.random(n) - 0.5) * w[i]
    cum = np.cumsum(gen.prob[i], axis=1)
    cum /= cum[:, -1:]  # draw within the grid (beyond-grid tail mass <~ 0.2%)
    k = (cum < rng.random(n)[:, None]).sum(axis=1)
    k = np.minimum(k, gen.prob.shape[1] - 1)
    dm = gen.edges[k] + rng.random(n) * gen.dm_step
    return z, dm


def _sim_pits(
    score: LikelihoodTable, gen: LikelihoodTable, n: int, m: int, rng: np.random.Generator
) -> np.ndarray:
    out = np.empty((m, n))
    for r in range(m):
        z, dm = simulate_sample(gen, n, rng)
        out[r] = pit(score, dm, z)
    return out


def null_rule(
    null_pits: np.ndarray,
    qs: Iterable[float] = (0.68, 0.95),
    band: float = 0.95,
    ks_p: float = 0.01,
) -> CalibrationRule:
    """Bands = central ``band`` quantile range of each coverage's null distribution."""
    lo, hi = (1 - band) / 2, (1 + band) / 2
    bands = {}
    for q in qs:
        cov = np.array([coverage(row, q) for row in null_pits])
        bands[q] = (float(np.quantile(cov, lo)), float(np.quantile(cov, hi)))
    return CalibrationRule(bands=bands, ks_p=ks_p)


def _fail_rate(rule: CalibrationRule, pits: np.ndarray) -> float:
    return float(np.mean([not rule.evaluate(row)["passed"] for row in pits]))


def run_controls(
    n: int,
    m: int = 2000,
    seed: int = 0,
    max_false_fail: float = 0.07,
    table_kw: dict | None = None,
) -> dict:
    """C0 (null false-fail rate) and C1 (power) as frozen in plan 99.

    C0: bands from m null samples; false-fail measured on an independent m. If it exceeds
    ``max_false_fail`` the band is widened from the central 95% to the central 99% of the
    null (the plan's "97.5%" wording read as the next standard step; logged in findings).
    C1: data generated with host mean x2 (lmean + log10 2) or F x2, scored unmodified.
    """
    kw = table_kw or {}
    rng = np.random.default_rng(seed)
    base = build_table(HOFFMANN_EMIN25, **kw)
    null_a = _sim_pits(base, base, n, m, rng)
    null_b = _sim_pits(base, base, n, m, rng)
    rule = null_rule(null_a, band=0.95)
    ff = _fail_rate(rule, null_b)
    widened = False
    if ff > max_false_fail:
        rule = null_rule(null_a, band=0.99)
        ff = _fail_rate(rule, null_b)
        widened = True
    cases = {
        "host_mean_x2": E2Params(lmean=HOFFMANN_EMIN25.lmean + np.log10(2.0)),
        "F_x2": E2Params(F=2 * HOFFMANN_EMIN25.F),
    }
    power: dict[str, dict[str, Any]] = {}
    for name, gp in cases.items():
        gen = build_table(gp, **kw)
        pits = _sim_pits(base, gen, n, m, rng)
        verdicts = [rule.evaluate(row)["verdict"] for row in pits]
        power[name] = {
            "power": float(np.mean([v != "CALIBRATED" for v in verdicts])),
            "verdict_counts": {v: verdicts.count(v) for v in sorted(set(verdicts))},
            "median_coverage_68": float(np.median([coverage(r, 0.68) for r in pits])),
            "median_coverage_95": float(np.median([coverage(r, 0.95) for r in pits])),
        }
    null_med = {q: float(np.median([coverage(r, q) for r in null_b])) for q in (0.68, 0.95)}
    return {
        "n": n,
        "m": m,
        "seed": seed,
        "table": {
            "zmin": kw.get("zmin", 0.01),
            "zmax": kw.get("zmax", 4.0),
            "nz": kw.get("nz", 600),
            "dm_step": kw.get("dm_step", 2.0),
        },
        "rule": {
            "bands": {str(q): list(b) for q, b in rule.bands.items()},
            "ks_p": rule.ks_p,
            "band_quantile": 0.99 if widened else 0.95,
            "widened": widened,
        },
        "C0": {
            "false_fail_rate": ff,
            "gate": max_false_fail,
            "passed": ff <= max_false_fail,
            "null_median_coverage": {str(q): v for q, v in null_med.items()},
            "null_ks_reject_rate": float(
                np.mean([float(stats.kstest(r, "uniform").pvalue) <= rule.ks_p for r in null_b])
            ),
        },
        "C1": {
            **power,
            "required_power_host": 0.8,
            "passed": power["host_mean_x2"]["power"] >= 0.8,
        },
    }


# --- sample assembly (plan 99 step 0) -------------------------------------------------


def _key(name: str) -> str:
    """'FRB20240210A' -> '20240210A'."""
    return re.sub(r"^FRB", "", name.strip())


def _date(name: str) -> str:
    m = re.match(r"(\d{8})", _key(name))
    if m is None:
        raise ValueError(f"no date in FRB name {name!r}")
    return m.group(1)


def _dm(js: dict) -> float | None:
    v = js.get("DM")
    if isinstance(v, dict):
        v = v.get("value")
    return None if v is None else float(v)


def _n_nonnull(js: dict) -> int:
    return sum(v not in (None, "", [], {}) for v in js.values())


def dedupe(jsons: list[dict], dm_tol: float = 1.0) -> tuple[list[dict], list[dict]]:
    """Group by UTC date; merge only files whose DMs agree within ``dm_tol`` (two genuine
    bursts can share a date). Keep the file with more non-null fields. Each dict needs a
    ``_file`` key holding the burst name."""
    groups: dict[str, list[dict]] = {}
    for js in jsons:
        groups.setdefault(_date(js["_file"]), []).append(js)
    kept, merges = [], []
    for date, g in sorted(groups.items()):
        clusters: list[list[dict]] = []
        for js in g:
            for c in clusters:
                a, b = _dm(js), _dm(c[0])
                if a is not None and b is not None and abs(a - b) < dm_tol:
                    c.append(js)
                    break
            else:
                clusters.append([js])
        for c in clusters:
            best = max(c, key=lambda j: (_n_nonnull(j), j["_file"]))
            kept.append(best)
            if len(c) > 1:
                merges.append(
                    {
                        "date": date,
                        "kept": best["_file"],
                        "dropped": sorted(j["_file"] for j in c if j is not best),
                    }
                )
    return kept, merges


def provenance_side(
    telescope: str,
    name: str,
    dsa_fit: Iterable[str],
    dsa_fit_z_used: Iterable[str],
    askap_cutoff: str,
) -> tuple[str, str]:
    """Production / certification / no_model for one burst against HoffmannEmin25."""
    k = _key(name)
    dsa_fit, dsa_z = set(dsa_fit), set(dsa_fit_z_used)
    if telescope == "DSA":
        if k in dsa_fit:
            return "production", "HoffmannEmin25 DSA table" + (
                " (z used)" if k in dsa_z else " (DM only)"
            )
        return "certification", "DSA, not in HoffmannEmin25 fit"
    if telescope == "ASKAP":
        if _date(name) <= askap_cutoff:
            return "production", "CRAFT <= 2023 (HoffmannEmin25 Sec. 3.3)"
        return "certification", "ASKAP 2024+ (after HoffmannEmin25 sample)"
    if telescope in ("CHIME", "MeerKAT"):
        return "certification", f"{telescope} excluded from fit (HoffmannEmin25 Sec. 2.1)"
    return "no_model", f"{telescope or 'unknown'}: no zdm survey model"
