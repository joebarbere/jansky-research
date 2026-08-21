"""How often are long-period transients actually on? (plan 90)

Rose et al. (2026, Nature Astron. 10, 1166) report that the bursts of ASKAP
J174508.9-505149 "turn off for several hours at a time" -- qualitatively, for one source.
This module turns the class version of that statement into a number, using the VAST snapshot
table the ``lptv`` slice already committed (``results/lptv_vast_epochs.csv``): per-epoch MJD,
snapshot ``duration_s``, forced Stokes-I/V flux and the per-epoch noise.

Two things decide whether the answer means anything, and both are enforced here rather than
left to the write-up:

**Only a product is identifiable.** A snapshot of length *T* on a source of period *P* whose
emitting phase has width *w* catches a pulse with probability ~min(1, (w + T)/P) *when the
source is active at all*. A detection rate therefore constrains

    p = f_active * (w + T) / P

and not its factors. :func:`duty_constraint` returns *p*; splitting it needs phase
information from an ephemeris, which most of these sources do not have.

**A null divides by sensitivity, not by epoch count.** This is the ``frblens`` lesson (that
slice's limit was 4x too tight until per-source efficiency entered the denominator). A pulse
of flux *S* is not detectable in every epoch: each has its own noise. The denominator here is
therefore the summed per-epoch detection efficiency at an assumed pulse flux, so epochs too
shallow to have seen the pulse contribute ~0 instead of padding the exposure. Every limit is
consequently a *function of the assumed pulse flux*, and is reported that way.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = [
    "Ephemeris",
    "EpochRow",
    "PhaseResolved",
    "SourceConstraint",
    "DETECT_THRESHOLD_SIGMA",
    "MIN_EFFICIENCY",
    "duty_constraint",
    "epoch_efficiency",
    "PhaseSampling",
    "load_epochs",
    "phase_resolved_activity",
    "phase_sampling",
    "read_period_precision",
    "read_periods",
]

#: Detection threshold used by the ``lptv`` VAST sweep, in sigma. Kept explicit because the
#: efficiency (and therefore every limit) depends on it.
DETECT_THRESHOLD_SIGMA = 5.0

#: 95% upper limit on the mean of a Poisson process with zero events (Gehrels 1986).
POISSON_ZERO_95 = 2.996

#: Efficiencies below this are treated as zero. A Gaussian tail is never exactly 0, so a
#: 0.5 mJy pulse in a 100 mJy epoch still scores ~3e-7 -- and summing enough such epochs
#: manufactures sensitivity that does not exist (10^6 of them would look like a third of an
#: epoch). An epoch that could not plausibly have seen the pulse must contribute nothing.
MIN_EFFICIENCY = 1e-3


@dataclass(frozen=True)
class EpochRow:
    """One forced-photometry snapshot from the committed VAST sweep."""

    name: str
    obs_id: str
    epoch_mjd: float
    duration_s: float
    v_mjy: float
    e_v: float


@dataclass(frozen=True)
class SourceConstraint:
    """Per-source duty-cycle constraint at one assumed pulse flux."""

    name: str
    n_epochs: int
    total_exposure_s: float
    effective_epochs: float
    n_detections: int
    assumed_pulse_mjy: float
    p_point: float | None
    p_upper_95: float
    identifiable_factors: bool


def load_epochs(path: str | Path) -> list[EpochRow]:
    """Read the committed VAST epoch table, keeping only rows that carry a measurement.

    Rows with an empty ``v_mjy``/``e_v`` are epochs the archive never released (the table
    marks them in ``note``); they measured nothing and must not enter any denominator.
    """
    out: list[EpochRow] = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                v, e = float(row["v_mjy"]), float(row["e_v"])
                dur = float(row["duration_s"])
                mjd = float(row["epoch_mjd"])
            except (TypeError, ValueError):
                continue
            if not math.isfinite(v) or not math.isfinite(e) or e <= 0:
                continue
            out.append(
                EpochRow(
                    name=row["name"],
                    obs_id=row["obs_id"],
                    epoch_mjd=mjd,
                    duration_s=dur,
                    v_mjy=v,
                    e_v=e,
                )
            )
    return out


def read_period_precision(path: str | Path) -> dict[str, float]:
    """``{source: quoted precision in seconds}`` from the catalogue's decimal places.

    A proxy for, NOT the same as, the published uncertainty: a period written 1318.1957 is
    quoted to 1e-4 s, but its paper may still quote +/-5e-4. Use this to decide whether a
    phase test is even worth attempting; the GATE-0 ephemeris audit must replace it with the
    real uncertainties before any phase-resolved claim.
    """
    out: dict[str, float] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            txt = (row.get("period_s") or "").strip()
            if not txt:
                continue
            try:
                float(txt)
            except ValueError:
                continue
            decimals = len(txt.split(".")[1]) if "." in txt else 0
            out[row["name"]] = 10.0**-decimals
    return out


def read_periods(path: str | Path) -> dict[str, float]:
    """``{source name: period in seconds}`` from the vendored LPT catalogue."""
    periods: dict[str, float] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                periods[row["name"]] = float(row["period_s"])
            except (TypeError, ValueError, KeyError):
                continue
    return periods


def _phi(x: float) -> float:
    """Standard normal CDF (no SciPy dependency for one function)."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def epoch_efficiency(
    sigma_mjy: float,
    pulse_mjy: float,
    *,
    threshold_sigma: float = DETECT_THRESHOLD_SIGMA,
    min_efficiency: float = MIN_EFFICIENCY,
) -> float:
    """Probability that a pulse of ``pulse_mjy`` clears the threshold in this epoch.

    Gaussian noise of width ``sigma_mjy`` on a forced measurement: the measured value exceeds
    ``threshold_sigma * sigma`` with probability ``Phi(S/sigma - threshold)``. An epoch whose
    noise is far above the pulse flux returns ~0 and so contributes nothing to the exposure --
    which is the entire point.
    """
    if sigma_mjy <= 0 or not math.isfinite(sigma_mjy):
        return 0.0
    eff = _phi(pulse_mjy / sigma_mjy - threshold_sigma)
    return eff if eff >= min_efficiency else 0.0


def duty_constraint(
    rows: list[EpochRow],
    *,
    pulse_mjy: float,
    threshold_sigma: float = DETECT_THRESHOLD_SIGMA,
    min_efficiency: float = MIN_EFFICIENCY,
    identifiable_factors: bool = False,
) -> SourceConstraint:
    """Constrain ``p = f_active * (w + T) / P`` for one source at an assumed pulse flux.

    ``rows`` must all belong to one source. Detections are counted as epochs whose measured
    ``|V|/sigma`` clears ``threshold_sigma`` -- the same criterion the ``lptv`` sweep applied.

    The effective exposure is ``sum(epoch_efficiency)``, not ``len(rows)``. With ``k``
    detections the point estimate is ``k / effective_epochs``; with zero it is ``None`` and
    only the 95% upper limit ``2.996 / effective_epochs`` is reported. If no epoch could have
    seen a pulse this faint the effective exposure is ~0, the limit is infinite, and that is
    the honest answer -- not a limit of 1.
    """
    if not rows:
        raise ValueError("no epochs supplied")
    names = {r.name for r in rows}
    if len(names) != 1:
        raise ValueError(f"rows span more than one source: {sorted(names)}")

    eff = np.array(
        [epoch_efficiency(r.e_v, pulse_mjy, threshold_sigma=threshold_sigma) for r in rows]
    )
    effective = float(eff.sum())
    k = sum(1 for r in rows if abs(r.v_mjy) / r.e_v >= threshold_sigma)

    if effective <= 0:
        p_point: float | None = None
        p_upper = math.inf
    else:
        p_point = k / effective if k else None
        p_upper = (POISSON_ZERO_95 if k == 0 else POISSON_ZERO_95 + k) / effective

    return SourceConstraint(
        name=rows[0].name,
        n_epochs=len(rows),
        total_exposure_s=float(sum(r.duration_s for r in rows)),
        effective_epochs=effective,
        n_detections=k,
        assumed_pulse_mjy=pulse_mjy,
        p_point=p_point,
        p_upper_95=p_upper,
        identifiable_factors=identifiable_factors,
    )


@dataclass(frozen=True)
class PhaseSampling:
    """Whether a source's snapshots sample its pulse phase uniformly."""

    name: str
    n_epochs: int
    period_s: float
    baseline_days: float
    rayleigh_z: float
    rayleigh_p: float
    kuiper_v: float
    required_period_precision_s: float
    coherent_assumption_needed: bool


def phase_sampling(rows: list[EpochRow], period_s: float) -> PhaseSampling:
    """Test whether snapshots sample pulse phase uniformly, for one source.

    The binomial constraint in :func:`duty_constraint` assumes each snapshot is an
    independent draw on pulse phase. VAST observes on a roughly fortnightly cadence, which
    is not random with respect to any LPT period, so that assumption has to be checked
    rather than asserted.

    Phases are ``(mjd * 86400 / period) mod 1``. **The zero point is arbitrary and does not
    matter**: clustering is invariant under a phase shift, so no ephemeris is needed to test
    *uniformity* (an ephemeris is needed only to assign physically meaningful phase).

    Two statistics, because they fail differently: Rayleigh Z catches a single concentration,
    Kuiper V is sensitive to any departure from uniformity including multimodal ones.

    ``required_period_precision_s`` is the period accuracy needed for phase to stay coherent
    to 0.1 cycles across the observing baseline: ``0.1 * P^2 / T_baseline``. If the catalogued
    period is less precise than this, phases smear and **this test cannot detect clustering
    that may still be present** -- a null result here is then uninformative, not reassuring.
    """
    if period_s <= 0 or not math.isfinite(period_s):
        raise ValueError(f"period must be positive and finite, got {period_s}")
    if len(rows) < 2:
        raise ValueError("need at least two epochs to test phase sampling")

    mjd = np.array([r.epoch_mjd for r in rows], dtype=float)
    baseline_days = float(mjd.max() - mjd.min())
    # Phase relative to the first epoch, not to MJD 0. At MJD ~59000 a 1 h period is ~1.2e6
    # cycles, where a float64 keeps only ~1e-10 of a cycle; referencing to the first epoch
    # drops the magnitude by orders of magnitude and makes the zero-point invariance exact
    # rather than approximate. The zero point is arbitrary either way.
    phase = np.mod((mjd - mjd.min()) * 86400.0 / period_s, 1.0)
    n = len(phase)

    ang = 2.0 * math.pi * phase
    c, s = float(np.cos(ang).sum()), float(np.sin(ang).sum())
    r_bar = math.hypot(c, s) / n
    z = n * r_bar**2
    # Standard small-sample correction (Mardia & Jupp); p ~ exp(-Z) to leading order.
    p_rayleigh = math.exp(-z) * (1.0 + (2.0 * z - z**2) / (4.0 * n))
    p_rayleigh = min(max(p_rayleigh, 0.0), 1.0)

    srt = np.sort(phase)
    i = np.arange(1, n + 1, dtype=float)
    d_plus = float((i / n - srt).max())
    d_minus = float((srt - (i - 1.0) / n).max())
    kuiper_v = d_plus + d_minus

    required = 0.1 * period_s**2 / (baseline_days * 86400.0) if baseline_days > 0 else math.inf

    return PhaseSampling(
        name=rows[0].name,
        n_epochs=n,
        period_s=period_s,
        baseline_days=baseline_days,
        rayleigh_z=z,
        rayleigh_p=p_rayleigh,
        kuiper_v=kuiper_v,
        required_period_precision_s=required,
        coherent_assumption_needed=True,
    )


@dataclass(frozen=True)
class Ephemeris:
    """A published timing solution, with the uncertainties that decide if it is usable."""

    name: str
    period_s: float
    sigma_period_s: float | None
    pepoch_mjd: float
    pulse_width_s: float
    reference: str

    def phase_uncertainty_at(self, mjd: float) -> float:
        """Accumulated phase error (cycles) at ``mjd``, from the period uncertainty alone.

        ``(t - PEPOCH) * sigma_P / P^2``. Returns inf when the paper reports no uncertainty:
        an unreported error is not a zero error, and treating it as one is how a phase-folded
        claim goes quietly wrong years after the reference epoch.
        """
        if self.sigma_period_s is None:
            return math.inf
        dt_s = abs(mjd - self.pepoch_mjd) * 86400.0
        return dt_s * self.sigma_period_s / self.period_s**2


@dataclass(frozen=True)
class PhaseResolved:
    """f_active separated from the in-period duty cycle, for one source."""

    name: str
    n_on_window: int
    effective_on_window: float
    n_detections_in_window: int
    n_detections_outside: int
    window_fraction: float
    f_active_point: float | None
    f_active_upper_95: float
    max_phase_uncertainty: float
    usable: bool


def phase_resolved_activity(
    rows: list[EpochRow],
    eph: Ephemeris,
    *,
    pulse_phase: float = 0.0,
    pulse_mjy: float,
    threshold_sigma: float = DETECT_THRESHOLD_SIGMA,
    min_efficiency: float = MIN_EFFICIENCY,
    max_phase_uncertainty: float = 0.1,
) -> PhaseResolved:
    """Split ``p`` into f_active and the in-period duty cycle, where an ephemeris allows it.

    :func:`duty_constraint` can only constrain the product, because a non-detection is
    ambiguous between "the source was off" and "the snapshot missed the pulse". With physical
    phase that ambiguity resolves: restrict to snapshots whose phase coverage overlaps the
    emitting window, and the detection rate *within that subset* estimates f_active directly,
    while ``(w + T)/P`` is computed from published quantities rather than fitted.

    A snapshot of length T covers a phase range T/P wide, so the window it must overlap is
    ``(w + T)/P`` wide -- the same combination that appears in the product.

    ``usable`` is False when accumulated phase error exceeds ``max_phase_uncertainty`` at any
    epoch: past that, snapshots cannot be assigned to the window at all, and the split is not
    available regardless of how many epochs exist. Detections *outside* the window are
    reported rather than discarded -- they are evidence the ephemeris is wrong, not noise.
    """
    if not rows:
        raise ValueError("no epochs supplied")
    names = {r.name for r in rows}
    if len(names) != 1:
        raise ValueError(f"rows span more than one source: {sorted(names)}")

    worst = max(eph.phase_uncertainty_at(r.epoch_mjd) for r in rows)
    usable = worst <= max_phase_uncertainty

    in_window: list[EpochRow] = []
    det_in = det_out = 0
    widths = []
    for r in rows:
        span = (eph.pulse_width_s + r.duration_s) / eph.period_s
        widths.append(span)
        start = (r.epoch_mjd - eph.pepoch_mjd) * 86400.0 / eph.period_s
        # phase of the snapshot's midpoint relative to the pulse
        mid = math.fmod(start + 0.5 * r.duration_s / eph.period_s - pulse_phase, 1.0)
        if mid < 0:
            mid += 1.0
        dist = min(mid, 1.0 - mid)  # circular distance to the pulse phase
        hit = dist <= 0.5 * span
        detected = abs(r.v_mjy) / r.e_v >= threshold_sigma
        if hit:
            in_window.append(r)
            det_in += int(detected)
        elif detected:
            det_out += 1

    eff = sum(
        epoch_efficiency(
            r.e_v, pulse_mjy, threshold_sigma=threshold_sigma, min_efficiency=min_efficiency
        )
        for r in in_window
    )
    if eff <= 0:
        f_point: float | None = None
        f_upper = math.inf
    else:
        f_point = det_in / eff if det_in else None
        f_upper = (POISSON_ZERO_95 if det_in == 0 else POISSON_ZERO_95 + det_in) / eff

    return PhaseResolved(
        name=rows[0].name,
        n_on_window=len(in_window),
        effective_on_window=float(eff),
        n_detections_in_window=det_in,
        n_detections_outside=det_out,
        window_fraction=float(np.mean(widths)) if widths else 0.0,
        f_active_point=f_point,
        f_active_upper_95=f_upper,
        max_phase_uncertainty=worst,
        usable=usable,
    )
