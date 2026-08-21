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
    "EpochRow",
    "SourceConstraint",
    "DETECT_THRESHOLD_SIGMA",
    "MIN_EFFICIENCY",
    "duty_constraint",
    "epoch_efficiency",
    "load_epochs",
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
