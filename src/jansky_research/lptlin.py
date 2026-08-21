"""Linear polarization at LPT pulse positions (plan 92).

Rose et al. (2026, Nature Astron. 10, 1166) describe the bursts of ASKAP J174508.9-505149 as
**elliptically** polarized -- carrying linear as well as circular power. The ``lptv`` census
is Stokes-V only and says so explicitly: its non-detections are "limits on circularly
polarized flux, not on fractional polarization". This module supplies the linear counterpart
where the archive allows.

**Why this survives the failure that killed plan 91.** ``lptspec`` tried to read a spectral
index from ``taylor.1`` and found the MFS model invalid for a transient: a pulse present for
part of the synthesis cannot be represented by one constant-flux power-law source, and the
deconvolution absorbs the mismatch. Here we read ``taylor.0`` -- the band-averaged flux -- in
Q and U, exactly as ``lptv`` already does for I and V. A pulse on for a fraction ``f`` of the
integration is diluted to ``f*I``, ``f*Q``, ``f*U``, ``f*V`` alike, so **dilution cancels in
any ratio of Stokes parameters from the same image**. Fractional polarization is therefore
recoverable even though the absolute flux is not, which is the same reason ``lptv`` can quote
|V|/I for a pulse it only partly integrates.

The cancellation is exact for the raw ratio and **approximate once debiasing is applied**,
because the noise does not dilute with the signal: at ``f = 0.5`` the recovered fraction is
within 0.02% of the undiluted value, at ``f = 0.2`` within 0.12%, but by ``f = 0.05`` -- where
the pulse is close to the noise -- it is 2% low, biased *down* by the debiasing term. So the
fraction is robust to dilution while the pulse stays well above the noise, and becomes a mild
underestimate when it does not. That direction is the safe one; it cannot manufacture
polarization.

Two corrections the ratio still needs:

* **Ricean bias.** ``L = sqrt(Q^2 + U^2)`` is a positive-definite combination of two noisy
  quantities, so it is biased high -- badly so near the noise floor, where it is positive even
  for an unpolarized source. :func:`debias_linear` removes the first-order term.
* **Instrumental leakage.** ASKAP leaks I into Q/U as it does into V, so a fractional
  polarization at or below the leakage floor is not a measurement of the source.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "LEAKAGE_FRAC_QU",
    "LinearPol",
    "debias_linear",
    "linear_fraction",
    "polarization_angle",
    "total_polarization",
]

#: ASKAP on-axis Stokes I -> Q/U leakage floor. Taken equal to the I->V figure the ``lptv``
#: slice uses (``lptv.LEAKAGE_FRAC = 0.006``) because ASKAP's published on-axis leakage is
#: quoted at that level for the linear terms too; a measured per-field value would be better
#: and this is deliberately conservative rather than precise.
LEAKAGE_FRAC_QU = 0.006


@dataclass(frozen=True)
class LinearPol:
    """Linear polarization of one pulse, with the two corrections applied."""

    q_mjy: float
    u_mjy: float
    i_mjy: float
    sigma_qu_mjy: float
    l_raw_mjy: float
    l_debiased_mjy: float
    linear_fraction: float
    angle_deg: float
    significance: float
    above_leakage: bool
    detected: bool


def debias_linear(l_raw: float, sigma_qu: float) -> float:
    """Remove the first-order Ricean bias from ``L = sqrt(Q^2 + U^2)``.

    ``L`` is positive-definite: for an unpolarized source it still returns ~1.25 sigma rather
    than zero, so an uncorrected ``L/I`` invents polarization for every faint source. The
    standard estimator subtracts the noise in quadrature, ``sqrt(L^2 - sigma^2)``, and returns
    zero where that would be imaginary -- which is the honest answer, not a small number.
    """
    if sigma_qu <= 0:
        raise ValueError("sigma_qu must be positive")
    if l_raw <= sigma_qu:
        return 0.0
    return math.sqrt(l_raw**2 - sigma_qu**2)


def total_polarization(l_mjy: float, v_mjy: float, i_mjy: float) -> float:
    """Total polarized fraction ``sqrt(L^2 + V^2)/|I|``, which must not exceed 1.

    A cheap invariant worth checking on every measurement: Stokes parameters obey
    ``Q^2 + U^2 + V^2 <= I^2`` for any physical signal, so a total fraction above 1 means the
    photometry, the leakage, or the image is wrong -- not that the source is exotic.
    """
    if i_mjy == 0:
        raise ValueError("Stokes I is zero; a fractional polarization is undefined")
    return math.hypot(l_mjy, v_mjy) / abs(i_mjy)


def polarization_angle(q: float, u: float) -> float:
    """Electric-vector position angle in degrees, ``0.5*atan2(U, Q)``, wrapped to [0, 180)."""
    ang = 0.5 * math.degrees(math.atan2(u, q))
    return ang % 180.0


def linear_fraction(
    q_mjy: float,
    u_mjy: float,
    i_mjy: float,
    sigma_qu_mjy: float,
    *,
    det_sigma: float = 5.0,
    leakage_frac: float = LEAKAGE_FRAC_QU,
) -> LinearPol:
    """Debiased linear fraction ``L/I`` at one position, with the leakage veto applied.

    ``detected`` requires both that the debiased ``L`` clears ``det_sigma`` and that it
    exceeds the leakage floor ``leakage_frac * |I|`` -- the same two-part criterion ``lptv``
    applies to V, because a significant ``L`` that sits at the leakage level is measuring the
    instrument rather than the source.
    """
    if sigma_qu_mjy <= 0:
        raise ValueError("sigma_qu must be positive")
    if i_mjy == 0:
        raise ValueError("Stokes I is zero; a fractional polarization is undefined")
    l_raw = math.hypot(q_mjy, u_mjy)
    l_deb = debias_linear(l_raw, sigma_qu_mjy)
    sig = l_deb / sigma_qu_mjy
    above = l_deb > leakage_frac * abs(i_mjy)
    return LinearPol(
        q_mjy=q_mjy,
        u_mjy=u_mjy,
        i_mjy=i_mjy,
        sigma_qu_mjy=sigma_qu_mjy,
        l_raw_mjy=l_raw,
        l_debiased_mjy=l_deb,
        linear_fraction=l_deb / abs(i_mjy),
        angle_deg=polarization_angle(q_mjy, u_mjy),
        significance=sig,
        above_leakage=above,
        detected=bool(sig >= det_sigma and above),
    )
