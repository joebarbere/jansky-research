"""In-band spectra of LPT pulses from ASKAP Taylor-term images (plan 91).

Rose et al. (2026, Nature Astron. 10, 1166) report that ASKAP J174508.9-505149's bursts
"drift in emission frequency". Frequency structure inside a single pulse is otherwise
unmeasured for the class from imaging surveys, because everyone quotes band-averaged flux.

ASKAP images continuum with a multi-frequency-synthesis Taylor expansion, so a spectral index
is available from images the archive already holds:

    I(nu) = T0 + T1 * x + ...,   x = (nu - nu0) / nu0,   alpha ~ T1 / T0

with ``taylor.0`` and ``taylor.1`` served side by side on CASDA. No sub-band re-imaging is
needed. What *is* needed is honesty about when that ratio means anything.

**The ratio is the whole problem.** ``alpha = T1/T0`` divides one noisy quantity by another.
Two consequences the estimator cannot escape:

1. **The lever arm is short.** RACS spans ~288 MHz at ~887 MHz, so ``x`` only ranges about
   +/-0.16. Fitting a slope over that range amplifies noise: the Taylor-1 image carries
   ``sigma_1 ~ sigma_0 / rms(x)``, roughly an order of magnitude worse than Taylor-0. The
   resulting ``sigma_alpha ~ sigma_1 / T0`` therefore scales as ``1 / (S/N)`` with a large
   constant, which is why published guidance puts usable alpha at S/N of order 50-100 rather
   than the ~5 that suffices for a detection.
2. **A noisy denominator biases the ratio**, it does not merely broaden it. E[T1/T0] is not
   T1/T0 when T0 scatters. :func:`injection_recovery` measures that bias rather than assuming
   it away, because a bias is invisible to an error bar.

Nothing here fetches data; :func:`alpha_uncertainty` and the injection harness are pure and
tested, and they decide whether fetching is worth doing at all.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = [
    "RACS_BANDWIDTH_HZ",
    "RACS_NU0_HZ",
    "AlphaMeasurement",
    "RecoveryStats",
    "MT_MFS_PENALTY",
    "alpha_uncertainty",
    "realistic_sigma_alpha",
    "injection_recovery",
    "taylor1_noise_ratio",
    "fetch_taylor_cutout",
    "taylor_alpha",
    "taylor_science_mask",
    "usable_snr_threshold",
]

#: RACS-low reference frequency and bandwidth (McConnell et al. 2020).
RACS_NU0_HZ = 887.5e6
RACS_BANDWIDTH_HZ = 288.0e6

#: How much worse real MT-MFS alpha is than the idealised Gaussian calculation here.
#:
#: Rashid et al. (2024, arXiv:2405.18978) simulated point sources through the MT-MFS
#: algorithm -- the same Taylor-expansion synthesis ASKAPsoft uses -- and found in-band
#: indices "for SNR <~ 100 that have errors >~ 0.2, making them unreliable". The idealised
#: formula here gives 0.107 at S/N 100, so reality is about twice as noisy: deconvolution
#: error, beam and bandpass systematics, and a spectrum that is not exactly a power law all
#: enter. Quoting the idealised number alone would understate every uncertainty by ~2x.
#:
#: This factor is a prior from a uGMRT simulation, NOT an ASKAP calibration. The plan's
#: injection-recovery study on real RACS cutouts must replace it with a measured value.
MT_MFS_PENALTY = 2.0


@dataclass(frozen=True)
class AlphaMeasurement:
    """One in-band spectral index with the uncertainty that decides whether to quote it."""

    alpha: float
    sigma_alpha: float
    snr_taylor0: float
    usable: bool


@dataclass(frozen=True)
class RecoveryStats:
    """Injection-recovery outcome at one signal-to-noise ratio."""

    snr: float
    alpha_true: float
    n_trials: int
    alpha_mean: float
    alpha_median: float
    alpha_std: float
    bias: float
    frac_within_0p3: float


def taylor1_noise_ratio(
    bandwidth_hz: float = RACS_BANDWIDTH_HZ, nu0_hz: float = RACS_NU0_HZ
) -> float:
    """``sigma_1 / sigma_0``: how much noisier the Taylor-1 image is than Taylor-0.

    Fitting ``I = T0 + T1 x`` over channels spread across the band gives
    ``sigma_1 = sigma_0 / rms(x)``, and for a band uniformly covering
    ``x in [-B/2nu0, +B/2nu0]`` the rms is ``(B/nu0)/sqrt(12)``. For RACS-low that is
    ~0.094, so the Taylor-1 image is about 11x noisier -- the reason a 5-sigma detection is
    nowhere near enough for a spectral index.
    """
    frac = bandwidth_hz / nu0_hz
    rms_x = frac / math.sqrt(12.0)
    return 1.0 / rms_x


def taylor_alpha(t0: float, t1: float) -> float:
    """In-band spectral index from the two Taylor coefficients (``alpha = T1/T0``)."""
    if t0 == 0:
        raise ValueError("taylor.0 flux is zero; alpha is undefined")
    return t1 / t0


def alpha_uncertainty(
    t0: float,
    sigma0: float,
    t1: float = 0.0,
    *,
    bandwidth_hz: float = RACS_BANDWIDTH_HZ,
    nu0_hz: float = RACS_NU0_HZ,
) -> float:
    """First-order uncertainty on ``alpha = T1/T0``.

    Standard ratio propagation, ``(sigma_alpha/alpha)^2 = (sigma_1/T1)^2 + (sigma_0/T0)^2``,
    rearranged so it stays finite as ``T1 -> 0``:

        sigma_alpha^2 = (sigma_1/T0)^2 + (alpha * sigma_0/T0)^2

    The first term dominates for the modest indices these sources show, giving
    ``sigma_alpha ~ 11 / (S/N)`` for RACS -- i.e. S/N ~ 110 for ``sigma_alpha = 0.1``.
    Valid only while ``T0 >> sigma0``; when it is not, the linearisation fails and
    :func:`injection_recovery` is the honest tool.
    """
    if t0 <= 0 or sigma0 <= 0:
        raise ValueError("need positive taylor.0 flux and noise")
    sigma1 = sigma0 * taylor1_noise_ratio(bandwidth_hz, nu0_hz)
    alpha = t1 / t0
    return math.hypot(sigma1 / t0, alpha * sigma0 / t0)


def usable_snr_threshold(
    target_sigma_alpha: float = 0.3,
    *,
    bandwidth_hz: float = RACS_BANDWIDTH_HZ,
    nu0_hz: float = RACS_NU0_HZ,
    realistic: bool = True,
) -> float:
    """S/N in Taylor-0 needed to reach ``target_sigma_alpha`` (ignoring the alpha term).

    ``realistic=True`` (the default) applies :data:`MT_MFS_PENALTY`, so the answer is what
    the published MT-MFS behaviour implies rather than what the idealised Gaussian
    calculation promises. Pass ``realistic=False`` for the idealised floor.
    """
    if target_sigma_alpha <= 0:
        raise ValueError("target_sigma_alpha must be positive")
    penalty = MT_MFS_PENALTY if realistic else 1.0
    return penalty * taylor1_noise_ratio(bandwidth_hz, nu0_hz) / target_sigma_alpha


def realistic_sigma_alpha(
    t0: float,
    sigma0: float,
    t1: float = 0.0,
    *,
    bandwidth_hz: float = RACS_BANDWIDTH_HZ,
    nu0_hz: float = RACS_NU0_HZ,
) -> float:
    """:func:`alpha_uncertainty` scaled by the published MT-MFS penalty.

    Use this to decide whether a pulse is worth fetching; use the idealised value only to
    state the floor that no processing can beat.
    """
    return MT_MFS_PENALTY * alpha_uncertainty(
        t0, sigma0, t1, bandwidth_hz=bandwidth_hz, nu0_hz=nu0_hz
    )


def injection_recovery(
    snr: float,
    alpha_true: float,
    *,
    n_trials: int = 20000,
    seed: int = 0,
    bandwidth_hz: float = RACS_BANDWIDTH_HZ,
    nu0_hz: float = RACS_NU0_HZ,
) -> RecoveryStats:
    """Inject a source of known ``alpha_true`` at ``snr`` and recover it ``n_trials`` times.

    Measures **bias as well as scatter**. A ratio with a noisy denominator is biased, and an
    error bar cannot show that; this is the same shape as the repo's standing lesson that a
    resampling test bounds variance and is silent on bias.

    ``frac_within_0p3`` is the fraction of trials recovering alpha to within 0.3 -- roughly
    the precision at which an in-band index distinguishes the physically interesting cases
    (steep aged plasma vs flat coherent emission).
    """
    if snr <= 0:
        raise ValueError("snr must be positive")
    rng = np.random.default_rng(seed)
    sigma0 = 1.0
    t0_true = snr * sigma0
    t1_true = alpha_true * t0_true
    sigma1 = sigma0 * taylor1_noise_ratio(bandwidth_hz, nu0_hz)

    t0 = t0_true + rng.normal(0.0, sigma0, n_trials)
    t1 = t1_true + rng.normal(0.0, sigma1, n_trials)
    # A non-detection in Taylor-0 has no defined index; drop those rather than let a
    # near-zero denominator manufacture a spectacular alpha.
    ok = t0 > 5.0 * sigma0
    rec = t1[ok] / t0[ok]

    if rec.size == 0:
        return RecoveryStats(snr, alpha_true, 0, math.nan, math.nan, math.nan, math.nan, 0.0)
    return RecoveryStats(
        snr=snr,
        alpha_true=alpha_true,
        n_trials=int(rec.size),
        alpha_mean=float(rec.mean()),
        alpha_median=float(np.median(rec)),
        alpha_std=float(rec.std(ddof=1)),
        bias=float(rec.mean() - alpha_true),
        frac_within_0p3=float(np.mean(np.abs(rec - alpha_true) <= 0.3)),
    )


def measure(
    t0: float,
    sigma0: float,
    t1: float,
    *,
    target_sigma_alpha: float = 0.3,
    bandwidth_hz: float = RACS_BANDWIDTH_HZ,
    nu0_hz: float = RACS_NU0_HZ,
) -> AlphaMeasurement:
    """Alpha, its uncertainty, and whether it clears the usability bar."""
    alpha = taylor_alpha(t0, t1)
    sig = alpha_uncertainty(t0, sigma0, t1, bandwidth_hz=bandwidth_hz, nu0_hz=nu0_hz)
    return AlphaMeasurement(
        alpha=alpha,
        sigma_alpha=sig,
        snr_taylor0=t0 / sigma0,
        usable=sig <= target_sigma_alpha,
    )


def taylor_science_mask(table, term: int, sbid: str, stokes: str = "i"):
    """Select the restored, convolved Taylor-``term`` science image for one SBID.

    Two filters the general RACS cutout helper does not apply, both essential here:

    * **the SBID**, because a pulse exists in one observation and the same position has
      images from many epochs -- fetching "an image here" would measure a random epoch in
      which the source is almost certainly off;
    * **the Taylor term**, since ``taylor.0`` and ``taylor.1`` sit side by side and differ
      only in that field.

    ``noiseMap``/``meanMap`` products carry the same substrings, so ``restored`` and ``conv``
    are required as in :func:`stokesv._racs_science_mask`.
    """
    import numpy as _np

    fn = _np.array([str(x) for x in table["filename"]])
    want = f"taylor.{term}"
    return _np.array(
        [
            f.startswith(f"image.{stokes}.")
            and f"SB{sbid}" in f
            and want in f
            and "restored" in f
            and "conv" in f
            for f in fn
        ]
    )


def fetch_taylor_cutout(
    ra: float,
    dec: float,
    sbid: str,
    term: int,
    *,
    stokes: str = "i",
    radius_deg: float = 0.03,
    casda=None,
    username: str | None = None,
    pw_path: str = "~/.casda_pw",
    retries: int = 3,
):  # pragma: no cover - network
    """Stage and read one Taylor-term cutout from CASDA -> ``(image_mJy, wcs, casda)``.

    Mirrors ``stokesv.fetch_racs_cutout`` (same OPAL login, SODA staging and retry-on-401)
    but selects by SBID and Taylor term via :func:`taylor_science_mask`.
    """
    import os
    import tempfile

    import astropy.units as _u
    import numpy as _np
    import requests
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.wcs import WCS
    from astroquery.casda import Casda

    from .stokesv import _casda_session

    username = username or os.environ.get("CASDA_USERNAME")
    if not username:
        raise RuntimeError("set CASDA_USERNAME (OPAL email) for the Taylor-term fetch")
    coord = SkyCoord(ra * _u.deg, dec * _u.deg)
    for _ in range(retries):
        try:
            if casda is None:
                casda = _casda_session(username, pw_path)
            table = Casda.query_region(coord, radius=0.1 * _u.deg)
            mask = taylor_science_mask(table, term, sbid, stokes)
            if not mask.any():
                return None
            urls = casda.cutout(table[mask][:1], coordinates=coord, radius=radius_deg * _u.deg)
            furl = next(u for u in urls if u.endswith(".fits"))
            raw = requests.get(furl, timeout=300).content
            with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as fh:
                fh.write(raw)
                path = fh.name
            with fits.open(path) as hd:
                data = _np.squeeze(_np.asarray(hd[0].data, float))
                wcs = WCS(hd[0].header).celestial
            os.unlink(path)
            return data * 1000.0, wcs, casda  # Jy/beam -> mJy/beam
        except Exception:
            casda = None
    return None
