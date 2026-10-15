"""RACS Stokes-V two-epoch census of the Gaia ultracool-dwarf sample (plan 59, slice ``ucdv``).

Forced Stokes-V photometry at the proper-motion-propagated positions of the Reyl\\'e et al. (2018)
Gaia DR2 ultracool dwarfs, in RACS-low2 (887.5 MHz) and RACS-low3 (943.5 MHz). Ultracool-dwarf
radio emission is coherent and often deeply circularly polarized, so V at ~0.9 GHz probes a
different mechanism and band from the closed 3 GHz VLASS Stokes-I search (arXiv:2506.21169).
The expected yield is near zero; the per-target two-epoch V limit table is the product.

Every control below was frozen in ``plans/59-racs-v-ucd-census.md`` before any real target cutout
was opened, and the constants here are those frozen values -- change them only with a dated plan
amendment that says why. Thresholds are production-side; anything that calibrates is certified on
the disjoint set (Gaia ``source_id`` parity, :func:`assign_set`).

This module holds the pure, tested pieces (NumPy/SciPy/astropy, synthetic fixture). Photometry
reuses :func:`jansky_research.stokesv.measure_circular_pol` in its genuinely forced mode only
(``search_arcsec=0``): the 12-arcsec search mode is a noise maximum on blank sky, the
``stokesv_discovery`` lesson.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

from .stokesv import measure_circular_pol

__all__ = [
    "ALPHA_MAX",
    "BRIGHT_GUARD_ARCMIN",
    "BRIGHT_GUARD_MJY",
    "BUDGET_HOURS",
    "CUTOUT_RADIUS_DEG",
    "GAIA_DR2_EPOCH",
    "NU_LOW2_GHZ",
    "NU_LOW3_GHZ",
    "SCALE_TERM",
    "SEED",
    "SNR_DETECT",
    "Beam",
    "apply_detection_rule",
    "assign_set",
    "forced_measure",
    "injection_check",
    "inject_point_source",
    "luminosity_limit",
    "null_check",
    "null_positions",
    "recover_known_check",
    "scope_gate",
    "synthetic_cutout",
    "two_epoch_compare",
    "upper_limit",
    "validate_synthetic",
]

# ---- frozen constants (plans/59, 2026-10-15) ---------------------------------------------------
GAIA_DR2_EPOCH = 2015.5  # Reyle+2018 positions are Gaia DR2, reference epoch J2015.5
NU_LOW2_GHZ = 0.8875
NU_LOW3_GHZ = 0.9435
SEED = 59
SNR_DETECT = 5.0  # |V|/sigma_V at the forced pixel
BRIGHT_GUARD_MJY = 100.0  # no Stokes-I source brighter than this ...
BRIGHT_GUARD_ARCMIN = 3.0  # ... within this radius
SCALE_TERM = 0.10  # per-epoch flux-scale term on any two-epoch V change
ALPHA_MAX = 3.0  # |spectral index| a low2/low3 difference must exceed to be more than spectral
BUDGET_HOURS = 120.0  # S0: projected census wall-clock budget
NULL_OFFSET_ARCMIN = 2.0
NULL_PER_TARGET = 5  # nulls per target-epoch, bearings 360/5 apart (> 2 beams between them)
NULL_SIGN_CI = 0.99  # amended 2026-10-15 (pre-data): 0.95 + one null gave 6.6% false-fail at n~80
CUTOUT_RADIUS_DEG = 0.06  # 2' null + 90" annulus = 210"; the 3' bright guard also needs it
NULL_STD_RANGE = (0.85, 1.20)
NULL_EXCESS_ALLOWED = 3  # n(|V|/sigma >= 5) may exceed the Gaussian expectation by at most this
INJ_LEVELS_SIGMA = (3.0, 5.0, 7.0, 10.0, 20.0)
INJ_COMPLETENESS_AT_7 = 0.90
INJ_FLUX_RATIO_RANGE = (0.95, 1.05)  # median recovered/injected at >= 10 sigma
KNOWN_RATIO_RANGE = (0.85, 1.15)  # C3(a) median forced/catalogue V
LIMIT_NSIGMA = 3.0
CUTOUT_RMS_MAX_FACTOR = (
    5.0  # a cutout whose annulus RMS exceeds this x the epoch median is unusable
)

PC_CM = 3.0856775814913673e18
MJY_CGS = 1e-26  # 1 mJy in erg s^-1 cm^-2 Hz^-1


@dataclass(frozen=True)
class Beam:
    """Restoring beam: FWHM major/minor axes (arcsec) and position angle (deg, E of N)."""

    bmaj_arcsec: float
    bmin_arcsec: float
    bpa_deg: float = 0.0


def assign_set(source_id: np.ndarray) -> np.ndarray:
    """Provenance split by Gaia ``source_id`` parity: even -> ``"A"`` (calibrates), odd -> ``"B"``.

    A deterministic, data-blind rule fixed in the plan: no flux or position enters it, so set B
    can certify what set A was used to tune.
    """
    sid = np.asarray(source_id, dtype=np.int64)
    return np.where(sid % 2 == 0, "A", "B")


def null_positions(
    ra: np.ndarray,
    dec: np.ndarray,
    *,
    offset_arcmin: float = NULL_OFFSET_ARCMIN,
    n_per_target: int = NULL_PER_TARGET,
    seed: int = SEED,
) -> tuple[np.ndarray, np.ndarray]:
    """Blank-sky control positions: ``n_per_target`` points ``offset_arcmin`` from each target.

    Bearings are ``360/n_per_target`` apart with one random rotation per target, so at 2 arcmin
    and n = 5 neighbouring nulls are ~141 arcsec apart -- several beams, i.e. independent. The
    offset is a true angular separation (astropy ``directional_offset_by``). Output arrays have
    ``len(ra) * n_per_target`` entries, target-major (target k's nulls are rows k*n .. k*n+n-1).
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    ra = np.atleast_1d(np.asarray(ra, float))
    dec = np.atleast_1d(np.asarray(dec, float))
    rng = np.random.default_rng(seed)
    rot = rng.uniform(0.0, 360.0, ra.size)
    pa = (rot[:, None] + 360.0 * np.arange(n_per_target)[None, :] / n_per_target) % 360.0
    c = SkyCoord(np.repeat(ra, n_per_target) * u.deg, np.repeat(dec, n_per_target) * u.deg)
    off = c.directional_offset_by(pa.ravel() * u.deg, offset_arcmin * u.arcmin)
    return np.asarray(off.ra.deg), np.asarray(off.dec.deg)


def inject_point_source(
    image: np.ndarray,
    wcs,
    ra: float,
    dec: float,
    flux_mjy: float,
    beam: Beam,
) -> np.ndarray:
    """Return a copy of ``image`` (mJy/beam) with a point source of ``flux_mjy`` at ``(ra, dec)``.

    A point source's peak in Jy/beam equals its flux, so the added component is the restoring-beam
    Gaussian (FWHM ``bmaj``/``bmin``, ``bpa`` east of north) scaled to peak ``flux_mjy`` at the
    exact sub-pixel position. Assumes a north-up, east-left image (RACS cutouts are).
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.wcs.utils import proj_plane_pixel_scales

    img = np.array(image, dtype=float, copy=True)
    px, py = wcs.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    scale = float(np.mean(proj_plane_pixel_scales(wcs)) * 3600.0)  # arcsec/pixel
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    # sky offsets in arcsec: +x pixel is west in a north-up east-left image, so east = -dx
    d_east = -(xx - float(px)) * scale
    d_north = (yy - float(py)) * scale
    theta = np.radians(beam.bpa_deg)
    # rotate into the beam frame: major axis at PA theta (east of north)
    u_maj = d_north * np.cos(theta) + d_east * np.sin(theta)
    u_min = -d_north * np.sin(theta) + d_east * np.cos(theta)
    k = 4.0 * np.log(2.0)
    g = np.exp(-k * ((u_maj / beam.bmaj_arcsec) ** 2 + (u_min / beam.bmin_arcsec) ** 2))
    return img + flux_mjy * g


def forced_measure(image_i, image_v, wcs, ra: float, dec: float) -> dict[str, float]:
    """Forced I/V photometry at the pixel containing ``(ra, dec)`` -- never the search mode."""
    return measure_circular_pol(image_i, image_v, wcs, ra, dec, search_arcsec=0.0)


def null_check(v: np.ndarray, sigma_v: np.ndarray) -> dict:
    """C1: is forced V on blank sky signed, unit-variance and free of a 5-sigma excess?

    Pass criteria (frozen): the fraction of V > 0 lies inside the binomial ``NULL_SIGN_CI``
    interval of 0.5;
    ``std(V/sigma_V)`` is in ``NULL_STD_RANGE``; the count with ``|V|/sigma_V >= 5`` is at most the
    Gaussian two-sided expectation plus ``NULL_EXCESS_ALLOWED``. A positive-only null means the
    measurement is searching rather than forced, which stops the run.
    """
    v = np.asarray(v, float)
    s = np.asarray(sigma_v, float)
    ok = np.isfinite(v) & np.isfinite(s) & (s > 0)
    z = v[ok] / s[ok]
    n = int(z.size)
    if n == 0:
        return {"n": 0, "passed": False}
    lo, hi = stats.binom.interval(NULL_SIGN_CI, n, 0.5)
    frac_pos = float(np.mean(z > 0))
    std_z = float(np.std(z, ddof=1)) if n > 1 else float("nan")
    n5 = int(np.sum(np.abs(z) >= SNR_DETECT))
    expected5 = float(n * 2.0 * stats.norm.sf(SNR_DETECT))
    sign_ok = bool(lo / n <= frac_pos <= hi / n)
    std_ok = bool(NULL_STD_RANGE[0] <= std_z <= NULL_STD_RANGE[1])
    excess_ok = bool(n5 <= expected5 + NULL_EXCESS_ALLOWED)
    return {
        "n": n,
        "frac_positive": frac_pos,
        "frac_positive_ci": (float(lo / n), float(hi / n)),
        "std_z": std_z,
        "n_ge5": n5,
        "expected_ge5": expected5,
        "sign_ok": sign_ok,
        "std_ok": std_ok,
        "excess_ok": excess_ok,
        "passed": sign_ok and std_ok and excess_ok,
    }


def injection_check(
    level_sigma: np.ndarray,
    injected_mjy: np.ndarray,
    measured_v_mjy: np.ndarray,
    sigma_v: np.ndarray,
) -> dict:
    """C2: completeness and flux recovery per injection level.

    A recovery is ``|V|/sigma_V >= SNR_DETECT`` at the forced pixel. Pass (frozen): completeness
    >= ``INJ_COMPLETENESS_AT_7`` at 7 sigma, and the median recovered/injected flux at >= 10 sigma
    in ``INJ_FLUX_RATIO_RANGE``. The completeness curve is reported, not tuned.
    """
    lev = np.asarray(level_sigma, float)
    inj = np.asarray(injected_mjy, float)
    meas = np.asarray(measured_v_mjy, float)
    s = np.asarray(sigma_v, float)
    rec = np.abs(meas) / s >= SNR_DETECT
    curve = {}
    for L in np.unique(lev):
        m = lev == L
        curve[float(L)] = float(np.mean(rec[m]))
    bright = lev >= 10.0
    ratio = float(np.median(meas[bright] / inj[bright])) if bright.any() else float("nan")
    c7 = curve.get(7.0, float("nan"))
    comp_ok = bool(np.isfinite(c7) and c7 >= INJ_COMPLETENESS_AT_7)
    flux_ok = bool(INJ_FLUX_RATIO_RANGE[0] <= ratio <= INJ_FLUX_RATIO_RANGE[1])
    return {
        "completeness": curve,
        "median_ratio_ge10": ratio,
        "completeness_ok": comp_ok,
        "flux_ok": flux_ok,
        "passed": comp_ok and flux_ok,
    }


def recover_known_check(forced_v_mjy: np.ndarray, catalogue_v_mjy: np.ndarray) -> dict:
    """C3(a): forced low2 V against the RACS-low2 V-catalogue flux for the same stars.

    Signs must agree (a handedness flip is a failure for that star, counted, not dropped); the
    pass rule is the median of forced/catalogue over sign-agreeing stars in ``KNOWN_RATIO_RANGE``.
    """
    f = np.asarray(forced_v_mjy, float)
    c = np.asarray(catalogue_v_mjy, float)
    ok = np.isfinite(f) & np.isfinite(c) & (c != 0)
    agree = ok & (np.sign(f) == np.sign(c))
    ratio = float(np.median(f[agree] / c[agree])) if agree.any() else float("nan")
    passed = bool(KNOWN_RATIO_RANGE[0] <= ratio <= KNOWN_RATIO_RANGE[1])
    return {
        "n": int(ok.sum()),
        "n_sign_flip": int((ok & ~agree).sum()),
        "median_ratio": ratio,
        "passed": passed,
    }


def apply_detection_rule(
    v_mjy: np.ndarray,
    sigma_v: np.ndarray,
    i_mjy: np.ndarray,
    sigma_i: np.ndarray,
    leak_floor: np.ndarray | float,
    bright_neighbour: np.ndarray,
) -> dict[str, np.ndarray]:
    """The frozen candidate rule, per target-epoch.

    Candidate iff ``|V|/sigma_V >= SNR_DETECT``, AND ``|V|/I`` exceeds the local leakage floor
    whenever Stokes I is itself detected (``I/sigma_I >= SNR_DETECT``), AND no Stokes-I source
    brighter than ``BRIGHT_GUARD_MJY`` within ``BRIGHT_GUARD_ARCMIN`` (``bright_neighbour``).
    A candidate is not a detection: it still needs ``classify_emitter`` vetting and inspection.
    Returns the boolean ``candidate`` mask and the per-criterion masks behind it.
    """
    v = np.asarray(v_mjy, float)
    sv = np.asarray(sigma_v, float)
    i = np.asarray(i_mjy, float)
    si = np.asarray(sigma_i, float)
    floor = np.broadcast_to(np.asarray(leak_floor, float), v.shape)
    v_sig = np.abs(v) / sv >= SNR_DETECT
    i_det = i / si >= SNR_DETECT
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(i > 0, np.abs(v) / i, np.inf)
    leak_ok = ~i_det | (frac > floor)
    guard_ok = ~np.asarray(bright_neighbour, bool)
    return {
        "candidate": v_sig & leak_ok & guard_ok,
        "v_significant": v_sig,
        "clears_leakage": leak_ok,
        "clear_of_bright": guard_ok,
    }


def upper_limit(v_mjy: np.ndarray, sigma_v: np.ndarray, nsigma: float = LIMIT_NSIGMA) -> np.ndarray:
    """Conservative upper limit on |V|: ``|V_meas| + nsigma * sigma_V`` (the plan's definition)."""
    return np.abs(np.asarray(v_mjy, float)) + nsigma * np.asarray(sigma_v, float)


def luminosity_limit(s_mjy: np.ndarray, dist_pc: np.ndarray) -> np.ndarray:
    """Spectral luminosity ``4 pi d^2 S`` in erg s^-1 Hz^-1 (isotropic; no beaming correction)."""
    d = np.asarray(dist_pc, float) * PC_CM
    return 4.0 * np.pi * d**2 * np.asarray(s_mjy, float) * MJY_CGS


def two_epoch_compare(
    v1: np.ndarray,
    s1: np.ndarray,
    v2: np.ndarray,
    s2: np.ndarray,
    *,
    nu1_ghz: float = NU_LOW2_GHZ,
    nu2_ghz: float = NU_LOW3_GHZ,
    scale: float = SCALE_TERM,
    alpha_max: float = ALPHA_MAX,
) -> dict[str, np.ndarray]:
    """Two-epoch V comparison with a per-epoch flux-scale term and a spectral allowance.

    ``z_diff`` is the V difference over noise plus ``scale`` x each epoch's V, so two epochs
    whose fluxes differ by a calibration factor do not read as a change (the ``stokesv_discovery``
    V/I lesson). Because low2 and low3 differ in frequency, a change is only ``beyond_spectral`` if
    the epoch ratio also falls outside what any ``|alpha| <= alpha_max`` predicts, or V flips sign.
    """
    v1 = np.asarray(v1, float)
    v2 = np.asarray(v2, float)
    s1 = np.asarray(s1, float)
    s2 = np.asarray(s2, float)
    err = np.sqrt(s1**2 + s2**2 + (scale * v1) ** 2 + (scale * v2) ** 2)
    z = (v2 - v1) / err
    differs = np.abs(z) >= SNR_DETECT
    r = (nu2_ghz / nu1_ghz) ** alpha_max
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = v2 / v1
    flip = np.sign(v1) * np.sign(v2) < 0
    outside = ~((ratio >= 1.0 / r) & (ratio <= r))
    return {
        "z_diff": z,
        "differs": differs,
        "beyond_spectral": differs & (flip | outside),
    }


def scope_gate(
    sec_per_pair: float,
    fail_rate: float,
    dist_pc: np.ndarray,
    *,
    n_epochs: int = 2,
    budget_hours: float = BUDGET_HOURS,
) -> dict:
    """S0: full sample, or the largest distance cut whose projected wall-clock fits the budget.

    Projected hours = targets x epochs x ``sec_per_pair`` / (1 - ``fail_rate``) / 3600 (failed
    slots are retried, so each success costs 1/(1-f) attempts). The cut uses distances only --
    by construction no V or I value can influence the scope.
    """
    if not 0.0 <= fail_rate < 1.0:
        raise ValueError("fail_rate must be in [0, 1)")
    d = np.sort(np.asarray(dist_pc, float)[np.isfinite(dist_pc)])
    per_target_h = n_epochs * sec_per_pair / (1.0 - fail_rate) / 3600.0
    full_h = d.size * per_target_h
    if full_h <= budget_hours:
        return {"full_sample": True, "d_max_pc": float("inf"), "n": int(d.size), "hours": full_h}
    n_fit = int(np.floor(budget_hours / per_target_h))
    if n_fit == 0:
        return {"full_sample": False, "d_max_pc": 0.0, "n": 0, "hours": 0.0}
    # keep every target strictly closer than the first one that would not fit (no partial shells)
    d_cut = float(d[n_fit]) if n_fit < d.size else float("inf")
    n = int(np.sum(d < d_cut))
    return {"full_sample": False, "d_max_pc": d_cut, "n": n, "hours": n * per_target_h}


def synthetic_cutout(
    *,
    n_pix: int = 180,
    pix_arcsec: float = 2.5,
    rms_mjy: float = 0.25,
    beam: Beam | None = None,
    ra0: float = 150.0,
    dec0: float = -30.0,
    seed: int = SEED,
):
    """Offline fixture: a north-up TAN cutout of beam-correlated Gaussian noise, I and V.

    White noise is smoothed with the restoring beam and rescaled to ``rms_mjy`` so that adjacent
    pixels are correlated as in a real restored image. Returns ``(image_i, image_v, wcs, beam)``.
    """
    from astropy.wcs import WCS
    from scipy.ndimage import gaussian_filter

    beam = beam or Beam(15.0, 12.0, 30.0)  # offline fixture only; real beams come from headers
    rng = np.random.default_rng(seed)
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [n_pix / 2.0 + 0.5, n_pix / 2.0 + 0.5]
    w.wcs.cdelt = [-pix_arcsec / 3600.0, pix_arcsec / 3600.0]
    sig_pix = np.sqrt(beam.bmaj_arcsec * beam.bmin_arcsec) / 2.3548 / pix_arcsec

    def field() -> np.ndarray:
        f = gaussian_filter(rng.normal(size=(n_pix, n_pix)), sig_pix)
        return f / f.std() * rms_mjy

    return field(), field(), w, beam


def validate_synthetic(*, n_cutouts: int = 80, seed: int = SEED) -> dict:
    """Run C1 (null) and C2 (injection) end to end on synthetic cutouts -- the offline leg.

    Mirrors the real geometry: a 0.06-degree cutout centred on the target, ``NULL_PER_TARGET``
    forced nulls 2 arcmin out (80 cutouts x 5 = 400, the C1 minimum), and one injection per level
    at the target, each level in its own independent noise realisation. It validates the
    machinery -- forced photometry, beam injection, both pass rules -- not RACS.
    """
    rng = np.random.default_rng(seed)
    v_null, s_null = [], []
    lev, inj, meas, sig = [], [], [], []
    for k in range(n_cutouts):
        img_i, img_v, w, beam = synthetic_cutout(seed=seed * 1000 + k)
        n_pix = img_i.shape[0]
        ra_c, dec_c = (float(x) for x in w.pixel_to_world_values(n_pix / 2.0, n_pix / 2.0))
        ra_n, dec_n = null_positions(np.array([ra_c]), np.array([dec_c]), seed=seed + k)
        for rn, dn in zip(ra_n, dec_n, strict=True):
            m = forced_measure(img_i, img_v, w, float(rn), float(dn))
            v_null.append(m["v_peak"])
            s_null.append(m["v_rms"])
        for j, L in enumerate(INJ_LEVELS_SIGMA):
            # an independent noise realisation per level, so the curve's points are independent
            fi, fv, _, _ = synthetic_cutout(seed=seed * 1000 + k + 10_000 * (j + 1))
            sigma = float(np.std(fv))
            sign = rng.choice([-1.0, 1.0])
            flux = sign * L * sigma
            iv = inject_point_source(fv, w, ra_c, dec_c, flux, beam)
            ii = inject_point_source(fi, w, ra_c, dec_c, abs(flux) * 2.0, beam)
            mm = forced_measure(ii, iv, w, ra_c, dec_c)
            lev.append(L)
            inj.append(flux)
            meas.append(mm["v_peak"])
            sig.append(mm["v_rms"])
    return {
        "null": null_check(np.array(v_null), np.array(s_null)),
        "injection": injection_check(np.array(lev), np.array(inj), np.array(meas), np.array(sig)),
    }
