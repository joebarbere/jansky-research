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
from typing import Any

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
    "analyze_census",
    "assign_set",
    "cutout_leakage_floor",
    "poisson_upper",
    "snapshot_fraction_limit",
    "synthetic_census",
    "forced_measure",
    "in_footprint",
    "parse_polygon",
    "pick_tile",
    "run_s0",
    "run_c3",
    "run_c3b",
    "w0623_position",
    "run",
    "run_census",
    "backoff_seconds",
    "retry_policy",
    "census_leak_floors",
    "followup_lightcurve",
    "read_done",
    "census_targets",
    "bright_neighbour",
    "fetch_epoch_pair",
    "injection_check",
    "inject_point_source",
    "load_cutout",
    "controls_c1_c2",
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
    """Provenance split: ``"A"`` (calibrates) or ``"B"`` (certifies) from a hash of the Gaia id.

    The low bit of a Gaia DR2 ``source_id`` is not random -- all 14,915 Reyle+2018 ids are even --
    so the parity rule first frozen in the plan would have put every target in A (amended
    2026-10-15, pre-data). The first byte of SHA-256 of the decimal id is uniform and still
    data-blind: no flux, position or distance enters it.
    """
    import hashlib

    sid = np.atleast_1d(np.asarray(source_id, dtype=np.int64))
    bit = np.array([hashlib.sha256(str(int(x)).encode()).digest()[0] & 1 for x in sid], dtype=int)
    return np.where(bit == 0, "A", "B")


def parse_polygon(s_region: str) -> np.ndarray:
    """Vertices ``(n, 2)`` in degrees from an ObsCore ``POLYGON ICRS ra1 dec1 ra2 dec2 ...``."""
    tok = str(s_region).split()
    if not tok or tok[0].upper() != "POLYGON":
        raise ValueError(f"not a POLYGON s_region: {s_region[:40]!r}")
    nums = [float(x) for x in tok[1:] if x.upper() not in ("ICRS", "J2000", "FK5")]
    return np.asarray(nums, float).reshape(-1, 2)


def _unit(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def in_footprint(ra: np.ndarray, dec: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    """Points inside a convex spherical polygon (great-circle edges), either vertex winding.

    Inside iff on the polygon's side of every edge's great circle, where "the polygon's side" is
    fixed by the vertex centroid -- without that, the antipodal region also passes (all signs
    flipped together). Exact on the sphere, so it holds near the poles and across RA = 0.
    """
    p = _unit(np.atleast_1d(np.asarray(ra, float)), np.atleast_1d(np.asarray(dec, float)))
    v = _unit(vertices[:, 0], vertices[:, 1])
    normals = np.cross(v, np.roll(v, -1, axis=0))  # one per edge
    orient = np.sign(v.mean(axis=0) @ normals.T)  # the interior's sign for each edge
    return np.all((p @ normals.T) * orient >= 0, axis=1)


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


# ---- network (real leg; not unit-tested) -------------------------------------------------------
CASDA_TAP_URL = "https://casda.csiro.au/casda_vo_tools/tap"
EPOCH_LAMBDA_M = {"low2": 0.338, "low3": 0.318}  # ObsCore (em_min+em_max)/2, rounded to 1 mm


def _stokes_of(header) -> str:
    """``"i"`` or ``"v"`` from a FITS header's STOKES axis (CASDA cutout names do not carry it)."""
    for n in range(1, int(header.get("NAXIS", 0)) + 1):
        if str(header.get(f"CTYPE{n}", "")).strip().upper() == "STOKES":
            code = int(round(float(header[f"CRVAL{n}"])))
            return {1: "i", 4: "v"}.get(code, f"stokes{code}")
    raise RuntimeError("cutout header has no STOKES axis")


def pick_tile(ra: float, dec: float, regions: list[str]) -> int:
    """Index of the tile whose vertex centroid is nearest ``(ra, dec)`` among those containing it.

    Leakage rises towards tile edges, so the most central tile is preferred. Positions only: no
    pixel value can influence which observation a target is measured in. Returns -1 if none.
    """
    best, best_cos = -1, -2.0
    p = _unit(np.array([ra]), np.array([dec]))[0]
    for k, reg in enumerate(regions):
        v = parse_polygon(reg)
        if not in_footprint(np.array([ra]), np.array([dec]), v)[0]:
            continue
        c = _unit(v[:, 0], v[:, 1]).mean(axis=0)
        cosang = float(p @ (c / np.linalg.norm(c)))
        if cosang > best_cos:
            best, best_cos = k, cosang
    return best


JOB_TIMEOUT_S = 600.0  # a CASDA cutout job not done in 10 min is abandoned (aborted) and retried
JOB_RETRIES = 2  # extra in-slot attempts after a job ends ERROR / HTTP 5xx (transient on CASDA)
JOB_RETRY_WAIT_S = 20.0


def _cutout_urls(
    casda, rows, coord, radius_deg: float, poll: float, timeout: float = JOB_TIMEOUT_S
):  # pragma: no cover - network
    """astroquery's ``Casda.cutout`` with a deadline.

    ``Casda._run_job`` polls until the job leaves PENDING/QUEUED/EXECUTING and has no timeout; on
    2026-10-16 CASDA's SODA queue left jobs PENDING indefinitely and the census hung silently. This
    runs the same steps (datalink tokens -> SODA job -> POS -> RUN -> poll) but aborts the job and
    raises ``TimeoutError`` after ``timeout`` s, so the caller's backoff takes over.
    """
    import time
    from urllib.parse import unquote

    import astropy.units as u

    job = casda._create_job(rows, "cutout_service", False)
    spec = casda._args_to_payload(radius=radius_deg * u.deg, coordinates=coord)
    casda._add_cutout_params(job, False, spec)
    casda._request("POST", job + "/phase", data={"phase": "RUN"}, cache=False)
    t0 = time.time()
    while True:
        status = casda._read_job_status(casda._get_job_details_xml(job), False)
        if status not in ("EXECUTING", "QUEUED", "PENDING", "SUSPENDED"):
            break
        if time.time() - t0 > timeout:
            try:
                casda._request("POST", job + "/phase", data={"phase": "ABORT"}, cache=False)
            except Exception:
                pass
            raise TimeoutError(f"CASDA cutout job still {status} after {timeout:.0f} s: {job}")
        time.sleep(poll)
    if status != "COMPLETED":
        msg = ""
        try:
            summ = casda._get_job_details_xml(job).find("uws:errorSummary", casda._uws_ns)
            if summ is not None:
                msg = " ".join(t.strip() for t in summ.itertext() if t.strip())
        except Exception:
            pass
        raise ValueError(f"CASDA cutout job ended {status} ({msg or 'no errorSummary'}): {job}")
    details = casda._get_job_details_xml(job)
    results = details.find("uws:results", casda._uws_ns).findall("uws:result", casda._uws_ns)
    return [unquote(r.get("{http://www.w3.org/1999/xlink}href")) for r in results]


def fetch_epoch_pair(
    ra: float,
    dec: float,
    epoch: str,
    casda,
    *,
    radius_deg: float = CUTOUT_RADIUS_DEG,
    poll_interval: float = 20.0,
    obs_id: str | None = None,
    any_project: bool = False,
):  # pragma: no cover - network
    """One SODA job staging a same-observation I+V cutout pair for ``epoch`` at ``(ra, dec)``.

    ``any_project`` (with ``obs_id``; follow-up only, never the census): drop the RACS collection
    and band filters so a VAST/EMU/other ASKAP observation of the position can be fetched.

    Queries ObsCore for RACS taylor.0 restored conv I and V products of that epoch whose footprint
    contains the position, keeps observations that have both, picks the most central tile
    (:func:`pick_tile`), and stages I and V in a single cutout job (astroquery's default is one job
    per product). ``obs_id`` pins the observation instead (C3 compares against a catalogue measured
    in one specific SBID). Returns ``dict(i, v, wcs, beam, mjd, obs_id, files, t_tap, t_job, t_dl)`` with
    images in mJy/beam and timings in seconds, or ``None`` if the epoch has no I+V pair here.
    """
    import tempfile
    import time
    from pathlib import Path

    import astropy.units as u
    import pyvo
    import requests
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.wcs import WCS

    t0 = time.time()
    if any_project:
        if obs_id is None:
            raise ValueError("any_project needs obs_id")
        where = f"obs_id = '{obs_id}'"
    else:
        lam = EPOCH_LAMBDA_M[epoch]
        where = (
            "obs_collection = 'The Rapid ASKAP Continuum Survey' "
            f"AND em_min > {lam - 0.01:.3f} AND em_max < {lam + 0.01:.3f}"
        )
    q = f"""
    SELECT obs_id, filename, t_min, em_min, access_url, s_region FROM ivoa.obscore
    WHERE {where}
    AND (filename LIKE 'image.i.%' OR filename LIKE 'image.v.%')
    AND filename LIKE '%.taylor.0.restored.conv.fits'
    AND 1 = CONTAINS(POINT('ICRS', {ra:.8f}, {dec:.8f}), s_region)
    """
    tab = pyvo.dal.TAPService(CASDA_TAP_URL).search(q).to_table()
    t_tap = time.time() - t0
    groups: dict[str, dict[str, int]] = {}
    for k, row in enumerate(tab):
        groups.setdefault(str(row["obs_id"]), {})[str(row["filename"]).split(".")[1]] = k
    complete = [
        g for o, g in groups.items() if "i" in g and "v" in g and (obs_id is None or o == obs_id)
    ]
    if not complete:
        return None
    k = pick_tile(ra, dec, [str(tab["s_region"][g["i"]]) for g in complete])
    if k < 0:
        return None
    rows = tab[[complete[k]["i"], complete[k]["v"]]]
    t1 = time.time()
    for attempt in range(JOB_RETRIES + 1):
        try:
            urls = _cutout_urls(
                casda,
                rows["access_url",],
                SkyCoord(ra * u.deg, dec * u.deg),
                radius_deg,
                poll_interval,
            )
            break
        except (ValueError, requests.HTTPError) as exc:
            # A job that CASDA ends ERROR, or a 5xx, usually succeeds when simply resubmitted
            # (2026-10-22: 77 of the last 300 slot errors were "ended ERROR"). Stuck jobs
            # (TimeoutError) and auth errors are not retried here.
            transient = isinstance(exc, ValueError) or (
                exc.response is not None and exc.response.status_code >= 500
            )
            if not transient or attempt == JOB_RETRIES:
                raise
            time.sleep(JOB_RETRY_WAIT_S)
    t_job = time.time() - t1
    t2 = time.time()
    out: dict[str, Any] = {}
    for furl in (x for x in urls if x.endswith(".fits")):
        raw = requests.get(furl, timeout=200).content
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "c.fits"
            path.write_bytes(raw)
            with fits.open(path) as hd:
                hdr = hd[0].header
                data = np.squeeze(np.asarray(hd[0].data, float)) * 1000.0  # Jy -> mJy/beam
                stokes = _stokes_of(hdr)
                if stokes in out:
                    raise RuntimeError(f"two Stokes-{stokes} cutouts in one job: {urls}")
                out[stokes] = data
                out["wcs"] = WCS(hdr).celestial
                out["beam"] = Beam(
                    float(hdr["BMAJ"]) * 3600.0, float(hdr["BMIN"]) * 3600.0, float(hdr["BPA"])
                )
                if "MJD-OBS" in hdr:
                    out["mjd_hdr"] = float(hdr["MJD-OBS"])
                elif "DATE-OBS" in hdr:  # early VAST pilot products
                    from astropy.time import Time

                    out["mjd_hdr"] = float(Time(hdr["DATE-OBS"], scale="utc").mjd)
    t_dl = time.time() - t2
    if "i" not in out or "v" not in out:
        raise RuntimeError(f"cutout job returned {sorted(out)} for {urls}")
    if np.shape(out["i"]) != np.shape(out["v"]):
        raise RuntimeError("I/V cutout grids differ")
    mjd = float(rows["t_min"][0])
    if not np.isfinite(mjd):  # early pilot products carry no t_min in ObsCore
        mjd = float(out.pop("mjd_hdr", np.nan))
    out.update(
        mjd=mjd,
        freq_mhz=299.792458 / float(rows["em_min"][0]) if rows["em_min"][0] else np.nan,
        obs_id=str(rows["obs_id"][0]),
        files=[str(x) for x in rows["filename"]],
        t_tap=t_tap,
        t_job=t_job,
        t_dl=t_dl,
    )
    return out


def run_s0(
    *,
    username: str,
    pw_path: str = "~/.casda_pw",
    n_targets: int = 40,
    poll_interval: float = 3.0,
    table_path: str = "data/ucdv/reyle2018_ucds.ecsv",
    products_path: str = "data/ucdv/racs_v_products.ecsv",
    cutout_dir: str = "data/ucdv/cutouts",
    out_json: str = "results/ucdv_s0_throughput.json",
):  # pragma: no cover - network
    """S0: time single-job I+V pair fetches for ``n_targets`` random set-B targets, both epochs.

    Saves each pair (images, header WCS, beam, MJD) under ``cutout_dir`` for C1-C3, but computes
    and records **timing and failure counts only** -- no pixel statistic is evaluated here, per the
    plan. Positions are propagated to each epoch's MJD before the cutout is centred.
    """
    import json
    import time
    from pathlib import Path

    from astropy.table import Table

    from .stokesv import _casda_session
    from .stokesv_discovery import epoch_position

    ucds = Table.read(table_path)
    prods = Table.read(products_path)
    lam = np.round((np.asarray(prods["em_min"]) + np.asarray(prods["em_max"])) / 2, 3)
    ra0 = np.asarray(ucds["RAICRS"], float)
    de0 = np.asarray(ucds["DEICRS"], float)
    both = np.ones(len(ucds), bool)
    for ep in ("low2", "low3"):
        hit = np.zeros(len(ucds), bool)
        for reg in prods["s_region"][lam == EPOCH_LAMBDA_M[ep]]:
            hit |= in_footprint(ra0, de0, parse_polygon(str(reg)))
        both &= hit
    sets = assign_set(np.asarray(ucds["Source"]))
    pool = np.flatnonzero(both & (sets == "B"))
    pick = np.sort(np.random.default_rng(SEED).choice(pool, n_targets, replace=False))
    Path(cutout_dir).mkdir(parents=True, exist_ok=True)
    casda = _casda_session(username, pw_path)
    rows = []
    t_all = time.time()
    for idx in pick:
        r = ucds[int(idx)]
        for ep in ("low2", "low3"):
            rec: dict[str, Any] = {"source_id": int(r["Source"]), "epoch": ep}
            t = time.time()
            try:
                # centre the cutout on the position propagated to this epoch (a 2-pass fetch would
                # cost a second job; the tile's MJD spread within an epoch is < 1 yr, ~1" for UCDs)
                mjd_mid = 59750.0 if ep == "low2" else 60330.0
                ra, de = epoch_position(
                    np.array([r["RAICRS"]]),
                    np.array([r["DEICRS"]]),
                    np.array([r["pmRA"]]),
                    np.array([r["pmDE"]]),
                    GAIA_DR2_EPOCH,
                    mjd_mid,
                )
                out = fetch_epoch_pair(
                    float(ra[0]), float(de[0]), ep, casda, poll_interval=poll_interval
                )
                if out is None:
                    rec["status"] = "no_pair"
                else:
                    rec.update(
                        status="ok",
                        obs_id=out["obs_id"],
                        mjd=out["mjd"],
                        t_tap=out["t_tap"],
                        t_job=out["t_job"],
                        t_dl=out["t_dl"],
                    )
                    np.savez_compressed(
                        Path(cutout_dir) / f"{int(r['Source'])}_{ep}.npz",
                        i=out["i"],
                        v=out["v"],
                        wcs_header=out["wcs"].to_header_string(),
                        beam=np.array(
                            [out["beam"].bmaj_arcsec, out["beam"].bmin_arcsec, out["beam"].bpa_deg]
                        ),
                        mjd=out["mjd"],
                        obs_id=out["obs_id"],
                    )
            except Exception as exc:  # recorded, not raised: the failure rate is the measurement
                rec.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
                casda = _casda_session(username, pw_path)
            rec["t_total"] = time.time() - t
            rows.append(rec)
            print(rec["source_id"], ep, rec["status"], f"{rec['t_total']:.1f}s", flush=True)
    ok = [x for x in rows if x["status"] == "ok"]
    attempts = len(rows)
    summary = {
        "source": "real CASDA S0 throughput run (timings and failure counts only; no pixel stats)",
        "date": time.strftime("%Y-%m-%d"),
        "mode": f"single job per I+V pair, poll_interval={poll_interval}s",
        "n_targets": n_targets,
        "n_attempts": attempts,
        "n_ok": len(ok),
        "n_no_pair": sum(x["status"] == "no_pair" for x in rows),
        "n_error": sum(x["status"] == "error" for x in rows),
        "fail_rate": 1.0 - len(ok) / attempts,
        "sec_per_attempt_mean": float(np.mean([x["t_total"] for x in rows])),
        "sec_per_ok_median": float(np.median([x["t_total"] for x in ok])) if ok else None,
        "wall_clock_s": time.time() - t_all,
        "rows": rows,
    }
    Path(out_json).write_text(json.dumps(summary, indent=1))
    return summary


def load_cutout(path) -> dict[str, Any]:
    """Read a saved I+V cutout pair (``run_s0``/census ``.npz``) back into images, WCS and beam."""
    from astropy.io import fits
    from astropy.wcs import WCS

    z = np.load(path, allow_pickle=False)
    b = z["beam"]
    return {
        "i": z["i"],
        "v": z["v"],
        "wcs": WCS(fits.Header.fromstring(str(z["wcs_header"]))),
        "beam": Beam(float(b[0]), float(b[1]), float(b[2])),
        "mjd": float(z["mjd"]),
        "obs_id": str(z["obs_id"]),
    }


def controls_c1_c2(
    cutouts: list[dict[str, Any]],
    targets: list[tuple[float, float, float, float]],
    *,
    seed: int = SEED,
) -> dict:
    """C1 and C2 on real cutouts. ``targets[k]`` is ``(ra, dec, pmra, pmdec)`` at the Gaia epoch.

    For each cutout the target is propagated to that cutout's MJD; C1 reads forced I/V at the
    ``NULL_PER_TARGET`` null positions; C2 injects beam-shaped V sources (header beam) at
    ``INJ_LEVELS_SIGMA`` x the local annulus sigma_V, at the target pixel and at each null site,
    one injection per site per level in a copy of the image. Returns both rule outcomes plus the
    target-pixel-only completeness curve, reported separately.
    """
    from .stokesv_discovery import epoch_position

    rng = np.random.default_rng(seed)
    nv, ns, ni, nsi = [], [], [], []
    lev, inj, meas, sig, at_target = [], [], [], [], []
    for k, (c, (ra0, de0, pmra, pmde)) in enumerate(zip(cutouts, targets, strict=True)):
        ra, de = epoch_position(
            np.array([ra0]),
            np.array([de0]),
            np.array([pmra]),
            np.array([pmde]),
            GAIA_DR2_EPOCH,
            c["mjd"],
        )
        rt, dt = float(ra[0]), float(de[0])
        rn, dn = null_positions(np.array([rt]), np.array([dt]), seed=seed + k)
        sites = [(rt, dt, True)] + [
            (float(a), float(b), False) for a, b in zip(rn, dn, strict=True)
        ]
        for a, b, is_t in sites:
            m = forced_measure(c["i"], c["v"], c["wcs"], a, b)
            if not is_t:
                nv.append(m["v_peak"])
                ns.append(m["v_rms"])
                ni.append(m["i_peak"])
                nsi.append(m["i_rms"])
            sigma = m["v_rms"]
            if not np.isfinite(sigma):
                continue
            for L in INJ_LEVELS_SIGMA:
                flux = rng.choice([-1.0, 1.0]) * L * sigma
                iv = inject_point_source(c["v"], c["wcs"], a, b, flux, c["beam"])
                mm = forced_measure(c["i"], iv, c["wcs"], a, b)
                lev.append(L)
                inj.append(flux)
                meas.append(mm["v_peak"])  # the underlying sky is part of it, as for a real source
                sig.append(mm["v_rms"])
                at_target.append(is_t)
    lev_a, inj_a, meas_a, sig_a = map(np.asarray, (lev, inj, meas, sig))
    t = np.asarray(at_target, bool)
    return {
        "c1_null_v": null_check(np.asarray(nv), np.asarray(ns)),
        "c1_null_i": null_check(np.asarray(ni), np.asarray(nsi)),
        "c2_injection": injection_check(lev_a, inj_a, meas_a, sig_a),
        "c2_target_pixel_only": injection_check(lev_a[t], inj_a[t], meas_a[t], sig_a[t]),
        "n_injections_per_level": int(len(lev) / len(INJ_LEVELS_SIGMA)),
        "n_target_pixel_per_level": int(t.sum() / len(INJ_LEVELS_SIGMA)),
    }


def run_c3(
    *,
    username: str,
    pw_path: str = "~/.casda_pw",
    catalogue_path: str = "data/ucdv/racs_low2_v_catalogue.ecsv",
    cutout_dir: str = "data/ucdv/c3",
    out_json: str = "results/ucdv_controls_c3.json",
    poll_interval: float = 3.0,
):  # pragma: no cover - network
    """C3(a): forced low2 V at every star of the RACS-low2 V catalogue, in that star's own SBID.

    The catalogue (CASDA ``AS110.racs_low2_v_v01``, 67 rows with ``source_type == 'star'``) is
    measured in one specific observation per source, and stars flare, so the cutout is pinned to
    ``ASKAP-<sbid>`` rather than the most central tile. Positions are the catalogue's fitted
    ``ra``/``dec`` in that image (not ``dec_corr``), so the forced pixel is the catalogued one.
    """
    import json
    import time
    from pathlib import Path

    from astropy.table import Table

    from .stokesv import _casda_session

    cat = Table.read(catalogue_path)
    stars = cat[np.array([str(x) == "star" for x in cat["source_type"]])]
    Path(cutout_dir).mkdir(parents=True, exist_ok=True)
    casda = _casda_session(username, pw_path)
    rows = []
    for r in stars:
        rec: dict[str, Any] = {"name": str(r["name"]), "cat_v": float(r["peak_flux_v"])}
        try:
            out = fetch_epoch_pair(
                float(r["ra"]),
                float(r["dec"]),
                "low2",
                casda,
                poll_interval=poll_interval,
                obs_id=f"ASKAP-{int(r['sbid'])}",
            )
            if out is None:
                rec["status"] = "no_pair_in_sbid"
            else:
                m = forced_measure(out["i"], out["v"], out["wcs"], float(r["ra"]), float(r["dec"]))
                rec.update(status="ok", forced_v=m["v_peak"], forced_sigma_v=m["v_rms"])
                rec["cat_noise_v"] = float(r["noise_v"])
        except Exception as exc:
            rec.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
            casda = _casda_session(username, pw_path)
        rows.append(rec)
        print(rec["name"], rec["status"], flush=True)
    ok = [x for x in rows if x["status"] == "ok"]
    res = recover_known_check(
        np.array([x["forced_v"] for x in ok]), np.array([x["cat_v"] for x in ok])
    )
    summary = {
        "source": "real RACS-low2 cutouts pinned to catalogue SBIDs; CASDA AS110.racs_low2_v_v01",
        "date": time.strftime("%Y-%m-%d"),
        "n_stars": len(rows),
        "n_ok": len(ok),
        "c3a": res,
        "rows": rows,
    }
    Path(out_json).write_text(json.dumps(summary, indent=1, default=float))
    return summary


def bright_neighbour(image_i, wcs, ra: float, dec: float) -> bool:
    """True if any Stokes-I pixel above ``BRIGHT_GUARD_MJY`` lies within ``BRIGHT_GUARD_ARCMIN``."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.wcs.utils import proj_plane_pixel_scales

    img = np.asarray(image_i, float)
    px, py = wcs.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    scale = float(np.mean(proj_plane_pixel_scales(wcs)) * 3600.0)
    yy, xx = np.mgrid[0 : img.shape[0], 0 : img.shape[1]]
    near = np.hypot(xx - float(px), yy - float(py)) * scale <= BRIGHT_GUARD_ARCMIN * 60.0
    vals = img[near & np.isfinite(img)]
    return bool(vals.size and np.max(vals) > BRIGHT_GUARD_MJY)


def census_targets(
    table_path: str = "data/ucdv/reyle2018_ucds.ecsv",
    products_path: str = "data/ucdv/racs_v_products.ecsv",
    d_max_pc: float = 107.77952446412739,
):
    """The census list: UCDs in both low2 and low3 V footprints with ``1000/Plx < d_max_pc``.

    ``d_max_pc`` is the S0 outcome (results/ucdv_s0_throughput.json), not a tunable.
    """
    from astropy.table import Table

    ucds = Table.read(table_path)
    prods = Table.read(products_path)
    lam = np.round((np.asarray(prods["em_min"]) + np.asarray(prods["em_max"])) / 2, 3)
    ra0 = np.asarray(ucds["RAICRS"], float)
    de0 = np.asarray(ucds["DEICRS"], float)
    keep = 1000.0 / np.asarray(ucds["Plx"], float) < d_max_pc
    for ep in ("low2", "low3"):
        hit = np.zeros(len(ucds), bool)
        for reg in prods["s_region"][lam == EPOCH_LAMBDA_M[ep]]:
            hit |= in_footprint(ra0, de0, parse_polygon(str(reg)))
        keep &= hit
    return ucds[keep]


def backoff_seconds(consecutive_failures: int, *, base: float = 60.0, cap: float = 1800.0) -> float:
    """Exponential retry wait: ``base * 2**(k-1)`` seconds after ``k`` consecutive failures, capped."""
    k = max(int(consecutive_failures), 1)
    return float(min(cap, base * 2.0 ** (k - 1)))


def retry_policy(consecutive_failures: int, error: str = "") -> tuple[float, bool]:
    """(wait seconds, log in again?) after a failed slot.

    An isolated failure is usually one bad CASDA job: wait 30 s and keep the session. Only a run of
    3+ consecutive failures (an outage) or an authentication error escalates to the exponential
    backoff with a fresh login, which keeps the 2026-10-16 protection against login storms while
    no longer paying minutes of backoff, and a login, for every stray job error.
    """
    # "'CasdaClass' object has no attribute '_auth'": a failed login left the session without
    # credentials (2026-10-18/21/23); only a fresh login fixes it.
    auth = "401" in error or "403" in error or "authenticat" in error.lower() or "_auth" in error
    k = max(int(consecutive_failures), 1)
    if k < 3 and not auth:
        return 30.0, False
    return backoff_seconds(max(k - 2, 1)), True


def read_done(path) -> set[tuple[int, str]]:
    """(source_id, epoch) slots already in the census JSONL. A truncated last line (the process
    killed mid-write) is ignored rather than crashing the resume; any other bad line raises."""
    import json
    from pathlib import Path

    done: set[tuple[int, str]] = set()
    p = Path(path)
    if not p.exists():
        return done
    lines = p.read_text().splitlines()
    for k, line in enumerate(lines):
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            if k == len(lines) - 1:
                continue
            raise
        done.add((int(r["source_id"]), str(r["epoch"])))
    return done


CENSUS_PASSES = 4  # a pass retries every slot that errored in the previous one


def run_census(
    *,
    username: str,
    pw_path: str = "~/.casda_pw",
    out_jsonl: str = "data/ucdv/census.jsonl",
    cutout_dir: str = "data/ucdv/census_cutouts",
    poll_interval: float = 3.0,
    limit: int | None = None,
):  # pragma: no cover - network, ~120 h
    """The census: forced I/V at every target, both epochs; resumable, one JSON line per slot.

    For each target-epoch: fetch the I+V pair (most central tile), propagate the Gaia DR2 position
    to that cutout's MJD, measure forced I/V there, the 5 nulls (ongoing C1), and the bright-source
    guard; save the cutout (float32) for vetting. Slots already in ``out_jsonl`` are skipped, so
    the run can be stopped and restarted. Detection, limits and comparison run afterwards on the
    JSONL, never inside this loop.
    """
    import time
    from pathlib import Path

    from .stokesv import _casda_session

    tg = census_targets()
    out = Path(out_jsonl)
    if out.exists() and out.stat().st_size and not out.read_bytes().endswith(b"\n"):
        with out.open("ab") as fh0:  # killed mid-write: end the partial line so appends stay valid
            fh0.write(b"\n")
    Path(cutout_dir).mkdir(parents=True, exist_ok=True)
    casda = _casda_session(username, pw_path)
    n = 0
    for pass_no in range(CENSUS_PASSES):
        done = read_done(out)
        errors_before = n
        attempted, failed = _census_pass(
            tg, done, out, cutout_dir, casda, username, pw_path, poll_interval, limit, n
        )
        n = attempted
        print(
            time.strftime("%Y-%m-%d %H:%M:%S"),
            f"census pass {pass_no + 1}: {attempted - errors_before} slots attempted, "
            f"{failed} failed",
            flush=True,
        )
        if failed == 0 or (limit is not None and n >= limit):
            break
    return n


def _census_pass(
    tg, done, out, cutout_dir, casda, username, pw_path, poll_interval, limit, n
):  # pragma: no cover - network
    """One pass over every slot not in ``done``. Returns (slots attempted so far, failures)."""
    import json
    import time
    from pathlib import Path

    from .stokesv import _casda_session
    from .stokesv_discovery import epoch_position

    fails = 0  # consecutive failed slots, drives the backoff
    failed = 0
    with out.open("a") as fh:
        for r in tg:
            sid = int(r["Source"])
            for ep in ("low2", "low3"):
                if (sid, ep) in done:
                    continue
                if limit is not None and n >= limit:
                    return n, failed
                n += 1
                rec: dict[str, Any] = {"source_id": sid, "epoch": ep, "t_start": time.time()}
                pos = (np.array([r["RAICRS"]]), np.array([r["DEICRS"]]))
                pm = (np.array([r["pmRA"]]), np.array([r["pmDE"]]))
                try:
                    mid = 59750.0 if ep == "low2" else 60330.0
                    ra_c, de_c = epoch_position(*pos, *pm, GAIA_DR2_EPOCH, mid)
                    pair = fetch_epoch_pair(
                        float(ra_c[0]), float(de_c[0]), ep, casda, poll_interval=poll_interval
                    )
                    if pair is None:
                        rec["status"] = "no_pair"
                    else:
                        ra_t, de_t = epoch_position(*pos, *pm, GAIA_DR2_EPOCH, pair["mjd"])
                        rt, dt = float(ra_t[0]), float(de_t[0])
                        m = forced_measure(pair["i"], pair["v"], pair["wcs"], rt, dt)
                        rn, dn = null_positions(
                            np.array([rt]), np.array([dt]), seed=SEED + sid % 997
                        )
                        nulls = [
                            forced_measure(pair["i"], pair["v"], pair["wcs"], float(a), float(b))
                            for a, b in zip(rn, dn, strict=True)
                        ]
                        rec.update(
                            status="ok",
                            obs_id=pair["obs_id"],
                            mjd=pair["mjd"],
                            ra=rt,
                            dec=dt,
                            v=m["v_peak"],
                            sigma_v=m["v_rms"],
                            i=m["i_peak"],
                            sigma_i=m["i_rms"],
                            bright=bright_neighbour(pair["i"], pair["wcs"], rt, dt),
                            null_v=[x["v_peak"] for x in nulls],
                            null_sigma_v=[x["v_rms"] for x in nulls],
                            beam=[pair["beam"].bmaj_arcsec, pair["beam"].bmin_arcsec],
                        )
                        np.savez_compressed(
                            Path(cutout_dir) / f"{sid}_{ep}.npz",
                            i=pair["i"].astype(np.float32),
                            v=pair["v"].astype(np.float32),
                            wcs_header=pair["wcs"].to_header_string(),
                            beam=np.array(
                                [
                                    pair["beam"].bmaj_arcsec,
                                    pair["beam"].bmin_arcsec,
                                    pair["beam"].bpa_deg,
                                ]
                            ),
                            mjd=pair["mjd"],
                            obs_id=pair["obs_id"],
                        )
                except Exception as exc:  # retried in the next pass
                    rec.update(status="error", error=f"{type(exc).__name__}: {exc}"[:500])
                    fails += 1
                    failed += 1
                    # retry_policy: isolated job errors wait 30 s on the same session; a run of
                    # failures (an outage) backs off exponentially with a fresh login, so an
                    # outage still costs tens of logins, not hundreds (2026-10-16).
                    wait, relogin = retry_policy(fails, rec["error"])
                    print(
                        time.strftime("%Y-%m-%d %H:%M:%S"),
                        f"slot error #{fails} {sid}/{ep} ({rec['error']}); sleeping {wait:.0f} s",
                        flush=True,
                    )
                    time.sleep(wait)
                    if relogin:
                        try:
                            casda = _casda_session(username, pw_path)
                        except Exception as login_exc:
                            print(time.strftime("%Y-%m-%d %H:%M:%S"), "login failed:", login_exc)
                else:
                    fails = 0
                rec["t_total"] = time.time() - rec["t_start"]
                if rec["status"] != "error":
                    fh.write(json.dumps(rec, default=float) + "\n")
                    fh.flush()
    return n, failed


# C3(b) target, from Rose et al. (arXiv:2306.15219, full text): RACS-mid position at MJD 59216
# and proper motion (Sec. II.4), distance 11.44 +/- 0.38 pc. A T8 dwarf -- not in Gaia, so not in
# the Reyle sample. Reported, not pass/fail (plan 59).
W0623_RA_DEG = (6 + 23 / 60 + 9.28 / 3600) * 15.0
W0623_DEC_DEG = -(4 + 56 / 60 + 22.8 / 3600)
W0623_REF_MJD = 59216.0
W0623_PM_MAS_YR = (-930.0, 170.0)  # (mu_alpha*, mu_delta)


def followup_lightcurve(
    ra: float, dec: float, pmra: float, pmdec: float, *, username: str, pw_path: str = "~/.casda_pw"
):  # pragma: no cover - network
    """Post-hoc follow-up of a census candidate: forced I/V at every ASKAP observation of the
    position that has I and V images (any project), at the position propagated from Gaia DR2
    (``ra``, ``dec`` at epoch 2015.5, ``pmra``/``pmdec`` in mas/yr) to each observation's MJD.

    Not part of the frozen census: everything it returns is post hoc and is labelled so.
    """
    import pyvo

    from .stokesv import _casda_session
    from .stokesv_discovery import epoch_position

    q = f"""
    SELECT DISTINCT obs_id FROM ivoa.obscore
    WHERE filename LIKE 'image.v.%' AND filename LIKE '%.taylor.0.restored.conv.fits'
    AND 1 = CONTAINS(POINT('ICRS', {ra:.8f}, {dec:.8f}), s_region)
    """
    obs = sorted(str(x) for x in pyvo.dal.TAPService(CASDA_TAP_URL).search(q).to_table()["obs_id"])
    casda = _casda_session(username, pw_path)
    out = []
    for oid in obs:
        rec: dict[str, Any] = {"obs_id": oid}
        try:
            pair = None
            for attempt in range(3):
                try:
                    pair = fetch_epoch_pair(ra, dec, "low2", casda, obs_id=oid, any_project=True)
                    break
                except Exception:
                    if attempt == 2:
                        raise
            if pair is None:
                rec["status"] = "no_pair"
            else:
                pos = (np.array([ra]), np.array([dec]))
                pm = (np.array([pmra]), np.array([pmdec]))
                rt, dt = (
                    float(x[0]) for x in epoch_position(*pos, *pm, GAIA_DR2_EPOCH, pair["mjd"])
                )
                m = forced_measure(pair["i"], pair["v"], pair["wcs"], rt, dt)
                rec.update(
                    status="ok",
                    mjd=pair["mjd"],
                    freq_mhz=pair["freq_mhz"],
                    v=m["v_peak"],
                    sigma_v=m["v_rms"],
                    i=m["i_peak"],
                    sigma_i=m["i_rms"],
                    bright=bright_neighbour(pair["i"], pair["wcs"], rt, dt),
                    files=pair["files"],
                )
        except Exception as exc:
            rec.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
        print(rec, flush=True)
        out.append(rec)
    return out


def w0623_position(mjd: float) -> tuple[float, float]:
    """WISE J0623-0456 propagated from its RACS-mid position to ``mjd`` with its proper motion."""
    from .stokesv_discovery import epoch_position

    ref_yr = (W0623_REF_MJD - 51544.5) / 365.25 + 2000.0
    ra, dec = epoch_position(
        np.array([W0623_RA_DEG]),
        np.array([W0623_DEC_DEG]),
        np.array([W0623_PM_MAS_YR[0]]),
        np.array([W0623_PM_MAS_YR[1]]),
        ref_yr,
        mjd,
    )
    return float(ra[0]), float(dec[0])


def run_c3b(*, username: str, pw_path: str = "~/.casda_pw"):  # pragma: no cover - network
    """C3(b): forced V for WISE J0623-0456 in low2 and low3 (reported, not a gate)."""
    from .stokesv import _casda_session

    casda = _casda_session(username, pw_path)
    out: dict[str, dict[str, Any]] = {}
    for ep, mid in (("low2", 59750.0), ("low3", 60330.0)):
        ra, dec = w0623_position(mid)
        pair = fetch_epoch_pair(ra, dec, ep, casda, poll_interval=3.0)
        if pair is None:
            out[ep] = {"status": "no_pair"}
            continue
        ra, dec = w0623_position(pair["mjd"])
        m = forced_measure(pair["i"], pair["v"], pair["wcs"], ra, dec)
        out[ep] = {"status": "ok", "mjd": pair["mjd"], "obs_id": pair["obs_id"], **m}
    return out


L_GRID = np.logspace(13.0, 17.0, 41)  # erg s^-1 Hz^-1, for the snapshot-fraction curve
L_QUOTE = (3e14, 1e15, 1e16)  # luminosities the paper quotes a fraction limit at
SENS_NSIGMA = (5.0, 7.0)  # detectability for the fraction limit: the rule's 5 sigma, and C2's 7


def _sensitivity_and_detections(res: dict[str, Any], nsigma: float):
    """(L each sane slot could have detected at ``nsigma`` sigma_V, L of held/candidate slots).

    GATE-2 (2026-10-26): the denominator must not use the measured V (``|V|+3 sigma`` drops
    positive fluctuations and the detection's own slot), and a held detection that vetting
    confirms belongs in the numerator at its measured L, not at its upper limit.
    """
    s = res["sane"]
    sens = luminosity_limit(nsigma * res["sigma_v"][s], res["dist_pc"][s])
    det = res["candidate"] | res["held"]
    l_det = luminosity_limit(np.abs(res["v"][det]), res["dist_pc"][det])
    return sens, l_det


def _ltag(L: float) -> str:
    e = int(np.floor(np.log10(L)))
    return f"{round(L / 10**e):d}e{e}"


def _metrics(res: dict[str, Any], prefix: str) -> dict[str, Any]:
    """Flatten an ``analyze_census`` result into the namespaced metrics the paper cites."""
    m: dict[str, Any] = {
        f"{prefix}_n_targets": res["n_targets"],
        f"{prefix}_n_slots": res["n_slots"],
        f"{prefix}_n_sane": res["n_sane"],
        f"{prefix}_n_candidates": res["n_candidates"],
        f"{prefix}_n_held": res["n_held"],
        f"{prefix}_n_pairs": res["n_pairs"],
        f"{prefix}_n_beyond_spectral": int(np.sum(res["two_epoch"]["beyond_spectral"])),
        f"{prefix}_median_sigma_v_mjy": round(float(np.nanmedian(res["sigma_v"])), 3),
        f"{prefix}_median_limit_mjy": round(float(np.nanmedian(res["upper_limit_mjy"])), 3),
        f"{prefix}_null_frac_positive": round(
            float(res["null_monitor"].get("frac_positive", np.nan)), 3
        ),
        f"{prefix}_null_std": round(float(res["null_monitor"].get("std_z", np.nan)), 3),
        f"{prefix}_null_passed": bool(res["null_monitor"].get("passed", False)),
    }
    for nsig in SENS_NSIGMA:
        sens, l_det = _sensitivity_and_detections(res, nsig)
        frac = snapshot_fraction_limit(sens, np.array(L_QUOTE), l_det)
        k = "" if nsig == SENS_NSIGMA[0] else f"{int(nsig)}"
        for L, n, nd, f in zip(
            L_QUOTE, frac["n_eff"], frac["n_det"], frac["frac_upper"], strict=True
        ):
            tag = _ltag(L)
            m[f"{prefix}_neff{k}_{tag}"] = int(n)
            m[f"{prefix}_ndet{k}_{tag}"] = int(nd)
            ok = np.isfinite(f) and f <= 1.0  # a "fraction limit" above 1 constrains nothing
            m[f"{prefix}_frac{k}_upper_{tag}"] = round(float(f), 5) if ok else None
    _, l_det = _sensitivity_and_detections(res, SENS_NSIGMA[0])
    m[f"{prefix}_n_detections"] = int(l_det.size)
    m[f"{prefix}_det_l"] = [float(f"{x:.3g}") for x in l_det]
    return m


def _control_metrics(results_dir) -> dict[str, Any]:
    """S0/C1/C2/C3 numbers from their committed real results files, as ``real_*`` metrics."""
    import json
    from pathlib import Path

    rd = Path(results_dir)
    m: dict[str, Any] = {}
    s0p, cp, c3p = (
        rd / f"ucdv_{n}.json" for n in ("s0_throughput", "controls_c1c2", "controls_c3")
    )
    if s0p.exists():
        s0 = json.loads(s0p.read_text())
        m.update(
            real_s0_n_ok=s0["n_ok"],
            real_s0_n_attempts=s0["n_attempts"],
            real_s0_sec_per_pair=round(float(s0["sec_per_attempt_mean"]), 1),
        )
    if cp.exists():
        c = json.loads(cp.read_text())
        nv, inj = c["c1_null_v"], c["c2_injection"]
        m.update(
            real_c1_n=nv["n"],
            real_c1_frac_pos=round(float(nv["frac_positive"]), 3),
            real_c1_std=round(float(nv["std_z"]), 2),
            real_c2_n_per_level=c["n_injections_per_level"],
            real_c2_comp7=round(float(inj["completeness"]["7.0"]), 3),
            real_c2_ratio=round(float(inj["median_ratio_ge10"]), 3),
        )
    scp = rd / "ucdv_scope.json"
    if scp.exists():
        sc = json.loads(scp.read_text())
        m.update(
            real_n_reyle=sc["n_reyle"],
            real_n_both_epochs=sc["n_both_epochs"],
            real_d_max_pc=sc["d_max_pc"],
            real_n_census=sc["n_census"],
            real_full_hours=sc["projected_full_hours"],
        )
    if c3p.exists():
        c3 = json.loads(c3p.read_text())
        m.update(
            real_c3_n=c3["c3a"]["n"],
            real_c3_ratio=round(float(c3["c3a"]["median_ratio"]), 3),
            real_c3_flips=c3["c3a"]["n_sign_flip"],
        )
    return m


def _tex_sci(x: float, digits: int = 2) -> str:
    """``4.23e14`` -> ``4.2\\times10^{14}`` for a macro (math mode in the paper)."""
    if not np.isfinite(x) or x == 0:
        return "--"
    e = int(np.floor(np.log10(abs(x))))
    return f"{x / 10**e:.{digits - 1}f}\\times10^{{{e}}}"


def _detection_metrics(results_dir, real: dict[str, Any], crow: list[dict[str, Any]]) -> dict:
    """Numbers for the note's results: slot accounting, each held/candidate slot, and the
    post-hoc follow-up light curve (``results/ucdv_followup_j1123.json``) of the detection."""
    import json
    from pathlib import Path
    from statistics import NormalDist

    from astropy.time import Time

    m: dict[str, Any] = {
        "real_n_missing": int(2 * real["n_targets"] - real["n_slots"]),
        "real_n_insane": int(real["n_slots"] - real["n_sane"]),
    }
    det = real["candidate"] | real["held"]
    if int(det.sum()) != 1:  # the note is written for exactly one vetted detection
        return m
    k = int(np.flatnonzero(det)[0])
    sid, ep = int(real["source_id"][k]), str(real["epoch"][k])
    row = next(r for r in crow if int(r["source_id"]) == sid and r["epoch"] == ep)
    other = [r for r in crow if int(r["source_id"]) == sid and r["epoch"] != ep]
    v, sv, i, si = (float(row[x]) for x in ("v", "sigma_v", "i", "sigma_i"))
    frac, frac_err = v / i, abs(v / i) * float(np.hypot(sv / v, si / i))
    d = float(real["dist_pc"][k])
    m.update(
        real_det_gaia=str(sid),
        real_det_epoch=ep,
        real_det_obs=str(row["obs_id"]).replace("ASKAP-", ""),
        real_det_date=Time(float(row["mjd"]), format="mjd").datetime.strftime("%Y %B %-d"),
        real_det_v=round(v, 2),
        real_det_sigma_v=round(sv, 2),
        real_det_snr=round(v / sv, 1),
        real_det_i=round(i, 2),
        real_det_isnr=round(i / si, 1),
        real_det_frac=round(frac, 1),
        real_det_frac_err=round(frac_err, 1),
        real_det_dist_pc=round(d, 1),
        real_det_l_tex=_tex_sci(float(luminosity_limit(np.array([abs(v)]), np.array([d]))[0])),
    )
    if other:
        o = other[0]
        m["real_det_other_epoch"] = str(o["epoch"])
        m["real_det_other_snr"] = round(float(o["v"]) / float(o["sigma_v"]), 1)
        lo, hi = (o, row) if ep == "low3" else (row, o)
        cmp = two_epoch_compare(
            np.array([float(lo["v"])]),
            np.array([float(lo["sigma_v"])]),
            np.array([float(hi["v"])]),
            np.array([float(hi["sigma_v"])]),
        )
        m["real_det_pair_z"] = round(float(cmp["z_diff"][0]), 1)
    # the frozen two-epoch rule's smallest detectable change at the median sigma_V
    m["real_two_epoch_floor_mjy"] = round(
        SNR_DETECT * float(np.sqrt(2.0)) * float(np.nanmedian(real["sigma_v"][real["sane"]])), 2
    )
    fp = Path(results_dir) / "ucdv_followup_j1123.json"
    if not fp.exists():
        return m
    fu = json.loads(fp.read_text())
    eps = fu["epochs"]
    meas = [
        r for r in eps if r["status"] == "ok" and np.isfinite(r["v"]) and np.isfinite(r["sigma_v"])
    ]
    indep = [r for r in meas if str(r["obs_id"]) != str(row["obs_id"])]
    z = np.array([r["v"] / r["sigma_v"] for r in indep])
    best = indep[int(np.argmax(z))]
    p1 = 1.0 - NormalDist().cdf(float(z.max()))
    sub_ = z[np.abs(z) < 2.5]
    # the deepest same-band observation without a detection: the matched non-detection
    best_freq = float(best.get("freq_mhz", np.nan))
    same = [
        r
        for r in indep
        if r is not best
        and abs(float(r.get("freq_mhz", np.nan)) - best_freq) < 1.0
        and r["v"] / r["sigma_v"] < 2.5
    ]
    if same:
        mt = min(same, key=lambda r: r["sigma_v"])
        m.update(
            real_fu_match_obs=str(mt["obs_id"]).replace("ASKAP-", ""),
            real_fu_match_year=int(Time(float(mt["mjd"]), format="mjd").datetime.year),
            real_fu_match_v=f"{float(mt['v']):.3f}",
            real_fu_match_sigma_v=f"{float(mt['sigma_v']):.3f}",
            real_fu_match_z=round(
                (float(best["v"]) - float(mt["v"]))
                / float(np.hypot(best["sigma_v"], mt["sigma_v"])),
                1,
            ),
        )
    # epochs >= 2.5 sigma in V whose V exceeds I by > 2 sigma: not physical on their own
    viol = [
        r
        for r in indep
        if r["v"] / r["sigma_v"] >= 2.5
        and (r["v"] - r["i"]) / float(np.hypot(r["sigma_v"], r["sigma_i"])) > 2.0
    ]
    m["real_fu_n_v_gt_i"] = len(viol)
    if viol:
        m["real_fu_v_gt_i_obs"] = ", ".join(str(r["obs_id"]).replace("ASKAP-", "") for r in viol)
    m.update(
        real_fu_n_obs=len(eps),
        real_fu_n_nan=sum(
            r["status"] == "ok" and not (np.isfinite(r["v"]) and np.isfinite(r["sigma_v"]))
            for r in eps
        ),
        real_fu_n_refused=sum(r["status"] == "error" for r in eps),
        real_fu_n_meas=len(meas),
        real_fu_n_indep=len(indep),
        real_fu_first_year=int(Time(min(r["mjd"] for r in meas), format="mjd").datetime.year),
        real_fu_last_year=int(Time(max(r["mjd"] for r in meas), format="mjd").datetime.year),
        real_fu_n_pos=int(np.sum(z >= 2.5)),
        real_fu_n_neg=int(np.sum(z <= -2.5)),
        real_fu_best_obs=str(best["obs_id"]).replace("ASKAP-", ""),
        real_fu_best_date=Time(float(best["mjd"]), format="mjd").datetime.strftime("%Y %B"),
        real_fu_best_v=f"{float(best['v']):.3f}",
        real_fu_best_sigma_v=f"{float(best['sigma_v']):.3f}",
        real_fu_best_snr=round(float(z.max()), 1),
        real_fu_best_isnr=round(float(best["i"]) / float(best["sigma_i"]), 1),
        real_fu_best_p_tex=_tex_sci(p1, 1),
        real_fu_best_p_trials_tex=_tex_sci(min(1.0, p1 * len(indep)), 1),
        real_fu_sub_n=int(sub_.size),
        real_fu_sub_mean=round(float(sub_.mean()), 2),
    )
    return m


def beam_centroid(stamp: np.ndarray, pix_arcsec: float, beam, star_offset_pix) -> dict:
    """Peak and centroid of a source near the stamp centre, fitted with the restoring beam's
    shape (amplitude and position free) inside +-15 arcsec; offsets are from the star's
    propagated position, in arcsec (+x = pixel x). Returns NaNs if the fit fails."""
    from scipy.optimize import curve_fit

    im = np.asarray(stamp, float)
    n = im.shape[0]
    c = (n - 1) / 2.0
    bmaj, bmin, bpa = beam
    sx, sy, th = bmin / 2.3548 / pix_arcsec, bmaj / 2.3548 / pix_arcsec, np.radians(bpa)
    yy, xx = np.mgrid[0:n, 0:n]
    r = int(round(15.0 / pix_arcsec))
    m = (np.abs(xx - c) <= r) & (np.abs(yy - c) <= r) & np.isfinite(im)

    def model(xy, a, x0, y0):
        x, y = xy
        ct, st = np.cos(th), np.sin(th)
        xr, yr = ct * (x - x0) + st * (y - y0), -st * (x - x0) + ct * (y - y0)
        return a * np.exp(-0.5 * ((xr / sx) ** 2 + (yr / sy) ** 2))

    try:
        popt, _ = curve_fit(model, (xx[m], yy[m]), im[m], p0=(float(im[m].max()), c, c))
    except Exception:
        return {"peak": np.nan, "dx": np.nan, "dy": np.nan, "r": np.nan}
    dx = (popt[1] - (c + star_offset_pix[0])) * pix_arcsec
    dy = (popt[2] - (c + star_offset_pix[1])) * pix_arcsec
    return {"peak": float(popt[0]), "dx": float(dx), "dy": float(dy), "r": float(np.hypot(dx, dy))}


def blank_sky_tail(image_i, image_v, wcs, ra: float, dec: float, *, step_arcsec: float = 16.0):
    """Forced V/sigma_V on a grid of blank positions >= ``step_arcsec`` apart across a cutout,
    skipping 30 arcsec around (ra, dec), the 45 arcsec border and any position with I > 5 sigma_I.
    The cutout's own empirical noise tail, for a detection made in it."""
    pix = float(abs(wcs.proj_plane_pixel_scales()[0].to_value("arcsec")))
    ny, nx = np.shape(image_v)
    step, edge = int(np.ceil(step_arcsec / pix)), int(np.ceil(45.0 / pix))
    cx, cy = (float(q) for q in wcs.world_to_pixel_values(ra, dec))
    z = []
    for y in range(edge, ny - edge, step):
        for x in range(edge, nx - edge, step):
            if np.hypot(x - cx, y - cy) * pix < 30.0:
                continue
            r, d = wcs.pixel_to_world_values(x, y)
            m = forced_measure(image_i, image_v, wcs, float(r), float(d))
            if not (np.isfinite(m["v_rms"]) and m["v_rms"] > 0) or m["i_peak"] > 5 * m["i_rms"]:
                continue
            z.append(m["v_peak"] / m["v_rms"])
    za = np.asarray(z)
    return {
        "n": int(za.size),
        "frac_positive": float(np.mean(za > 0)) if za.size else np.nan,
        "std": float(za.std()) if za.size else np.nan,
        "max_abs": float(np.abs(za).max()) if za.size else np.nan,
        "n_ge3": int(np.sum(np.abs(za) >= 3.0)),
    }


def completeness_weighted_neff(
    res: dict[str, Any], comp: dict[str, float], l_quote: float
) -> float:
    """Effective number of sane snapshots that could detect ``l_quote``, each weighted by the C2
    injection completeness at its own L/(L per sigma_V) (linear interpolation in sigma, 0 below
    the lowest level, 1 above the highest). A cross-check on the 5 sigma step definition."""
    s = res["sane"]
    per_sigma = luminosity_limit(res["sigma_v"][s], res["dist_pc"][s])
    snr = l_quote / per_sigma
    lv = np.array(sorted(float(k) for k in comp))
    cv = np.array([float(comp[k]) for k in sorted(comp, key=float)])
    w = np.interp(snr, lv, cv, left=0.0, right=float(cv[-1]))
    return float(np.nansum(w))


def _tail_and_astrometry_metrics(results_dir, real, crow) -> dict[str, Any]:
    """Referee round 1 (2026-10-26): heavy-tail check of the census nulls, the expected number of
    chance >= 5 sigma slots, a completeness-weighted fraction limit, and (from the committed
    stamps file) beam-fitted centroids and each detection image's own blank-sky tail."""
    import json
    from pathlib import Path
    from statistics import NormalDist

    m: dict[str, Any] = {}
    nz = np.concatenate(
        [
            np.asarray(r.get("null_v", []), float) / np.asarray(r.get("null_sigma_v", []), float)
            for r in crow
            if r.get("status") == "ok"
        ]
    )
    nz = nz[np.isfinite(nz)]
    nd = NormalDist()
    for k in (3, 4, 5):
        obs = int(np.sum(np.abs(nz) >= k))
        exp = 2.0 * (1.0 - nd.cdf(k)) * nz.size
        m[f"real_null_n_ge{k}"] = obs
        m[f"real_null_exp_ge{k}"] = round(exp, 1 if exp >= 1 else 2)
        m[f"real_null_ratio_ge{k}"] = round(obs / exp, 1) if exp > 0 else None
    m["real_chance_ge5"] = round(m["real_null_n_ge5"] / nz.size * real["n_slots"], 1)
    c1p = Path(results_dir) / "ucdv_controls_c1c2.json"
    if c1p.exists():
        comp = json.loads(c1p.read_text())["c2_injection"]["completeness"]
        neff = completeness_weighted_neff(real, comp, 1e15)
        m["real_neff_cw_1e15"] = int(round(neff))
        m["real_frac_cw_upper_1e15"] = round(poisson_upper(0) / neff, 4) if neff > 0 else None
    sp = Path(results_dir) / "ucdv_followup_j1123_stamps.json"
    if sp.exists():
        stamps = json.loads(sp.read_text())
        tags = {"ASKAP-56518": "racs", "ASKAP-82408": "emu"}
        for oid, st in stamps["stamps"].items():
            t = tags.get(oid, oid)
            for stokes in ("v", "i"):
                cfit = beam_centroid(
                    st[f"{stokes}_over_sigma"],
                    st["pix_arcsec"],
                    st["beam_arcsec_deg"],
                    st["star_offset_pix"],
                )
                m[f"real_cen_{stokes}_{t}"] = round(cfit["r"], 1)
            if "blank_tail" in st:
                bt = st["blank_tail"]
                m[f"real_tail_n_{t}"] = bt["n"]
                m[f"real_tail_max_{t}"] = round(bt["max_abs"], 1)
        cat = stamps.get("racs_low3_catalogue")
        if cat:
            m["real_cat_name"] = cat["name"]
            m["real_cat_offset"] = round(float(cat["offset_arcsec"]), 1)
            m["real_cat_e_dec"] = round(float(cat["e_dec_arcsec"]), 1)
    return m


def _c3b_metrics(results_dir) -> dict[str, Any]:
    import json
    from pathlib import Path

    fp = Path(results_dir) / "ucdv_controls_c3b.json"
    if not fp.exists():
        return {}
    ep = json.loads(fp.read_text())["epochs"]
    m: dict[str, Any] = {}
    for e in ("low2", "low3"):
        r = ep.get(e, {})
        if r.get("status") == "ok":
            v, sv, i = float(r["v_peak"]), float(r["v_rms"]), float(r["i_peak"])
            m[f"real_c3b_{e}_v"] = round(v, 2)
            m[f"real_c3b_{e}_snr"] = round(v / sv, 1)
            m[f"real_c3b_{e}_isnr"] = round(i / float(r["i_rms"]), 1)
            m[f"real_c3b_{e}_frac"] = round(abs(v) / i, 1) if i > 0 else None
    return m


def run(out: str = ".", *, offline: bool = True) -> dict[str, Any]:
    """Slice entry point. Offline: synthetic census with planted emitters (recover-a-known).

    Real (``offline=False``): the committed census rows in ``data/ucdv/census.jsonl``. Offline and
    real metrics are namespaced (``syn_*`` / ``real_*`` -> ``\\ucdvSyn*`` / ``\\ucdvReal*``), the
    results JSON and macros merge through ``report.write_results`` / ``preserve_live_macros``,
    and an offline run must be pointed at a scratch ``out`` (CLAUDE.md).
    """
    import json
    from pathlib import Path

    from .report import write_results

    rows, dist, emit = synthetic_census()
    res = analyze_census(rows, dist)
    found = set(res["source_id"][res["candidate"]].tolist())
    metrics: dict[str, Any] = {"source": "synthetic census with planted emitters"}
    metrics.update(_metrics(res, "syn"))
    metrics["syn_n_planted"] = len(emit)
    metrics["syn_n_recovered"] = len(found & emit)
    metrics["syn_n_false"] = len(found - emit)
    real = None
    if not offline:  # pragma: no cover - needs the committed real control results
        metrics.update(_control_metrics(Path(out) / "results"))
    if not offline:  # pragma: no cover - needs the real census
        cpath = Path(out) / "data" / "ucdv" / "census.jsonl"
        crow = [json.loads(x) for x in cpath.read_text().splitlines()]
        tg = census_targets()
        dist_r = {int(s): 1000.0 / float(p) for s, p in zip(tg["Source"], tg["Plx"], strict=True)}
        floors = census_leak_floors(crow, Path(out) / "data" / "ucdv" / "census_cutouts")
        real = analyze_census(crow, dist_r, floors)
        metrics.update(_metrics(real, "real"))
        metrics["real_n_floor_computed"] = len(floors)
        metrics["real_held"] = [
            {"source_id": int(s), "epoch": str(e)}
            for s, e in zip(
                real["source_id"][real["held"]], real["epoch"][real["held"]], strict=True
            )
        ]
        metrics["real_candidates"] = [
            {"source_id": int(s), "epoch": str(e)}
            for s, e in zip(
                real["source_id"][real["candidate"]], real["epoch"][real["candidate"]], strict=True
            )
        ]
        metrics["real_null_n"] = int(real["null_monitor"].get("n", 0))
        metrics.update(_detection_metrics(Path(out) / "results", real, crow))
        metrics.update(_c3b_metrics(Path(out) / "results"))
        metrics.update(_tail_and_astrometry_metrics(Path(out) / "results", real, crow))
        metrics["source"] = "synthetic recover-a-known + real RACS-low2/low3 forced-V census"
    op = Path(out)
    write_results(metrics, op / "results" / "ucdv_metrics.json")
    fu = stamps = None
    fp = op / "results" / "ucdv_followup_j1123.json"
    sp = op / "results" / "ucdv_followup_j1123_stamps.json"
    if real is not None and fp.exists():
        fu = json.loads(fp.read_text())["epochs"]
    if real is not None and sp.exists():
        stamps = json.loads(sp.read_text())
    _figure(real if real is not None else res, op / "papers" / "ucdv" / "figures", fu, stamps)
    _write_macros(metrics, op / "papers" / "ucdv" / "generated" / "macros.tex")
    return metrics


def _figure(res: dict[str, Any], out_dir, followup=None, stamps=None) -> None:
    """(a) 95% limit on the fraction of snapshots above L_nu,V (5 sigma_V sensitivity);
    (b) with ``followup`` (the post-hoc light curve), V/sigma_V of the detected star in every
    ASKAP observation; otherwise (b) is the per-slot limit against distance. With ``stamps``
    (``results/ucdv_followup_j1123_stamps.json``), (c)-(d) are the Stokes-V images of the two
    detections, with the propagated position and the beam."""
    from pathlib import Path

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    od = Path(out_dir)
    od.mkdir(parents=True, exist_ok=True)
    s = res["sane"]
    sens, l_det = _sensitivity_and_detections(res, SENS_NSIGMA[0])
    frac = snapshot_fraction_limit(sens, L_GRID, l_det)
    if stamps:
        fig, axs = plt.subplots(2, 2, figsize=(9, 7.6))
        ax = axs.ravel()
    else:
        fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
        ax = axs.ravel()
    fin = np.isfinite(frac["frac_upper"]) & (frac["frac_upper"] <= 1.0)
    ax[0].plot(frac["l_grid"][fin], frac["frac_upper"][fin], color="k")
    for L in l_det:
        ax[0].axvline(L, color="tab:red", ls=":", lw=1)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel(r"$L_{\nu,V}$ (erg s$^{-1}$ Hz$^{-1}$)")
    ax[0].set_ylabel("95% upper limit on snapshot fraction")
    ax[0].set_title("(a)", loc="left", fontsize=9)
    from astropy.time import Time

    if followup:
        ok = [
            r
            for r in followup
            if r.get("status") == "ok" and np.isfinite(r["v"]) and np.isfinite(r["sigma_v"])
        ]
        t = Time([r["mjd"] for r in ok], format="mjd").decimalyear
        z = np.array([r["v"] / r["sigma_v"] for r in ok])
        census = np.array(["RACS_" in " ".join(r.get("files", [])) for r in ok])
        ax[1].axhline(0, color="0.6", lw=0.8)
        ax[1].axhspan(-2.5, 2.5, color="0.92")
        ax[1].scatter(t[~census], z[~census], s=14, color="0.2", label="other ASKAP projects")
        ax[1].scatter(t[census], z[census], s=22, marker="s", color="tab:red", label="RACS")
        ax[1].set_xlabel("year")
        ax[1].set_ylabel(r"$V/\sigma_V$ at the star")
        ax[1].legend(fontsize=7, loc="upper left")
    else:
        ax[1].scatter(res["dist_pc"][s], res["l_upper"][s], s=2, alpha=0.4, color="0.3")
        ax[1].set_xlabel("distance (pc)")
        ax[1].set_ylabel(r"$L_{\nu,V}$ upper limit (erg s$^{-1}$ Hz$^{-1}$)")
        ax[1].set_yscale("log")
    ax[1].set_title("(b)", loc="left", fontsize=9)
    if stamps:
        from matplotlib.patches import Ellipse

        labels = {"ASKAP-56518": "RACS-low3 (census)", "ASKAP-82408": "EMU (follow-up)"}
        for a, (oid, st), tag in zip(
            ax[2:], stamps["stamps"].items(), ("(c)", "(d)"), strict=False
        ):
            v = np.asarray(st["v_over_sigma"], float)
            hw = (v.shape[0] - 1) / 2 * st["pix_arcsec"]
            ext = (hw, -hw, -hw, hw)  # east to the left
            im = a.imshow(v, origin="lower", cmap="RdBu_r", vmin=-6, vmax=6, extent=ext)
            dx, dy = (o * st["pix_arcsec"] for o in st["star_offset_pix"])
            a.plot(-dx, dy, marker="+", color="k", ms=14, mew=1.4)
            bmaj, bmin, bpa = st["beam_arcsec_deg"]
            a.add_patch(
                Ellipse((0.75 * hw, -0.75 * hw), bmin, bmaj, angle=-bpa, fill=False, color="k")
            )
            yr = Time(st["mjd"], format="mjd").datetime.strftime("%Y %b")
            a.set_title(f"{tag} Stokes V, {labels.get(oid, oid)}, {yr}", loc="left", fontsize=8)
            a.set_xlabel(r"$\Delta\alpha$ (arcsec)")
            a.set_ylabel(r"$\Delta\delta$ (arcsec)")
            cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.02)
            if tag == "(d)":
                cb.set_label(r"$V/\sigma_V$")
    fig.tight_layout(w_pad=2.5)
    fig.savefig(od / "ucdv_limits.pdf")
    plt.close(fig)


_MACROS = {
    "SynNTargets": "syn_n_targets",
    "SynNPlanted": "syn_n_planted",
    "SynNRecovered": "syn_n_recovered",
    "SynNFalse": "syn_n_false",
    "RealNReyle": "real_n_reyle",
    "RealNBoth": "real_n_both_epochs",
    "RealDmax": "real_d_max_pc",
    "RealNCensus": "real_n_census",
    "RealFullHours": "real_full_hours",
    "RealSzeroOk": "real_s0_n_ok",
    "RealSzeroN": "real_s0_n_attempts",
    "RealSzeroSec": "real_s0_sec_per_pair",
    "RealConeN": "real_c1_n",
    "RealConePos": "real_c1_frac_pos",
    "RealConeStd": "real_c1_std",
    "RealCtwoN": "real_c2_n_per_level",
    "RealCtwoComp": "real_c2_comp7",
    "RealCtwoRatio": "real_c2_ratio",
    "RealCthreeN": "real_c3_n",
    "RealCthreeRatio": "real_c3_ratio",
    "RealCthreeFlips": "real_c3_flips",
    "RealNTargets": "real_n_targets",
    "RealNSlots": "real_n_slots",
    "RealNSane": "real_n_sane",
    "RealNCand": "real_n_candidates",
    "RealNHeld": "real_n_held",
    "RealNPairs": "real_n_pairs",
    "RealNBeyondSpec": "real_n_beyond_spectral",
    "RealSigmaV": "real_median_sigma_v_mjy",
    "RealLimit": "real_median_limit_mjy",
    "RealNullPos": "real_null_frac_positive",
    "RealNullStd": "real_null_std",
    "RealNeffThreeFourteen": "real_neff_3e14",
    "RealFracThreeFourteen": "real_frac_upper_3e14",
    "RealNdetThreeFourteen": "real_ndet_3e14",
    "RealFracSevenFifteen": "real_frac7_upper_1e15",
    "RealNeffSevenFifteen": "real_neff7_1e15",
    "RealNDetections": "real_n_detections",
    "RealNullN": "real_null_n",
    "RealNMissing": "real_n_missing",
    "RealNInsane": "real_n_insane",
    "RealDetGaia": "real_det_gaia",
    "RealDetObs": "real_det_obs",
    "RealDetDate": "real_det_date",
    "RealDetV": "real_det_v",
    "RealDetSigmaV": "real_det_sigma_v",
    "RealDetSnr": "real_det_snr",
    "RealDetI": "real_det_i",
    "RealDetIsnr": "real_det_isnr",
    "RealDetFrac": "real_det_frac",
    "RealDetFracErr": "real_det_frac_err",
    "RealDetDist": "real_det_dist_pc",
    "RealDetL": "real_det_l_tex",
    "RealDetOtherSnr": "real_det_other_snr",
    "RealFuNObs": "real_fu_n_obs",
    "RealFuNRefused": "real_fu_n_refused",
    "RealFuNMeas": "real_fu_n_meas",
    "RealFuNIndep": "real_fu_n_indep",
    "RealFuFirstYear": "real_fu_first_year",
    "RealFuLastYear": "real_fu_last_year",
    "RealFuNPos": "real_fu_n_pos",
    "RealFuNNeg": "real_fu_n_neg",
    "RealFuBestObs": "real_fu_best_obs",
    "RealFuBestDate": "real_fu_best_date",
    "RealFuBestV": "real_fu_best_v",
    "RealFuBestSigmaV": "real_fu_best_sigma_v",
    "RealFuBestSnr": "real_fu_best_snr",
    "RealFuBestIsnr": "real_fu_best_isnr",
    "RealFuBestP": "real_fu_best_p_tex",
    "RealFuBestPTrials": "real_fu_best_p_trials_tex",
    "RealFuSubN": "real_fu_sub_n",
    "RealFuSubMean": "real_fu_sub_mean",
    "RealCthreebLowtwoSnr": "real_c3b_low2_snr",
    "RealCthreebLowthreeSnr": "real_c3b_low3_snr",
    "RealNullNGeFour": "real_null_n_ge4",
    "RealNullExpGeFour": "real_null_exp_ge4",
    "RealNullRatioGeFour": "real_null_ratio_ge4",
    "RealNullNGeFive": "real_null_n_ge5",
    "RealChanceGeFive": "real_chance_ge5",
    "RealNeffCwFifteen": "real_neff_cw_1e15",
    "RealFracCwFifteen": "real_frac_cw_upper_1e15",
    "RealCenVRacs": "real_cen_v_racs",
    "RealCenVEmu": "real_cen_v_emu",
    "RealCenIRacs": "real_cen_i_racs",
    "RealCenIEmu": "real_cen_i_emu",
    "RealTailNRacs": "real_tail_n_racs",
    "RealTailMaxRacs": "real_tail_max_racs",
    "RealTailNEmu": "real_tail_n_emu",
    "RealTailMaxEmu": "real_tail_max_emu",
    "RealCatName": "real_cat_name",
    "RealCatOffset": "real_cat_offset",
    "RealCatEDec": "real_cat_e_dec",
    "RealDetPairZ": "real_det_pair_z",
    "RealTwoEpochFloor": "real_two_epoch_floor_mjy",
    "RealFuMatchObs": "real_fu_match_obs",
    "RealFuMatchYear": "real_fu_match_year",
    "RealFuMatchV": "real_fu_match_v",
    "RealFuMatchSigmaV": "real_fu_match_sigma_v",
    "RealFuMatchZ": "real_fu_match_z",
    "RealFuNVgtI": "real_fu_n_v_gt_i",
    "RealFuVgtIObs": "real_fu_v_gt_i_obs",
    "RealFuNNan": "real_fu_n_nan",
    "RealCthreebLowtwoV": "real_c3b_low2_v",
    "RealCthreebLowthreeV": "real_c3b_low3_v",
    "RealCthreebLowtwoFrac": "real_c3b_low2_frac",
    "RealCthreebLowthreeFrac": "real_c3b_low3_frac",
    "RealNeffFifteen": "real_neff_1e15",
    "RealFracFifteen": "real_frac_upper_1e15",
    "RealNeffSixteen": "real_neff_1e16",
    "RealFracSixteen": "real_frac_upper_1e16",
}


def _write_macros(m: dict[str, Any], path) -> None:
    from pathlib import Path

    from .report import preserve_live_macros

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.ucdv._write_macros -- do not edit.",
        rf"\newcommand{{\ucdvSource}}{{{m['source']}}}",
    ]
    lines += [rf"\newcommand{{\ucdv{name}}}{{{_fmt(key)}}}" for name, key in _MACROS.items()]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import os

    ap = argparse.ArgumentParser(prog="python -m jansky_research.ucdv")
    ap.add_argument("command", choices=["census", "run"])
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=None, help="run: output root (use a scratch dir offline)")
    ap.add_argument("--offline", action="store_true")
    a = ap.parse_args(argv)
    if a.command == "run":
        if a.out is None:
            raise SystemExit("run: pass --out (offline runs must not write into the repo root)")
        import json

        print(json.dumps(run(a.out, offline=a.offline), indent=2, default=str))
        return 0
    user = os.environ.get("CASDA_USERNAME")
    if not user:
        raise SystemExit("set CASDA_USERNAME")
    n = run_census(username=user, limit=a.limit)
    print(f"census: {n} slots attempted this run")
    return 0


# ---- analysis of the census (runs after the fetch, never inside it) -----------------------------
BRIGHT_LEAK_SNR = 50.0  # field sources used for the leakage floor: I/sigma_I above this
MIN_LEAK_SOURCES = 3


def cutout_leakage_floor(image_i, image_v, *, snr: float = BRIGHT_LEAK_SNR) -> float:
    """Leakage floor from the cutout's own bright Stokes-I peaks (``stokesv.leakage_floor``).

    Local maxima with ``I > snr * sigma_I`` (sigma from a robust MAD of the whole cutout) give
    |V/I|; fewer than ``MIN_LEAK_SOURCES`` such peaks returns NaN, which the detection rule treats
    as "floor undetermined" -- the candidate is held for vetting, never passed or dropped silently.
    """
    from scipy.ndimage import maximum_filter

    from .stokesv import leakage_floor

    i = np.asarray(image_i, float)
    v = np.asarray(image_v, float)
    ok = np.isfinite(i)
    if not ok.any():
        return float("nan")
    sig = 1.4826 * float(np.nanmedian(np.abs(i[ok] - np.nanmedian(i[ok]))))
    peaks = (i == maximum_filter(np.nan_to_num(i, nan=-np.inf), size=5)) & (i > snr * sig)
    if int(peaks.sum()) < MIN_LEAK_SOURCES:
        return float("nan")
    return leakage_floor(np.abs(v[peaks]) / i[peaks])


def poisson_upper(n: int, cl: float = 0.95) -> float:
    """One-sided Poisson upper limit on the mean for ``n`` observed events (2.996 for n = 0)."""
    return float(stats.chi2.ppf(cl, 2 * (int(n) + 1)) / 2.0)


def snapshot_fraction_limit(
    l_limit: np.ndarray, l_grid: np.ndarray, l_detected: np.ndarray | None = None, cl: float = 0.95
) -> dict[str, np.ndarray]:
    """Upper limit on the fraction of target-snapshots above luminosity ``L``, for each ``L``.

    Only snapshots whose own limit reaches ``L`` could have seen a source that bright, so the
    denominator is ``N_eff(L) = #(l_limit <= L)`` -- not the sample size (the ``frblens`` lesson:
    a null divides by sensitivity). Detections with ``L_det >= L`` enter the Poisson numerator.
    """
    lim = np.asarray(l_limit, float)
    lim = lim[np.isfinite(lim)]
    det = np.asarray([] if l_detected is None else l_detected, float)
    grid = np.asarray(l_grid, float)
    n_eff = np.array([int(np.sum(lim <= L)) for L in grid])
    n_det = np.array([int(np.sum(det >= L)) for L in grid])
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.array(
            [
                poisson_upper(k, cl) / n if n > 0 else np.inf
                for k, n in zip(n_det, n_eff, strict=True)
            ]
        )
    return {"l_grid": grid, "n_eff": n_eff, "n_det": n_det, "frac_upper": frac}


def analyze_census(
    rows: list[dict[str, Any]],
    dist_pc: dict[int, float],
    leak_floor: dict[tuple[int, str], float] | None = None,
) -> dict[str, Any]:
    """Apply the frozen rules to census rows (one per target-epoch, as ``run_census`` writes them).

    Slot sanity: rows with undefined sigma_V, or sigma_V > ``CUTOUT_RMS_MAX_FACTOR`` x that
    epoch's median, get no limit. Candidates follow :func:`apply_detection_rule`; when Stokes I
    is detected and the leakage floor is undetermined the slot is ``held`` for vetting. Returns
    per-slot arrays, the per-target two-epoch comparison, the ongoing C1 null monitor, and counts.
    """
    leak_floor = leak_floor or {}
    ok = [r for r in rows if r.get("status") == "ok"]
    sid = np.array([int(r["source_id"]) for r in ok], dtype=np.int64)
    ep = np.array([str(r["epoch"]) for r in ok])
    v = np.array([float(r["v"]) for r in ok])
    sv = np.array([float(r["sigma_v"]) for r in ok])
    i = np.array([float(r["i"]) for r in ok])
    si = np.array([float(r["sigma_i"]) for r in ok])
    bright = np.array([bool(r["bright"]) for r in ok])
    sane = np.isfinite(sv) & (sv > 0)
    for e in np.unique(ep):
        m = ep == e
        med = float(np.nanmedian(sv[m & sane])) if (m & sane).any() else np.nan
        sane &= ~(m & (sv > CUTOUT_RMS_MAX_FACTOR * med))
    floor = np.array(
        [leak_floor.get((int(a), str(b)), np.nan) for a, b in zip(sid, ep, strict=True)]
    )
    rule = apply_detection_rule(v, sv, i, si, np.nan_to_num(floor, nan=np.inf), bright)
    i_det = i / si >= SNR_DETECT
    held = rule["v_significant"] & i_det & ~np.isfinite(floor) & rule["clear_of_bright"]
    cand = rule["candidate"] & sane
    ul = np.where(sane, upper_limit(v, sv), np.nan)
    d = np.array([dist_pc.get(int(s), np.nan) for s in sid])
    l_ul = luminosity_limit(ul, d)
    nv = (
        np.concatenate([np.asarray(r.get("null_v", []), float) for r in ok]) if ok else np.array([])
    )
    ns = (
        np.concatenate([np.asarray(r.get("null_sigma_v", []), float) for r in ok])
        if ok
        else np.array([])
    )
    # two-epoch comparison for targets with both epochs sane
    pairs = []
    for s in np.unique(sid):
        a = np.flatnonzero((sid == s) & (ep == "low2") & sane)
        b = np.flatnonzero((sid == s) & (ep == "low3") & sane)
        if a.size and b.size:
            pairs.append((int(s), int(a[0]), int(b[0])))
    if pairs:
        ia = np.array([p[1] for p in pairs])
        ib = np.array([p[2] for p in pairs])
        cmp = two_epoch_compare(v[ia], sv[ia], v[ib], sv[ib])
    else:
        cmp = {
            "z_diff": np.array([]),
            "differs": np.array([], bool),
            "beyond_spectral": np.array([]),
        }
    return {
        "source_id": sid,
        "epoch": ep,
        "v": v,
        "sigma_v": sv,
        "sane": sane,
        "candidate": cand,
        "held": held & sane,
        "upper_limit_mjy": ul,
        "dist_pc": d,
        "l_upper": l_ul,
        "pair_source_id": np.array([p[0] for p in pairs], dtype=np.int64),
        "two_epoch": cmp,
        "null_monitor": null_check(nv, ns) if nv.size else {"n": 0, "passed": False},
        "n_targets": int(np.unique(sid).size),
        "n_slots": int(len(ok)),
        "n_sane": int(sane.sum()),
        "n_candidates": int(cand.sum()),
        "n_held": int((held & sane).sum()),
        "n_pairs": len(pairs),
    }


def census_leak_floors(rows: list[dict[str, Any]], cutout_dir) -> dict[tuple[int, str], float]:
    """Per-slot leakage floor (:func:`cutout_leakage_floor`) from the saved census cutouts.

    Computed only where the frozen detection rule consults it -- |V|/sigma_V >= SNR_DETECT with
    Stokes I detected (I/sigma_I >= SNR_DETECT); elsewhere the floor cannot change the outcome.
    A slot whose cutout file is missing gets no entry (the rule then holds it, never passes it).
    """
    from pathlib import Path

    out: dict[tuple[int, str], float] = {}
    for r in rows:
        if r.get("status") != "ok":
            continue
        v, sv, i, si = (float(r[k]) for k in ("v", "sigma_v", "i", "sigma_i"))
        if not (sv > 0 and si > 0 and abs(v) / sv >= SNR_DETECT and i / si >= SNR_DETECT):
            continue
        f = Path(cutout_dir) / f"{int(r['source_id'])}_{r['epoch']}.npz"
        if f.is_file():
            c = load_cutout(f)
            out[(int(r["source_id"]), str(r["epoch"]))] = cutout_leakage_floor(c["i"], c["v"])
    return out


def synthetic_census(
    n_targets: int = 2000, *, n_emitters: int = 6, seed: int = SEED
) -> tuple[list[dict[str, Any]], dict[int, float], set[int]]:
    """Offline fixture shaped like ``census.jsonl``: noise-only targets plus planted emitters.

    Noise per slot is drawn around the measured RACS level (median sigma_V 0.15 mJy/beam, S0);
    ``n_emitters`` targets carry a strongly circularly polarised source (|V| = 10-20 sigma, I below
    detection) in one or both epochs. Returns ``(rows, dist_pc, emitter_ids)``.
    """
    rng = np.random.default_rng(seed)
    ids = rng.choice(np.arange(10**15, 10**15 + 10**8, 2), n_targets, replace=False).astype(
        np.int64
    )
    dist = dict(zip(ids.tolist(), rng.uniform(10.0, 107.8, n_targets).tolist(), strict=True))
    emit = set(rng.choice(ids, n_emitters, replace=False).tolist())
    rows = []
    for s in ids.tolist():
        both = rng.random() < 0.5
        for k, e in enumerate(("low2", "low3")):
            sv = float(rng.lognormal(np.log(0.15), 0.2))
            si = sv * 1.3
            vt = 0.0
            if s in emit and (both or k == 0):
                vt = float(rng.choice([-1.0, 1.0]) * rng.uniform(10, 20) * sv)
            rows.append(
                {
                    "source_id": s,
                    "epoch": e,
                    "status": "ok",
                    "v": vt + rng.normal(0, sv),
                    "sigma_v": sv,
                    "i": abs(vt) * 0.3 + rng.normal(0, si),
                    "sigma_i": si,
                    "bright": bool(rng.random() < 0.01),
                    "null_v": rng.normal(0, sv, NULL_PER_TARGET).tolist(),
                    "null_sigma_v": [sv] * NULL_PER_TARGET,
                }
            )
    return rows, dist, emit


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
