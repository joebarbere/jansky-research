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
    "assign_set",
    "forced_measure",
    "in_footprint",
    "parse_polygon",
    "pick_tile",
    "run_s0",
    "run_c3",
    "run_census",
    "backoff_seconds",
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
        raise ValueError(f"CASDA cutout job ended {status}: {job}")
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
):  # pragma: no cover - network
    """One SODA job staging a same-observation I+V cutout pair for ``epoch`` at ``(ra, dec)``.

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

    lam = EPOCH_LAMBDA_M[epoch]
    t0 = time.time()
    q = f"""
    SELECT obs_id, filename, t_min, access_url, s_region FROM ivoa.obscore
    WHERE obs_collection = 'The Rapid ASKAP Continuum Survey'
    AND (filename LIKE 'image.i.%' OR filename LIKE 'image.v.%')
    AND filename LIKE '%.taylor.0.restored.conv.fits'
    AND em_min > {lam - 0.01:.3f} AND em_max < {lam + 0.01:.3f}
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
    urls = _cutout_urls(
        casda, rows["access_url",], SkyCoord(ra * u.deg, dec * u.deg), radius_deg, poll_interval
    )
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
    t_dl = time.time() - t2
    if "i" not in out or "v" not in out:
        raise RuntimeError(f"cutout job returned {sorted(out)} for {urls}")
    if np.shape(out["i"]) != np.shape(out["v"]):
        raise RuntimeError("I/V cutout grids differ")
    out.update(
        mjd=float(rows["t_min"][0]),
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
    import json
    import time
    from pathlib import Path

    from .stokesv import _casda_session
    from .stokesv_discovery import epoch_position

    tg = census_targets()
    done: set[tuple[int, str]] = set()
    out = Path(out_jsonl)
    if out.exists():
        for line in out.read_text().splitlines():
            r = json.loads(line)
            done.add((int(r["source_id"]), str(r["epoch"])))
    Path(cutout_dir).mkdir(parents=True, exist_ok=True)
    casda = _casda_session(username, pw_path)
    n = 0
    fails = 0  # consecutive failed slots, drives the backoff
    with out.open("a") as fh:
        for r in tg:
            sid = int(r["Source"])
            for ep in ("low2", "low3"):
                if (sid, ep) in done:
                    continue
                if limit is not None and n >= limit:
                    return n
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
                except Exception as exc:  # recorded and retried on the next resume
                    rec.update(status="error", error=f"{type(exc).__name__}: {exc}"[:300])
                    fails += 1
                    # Back off exponentially (cap 30 min) so an outage costs tens of logins,
                    # not hundreds: on 2026-10-16 a fixed 60 s retry made ~285 failed OPAL logins
                    # over 11.5 h while CASDA authentication was refusing.
                    wait = backoff_seconds(fails)
                    print(
                        time.strftime("%Y-%m-%d %H:%M:%S"),
                        f"slot error #{fails} ({rec['error'][:120]}); sleeping {wait:.0f} s",
                        flush=True,
                    )
                    time.sleep(wait)
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
    return n


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import os

    ap = argparse.ArgumentParser(prog="python -m jansky_research.ucdv")
    ap.add_argument("command", choices=["census"])
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args(argv)
    user = os.environ.get("CASDA_USERNAME")
    if not user:
        raise SystemExit("set CASDA_USERNAME")
    n = run_census(username=user, limit=a.limit)
    print(f"census: {n} slots attempted this run")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
