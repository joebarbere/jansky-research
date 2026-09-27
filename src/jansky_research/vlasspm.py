"""Blind radio proper motions across the VLASS epochs --- movers found without Gaia (plan 64).

Every published radio proper motion starts from a known object (the field is Gaia-anchored:
De et al. 2024, arXiv:2409.18466; Driessen et al. 2023, arXiv:2306.08059). This module asks the
blind question: which VLASS components *move* between epochs, found from the radio catalogues
alone, with optical/IR counterparts excluded only at the very end?

The method, and the simplification that sizes it
------------------------------------------------
A source moving at ``mu`` between epochs a time ``dt`` apart shifts by ``mu * dt``. At the plan's
upper rate of 5"/yr and VLASS's ~3.5-yr E1->E2 baseline that is < 20", so linkage is a
fixed-radius neighbour search, not an all-pairs problem: a KD-tree on unit vectors does E1 x E2
over the whole survey in minutes on a CPU. The compute that matters is in the *null* (many
position-scrambled re-linkages) and the *completeness* (many injection trials), which are
embarrassingly repeatable and checkpointed.

1. **Orphans.** A mover leaves its E1 position empty in E2. Each epoch's components with no
   counterpart in the other epoch within the static-match radius are orphans; only orphans are
   linked. Static sources never enter the search.
2. **E1 x E2 linkage.** Orphan pairs whose implied rate ``sep / dt`` lies in [mu_min, mu_max].
   These are dominated by chance: variable sources near the detection limit orphan themselves.
3. **E3 collinearity.** The pair predicts where the source must be at each E3 component's own
   epoch; an E3 orphan within ``tol`` of that prediction confirms a three-point straight line at a
   constant rate. No candidate survives on two epochs alone.
4. **Flux consistency** across the three detections (a cheap cut against unrelated orphans).
5. **The null.** Rotating the later epochs' orphans in RA by a few arcminutes destroys every real
   mover while preserving the orphan density, so re-running steps 2-3 on scrambles gives the
   expected number of chance pairs and chance triplets (Poisson), per rate bin.
6. **Completeness.** Movers injected into the real catalogues, run through the same pipeline.
7. **Limit.** At zero survivors, the surface-density upper limit divides by the
   *completeness-weighted* area, not the raw area (the ``frblens`` lesson: a null divides by
   sensitivity, not sample size).

Only after 1-6 are Gaia/CatWISE counterparts removed; the survivors are optically dark movers
(Y dwarfs, high-velocity pulsars) --- or, at zero yield, the limit is the result.

Pure NumPy/SciPy; the offline fixture (:func:`synthetic_epochs`) plants movers, statics and
variable-source orphans with the real astrometric floor so every step is testable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

__all__ = [
    "EpochCatalog",
    "ErrorModel",
    "SizeNoiseModel",
    "choose_threshold",
    "compact_keep",
    "compactness",
    "compactness_keep_table",
    "injection_compactness",
    "abs_cov",
    "calibration_table",
    "chi2_2d",
    "condon_cov",
    "condon_position_errors",
    "condon_rho",
    "cov_axes",
    "ellipse_cov",
    "fit_error_model",
    "match_statics",
    "static_pair_terms",
    "collinearity_test",
    "banded_floor",
    "calibrate_floors",
    "calibrate_floors_banded",
    "completeness",
    "epoch_triples",
    "flux_consistent",
    "inject_movers",
    "isolated_mask",
    "link_pairs",
    "orphan_masks",
    "run",
    "scramble_null",
    "search",
    "search_multi",
    "surface_density_limit",
    "synthetic_epochs",
    "synthetic_statics",
    "chance_within",
    "comoving",
    "fit_track",
    "parallax_factors",
    "parallax_residual_factor",
    "predict_track",
    "tangent_offsets_arcsec",
    "triplet_significance",
]

ARCSEC = 1.0 / 3600.0
MU_MIN_ARCSEC_YR = 0.3  # plan 64 annulus; its low edge sits at the E1 astrometric floor
MU_MAX_ARCSEC_YR = 5.0
STATIC_MATCH_ARCSEC = 2.5  # one VLASS beam: a counterpart this close means "did not move"
E3_TOL_SIGMA = 3.0  # E3 must lie within this many combined sigma of the predicted position
FLUX_RATIO_MAX: float | None = None  # max/min peak flux across detections; None = no cut.
# Off by default since the first real run: radio stars FLARE (UV Ceti: 1.9, 11.0, 1.2 mJy across
# E2-E4, ratio 9.5), so a flux-consistency cut rejects exactly the population being searched for.
ISOLATION_ARCSEC = 30.0  # an orphan must have no other component this close IN ITS OWN EPOCH.
# The first real run's candidates had a same-epoch neighbour within 30" 86% of the time vs 25%
# for a random component: extended sources decomposed differently per epoch produce spatially
# correlated orphans that line up, and an arcminute RA-scramble null cannot see them.
POISSON95_ZERO = 2.996  # 95% one-sided upper limit on a Poisson mean when 0 events are seen
DAYS_PER_YEAR = 365.25


@dataclass
class EpochCatalog:
    """One epoch's component list. ``t_yr`` is per component (VLASS tiles span years).

    The shape-aware fields are optional; when absent the catalogue behaves exactly as the
    isotropic ``pos_err`` model of runs 1-3:

    - ``cov`` (n, 3): measurement covariance ``[xx, yy, xy]`` in arcsec^2 on the sky
      (x = RA*cos(dec), east; y = Dec, north) --- the Condon (1997) error ellipse plus the
      astrometric floor squared. Default: ``pos_err^2`` on both axes, no correlation.
    - ``shape`` (n, 3): the deconvolved source's FWHM^2 covariance (arcsec^2); zeros = point.
    - ``beam`` (n, 3): the restoring beam's FWHM^2 covariance (arcsec^2); zeros = unknown.
    - ``rms`` (n,): local noise, mJy/beam (NaN = unknown). ``floor`` (n,): floor, arcsec.
    """

    ra: np.ndarray  # deg
    dec: np.ndarray  # deg
    t_yr: np.ndarray  # decimal year of the observation of each component
    flux: np.ndarray  # peak flux, mJy/beam (scale-corrected)
    pos_err: np.ndarray  # 1-sigma per-axis positional error, arcsec (isotropic equivalent)
    ident: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    cov: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    shape: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    beam: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    rms: np.ndarray = field(default_factory=lambda: np.zeros(0))
    floor: np.ndarray = field(default_factory=lambda: np.zeros(0))

    def __post_init__(self) -> None:
        for k in ("ra", "dec", "t_yr", "flux", "pos_err"):
            setattr(self, k, np.asarray(getattr(self, k), dtype=float))
        n = self.ra.size
        self.ident = np.asarray(self.ident)
        if self.ident.size == 0:
            self.ident = np.full(n, -1, dtype=np.int64)
        self.ident = np.asarray(self.ident)
        if np.asarray(self.cov).size == 0:
            e2 = self.pos_err**2
            self.cov = np.column_stack([e2, e2, np.zeros(n)])
        for k in ("shape", "beam"):
            if np.asarray(getattr(self, k)).size == 0:
                setattr(self, k, np.zeros((n, 3)))
        if np.asarray(self.rms).size == 0:
            self.rms = np.full(n, np.nan)
        if np.asarray(self.floor).size == 0:
            self.floor = np.zeros(n)
        for k in ("cov", "shape", "beam"):
            setattr(self, k, np.asarray(getattr(self, k), float).reshape(n, 3))
        self.rms = np.asarray(self.rms, float)
        self.floor = np.asarray(self.floor, float)

    def __len__(self) -> int:
        return int(self.ra.size)

    @property
    def has_shape(self) -> bool:
        """True when the catalogue carries per-component noise and beams (the real epochs)."""
        return bool(np.isfinite(self.rms).any() and np.any(self.beam[:, 0] > 0))

    def subset(self, mask: np.ndarray) -> EpochCatalog:
        return EpochCatalog(*(getattr(self, k)[mask] for k in _FIELDS))

    @staticmethod
    def concat(parts: list[EpochCatalog]) -> EpochCatalog:
        return EpochCatalog(
            *(np.concatenate([getattr(p, k) for p in parts]) for k in _FIELDS),
        )


_FIELDS = (
    "ra",
    "dec",
    "t_yr",
    "flux",
    "pos_err",
    "ident",
    "cov",
    "shape",
    "beam",
    "rms",
    "floor",
)


# ------------------------------------------------------------------ shape-aware error model
#
# Condon (1997, PASP 109, 166), eq. 21 with the effective S/N of Condon et al. (1998, AJ 115,
# 1693) eq. 26 (exponents cross-checked against Prandoni et al. 2000, A&AS 146, 41, App. A, and
# against PyBDSF's ``functions.get_errors``, which produced every VLASS catalogue used here):
#
#   rho^2 = (theta_M theta_m / 4 theta_N^2) (1 + theta_N^2/theta_M^2)^aM
#                                           (1 + theta_N^2/theta_m^2)^am  (A/mu)^2
#   8 ln2 sigma^2(x0) / theta_M^2 = 2 / rho^2   with (aM, am) = (5/2, 1/2)   [along the major axis]
#   8 ln2 sigma^2(y0) / theta_m^2 = 2 / rho^2   with (aM, am) = (1/2, 5/2)   [along the minor axis]
#
# theta_M, theta_m are the FITTED (beam-convolved) FWHM axes, theta_N^2 = bmaj * bmin is the
# noise correlation area (PyBDSF's choice), A the peak and mu the local rms. For a point source
# in a round beam (theta_M = theta_m = theta_N) rho^2 = 2 SNR^2 and sigma(x0) = sigma_beam / SNR
# (sigma_beam = FWHM / sqrt(8 ln2)); a resolved source (theta_M > beam) gets a LARGER error at the
# same peak S/N. All four VLASS catalogues' E_RA/E_DEC equal this to 1 part in 10^3.

FWHM_PER_SIGMA = float(np.sqrt(8.0 * np.log(2.0)))


def condon_rho(maj, min_, bmaj, bmin, snr, alpha_major: float, alpha_minor: float) -> np.ndarray:
    """Condon's effective signal-to-noise ``rho`` (not squared) for one parameter's exponents."""
    maj, min_, snr = (np.asarray(x, float) for x in (maj, min_, snr))
    th2 = np.asarray(bmaj, float) * np.asarray(bmin, float)
    rho2 = (
        maj
        * min_
        / (4.0 * th2)
        * (1.0 + th2 / maj**2) ** alpha_major
        * (1.0 + th2 / min_**2) ** alpha_minor
        * snr**2
    )
    return np.sqrt(rho2)


def condon_position_errors(maj, min_, bmaj, bmin, snr) -> tuple[np.ndarray, np.ndarray]:
    """1-sigma centroid errors (arcsec) ALONG the fitted major and minor axes (Condon 1997)."""
    k = np.sqrt(2.0) / FWHM_PER_SIGMA
    s_major = np.asarray(maj, float) * k / condon_rho(maj, min_, bmaj, bmin, snr, 2.5, 0.5)
    s_minor = np.asarray(min_, float) * k / condon_rho(maj, min_, bmaj, bmin, snr, 0.5, 2.5)
    return s_major, s_minor


def ellipse_cov(a, b, pa_deg) -> np.ndarray:
    """(n, 3) sky covariance ``[xx, yy, xy]`` of an ellipse with 1-sigma (or FWHM) semi-values
    ``a`` along the major axis and ``b`` along the minor, major axis at ``pa_deg`` east of north.

    x = east (RA*cos dec), y = north. The major axis points along (sin PA, cos PA).
    """
    a2, b2 = np.asarray(a, float) ** 2, np.asarray(b, float) ** 2
    p = np.radians(np.asarray(pa_deg, float))
    s, c = np.sin(p), np.cos(p)
    return np.column_stack(
        np.broadcast_arrays(a2 * s * s + b2 * c * c, a2 * c * c + b2 * s * s, (a2 - b2) * s * c)
    )


def cov_axes(cov: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse of :func:`ellipse_cov`: (major, minor, pa_deg) of an (n, 3) covariance."""
    cov = np.atleast_2d(cov)
    xx, yy, xy = cov[:, 0], cov[:, 1], cov[:, 2]
    tr, d = 0.5 * (xx + yy), np.sqrt((0.5 * (yy - xx)) ** 2 + xy**2)
    lam1, lam2 = tr + d, np.maximum(tr - d, 0.0)
    pa = np.degrees(0.5 * np.arctan2(2.0 * xy, yy - xx))  # of the major axis, E of N
    return np.sqrt(lam1), np.sqrt(lam2), pa


def condon_cov(maj, min_, pa_deg, bmaj, bmin, snr) -> np.ndarray:
    """(n, 3) Condon (1997) position covariance on the sky, rotated by the fitted PA."""
    s_major, s_minor = condon_position_errors(maj, min_, bmaj, bmin, snr)
    return ellipse_cov(s_major, s_minor, pa_deg)


def abs_cov(cov: np.ndarray) -> np.ndarray:
    """Matrix absolute value |M| of symmetric 2x2 matrices stored as (n, 3) ``[xx, yy, xy]``.

    Used for the beam-change term: a beam that got LONGER in one direction and SHORTER in the
    other moves a centroid along both, so the variance is |B_a - B_b|, never negative.
    """
    cov = np.atleast_2d(np.asarray(cov, float))
    xx, yy, xy = cov[:, 0], cov[:, 1], cov[:, 2]
    tr, d = 0.5 * (xx + yy), np.sqrt((0.5 * (xx - yy)) ** 2 + xy**2)
    l1, l2 = np.abs(tr + d), np.abs(tr - d)
    mean, half = 0.5 * (l1 + l2), 0.5 * (l1 - l2)
    with np.errstate(divide="ignore", invalid="ignore"):
        ux = np.where(d > 0, (0.5 * (xx - yy)) / d, 0.0)
        uxy = np.where(d > 0, xy / d, 0.0)
    return np.column_stack([mean + half * ux, mean - half * ux, half * uxy])


def chi2_2d(dx, dy, cov: np.ndarray) -> np.ndarray:
    """Mahalanobis distance squared of offsets (dx, dy) under (n, 3) covariances (chi^2, 2 dof)."""
    cov = np.atleast_2d(cov)
    xx, yy, xy = cov[:, 0], cov[:, 1], cov[:, 2]
    det = xx * yy - xy**2
    dx, dy = np.asarray(dx, float), np.asarray(dy, float)
    return (yy * dx * dx + xx * dy * dy - 2.0 * xy * dx * dy) / det


@dataclass(frozen=True)
class ErrorModel:
    """Systematic terms added to each component's measurement covariance.

    ``sys_cov = k_struct^2 * shape + q_beam^2 * |beam - beam_ref|``

    - ``k_struct`` (dimensionless): a resolved source's fitted centroid wanders by this fraction
      of its deconvolved FWHM, along each of its own axes (structure the single-Gaussian fit
      does not describe, redistributed differently by each epoch's imaging).
    - ``q_beam`` (dimensionless): extra variance per unit beam CHANGE (FWHM^2, matrix |.|)
      relative to ``beam_ref``, the mean beam over the epochs being compared.

    Both are measured from static sources (:func:`fit_error_model`); the zero model reproduces
    the run 1-3 behaviour (measurement covariance only).
    """

    k_struct: float = 0.0
    q_beam: float = 0.0

    def sys_cov(self, shape: np.ndarray, beam: np.ndarray, beam_ref: np.ndarray) -> np.ndarray:
        out = self.k_struct**2 * np.atleast_2d(shape)
        if self.q_beam:
            out = out + self.q_beam**2 * abs_cov(np.atleast_2d(beam) - np.atleast_2d(beam_ref))
        return out


NO_MODEL = ErrorModel()


def _xyz(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    r, d = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])


def _chord(arcsec: float) -> float:
    """Chord length on the unit sphere for an angular separation (for KD-tree radii)."""
    return float(2.0 * np.sin(np.radians(arcsec / 3600.0) / 2.0))


def tangent_offsets_arcsec(ra0, dec0, ra1, dec1) -> tuple[np.ndarray, np.ndarray]:
    """(dRA*cos(dec), dDec) in arcsec from (ra0, dec0) to (ra1, dec1); small-angle, RA-wrap safe."""
    dra = (np.asarray(ra1, float) - np.asarray(ra0, float) + 180.0) % 360.0 - 180.0
    cosd = np.cos(np.radians(0.5 * (np.asarray(dec0, float) + np.asarray(dec1, float))))
    return dra * cosd * 3600.0, (np.asarray(dec1, float) - np.asarray(dec0, float)) * 3600.0


def isolated_mask(cat: EpochCatalog, *, radius_arcsec: float = ISOLATION_ARCSEC) -> np.ndarray:
    """True for components with no OTHER component of the same epoch within ``radius_arcsec``."""
    if len(cat) == 0:
        return np.zeros(0, dtype=bool)
    xyz = _xyz(cat.ra, cat.dec)
    n = cKDTree(xyz).query_ball_point(xyz, r=_chord(radius_arcsec), return_length=True)
    return np.asarray(n) <= 1


def orphan_masks(
    a: EpochCatalog, b: EpochCatalog, *, radius_arcsec: float = STATIC_MATCH_ARCSEC
) -> tuple[np.ndarray, np.ndarray]:
    """Boolean masks of the components in ``a`` and ``b`` with NO counterpart within the radius."""
    ta, tb = cKDTree(_xyz(a.ra, a.dec)), cKDTree(_xyz(b.ra, b.dec))
    r = _chord(radius_arcsec)
    da, _ = tb.query(_xyz(a.ra, a.dec), k=1, distance_upper_bound=r)
    db, _ = ta.query(_xyz(b.ra, b.dec), k=1, distance_upper_bound=r)
    return ~np.isfinite(da), ~np.isfinite(db)


@dataclass
class Pairs:
    """Linked orphan pairs: indices into the two (orphan) catalogues plus the implied motion."""

    i: np.ndarray
    j: np.ndarray
    mu_ra: np.ndarray  # arcsec/yr (RA*cos dec)
    mu_dec: np.ndarray  # arcsec/yr
    dt: np.ndarray  # yr

    def __len__(self) -> int:
        return int(self.i.size)

    @property
    def mu(self) -> np.ndarray:
        return np.hypot(self.mu_ra, self.mu_dec)


def _pair_cov(a: EpochCatalog, ia, b: EpochCatalog, jb, model: ErrorModel) -> np.ndarray:
    """Covariance of the offset between detections ``a[ia]`` and ``b[jb]`` (arcsec^2)."""
    cov = a.cov[ia] + b.cov[jb]
    if model.k_struct or model.q_beam:
        ref = 0.5 * (a.beam[ia] + b.beam[jb])
        cov = cov + model.sys_cov(a.shape[ia], a.beam[ia], ref)
        cov = cov + model.sys_cov(b.shape[jb], b.beam[jb], ref)
    return cov


def link_pairs(
    a: EpochCatalog,
    b: EpochCatalog,
    *,
    mu_min: float = MU_MIN_ARCSEC_YR,
    mu_max: float = MU_MAX_ARCSEC_YR,
    n_sigma_min: float = 3.0,
    model: ErrorModel = NO_MODEL,
) -> Pairs:
    """All (a, b) pairs whose implied proper motion lies in [mu_min, mu_max] arcsec/yr.

    The shift must also be significant: its chi^2 (2 dof) under the pair's full covariance --- both
    measurement ellipses plus the ``model`` systematics --- must reach ``n_sigma_min^2``, so an
    astrometric-noise offset of a static source cannot enter even if it slipped past the orphan
    cut. With isotropic errors and no model this is ``sep >= n_sigma_min * hypot(err_a, err_b)``.
    Per-component epochs make ``dt`` pair-specific.
    """
    empty = Pairs(*(np.zeros(0, dtype=np.int64),) * 2, *(np.zeros(0),) * 3)
    if len(a) == 0 or len(b) == 0:
        return empty
    dt_max = float(np.max(b.t_yr) - np.min(a.t_yr))
    if dt_max <= 0:
        return empty
    tb = cKDTree(_xyz(b.ra, b.dec))
    hits = tb.query_ball_point(_xyz(a.ra, a.dec), r=_chord(mu_max * dt_max + 1.0))
    ii = np.repeat(np.arange(len(a)), [len(h) for h in hits])
    jj = np.fromiter((j for h in hits for j in h), dtype=np.int64, count=ii.size)
    if ii.size == 0:
        return empty
    dx, dy = tangent_offsets_arcsec(a.ra[ii], a.dec[ii], b.ra[jj], b.dec[jj])
    dt = b.t_yr[jj] - a.t_yr[ii]
    sep = np.hypot(dx, dy)
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = np.where(dt > 0, sep / dt, np.inf)
    z2 = chi2_2d(dx, dy, _pair_cov(a, ii, b, jj, model))
    keep = (dt > 0) & (mu >= mu_min) & (mu <= mu_max) & (z2 >= n_sigma_min**2)
    return Pairs(ii[keep], jj[keep], dx[keep] / dt[keep], dy[keep] / dt[keep], dt[keep])


@dataclass
class Triplets:
    """Pairs confirmed by an E3 detection on the predicted straight line."""

    i: np.ndarray
    j: np.ndarray
    k: np.ndarray
    mu_ra: np.ndarray
    mu_dec: np.ndarray
    resid_sigma: np.ndarray  # E3 miss distance in units of the combined positional error

    def __len__(self) -> int:
        return int(self.i.size)

    @property
    def mu(self) -> np.ndarray:
        return np.hypot(self.mu_ra, self.mu_dec)


def collinearity_test(
    a: EpochCatalog,
    b: EpochCatalog,
    c: EpochCatalog,
    pairs: Pairs,
    *,
    tol_sigma: float = E3_TOL_SIGMA,
    mu_max: float = MU_MAX_ARCSEC_YR,
    model: ErrorModel = NO_MODEL,
) -> Triplets:
    """Keep pairs that an E3 component continues on a straight line at the same rate.

    For each pair, E3 components near the E2 position are candidates; each is tested against the
    position the pair's motion predicts at *that component's own* epoch. The tolerance combines
    the E3 error with the prediction error (the pair's two positional errors propagated forward).
    Where several E3 components qualify, the best (smallest residual) is kept.

    Errors are full 2-D covariances: ``S = r^2 S1 + (1 + r)^2 S2 + S3`` with r = dt23/dt12 and
    ``S_e`` = measurement ellipse + ``model`` systematics (beam change relative to the mean beam
    of the three detections). The residual is the Mahalanobis distance, Rayleigh(1) if the model
    is right.
    """
    out: dict[str, list] = {k: [] for k in ("i", "j", "k", "mura", "mudec", "res")}
    if len(pairs) == 0 or len(c) == 0:
        return Triplets(*(np.zeros(0, dtype=np.int64),) * 3, *(np.zeros(0),) * 3)
    dt_max = max(float(np.max(c.t_yr) - np.min(b.t_yr)), 0.0)
    tc = cKDTree(_xyz(c.ra, c.dec))
    near = tc.query_ball_point(
        _xyz(b.ra[pairs.j], b.dec[pairs.j]), r=_chord(mu_max * dt_max + 10.0)
    )
    for n, ks in enumerate(near):
        if not ks:
            continue
        i, j = pairs.i[n], pairs.j[n]
        ks = np.asarray(ks)
        t3 = c.t_yr[ks]
        dt23, dt12 = t3 - b.t_yr[j], pairs.dt[n]
        # Predicted E3 offset from the E2 position, and the observed one.
        px, py = pairs.mu_ra[n] * dt23, pairs.mu_dec[n] * dt23
        ox, oy = tangent_offsets_arcsec(b.ra[j], b.dec[j], c.ra[ks], c.dec[ks])
        # Residual = -d1*r + d2*(1 + r) - d3 with r = dt23/dt12 (d = positional error per
        # epoch). E2's error enters TWICE -- through the fitted rate and as the anchor the
        # prediction starts from -- so it is not independent of the rate error; treating the
        # two as independent underestimated the tolerance and rejected real movers at ~3 sigma.
        r = dt23 / dt12
        s1, s2, s3 = a.cov[i][None, :], b.cov[j][None, :], c.cov[ks]
        if model.k_struct or model.q_beam:
            ref = (a.beam[i][None, :] + b.beam[j][None, :] + c.beam[ks]) / 3.0
            s1 = s1 + model.sys_cov(a.shape[i][None, :], a.beam[i][None, :], ref)
            s2 = s2 + model.sys_cov(b.shape[j][None, :], b.beam[j][None, :], ref)
            s3 = s3 + model.sys_cov(c.shape[ks], c.beam[ks], ref)
        cov = (r**2)[:, None] * s1 + ((1.0 + r) ** 2)[:, None] * s2 + s3
        res = np.sqrt(chi2_2d(ox - px, oy - py, cov))
        ok = (dt23 > 0) & (res <= tol_sigma)
        if ok.any():
            best = int(np.argmin(np.where(ok, res, np.inf)))
            out["i"].append(i)
            out["j"].append(j)
            out["k"].append(int(ks[best]))
            out["mura"].append(pairs.mu_ra[n])
            out["mudec"].append(pairs.mu_dec[n])
            out["res"].append(float(res[best]))
    return Triplets(
        np.asarray(out["i"], dtype=np.int64),
        np.asarray(out["j"], dtype=np.int64),
        np.asarray(out["k"], dtype=np.int64),
        np.asarray(out["mura"], float),
        np.asarray(out["mudec"], float),
        np.asarray(out["res"], float),
    )


def triplet_significance(
    a: EpochCatalog,
    i: int,
    b: EpochCatalog,
    j: int,
    c: EpochCatalog,
    k: int,
    model: ErrorModel = NO_MODEL,
) -> dict:
    """Significances of one detection triple under ``model``, for re-scoring old candidates.

    ``z_12`` and ``z_13``: the E1->E2 and E1->E3 displacements in units of their covariance (a
    static source gives Rayleigh(1); the search links a pair only when ``z_12 >= 3``).
    ``resid``: the E3 straight-line residual exactly as :func:`collinearity_test` computes it.
    """
    ia, jb, kc = np.array([i]), np.array([j]), np.array([k])
    dx12, dy12 = tangent_offsets_arcsec(a.ra[ia], a.dec[ia], b.ra[jb], b.dec[jb])
    dx13, dy13 = tangent_offsets_arcsec(a.ra[ia], a.dec[ia], c.ra[kc], c.dec[kc])
    z12 = float(np.sqrt(chi2_2d(dx12, dy12, _pair_cov(a, ia, b, jb, model)))[0])
    z13 = float(np.sqrt(chi2_2d(dx13, dy13, _pair_cov(a, ia, c, kc, model)))[0])
    dt12 = float(b.t_yr[j] - a.t_yr[i])
    pairs = Pairs(ia, jb, dx12 / dt12, dy12 / dt12, np.array([dt12]))
    t = collinearity_test(a, b, c.subset(kc), pairs, tol_sigma=np.inf, model=model)
    return {
        "z_12": z12,
        "z_13": z13,
        "resid": float(t.resid_sigma[0]) if len(t) else float("nan"),
    }


# ------------------------------------------------------- tracks with parallax (referee 1)


def parallax_factors(ra_deg, dec_deg, t_yr) -> tuple[np.ndarray, np.ndarray]:
    """Parallax factors (P_alpha*, P_delta): the apparent offset in arcsec per arcsec of
    parallax at decimal years ``t_yr``, from the Earth's barycentric position (astropy's
    built-in ephemeris). Standard form, e.g. Green (1985):
    ``P_a = X sin a - Y cos a``, ``P_d = X cos a sin d + Y sin a sin d - Z cos d``."""
    from astropy.coordinates import get_body_barycentric
    from astropy.time import Time

    t = np.atleast_1d(np.asarray(t_yr, float))
    xyz = get_body_barycentric("earth", Time(t, format="decimalyear")).xyz.to_value("au")
    x, y, z = xyz
    a, d = np.radians(ra_deg), np.radians(dec_deg)
    pa = x * np.sin(a) - y * np.cos(a)
    pd = x * np.cos(a) * np.sin(d) + y * np.sin(a) * np.sin(d) - z * np.cos(d)
    return np.asarray(pa, float), np.asarray(pd, float)


def fit_track(
    ra_deg,
    dec_deg,
    t_yr,
    cov: np.ndarray,
    *,
    parallax_arcsec: float | None = 0.0,
    t_ref: float | None = None,
) -> dict:
    """Weighted least-squares straight line (+ parallax) through a multi-epoch radio track.

    Each epoch is weighted by its full 2x2 positional covariance ``cov`` ((n, 3) arcsec^2).
    ``parallax_arcsec``: a fixed parallax (0 = none) or ``None`` to fit it. Returns the proper
    motion (arcsec/yr), its covariance and 1-sigma error on |mu|, the position at ``t_ref``,
    chi^2, degrees of freedom and, for a free fit, the parallax and its error.
    """
    ra, dec, t = (np.atleast_1d(np.asarray(x, float)) for x in (ra_deg, dec_deg, t_yr))
    cov = np.atleast_2d(cov)
    n = ra.size
    t_ref = float(np.mean(t)) if t_ref is None else float(t_ref)
    ra0, dec0 = float(ra[0]), float(dec[0])
    dx, dy = tangent_offsets_arcsec(ra0, dec0, ra, dec)
    pa, pd = parallax_factors(ra0, dec0, t)
    free = parallax_arcsec is None
    plx = 0.0 if parallax_arcsec is None else float(parallax_arcsec)
    npar = 5 if free else 4
    a_mat = np.zeros((2 * n, npar))
    b = np.zeros(2 * n)
    w = np.zeros((2 * n, 2 * n))
    for k in range(n):
        dt = t[k] - t_ref
        a_mat[2 * k, [0, 2]] = [1.0, dt]
        a_mat[2 * k + 1, [1, 3]] = [1.0, dt]
        if free:
            a_mat[2 * k, 4], a_mat[2 * k + 1, 4] = pa[k], pd[k]
            b[2 * k], b[2 * k + 1] = dx[k], dy[k]
        else:
            b[2 * k] = dx[k] - plx * pa[k]
            b[2 * k + 1] = dy[k] - plx * pd[k]
        c = np.array([[cov[k, 0], cov[k, 2]], [cov[k, 2], cov[k, 1]]])
        w[2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = np.linalg.inv(c)
    normal = a_mat.T @ w @ a_mat
    pcov = np.linalg.inv(normal)
    sol = pcov @ a_mat.T @ w @ b
    resid = b - a_mat @ sol
    chi2 = float(resid @ w @ resid)
    mu_ra, mu_dec = float(sol[2]), float(sol[3])
    mu = float(np.hypot(mu_ra, mu_dec))
    grad = np.array([mu_ra, mu_dec]) / max(mu, 1e-12)
    mu_err = float(np.sqrt(grad @ pcov[2:4, 2:4] @ grad))
    out = {
        "n_epochs": int(n),
        "t_ref": t_ref,
        "mu_ra": mu_ra,
        "mu_dec": mu_dec,
        "mu": mu,
        "mu_err": mu_err,
        "mu_cov": pcov[2:4, 2:4].tolist(),
        "x0": float(sol[0]),
        "y0": float(sol[1]),
        "pos_cov_ref": pcov[0:2, 0:2].tolist(),
        "chi2": chi2,
        "dof": int(2 * n - npar),
        "parallax_fixed_arcsec": None if free else plx,
    }
    if free:
        out["parallax_arcsec"] = float(sol[4])
        out["parallax_err_arcsec"] = float(np.sqrt(pcov[4, 4]))
    return out


def predict_track(fit: dict, ra0: float, dec0: float, t: float) -> tuple[float, float, np.ndarray]:
    """Position (deg) and 2x2 covariance (arcsec^2) of a :func:`fit_track` line at epoch ``t``
    (no parallax term: used to place a counterpart at its catalogue epoch)."""
    dt = t - fit["t_ref"]
    j = np.array([[1.0, 0.0, dt, 0.0], [0.0, 1.0, 0.0, dt]])
    full = np.zeros((4, 4))
    full[:2, :2] = np.asarray(fit["pos_cov_ref"])
    full[2:, 2:] = np.asarray(fit["mu_cov"])
    # the position-rate covariance is not stored; t_ref at the mean epoch makes it ~0
    x = fit["x0"] + fit["mu_ra"] * dt
    y = fit["y0"] + fit["mu_dec"] * dt
    dec = dec0 + y / 3600.0
    ra = ra0 + x / 3600.0 / np.cos(np.radians(dec0))
    return float(ra % 360.0), float(dec), j @ full @ j.T


def parallax_residual_factor(ra_deg, dec_deg, t1, t2, t3) -> np.ndarray:
    """|E3 straight-line residual| per arcsec of parallax for detections at t1 < t2 < t3: the
    search's line through E1 and E2 ignores parallax, so a nearby source misses its E3 position
    by ``parallax x`` this factor (arcsec/arcsec)."""
    t1, t2, t3 = (np.asarray(x, float) for x in (t1, t2, t3))
    r = (t3 - t2) / (t2 - t1)
    p1, p2, p3 = (np.column_stack(parallax_factors(ra_deg, dec_deg, t)) for t in (t1, t2, t3))
    res = -p1 * r[:, None] + p2 * (1.0 + r)[:, None] - p3
    return np.hypot(res[:, 0], res[:, 1])


# ---------------------------------------------------------- counterparts (referee 1)

COUNTERPART_RADIUS_ARCSEC = 5.0  # the vet_counterparts search radius
COMOVING_NSIGMA = 3.0
COMOVING_POS_FLOOR_ARCSEC = 0.5  # catalogue astrometry + epoch-propagation floor
COMOVING_PM_FRAC = 0.3  # |mu_cat - mu_radio| < 30% of |mu_radio| also counts as co-moving


def comoving(
    pred_ra: float,
    pred_dec: float,
    pred_cov: np.ndarray,
    mu_radio: tuple[float, float],
    cat_ra: np.ndarray,
    cat_dec: np.ndarray,
    cat_pmra: np.ndarray | None = None,
    cat_pmdec: np.ndarray | None = None,
) -> np.ndarray:
    """Which catalogue sources are CO-MOVING counterparts of a radio track.

    A source counts if (a) its position at its own catalogue epoch lies within
    ``COMOVING_NSIGMA`` of where the radio track predicts the mover at that epoch (prediction
    covariance + a ``COMOVING_POS_FLOOR_ARCSEC`` floor), or (b) its catalogue proper motion
    agrees with the radio rate to ``COMOVING_PM_FRAC`` of |mu|. An unrelated source that merely
    falls inside the 5" search circle does neither, so it no longer makes a mover "not dark".
    """
    cat_ra, cat_dec = np.atleast_1d(cat_ra).astype(float), np.atleast_1d(cat_dec).astype(float)
    dx, dy = tangent_offsets_arcsec(pred_ra, pred_dec, cat_ra, cat_dec)
    c = np.asarray(pred_cov, float)
    f2 = COMOVING_POS_FLOOR_ARCSEC**2
    cov3 = np.array([[c[0, 0] + f2, c[1, 1] + f2, c[0, 1]]])
    pos_ok = chi2_2d(dx, dy, np.repeat(cov3, dx.size, axis=0)) <= COMOVING_NSIGMA**2
    pm_ok = np.zeros(dx.size, dtype=bool)
    if cat_pmra is not None and cat_pmdec is not None:
        pmra = np.atleast_1d(np.asarray(cat_pmra, float))
        pmde = np.atleast_1d(np.asarray(cat_pmdec, float))
        mu = float(np.hypot(*mu_radio))
        diff = np.hypot(pmra - mu_radio[0], pmde - mu_radio[1])
        pm_ok = np.isfinite(diff) & (diff < COMOVING_PM_FRAC * mu)
    return pos_ok | pm_ok


def chance_within(distances_arcsec: list, radii) -> list[float]:
    """P(at least one catalogue source within r) for each r, from nearest-source distances at
    random positions (``inf`` where none was found within the query radius)."""
    d = np.asarray(distances_arcsec, float)
    return [float(np.mean(d <= r)) for r in radii]


def flux_consistent(
    fa: np.ndarray, fb: np.ndarray, fc: np.ndarray, *, max_ratio: float | None = FLUX_RATIO_MAX
) -> np.ndarray:
    """True where max/min of the three peak fluxes is within ``max_ratio`` (all True if None)."""
    f = np.column_stack([fa, fb, fc]).astype(float)
    if max_ratio is None:
        return np.ones(f.shape[0], dtype=bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = f.max(axis=1) / f.min(axis=1)
    return np.isfinite(r) & (r <= max_ratio)


# ----------------------------------------------------------------------------- compactness
#
# Every candidate in runs 3 and 4 other than UV Ceti was a resolved static source. Stars and
# pulsars are unresolved at 2.5", so a mover should be compact in its detections. The metric is
# the catalogue's own deconvolved major axis over the restoring beam's major axis, per detection
# (0 where PyBDSF could not deconvolve, i.e. the fit is no larger than the beam). A faint point
# source does NOT come out at 0 -- its deconvolved size is noisy -- so the threshold must be
# calibrated on injections whose sizes carry that noise (:class:`SizeNoiseModel`).

COMPACT_RULES = ("all", "2of3")


def compactness(cat: EpochCatalog) -> np.ndarray:
    """Deconvolved major FWHM / beam major FWHM per component (0 = unresolved or no beam)."""
    dc = cov_axes(cat.shape)[0]
    bm = cov_axes(cat.beam)[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(bm > 0, dc / bm, 0.0)


def compact_keep(c: np.ndarray, threshold: float | None, rule: str = "all") -> np.ndarray:
    """Which detection triples pass the cut. ``c`` is (n, 3) compactness; ``rule`` is ``"all"``
    (every detection <= threshold) or ``"2of3"`` (at least two). ``threshold=None``: no cut."""
    c = np.atleast_2d(np.asarray(c, float))
    if threshold is None:
        return np.ones(c.shape[0], dtype=bool)
    n_ok = np.sum(c <= threshold, axis=1)
    if rule == "all":
        return n_ok == c.shape[1]
    if rule == "2of3":
        return n_ok >= 2
    raise ValueError(f"unknown compactness rule {rule!r}; use one of {COMPACT_RULES}")


@dataclass
class SizeNoiseModel:
    """Empirical deconvolved sizes of POINT sources, by epoch and peak S/N.

    ``samples`` is (m, 4): ``[epoch_index, snr, dcmaj/bmaj, dcmin/bmaj]`` for detections of a
    reference population known to be point-like at 2.5" (e.g. nearby Gaia stars). :meth:`draw`
    resamples the ratio pair from the same epoch and S/N bin, pooling epochs when that bin holds
    fewer than ``min_count`` detections, and the nearest populated bin when the S/N is outside
    the reference range.
    """

    samples: np.ndarray
    snr_edges: tuple[float, ...] = (0.0, 7.0, 10.0, 15.0, 25.0, 50.0, np.inf)
    min_count: int = 30

    def _bin(self, snr: np.ndarray) -> np.ndarray:
        return np.clip(np.digitize(snr, self.snr_edges) - 1, 0, len(self.snr_edges) - 2)

    def pool(self, epoch: int, snr_bin: int) -> np.ndarray:
        """Row indices into ``samples`` used for one (epoch, S/N bin)."""
        sb = self._bin(self.samples[:, 1])
        own = np.flatnonzero((self.samples[:, 0] == epoch) & (sb == snr_bin))
        if own.size >= self.min_count:
            return own
        pooled = np.flatnonzero(sb == snr_bin)
        if pooled.size >= self.min_count:
            return pooled
        filled = [b for b in range(len(self.snr_edges) - 1) if np.sum(sb == b) >= 1]
        if not filled:
            raise ValueError("empty size-noise reference sample")
        nearest = min(filled, key=lambda b: abs(b - snr_bin))
        return np.flatnonzero(sb == nearest)

    def draw(
        self, snr: np.ndarray, epoch: int, rng: np.random.Generator
    ) -> tuple[np.ndarray, np.ndarray]:
        snr = np.asarray(snr, float)
        rmaj, rmin = np.zeros(snr.size), np.zeros(snr.size)
        bins = self._bin(snr)
        for b in np.unique(bins):
            sel = np.flatnonzero(bins == b)
            rows = self.pool(epoch, int(b))
            pick = rows[rng.integers(0, rows.size, sel.size)]
            rmaj[sel], rmin[sel] = self.samples[pick, 2], self.samples[pick, 3]
        return rmaj, np.minimum(rmin, rmaj)


@dataclass
class SearchResult:
    n_orphans: tuple[int, int, int]
    pairs: Pairs
    triplets: Triplets
    candidates: Triplets  # triplets that also pass the flux cut
    orphans: tuple[EpochCatalog, EpochCatalog, EpochCatalog]


def search(
    e1: EpochCatalog,
    e2: EpochCatalog,
    e3: EpochCatalog,
    *,
    mu_min: float = MU_MIN_ARCSEC_YR,
    mu_max: float = MU_MAX_ARCSEC_YR,
    isolation_arcsec: float | None = ISOLATION_ARCSEC,
    max_flux_ratio: float | None = FLUX_RATIO_MAX,
    model: ErrorModel = NO_MODEL,
    compact_max: float | None = None,
    compact_rule: str = "all",
) -> SearchResult:
    """Steps 1-4: isolated orphans, E1 x E2 linkage, E3 collinearity, (optional) flux cut, and
    (optional) compactness cut on the three detections (:func:`compact_keep`). The cut only
    removes: ``triplets`` is unchanged, ``candidates`` is the subset that passes."""
    o1, o2 = orphan_masks(e1, e2)
    # An E3 orphan must be unmatched in BOTH earlier epochs (a mover's E3 position is new sky).
    o3a, _ = orphan_masks(e3, e1)
    o3b, _ = orphan_masks(e3, e2)
    o3 = o3a & o3b
    if isolation_arcsec is not None:
        o1 &= isolated_mask(e1, radius_arcsec=isolation_arcsec)
        o2 &= isolated_mask(e2, radius_arcsec=isolation_arcsec)
        o3 &= isolated_mask(e3, radius_arcsec=isolation_arcsec)
    a, b, c = e1.subset(o1), e2.subset(o2), e3.subset(o3)
    pairs = link_pairs(a, b, mu_min=mu_min, mu_max=mu_max, model=model)
    trip = collinearity_test(a, b, c, pairs, mu_max=mu_max, model=model)
    fok = flux_consistent(a.flux[trip.i], b.flux[trip.j], c.flux[trip.k], max_ratio=max_flux_ratio)
    if compact_max is not None:
        cvals = np.column_stack(
            [compactness(a)[trip.i], compactness(b)[trip.j], compactness(c)[trip.k]]
        )
        fok &= compact_keep(cvals, compact_max, compact_rule)
    cand = Triplets(
        trip.i[fok],
        trip.j[fok],
        trip.k[fok],
        trip.mu_ra[fok],
        trip.mu_dec[fok],
        trip.resid_sigma[fok],
    )
    return SearchResult((len(a), len(b), len(c)), pairs, trip, cand, (a, b, c))


def _rotate_ra(cat: EpochCatalog, shift_deg: float) -> EpochCatalog:
    from dataclasses import replace

    return replace(cat, ra=(cat.ra + shift_deg) % 360.0)


def scramble_null(
    orphans: tuple[EpochCatalog, EpochCatalog, EpochCatalog],
    *,
    n_reps: int = 20,
    shift_arcmin: tuple[float, float] = (3.0, 30.0),
    seed: int = 0,
    mu_min: float = MU_MIN_ARCSEC_YR,
    mu_max: float = MU_MAX_ARCSEC_YR,
    model: ErrorModel = NO_MODEL,
    compact_max: float | None = None,
    compact_rule: str = "all",
) -> dict:
    """Chance pairs/triplets/candidates per scramble, from RA-rotated E2 and E3 orphans.

    E2 and E3 are rotated by *independent* random shifts, so a real mover can survive neither
    the pairing nor the E3 line; the orphan density (and its Dec dependence) is preserved, and
    every component keeps its own error ellipse, shape and beam (the same ``model`` is applied).
    """
    a, b, c = orphans
    rng = np.random.default_rng(seed)
    k_a, k_b, k_c = compactness(a), compactness(b), compactness(c)
    n_pairs, n_trip, n_cand = [], [], []
    for _ in range(n_reps):
        s2, s3 = rng.uniform(*shift_arcmin, size=2) / 60.0 * rng.choice([-1, 1], size=2)
        bb, cc = _rotate_ra(b, float(s2)), _rotate_ra(c, float(s3))
        p = link_pairs(a, bb, mu_min=mu_min, mu_max=mu_max, model=model)
        t = collinearity_test(a, bb, cc, p, mu_max=mu_max, model=model)
        f = flux_consistent(a.flux[t.i], bb.flux[t.j], cc.flux[t.k])
        f &= compact_keep(
            np.column_stack([k_a[t.i], k_b[t.j], k_c[t.k]]), compact_max, compact_rule
        )
        n_pairs.append(len(p))
        n_trip.append(len(t))
        n_cand.append(int(f.sum()))
    return {
        "n_reps": n_reps,
        "pairs_mean": float(np.mean(n_pairs)),
        "triplets_mean": float(np.mean(n_trip)),
        "candidates_mean": float(np.mean(n_cand)),
        "candidates_per_rep": n_cand,
    }


def _mover_track(ra0, dec0, mu_ra, mu_dec, t0, t):
    """Positions at epochs ``t`` of sources at (ra0, dec0) at ``t0`` moving at (mu_ra, mu_dec)."""
    dt = np.asarray(t, float) - np.asarray(t0, float)
    dec = dec0 + mu_dec * dt * ARCSEC
    ra = (ra0 + mu_ra * dt * ARCSEC / np.cos(np.radians(dec0))) % 360.0
    return ra, dec


INJ_FLUX_SPREAD = (0.8, 1.25)  # per-epoch peak = flux x U(0.8, 1.25): a mildly variable source
PLACE_ISOLATED_ARCSEC = (60.0, 120.0)  # runs 1-3: into empty sky, never near a real component
PLACE_RANDOM_ARCSEC = 900.0  # "realistic": uniform over a 15' disc round a host = random sky


def _injected_errors(
    cat: EpochCatalog,
    near: np.ndarray,
    flux_mjy: np.ndarray,
    size_arcsec: float | tuple[np.ndarray, np.ndarray],
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Measurement covariance etc. for injected detections observed like component ``near``.

    The injected source is a Gaussian of deconvolved FWHM ``size_arcsec`` (a scalar: circular;
    or a ``(dcmaj, dcmin)`` pair of arrays: the MEASURED noisy size of each detection, at a random
    PA) at peak ``flux_mjy``, imaged with the neighbour's beam and noise and given the
    neighbour's floor, so its Condon ellipse is what the real catalogue would have quoted.
    """
    beam = cat.beam[near]
    bmaj, bmin, bpa = cov_axes(beam)
    if isinstance(size_arcsec, tuple):
        dcmaj, dcmin = size_arcsec
        pa = (rng or np.random.default_rng(0)).uniform(-90.0, 90.0, near.size)
    else:
        dcmaj = dcmin = np.full(near.size, float(size_arcsec))
        pa = np.zeros(near.size)
    maj, mn = np.hypot(bmaj, dcmaj), np.hypot(bmin, dcmin)
    snr = flux_mjy / cat.rms[near]
    cov = condon_cov(maj, mn, bpa, bmaj, bmin, snr) + (cat.floor[near] ** 2)[:, None] * np.array(
        [1.0, 1.0, 0.0]
    )
    shape = ellipse_cov(dcmaj, dcmin, pa)
    return cov, shape, beam, cat.rms[near], cat.floor[near]


def _inject(
    *cats: EpochCatalog,
    n: int,
    mu_range: tuple[float, float] = (MU_MIN_ARCSEC_YR, MU_MAX_ARCSEC_YR),
    flux_mjy: float = 3.0,
    seed: int = 0,
    first_ident: int = 10_000_000,
    realistic_frac: float = 0.0,
    size_arcsec: float = 0.0,
    model: ErrorModel = NO_MODEL,
    size_noise: SizeNoiseModel | None = None,
) -> tuple[list[EpochCatalog], np.ndarray, np.ndarray, np.ndarray]:
    """:func:`inject_movers` plus the placement class and each mover's E1 neighbour distance.

    With ``size_noise`` each detection's MEASURED deconvolved size is drawn from the point-source
    reference sample at that detection's S/N and epoch (``size_arcsec`` is then ignored): the
    injection is a point source that the catalogue would nevertheless report with a noisy size.
    """
    e1 = cats[0]
    rng = np.random.default_rng(seed)
    host = rng.integers(0, len(e1), n)
    realistic = rng.random(n) < realistic_frac
    lo, hi = PLACE_ISOLATED_ARCSEC
    off = np.where(
        realistic,
        PLACE_RANDOM_ARCSEC * np.sqrt(rng.random(n)),  # uniform in area: a random sky point
        rng.uniform(lo, hi, n),
    )
    off = off * ARCSEC
    ang = rng.uniform(0, 2 * np.pi, n)
    ra0 = (e1.ra[host] + off * np.cos(ang) / np.cos(np.radians(e1.dec[host]))) % 360.0
    dec0 = np.clip(e1.dec[host] + off * np.sin(ang), -89.9, 89.9)
    mu = np.exp(rng.uniform(np.log(mu_range[0]), np.log(mu_range[1]), n))
    phi = rng.uniform(0, 2 * np.pi, n)
    mu_ra, mu_dec = mu * np.cos(phi), mu * np.sin(phi)
    ids = first_ident + np.arange(n)
    t1 = e1.t_yr[host]
    nn_d, _ = cKDTree(_xyz(e1.ra, e1.dec)).query(_xyz(ra0, dec0), k=1)
    nn_arcsec = np.degrees(2.0 * np.arcsin(nn_d / 2.0)) * 3600.0
    shaped = all(c.has_shape for c in cats)
    nears = [cKDTree(_xyz(c.ra, c.dec)).query(_xyz(ra0, dec0), k=1)[1] for c in cats]
    # The beam every epoch would give this sky position: the mean over epochs is the reference
    # for the beam-change term, exactly as in the collinearity test.
    beam_ref = np.mean([c.beam[k] for c, k in zip(cats, nears, strict=True)], axis=0)
    out = []
    for e, cat in enumerate(cats):
        near = nears[e]
        t = t1 if e == 0 else cat.t_yr[near]
        flux = flux_mjy * rng.uniform(*INJ_FLUX_SPREAD, n)
        if shaped:
            size: float | tuple[np.ndarray, np.ndarray] = size_arcsec
            if size_noise is not None:
                rmaj, rmin = size_noise.draw(flux / cat.rms[near], e, rng)
                bmaj_e = cov_axes(cat.beam[near])[0]
                size = (rmaj * bmaj_e, rmin * bmaj_e)
            cov, shape, beam, rms, floor = _injected_errors(cat, near, flux, size, rng)
            scatter = cov + model.sys_cov(shape, beam, beam_ref)
        else:  # the synthetic isotropic fixture: the epoch's median error, as in runs 1-3
            err0 = float(np.median(cat.pos_err))
            cov = np.column_stack([np.full(n, err0**2), np.full(n, err0**2), np.zeros(n)])
            shape, beam = np.zeros((n, 3)), np.zeros((n, 3))
            rms, floor = np.full(n, np.nan), np.zeros(n)
            scatter = cov
        # Draw the 2-D offset from the full covariance (Cholesky of each 2x2).
        l11 = np.sqrt(scatter[:, 0])
        l21 = scatter[:, 2] / l11
        l22 = np.sqrt(np.maximum(scatter[:, 1] - l21**2, 0.0))
        g1, g2 = rng.normal(0, 1, n), rng.normal(0, 1, n)
        ex, ey = l11 * g1, l21 * g1 + l22 * g2
        ra, dec = _mover_track(ra0, dec0, mu_ra, mu_dec, t1, t)
        ra = (ra + ex * ARCSEC / np.cos(np.radians(dec))) % 360.0
        dec = dec + ey * ARCSEC
        err = np.sqrt(0.5 * (cov[:, 0] + cov[:, 1]))
        inj = EpochCatalog(ra, dec, t, flux, err, ids, cov, shape, beam, rms, floor)
        out.append(EpochCatalog.concat([cat, inj]))
    return out, mu, realistic, nn_arcsec


def inject_movers(*cats: EpochCatalog, n: int, **kw) -> tuple:
    """Plant ``n`` movers into copies of every catalogue; return ``(*cats, mu)``.

    Each mover borrows the sky position and observation time of a random real component of the
    FIRST epoch (so it inherits the real tiling, baselines and Dec distribution) and gets a
    log-uniform rate in ``mu_range`` with a random direction. In later epochs it is observed at
    the time of the nearest real component (that epoch's tile date). Injected rows carry
    ``ident >= first_ident``.

    Placement: by default 60-120" from the host (into empty sky --- runs 1-3, which never paid
    the isolation cost); a fraction ``realistic_frac`` is instead placed uniformly over a 15'
    disc round the host, i.e. at a random sky position, so its distance to the nearest real
    component follows the real neighbour-distance distribution and the isolation cut and
    linkage confusion act on it as on a real mover.

    Scatter: on catalogues with shape information each detection gets the error the catalogue
    would have quoted for it (Condon ellipse at the local beam and noise for the assumed
    ``flux_mjy`` and deconvolved ``size_arcsec``, plus the local floor) plus the ``model``
    systematics, and positions are DRAWN from that full covariance. On the isotropic synthetic
    fixture: the epoch's median ``pos_err``.
    """
    out, mu, _, _ = _inject(*cats, n=n, **kw)
    return (*out, mu)


def epoch_triples(n_epochs: int) -> list[tuple[int, int, int]]:
    """Every time-ordered triple of epoch indices (E1-E2-E3, E1-E2-E4, ... for four epochs)."""
    from itertools import combinations

    return list(combinations(range(n_epochs), 3))


def search_multi(
    cats: list[EpochCatalog],
    triples: list[tuple[int, int, int]] | None = None,
    **kw,
) -> dict[tuple[int, int, int], SearchResult]:
    """:func:`search` on each epoch triple. A source absent from one epoch (a flare star below
    the limit, a tile gap) can still be found through the triples that skip that epoch."""
    triples = epoch_triples(len(cats)) if triples is None else triples
    return {t: search(cats[t[0]], cats[t[1]], cats[t[2]], **kw) for t in triples}


def _recovered_idents(res: SearchResult, first: int) -> np.ndarray:
    a, b, c = res.orphans
    cand = res.candidates
    # A recovery = a candidate whose three components are the SAME injected source.
    same = (
        (a.ident[cand.i] >= first)
        & (a.ident[cand.i] == b.ident[cand.j])
        & (b.ident[cand.j] == c.ident[cand.k])
    )
    return a.ident[cand.i][same] - first


def completeness(
    *cats: EpochCatalog,
    n: int = 2000,
    bins: np.ndarray | None = None,
    seed: int = 0,
    flux_mjy: float = 3.0,
    triples: list[tuple[int, int, int]] | None = None,
    realistic_frac: float = 0.0,
    size_arcsec: float = 0.0,
    model: ErrorModel = NO_MODEL,
    size_noise: SizeNoiseModel | None = None,
    compact_max: float | None = None,
    compact_rule: str = "all",
    mu_range: tuple[float, float] | None = None,
) -> dict:
    """Fraction of injected movers recovered by ANY epoch triple, overall and per log-rate bin.

    Rates are drawn log-uniformly over ``mu_range`` (default: the bin range), so an average over
    equal-width log bins is a log-uniform (scale-free) rate prior.

    With ``realistic_frac > 0`` the result also splits by placement class: ``isolated``
    (60-120" from a real component, the runs 1-3 injections) and ``realistic`` (random sky,
    paying the isolation cut and linkage confusion). The realistic class is the completeness
    for real movers, which do not avoid other sources.
    """
    bins = np.geomspace(MU_MIN_ARCSEC_YR, MU_MAX_ARCSEC_YR, 6) if bins is None else bins
    first = 10_000_000
    inj, mu, realistic, nn = _inject(
        *cats,
        n=n,
        mu_range=(float(bins[0]), float(bins[-1])) if mu_range is None else mu_range,
        seed=seed,
        flux_mjy=flux_mjy,
        first_ident=first,
        realistic_frac=realistic_frac,
        size_arcsec=size_arcsec,
        model=model,
        size_noise=size_noise,
    )
    found = np.zeros(n, dtype=bool)
    per_triple = {}
    kw = {"model": model, "compact_max": compact_max, "compact_rule": compact_rule}
    for t, res in search_multi(list(inj), triples, **kw).items():
        rec = _recovered_idents(res, first)
        found[rec] = True
        per_triple["-".join(f"E{x + 1}" for x in t)] = float(np.unique(rec).size / n)
    idx = np.digitize(mu, bins) - 1

    def _per_bin(sel: np.ndarray) -> list[float]:
        return [
            float(found[sel & (idx == k)].mean()) if np.any(sel & (idx == k)) else float("nan")
            for k in range(len(bins) - 1)
        ]

    out = {
        "n_injected": n,
        "flux_mjy": flux_mjy,
        "size_arcsec": size_arcsec,
        "model": {"k_struct": model.k_struct, "q_beam": model.q_beam},
        "size_noise": size_noise is not None,
        "compact_max": compact_max,
        "compact_rule": compact_rule,
        "overall": float(found.mean()),
        "bin_edges": [float(x) for x in bins],
        "per_bin": _per_bin(np.ones(n, bool)),
        "per_triple": per_triple,
    }
    out["n_per_bin"] = [int(np.sum(idx == k)) for k in range(len(bins) - 1)]
    if realistic_frac > 0:
        out["realistic_frac"] = realistic_frac
        for name, sel in (("isolated", ~realistic), ("realistic", realistic)):
            out[name] = {
                "n": int(sel.sum()),
                "overall": float(found[sel].mean()) if sel.any() else float("nan"),
                "per_bin": _per_bin(sel),
                "frac_neighbour_within_isolation": float(np.mean(nn[sel] < ISOLATION_ARCSEC))
                if sel.any()
                else float("nan"),
            }
    return out


def injection_compactness(
    *cats: EpochCatalog,
    n: int,
    flux_mjy: float,
    size_noise: SizeNoiseModel,
    model: ErrorModel = NO_MODEL,
    seed: int = 0,
    realistic_frac: float = 1.0,
    triples: list[tuple[int, int, int]] | None = None,
) -> dict:
    """Point-source movers with realistic size noise, searched WITHOUT the compactness cut.

    Returns per injection: ``mu``, ``recovered``, ``snr`` (lowest detection S/N of its first
    recovering triple) and ``c`` (n_triples, n, 3) compactness of each recovering triple's
    detections (NaN where that triple did not recover it) -- everything
    :func:`compactness_keep_table` needs to evaluate any threshold and rule afterwards.
    """
    first = 10_000_000
    inj, mu, _, _ = _inject(
        *cats,
        n=n,
        seed=seed,
        flux_mjy=flux_mjy,
        first_ident=first,
        realistic_frac=realistic_frac,
        model=model,
        size_noise=size_noise,
    )
    res = search_multi(list(inj), triples, model=model)
    c = np.full((len(res), n, 3), np.nan)
    snr = np.full(n, np.nan)
    for t_i, r in enumerate(res.values()):
        a, b, cc = r.orphans
        cand = r.candidates
        same = (
            (a.ident[cand.i] >= first)
            & (a.ident[cand.i] == b.ident[cand.j])
            & (b.ident[cand.j] == cc.ident[cand.k])
        )
        idx = a.ident[cand.i][same] - first
        ii, jj, kk = cand.i[same], cand.j[same], cand.k[same]
        c[t_i, idx] = np.column_stack([compactness(a)[ii], compactness(b)[jj], compactness(cc)[kk]])
        s3 = np.min(
            np.column_stack(
                [a.flux[ii] / a.rms[ii], b.flux[jj] / b.rms[jj], cc.flux[kk] / cc.rms[kk]]
            ),
            axis=1,
        )
        fill = np.isnan(snr[idx])
        snr[idx[fill]] = s3[fill]
    return {"mu": mu, "recovered": ~np.all(np.isnan(c[..., 0]), axis=0), "snr": snr, "c": c}


def _kept(c: np.ndarray, threshold: float, rule: str) -> np.ndarray:
    """Injection kept if ANY triple that recovered it passes the cut."""
    out = np.zeros(c.shape[1], dtype=bool)
    for t in range(c.shape[0]):
        ok = ~np.isnan(c[t, :, 0])
        out[ok] |= compact_keep(c[t, ok], threshold, rule)
    return out


def compactness_keep_table(
    inj: dict,
    thresholds,
    *,
    rules: tuple[str, ...] = COMPACT_RULES,
    snr_edges=(0.0, 7.0, 10.0, 15.0, 25.0, np.inf),
    mu_edges=None,
) -> dict:
    """Fraction of RECOVERED injected point movers the cut keeps, per rule and threshold:
    overall, per S/N bin and per rate bin, with the counts behind each fraction."""
    mu_edges = np.geomspace(MU_MIN_ARCSEC_YR, MU_MAX_ARCSEC_YR, 6) if mu_edges is None else mu_edges
    rec = inj["recovered"]
    sb = np.digitize(inj["snr"], snr_edges) - 1
    mb = np.digitize(inj["mu"], mu_edges) - 1
    out: dict = {
        "n_recovered": int(rec.sum()),
        "snr_edges": [float(x) for x in snr_edges],
        "mu_edges": [float(x) for x in mu_edges],
        "snr_counts": [int(np.sum(rec & (sb == k))) for k in range(len(snr_edges) - 1)],
        "mu_counts": [int(np.sum(rec & (mb == k))) for k in range(len(mu_edges) - 1)],
        "rules": {},
    }
    for rule in rules:
        rows = []
        for thr in thresholds:
            kept = _kept(inj["c"], float(thr), rule)

            def frac(sel: np.ndarray, kept: np.ndarray = kept) -> float:
                m = rec & sel
                return float(kept[m].mean()) if m.any() else float("nan")

            rows.append(
                {
                    "threshold": float(thr),
                    "overall": frac(np.ones(rec.size, bool)),
                    "per_snr": [frac(sb == k) for k in range(len(snr_edges) - 1)],
                    "per_mu": [frac(mb == k) for k in range(len(mu_edges) - 1)],
                }
            )
        out["rules"][rule] = rows
    return out


def choose_threshold(
    table: dict, rule: str, *, keep_min: float = 0.95, min_count: int = 200
) -> float | None:
    """The pre-stated criterion: the SMALLEST (most aggressive) threshold on the grid at which the
    cut keeps >= ``keep_min`` of recovered injected point movers overall AND in every S/N bin
    holding >= ``min_count`` recovered injections. None if no threshold on the grid qualifies."""
    ok_bins = [k for k, nk in enumerate(table["snr_counts"]) if nk >= min_count]
    for row in table["rules"][rule]:  # thresholds ascending
        if row["overall"] >= keep_min and all(row["per_snr"][k] >= keep_min for k in ok_bins):
            return float(row["threshold"])
    return None


def calibrate_floors(
    cats: list[EpochCatalog],
    *,
    flux_min_mjy: float = 10.0,
    radius_arcsec: float = STATIC_MATCH_ARCSEC,
) -> dict:
    """Per-epoch astrometric floor, measured from bright static sources, not assumed.

    For every epoch pair, bright components (peak >= ``flux_min_mjy``, where catalogue fitting
    errors are negligible) are matched within ``radius_arcsec`` and the robust per-axis scatter
    of their offsets gives ``s_ab^2 = f_a^2 + f_b^2``. With three or more epochs the per-epoch
    ``f`` follow by least squares. Proper motions of bright extragalactic sources are zero, so
    the scatter is astrometry.
    """
    k = len(cats)
    rows, rhs, pair_sigma = [], [], {}
    for a in range(k):
        for b in range(a + 1, k):
            ca = cats[a].subset(cats[a].flux >= flux_min_mjy)
            cb = cats[b].subset(cats[b].flux >= flux_min_mjy)
            if len(ca) < 10 or len(cb) < 10:
                continue
            d, j = cKDTree(_xyz(cb.ra, cb.dec)).query(
                _xyz(ca.ra, ca.dec), k=1, distance_upper_bound=_chord(radius_arcsec)
            )
            ok = np.isfinite(d)
            if ok.sum() < 10:
                continue
            dx, dy = tangent_offsets_arcsec(ca.ra[ok], ca.dec[ok], cb.ra[j[ok]], cb.dec[j[ok]])
            mad = np.concatenate([dx - np.median(dx), dy - np.median(dy)])
            s = float(1.4826 * np.median(np.abs(mad)))
            pair_sigma[f"E{a + 1}-E{b + 1}"] = s
            row = np.zeros(k)
            row[[a, b]] = 1.0
            rows.append(row)
            rhs.append(s**2)
    if len(rows) < k:
        return {"pair_sigma_arcsec": pair_sigma, "floor_arcsec": None}
    sol, *_ = np.linalg.lstsq(np.asarray(rows), np.asarray(rhs), rcond=None)
    return {
        "pair_sigma_arcsec": pair_sigma,
        "floor_arcsec": [float(np.sqrt(max(x, 0.0))) for x in sol],
    }


DEC_BAND_EDGES = (-40.0, -20.0, 0.0, 30.0, 90.0)  # E1 astrometry degrades south of Dec -20


def calibrate_floors_banded(
    cats: list[EpochCatalog],
    *,
    dec_edges: tuple[float, ...] = DEC_BAND_EDGES,
    flux_min_mjy: float = 10.0,
) -> dict:
    """:func:`calibrate_floors` separately in each declination band.

    VLASS epoch-1 astrometry is documented to be worse in the south (~1" below Dec -20 vs ~0.5"
    north; Memo 22). A single all-sky floor then understates southern errors, and a
    faint southern offset looks significant. Returns ``{"edges": [...], "bands": [{lo, hi,
    floor_arcsec, pair_sigma_arcsec}, ...]}``; a band with too few bright matches gets
    ``floor_arcsec = None`` and callers fall back to the all-sky floor there.
    """
    bands = []
    for lo, hi in zip(dec_edges[:-1], dec_edges[1:], strict=True):
        sub = [c.subset((c.dec >= lo) & (c.dec < hi)) for c in cats]
        cal = calibrate_floors(sub, flux_min_mjy=flux_min_mjy)
        bands.append({"lo": lo, "hi": hi, **cal})
    return {"edges": list(dec_edges), "bands": bands}


def banded_floor(dec: np.ndarray, epoch: int, banded: dict, fallback: float) -> np.ndarray:
    """Per-component floor for epoch index ``epoch`` from :func:`calibrate_floors_banded`."""
    dec = np.asarray(dec, float)
    out = np.full(dec.size, float(fallback))
    for band in banded["bands"]:
        f = band.get("floor_arcsec")
        if f is None:
            continue
        sel = (dec >= band["lo"]) & (dec < band["hi"])
        out[sel] = f[epoch]
    return out


# ------------------------------------------------------- calibration against static sources

CALIB_MATCH_ARCSEC = 5.0  # > the 2.5" orphan radius, so the tail that MAKES orphans is measured
RAYLEIGH_TAIL_3SIGMA = float(np.exp(-4.5))  # P(chi_2 > 3) = 0.0111
RAYLEIGH_MEDIAN = float(np.sqrt(2.0 * np.log(2.0)))  # 1.177


def match_statics(
    a: EpochCatalog,
    b: EpochCatalog,
    *,
    radius_arcsec: float = CALIB_MATCH_ARCSEC,
    isolation_arcsec: float | None = ISOLATION_ARCSEC,
) -> dict:
    """Static sources seen in both epochs: isolated in each (as the search's orphans must be),
    nearest match within ``radius_arcsec``. Returns indices and the (b - a) offsets, arcsec.

    The match radius is deliberately wider than the 2.5" orphan radius: a static source whose
    two positions differ by more than 2.5" is exactly the population that becomes a pair of
    orphans, so a calibration truncated at 2.5" could not see the tail that matters.
    """
    ma = (
        np.ones(len(a), bool)
        if isolation_arcsec is None
        else isolated_mask(a, radius_arcsec=isolation_arcsec)
    )
    mb = (
        np.ones(len(b), bool)
        if isolation_arcsec is None
        else isolated_mask(b, radius_arcsec=isolation_arcsec)
    )
    ia, ib = np.flatnonzero(ma), np.flatnonzero(mb)
    if ia.size == 0 or ib.size == 0:
        z = np.zeros(0, dtype=np.int64)
        return {"i": z, "j": z, "dx": np.zeros(0), "dy": np.zeros(0)}
    d, k = cKDTree(_xyz(b.ra[ib], b.dec[ib])).query(
        _xyz(a.ra[ia], a.dec[ia]), k=1, distance_upper_bound=_chord(radius_arcsec)
    )
    ok = np.isfinite(d)
    i, j = ia[ok], ib[k[ok]]
    dx, dy = tangent_offsets_arcsec(a.ra[i], a.dec[i], b.ra[j], b.dec[j])
    return {"i": i, "j": j, "dx": dx, "dy": dy}


def static_pair_terms(a: EpochCatalog, b: EpochCatalog, match: dict) -> dict:
    """Per matched static: offsets, measurement covariance, and the two systematic templates.

    ``shape`` = sum of the two detections' deconvolved FWHM^2 covariances (k_struct^2 scales it);
    ``dbeam`` = |B_a - Bref| + |B_b - Bref| with Bref the pair's mean beam (= |B_a - B_b|),
    which q_beam^2 scales. Also the binning variables: ``snr`` (lower of the two), ``size``
    (quadrature sum of the two deconvolved major axes, arcsec), ``beam_change`` (|B_a - B_b|
    Frobenius norm over the mean beam area, dimensionless) and ``dec``.
    """
    i, j = match["i"], match["j"]
    ref = 0.5 * (a.beam[i] + b.beam[j])
    dbeam = abs_cov(a.beam[i] - ref) + abs_cov(b.beam[j] - ref)
    diff = a.beam[i] - b.beam[j]
    area = 0.5 * (
        np.prod(cov_axes(a.beam[i])[:2], axis=0) + np.prod(cov_axes(b.beam[j])[:2], axis=0)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        change = np.sqrt(diff[:, 0] ** 2 + diff[:, 1] ** 2 + 2 * diff[:, 2] ** 2) / area
        snr = np.minimum(a.flux[i] / a.rms[i], b.flux[j] / b.rms[j])
    return {
        "dx": match["dx"],
        "dy": match["dy"],
        "cov": a.cov[i] + b.cov[j],
        "shape": a.shape[i] + b.shape[j],
        "dbeam": dbeam,
        "snr": snr,
        "size": np.hypot(cov_axes(a.shape[i])[0], cov_axes(b.shape[j])[0]),
        "beam_change": np.nan_to_num(change),
        "dec": a.dec[i],
    }


def static_z(terms: dict, model: ErrorModel) -> np.ndarray:
    """Normalised offset (Mahalanobis, arcsec/arcsec) of each static pair under ``model``."""
    cov = terms["cov"] + model.k_struct**2 * terms["shape"] + model.q_beam**2 * terms["dbeam"]
    return np.sqrt(chi2_2d(terms["dx"], terms["dy"], cov))


def _subset_terms(terms: dict, sel: np.ndarray) -> dict:
    return {k: x[sel] for k, x in terms.items()}


SIZE_BIN_EDGES = (0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 10.0, np.inf)
SNR_BIN_EDGES = (0.0, 7.0, 15.0, 50.0, np.inf)
BEAM_BIN_EDGES = (0.0, 0.1, 0.2, 0.4, 0.8, np.inf)


def fit_error_model(
    terms: dict,
    *,
    k_grid: np.ndarray | None = None,
    q_grid: np.ndarray | None = None,
    min_count: int = 200,
) -> tuple[ErrorModel, dict]:
    """Grid-fit (k_struct, q_beam) so statics' 3-sigma tail matches Rayleigh in every bin.

    The target is the tail, not the core, because the search's decisions are 3-sigma cuts
    (linkage significance, E3 tolerance). Bins: size x S/N and beam change x size; the loss is
    the sum over bins with >= ``min_count`` pairs of log(f_tail / 0.0111)^2, so every regime
    counts equally rather than the most populous one. Returns the model and the loss surface.
    """
    k_grid = np.linspace(0.0, 0.3, 31) if k_grid is None else np.asarray(k_grid, float)
    q_grid = np.linspace(0.0, 0.2, 21) if q_grid is None else np.asarray(q_grid, float)
    labels = []
    si = np.digitize(terms["size"], SIZE_BIN_EDGES) - 1
    ni = np.digitize(terms["snr"], SNR_BIN_EDGES) - 1
    bi = np.digitize(terms["beam_change"], BEAM_BIN_EDGES) - 1
    labels.append(si * 10 + ni)
    labels.append(100 + si * 10 + bi)
    groups = []
    for lab in labels:
        u, inv, cnt = np.unique(lab, return_inverse=True, return_counts=True)
        groups.append((inv, cnt >= min_count, u.size))
    loss = np.full((k_grid.size, q_grid.size), np.nan)
    for a_, k in enumerate(k_grid):
        for b_, q in enumerate(q_grid):
            tail = static_z(terms, ErrorModel(float(k), float(q))) > 3.0
            tot = 0.0
            for inv, ok, nb in groups:
                f = np.bincount(inv, weights=tail, minlength=nb) / np.bincount(inv, minlength=nb)
                f = np.maximum(f, 1e-4)
                tot += float(np.sum(np.log(f[ok] / RAYLEIGH_TAIL_3SIGMA) ** 2))
            loss[a_, b_] = tot
    ka, qb = np.unravel_index(np.nanargmin(loss), loss.shape)
    best = (float(k_grid[ka]), float(q_grid[qb]))
    return ErrorModel(*best), {
        "k_grid": k_grid.tolist(),
        "q_grid": q_grid.tolist(),
        "loss": loss.tolist(),
        "k_struct": best[0],
        "q_beam": best[1],
    }


def calibration_table(terms: dict, model: ErrorModel, *, by: str, edges) -> list[dict]:
    """Per-bin calibration of ``model`` on statics: median z (Rayleigh 1.177), tail fractions
    above 3 and 5 sigma (Rayleigh 0.0111 and 3.7e-6), and the count."""
    z = static_z(terms, model)
    x = terms[by]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        sel = (x >= lo) & (x < hi)
        if not sel.any():
            continue
        rows.append(
            {
                "lo": float(lo),
                "hi": float(hi),
                "n": int(sel.sum()),
                "median_z": float(np.median(z[sel])),
                "f_gt3": float(np.mean(z[sel] > 3.0)),
                "f_gt5": float(np.mean(z[sel] > 5.0)),
            }
        )
    return rows


def surface_density_limit(n_found: int, area_deg2: float, completeness_frac: float) -> float:
    """95% upper limit on movers per deg^2 at zero survivors (or the Poisson-mean estimate).

    Divides by the completeness-weighted area: a search that could only have seen a fraction of
    the movers in a region constrains that fraction, not the whole area.
    """
    eff = area_deg2 * completeness_frac
    if eff <= 0:
        return float("inf")
    return (POISSON95_ZERO if n_found == 0 else float(n_found)) / eff


def synthetic_epochs(
    *,
    n_static: int = 20_000,
    n_movers: int = 25,
    n_variable: int = 3_000,
    area_side_deg: float = 10.0,
    pos_err: tuple[float, ...] = (0.5, 0.3, 0.15),
    epochs_t: tuple[float, ...] = (2018.5, 2021.5, 2024.0),
    seed: int = 0,
) -> tuple:
    """Offline fixture: statics, planted movers, and variable-source orphans in three epochs.

    Positional errors default to the real per-epoch floors (E1 ~0.5", E2 ~0.3", E3 ~0.15";
    VLASS Memo 22). Variable sources appear in a random subset of epochs at random positions,
    which is what produces chance orphan pairs in real data. Returns ``(*epochs, mu)`` -- one
    catalogue per entry of ``epochs_t`` plus the planted movers' true rates (arcsec/yr); movers
    carry ``ident`` 0..n_movers-1.
    """
    if len(pos_err) != len(epochs_t):
        raise ValueError("pos_err and epochs_t must have one entry per epoch")
    rng = np.random.default_rng(seed)
    half = area_side_deg / 2

    def _sky(n):
        return rng.uniform(180 - half, 180 + half, n), rng.uniform(-half, half, n)

    sra, sdec = _sky(n_static)
    sflux = rng.lognormal(np.log(3.0), 0.8, n_static)
    mra, mdec = _sky(n_movers)
    mu = np.exp(rng.uniform(np.log(0.5), np.log(4.5), n_movers))
    phi = rng.uniform(0, 2 * np.pi, n_movers)
    mflux = rng.uniform(2.0, 6.0, n_movers)
    vra, vdec = _sky(n_variable)
    vflux = rng.lognormal(np.log(1.2), 0.3, n_variable)
    seen = rng.random((n_variable, len(epochs_t))) < 0.5

    cats = []
    for e, t in enumerate(epochs_t):
        err = pos_err[e]
        tt = t + rng.uniform(-0.3, 0.3, n_static + n_movers + n_variable)
        ra_m, dec_m = _mover_track(
            mra,
            mdec,
            mu * np.cos(phi),
            mu * np.sin(phi),
            epochs_t[0],
            tt[n_static : n_static + n_movers],
        )
        ra = np.concatenate([sra, ra_m, vra])
        dec = np.concatenate([sdec, dec_m, vdec])
        flux = np.concatenate([sflux, mflux, vflux])
        ident = np.concatenate(
            [np.full(n_static, -1), np.arange(n_movers), np.full(n_variable, -2)]
        ).astype(np.int64)
        keep = np.concatenate([np.ones(n_static + n_movers, bool), seen[:, e]])
        noise_ra = rng.normal(0, err, ra.size) * ARCSEC / np.cos(np.radians(dec))
        noise_dec = rng.normal(0, err, ra.size) * ARCSEC
        cat = EpochCatalog(
            (ra + noise_ra) % 360.0, dec + noise_dec, tt, flux, np.full(ra.size, err), ident
        )
        cats.append(cat.subset(keep))
    return (*cats, mu)


def synthetic_statics(
    *,
    n: int = 20_000,
    n_epochs: int = 2,
    k_struct: float = 0.1,
    q_beam: float = 0.03,
    floor_arcsec: float = 0.1,
    resolved_frac: float = 0.6,
    area_side_deg: float = 10.0,
    seed: int = 0,
) -> list[EpochCatalog]:
    """Offline calibration fixture: static sources with shapes, beams and S/N, whose positions
    scatter by the Condon ellipse + floor PLUS a structure term ``k_struct^2 * shape`` and a
    beam-change term ``q_beam^2 * |B_e - mean B|`` --- i.e. more than the catalogue would quote
    for resolved sources. :func:`fit_error_model` must recover (k_struct, q_beam) from it, and a
    measurement-only error model must visibly fail on the resolved ones.
    """
    rng = np.random.default_rng(seed)
    half = area_side_deg / 2
    ra = rng.uniform(180 - half, 180 + half, n)
    dec = rng.uniform(-half, half, n)
    theta = np.where(rng.random(n) < resolved_frac, rng.lognormal(np.log(2.0), 0.6, n), 0.0)
    dcmin = theta * rng.uniform(0.3, 1.0, n)
    dcpa = rng.uniform(0, 180, n)
    shape = ellipse_cov(theta, dcmin, dcpa)
    snr = np.exp(rng.uniform(np.log(5.0), np.log(200.0), n))
    rms = np.full(n, 0.14)
    beams = []
    for _ in range(n_epochs):
        bmaj = rng.uniform(2.4, 4.8, n)
        bmin = rng.uniform(1.7, 2.4, n)
        beams.append((bmaj, bmin, rng.uniform(-90, 90, n)))
    bcov = [ellipse_cov(*b) for b in beams]
    bref = np.mean(bcov, axis=0)
    truth = ErrorModel(k_struct, q_beam)
    cats = []
    for e, (bmaj, bmin, bpa) in enumerate(beams):
        maj, mn = np.hypot(bmaj, theta), np.hypot(bmin, dcmin)
        pa = np.where(theta > 0, dcpa, bpa)
        cov = condon_cov(maj, mn, pa, bmaj, bmin, snr)
        cov[:, :2] += floor_arcsec**2
        tot = cov + truth.sys_cov(shape, bcov[e], bref)
        l11 = np.sqrt(tot[:, 0])
        l21 = tot[:, 2] / l11
        l22 = np.sqrt(np.maximum(tot[:, 1] - l21**2, 0.0))
        g1, g2 = rng.normal(0, 1, n), rng.normal(0, 1, n)
        ex, ey = l11 * g1, l21 * g1 + l22 * g2
        cats.append(
            EpochCatalog(
                (ra + ex * ARCSEC / np.cos(np.radians(dec))) % 360.0,
                dec + ey * ARCSEC,
                np.full(n, 2018.5 + 3.0 * e),
                snr * rms,
                np.sqrt(0.5 * (cov[:, 0] + cov[:, 1])),
                np.arange(n),
                cov,
                shape,
                bcov[e],
                rms,
                np.full(n, floor_arcsec),
            )
        )
    return cats


# ------------------------------------------------------------------------------------------ run


# ----------------------------------------------------------------------------- paper


def _fmt_int(x) -> str:
    return f"{int(round(float(x))):,}".replace(",", "{,}")


def _fmt_sci(x: float, digits: int = 2) -> str:
    mant, exp = f"{x:.{digits - 1}e}".split("e")
    return rf"{mant}\times10^{{{int(exp)}}}"


def _syn_macro_values(m: dict) -> dict[str, str]:
    """Offline-fixture numbers for the vpmSyn* namespace (empty for a real metrics dict)."""
    if m.get("is_real") or "syn_n_movers" not in m:
        return {}
    out = {
        "NMovers": str(m["syn_n_movers"]),
        "NRecovered": str(m["syn_n_recovered"]),
        "NFalse": str(m["syn_n_false"]),
        "Completeness": f"{m['syn_completeness']:.2f}",
        "NullCand": f"{m['syn_null_candidates_mean']:.2f}",
    }
    if m.get("syn_missed_mu") is not None:
        out |= {
            "NMissed": str(len(m["syn_missed_mu"])),
            "MissedMuMax": f"{max(m['syn_missed_mu'], default=0.0):.1f}",
            "NMissedNotIsolated": str(m["syn_n_missed_not_isolated"]),
        }
    if m.get("syn_calib_k_struct_fit") is not None:
        out |= {
            "KTrue": f"{m['syn_calib_k_struct_true']:.2f}",
            "KFit": f"{m['syn_calib_k_struct_fit']:.2f}",
            "TailOld": f"{m['syn_calib_resolved_f_gt3_measurement_only']:.2f}",
            "TailFit": f"{m['syn_calib_resolved_f_gt3_fitted']:.3f}",
        }
    return out


def _real_macro_values(m: dict, vet: dict | None) -> dict[str, str]:
    """Every data-derived number the paper quotes, from the committed real evidence."""
    if not m.get("is_real") or "real_run5_compactness" not in m:
        return {}
    r5 = m["real_run5_compactness"]
    cal = m["real_calibration"]
    val = cal["validation"]
    ep = m["real_epochs"]
    uv = m["real_uvcet"]["hits"][0]
    floors = [f for b in m["real_floors_banded"]["bands"] for f in (b.get("floor_arcsec") or [])]
    dec_rows = {(r["lo"], r["hi"]): r for r in val["fitted_model"]["dec"]}
    south = dec_rows[(-90.0, -35.0)]
    tab = r5["calibration"]["table"]
    rule_rows = {r["threshold"]: r for r in tab["rules"][r5["rule"]]}
    kept = rule_rows[r5["threshold"]]
    slow_bin = 1  # 0.53-0.92 arcsec/yr
    uv5 = next(c for c in r5["run4_candidates"] if c["triple"] == uv["triple"] and c["passes"])
    comp = r5["completeness"]
    lim = r5["limits"]
    null = max(v["candidates_mean"] for v in r5["null"].values())
    edges = comp["3mJy_cut"]["bin_edges"]
    n_static_pairs = cal["n_pairs_fit"] + cal["n_pairs_validate"]
    stars = r5["size_noise_reference"]
    out = {
        "NEOne": _fmt_int(ep["E1"]["n_clean"]),
        "NETwo": _fmt_int(ep["E2"]["n_clean"]),
        "NEThree": _fmt_int(ep["E3"]["n_clean"]),
        "NEFour": _fmt_int(ep["E4"]["n_clean"]),
        "TStart": f"{ep['E1']['t_min']:.1f}",
        "TEnd": f"{ep['E4']['t_max']:.1f}",
        "AreaDeg": _fmt_int(r5["area_deg2"]),
        "FloorMin": f"{min(floors):.2f}",
        "FloorMax": f"{max(floors):.2f}",
        "NStaticPairs": f"{n_static_pairs / 1e6:.1f}",
        "KStruct": f"{cal['model']['k_struct']:.2f}",
        "TailLegacy": f"{val['legacy_runs1to3']['all']['f_gt3']:.3f}",
        "TailMeas": f"{val['measurement_only']['all']['f_gt3']:.3f}",
        "TailFit": f"{val['fitted_model']['all']['f_gt3']:.4f}",
        "TailFitFive": f"{val['fitted_model']['all']['f_gt5']:.4f}",
        "MedZFit": f"{val['fitted_model']['all']['median_z']:.2f}",
        "TailSouth": f"{south['f_gt3']:.3f}",
        "TailSouthRatio": f"{south['f_gt3'] / RAYLEIGH_TAIL_3SIGMA:.1f}",
        "UVMu": f"{uv['mu']:.2f}",
        "UVGaiaMu": f"{uv['gaia_mu']:.2f}",
        "UVResid": f"{uv['resid_sigma']:.2f}",
        "UVCompA": f"{uv5['compactness'][0]:.2f}",
        "UVCompB": f"{uv5['compactness'][1]:.2f}",
        "UVCompC": f"{uv5['compactness'][2]:.2f}",
        "NCandBefore": str(len(r5["run4_candidates"])),
        "NCandAfter": str(r5["n_candidates"]),
        "NStaticVetted": str(len(r5["run4_candidates"]) - r5["n_candidates"]),
        "NullCand": f"{null:.2f}",
        "NStars": str(sum(stars["n_per_epoch"])),
        "NStarChance": str(stars["n_chance_expected"]),
        "CutThr": f"{r5['threshold']:.1f}",
        "CutThrAll": f"{r5['calibration']['choice']['all']['threshold']:.1f}",
        "NInjCal": _fmt_int(r5["calibration"]["n_injected"]),
        "NRecCal": _fmt_int(tab["n_recovered"]),
        "KeepAll": f"{kept['overall']:.3f}",
        "KeepLowSNR": f"{kept['per_snr'][0]:.3f}",
        "KeepSlow": f"{kept['per_mu'][slow_bin]:.2f}",
        "StaticRej": f"{100 * r5['calibration']['choice'][r5['rule']]['static_rejected_fraction']:.0f}",
        "StaticRejAll": f"{100 * r5['calibration']['choice']['all']['static_rejected_fraction']:.1f}",
        "NInj": _fmt_int(comp["3mJy_cut"]["n_injected"]),
        "CompThreeCut": f"{lim['3mJy_cut']['completeness_mean_0p92_5']:.3f}",
        "CompThreeNoCut": f"{lim['3mJy_nocut']['completeness_mean_0p92_5']:.3f}",
        "CompOneFiveCut": f"{lim['1.5mJy_cut']['completeness_mean_0p92_5']:.3f}",
        "CompOneFiveNoCut": f"{lim['1.5mJy_nocut']['completeness_mean_0p92_5']:.3f}",
        "CompSlowThreeCut": f"{comp['3mJy_cut']['per_bin'][slow_bin]:.2f}",
        "CompFloorThreeCut": f"{comp['3mJy_cut']['per_bin'][0]:.2f}",
        "RateLo": f"{edges[1]:.2f}",
        "RateMid": f"{edges[2]:.2f}",
        "RateMin": f"{edges[0]:.1f}",
        "RateMax": f"{edges[-1]:.0f}",
        "LimThree": _fmt_sci(lim["3mJy_cut"]["limit_per_deg2_95"]),
        "LimOneFive": _fmt_sci(lim["1.5mJy_cut"]["limit_per_deg2_95"]),
        "AllSky": f"{lim['3mJy_cut']['limit_per_deg2_95'] * 4 * np.pi * (180 / np.pi) ** 2:.1f}",
    }
    # How close the removed candidates came: the threshold at which each would have passed the
    # chosen rule (for 2 of 3, its second-smallest compactness), smallest two first.
    need = sorted(
        (
            (
                sorted(c["compactness"])[len(c["compactness"]) - 2]
                if r5["rule"] == "2of3"
                else max(c["compactness"])
            ),
            c["key"],
        )
        for c in r5["run4_candidates"]
        if not c["passes"]
    )
    for tag, (thr, key) in zip(("A", "B"), need[:2], strict=False):
        ra_part, dec_part = key[1:].split("_")
        out[f"Margin{tag}Name"] = f"J{float(ra_part):.1f}${dec_part[0]}${abs(float(dec_part)):.1f}"
        out[f"Margin{tag}Thr"] = f"{thr:.2f}"
    out["FluxThree"] = f"{comp['3mJy_cut']['flux_mjy']:g}"
    out["FluxOneFive"] = f"{comp['1.5mJy_cut']['flux_mjy']:g}"
    if vet:
        out["NVetted"] = str(vet.get("n_candidates"))
        out["NNewMovers"] = str(vet.get("n_new_movers"))
        cands = (vet.get("candidates") or {}).values()
        out["NSouthVetted"] = str(
            sum(1 for c in cands if c["dec"] < -35.0 and c["verdict"].startswith("static"))
        )
    return out


REAL_MACRO_NAMES = (
    "NEOne",
    "NETwo",
    "NEThree",
    "NEFour",
    "TStart",
    "TEnd",
    "AreaDeg",
    "FloorMin",
    "FloorMax",
    "NStaticPairs",
    "KStruct",
    "TailLegacy",
    "TailMeas",
    "TailFit",
    "TailFitFive",
    "MedZFit",
    "TailSouth",
    "TailSouthRatio",
    "UVMu",
    "UVGaiaMu",
    "UVResid",
    "UVCompA",
    "UVCompB",
    "UVCompC",
    "NCandBefore",
    "NCandAfter",
    "NStaticVetted",
    "NullCand",
    "NStars",
    "NStarChance",
    "CutThr",
    "CutThrAll",
    "NInjCal",
    "NRecCal",
    "KeepAll",
    "KeepLowSNR",
    "KeepSlow",
    "StaticRej",
    "StaticRejAll",
    "NInj",
    "CompThreeCut",
    "CompThreeNoCut",
    "CompOneFiveCut",
    "CompOneFiveNoCut",
    "CompSlowThreeCut",
    "CompFloorThreeCut",
    "RateLo",
    "RateMid",
    "RateMin",
    "RateMax",
    "LimThree",
    "LimOneFive",
    "AllSky",
    "NVetted",
    "NNewMovers",
    "NSouthVetted",
    "MarginAName",
    "MarginAThr",
    "MarginBName",
    "MarginBThr",
    "FluxThree",
    "FluxOneFive",
)
SYN_MACRO_NAMES = (
    "NMovers",
    "NRecovered",
    "NFalse",
    "Completeness",
    "NullCand",
    "NMissed",
    "MissedMuMax",
    "NMissedNotIsolated",
    "KTrue",
    "KFit",
    "TailOld",
    "TailFit",
)


def _ref1_macro_values(m: dict, ref1: dict | None) -> dict[str, str]:
    """Referee-round-1 numbers (results/vlasspm_referee1.json): the UV Ceti track fit, the
    chance-coincidence rates, fine-bin completeness, tail bound, parallax floor, worst-bin
    limits, and the flux and sky domains."""
    if not ref1 or not m.get("is_real"):
        return {}
    uv = ref1["uvcet"]
    g = uv["gaia"]
    fx, fn, ff = (uv["fits"][k] for k in ("fixed_gaia_parallax", "no_parallax", "free_parallax"))
    sysmu = g["system"]["mu_arcsec_yr"]
    ch = ref1["chance_coincidence"]
    radii = ch["radii_arcsec"]
    co_r = COMOVING_NSIGMA * COMOVING_POS_FLOOR_ARCSEC
    fine = ref1["fine_completeness_3mJy_cut"]
    fe = fine["bin_edges"]
    tb = ref1["tail_completeness_bound"]
    lim = ref1["limits"]
    per = m["real_per_triple"]
    area4 = max(v["area_deg2"] for k, v in per.items() if "E4" in k)
    above = [p for p, lo in zip(fine["per_bin"], fe[:-1], strict=True) if lo >= 1.1 - 1e-9]
    return {
        "UVFitMu": f"{fx['mu']:.2f}",
        "UVFitErr": f"{fx['mu_err']:.2f}",
        "UVFitChi": f"{fx['chi2']:.1f}",
        "UVNoPlxMu": f"{fn['mu']:.2f}",
        "UVFreePlx": f"{1000 * ff['parallax_arcsec']:.0f}",
        "UVFreePlxErr": f"{1000 * ff['parallax_err_arcsec']:.0f}",
        "UVGaiaPlx": f"{1000 * uv['parallax_used_arcsec']:.0f}",
        "SysMu": f"{sysmu:.2f}",
        "SysSigma": f"{abs(fx['mu'] - sysmu) / fx['mu_err']:.1f}",
        "BLMu": f"{g['components']['BL Cet']['mu_arcsec_yr']:.2f}",
        "RuweMin": f"{min(c['ruwe'] for c in g['components'].values()):.0f}",
        "RuweMax": f"{max(c['ruwe'] for c in g['components'].values()):.0f}",
        "UVPlxResid": f"{uv['e3_parallax_residual_arcsec']:.1f}",
        "CompSep": f"{np.hypot(*tangent_offsets_arcsec(*(g['components'][n][k] for n in ('UV Cet', 'BL Cet') for k in ('ra', 'dec')))):.1f}",
        "ChanceFive": f"{ch['p_either'][radii.index(COUNTERPART_RADIUS_ARCSEC)]:.2f}",
        "ChanceComoving": f"{ch['p_either'][radii.index(co_r)]:.2f}",
        "ComovingRadius": f"{co_r:.1f}",
        "NChance": str(ch["n_positions"]),
        "FineEdgeLo": f"{fe[1]:.2f}",
        "FineEdgeMid": f"{fe[2]:.1f}",
        "FineEdgeHi": f"{fe[3]:.1f}",
        "FineCompA": f"{fine['per_bin'][1]:.2f}",
        "FineCompB": f"{fine['per_bin'][2]:.2f}",
        "FineCompAboveMin": f"{min(above):.2f}",
        "LimWorstBin": _fmt_sci(lim["worst_bin_limit_per_deg2_95"]),
        "LimEdge": _fmt_sci(lim["worst_fine_bin_limit_per_deg2_95"]),
        "TailBoundThree": f"{100 * tb['3mJy']['max_excess']:.1f}",
        "TailBoundOneFive": f"{100 * tb['1.5mJy']['max_excess']:.1f}",
        "PlxDistMin": f"{ref1['parallax_floor']['distance_min_pc']:.0f}",
        "AreaFour": _fmt_int(area4),
        "FluxLoThree": f"{3.0 * INJ_FLUX_SPREAD[0]:.1f}",
        "FluxHiThree": f"{3.0 * INJ_FLUX_SPREAD[1]:.2f}",
        "FluxLoOneFive": f"{1.5 * INJ_FLUX_SPREAD[0]:.1f}",
        "FluxHiOneFive": f"{1.5 * INJ_FLUX_SPREAD[1]:.1f}",
    }


REF1_MACRO_NAMES = (
    "UVFitMu",
    "UVFitErr",
    "UVFitChi",
    "UVNoPlxMu",
    "UVFreePlx",
    "UVFreePlxErr",
    "UVGaiaPlx",
    "SysMu",
    "SysSigma",
    "BLMu",
    "RuweMin",
    "RuweMax",
    "UVPlxResid",
    "CompSep",
    "ChanceFive",
    "ChanceComoving",
    "ComovingRadius",
    "NChance",
    "FineEdgeLo",
    "FineEdgeMid",
    "FineEdgeHi",
    "FineCompA",
    "FineCompB",
    "FineCompAboveMin",
    "LimWorstBin",
    "LimEdge",
    "TailBoundThree",
    "TailBoundOneFive",
    "PlxDistMin",
    "AreaFour",
    "FluxLoThree",
    "FluxHiThree",
    "FluxLoOneFive",
    "FluxHiOneFive",
)


def _write_macros(m: dict, path, vet: dict | None = None, ref1: dict | None = None) -> None:
    """Both namespaces always emitted; the inactive one as placeholders, merged by
    :func:`report.preserve_live_macros` so neither leg can blank or overwrite the other."""
    from .report import MACRO_PLACEHOLDER, preserve_live_macros

    real = _real_macro_values(m, vet) | _ref1_macro_values(m, ref1)
    syn = _syn_macro_values(m)
    names_real, names_syn = REAL_MACRO_NAMES + REF1_MACRO_NAMES, SYN_MACRO_NAMES
    lines = [
        "% Auto-generated by jansky_research.vlasspm._write_macros -- do not edit.",
        "% vpmReal* come from results/vlasspm_metrics.json + vlasspm_vetting.json (real leg);",
        "% vpmSyn* from the offline fixture. The inactive namespace is written as placeholders",
        "% and preserve_live_macros keeps the other leg's live values.",
        rf"\newcommand{{\vpmSource}}{{{m['source']}}}",
    ]
    lines += [rf"\newcommand{{\vpmSyn{k}}}{{{syn.get(k, MACRO_PLACEHOLDER)}}}" for k in names_syn]
    lines += [
        rf"\newcommand{{\vpmReal{k}}}{{{real.get(k, MACRO_PLACEHOLDER)}}}" for k in names_real
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    text = preserve_live_macros("\n".join(lines) + "\n", p)
    p.write_text(text)


def _paper_figure(m: dict, path, ref1: dict | None = None) -> Path:
    """Completeness per rate bin (steps), 1.5 and 3 mJy, with and without the compactness cut;
    the fine bins near the 0.92"/yr edge (3 mJy, cut) as points; UV Ceti's rate marked."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker

    r5 = m["real_run5_compactness"]
    comp = r5["completeness"]
    edges = np.asarray(comp["3mJy_cut"]["bin_edges"])
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    styles = {
        "3mJy_nocut": ("3 mJy, no cut", "C0", "--"),
        "3mJy_cut": ("3 mJy, cut", "C0", "-"),
        "1.5mJy_nocut": ("1.5 mJy, no cut", "C1", "--"),
        "1.5mJy_cut": ("1.5 mJy, cut", "C1", "-"),
    }
    for key, (lab, col, ls) in styles.items():
        y = np.asarray(comp[key]["per_bin"])
        ax.stairs(y, edges, color=col, ls=ls, lw=1.1, label=lab, baseline=None)
    if ref1:
        fine = ref1["fine_completeness_3mJy_cut"]
        fe = np.asarray(fine["bin_edges"])
        fc = np.sqrt(fe[:-1] * fe[1:])
        ax.errorbar(
            fc,
            fine["per_bin"],
            xerr=[fc - fe[:-1], fe[1:] - fc],
            fmt="o",
            color="k",
            ms=2.5,
            lw=0.7,
            label="3 mJy, cut, fine bins",
        )
    mu_uv = m["real_uvcet"]["hits"][0]["mu"]
    if ref1:
        mu_uv = ref1["uvcet"]["fits"]["fixed_gaia_parallax"]["mu"]
    ax.axvline(mu_uv, color="k", lw=0.8, ls=":")
    ax.text(mu_uv * 0.96, 0.05, "UV Cet", rotation=90, ha="right", va="bottom", fontsize=7)
    ax.set_xscale("log")
    ax.set_xlim(edges[0], edges[-1])
    ticks = [0.3, 0.5, 1, 2, 5]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel(r"proper motion (arcsec yr$^{-1}$)")
    ax.set_ylabel("completeness")
    ax.legend(fontsize=6, loc="upper left", frameon=False)
    fig.tight_layout()
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p)
    plt.close(fig)
    return p


def write_real_paper(out: str | Path = ".") -> dict:
    """Paper macros + figure from the committed real evidence under ``out/results``."""
    import json

    op = Path(out)
    m = json.loads((op / "results" / "vlasspm_metrics.json").read_text())
    vp = op / "results" / "vlasspm_vetting.json"
    vet = json.loads(vp.read_text()) if vp.exists() else None
    rp = op / "results" / "vlasspm_referee1.json"
    ref1 = json.loads(rp.read_text()) if rp.exists() else None
    if not m.get("is_real"):
        raise ValueError("results/vlasspm_metrics.json is not real evidence; refusing to build")
    _write_macros(m, op / "papers" / "vlasspm" / "generated" / "macros.tex", vet, ref1)
    _paper_figure(m, op / "papers" / "vlasspm" / "figures" / "vlasspm_completeness.pdf", ref1)
    return m


def run(out: str = ".", *, offline: bool = True, n_null: int = 20, n_inject: int = 2000) -> dict:
    """Offline: the synthetic recover-a-known (planted movers found; null predicts chance count).

    The real leg is ``scripts/vlasspm_real.py`` (full-sky catalogues, checkpointed, run detached);
    it ends by calling :func:`write_real_paper`, which ``run(offline=False)`` also calls to
    regenerate the paper's macros and figure from the committed results. Both legs write
    ``papers/vlasspm/generated/macros.tex`` under ``out``, merged so neither blanks the other.
    """
    if not offline:  # pragma: no cover - reads the committed real evidence
        return write_real_paper(out)
    e1, e2, e3, mu_true = synthetic_epochs()
    res = search(e1, e2, e3)
    a, b, c = res.orphans
    cand = res.candidates
    same = (
        (a.ident[cand.i] >= 0)
        & (a.ident[cand.i] == b.ident[cand.j])
        & (b.ident[cand.j] == c.ident[cand.k])
    )
    n_movers = int(mu_true.size)
    found_ids = set(np.unique(a.ident[cand.i][same]).tolist())
    missed = np.array([k not in found_ids for k in range(n_movers)])
    iso_all = np.ones(n_movers, dtype=bool)
    for e in (e1, e2, e3):
        iso_ids = set(e.ident[isolated_mask(e)].tolist())
        iso_all &= np.array([k in iso_ids for k in range(n_movers)])
    null = scramble_null(res.orphans, n_reps=n_null)
    comp = completeness(e1, e2, e3, n=n_inject)
    sa, sb = synthetic_statics()
    terms = static_pair_terms(sa, sb, match_statics(sa, sb))
    fitted, _ = fit_error_model(terms)
    old_tail = calibration_table(terms, NO_MODEL, by="size", edges=(3.0, np.inf))
    new_tail = calibration_table(terms, fitted, by="size", edges=(3.0, np.inf))
    metrics = {
        "source": "synthetic: planted movers + variable-source orphans (offline fixture)",
        "is_real": False,
        "syn_n_movers": n_movers,
        "syn_n_recovered": int(np.unique(a.ident[cand.i][same]).size),
        "syn_n_false": int((~same).sum()),
        "syn_missed_mu": [float(x) for x in np.sort(mu_true[missed])],
        "syn_n_missed_not_isolated": int(np.sum(missed & ~iso_all)),
        "syn_n_missed_isolated_slow": int(np.sum(missed & iso_all & (mu_true < 1.0))),
        "syn_n_missed_isolated_fast": int(np.sum(missed & iso_all & (mu_true >= 1.0))),
        "syn_n_orphans": list(res.n_orphans),
        "syn_n_pairs": len(res.pairs),
        "syn_n_triplets": len(res.triplets),
        "syn_null_candidates_mean": null["candidates_mean"],
        "syn_null_pairs_mean": null["pairs_mean"],
        "syn_completeness": comp["overall"],
        "syn_completeness_per_bin": comp["per_bin"],
        "syn_bin_edges": comp["bin_edges"],
        "syn_calib_k_struct_true": 0.1,
        "syn_calib_q_beam_true": 0.03,
        "syn_calib_k_struct_fit": fitted.k_struct,
        "syn_calib_q_beam_fit": fitted.q_beam,
        "syn_calib_resolved_f_gt3_measurement_only": old_tail[0]["f_gt3"] if old_tail else None,
        "syn_calib_resolved_f_gt3_fitted": new_tail[0]["f_gt3"] if new_tail else None,
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "vlasspm_metrics.json")
    _write_macros(metrics, op / "papers" / "vlasspm" / "generated" / "macros.tex")
    return metrics


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Blind VLASS proper-motion search (plan 64).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
