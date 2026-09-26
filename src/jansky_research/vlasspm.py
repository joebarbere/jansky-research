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
    "collinearity_test",
    "completeness",
    "flux_consistent",
    "inject_movers",
    "link_pairs",
    "orphan_masks",
    "run",
    "scramble_null",
    "search",
    "surface_density_limit",
    "synthetic_epochs",
    "tangent_offsets_arcsec",
]

ARCSEC = 1.0 / 3600.0
MU_MIN_ARCSEC_YR = 0.3  # plan 64 annulus; its low edge sits at the E1 astrometric floor
MU_MAX_ARCSEC_YR = 5.0
STATIC_MATCH_ARCSEC = 2.5  # one VLASS beam: a counterpart this close means "did not move"
E3_TOL_SIGMA = 3.0  # E3 must lie within this many combined sigma of the predicted position
FLUX_RATIO_MAX = 3.0  # max/min peak flux across the three detections
POISSON95_ZERO = 2.996  # 95% one-sided upper limit on a Poisson mean when 0 events are seen
DAYS_PER_YEAR = 365.25


@dataclass
class EpochCatalog:
    """One epoch's component list. ``t_yr`` is per component (VLASS tiles span years)."""

    ra: np.ndarray  # deg
    dec: np.ndarray  # deg
    t_yr: np.ndarray  # decimal year of the observation of each component
    flux: np.ndarray  # peak flux, mJy/beam (scale-corrected)
    pos_err: np.ndarray  # 1-sigma positional error, arcsec (catalogue + astrometric floor)
    ident: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))

    def __post_init__(self) -> None:
        for k in ("ra", "dec", "t_yr", "flux", "pos_err"):
            setattr(self, k, np.asarray(getattr(self, k), dtype=float))
        if self.ident.size == 0:
            self.ident = np.full(self.ra.size, -1, dtype=np.int64)

    def __len__(self) -> int:
        return int(self.ra.size)

    def subset(self, mask: np.ndarray) -> EpochCatalog:
        return EpochCatalog(
            self.ra[mask],
            self.dec[mask],
            self.t_yr[mask],
            self.flux[mask],
            self.pos_err[mask],
            self.ident[mask],
        )

    @staticmethod
    def concat(parts: list[EpochCatalog]) -> EpochCatalog:
        return EpochCatalog(
            *(np.concatenate([getattr(p, k) for p in parts]) for k in _FIELDS),
        )


_FIELDS = ("ra", "dec", "t_yr", "flux", "pos_err", "ident")


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


def link_pairs(
    a: EpochCatalog,
    b: EpochCatalog,
    *,
    mu_min: float = MU_MIN_ARCSEC_YR,
    mu_max: float = MU_MAX_ARCSEC_YR,
    n_sigma_min: float = 3.0,
) -> Pairs:
    """All (a, b) pairs whose implied proper motion lies in [mu_min, mu_max] arcsec/yr.

    The shift must also be significant: ``sep >= n_sigma_min * hypot(err_a, err_b)``, so an
    astrometric-noise offset of a static source cannot enter even if it slipped past the orphan
    cut. Per-component epochs make ``dt`` pair-specific.
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
    sig = np.hypot(a.pos_err[ii], b.pos_err[jj])
    keep = (dt > 0) & (mu >= mu_min) & (mu <= mu_max) & (sep >= n_sigma_min * sig)
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
) -> Triplets:
    """Keep pairs that an E3 component continues on a straight line at the same rate.

    For each pair, E3 components near the E2 position are candidates; each is tested against the
    position the pair's motion predicts at *that component's own* epoch. The tolerance combines
    the E3 error with the prediction error (the pair's two positional errors propagated forward).
    Where several E3 components qualify, the best (smallest residual) is kept.
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
        s = np.sqrt((a.pos_err[i] * r) ** 2 + (b.pos_err[j] * (1.0 + r)) ** 2 + c.pos_err[ks] ** 2)
        res = np.hypot(ox - px, oy - py) / s
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


def flux_consistent(
    fa: np.ndarray, fb: np.ndarray, fc: np.ndarray, *, max_ratio: float = FLUX_RATIO_MAX
) -> np.ndarray:
    """True where max/min of the three peak fluxes is within ``max_ratio``."""
    f = np.column_stack([fa, fb, fc]).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = f.max(axis=1) / f.min(axis=1)
    return np.isfinite(r) & (r <= max_ratio)


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
) -> SearchResult:
    """Steps 1-4: orphans, E1 x E2 linkage, E3 collinearity, flux consistency."""
    o1, o2 = orphan_masks(e1, e2)
    # An E3 orphan must be unmatched in BOTH earlier epochs (a mover's E3 position is new sky).
    o3a, _ = orphan_masks(e3, e1)
    o3b, _ = orphan_masks(e3, e2)
    a, b, c = e1.subset(o1), e2.subset(o2), e3.subset(o3a & o3b)
    pairs = link_pairs(a, b, mu_min=mu_min, mu_max=mu_max)
    trip = collinearity_test(a, b, c, pairs, mu_max=mu_max)
    fok = flux_consistent(a.flux[trip.i], b.flux[trip.j], c.flux[trip.k])
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
    return EpochCatalog(
        (cat.ra + shift_deg) % 360.0, cat.dec, cat.t_yr, cat.flux, cat.pos_err, cat.ident
    )


def scramble_null(
    orphans: tuple[EpochCatalog, EpochCatalog, EpochCatalog],
    *,
    n_reps: int = 20,
    shift_arcmin: tuple[float, float] = (3.0, 30.0),
    seed: int = 0,
    mu_min: float = MU_MIN_ARCSEC_YR,
    mu_max: float = MU_MAX_ARCSEC_YR,
) -> dict:
    """Chance pairs/triplets/candidates per scramble, from RA-rotated E2 and E3 orphans.

    E2 and E3 are rotated by *independent* random shifts, so a real mover can survive neither
    the pairing nor the E3 line; the orphan density (and its Dec dependence) is preserved.
    """
    a, b, c = orphans
    rng = np.random.default_rng(seed)
    n_pairs, n_trip, n_cand = [], [], []
    for _ in range(n_reps):
        s2, s3 = rng.uniform(*shift_arcmin, size=2) / 60.0 * rng.choice([-1, 1], size=2)
        bb, cc = _rotate_ra(b, float(s2)), _rotate_ra(c, float(s3))
        p = link_pairs(a, bb, mu_min=mu_min, mu_max=mu_max)
        t = collinearity_test(a, bb, cc, p, mu_max=mu_max)
        f = flux_consistent(a.flux[t.i], bb.flux[t.j], cc.flux[t.k])
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


def inject_movers(
    e1: EpochCatalog,
    e2: EpochCatalog,
    e3: EpochCatalog,
    *,
    n: int,
    mu_range: tuple[float, float] = (MU_MIN_ARCSEC_YR, MU_MAX_ARCSEC_YR),
    flux_mjy: float = 3.0,
    seed: int = 0,
    first_ident: int = 10_000_000,
) -> tuple[EpochCatalog, EpochCatalog, EpochCatalog, np.ndarray]:
    """Plant ``n`` movers into copies of the three catalogues; return them + injected rates.

    Each mover borrows the sky position and per-epoch observation times of a random real E1
    component (so it inherits the real tiling, baselines and Dec distribution), is offset 60-120"
    from it (into empty sky, not onto a static source), and gets a log-uniform rate in
    ``mu_range`` with a random direction. Detected positions scatter by that epoch's typical
    positional error. Injected rows carry ``ident >= first_ident``.
    """
    rng = np.random.default_rng(seed)
    host = rng.integers(0, len(e1), n)
    off = rng.uniform(60.0, 120.0, n) * ARCSEC
    ang = rng.uniform(0, 2 * np.pi, n)
    ra0 = (e1.ra[host] + off * np.cos(ang) / np.cos(np.radians(e1.dec[host]))) % 360.0
    dec0 = np.clip(e1.dec[host] + off * np.sin(ang), -89.9, 89.9)
    mu = np.exp(rng.uniform(np.log(mu_range[0]), np.log(mu_range[1]), n))
    phi = rng.uniform(0, 2 * np.pi, n)
    mu_ra, mu_dec = mu * np.cos(phi), mu * np.sin(phi)
    ids = first_ident + np.arange(n)

    def _epoch_times(cat: EpochCatalog) -> np.ndarray:
        # The time the injected source is observed in this epoch: the nearest real component's.
        _, k = cKDTree(_xyz(cat.ra, cat.dec)).query(_xyz(ra0, dec0), k=1)
        return cat.t_yr[k]

    t1 = e1.t_yr[host]
    out = []
    for cat in (e1, e2, e3):
        t = t1 if cat is e1 else _epoch_times(cat)
        err = np.full(n, float(np.median(cat.pos_err)))
        ra, dec = _mover_track(ra0, dec0, mu_ra, mu_dec, t1, t)
        ra = (ra + rng.normal(0, 1, n) * err * ARCSEC / np.cos(np.radians(dec))) % 360.0
        dec = dec + rng.normal(0, 1, n) * err * ARCSEC
        flux = flux_mjy * rng.uniform(0.8, 1.25, n)
        out.append(EpochCatalog.concat([cat, EpochCatalog(ra, dec, t, flux, err, ids)]))
    return out[0], out[1], out[2], mu


def completeness(
    e1: EpochCatalog,
    e2: EpochCatalog,
    e3: EpochCatalog,
    *,
    n: int = 2000,
    bins: np.ndarray | None = None,
    seed: int = 0,
    flux_mjy: float = 3.0,
) -> dict:
    """Fraction of injected movers recovered as candidates, overall and per log-rate bin."""
    bins = np.geomspace(MU_MIN_ARCSEC_YR, MU_MAX_ARCSEC_YR, 6) if bins is None else bins
    first = 10_000_000
    i1, i2, i3, mu = inject_movers(e1, e2, e3, n=n, seed=seed, flux_mjy=flux_mjy, first_ident=first)
    res = search(i1, i2, i3)
    a, b, c = res.orphans
    cand = res.candidates
    # A recovery = a candidate whose three components are the SAME injected source.
    same = (
        (a.ident[cand.i] >= first)
        & (a.ident[cand.i] == b.ident[cand.j])
        & (b.ident[cand.j] == c.ident[cand.k])
    )
    found = np.zeros(n, dtype=bool)
    found[a.ident[cand.i][same] - first] = True
    idx = np.digitize(mu, bins) - 1
    per_bin = [
        float(found[idx == k].mean()) if np.any(idx == k) else float("nan")
        for k in range(len(bins) - 1)
    ]
    return {
        "n_injected": n,
        "flux_mjy": flux_mjy,
        "overall": float(found.mean()),
        "bin_edges": [float(x) for x in bins],
        "per_bin": per_bin,
    }


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
    pos_err: tuple[float, float, float] = (0.5, 0.3, 0.15),
    seed: int = 0,
) -> tuple[EpochCatalog, EpochCatalog, EpochCatalog, np.ndarray]:
    """Offline fixture: statics, planted movers, and variable-source orphans in three epochs.

    Positional errors default to the real per-epoch floors (E1 ~0.5", E2 ~0.3", E3 ~0.15";
    VLASS Memo 22). Variable sources appear in a random subset of epochs at random positions,
    which is what produces chance orphan pairs in real data. Returns the epochs and the planted
    movers' true rates (arcsec/yr); movers carry ``ident`` 0..n_movers-1.
    """
    rng = np.random.default_rng(seed)
    half = area_side_deg / 2
    epochs_t = (2018.5, 2021.5, 2024.0)

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
    seen = rng.random((n_variable, 3)) < 0.5

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
    return cats[0], cats[1], cats[2], mu


# ------------------------------------------------------------------------------------------ run


def run(out: str = ".", *, offline: bool = True, n_null: int = 20, n_inject: int = 2000) -> dict:
    """Offline: the synthetic recover-a-known (planted movers found; null predicts chance count).

    The real leg is ``scripts/vlasspm_real.py`` (full-sky catalogues, checkpointed, run detached).
    """
    if not offline:  # pragma: no cover - the real leg is a separate, long, resumable job
        raise SystemExit("the real leg is scripts/vlasspm_real.py; run() is the offline fixture")
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
    null = scramble_null(res.orphans, n_reps=n_null)
    comp = completeness(e1, e2, e3, n=n_inject)
    metrics = {
        "source": "synthetic: planted movers + variable-source orphans (offline fixture)",
        "is_real": False,
        "syn_n_movers": n_movers,
        "syn_n_recovered": int(np.unique(a.ident[cand.i][same]).size),
        "syn_n_false": int((~same).sum()),
        "syn_n_orphans": list(res.n_orphans),
        "syn_n_pairs": len(res.pairs),
        "syn_n_triplets": len(res.triplets),
        "syn_null_candidates_mean": null["candidates_mean"],
        "syn_null_pairs_mean": null["pairs_mean"],
        "syn_completeness": comp["overall"],
        "syn_completeness_per_bin": comp["per_bin"],
        "syn_bin_edges": comp["bin_edges"],
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "vlasspm_metrics.json")
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
