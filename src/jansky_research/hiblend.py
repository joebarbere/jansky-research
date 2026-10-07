"""Does single-dish beam blending inflate HI fluxes? A two-beam test, FASHI x ALFALFA (plan 97).

FAST (FASHI) and Arecibo/ALFA (ALFALFA) observed the same galaxies with beams of different size.
Blending -- a neighbour's HI counted in the target's flux -- depends on the beam; the galaxies'
physics does not. For a target with neighbours at separation ``s`` the first-order model is

    S_obs(theta) = S_c + sum_n S_n * B_theta(s_n) * V_n,    B_theta(s) = exp(-4 ln2 s^2 / theta^2)

with ``V_n = 1`` when the neighbour's profile overlaps the target's velocity window. The
predicted two-survey log ratio ``R_pred = log10(S_obs(ALFA) / S_obs(FAST))`` is largest for
neighbours near 2 arcmin. The test statistic is the slope ``beta`` of the observed (calibrated)
log flux ratio against ``R_pred``: blending predicts beta > 0 (about 1 if the model is right),
a physical HI excess predicts beta = 0, because both surveys see the same gas.

Everything that could bias the answer was frozen in ``plans/97-hiblend-fashi-alfalfa.md`` before
the first real cross-match: the beams, the match tolerances, the samples and the five controls
(C0 power, C1 planted truth, C2 spectral negative control, C3 calibration on a disjoint isolated
sample, C4 match reliability). Changes are logged in ``survey/hiblend-findings.md``.

Pure NumPy/SciPy; the offline fixture (:func:`synthetic_field`) plants blends with known
strength so every step is testable without network.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Beams, with locators (survey/hiblend-findings.md, step 1).
FAST_FWHM_ARCMIN = 2.9  # FASHI DR2, arXiv:2606.31539 p.4 and Table 1: 2'.9 at 1420 MHz
ALFA_FWHM_ARCMIN = 3.5  # Giovanelli+2005 (AJ 130, 2598) p.22: "averaging 3.5' ... at half power"

# Frozen analysis choices (plan 97).
MATCH_RADIUS_ARCMIN = 1.5
MATCH_DV_KMS = 100.0
NEIGHBOUR_RMIN_ARCMIN = 1.0  # below this both surveys merge a pair into one source
NEIGHBOUR_RMAX_ARCMIN = 6.0  # beam response < 0.003 beyond 5', so 6' closes the sum
ISOLATION_ARCMIN = 6.0  # calibration sample: no neighbour of either catalogue within this
NULL_DV_KMS = 600.0  # C2: neighbours this far off in velocity cannot blend spectrally
SHIFT_ARCMIN = 10.0  # C4: catalogue shift used to measure the chance-match rate


def beam_response(s_arcmin: np.ndarray, fwhm_arcmin: float) -> np.ndarray:
    """Gaussian main-beam response at separation ``s`` (1 at the centre, 0.5 at FWHM/2)."""
    s = np.asarray(s_arcmin, float)
    return np.exp(-4.0 * np.log(2.0) * s**2 / fwhm_arcmin**2)


def _unit(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra, dec = np.radians(np.asarray(ra_deg, float)), np.radians(np.asarray(dec_deg, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def _chord(arcmin: float) -> float:
    return float(2.0 * np.sin(np.radians(arcmin / 60.0) / 2.0))


def _sep_arcmin(u1: np.ndarray, u2: np.ndarray) -> np.ndarray:
    c = np.linalg.norm(u1 - u2, axis=-1)
    return np.degrees(2.0 * np.arcsin(np.clip(c / 2.0, 0.0, 1.0))) * 60.0


def crossmatch(
    ra1: np.ndarray,
    dec1: np.ndarray,
    v1: np.ndarray,
    ra2: np.ndarray,
    dec2: np.ndarray,
    v2: np.ndarray,
    *,
    radius_arcmin: float = MATCH_RADIUS_ARCMIN,
    dv_kms: float = MATCH_DV_KMS,
) -> tuple[np.ndarray, np.ndarray]:
    """Mutual-nearest positional match within ``radius`` and ``|dv| <= dv_kms``.

    Returns index arrays ``(i1, i2)``. A pair is kept only if each is the other's nearest
    velocity-consistent counterpart, so a target with a close neighbour is not matched to it.
    """
    from scipy.spatial import cKDTree

    u1, u2 = _unit(ra1, dec1), _unit(ra2, dec2)
    v1, v2 = np.asarray(v1, float), np.asarray(v2, float)

    def nearest(ua, va, ub, vb) -> np.ndarray:
        best = np.full(len(ua), -1, int)
        hits = cKDTree(ub).query_ball_point(ua, r=_chord(radius_arcmin))
        for i, h in enumerate(hits):
            if not h:
                continue
            h = np.asarray(h, int)
            h = h[np.abs(vb[h] - va[i]) <= dv_kms]
            if h.size:
                best[i] = int(h[np.argmin(_sep_arcmin(ub[h], ua[i]))])
        return best

    b12 = nearest(u1, v1, u2, v2)
    b21 = nearest(u2, v2, u1, v1)
    i1 = np.flatnonzero(b12 >= 0)
    i2 = b12[i1]
    keep = b21[i2] == i1
    return i1[keep], i2[keep]


@dataclass
class Neighbours:
    """Target-neighbour pairs within ``NEIGHBOUR_RMAX_ARCMIN`` (flat arrays, one row per pair)."""

    target: np.ndarray  # index into the target arrays
    sep_arcmin: np.ndarray
    dv_kms: np.ndarray
    flux: np.ndarray  # neighbour HI flux, Jy km/s
    overlap: np.ndarray  # bool: profiles overlap in velocity (V_n = 1)


def find_neighbours(
    t_ra: np.ndarray,
    t_dec: np.ndarray,
    t_v: np.ndarray,
    t_w50: np.ndarray,
    n_ra: np.ndarray,
    n_dec: np.ndarray,
    n_v: np.ndarray,
    n_w50: np.ndarray,
    n_flux: np.ndarray,
    *,
    self_index: np.ndarray | None = None,
    rmax_arcmin: float = NEIGHBOUR_RMAX_ARCMIN,
) -> Neighbours:
    """Every neighbour within ``rmax`` of every target, in any velocity.

    ``self_index[i]`` is the target's own row in the neighbour catalogue (or -1), excluded so a
    target is never its own neighbour.
    """
    from scipy.spatial import cKDTree

    ut, un = _unit(t_ra, t_dec), _unit(n_ra, n_dec)
    hits = cKDTree(un).query_ball_point(ut, r=_chord(rmax_arcmin))
    rows_t, rows_n = [], []
    for i, h in enumerate(hits):
        for j in h:
            if self_index is not None and j == self_index[i]:
                continue
            rows_t.append(i)
            rows_n.append(j)
    it, jn = np.asarray(rows_t, int), np.asarray(rows_n, int)
    tv, nv = np.asarray(t_v, float), np.asarray(n_v, float)
    dv = nv[jn] - tv[it] if it.size else np.empty(0)
    half = (
        np.nan_to_num(np.asarray(t_w50, float)[it]) + np.nan_to_num(np.asarray(n_w50, float)[jn])
    ) / 2
    return Neighbours(
        target=it,
        sep_arcmin=_sep_arcmin(ut[it], un[jn]) if it.size else np.empty(0),
        dv_kms=dv,
        flux=np.asarray(n_flux, float)[jn] if it.size else np.empty(0),
        overlap=np.abs(dv) < half if it.size else np.zeros(0, bool),
    )


def predicted_log_ratio(
    flux_c: np.ndarray,
    nb: Neighbours,
    *,
    use: np.ndarray | None = None,
    force_overlap: bool = False,
    fwhm_a: float = ALFA_FWHM_ARCMIN,
    fwhm_f: float = FAST_FWHM_ARCMIN,
) -> np.ndarray:
    """R_pred per target: log10 of the ALFA-to-FAST ratio the beam model predicts.

    ``use`` selects which neighbour pairs enter the sum (default: all); ``force_overlap`` treats
    every used pair as spectrally overlapping, which C2 needs (what blending *would* predict for
    neighbours that in fact cannot blend).
    """
    sc = np.asarray(flux_c, float)
    sel = np.ones(nb.target.size, bool) if use is None else np.asarray(use, bool)
    v = np.ones(nb.target.size, bool) if force_overlap else nb.overlap
    w = sel & v
    add_a = np.zeros(sc.size)
    add_f = np.zeros(sc.size)
    np.add.at(add_a, nb.target[w], nb.flux[w] * beam_response(nb.sep_arcmin[w], fwhm_a))
    np.add.at(add_f, nb.target[w], nb.flux[w] * beam_response(nb.sep_arcmin[w], fwhm_f))
    return np.log10((sc + add_a) / (sc + add_f))


def target_flux_estimate(
    flux_a: np.ndarray,
    flux_f: np.ndarray,
    sig_a_dex: np.ndarray,
    sig_f_dex: np.ndarray,
    nb: Neighbours,
    *,
    subtract_blend: bool = True,
) -> np.ndarray:
    """Unblended target flux for R_pred, with noise uncorrelated with the observed log ratio.

    Using either survey's flux as S_c puts that survey's measurement noise into both the response
    (log S_A - log S_F) and the regressor (R_pred falls as S_c rises): the noise then manufactures
    a slope (beta came out 3.4 for a planted 1 on the synthetic field). The weighted geometric mean
    log S = w log S_A + (1 - w) log S_F with w = s_F^2 / (s_A^2 + s_F^2) has
    cov(e_A - e_F, w e_A + (1 - w) e_F) = 0, so it is uncorrelated with the response noise to first
    order. The beam model's own blend contribution is then subtracted, so S_c estimates the
    unblended flux the model's S_c denotes. ``subtract_blend=False`` returns the geometric mean
    alone -- the version the C3 covariates must use, because subtracting a blend that may not
    exist shifts the covariates with R_pred by construction and biases beta under the null.
    """
    sa2 = np.asarray(sig_a_dex, float) ** 2
    sf2 = np.asarray(sig_f_dex, float) ** 2
    w = np.where(sa2 + sf2 > 0, sf2 / (sa2 + sf2), 0.5)
    logs = w * np.log10(flux_a) + (1.0 - w) * np.log10(flux_f)
    est = 10**logs
    blend = np.zeros(est.size)
    ov = nb.overlap
    contrib = nb.flux[ov] * (
        w[nb.target[ov]] * beam_response(nb.sep_arcmin[ov], ALFA_FWHM_ARCMIN)
        + (1.0 - w[nb.target[ov]]) * beam_response(nb.sep_arcmin[ov], FAST_FWHM_ARCMIN)
    )
    np.add.at(blend, nb.target[ov], contrib)
    if not subtract_blend:
        return est
    return np.maximum(est - blend, 0.1 * est)


def sample_masks(n_targets: int, nb: Neighbours) -> dict[str, np.ndarray]:
    """The three frozen samples, as boolean masks over targets.

    - ``isolated`` (C3 calibration): no neighbour within ``ISOLATION_ARCMIN``, any velocity.
    - ``primary`` (beta): at least one velocity-overlapping neighbour at 1-6', and no neighbour
      inside 1' (an unresolved pair is merged by both surveys).
    - ``null`` (C2): no overlapping neighbour within 6', and at least one neighbour at 1-6'
      offset by more than ``NULL_DV_KMS``.
    """
    close = np.zeros(n_targets, bool)
    close[nb.target[nb.sep_arcmin < NEIGHBOUR_RMIN_ARCMIN]] = True
    any6 = np.zeros(n_targets, bool)
    any6[nb.target[nb.sep_arcmin <= ISOLATION_ARCMIN]] = True
    ring = (nb.sep_arcmin >= NEIGHBOUR_RMIN_ARCMIN) & (nb.sep_arcmin <= NEIGHBOUR_RMAX_ARCMIN)
    ov = np.zeros(n_targets, bool)
    ov[nb.target[ring & nb.overlap]] = True
    ov_any = np.zeros(n_targets, bool)
    ov_any[nb.target[nb.overlap & (nb.sep_arcmin <= NEIGHBOUR_RMAX_ARCMIN)]] = True
    far = np.zeros(n_targets, bool)
    far[nb.target[ring & (np.abs(nb.dv_kms) > NULL_DV_KMS)]] = True
    return {
        "isolated": ~any6,
        "primary": ov & ~close,
        "null": far & ~ov_any & ~close,
    }


def null_pairs(nb: Neighbours) -> np.ndarray:
    """The C2 neighbour pairs: in the 1-6' ring and offset by more than ``NULL_DV_KMS``."""
    ring = (nb.sep_arcmin >= NEIGHBOUR_RMIN_ARCMIN) & (nb.sep_arcmin <= NEIGHBOUR_RMAX_ARCMIN)
    return ring & (np.abs(nb.dv_kms) > NULL_DV_KMS)


def calibration_design(
    log_snr: np.ndarray, log_w50: np.ndarray, log_flux: np.ndarray, dec: np.ndarray
) -> np.ndarray:
    """C3 covariates for the survey-to-survey flux scale (S/N enters to second order, for the
    Eddington-bias offset FASHI DR2 reports at SNR <~ 20, arXiv:2606.31539 sec. 5.2)."""
    x = [
        np.ones(len(log_snr)),
        log_snr,
        log_snr**2,
        log_w50,
        log_flux,
        np.asarray(dec, float) / 90.0,
    ]
    return np.column_stack(x)


def calibration_design_v2(
    log_snr: np.ndarray,
    log_w50: np.ndarray,
    log_flux: np.ndarray,
    dec: np.ndarray,
    log_n15: np.ndarray,
    knots: np.ndarray,
) -> np.ndarray:
    """Plan 98 C3 covariates: a linear spline in log S/N (hinges at ``knots``, the calibration
    sample's log-S/N deciles) replaces plan 97's quadratic, plus log(1 + N15) (D3)."""
    hinges = [np.clip(log_snr - k, 0, None) for k in np.asarray(knots, float)]
    x = [np.ones(len(log_snr)), log_snr, *hinges, log_w50, log_flux, np.asarray(dec, float) / 90.0,
         np.asarray(log_n15, float)]  # fmt: skip
    return np.column_stack(x)


def snr_knots(log_snr: np.ndarray) -> np.ndarray:
    """Plan 98's frozen knots: the 10th..90th percentiles of log S/N in the calibration sample."""
    return np.nanquantile(log_snr, np.linspace(0.1, 0.9, 9))


def local_density(
    ra: np.ndarray,
    dec: np.ndarray,
    cat_ra: np.ndarray,
    cat_dec: np.ndarray,
    *,
    radius: float = 15.0,
) -> np.ndarray:
    """log10(1 + N), N = catalogued HI sources within ``radius`` arcmin, excluding the target
    itself (D3; a plan-98 calibration covariate)."""
    from scipy.spatial import cKDTree

    tree = cKDTree(_unit(cat_ra, cat_dec))
    hits = tree.query_ball_point(_unit(ra, dec), r=_chord(radius))
    n = np.array([len(h) for h in hits]) - 1
    return np.log10(1.0 + np.clip(n, 0, None))


def fit_calibration(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    """Least-squares coefficients of the log flux ratio on the C3 covariates."""
    ok = np.all(np.isfinite(design), axis=1) & np.isfinite(y)
    coef, *_ = np.linalg.lstsq(design[ok], y[ok], rcond=None)
    return coef


def target_clusters(
    ra: np.ndarray, dec: np.ndarray, *, link_arcmin: float = NEIGHBOUR_RMAX_ARCMIN
) -> np.ndarray:
    """Connected components of targets linked within ``link_arcmin`` (bootstrap units)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    u = _unit(ra, dec)
    pairs = cKDTree(u).query_pairs(r=_chord(link_arcmin), output_type="ndarray")
    n = len(u)
    if pairs.size == 0:
        return np.arange(n)
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _k, labels = connected_components(g, directed=False)
    return labels


def fit_beta(
    y: np.ndarray,
    r_pred: np.ndarray,
    clusters: np.ndarray,
    rng: np.random.Generator,
    *,
    n_boot: int = 1000,
) -> dict:
    """OLS slope of the calibrated log ratio on R_pred, with a cluster bootstrap error."""
    y, r, cl = np.asarray(y, float), np.asarray(r_pred, float), np.asarray(clusters)
    ok = np.isfinite(y) & np.isfinite(r)
    y, r, cl = y[ok], r[ok], cl[ok]
    if y.size < 3 or np.ptp(r) == 0:
        return {"n": int(y.size)}

    def slope(yy, rr):
        x = np.column_stack([np.ones(rr.size), rr])
        coef, *_ = np.linalg.lstsq(x, yy, rcond=None)
        return coef

    b0 = slope(y, r)
    labels, inv = np.unique(cl, return_inverse=True)
    members = [np.flatnonzero(inv == k) for k in range(labels.size)]
    boots = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, labels.size, labels.size)
        idx = np.concatenate([members[k] for k in pick])
        boots[b] = slope(y[idx], r[idx])[1] if np.ptp(r[idx]) > 0 else np.nan
    se = float(np.nanstd(boots, ddof=1))
    return {
        "n": int(y.size),
        "n_clusters": int(labels.size),
        "beta": round(float(b0[1]), 4),
        "beta_se": round(se, 4),
        "beta_sigma": round(float(b0[1]) / se, 2) if se > 0 else None,
        "intercept": round(float(b0[0]), 5),
        "r_pred_median": round(float(np.median(r)), 5),
        "r_pred_max": round(float(np.max(r)), 5),
    }


def blended_fluxes(
    flux_true: np.ndarray, nb: Neighbours, *, strength: float = 1.0, use: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Noise-free ALFA and FAST fluxes a target would have if blending acts with ``strength``
    times the beam model (0 = no blending, 1 = the model). Used by C0/C1 injections."""
    sel = np.ones(nb.target.size, bool) if use is None else np.asarray(use, bool)
    w = sel & nb.overlap
    add_a = np.zeros(len(flux_true))
    add_f = np.zeros(len(flux_true))
    np.add.at(add_a, nb.target[w], nb.flux[w] * beam_response(nb.sep_arcmin[w], ALFA_FWHM_ARCMIN))
    np.add.at(add_f, nb.target[w], nb.flux[w] * beam_response(nb.sep_arcmin[w], FAST_FWHM_ARCMIN))
    return flux_true + strength * add_a, flux_true + strength * add_f


def synthetic_field(
    n_targets: int = 3000,
    *,
    pair_frac: float = 0.35,
    strength: float = 1.0,
    scatter_dex: float = 0.08,
    seed: int = 0,
    snr_hinge_dex: float = 0.0,
) -> dict:
    """Offline fixture: targets, some with an HI-detected neighbour at 0.5-6', two surveys that
    blend per the beam model scaled by ``strength``, a flux-scale term that depends on S/N (the
    Eddington-like offset C3 must remove), and independent log-normal measurement scatter."""
    rng = np.random.default_rng(seed)
    ra = rng.uniform(150, 210, n_targets)
    dec = rng.uniform(0, 30, n_targets)
    v = rng.uniform(3000, 12000, n_targets)
    w50 = rng.uniform(80, 300, n_targets)
    flux = 10 ** rng.uniform(-0.5, 1.0, n_targets)
    has = rng.uniform(size=n_targets) < pair_frac
    k = int(has.sum())
    sep = rng.uniform(0.5, 6.0, k)
    ang = rng.uniform(0, 2 * np.pi, k)
    t_idx = np.flatnonzero(has)
    n_ra = ra[t_idx] + sep / 60 * np.cos(ang) / np.cos(np.radians(dec[t_idx]))
    n_dec = dec[t_idx] + sep / 60 * np.sin(ang)
    fast_v = rng.uniform(size=k) < 0.8  # most companions share the velocity window
    n_v = v[t_idx] + np.where(
        fast_v, rng.normal(0, 60, k), rng.choice([-1, 1], k) * rng.uniform(700, 1500, k)
    )
    n_w50 = rng.uniform(80, 300, k)
    n_flux = flux[t_idx] * 10 ** rng.uniform(-0.7, 0.7, k)
    nb = find_neighbours(ra, dec, v, w50, n_ra, n_dec, n_v, n_w50, n_flux)
    sa, sf = blended_fluxes(flux, nb, strength=strength)
    snr = 5 + 60 * (flux / flux.max())
    offset = 0.05 * np.exp(-(snr - 5) / 8)  # ALFALFA high at low S/N
    # optional plan-98 plant: a hinge at the top S/N quintile, the shape D0 found on real data
    ls = np.log10(snr)
    offset = offset + snr_hinge_dex * np.clip(ls - np.quantile(ls, 0.8), 0, None) / np.ptp(ls)
    obs_a = sa * 10 ** (offset + rng.normal(0, scatter_dex, n_targets))
    obs_f = sf * 10 ** rng.normal(0, scatter_dex, n_targets)
    return {
        "ra": ra, "dec": dec, "v": v, "w50": w50, "flux_true": flux, "snr": snr,
        "flux_a": obs_a, "flux_f": obs_f, "nb": nb,
        "flux_err_f": flux / snr,
        "sig_a_dex": np.full(n_targets, scatter_dex), "sig_f_dex": np.full(n_targets, scatter_dex),
        "log_n15": local_density(ra, dec, np.concatenate([ra, n_ra]), np.concatenate([dec, n_dec])),
    }  # fmt: skip


def _design(field: dict, log_snr: np.ndarray, s_g: np.ndarray, calib: np.ndarray) -> np.ndarray:
    """The C3 design matrix: plan 97's form, or plan 98's when ``field["calib"] == "v2"`` (knots
    from the rows in ``calib``, the sample the calibration is fitted on)."""
    if field.get("calib", "v1") == "v2":
        return calibration_design_v2(
            log_snr, np.log10(field["w50"]), np.log10(s_g), field["dec"], field["log_n15"],
            snr_knots(log_snr[calib]),
        )  # fmt: skip
    return calibration_design(log_snr, np.log10(field["w50"]), np.log10(s_g), field["dec"])


def _prepare(field: dict) -> dict:
    """Everything step 3 computes before fitting: samples, C3 calibration, residuals, R_pred."""
    nb: Neighbours = field["nb"]
    n = len(field["flux_f"])
    masks = sample_masks(n, nb)
    y = np.log10(field["flux_a"] / field["flux_f"])
    s_c = target_flux_estimate(
        field["flux_a"], field["flux_f"], field["sig_a_dex"], field["sig_f_dex"], nb
    )
    # C3 covariates assume nothing about blending: the noise-decorrelated geometric-mean flux, NOT
    # blend-subtracted (subtracting a blend that may not exist shifts the covariates with R_pred
    # and biased beta to -0.32 +/- 0.08 under the null), and an S/N from the flux ERROR (a noise
    # level) -- dividing by a measured flux puts that survey's noise back into the covariate.
    # If blending is real the covariates carry some of it and the calibration absorbs a little:
    # that biases beta toward 0, the conservative direction, and C1 measures by how much.
    s_g = target_flux_estimate(
        field["flux_a"], field["flux_f"], field["sig_a_dex"], field["sig_f_dex"], nb,
        subtract_blend=False,
    )  # fmt: skip
    snr_g = s_g / field["flux_err_f"]
    design = _design(field, np.log10(snr_g), s_g, masks["isolated"])
    coef = fit_calibration(y[masks["isolated"]], design[masks["isolated"]])
    return {
        "n": n, "masks": masks, "y": y, "s_c": s_c, "s_g": s_g, "snr_g": snr_g,
        "design": design, "coef": coef, "resid": y - design @ coef,
        "clusters": target_clusters(field["ra"], field["dec"]),
        "r_pred": predicted_log_ratio(s_c, nb),
        "r_null": predicted_log_ratio(s_c, nb, use=null_pairs(nb), force_overlap=True),
    }  # fmt: skip


def analyse(field: dict, rng: np.random.Generator, *, n_boot: int = 500) -> dict:
    """The full frozen analysis on a matched catalogue: C3 calibration on isolated targets,
    beta on the primary sample, beta_null on the C2 sample."""
    d = _prepare(field)
    p, q = d["masks"]["primary"], d["masks"]["null"]
    return {
        "n_targets": d["n"],
        "n_isolated": int(d["masks"]["isolated"].sum()),
        "n_primary": int(p.sum()),
        "n_null": int(q.sum()),
        "calibration_coef": [round(float(c), 5) for c in d["coef"]],
        "primary": fit_beta(d["resid"][p], d["r_pred"][p], d["clusters"][p], rng, n_boot=n_boot),
        "null": fit_beta(d["resid"][q], d["r_null"][q], d["clusters"][q], rng, n_boot=n_boot),
    }


def fit_ols(
    y: np.ndarray, x: np.ndarray, clusters: np.ndarray, rng: np.random.Generator, *,
    weights: np.ndarray | None = None, n_boot: int = 500,
) -> dict:  # fmt: skip
    """(Weighted) least squares of y on [1, x...] with cluster-bootstrap standard errors."""
    x = np.column_stack([np.ones(len(y)), np.asarray(x, float).reshape(len(y), -1)])
    w = np.ones(len(y)) if weights is None else np.asarray(weights, float)
    ok = np.isfinite(y) & np.all(np.isfinite(x), axis=1) & (w > 0)
    y, x, w, cl = np.asarray(y, float)[ok], x[ok], w[ok], np.asarray(clusters)[ok]
    if y.size <= x.shape[1] + 2:
        return {"n": int(y.size)}

    def solve(idx):
        sw = np.sqrt(w[idx])
        coef, *_ = np.linalg.lstsq(x[idx] * sw[:, None], y[idx] * sw, rcond=None)
        return coef

    c0 = solve(np.arange(y.size))
    labels, inv = np.unique(cl, return_inverse=True)
    members = [np.flatnonzero(inv == k) for k in range(labels.size)]
    boots = np.array([
        solve(np.concatenate([members[k] for k in rng.integers(0, labels.size, labels.size)]))
        for _ in range(n_boot)
    ])  # fmt: skip
    se = boots.std(axis=0, ddof=1)
    return {
        "n": int(y.size),
        "coef": [round(float(v), 5) for v in c0],
        "se": [round(float(v), 5) for v in se],
        "sigma": [round(float(v / s), 2) if s > 0 else None for v, s in zip(c0, se, strict=True)],
    }


def _binned_median(v: np.ndarray, key: np.ndarray, n_bins: int) -> list[dict]:
    edges = np.nanquantile(key, np.linspace(0, 1, n_bins + 1))
    out = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        s = (key >= lo) & (key <= hi) & np.isfinite(v)
        vv = v[s]
        se = 1.2533 * vv.std(ddof=1) / np.sqrt(vv.size) if vv.size > 1 else np.nan
        out.append({"lo": round(float(lo), 3), "hi": round(float(hi), 3), "n": int(vv.size),
                    "median": round(float(np.median(vv)), 5) if vv.size else None,
                    "se": round(float(se), 5), "sigma": round(float(np.median(vv) / se), 2) if vv.size > 1 and se > 0 else None})  # fmt: skip
    return out


def diagnostics(
    field: dict,
    rng: np.random.Generator,
    *,
    n_cat_ra: np.ndarray,
    n_cat_dec: np.ndarray,
    n_boot: int = 500,
) -> dict:
    """Post-hoc D0-D3 of the C2 failure (survey/hiblend-findings.md step 4)."""
    d = _prepare(field)
    m, resid, cl = d["masks"], d["resid"], d["clusters"]
    iso, q, p = np.flatnonzero(m["isolated"]), m["null"], m["primary"]
    log_snr = np.log10(d["snr_g"])
    out: dict = {}
    # D0: calibrate on half the isolated sample, inspect the other half by S/N
    perm = rng.permutation(iso)
    a, b = perm[: perm.size // 2], perm[perm.size // 2 :]
    coef_a = fit_calibration(d["y"][a], d["design"][a])
    res_b = d["y"][b] - d["design"][b] @ coef_a
    out["D0_heldout_isolated_by_snr"] = _binned_median(res_b, log_snr[b], 5)
    # D1: null reweighted to the primary S/N distribution; and by null S/N tercile
    edges = np.nanquantile(log_snr[p], np.linspace(0, 1, 11))
    hp = np.histogram(log_snr[p], edges)[0] / p.sum()
    hq = np.histogram(log_snr[q], edges)[0] / max(q.sum(), 1)
    k = np.clip(np.digitize(log_snr, edges) - 1, 0, 9)
    wq = np.where(hq[k] > 0, hp[k] / np.where(hq[k] > 0, hq[k], 1), 0.0)
    out["D1_null_snr_matched"] = fit_ols(
        resid[q], d["r_null"][q], cl[q], rng, weights=wq[q], n_boot=n_boot
    )
    terc = np.nanquantile(log_snr[q], [0, 1 / 3, 2 / 3, 1])
    out["D1_null_by_snr_tercile"] = []
    for lo, hi in zip(terc[:-1], terc[1:], strict=True):
        s = q & (log_snr >= lo) & (log_snr <= hi)
        out["D1_null_by_snr_tercile"].append({"log_snr": [round(float(lo), 3), round(float(hi), 3)],
                                              **fit_beta(resid[s], d["r_null"][s], cl[s], rng, n_boot=n_boot)})  # fmt: skip
    # D2: target flux vs brightest null-neighbour flux
    nb = field["nb"]
    npairs = null_pairs(nb)
    snmax = np.full(d["n"], np.nan)
    np.fmax.at(snmax, nb.target[npairs], nb.flux[npairs])
    x2 = np.column_stack([d["r_null"], np.log10(d["s_c"]), np.log10(snmax)])
    out["D2_null_terms"] = {"terms": ["const", "R_null", "log S_c", "log S_n,max"],
                            **fit_ols(resid[q], x2[q], cl[q], rng, n_boot=n_boot)}  # fmt: skip
    # D3: local density within 15' of every catalogued HI source
    ld = local_density(field["ra"], field["dec"], n_cat_ra, n_cat_dec)
    n15 = np.rint(10**ld - 1.0)
    out["D3_null_density"] = {"terms": ["const", "R_null", "log(1+N15)"],
                              **fit_ols(resid[q], np.column_stack([d["r_null"], ld])[q], cl[q], rng, n_boot=n_boot)}  # fmt: skip
    out["D3_isolated_density"] = {"terms": ["const", "log(1+N15)"],
                                  **fit_ols(resid[iso], ld[iso], cl[iso], rng, n_boot=n_boot)}  # fmt: skip
    out["N15_median"] = {"isolated": float(np.median(n15[iso])), "null": float(np.median(n15[q])),
                         "primary": float(np.median(n15[p]))}  # fmt: skip
    return out


def power(
    n_targets: int, *, strength: float = 1.0, n_real: int = 20, seed: int = 0, n_boot: int = 200
) -> dict:
    """C0: fraction of synthetic realizations where beta > 0 at >= 3 sigma."""
    rng = np.random.default_rng(seed)
    hits, betas = 0, []
    for k in range(n_real):
        res = analyse(
            synthetic_field(n_targets, strength=strength, seed=seed + 1000 + k), rng, n_boot=n_boot
        )
        b = res["primary"]
        betas.append(b.get("beta", np.nan))
        if (b.get("beta_sigma") or 0) >= 3:
            hits += 1
    return {
        "n_real": n_real,
        "strength": strength,
        "detect_frac": round(hits / n_real, 3),
        "beta_mean": round(float(np.nanmean(betas)), 4),
        "beta_std": round(float(np.nanstd(betas, ddof=1)), 4),
    }


# ---------------------------------------------------------------------------------------------
# The real leg: catalogues -> matched targets + neighbours -> gated controls -> beta
# ---------------------------------------------------------------------------------------------


def build_field(fashi: dict, alfalfa: dict, *, shift_dec_arcmin: float = 0.0) -> dict:
    """Matched targets (FASHI x ALFALFA code 1) with their neighbours from both catalogues.

    ``fashi`` keys: ra, dec, cz, w50, flux, flux_err (Jy km/s). ``alfalfa`` keys: ra, dec, v,
    w50, flux, flux_err, code. The neighbour catalogue is every FASHI source plus every ALFALFA
    detection (code 1 or 2) with no FASHI counterpart, so no galaxy is counted twice; a
    neighbour's flux comes from FASHI where it has one (the narrower beam), else ALFALFA.
    ``shift_dec_arcmin`` displaces ALFALFA for the C4 chance-match test.
    """
    a_dec = np.asarray(alfalfa["dec"], float) + shift_dec_arcmin / 60.0
    code1 = np.asarray(alfalfa["code"]) == 1
    ia1 = np.flatnonzero(code1)
    i_f, j = crossmatch(
        fashi["ra"], fashi["dec"], fashi["cz"], alfalfa["ra"][ia1], a_dec[ia1], alfalfa["v"][ia1]
    )
    i_a = ia1[j]
    # ALFALFA detections (any code) with a FASHI counterpart are not separate neighbours
    _f_any, a_any = crossmatch(
        fashi["ra"], fashi["dec"], fashi["cz"], alfalfa["ra"], a_dec, alfalfa["v"]
    )
    a_only = np.setdiff1d(np.arange(len(a_dec)), a_any)
    n_ra = np.concatenate([fashi["ra"], alfalfa["ra"][a_only]])
    n_dec = np.concatenate([fashi["dec"], a_dec[a_only]])
    n_v = np.concatenate([fashi["cz"], alfalfa["v"][a_only]])
    n_w50 = np.concatenate([fashi["w50"], alfalfa["w50"][a_only]])
    n_flux = np.concatenate([fashi["flux"], alfalfa["flux"][a_only]])
    nb = find_neighbours(
        fashi["ra"][i_f], fashi["dec"][i_f], fashi["cz"][i_f], fashi["w50"][i_f],
        n_ra, n_dec, n_v, n_w50, n_flux, self_index=i_f,
    )  # fmt: skip
    ln10 = np.log(10.0)
    fa, fa_e = np.asarray(alfalfa["flux"], float)[i_a], np.asarray(alfalfa["flux_err"], float)[i_a]
    ff, ff_e = np.asarray(fashi["flux"], float)[i_f], np.asarray(fashi["flux_err"], float)[i_f]
    return {
        "ra": np.asarray(fashi["ra"], float)[i_f],
        "dec": np.asarray(fashi["dec"], float)[i_f],
        "v": np.asarray(fashi["cz"], float)[i_f],
        "w50": np.asarray(fashi["w50"], float)[i_f],
        "flux_a": fa, "flux_f": ff, "flux_err_f": ff_e,
        "sig_a_dex": fa_e / (fa * ln10), "sig_f_dex": ff_e / (ff * ln10),
        "nb": nb, "i_fashi": i_f, "i_alfalfa": i_a,
        "n_alfalfa_code1": int(code1.sum()), "n_alfalfa_only_neighbours": int(a_only.size),
        "n_cat_ra": n_ra, "n_cat_dec": n_dec,
        "log_n15": local_density(fashi["ra"][i_f], fashi["dec"][i_f], n_ra, n_dec),
    }  # fmt: skip


def _subset(field: dict, keep: np.ndarray, nb: Neighbours) -> dict:
    out = {
        k: (v[keep] if isinstance(v, np.ndarray) and v.shape[:1] == keep.shape else v)
        for k, v in field.items()
    }
    out["nb"] = nb
    return out


def injection_field(field: dict, rng: np.random.Generator, *, strength: float) -> dict:
    """C1 on real data: split the real ISOLATED targets in half. One half stays as the
    calibration sample; the other gets one synthetic neighbour each, drawn from the real primary
    sample's (separation, neighbour/target flux ratio) pairs, and its real measured fluxes are
    blended per the model at ``strength``. Beta must come back near ``strength``."""
    masks = sample_masks(len(field["flux_f"]), field["nb"])
    iso = np.flatnonzero(masks["isolated"])
    prim = masks["primary"]
    nb = field["nb"]
    ring = prim[nb.target] & nb.overlap & (nb.sep_arcmin >= NEIGHBOUR_RMIN_ARCMIN)
    sep_pool = nb.sep_arcmin[ring]
    ratio_pool = nb.flux[ring] / field["flux_f"][nb.target[ring]]
    rng.shuffle(iso)
    half = iso.size // 2
    calib, inj = iso[:half], iso[half:]
    keep = np.concatenate([calib, inj])
    sub = _subset(field, np.isin(np.arange(len(field["flux_f"])), keep), nb)
    # re-index: positions of `inj` inside the subset
    order = np.flatnonzero(np.isin(np.arange(len(field["flux_f"])), keep))
    pos = np.searchsorted(order, inj)
    pick = rng.integers(0, sep_pool.size, inj.size)
    s_c = sub["flux_f"][pos]
    syn = Neighbours(
        target=pos, sep_arcmin=sep_pool[pick], dv_kms=np.zeros(inj.size),
        flux=ratio_pool[pick] * s_c, overlap=np.ones(inj.size, bool),
    )  # fmt: skip
    add_a = syn.flux * beam_response(syn.sep_arcmin, ALFA_FWHM_ARCMIN)
    add_f = syn.flux * beam_response(syn.sep_arcmin, FAST_FWHM_ARCMIN)
    sub["flux_a"] = sub["flux_a"].copy()
    sub["flux_f"] = sub["flux_f"].copy()
    sub["flux_a"][pos] += strength * add_a
    sub["flux_f"][pos] += strength * add_f
    sub["nb"] = syn
    return sub


def geometry_power(
    field: dict,
    rng: np.random.Generator,
    *,
    n_real: int = 20,
    strength: float = 1.0,
    n_boot: int = 200,
) -> dict:
    """C0 on the real geometry: the real targets, neighbours, neighbour fluxes and flux errors,
    with SYNTHETIC target fluxes (the geometric-mean flux as truth) blended per the model and
    re-noised with each survey's own errors. The real flux ratios are never used."""
    s_true = target_flux_estimate(
        field["flux_a"], field["flux_f"], field["sig_a_dex"], field["sig_f_dex"], field["nb"],
        subtract_blend=False,
    )  # fmt: skip
    hits, betas = 0, []
    for _k in range(n_real):
        sa, sf = blended_fluxes(s_true, field["nb"], strength=strength)
        sim = dict(field)
        sim["flux_a"] = sa * 10 ** rng.normal(0, np.nan_to_num(field["sig_a_dex"], nan=0.1))
        sim["flux_f"] = sf * 10 ** rng.normal(0, np.nan_to_num(field["sig_f_dex"], nan=0.1))
        b = analyse(sim, rng, n_boot=n_boot)["primary"]
        betas.append(b.get("beta", np.nan))
        hits += int((b.get("beta_sigma") or 0) >= 3)
    return {
        "n_real": n_real, "strength": strength, "detect_frac": round(hits / n_real, 3),
        "beta_mean": round(float(np.nanmean(betas)), 4), "beta_std": round(float(np.nanstd(betas, ddof=1)), 4),
    }  # fmt: skip


def run_gated(
    field: dict,
    shifted: dict,
    rng: np.random.Generator,
    *,
    n_power: int = 20,
    n_inj: int = 10,
    n_boot: int = 1000,
) -> dict:
    """The frozen order: C0 power -> C1 planted truth -> C4 match reliability -> (C3 + C2 + beta).
    A failed gate stops the run and records why; beta is computed only if all gates pass."""
    out: dict = {"n_targets": len(field["flux_f"]), "gates": {}}
    c0 = geometry_power(field, rng, n_real=n_power, n_boot=200)
    c0_null = geometry_power(field, rng, n_real=n_power, strength=0.0, n_boot=200)
    out["gates"]["C0_power"] = {
        "planted_1": c0,
        "planted_0": c0_null,
        "pass": c0["detect_frac"] >= 0.8,
    }
    if not out["gates"]["C0_power"]["pass"]:
        out["stopped_at"] = "C0: underpowered (detect_frac < 0.8); no real result is quoted"
        return out
    inj1 = [
        analyse(injection_field(field, rng, strength=1.0), rng, n_boot=300)["primary"]
        for _ in range(n_inj)
    ]
    inj0 = [
        analyse(injection_field(field, rng, strength=0.0), rng, n_boot=300)["primary"]
        for _ in range(n_inj)
    ]
    b1 = np.array([r["beta"] for r in inj1])
    b0 = np.array([r["beta"] for r in inj0])
    se0 = np.array([r["beta_se"] for r in inj0])
    c1_pass = bool(abs(b1.mean() - 1.0) <= 0.3 and abs(b0.mean()) < 2 * np.mean(se0))
    out["gates"]["C1_planted"] = {
        "n_inj": n_inj, "beta_planted_1_mean": round(float(b1.mean()), 4), "beta_planted_1_sd": round(float(b1.std(ddof=1)), 4),
        "beta_planted_0_mean": round(float(b0.mean()), 4), "beta_planted_0_sd": round(float(b0.std(ddof=1)), 4),
        "pass": c1_pass,
    }  # fmt: skip
    if not c1_pass:
        out["stopped_at"] = "C1: planted truth not recovered (|beta1-1| > 0.3 or beta0 != 0)"
        return out
    chance = len(shifted["flux_f"]) / max(len(field["flux_f"]), 1)
    out["gates"]["C4_match"] = {"n_matched": len(field["flux_f"]), "n_matched_shifted": len(shifted["flux_f"]),
                               "chance_rate": round(chance, 5), "pass": chance < 0.01}  # fmt: skip
    if chance >= 0.01:
        out["stopped_at"] = "C4: chance-match rate >= 1%"
        return out
    res = analyse(field, rng, n_boot=n_boot)
    nul = res["null"]
    c2_pass = bool(nul.get("beta_se") and abs(nul["beta"]) < 2 * nul["beta_se"])
    out["gates"]["C2_spectral_null"] = {**nul, "pass": c2_pass}
    out["C3_calibration_coef"] = res["calibration_coef"]
    out["samples"] = {k: res[k] for k in ("n_isolated", "n_primary", "n_null")}
    out["primary"] = res["primary"]
    b = res["primary"]
    if not c2_pass:
        out["outcome"] = "ambiguous: C2 spectral null control is non-zero"
    elif (b.get("beta_sigma") or 0) >= 3:
        out["outcome"] = "blending supported: beta > 0 at >= 3 sigma, C1 and C2 pass"
    else:
        out["outcome"] = "blending not supported: beta consistent with 0"
        out["beta_upper_95"] = round(b["beta"] + 1.645 * b["beta_se"], 4)
    return out


def heldout_calibration_check(
    field: dict, rng: np.random.Generator, *, n_bins: int = 10, n_boot: int = 200,
    max_abs_dex: float = 0.010, min_p: float = 0.01, strip_deg: float = 2.0,
) -> dict:  # fmt: skip
    """C3' (plan 98): fit the calibration on isolated targets in even RA strips, test on odd,
    and the reverse. Each direction passes if every held-out S/N-decile median residual is under
    ``max_abs_dex`` and the chi^2 of the medians against 0 has p > ``min_p``."""
    from scipy.stats import chi2

    d = _prepare(field)
    iso = d["masks"]["isolated"]
    log_snr = np.log10(d["snr_g"])
    parity = np.floor(np.asarray(field["ra"], float) / strip_deg).astype(int) % 2
    out: dict = {"form": field.get("calib", "v1"), "directions": []}
    for train in (0, 1):
        tr, te = iso & (parity == train), iso & (parity != train)
        design = _design(field, log_snr, d["s_g"], tr)
        coef = fit_calibration(d["y"][tr], design[tr])
        res, key = (d["y"] - design @ coef)[te], log_snr[te]
        ok = np.isfinite(res) & np.isfinite(key)
        res, key = res[ok], key[ok]
        edges = np.quantile(key, np.linspace(0, 1, n_bins + 1))
        k = np.clip(np.searchsorted(edges, key, side="right") - 1, 0, n_bins - 1)
        bins = []
        for b in range(n_bins):
            v = res[k == b]
            med = float(np.median(v))
            se = float(np.std([np.median(rng.choice(v, v.size)) for _ in range(n_boot)], ddof=1))
            bins.append({"lo": round(float(edges[b]), 3), "hi": round(float(edges[b + 1]), 3),
                         "n": int(v.size), "median": round(med, 5), "se": round(se, 5)})  # fmt: skip
        meds = np.array([b["median"] for b in bins])
        ses = np.array([b["se"] for b in bins])
        chi = float(np.sum((meds / ses) ** 2))
        p = float(chi2.sf(chi, n_bins))
        max_abs = float(np.max(np.abs(meds)))
        out["directions"].append({
            "train_parity": train, "n_train": int(tr.sum()), "n_test": int(ok.sum()),
            "bins": bins, "max_abs_median": round(max_abs, 5), "chi2": round(chi, 2),
            "p": round(p, 5), "pass": bool(max_abs < max_abs_dex and p > min_p),
        })  # fmt: skip
    out["pass"] = all(x["pass"] for x in out["directions"])
    return out


def _by_tercile(resid, r, log_snr, sel, clusters, rng, n_boot) -> list[dict]:
    terc = np.nanquantile(log_snr[sel], [0, 1 / 3, 2 / 3, 1])
    rows = []
    for lo, hi in zip(terc[:-1], terc[1:], strict=True):
        s = sel & (log_snr >= lo) & (log_snr <= hi)
        rows.append({"log_snr": [round(float(lo), 3), round(float(hi), 3)],
                     **fit_beta(resid[s], r[s], clusters[s], rng, n_boot=n_boot)})  # fmt: skip
    return rows


def _within_2sigma(b: dict) -> bool:
    return bool(b.get("beta_se") and abs(b["beta"]) < 2 * b["beta_se"])


def run_gated_v2(
    field: dict,
    shifted: dict,
    rng: np.random.Generator,
    *,
    n_power: int = 20,
    n_inj: int = 10,
    n_boot: int = 1000,
) -> dict:
    """Plan 98's frozen order: C3' held-out calibration -> C0 -> C1 -> C4 -> C2 (overall AND per
    null S/N tercile) -> beta. ``field`` must carry ``log_n15``; the v2 calibration is used
    throughout, including inside the C0/C1 injections."""
    field = {**field, "calib": "v2"}
    out: dict = {"n_targets": len(field["flux_f"]), "calibration": "v2", "gates": {}}
    c3 = heldout_calibration_check(field, rng)
    c3_v1 = heldout_calibration_check({**field, "calib": "v1"}, rng)
    out["gates"]["C3prime_heldout"] = c3
    out["C3prime_plan97_form_for_comparison"] = {
        "pass": c3_v1["pass"],
        "max_abs_median": [x["max_abs_median"] for x in c3_v1["directions"]],
        "p": [x["p"] for x in c3_v1["directions"]],
    }
    if not c3["pass"]:
        out["stopped_at"] = "C3': flux scale cannot be calibrated to the precision the test needs"
        out["outcome"] = "stopped: calibration is the limiting systematic; no beta is quoted"
        return out
    gated = run_gated(field, {**shifted, "calib": "v2"}, rng, n_power=n_power, n_inj=n_inj,
                      n_boot=n_boot)  # fmt: skip
    out["gates"].update(gated["gates"])
    if "stopped_at" in gated:
        out["stopped_at"] = gated["stopped_at"]
        out["outcome"] = "ambiguous: " + gated["stopped_at"]
        return out
    d = _prepare(field)
    m = d["masks"]
    log_snr = np.log10(d["snr_g"])
    terc = _by_tercile(d["resid"], d["r_null"], log_snr, m["null"], d["clusters"], rng, n_boot)
    c2 = out["gates"]["C2_spectral_null"]
    c2["by_snr_tercile"] = terc
    c2["pass_overall"] = c2["pass"]
    c2["pass"] = bool(c2["pass_overall"] and all(_within_2sigma(t) for t in terc))
    out["C3_calibration_coef"] = gated["C3_calibration_coef"]
    out["samples"] = gated["samples"]
    out["primary"] = b = gated["primary"]
    out["secondary_primary_by_snr_tercile"] = _by_tercile(
        d["resid"], d["r_pred"], log_snr, m["primary"], d["clusters"], rng, n_boot
    )
    if not c2["pass"]:
        out["outcome"] = "ambiguous: C2 spectral null non-zero (overall or in an S/N tercile)"
    elif (b.get("beta_sigma") or 0) >= 3:
        out["outcome"] = ("blending supported after recalibration (a second attempt designed "
                          "after plan 97's samples were seen)")  # fmt: skip
    else:
        out["outcome"] = "blending not supported: beta consistent with 0"
        out["beta_upper_95"] = round(b["beta"] + 1.645 * b["beta_se"], 4)
    return out


def fetch_alfalfa() -> dict:  # pragma: no cover - network
    """ALFALFA alpha.100 (Haynes et al. 2018; VizieR J/ApJ/861/49/table2), HI centroids."""
    from .fashienv import _vizier

    v = _vizier(
        columns=["_RAJ2000", "_DEJ2000", "AGC", "Vhel", "W50", "HIflux", "e_HIflux", "SNR", "HI"]
    )
    v.ROW_LIMIT = -1
    t = v.get_catalogs("J/ApJ/861/49/table2")[0]
    out = {
        "ra": np.asarray(t["_RAJ2000"], float), "dec": np.asarray(t["_DEJ2000"], float),
        "v": np.asarray(t["Vhel"], float), "w50": np.asarray(t["W50"], float),
        "flux": np.asarray(t["HIflux"], float), "flux_err": np.asarray(t["e_HIflux"], float),
        "snr": np.asarray(t["SNR"], float), "code": np.asarray(t["HI"], int),
    }  # fmt: skip
    ok = np.isfinite(out["flux"]) & (out["flux"] > 0) & (out["flux"] < 999) & (out["flux_err"] > 0)
    return {k: val[ok] for k, val in out.items()}
