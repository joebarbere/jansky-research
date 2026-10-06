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
    obs_a = sa * 10 ** (offset + rng.normal(0, scatter_dex, n_targets))
    obs_f = sf * 10 ** rng.normal(0, scatter_dex, n_targets)
    return {
        "ra": ra, "dec": dec, "v": v, "w50": w50, "flux_true": flux, "snr": snr,
        "flux_a": obs_a, "flux_f": obs_f, "nb": nb,
        "flux_err_f": flux / snr,
        "sig_a_dex": np.full(n_targets, scatter_dex), "sig_f_dex": np.full(n_targets, scatter_dex),
    }  # fmt: skip


def analyse(field: dict, rng: np.random.Generator, *, n_boot: int = 500) -> dict:
    """The full frozen analysis on a matched catalogue: C3 calibration on isolated targets,
    beta on the primary sample, beta_null on the C2 sample."""
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
    design = calibration_design(
        np.log10(s_g / field["flux_err_f"]), np.log10(field["w50"]), np.log10(s_g), field["dec"]
    )
    coef = fit_calibration(y[masks["isolated"]], design[masks["isolated"]])
    resid = y - design @ coef
    clusters = target_clusters(field["ra"], field["dec"])
    r_pred = predicted_log_ratio(s_c, nb)
    r_null = predicted_log_ratio(s_c, nb, use=null_pairs(nb), force_overlap=True)
    p, q = masks["primary"], masks["null"]
    return {
        "n_targets": n,
        "n_isolated": int(masks["isolated"].sum()),
        "n_primary": int(p.sum()),
        "n_null": int(q.sum()),
        "calibration_coef": [round(float(c), 5) for c in coef],
        "primary": fit_beta(resid[p], r_pred[p], clusters[p], rng, n_boot=n_boot),
        "null": fit_beta(resid[q], r_null[q], clusters[q], rng, n_boot=n_boot),
    }


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
