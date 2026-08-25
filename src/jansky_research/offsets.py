"""Radio--optical position offsets of AGN: ICRF3 (VLBI) vs Gaia DR3 (optical).

The VLBI radio position of an AGN marks the synchrotron self-absorbed jet base; the Gaia optical
position can be pulled toward optical jet/host structure. The two therefore do not coincide, and the
*normalised* offset $X=\\sqrt{(\\Delta\\alpha^*/\\sigma_\\alpha)^2+(\\Delta\\delta/\\sigma_\\delta)^2}$
shows a heavy tail far beyond the Rayleigh expectation for pure Gaussian astrometric noise --- a
well-established result (Mignard et al. 2016; Petrov & Kovalev 2017; Kovalev et al. 2017; Lindegren
et al. 2018). This module **reproduces** that excess tail with a small tested tool and builds a
reproducible offset catalogue.

Catalogue-only: ICRF3 (Charlot et al. 2020) cross-matched to Gaia DR3. Composes
:mod:`jansky_research.spectra` (``crossmatch``) and the CDS XMatch service. Pure NumPy + a synthetic
offline fixture for tests.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "alignment_stats",
    "directional_split",
    "fetch_icrf3_gaia",
    "fetch_mojave_jets",
    "inflation_sweep",
    "jet_axis_angles",
    "mahalanobis_offset",
    "match_jets",
    "noise_core_scale",
    "normalised_offset",
    "offset_statistics",
    "radio_optical_offset",
    "run",
    "synthetic_alignment",
    "synthetic_field",
]

DEG_TO_MAS = 3.6e6  # degrees -> milliarcseconds


def radio_optical_offset(
    ra_r: np.ndarray, dec_r: np.ndarray, ra_o: np.ndarray, dec_o: np.ndarray
) -> dict[str, np.ndarray]:
    r"""Radio$\to$optical offset components, total separation (mas), and position angle (deg E of N).

    $\Delta\alpha^*=(\alpha_o-\alpha_r)\cos\delta_r$ and $\Delta\delta=\delta_o-\delta_r$ (both in mas);
    the separation is $\sqrt{\Delta\alpha^{*2}+\Delta\delta^2}$ and the position angle
    $\mathrm{PA}=\mathrm{atan2}(\Delta\alpha^*,\Delta\delta)$ runs from the radio position toward the
    optical, measured East of North.
    """
    ra_r = np.asarray(ra_r, float)
    dec_r = np.asarray(dec_r, float)
    cosd = np.cos(np.radians(dec_r))
    dra = (np.asarray(ra_o, float) - ra_r) * cosd * DEG_TO_MAS
    ddec = (np.asarray(dec_o, float) - dec_r) * DEG_TO_MAS
    offset = np.hypot(dra, ddec)
    pa = np.degrees(np.arctan2(dra, ddec)) % 360.0
    return {"dra_mas": dra, "ddec_mas": ddec, "offset_mas": offset, "pa_deg": pa}


def normalised_offset(
    dra_mas: np.ndarray, ddec_mas: np.ndarray, sig_a_mas: np.ndarray, sig_d_mas: np.ndarray
) -> np.ndarray:
    r"""Significance of an offset: $X=\sqrt{(\Delta\alpha^*/\sigma_\alpha)^2+(\Delta\delta/\sigma_\delta)^2}$.

    The per-axis errors combine the radio and optical formal uncertainties in quadrature
    ($\sigma^2=\sigma_\mathrm{radio}^2+\sigma_\mathrm{Gaia}^2$, supplied by the caller). For pure
    Gaussian astrometric noise $X$ follows a 2-D Rayleigh, so $P(X>x)=e^{-x^2/2}$; a heavier tail is
    real structure.
    """
    a = np.asarray(dra_mas, float) / np.asarray(sig_a_mas, float)
    d = np.asarray(ddec_mas, float) / np.asarray(sig_d_mas, float)
    return np.sqrt(a**2 + d**2)


def mahalanobis_offset(
    dra_mas: np.ndarray,
    ddec_mas: np.ndarray,
    radio: dict,
    optical: dict,
) -> np.ndarray:
    r"""Correlation-aware normalised offset: the Mahalanobis distance of the 2-D offset.

    The naive per-axis $X$ treats the RA/Dec errors as independent, but both catalogues ship
    a per-source correlation (ICRF3 ``Corr``, Gaia ``RADEcor``), and the naive statistic's
    pure-noise tail is always *heavier* than Rayleigh when it is ignored -- the excess was
    being quoted against the wrong null. Here the combined covariance
    :math:`C = C_\mathrm{radio} + C_\mathrm{Gaia}` (each
    :math:`[[\sigma_a^2, \rho\sigma_a\sigma_d],[\rho\sigma_a\sigma_d, \sigma_d^2]]`) whitens
    the offset exactly, so for pure Gaussian noise :math:`X` is exactly Rayleigh and
    :math:`P(X>x)=e^{-x^2/2}` is the correct null by construction. Missing correlations are
    treated as zero.
    """
    dra = np.asarray(dra_mas, float)
    ddec = np.asarray(ddec_mas, float)
    out = np.full(dra.shape, np.nan)
    c11 = np.zeros(dra.shape)
    c22 = np.zeros(dra.shape)
    c12 = np.zeros(dra.shape)
    for cat in (radio, optical):
        sa = np.asarray(cat["e_a"], float)
        sd = np.asarray(cat["e_d"], float)
        rho = np.nan_to_num(np.asarray(cat.get("rho", np.zeros_like(sa)), float))
        c11 = c11 + sa**2
        c22 = c22 + sd**2
        c12 = c12 + rho * sa * sd
    det = c11 * c22 - c12**2
    ok = np.isfinite(det) & (det > 0)
    q = c22 * dra**2 - 2.0 * c12 * dra * ddec + c11 * ddec**2
    out[ok] = np.sqrt(q[ok] / det[ok])
    return out


def noise_core_scale(x_norm: np.ndarray, *, x_max: float = 1.5) -> dict:
    r"""Fit a truncated-Rayleigh scale to the noise core of the $X$ distribution.

    If the formal errors were underestimated by a common factor $f$, the noise core of $X$
    would follow a Rayleigh of scale $f$, not 1. Maximum-likelihood fit of the scale to the
    sources with $X<$ ``x_max`` (where structural offsets contribute least); structural
    contamination biases the fitted scale *high*, so it is an upper limit on any error
    underestimation. Returns the fitted scale and the core count.
    """
    x = np.asarray(x_norm, float)
    x = x[np.isfinite(x) & (x > 0) & (x < x_max)]
    n = int(x.size)
    if n < 10:
        return {"scale": float("nan"), "n_core": n, "x_max": x_max}
    s2 = float(np.mean(x**2)) / 2.0  # untruncated MLE as the start

    def _neg_ll(s: float) -> float:
        z = x / s
        c = x_max / s
        trunc = 1.0 - np.exp(-(c**2) / 2.0)  # Rayleigh CDF at the truncation
        return -float(np.sum(np.log(z / s) - z**2 / 2.0) - n * np.log(trunc))

    from scipy.optimize import minimize_scalar

    res = minimize_scalar(_neg_ll, bounds=(0.3, 3.0), method="bounded")
    return {"scale": round(float(res.x), 3), "n_core": n, "x_max": x_max, "s_raw": round(s2, 3)}


def inflation_sweep(
    x_norm: np.ndarray, *, x_cut: float = 3.0, factors: tuple[float, ...] = (1.2, 1.5, 2.0, 3.0)
) -> dict:
    """The excess under a uniform error inflation $f$: $X$ scales as $1/f$, so the inflated
    tail fraction is simply $P(X_\\mathrm{obs} > f\\,x_\\mathrm{cut})$. Retires 'robust to any
    plausible error inflation' by computing it. Keys are the factors as strings."""
    x = np.asarray(x_norm, float)
    x = x[np.isfinite(x)]
    rayleigh = float(np.exp(-(x_cut**2) / 2.0))
    out = {}
    for f in factors:
        frac = float((x > f * x_cut).mean())
        out[str(f)] = {
            "frac_pct": round(100.0 * frac, 2),
            "excess": round(frac / rayleigh, 1),
        }
    return out


def directional_split(
    offset_pa_deg: np.ndarray, jet_pa_deg: np.ndarray, *, axis_cut_deg: float = 45.0
) -> dict:
    r"""The directional test with the correct null: among AXIS-ALIGNED offsets, do more point
    downstream than upstream?

    Comparing frac(within 45° of downstream) = 50% to a "random 25%" double-counts the axis
    alignment the previous sentence has already established. Under a sign-symmetric
    axis-aligned distribution the two directions are equally likely, so the null is a fair
    coin over the aligned subset. Returns the split and its two-sided binomial p.
    """
    from scipy import stats as _stats

    down, axis = jet_axis_angles(np.asarray(offset_pa_deg, float), np.asarray(jet_pa_deg, float))
    aligned = axis < axis_cut_deg
    n_aligned = int(aligned.sum())
    n_down = int((down[aligned] < 90.0).sum())
    n_up = n_aligned - n_down
    p = float(_stats.binomtest(n_down, n_aligned, 0.5).pvalue) if n_aligned else float("nan")
    return {
        "axis_cut_deg": axis_cut_deg,
        "n_axis_aligned": n_aligned,
        "n_downstream": n_down,
        "n_upstream": n_up,
        "binom_p_two_sided": p,
        "frac_axis_lt45": (
            round(float(np.mean(axis < axis_cut_deg)), 3) if axis.size else float("nan")
        ),
    }


def offset_statistics(
    x_norm: np.ndarray, offset_mas: np.ndarray | None = None, *, x_cut: float = 3.0
) -> dict:
    r"""Summarise the offset population: the $X>$ ``x_cut`` fraction vs its Rayleigh expectation.

    Returns the count, the median raw offset (mas, if given), the fraction with $X>$ ``x_cut``, the
    Rayleigh expectation $e^{-x_\mathrm{cut}^2/2}$, and their ratio (the *excess* --- the reproduced
    result: AGN show many more significant offsets than Gaussian errors allow).
    """
    x = np.asarray(x_norm, float)
    x = x[np.isfinite(x)]
    n = int(x.size)
    frac = float((x > x_cut).mean()) if n else 0.0
    rayleigh = float(np.exp(-(x_cut**2) / 2.0))
    med = (
        float(np.nanmedian(np.asarray(offset_mas, float)))
        if offset_mas is not None and np.asarray(offset_mas).size
        else float("nan")
    )
    return {
        "n": n,
        "median_offset_mas": med,
        "x_cut": x_cut,
        "frac_x_gt_cut": frac,
        "rayleigh_expectation": rayleigh,
        "excess_ratio": frac / rayleigh if rayleigh > 0 else float("nan"),
    }


def synthetic_field(
    n_sources: int = 4000,
    *,
    structured_fraction: float = 0.15,
    sigma_mas: float = 0.4,
    struct_scale_mas: float = 3.0,
    rho: float = 0.0,
    seed: int = 0,
) -> tuple[dict, dict, np.ndarray]:
    """Synthetic ICRF3-like radio + Gaia-like optical positions with an injected structural-offset tail.

    Most sources have a pure-Gaussian radio$\\to$optical offset (per-axis error distributions
    scaled so their mean combined per-axis sigma is ``sigma_mas``); a ``structured_fraction``
    minority carry an additional real offset (exponential magnitude, random direction) standing
    in for optical jet/host structure. ``rho`` injects a common RA/Dec error correlation into
    the noise (and is returned per catalogue), so the correlation-aware statistic is testable.
    Returns ``(radio, optical, is_structured)`` with ra/dec (deg) and per-axis errors
    ``e_a``/``e_d`` (mas, on $\\alpha\\cos\\delta$ and $\\delta$).
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(0.0, 360.0, n_sources)
    dec = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, n_sources)))  # uniform on the sphere
    cosd = np.cos(np.radians(dec))
    k = sigma_mas / 0.4  # 0.4 mas is the mean combined sigma of the unscaled draws
    e_radio = rng.uniform(0.05, 0.5, n_sources) * k
    e_gaia = rng.uniform(0.05, 0.8, n_sources) * k
    sig = np.hypot(e_radio, e_gaia)
    # correlated noise offset (mas): z2 correlated with z1 at rho
    z1 = rng.standard_normal(n_sources)
    z2 = rho * z1 + np.sqrt(max(1.0 - rho**2, 0.0)) * rng.standard_normal(n_sources)
    dra = sig * z1
    ddec = sig * z2
    is_struct = rng.random(n_sources) < structured_fraction
    mag = rng.exponential(struct_scale_mas, n_sources)
    ang = rng.uniform(0.0, 2 * np.pi, n_sources)
    dra = np.where(is_struct, dra + mag * np.sin(ang), dra)
    ddec = np.where(is_struct, ddec + mag * np.cos(ang), ddec)
    rho_arr = np.full(n_sources, rho)
    radio = {"ra": ra, "dec": dec, "e_a": e_radio, "e_d": e_radio, "rho": rho_arr}
    optical = {
        "ra": ra + (dra / DEG_TO_MAS) / cosd,
        "dec": dec + ddec / DEG_TO_MAS,
        "e_a": e_gaia,
        "e_d": e_gaia,
        "rho": rho_arr,
    }
    return radio, optical, is_struct


def jet_axis_angles(
    offset_pa_deg: np.ndarray, jet_pa_deg: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    r"""Angles between a radio→optical offset PA and the parsec-scale jet PA.

    Returns ``(downstream_deg, axis_deg)``: the *downstream* angle wrapped to $[0,180]$ (0 = the offset
    points along the jet toward the approaching/downstream side, 180 = anti-jet/upstream), and the
    *jet-axis* angle $\min(\Delta,180-\Delta)\in[0,90]$ (0 = aligned with the jet axis, 90 =
    perpendicular). Both PAs are degrees East of North.
    """
    d = np.abs(
        ((np.asarray(offset_pa_deg, float) - np.asarray(jet_pa_deg, float) + 180.0) % 360.0) - 180.0
    )
    return d, np.minimum(d, 180.0 - d)


def alignment_stats(
    offset_pa_deg: np.ndarray, jet_pa_deg: np.ndarray, x_norm: np.ndarray, *, x_cut: float = 2.0
) -> dict:
    r"""Test whether radio→optical offsets align with the jet (the Kovalev/Petrov/Plavin result).

    The **full matched sample** is the primary test: a Kolmogorov--Smirnov comparison of the jet-axis
    angle (:func:`jet_axis_angles`, folded to $[0,90]$) against the uniform distribution expected if
    offsets were randomly oriented (median $45°$, fraction within $30°$ of the axis $=1/3$). The
    fraction of offsets pointing *downstream* (within $45°$ of the jet direction, random $=1/4$) isolates
    the downstream component. As a *qualitative* consistency check the same statistics are reported for
    the significant ($X>$ ``x_cut``) subset, which is expected to align more tightly (weak offsets are
    astrometric noise with random PA). Returns the counts and these statistics.
    """
    from scipy import stats as _stats

    opa = np.asarray(offset_pa_deg, float)
    jpa = np.asarray(jet_pa_deg, float)
    x = np.asarray(x_norm, float)
    good = np.isfinite(opa) & np.isfinite(jpa) & np.isfinite(x)
    opa, jpa, x = opa[good], jpa[good], x[good]
    down, axis = jet_axis_angles(opa, jpa)
    nan = float("nan")
    ks_p = float(_stats.kstest(axis / 90.0, "uniform").pvalue) if axis.size >= 5 else nan
    sig = x > x_cut
    return {
        "n_jet": int(axis.size),
        "median_axis_deg": float(np.median(axis)) if axis.size else nan,
        "frac_axis_lt30": float(np.mean(axis < 30.0)) if axis.size else nan,
        "frac_down_lt45": float(np.mean(down < 45.0)) if axis.size else nan,
        "ks_p": ks_p,
        "x_cut": x_cut,
        "n_jet_signif": int(sig.sum()),
        "median_axis_signif_deg": float(np.median(axis[sig])) if sig.any() else nan,
        "frac_axis_signif": float(np.mean(axis[sig] < 30.0)) if sig.any() else nan,
        "_axis_deg": axis,  # for the figure
    }


def synthetic_alignment(
    *,
    n: int = 420,
    aligned_fraction: float = 0.6,
    downstream_fraction: float = 0.8,
    jet_scatter_deg: float = 18.0,
    seed: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Synthetic matched offset/jet sample with an injected jet-aligned, mostly-downstream population.

    A fraction ``aligned_fraction`` of sources are *structural*: their offset points along the jet (a
    ``downstream_fraction`` majority downstream, the rest anti-jet) with ``jet_scatter_deg`` of scatter,
    and they carry a significant ``X``; the remainder are astrometric noise with a random offset PA and
    low ``X``. Returns ``(offset_pa_deg, jet_pa_deg, x_norm)`` so :func:`alignment_stats` recovers the
    injected alignment offline (no network).
    """
    rng = np.random.default_rng(seed)
    jet_pa = rng.uniform(0.0, 360.0, n)
    aligned = rng.random(n) < aligned_fraction
    downstream = rng.random(n) < downstream_fraction
    along = np.where(downstream, jet_pa, jet_pa + 180.0) + rng.normal(0.0, jet_scatter_deg, n)
    offset_pa = np.where(aligned, along % 360.0, rng.uniform(0.0, 360.0, n))
    x = np.where(aligned, rng.uniform(2.5, 8.0, n), rng.uniform(0.3, 2.0, n))
    return offset_pa, jet_pa, x


def fetch_mojave_jets() -> dict:  # pragma: no cover - network
    """Per-source mean innermost jet position angle from MOJAVE XVIII (Lister et al. 2021).

    VizieR ``J/ApJ/923/30/mojave18``: ``PA`` is the flux-weighted innermost jet PA measured from the
    core toward the approaching (downstream) side, deg East of North; ``delPA`` is its range (jet
    wobble). Returns sky positions (deg), ``jet_pa`` (deg), and ``delpa`` (deg).
    """
    from astroquery.vizier import Vizier

    v = Vizier(columns=["_RA", "_DE", "PA", "delPA"], row_limit=-1)
    t = v.get_catalogs("J/ApJ/923/30")["J/ApJ/923/30/mojave18"]
    return {
        "ra": np.asarray(t["_RA"], float),
        "dec": np.asarray(t["_DE"], float),
        "jet_pa": np.asarray(t["PA"], float),
        "delpa": np.asarray(t["delPA"], float),
    }


def match_jets(
    radio_ra_deg: np.ndarray, radio_dec_deg: np.ndarray, jets: dict, *, max_arcsec: float = 1.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Positionally match radio sources to the jet catalogue → ``(mask, jet_pa, delpa)``.

    Nearest-neighbour match within ``max_arcsec`` (same AGN). ``mask`` selects matched radio sources;
    ``jet_pa``/``delpa`` are aligned to the radio array order. Pure astropy (no network).
    """
    import astropy.units as _u
    from astropy.coordinates import SkyCoord

    rc = SkyCoord(
        np.asarray(radio_ra_deg, float) * _u.deg, np.asarray(radio_dec_deg, float) * _u.deg
    )
    jc = SkyCoord(jets["ra"] * _u.deg, jets["dec"] * _u.deg)
    idx, sep, _ = rc.match_to_catalog_sky(jc)
    mask = sep.arcsec < max_arcsec
    return mask, np.asarray(jets["jet_pa"], float)[idx], np.asarray(jets["delpa"], float)[idx]


def fetch_icrf3_gaia(*, max_arcsec: float = 0.5) -> tuple[dict, dict]:  # pragma: no cover - network
    """ICRF3 S/X (VizieR ``J/A+A/644/A159/table10``) cross-matched to Gaia DR3 via CDS X-Match.

    Returns ``(radio, optical)`` dicts with ra/dec (deg), per-axis errors ``e_a``/``e_d`` (mas,
    on $\\alpha\\cos\\delta$ and $\\delta$), the per-source RA/Dec error correlation ``rho``
    (ICRF3 ``Corr``; Gaia ``RADEcor``), and the ICRF designation ``name``. ICRF3 stores the RA
    error in **time-seconds** (so $\\sigma_{\\alpha^*}=e_\\mathrm{RA}\\times15000\\cos\\delta$
    mas) and the Dec error in **arcsec** ($\\times1000$ mas); Gaia ``e_RA_ICRS``/``e_DE_ICRS``
    are already mas. Keeps the nearest Gaia match per ICRF3 source within ``max_arcsec``.
    The chance-coincidence rate is measured, not assumed: the same query is repeated with
    every ICRF3 declination shifted by +0.2° (a decoy field), and the match count is returned
    as ``n_chance``.
    """
    import numpy as _np
    from astropy import units as _u
    from astropy.coordinates import SkyCoord
    from astropy.table import Table
    from astroquery.vizier import Vizier
    from astroquery.xmatch import XMatch

    v = Vizier(columns=["ICRF", "RAICRS", "DEICRS", "e_RAICRS", "e_DEICRS", "Corr"])
    v.ROW_LIMIT = -1
    icrf = v.get_catalogs("J/A+A/644/A159/table10")[0]
    # RAICRS/DEICRS are sexagesimal (RA in hours) -> parse to decimal degrees.
    coo = SkyCoord(icrf["RAICRS"], icrf["DEICRS"], unit=(_u.hourangle, _u.deg))
    ra = _np.asarray(coo.ra.deg, float)
    dec = _np.asarray(coo.dec.deg, float)
    keep = _np.isfinite(ra) & _np.isfinite(_np.asarray(icrf["e_RAICRS"], float))
    corr = _np.asarray(icrf["Corr"], float)
    if hasattr(icrf["Corr"], "mask"):
        corr = _np.where(_np.asarray(icrf["Corr"].mask), _np.nan, corr)
    names = _np.asarray([str(x).strip() for x in icrf["ICRF"]], dtype=object)
    t1 = Table(
        {
            "icrf_id": _np.arange(int(keep.sum())),
            "RAdeg": ra[keep],
            "DEdeg": dec[keep],
            "e_ra_s": _np.asarray(icrf["e_RAICRS"], float)[keep],
            "e_de_as": _np.asarray(icrf["e_DEICRS"], float)[keep],
            "rho_icrf": corr[keep],
        }
    )
    kept_names = names[keep]
    xm = XMatch.query(
        cat1=t1,
        cat2="vizier:I/355/gaiadr3",
        max_distance=max_arcsec * _u.arcsec,
        colRA1="RAdeg",
        colDec1="DEdeg",
    )
    # keep the nearest Gaia match per ICRF3 source
    xm.sort("angDist")
    _, first = _np.unique(_np.asarray(xm["icrf_id"]), return_index=True)
    xm = xm[first]
    decr = _np.asarray(xm["DEdeg"], float)
    cosd = _np.cos(_np.radians(decr))
    rho_gaia = _np.asarray(xm["RADEcor"], float)
    if hasattr(xm["RADEcor"], "mask"):
        rho_gaia = _np.where(_np.asarray(xm["RADEcor"].mask), _np.nan, rho_gaia)
    radio: dict = {
        "name": kept_names[_np.asarray(xm["icrf_id"], int)],
        "ra": _np.asarray(xm["RAdeg"], float),
        "dec": decr,
        "e_a": _np.asarray(xm["e_ra_s"], float) * 15000.0 * cosd,
        "e_d": _np.asarray(xm["e_de_as"], float) * 1000.0,
        "rho": _np.asarray(xm["rho_icrf"], float),
    }
    optical = {  # Gaia: RAdeg2/DEdeg2 positions; e_RAdeg already on alpha*cos(dec), mas
        "ra": _np.asarray(xm["RAdeg2"], float),
        "dec": _np.asarray(xm["DEdeg2"], float),
        "e_a": _np.asarray(xm["e_RAdeg"], float),
        "e_d": _np.asarray(xm["e_DEdeg"], float),
        "rho": rho_gaia,
    }
    # decoy field: same sources, dec + 0.2 deg -> every match is chance
    decoy = Table(
        {
            "icrf_id": t1["icrf_id"],
            "RAdeg": t1["RAdeg"],
            "DEdeg": _np.clip(_np.asarray(t1["DEdeg"], float) + 0.2, -89.9, 89.9),
        }
    )
    xd = XMatch.query(
        cat1=decoy,
        cat2="vizier:I/355/gaiadr3",
        max_distance=max_arcsec * _u.arcsec,
        colRA1="RAdeg",
        colDec1="DEdeg",
    )
    radio["n_chance"] = int(len(_np.unique(_np.asarray(xd["icrf_id"])))) if len(xd) else 0
    radio["n_icrf3"] = int(keep.sum())
    return radio, optical


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Full slice: synthesise (or fetch) ICRF3×Gaia, compute offsets, reproduce the excess tail."""
    from pathlib import Path

    if offline:
        radio, optical, truth = synthetic_field()
        source = "synthetic"
    else:  # pragma: no cover - network
        radio, optical = fetch_icrf3_gaia()
        truth = None
        source = "ICRF3 x Gaia DR3"

    off = radio_optical_offset(radio["ra"], radio["dec"], optical["ra"], optical["dec"])
    sig_a = np.hypot(radio["e_a"], optical["e_a"])
    sig_d = np.hypot(radio["e_d"], optical["e_d"])
    x_naive = normalised_offset(off["dra_mas"], off["ddec_mas"], sig_a, sig_d)
    # the headline statistic is correlation-aware, so the Rayleigh null is exact
    x = mahalanobis_offset(off["dra_mas"], off["ddec_mas"], radio, optical)
    stats = offset_statistics(x, off["offset_mas"])
    stats_naive = offset_statistics(x_naive)
    core = noise_core_scale(x)
    infl = inflation_sweep(x)

    # tail fraction by formal-error quartile: the tail is a property of the population
    # convolved with the error distribution, and this shows how strongly
    sig_mean = 0.5 * (sig_a + sig_d)
    qs = np.nanpercentile(sig_mean, [25, 50, 75])
    frac_by_err = []
    edges = [-np.inf, *qs, np.inf]
    for lo, hi in zip(edges[:-1], edges[1:], strict=False):
        m = (sig_mean > lo) & (sig_mean <= hi) & np.isfinite(x)
        frac_by_err.append(round(float((x[m] > 3.0).mean()), 3) if m.any() else None)

    # declination-band jackknife on the tail fraction: zonal systematics are the term the
    # binomial SE cannot see
    finite = np.isfinite(x)
    frac_all = float((x[finite] > 3.0).mean())
    bands = np.nanpercentile(radio["dec"][finite], [20, 40, 60, 80])
    band_id = np.digitize(radio["dec"][finite], bands)
    vals = []
    for b in range(5):
        keep = band_id != b
        vals.append(float((x[finite][keep] > 3.0).mean()))
    vals_a = np.asarray(vals)
    jk_se = float(np.sqrt(4.0 / 5.0 * np.sum((vals_a - vals_a.mean()) ** 2)))
    binom_se = float(np.sqrt(frac_all * (1.0 - frac_all) / max(int(finite.sum()), 1)))

    # raw-offset outliers: nothing structural produces a 500 mas photocentre shift
    om = np.asarray(off["offset_mas"], float)
    n_gt50 = int(np.sum(om > 50.0))
    n_gt100 = int(np.sum(om > 100.0))
    m50 = finite & (om <= 50.0)
    frac_le50 = float((x[m50] > 3.0).mean()) if m50.any() else float("nan")

    # jet-alignment test: does the offset DIRECTION point along the parsec-scale jet?
    delpa_matched: np.ndarray | None = None
    if offline:
        a_off_pa, a_jet_pa, a_x = synthetic_alignment()
        align = alignment_stats(a_off_pa, a_jet_pa, a_x)
        direction = directional_split(a_off_pa, a_jet_pa)
        align_goodjet = None
        jet_source = "synthetic"
        jet_pa_full = None
        mask = None
    else:  # pragma: no cover - network
        jets = fetch_mojave_jets()
        mask, jet_pa, delpa = match_jets(radio["ra"], radio["dec"], jets)
        align = alignment_stats(off["pa_deg"][mask], jet_pa[mask], x[mask])
        direction = directional_split(off["pa_deg"][mask], jet_pa[mask])
        # robustness: drop jets with poorly defined direction (wobble >= 45 deg)
        goodjet = mask & (delpa < 45.0)
        align_goodjet = alignment_stats(off["pa_deg"][goodjet], jet_pa[goodjet], x[goodjet])
        delpa_matched = delpa
        jet_source = "MOJAVE XVIII"
        jet_pa_full = jet_pa

    metrics = {
        "source": source,
        "n": stats["n"],
        "median_offset_mas": round(stats["median_offset_mas"], 3),
        "x_statistic": "mahalanobis (per-source RA/Dec correlations from both catalogues)",
        "x_cut": stats["x_cut"],
        "frac_x_gt3_pct": round(100.0 * stats["frac_x_gt_cut"], 2),
        "rayleigh_pct": round(100.0 * stats["rayleigh_expectation"], 2),
        "excess_ratio": round(stats["excess_ratio"], 1),
        "frac_x_gt3_pct_naive": round(100.0 * stats_naive["frac_x_gt_cut"], 2),
        "frac_x_gt3_binom_se_pct": round(100.0 * binom_se, 2),
        "frac_x_gt3_dec_jackknife_se_pct": round(100.0 * jk_se, 2),
        "frac_x_gt3_dec_bands_pct": [round(100.0 * v, 1) for v in vals],
        "noise_core_scale": core["scale"],
        "noise_core_n": core["n_core"],
        "noise_core_x_max": core["x_max"],
        "inflation_sweep": infl,
        "frac_x_gt3_by_err_quartile": frac_by_err,
        "n_offset_gt50mas": n_gt50,
        "n_offset_gt100mas": n_gt100,
        "frac_x_gt3_pct_le50mas": round(100.0 * frac_le50, 2),
        "jet_source": jet_source,
        "n_jet": align["n_jet"],
        "median_axis_deg": round(align["median_axis_deg"], 1),
        "frac_axis_lt30": round(align["frac_axis_lt30"], 3),
        "frac_axis_lt45": direction["frac_axis_lt45"],
        "frac_down_lt45": round(align["frac_down_lt45"], 3),
        "ks_p": align["ks_p"],
        "n_axis_aligned": direction["n_axis_aligned"],
        "n_downstream": direction["n_downstream"],
        "n_upstream": direction["n_upstream"],
        "downstream_binom_p": direction["binom_p_two_sided"],
        "n_jet_signif": align["n_jet_signif"],
        "median_axis_signif_deg": round(align["median_axis_signif_deg"], 1),
        "frac_axis_signif": round(align["frac_axis_signif"], 3),
    }
    for key in ("n_chance", "n_icrf3"):
        if key in radio:
            metrics[key] = radio[key]
    if align_goodjet is not None:  # pragma: no cover - real leg only
        metrics["goodjet_n"] = align_goodjet["n_jet"]
        metrics["goodjet_median_axis_deg"] = round(align_goodjet["median_axis_deg"], 1)
        metrics["goodjet_frac_axis_lt30"] = round(align_goodjet["frac_axis_lt30"], 3)
        metrics["goodjet_ks_p"] = align_goodjet["ks_p"]
    if (
        truth is not None
    ):  # synthetic: the excess tail should sit on the injected structured sources
        metrics["n_structured"] = int(truth.sum())
        metrics["frac_struct_in_tail"] = (
            round(float(truth[x > 3.0].mean()), 3) if (x > 3.0).any() else 0.0
        )

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "offsets_metrics.json")
    if radio.get("name") is not None:  # pragma: no cover - real leg only
        _write_sources(
            radio,
            optical,
            off,
            sig_a,
            sig_d,
            x,
            x_naive,
            mask,
            jet_pa_full,
            delpa_matched,
            op / "results" / "offsets_sources.csv",
        )
    _figure(x, align["_axis_deg"], core, op / "papers" / "offsets" / "figures")
    _write_macros(metrics, op / "papers" / "offsets" / "generated" / "macros.tex")
    return metrics


def _write_sources(
    radio, optical, off, sig_a, sig_d, x, x_naive, mask, jet_pa, delpa, path
) -> None:  # pragma: no cover - real leg only
    """The per-source matched catalogue -- the 'reproducible ICRF3xGaiaxMOJAVE matched
    catalogue' the abstract promises, previously claimed but never shipped."""
    import csv
    from pathlib import Path

    pt = Path(path)
    pt.parent.mkdir(parents=True, exist_ok=True)
    down = axis = None
    if jet_pa is not None:
        down, axis = jet_axis_angles(off["pa_deg"], jet_pa)
    with pt.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "icrf_name",
                "ra_deg",
                "dec_deg",
                "dra_mas",
                "ddec_mas",
                "offset_mas",
                "pa_deg",
                "sig_a_mas",
                "sig_d_mas",
                "rho_icrf3",
                "rho_gaia",
                "x_mahalanobis",
                "x_naive",
                "jet_pa_deg",
                "delpa_deg",
                "axis_angle_deg",
                "downstream_angle_deg",
            ]
        )
        for i in range(len(x)):
            has_jet = mask is not None and bool(mask[i])
            w.writerow(
                [
                    radio["name"][i],
                    f"{radio['ra'][i]:.8f}",
                    f"{radio['dec'][i]:.8f}",
                    f"{off['dra_mas'][i]:.4f}",
                    f"{off['ddec_mas'][i]:.4f}",
                    f"{off['offset_mas'][i]:.4f}",
                    f"{off['pa_deg'][i]:.2f}",
                    f"{sig_a[i]:.4f}",
                    f"{sig_d[i]:.4f}",
                    f"{radio['rho'][i]:.3f}" if np.isfinite(radio["rho"][i]) else "",
                    f"{optical['rho'][i]:.3f}" if np.isfinite(optical["rho"][i]) else "",
                    f"{x[i]:.3f}",
                    f"{x_naive[i]:.3f}",
                    f"{jet_pa[i]:.1f}" if has_jet else "",
                    f"{delpa[i]:.1f}" if has_jet and delpa is not None else "",
                    f"{axis[i]:.2f}" if has_jet and axis is not None else "",
                    f"{down[i]:.2f}" if has_jet and down is not None else "",
                ]
            )


def _figure(x_norm: np.ndarray, axis_deg: np.ndarray, core: dict, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.0))

    # Left: the offset-magnitude excess over Rayleigh (the existing result). The histogram is
    # normalised over the plotted range only, so the Rayleigh curve is scaled to the FITTED
    # noise core (amplitude = core count / Rayleigh mass below x_max, over in-range total):
    # comparing an all-mass Rayleigh to a range-truncated density misled at small X.
    x = np.asarray(x_norm, float)
    x = x[np.isfinite(x)]
    bins = np.linspace(0, 8, 40)
    n_inrange = int(np.sum(x <= 8.0))
    ax1.hist(x, bins=bins, density=True, color="0.6", label="ICRF3$\\times$Gaia")
    xs = 0.5 * (bins[:-1] + bins[1:])
    s = core.get("scale") or 1.0
    mass_core = 1.0 - np.exp(-(core.get("x_max", 1.5) ** 2) / (2.0 * s**2))
    amp = (core.get("n_core", 0) / mass_core) / max(n_inrange, 1)
    ax1.plot(
        xs,
        amp * (xs / s**2) * np.exp(-(xs**2) / (2.0 * s**2)),
        "r-",
        lw=1.2,
        label=f"Rayleigh, fitted noise core (scale {s:.2f})",
    )
    ax1.axvline(3.0, color="k", ls=":", lw=0.8, label="$X=3$")
    ax1.set(xlabel="normalised offset $X$", ylabel="density", title="Offset significance")
    ax1.legend(fontsize=8)
    ax1.set_yscale("log")

    # Right: the offset DIRECTION vs the jet axis (the new result) -- a peak at 0 = aligned
    a = np.asarray(axis_deg, float)
    a = a[np.isfinite(a)]
    abins = np.linspace(0, 90, 19)
    ax2.hist(a, bins=abins, density=True, color="C0", label="offset vs jet")
    ax2.axhline(1.0 / 90.0, color="r", ls="-", lw=1.2, label="uniform (random)")
    ax2.set(
        xlabel="jet-axis angle (deg)",
        ylabel="density",
        title="Offset direction vs parsec-scale jet",
    )
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "xnorm.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _texp(p: float) -> str:
        """Format a (tiny) p-value as a LaTeX exponent, e.g. ``3.2\\times10^{-22}``."""
        if not np.isfinite(p) or p <= 0:
            return "--"
        e = int(np.floor(np.log10(p)))
        mant = p / 10.0**e
        return rf"{mant:.1f}\times10^{{{e}}}" if e < -2 else f"{p:.3f}"

    def _fmt(key: str) -> str:
        v = m.get(key)
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "--"
        return str(v)

    infl = m.get("inflation_sweep") or {}

    def _infl(f: str, field: str) -> str:
        v = (infl.get(f) or {}).get(field)
        return "--" if v is None else str(v)

    lines = [
        "% Auto-generated by jansky_research.offsets._write_macros — do not edit by hand.",
        rf"\newcommand{{\offSource}}{{{m['source']}}}",
        rf"\newcommand{{\offN}}{{{m['n']}}}",
        rf"\newcommand{{\offMedian}}{{{m['median_offset_mas']}}}",
        rf"\newcommand{{\offFracTail}}{{{m['frac_x_gt3_pct']}}}",
        rf"\newcommand{{\offRayleigh}}{{{m['rayleigh_pct']}}}",
        rf"\newcommand{{\offExcess}}{{{m['excess_ratio']}}}",
        # the naive (correlation-blind) statistic, kept as the contrast with the old headline
        rf"\newcommand{{\offFracTailNaive}}{{{_fmt('frac_x_gt3_pct_naive')}}}",
        # uncertainties: binomial and the declination-band jackknife that sees zonal terms
        rf"\newcommand{{\offFracTailSe}}{{{_fmt('frac_x_gt3_binom_se_pct')}}}",
        rf"\newcommand{{\offFracTailJkSe}}{{{_fmt('frac_x_gt3_dec_jackknife_se_pct')}}}",
        # the measured error model: noise-core Rayleigh scale + the inflation sweep
        rf"\newcommand{{\offCoreScale}}{{{_fmt('noise_core_scale')}}}",
        rf"\newcommand{{\offCoreN}}{{{_fmt('noise_core_n')}}}",
        rf"\newcommand{{\offInflOneFive}}{{{_infl('1.5', 'excess')}}}",
        rf"\newcommand{{\offInflTwo}}{{{_infl('2.0', 'excess')}}}",
        rf"\newcommand{{\offInflThree}}{{{_infl('3.0', 'excess')}}}",
        # chance matches and raw-offset outliers
        rf"\newcommand{{\offNchance}}{{{_fmt('n_chance')}}}",
        rf"\newcommand{{\offNicrf}}{{{_fmt('n_icrf3')}}}",
        rf"\newcommand{{\offNgtFifty}}{{{_fmt('n_offset_gt50mas')}}}",
        rf"\newcommand{{\offFracTailLeFifty}}{{{_fmt('frac_x_gt3_pct_le50mas')}}}",
        rf"\newcommand{{\offJetSource}}{{{m['jet_source']}}}",
        rf"\newcommand{{\offJetN}}{{{m['n_jet']}}}",
        rf"\newcommand{{\offJetMedAxis}}{{{m['median_axis_deg']}}}",
        rf"\newcommand{{\offJetFracAxis}}{{{round(100.0 * m['frac_axis_lt30'])}}}",
        rf"\newcommand{{\offJetFracDown}}{{{round(100.0 * m['frac_down_lt45'])}}}",
        rf"\newcommand{{\offJetKsP}}{{{_texp(m['ks_p'])}}}",
        # the directional test with the correct (sign-symmetric) null
        rf"\newcommand{{\offNaligned}}{{{_fmt('n_axis_aligned')}}}",
        rf"\newcommand{{\offNdown}}{{{_fmt('n_downstream')}}}",
        rf"\newcommand{{\offNup}}{{{_fmt('n_upstream')}}}",
        rf"\newcommand{{\offDownP}}{{{_texp(m.get('downstream_binom_p', float('nan')))}}}",
        rf"\newcommand{{\offJetNsig}}{{{m['n_jet_signif']}}}",
        rf"\newcommand{{\offJetMedAxisSig}}{{{m['median_axis_signif_deg']}}}",
        rf"\newcommand{{\offJetFracAxisSig}}{{{round(100.0 * m['frac_axis_signif'])}}}",
        # the delPA<45 robustness variant, previously only in the findings file
        rf"\newcommand{{\offGoodjetN}}{{{_fmt('goodjet_n')}}}",
        rf"\newcommand{{\offGoodjetMedAxis}}{{{_fmt('goodjet_median_axis_deg')}}}",
        rf"\newcommand{{\offGoodjetKsP}}{{{_texp(m.get('goodjet_ks_p', float('nan')))}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and
    # would otherwise blank the other mode's macros with '--'. `make figures`
    # runs every slice offline in the repo root, so without this an offline
    # rebuild silently empties this paper. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Radio-optical offsets of AGN (ICRF3 x Gaia DR3).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
