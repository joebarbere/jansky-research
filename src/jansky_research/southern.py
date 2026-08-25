"""Southern peaked-spectrum (GPS/CSS) selection via GLEAM-X + RACS multi-band spectral curvature.

The northern :mod:`jansky_research.peaked` slice could only *bound* the low-frequency spectral index
(TGSS is shallow, so 150 MHz is an upper limit). In the south, GLEAM-X DR2 (Ross et al. 2024) measures
each source in **19 clean in-band sub-bands over 76--227 MHz**, so the low-frequency shape is *measured*.
Adding the three RACS bands (887.5, 1367.5, 1655.5 MHz) gives up to 22 flux points over a factor of
~22 in frequency --- enough to fit a real log-parabola SED and **measure the turnover frequency**
$\\nu_\\mathrm{pk}$, the upgrade this slice is built for.

This composes :mod:`jansky_research.spectra` (two-point index, cross-match) and reuses
:func:`jansky_research.peaked.classify_sed`, generalising ``peaked.peak_frequency`` to a weighted
N-point fit. Pure NumPy + a synthetic offline fixture for tests.
"""

from __future__ import annotations

import numpy as np

from .spectra import spectral_index

__all__ = [
    "GLEAMX_NU_GHZ",
    "RACS_NU_GHZ",
    "USS_THRESHOLD",
    "classify_curved",
    "fetch_gleamx",
    "fetch_racs_bands",
    "find_peaked_south",
    "fit_log_parabola",
    "run",
    "synthetic_field",
    "validate_callingham",
]

# GLEAM-X DR2 sub-band centres (MHz): the real ``Fint###`` columns of VizieR VIII/113/catalog2
# (the 092 band is omitted — its integrated-flux column is mislabelled ``Fpint092`` in VizieR).
_GLEAM_BANDS_MHZ = (
    76,
    84,
    99,
    107,
    115,
    122,
    130,
    143,
    151,
    158,
    166,
    174,
    181,
    189,
    197,
    204,
    212,
    220,
    227,
)
GLEAMX_NU_GHZ = np.array(_GLEAM_BANDS_MHZ) / 1000.0
# RACS band reference frequencies (GHz): RACS-low, RACS-mid, RACS-high.
RACS_NU_GHZ = np.array([0.8875, 1.3675, 1.6555])
# Ultra-steep-spectrum threshold (candidate high-z radio galaxies).
USS_THRESHOLD = -1.2
# GLEAM-X integrated-to-peak ratio above which a source counts as extended (compactness cut).
COMPACT_RATIO_MAX = 1.2

# RadioSED II (Kerrison et al. 2025) peaked-source surface density over Stripe 82, deg^-2 --
# the literature comparison point. The comparison is an indication, not a controlled ratio:
# the two samples differ in survey depth, frequency coverage, and selection method.
KERRISON_DENSITY_PER_DEG2 = 1.2

# VizieR catalogue IDs and the integrated-flux column per RACS band (RACS-low names it Fpk; mid/high
# use Fpeak; all three carry Ftot which we use for an integrated-flux SED consistent with GLEAM-X Fint).
VIZIER_GLEAMX = "VIII/113/catalog2"  # GLEAM-X DR2 (Ross+2024) — 20 sub-bands, fluxes in Jy
VIZIER_RACS = {
    "low": "J/other/PASA/38.58/galcut",  # RACS-low DR1 (Hale+2021), mJy
    "mid": "J/other/PASA/41.3/sourcesm",  # RACS-mid (Duchesne+2024), mJy
    "high": "J/other/PASA/42.38/sourcesh",  # RACS-high (Duchesne+2025), mJy
}


def fit_log_parabola(nu_ghz: np.ndarray, flux: np.ndarray, eflux: np.ndarray | None = None) -> dict:
    r"""Weighted log-parabola fit $\log_{10}S = a\,x^2 + b\,x + c$ with $x=\log_{10}\nu$.

    Generalises ``peaked.peak_frequency`` to N points with errors. The extremum is at
    $x_\mathrm{pk}=-b/2a$ and is a genuine (concave) peak when $a<0$; the turnover is
    $\nu_\mathrm{pk}=10^{x_\mathrm{pk}}$. ``is_peaked`` requires a concave fit whose turnover lies
    *within the sampled band* (a measured turnover, not an extrapolation). Returns a dict with
    ``nu_pk_ghz``, ``a`` (curvature), ``b``, ``c``, ``is_peaked``, ``chi2_red``, and ``n_points``
    (finite positive points used). Needs $\geq 4$ such points.
    """
    nu = np.asarray(nu_ghz, float)
    s = np.asarray(flux, float)
    good = np.isfinite(nu) & np.isfinite(s) & (s > 0) & (nu > 0)
    if eflux is not None:
        # With real per-band errors available, a sub-3-sigma point is not a detection: keeping
        # only its positive noise excursions (the log fit cannot take negatives) censors the
        # noise distribution asymmetrically and manufactures a rising low-frequency side --
        # the round-5 referee's second blocker. Require a detection to enter the fit.
        e = np.asarray(eflux, float)
        good &= np.isfinite(e) & (e > 0) & (s > 3.0 * e)
    nan = float("nan")
    out = {
        "nu_pk_ghz": nan,
        "a": nan,
        "b": nan,
        "c": nan,
        "is_peaked": False,
        "chi2_red": nan,
        "n_points": int(good.sum()),
    }
    if good.sum() < 4:
        return out
    x = np.log10(nu[good])
    y = np.log10(s[good])
    if eflux is not None:
        e = np.asarray(eflux, float)[good]
        sigma_y = np.where((e > 0) & np.isfinite(e), e / (s[good] * np.log(10.0)), np.nan)
        w = np.where(np.isfinite(sigma_y) & (sigma_y > 0), 1.0 / sigma_y, 1.0)
    else:
        w = np.ones_like(x)
        sigma_y = np.ones_like(x)
    a, b, c = np.polyfit(x, y, 2, w=w)
    if a == 0.0 or not np.isfinite(a):
        return out
    x_pk = -b / (2.0 * a)
    model = a * x**2 + b * x + c
    dof = max(x.size - 3, 1)
    chi2 = float(
        np.sum(((y - model) / np.where(np.isfinite(sigma_y) & (sigma_y > 0), sigma_y, 1.0)) ** 2)
    )
    in_band = float(x.min()) <= x_pk <= float(x.max())
    out.update(
        nu_pk_ghz=float(10.0**x_pk) if abs(x_pk) < 6.0 else nan,
        a=float(a),
        b=float(b),
        c=float(c),
        is_peaked=bool(a < 0.0 and in_band),
        chi2_red=chi2 / dof,
    )
    return out


def fitted_index(nu_ghz: np.ndarray, flux: np.ndarray, eflux: np.ndarray) -> tuple[float, float]:
    r"""Weighted in-band spectral index with its standard error, over $\ge3\sigma$ detections.

    A weighted linear fit of $\log_{10}S$ on $\log_{10}\nu$. This replaces the two-point
    76-vs-227 MHz ratio the rising-side gate used to rest on: 76 MHz is the least sensitive
    sub-band of the nineteen, so that gate was anchored on the noisiest measurement and flipped
    under a one-sub-band change of anchor (the round-5 referee decoded the committed figure's
    own candidates to show it). Returns ``(alpha, alpha_err)``; NaNs with $<3$ usable points.
    """
    nu = np.asarray(nu_ghz, float)
    s = np.asarray(flux, float)
    e = np.asarray(eflux, float)
    good = np.isfinite(nu) & np.isfinite(s) & np.isfinite(e) & (s > 3.0 * e) & (e > 0) & (nu > 0)
    if good.sum() < 3:
        return float("nan"), float("nan")
    x = np.log10(nu[good])
    y = np.log10(s[good])
    sy = e[good] / (s[good] * np.log(10.0))
    w = 1.0 / sy**2
    sw, sx, sxx = float(w.sum()), float((w * x).sum()), float((w * x * x).sum())
    sxy, sy_ = float((w * x * y).sum()), float((w * y).sum())
    den = sw * sxx - sx * sx
    if den <= 0:
        return float("nan"), float("nan")
    alpha = (sw * sxy - sx * sy_) / den
    return float(alpha), float(np.sqrt(sw / den))


def classify_curved(
    fit: dict,
    alpha_lo: float,
    alpha_hi: float,
    *,
    rise_min: float = -0.1,
    alpha_lo_err: float = 0.0,
) -> str:
    """Classify a southern SED: peaked / uss / steep / flat / inverted (reuses ``peaked.classify_sed``).

    ``peaked`` requires *both* a concave in-band log-parabola turnover **and** an optically-thick
    rising low-frequency side --- $\\alpha_\\mathrm{lo}-\\sigma_{\\alpha}>$ ``rise_min``, i.e.\\ the
    rising/flat side must be *significant*, not merely positive: with per-band errors a noise
    fluctuation on the faint low-frequency points can fake a positive two-point index (the
    round-5 blocker). Otherwise ``uss`` when both indices are ultra-steep ($<$ ``USS_THRESHOLD``;
    candidate high-z radio galaxy), else the coarse two-index class.
    """
    from .peaked import classify_sed

    if not (np.isfinite(alpha_lo) and np.isfinite(alpha_hi)):
        return "nan"
    err = alpha_lo_err if np.isfinite(alpha_lo_err) else 0.0
    if fit.get("is_peaked") and (alpha_lo - err) > rise_min:
        return "peaked"
    if alpha_lo < USS_THRESHOLD and alpha_hi < USS_THRESHOLD:
        return "uss"
    return classify_sed(alpha_lo, alpha_hi)


def find_peaked_south(
    gleamx: dict[str, np.ndarray], racs: dict[str, np.ndarray], *, radius_arcsec: float = 25.0
) -> dict[str, np.ndarray]:
    """Cross-match GLEAM-X (multi-band) × RACS (3-band), fit each SED, classify. Reuses ``spectra``.

    ``gleamx`` carries ``ra``/``dec`` and ``flux``/``eflux`` arrays of shape ``(n, len(GLEAMX_NU_GHZ))``;
    ``racs`` likewise with shape ``(m, 3)``. For each GLEAM-X source with a RACS match the full SED is
    fit with :func:`fit_log_parabola`; ``alpha_lo`` is the measured GLEAM-X in-band index and
    ``alpha_hi`` the RACS-low$\\to$high index. Returns per-matched-source arrays incl. ``nu_pk_ghz``,
    ``cls``, ``is_peaked``, ``is_uss``.
    """
    from .spectra import crossmatch

    gra, gdec = np.asarray(gleamx["ra"], float), np.asarray(gleamx["dec"], float)
    ig, ir, _ = crossmatch(gra, gdec, racs["ra"], racs["dec"], radius_arcsec)
    if ig.size == 0:
        return {k: np.array([]) for k in ("ra", "dec", "nu_pk_ghz", "cls", "is_peaked", "is_uss")}
    gflux = np.asarray(gleamx["flux"], float)
    geflux = np.asarray(gleamx["eflux"], float)
    rflux = np.asarray(racs["flux"], float)
    reflux = np.asarray(racs["eflux"], float)
    nu_all = np.concatenate([GLEAMX_NU_GHZ, RACS_NU_GHZ])
    # GLEAM compactness flag (resolution-artefact guard); if absent, treat all as compact.
    compact = gleamx.get("compact")
    compact_arr = np.ones(gra.size, bool) if compact is None else np.asarray(compact, bool)

    nu_pk, cls, is_pk, is_uss, a_lo, a_lo_e, a_hi, chi2, npts = [], [], [], [], [], [], [], [], []
    n_naive = n_rising = 0
    for k in range(ig.size):
        gi, rj = int(ig[k]), int(ir[k])
        flux = np.concatenate([gflux[gi], rflux[rj]])
        eflux = np.concatenate([geflux[gi], reflux[rj]])
        fit = fit_log_parabola(nu_all, flux, eflux)
        # the rising-side gate uses the FITTED in-band index with its error, not the two-point
        # 76-vs-227 ratio (which was anchored on the noisiest sub-band); alpha_hi stays the
        # RACS low->high two-point index (all three RACS bands are high-S/N)
        alpha_lo, alpha_lo_err = fitted_index(GLEAMX_NU_GHZ, gflux[gi], geflux[gi])
        alpha_hi, _ = spectral_index(rflux[rj, 0], RACS_NU_GHZ[0], rflux[rj, -1], RACS_NU_GHZ[-1])
        n_naive += int(bool(fit["is_peaked"]))
        c = classify_curved(fit, float(alpha_lo), float(alpha_hi), alpha_lo_err=float(alpha_lo_err))
        n_rising += int(c == "peaked")
        if not compact_arr[gi] and c in ("peaked", "inverted"):
            c = "extended"  # GLEAM-extended -> RACS flux loss fakes a turnover; reject
        nu_pk.append(fit["nu_pk_ghz"])
        cls.append(c)
        is_pk.append(c == "peaked")
        is_uss.append(c == "uss")
        a_lo.append(float(alpha_lo))
        a_lo_e.append(float(alpha_lo_err))
        a_hi.append(float(alpha_hi))
        chi2.append(float(fit["chi2_red"]))
        npts.append(int(fit["n_points"]))
    return {
        "ra": gra[ig],
        "dec": gdec[ig],
        "nu_pk_ghz": np.asarray(nu_pk),
        "alpha_lo": np.asarray(a_lo),
        "alpha_lo_err": np.asarray(a_lo_e),
        "alpha_hi": np.asarray(a_hi),
        "chi2_red": np.asarray(chi2),
        "n_fit_points": np.asarray(npts),
        "cls": np.asarray(cls, dtype=object).astype(str),
        "is_peaked": np.asarray(is_pk, bool),
        "is_uss": np.asarray(is_uss, bool),
        # the selection cascade, in the order the code applies it (naive concave-in-band ->
        # significant rising side -> compactness), plus the crossmatch's many-to-one count:
        # nearest-neighbour matching is not bijective, and a RACS source claimed by several
        # GLEAM-X sources is exactly the crowded/extended case the compactness cut hunts
        "n_peaked_naive": int(n_naive),
        "n_peaked_after_rising": int(n_rising),
        "n_racs_matched_multiply": int(ir.size - np.unique(ir).size),
    }


def synthetic_field(
    n_sources: int = 1200,
    *,
    peaked_fraction: float = 0.05,
    uss_fraction: float = 0.05,
    flatten_fraction: float = 0.08,
    extended_fraction: float = 0.05,
    rel_err: float = 0.08,
    noise_floor_mjy: float = 6.0,
    seed: int = 0,
) -> tuple[dict, dict, np.ndarray, np.ndarray]:
    """Synthetic GLEAM-X(19-band)+RACS(3-band) catalogues with injected classes and real-shaped noise.

    Beyond the injected peaked/uss/steep/flat mix, the fixture now carries what a validation
    must be able to fail on (round-5 referee):

    - a **noise floor** per GLEAM sub-band (strongest at 76 MHz, falling with frequency, like
      the real catalogue's ~7--8 mJy floors) added to the fluxes and recorded in ``eflux`` --
      faint sources' low-frequency points are noise, exactly the population the real gate
      selects on, and negatives are kept rather than clipped;
    - a **low-frequency-flattening contaminant** (a log-parabola peaking *below* the band):
      mild flattening that a merely-positive two-point index would select but a significant
      rising side must reject;
    - an injected **extended** population (``compact=False``, with RACS fluxes halved by
      resolution loss) so the compactness cut has something to catch — it was previously
      dead code offline.

    Injected turnovers span 0.1--0.7 GHz, covering both the GLEAM band and the 0.23--0.89 GHz
    sampling gap (previously all injections sat inside the gap).
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(20.0, 25.0, n_sources)
    dec = rng.uniform(-40.0, -35.0, n_sources)
    s_ref = 10.0 ** rng.uniform(1.0, 2.5, n_sources)  # ~10-300 mJy at 200 MHz
    nu_ref = 0.2  # GHz
    u = rng.random(n_sources)
    is_peaked = u < peaked_fraction
    is_uss = (~is_peaked) & (u < peaked_fraction + uss_fraction)
    is_flatten = (~is_peaked) & (~is_uss) & (u < peaked_fraction + uss_fraction + flatten_fraction)
    alpha = rng.uniform(-0.9, -0.5, n_sources)  # ordinary steep/flat
    alpha[is_uss] = rng.uniform(-1.6, -1.3, int(is_uss.sum()))  # ultra-steep
    nu_pk = rng.uniform(0.1, 0.7, n_sources)  # injected turnover (GHz) for peaked sources
    nu_pk_flat = rng.uniform(0.02, 0.05, n_sources)  # below-band peak -> in-band flattening
    curv = rng.uniform(0.6, 1.4, n_sources)
    is_extended = rng.random(n_sources) < extended_fraction

    def sed(nu):
        out = np.empty((n_sources, nu.size))
        for i in range(n_sources):
            if is_peaked[i]:
                out[i] = 10.0 ** (np.log10(s_ref[i]) - curv[i] * (np.log10(nu / nu_pk[i])) ** 2)
            elif is_flatten[i]:
                out[i] = 10.0 ** (np.log10(s_ref[i]) - 0.4 * (np.log10(nu / nu_pk_flat[i])) ** 2)
            else:
                out[i] = s_ref[i] * (nu / nu_ref) ** alpha[i]
        return out

    # GLEAM noise floor falls with frequency (76 MHz is the least sensitive sub-band)
    gfloor = noise_floor_mjy * (0.076 / GLEAMX_NU_GHZ)
    rfloor = np.full(RACS_NU_GHZ.size, 0.3)
    gtrue = sed(GLEAMX_NU_GHZ)
    rtrue = sed(RACS_NU_GHZ)
    rtrue[is_extended] *= 0.5  # resolution loss at RACS fakes a high-frequency turnover
    gflux = gtrue * rng.normal(1.0, rel_err, gtrue.shape) + rng.normal(
        0.0, gfloor[None, :], gtrue.shape
    )
    rflux = rtrue * rng.normal(1.0, rel_err, rtrue.shape) + rng.normal(
        0.0, rfloor[None, :], rtrue.shape
    )
    geflux = np.hypot(gfloor[None, :], rel_err * np.abs(gflux))
    reflux = np.hypot(rfloor[None, :], rel_err * np.abs(rflux))
    jit = lambda: rng.normal(0.0, 2.0 / 3600.0, n_sources)  # noqa: E731  (~2" jitter)
    gleamx = {
        "ra": ra + jit(),
        "dec": dec + jit(),
        "flux": gflux,
        "eflux": geflux,
        "compact": ~is_extended,
    }
    racs = {"ra": ra + jit(), "dec": dec + jit(), "flux": rflux, "eflux": reflux}
    return gleamx, racs, is_peaked, is_uss


def fetch_gleamx(center, radius_deg: float) -> dict:  # pragma: no cover - network
    """Cone-search GLEAM-X DR2 (VizieR ``VIII/113``); 19 in-band sub-band fluxes, **Jy → mJy**.

    GLEAM-X DR2 stores integrated fluxes in **Jy** (RACS is in mJy) — converted here ×1000 to a
    common mJy scale. The catalogue's own per-sub-band uncertainties (``e_Fint076..e_Fint227``)
    are fetched alongside — an earlier version claimed VizieR exposed no per-band errors and
    substituted a 10% proportional error, which is inverted relative to the true noise *floor*
    (~7--8 mJy at 76 MHz regardless of source brightness) and understated faint sources' errors
    by an order of magnitude, exactly where the peaked-candidate selection operates (the round-5
    referee's blocker). Fluxes are returned as measured, negatives included: the log-space fit
    can only use positive points, but which points those are is decided downstream against the
    real errors, not by a silent clip with a 10% error attached.
    """
    import numpy as _np
    from astropy import units as _u
    from astroquery.vizier import Vizier

    nb = len(_GLEAM_BANDS_MHZ)
    cols = [
        "RAJ2000",
        "DEJ2000",
        "Fintwide",
        "Fpwide",
        *[f"Fint{b:03d}" for b in _GLEAM_BANDS_MHZ],
        *[f"e_Fint{b:03d}" for b in _GLEAM_BANDS_MHZ],
    ]
    v = Vizier(columns=cols)
    v.ROW_LIMIT = -1
    res = v.query_region(center, radius=radius_deg * _u.deg, catalog=VIZIER_GLEAMX)
    if not res:
        z = _np.zeros((0, nb))
        return {
            "ra": _np.array([]),
            "dec": _np.array([]),
            "flux": z,
            "eflux": z,
            "compact": _np.array([], bool),
        }
    t = res[0]
    flux = _np.column_stack([_np.asarray(t[f"Fint{b:03d}"], float) for b in _GLEAM_BANDS_MHZ]) * 1e3
    eflux_cat = (
        _np.column_stack([_np.asarray(t[f"e_Fint{b:03d}"], float) for b in _GLEAM_BANDS_MHZ]) * 1e3
    )
    # compactness at GLEAM's ~2' resolution: integrated/peak ~ 1 for point sources. Non-compact (>1.2)
    # sources lose flux at RACS's 15-25" resolution, faking a high-frequency turnover (resolution
    # artefact) -- exclude them, the southern analogue of the peaked slice's alpha_high floor.
    fint_w = _np.asarray(t["Fintwide"], float)
    fp_w = _np.asarray(t["Fpwide"], float)
    with _np.errstate(divide="ignore", invalid="ignore"):
        ratio = fint_w / fp_w
    compact = _np.isfinite(ratio) & (ratio < COMPACT_RATIO_MAX)
    return {
        "ra": _np.asarray(t["RAJ2000"], float),
        "dec": _np.asarray(t["DEJ2000"], float),
        "flux": flux,
        "eflux": eflux_cat,
        "compact": compact,
    }


def fetch_racs_bands(
    center, radius_deg: float, *, match_arcsec: float = 15.0, rel_err: float = 0.1
) -> dict:  # pragma: no cover - network
    """Fetch RACS-low/mid/high integrated fluxes (``Ftot``, mJy) and assemble a unified ``(m, 3)`` table.

    Anchors on RACS-low and cross-matches RACS-mid and RACS-high onto its positions (NaN where a band
    has no match), giving a per-source 3-band flux array in ``RACS_NU_GHZ`` order.
    """
    import numpy as _np
    from astropy import units as _u
    from astroquery.vizier import Vizier

    from .spectra import crossmatch

    def _band(cid):
        v = Vizier(columns=["RAJ2000", "DEJ2000", "Ftot"])
        v.ROW_LIMIT = -1
        res = v.query_region(center, radius=radius_deg * _u.deg, catalog=cid)
        if not res:
            return _np.array([]), _np.array([]), _np.array([])
        t = res[0]
        return (
            _np.asarray(t["RAJ2000"], float),
            _np.asarray(t["DEJ2000"], float),
            _np.asarray(t["Ftot"], float),
        )

    lra, ldec, lf = _band(VIZIER_RACS["low"])
    if lra.size == 0:
        return {
            "ra": _np.array([]),
            "dec": _np.array([]),
            "flux": _np.zeros((0, 3)),
            "eflux": _np.zeros((0, 3)),
        }
    flux = _np.full((lra.size, 3), _np.nan)
    flux[:, 0] = lf
    for col, band in ((1, "mid"), (2, "high")):
        bra, bdec, bf = _band(VIZIER_RACS[band])
        if bra.size:
            il, ib, _ = crossmatch(lra, ldec, bra, bdec, match_arcsec)
            flux[il, col] = bf[ib]
    return {"ra": lra, "dec": ldec, "flux": flux, "eflux": rel_err * flux}


def validate_callingham(
    *, max_sources: int = 50, cone_deg: float = 0.06
) -> dict:  # pragma: no cover - network
    r"""Recover-a-known: do we *measure* the Callingham et al. (2017) published turnover $\nu_\mathrm{pk}$?

    The northern slice could only show purity against Callingham (it bounded, never measured, the
    turnover). Here we test the southern method's headline ability: for each Callingham peaked source
    (VizieR ``J/ApJ/836/174/pkfreq``, which gives a measured $\nu_\mathrm{pk}$) we fetch a small
    GLEAM-X$+$RACS cone, run :func:`find_peaked_south`, and compare *our* measured $\nu_\mathrm{pk}$ to
    theirs. Reports the recovery binned by published $\nu_\mathrm{pk}$ and the agreement (median
    $|\Delta\log_{10}\nu_\mathrm{pk}|$ and the fraction within 0.3 dex) for the recovered sources.
    Because Callingham is GLEAM-selected and MHz-peaked-dominated, sources peaking below the GLEAM band
    fall (correctly) outside this method's rising-side window --- the recovery is expected to climb with
    $\nu_\mathrm{pk}$.
    """
    import numpy as _np
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    from .spectra import crossmatch

    v = Vizier(columns=["RAJ2000", "DEJ2000", "nuPk"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs("J/ApJ/836/174/pkfreq")[0]
    ra = _np.asarray(t["RAJ2000"], float)
    dec = _np.asarray(t["DEJ2000"], float)
    nupk = _np.asarray(t["nuPk"], float)  # MHz, published
    # GLEAM-X footprint: Dec < +30, and RA in the GLEAM-X DR2 strip (20h40m-6h40m -> >310 or <100 deg)
    sel = _np.where((dec < 28) & _np.isfinite(nupk) & ((ra > 310) | (ra < 100)))[0]
    sel = sel[:max_sources]  # bound the per-source fetches; report how many are recovered

    pub: list[float] = []
    meas: list[float] = []
    pub_any: list[float] = []
    meas_any: list[float] = []
    disposition: dict[str, int] = {
        "no_gleamx_rows": 0,
        "no_match": 0,
        "fetch_error": 0,
        "recovered_peaked": 0,
        "not_recovered": 0,
    }
    tested_pub: list[float] = []
    for i in sel:
        center = SkyCoord(ra[i], dec[i], unit="deg")
        try:
            gleamx = fetch_gleamx(center, cone_deg)
            if gleamx["ra"].size == 0:
                disposition["no_gleamx_rows"] += 1
                continue
            racs = fetch_racs_bands(center, cone_deg)
            res = find_peaked_south(gleamx, racs, radius_arcsec=30.0)
        except Exception:
            # counted separately: a VizieR timeout is NOT an astrophysical non-recovery, and
            # folding the two together made the 38/50 rate unauditable (referee finding)
            disposition["fetch_error"] += 1
            continue
        if res["ra"].size == 0:
            disposition["no_match"] += 1
            continue
        _, jr, _ = crossmatch(_np.array([ra[i]]), _np.array([dec[i]]), res["ra"], res["dec"], 30.0)
        if jr.size == 0:
            disposition["no_match"] += 1
            continue
        ridx = int(jr[0])  # the matched source index in res for this Callingham position
        tested_pub.append(float(nupk[i]))
        # UNCONDITIONAL turnover comparison: every covered source with a finite fitted peak,
        # whether or not it passes the selection gate. The gated-only statistic routes any
        # catastrophically mis-measured turnover into "not recovered" and so cannot contain a
        # large error by construction (referee finding); both are reported.
        if _np.isfinite(res["nu_pk_ghz"][ridx]):
            pub_any.append(float(nupk[i]))
            meas_any.append(float(res["nu_pk_ghz"][ridx] * 1e3))
        if res["is_peaked"][ridx] and _np.isfinite(res["nu_pk_ghz"][ridx]):
            disposition["recovered_peaked"] += 1
            pub.append(float(nupk[i]))
            meas.append(float(res["nu_pk_ghz"][ridx] * 1e3))  # GHz -> MHz
        else:
            disposition["not_recovered"] += 1

    pub_a = _np.asarray(pub)
    meas_a = _np.asarray(meas)
    n = int(pub_a.size)
    n_covered = disposition["recovered_peaked"] + disposition["not_recovered"]
    bins = [(72, 250), (250, 500), (500, 2000)]
    # per-bin TESTED and RECOVERED counts, so a rate is expressible (a recovered-only count
    # cannot support "recovery climbs with nu_pk" -- another referee finding)
    tp = _np.asarray(tested_pub)
    recovery = {
        f"{lo}-{hi}MHz": {
            "tested": int(((tp >= lo) & (tp < hi)).sum()),
            "recovered": int(((pub_a >= lo) & (pub_a < hi)).sum()),
        }
        for lo, hi in bins
    }
    dlog = _np.abs(_np.log10(meas_a / pub_a)) if n else _np.array([])
    pa, ma = _np.asarray(pub_any), _np.asarray(meas_any)
    dlog_any = _np.abs(_np.log10(ma / pa)) if pa.size else _np.array([])
    return {
        "n_tried": int(sel.size),
        "n_with_coverage": n_covered,
        "disposition": disposition,
        "n_recovered_peaked": n,
        "recovered_by_published_nupk": recovery,
        "median_abs_dlog_nupk": float(_np.median(dlog)) if n else float("nan"),
        "frac_within_0p3dex": float((dlog < 0.3).mean()) if n else 0.0,
        "n_unconditional": int(pa.size),
        "median_abs_dlog_nupk_unconditional": float(_np.median(dlog_any))
        if pa.size
        else float("nan"),
        "frac_within_factor_two_of_tested": float(
            (_np.abs(_np.log10(ma / pa)) < _np.log10(2.0)).sum() / max(len(tested_pub), 1)
        )
        if pa.size
        else 0.0,
    }


def run(
    center=None,
    radius_deg: float = 3.0,
    out: str = ".",
    *,
    offline: bool = True,
    validate: bool = False,
) -> dict:
    """Full slice: synthesise (or fetch) GLEAM-X×RACS, fit curvature, classify, write artifacts.

    With ``validate`` (real data only), also runs the Callingham (2017) recover-a-known
    (:func:`validate_callingham`) and folds its headline numbers into the metrics and macros.
    """
    from pathlib import Path

    truth_pk: np.ndarray | None
    truth_uss: np.ndarray | None
    if offline or center is None:
        gleamx, racs, truth_pk, truth_uss = synthetic_field()
        source = "synthetic"
    else:  # pragma: no cover - network
        gleamx = fetch_gleamx(center, radius_deg)
        racs = fetch_racs_bands(center, radius_deg)
        truth_pk = truth_uss = None
        source = (
            f"GLEAM-X DR2 x RACS @ ({center.ra.deg:.1f},{center.dec.deg:.1f}) r={radius_deg}deg"
        )

    res = find_peaked_south(gleamx, racs)
    cls = res["cls"]
    area = (
        float("nan")
        if (offline or center is None)
        else float(2.0 * np.pi * (1.0 - np.cos(np.radians(radius_deg))) * (180.0 / np.pi) ** 2)
    )
    n_peaked = int(res["is_peaked"].sum())
    pk_chi2 = res["chi2_red"][res["is_peaked"]]
    metrics = {
        "source": source,
        "n_matched": int(cls.size),
        # the selection cascade, committed (the paper's "~16%" lived only in prose before)
        "n_peaked_naive": int(res["n_peaked_naive"]),
        "n_peaked_after_rising": int(res["n_peaked_after_rising"]),
        "n_peaked": n_peaked,
        "n_uss": int(res["is_uss"].sum()),
        "n_extended": int(np.sum(cls == "extended")),
        "n_racs_matched_multiply": int(res["n_racs_matched_multiply"]),
        "median_chi2_red_peaked": (
            round(float(np.nanmedian(pk_chi2)), 2) if np.isfinite(pk_chi2).any() else None
        ),
        "median_nu_pk_mhz": (
            round(1e3 * float(np.nanmedian(res["nu_pk_ghz"][res["is_peaked"]])), 1)
            if res["is_peaked"].any()
            else 0.0
        ),
        # A turnover between the GLEAM band top and RACS-low is interpolated across the
        # 0.23-0.89 GHz sampling gap, not sampled; the fraction is the one number that keeps
        # the title's "measuring" honest (referee finding).
        "peaked_in_gap_fraction": (
            round(
                float(
                    np.mean(
                        (res["nu_pk_ghz"][res["is_peaked"]] > GLEAMX_NU_GHZ[-1])
                        & (res["nu_pk_ghz"][res["is_peaked"]] < RACS_NU_GHZ[0])
                    )
                ),
                2,
            )
            if res["is_peaked"].any()
            else None
        ),
        "area_deg2": round(area, 2) if np.isfinite(area) else None,
        "peaked_density_per_deg2": (
            round(n_peaked / area, 3) if np.isfinite(area) and area > 0 else None
        ),
        # RadioSED II (Kerrison et al. 2025) Stripe 82 density, for the (uncontrolled: different
        # survey depths and selection) over-density indication the paper reports as such
        "kerrison_density_per_deg2": KERRISON_DENSITY_PER_DEG2,
        "density_ratio_vs_kerrison": (
            round(n_peaked / area / KERRISON_DENSITY_PER_DEG2, 2)
            if np.isfinite(area) and area > 0
            else None
        ),
    }
    if truth_pk is not None:
        from .spectra import crossmatch

        for key, truth in (("peaked", truth_pk), ("uss", truth_uss)):
            assert truth is not None  # both set together in the synthetic path
            mask = res["is_peaked"] if key == "peaked" else res["is_uss"]
            rt = np.flatnonzero(truth)
            ra_t = np.asarray(gleamx["ra"], float)[rt]
            dec_t = np.asarray(gleamx["dec"], float)[rt]
            if mask.any() and rt.size:
                i, _, _ = crossmatch(ra_t, dec_t, res["ra"][mask], res["dec"][mask], 5.0)
                metrics[f"n_{key}_recovered"] = int(i.size)
            else:
                metrics[f"n_{key}_recovered"] = 0
            metrics[f"n_injected_{key}"] = int(truth.sum())

    if validate and not offline and center is not None:  # pragma: no cover - network
        cal = validate_callingham()
        metrics["call_tried"] = int(cal["n_tried"])
        metrics["call_with_coverage"] = int(cal["n_with_coverage"])
        metrics["call_disposition"] = cal["disposition"]
        metrics["call_recovered"] = int(cal["n_recovered_peaked"])
        metrics["call_by_published_nupk"] = cal["recovered_by_published_nupk"]
        metrics["call_dlog"] = round(float(cal["median_abs_dlog_nupk"]), 3)
        metrics["call_within_pct"] = round(100.0 * float(cal["frac_within_0p3dex"]))
        metrics["call_n_unconditional"] = int(cal["n_unconditional"])
        metrics["call_dlog_unconditional"] = round(
            float(cal["median_abs_dlog_nupk_unconditional"]), 3
        )
        metrics["call_within_factor_two_of_tested_pct"] = round(
            100.0 * float(cal["frac_within_factor_two_of_tested"])
        )

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "southern_metrics.json")
    if not offline:  # pragma: no cover - the real catalogue is committed evidence
        _write_catalogue(res, op / "results" / "southern_candidates.csv", source)
    _figure(gleamx, racs, res, op / "papers" / "southern" / "figures")
    _write_macros(metrics, op / "papers" / "southern" / "generated" / "macros.tex")
    return metrics


def _write_catalogue(res: dict, path, source: str) -> None:
    """Commit the per-source catalogue: a scout-catalogue paper must publish its catalogue.

    The previous evidence file held nine scalars — no positions, fluxes, turnovers or classes
    for any of the 90 candidates, so every population claim was unauditable (referee finding).
    """
    import csv as _csv
    from pathlib import Path

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="") as fh:
        fh.write(f"# {source}\n")
        w = _csv.writer(fh)
        w.writerow(
            [
                "ra",
                "dec",
                "cls",
                "nu_pk_mhz",
                "alpha_lo",
                "alpha_lo_err",
                "alpha_hi",
                "chi2_red",
                "n_fit_points",
            ]
        )
        for i in range(res["ra"].size):
            nupk = res["nu_pk_ghz"][i]
            w.writerow(
                [
                    f"{res['ra'][i]:.5f}",
                    f"{res['dec'][i]:.5f}",
                    res["cls"][i],
                    f"{1e3 * nupk:.1f}" if np.isfinite(nupk) else "",
                    f"{res['alpha_lo'][i]:.3f}" if np.isfinite(res["alpha_lo"][i]) else "",
                    f"{res['alpha_lo_err'][i]:.3f}" if np.isfinite(res["alpha_lo_err"][i]) else "",
                    f"{res['alpha_hi'][i]:.3f}" if np.isfinite(res["alpha_hi"][i]) else "",
                    f"{res['chi2_red'][i]:.2f}" if np.isfinite(res["chi2_red"][i]) else "",
                    int(res["n_fit_points"][i]),
                ]
            )


def _figure(gleamx: dict, racs: dict, res: dict, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    nu = np.concatenate([GLEAMX_NU_GHZ, RACS_NU_GHZ])
    fig, ax = plt.subplots(figsize=(5, 4))
    pk = np.flatnonzero(res.get("is_peaked", np.zeros(0, bool)))[:6]
    from .spectra import crossmatch

    ig, ir, _ = crossmatch(gleamx["ra"], gleamx["dec"], racs["ra"], racs["dec"], 25.0)
    for k in pk:
        gi, rj = int(ig[k]), int(ir[k])
        sed = np.concatenate([gleamx["flux"][gi], racs["flux"][rj]])
        ax.plot(nu * 1e3, sed, ".-", lw=0.8, ms=3)
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="frequency (MHz)",
        ylabel="flux density (mJy)",
        title="Peaked-spectrum SEDs (GLEAM-X + RACS)",
    )
    fig.tight_layout()
    fig.savefig(out / "seds.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    """Emit BOTH namespaces, merged, never overwritten.

    Every count here is mode-dependent: the synthetic field and the real cone produce
    different values for the same quantity. Sharing one macro name meant an offline rebuild
    replaced 1545 real matches with the synthetic field's count -- the documented
    ``\tiiNEvents`` clobber, which ships a wrong number rather than a hole. So the counts are
    namespaced ``soSyn*``/``soReal*`` and only the active mode's are filled.

    The Callingham validation runs on the real path only, so its macros default to a ``--``
    placeholder rather than to 0; ``preserve_live_macros`` then guarantees a real value is
    never replaced by that placeholder. Without it, one ``make figures`` (which runs every
    slice offline in the repo root) turned the abstract's validation sentence into
    "of 0 tested, 0 are recovered".
    """
    from pathlib import Path

    from .report import preserve_live_macros

    real = not str(m.get("source", "")).startswith("synthetic")
    ns = "soReal" if real else "soSyn"
    other = "soSyn" if real else "soReal"

    def g(key: str, default: str = "--") -> str:
        v = m.get(key)
        return default if v is None else str(v)

    counts = (
        ("Nmatched", "n_matched"),
        ("NpeakedNaive", "n_peaked_naive"),
        ("NpeakedRising", "n_peaked_after_rising"),
        ("Npeaked", "n_peaked"),
        ("Nuss", "n_uss"),
        ("Nextended", "n_extended"),
        ("NracsMulti", "n_racs_matched_multiply"),
        ("MedianNupk", "median_nu_pk_mhz"),
        ("GapFrac", "peaked_in_gap_fraction"),
        ("MedianChisq", "median_chi2_red_peaked"),
        ("Area", "area_deg2"),
        ("Density", "peaked_density_per_deg2"),
        ("KerrisonDensity", "kerrison_density_per_deg2"),
        ("DensityRatio", "density_ratio_vs_kerrison"),
    )
    lines = [
        "% Auto-generated by jansky_research.southern._write_macros — do not edit by hand.",
        "% Synthetic (soSyn*) and real (soReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders, so synthetic counts can never masquerade",
        "% under soReal*. The file is MERGED, not overwritten (report.preserve_live_macros).",
        rf"\newcommand{{\soSource}}{{{m['source']}}}",
        rf"\newcommand{{\soUssThresh}}{{{USS_THRESHOLD}}}",
        rf"\newcommand{{\soCompactMax}}{{{COMPACT_RATIO_MAX}}}",
    ]
    for suffix, key in counts:
        lines.append(rf"\newcommand{{\{ns}{suffix}}}{{{g(key)}}}")
        lines.append(rf"\newcommand{{\{other}{suffix}}}{{--}}")
    # Real-path-only validation; a placeholder here is honest, a 0 is a wrong number.
    for suffix, key in (
        ("CallTried", "call_tried"),
        ("CallCovered", "call_with_coverage"),
        ("CallRecovered", "call_recovered"),
        ("CallDlog", "call_dlog"),
        ("CallWithin", "call_within_pct"),
        ("CallNuncond", "call_n_unconditional"),
        ("CallDlogUncond", "call_dlog_unconditional"),
        ("CallWithinTwoOfTested", "call_within_factor_two_of_tested_pct"),
    ):
        lines.append(rf"\newcommand{{\soReal{suffix}}}{{{g(key)}}}")
    # per-bin tested/recovered counts, so "recovery climbs with nu_pk" is macro-backed
    bins = m.get("call_by_published_nupk") or {}
    for name, key in (("Low", "72-250MHz"), ("Mid", "250-500MHz"), ("Hi", "500-2000MHz")):
        row = bins.get(key) or {}
        for part in ("tested", "recovered"):
            val = row.get(part)
            lines.append(
                rf"\newcommand{{\soRealCallBin{name}{part.capitalize()}}}"
                rf"{{{'--' if val is None else val}}}"
            )
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    from astropy.coordinates import SkyCoord

    p = argparse.ArgumentParser(description="Southern peaked-spectrum selection (GLEAM-X + RACS).")
    p.add_argument("--ra", type=float, help="field-centre RA (deg)")
    p.add_argument("--dec", type=float, help="field-centre Dec (deg)")
    p.add_argument("--radius", type=float, default=3.0, help="cone radius (deg)")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument(
        "--validate", action="store_true", help="also run the Callingham recover-a-known"
    )
    args = p.parse_args(argv)
    center = None if (args.offline or args.ra is None) else SkyCoord(args.ra, args.dec, unit="deg")
    metrics = run(center, args.radius, args.out, offline=args.offline, validate=args.validate)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
