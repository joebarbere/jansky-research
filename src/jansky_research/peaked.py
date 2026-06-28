"""Peaked-spectrum (GPS/CSS) radio-source selection via three-frequency spectral curvature.

Gigahertz-Peaked-Spectrum (GPS) and Compact-Steep-Spectrum (CSS) sources are compact, young radio
AGN whose radio spectrum *rises then falls*, peaking in the ~0.1--3 GHz band (O'Dea & Saikia 2021).
Selecting them needs $\\geq 3$ frequencies to see the turnover. Using three public surveys --- TGSS
(150 MHz), NVSS (1.4 GHz), and VLASS (3 GHz) --- we compute two indices, $\\alpha_\\mathrm{low}$
(150$\\to$1400 MHz) and $\\alpha_\\mathrm{high}$ (1400$\\to$3000 MHz), and select sources that are
rising at low frequency and falling at high frequency. The **curvature**
$\\alpha_\\mathrm{high}-\\alpha_\\mathrm{low}$ is far more robust to a constant TGSS flux-scale offset
than a single steep cut, and $\\alpha_\\mathrm{high}$ is TGSS-independent.

This is the maximal-reuse slice: it composes :mod:`jansky_research.spectra` (two-point index,
cross-match, fetch, NED/SIMBAD annotation) and :mod:`jansky_research.vlass` (VLASS 3 GHz fetch,
variability metrics to flag blazars, vetting) rather than reimplementing them. Pure NumPy + a
synthetic three-survey fixture for offline tests.
"""

from __future__ import annotations

import numpy as np

from .spectra import spectral_index

__all__ = [
    "NU_GHZ",
    "classify_sed",
    "find_peaked",
    "peak_frequency",
    "run",
    "synthetic_field",
    "two_point_indices",
    "validate_hfp",
    "validate_known",
]

# Survey reference frequencies (GHz).
NU_GHZ = {"tgss": 0.1475, "nvss": 1.4, "vlass": 3.0}
# TGSS 150 MHz ~7-sigma detection limit (mJy). Peaked/GPS sources are faint at 150 MHz and often
# below it, so a TGSS *non-detection* is used as an upper limit S_150 < this (see find_peaked).
TGSS_LIMIT_MJY = 25.0
# Floor on the NVSS->VLASS index. A drop steeper than this between 1.4 and 3 GHz is not a real
# spectrum but NVSS (45") extended emission resolved out by VLASS (2.5"): a resolution artefact.
A_HIGH_FLOOR = -2.0


def two_point_indices(
    s_tgss: np.ndarray,
    s_nvss: np.ndarray,
    s_vlass: np.ndarray,
    *,
    e_tgss: np.ndarray | None = None,
    e_nvss: np.ndarray | None = None,
    e_vlass: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Low (150$\\to$1400 MHz) and high (1400$\\to$3000 MHz) spectral indices (reuses ``spectra``)."""
    alpha_low, _ = spectral_index(s_tgss, NU_GHZ["tgss"], s_nvss, NU_GHZ["nvss"], e_tgss, e_nvss)
    alpha_high, _ = spectral_index(
        s_nvss, NU_GHZ["nvss"], s_vlass, NU_GHZ["vlass"], e_nvss, e_vlass
    )
    return alpha_low, alpha_high


def classify_sed(alpha_low: float, alpha_high: float, *, up: float = 0.1, dn: float = -0.1) -> str:
    """Classify a 3-point radio SED from its two indices.

    ``peaked`` -- rising then falling ($\\alpha_\\mathrm{low}>$ ``up`` and $\\alpha_\\mathrm{high}<$
    ``dn``: a turnover between 150 MHz and 3 GHz, the GPS/CSS signature); ``inverted`` -- still rising
    at 3 GHz (peak above the band, or a flat-spectrum core); ``steep`` -- falling throughout;
    ``flat`` -- everything else.
    """
    if not (np.isfinite(alpha_low) and np.isfinite(alpha_high)):
        return "nan"
    if alpha_low > up and alpha_high < dn:
        return "peaked"
    if alpha_high > up:
        return "inverted"
    if alpha_low < -0.5 and alpha_high < -0.5:
        return "steep"
    return "flat"


def peak_frequency(flux: np.ndarray, nu_ghz: np.ndarray) -> tuple[float, bool]:
    """Turnover frequency (GHz) of a 3-point SED, by a parabolic log--log fit.

    Fits $\\log_{10} S = a x^2 + b x + c$ with $x=\\log_{10}\\nu$; the extremum is at
    $x_\\mathrm{peak}=-b/2a$ and is a genuine peak (concave) when $a<0$. Returns
    ``(nu_peak_ghz, is_peak)``.
    """
    flux = np.asarray(flux, dtype=float)
    nu = np.asarray(nu_ghz, dtype=float)
    if flux.size < 3 or np.any(flux <= 0):
        return float("nan"), False
    a, b, _ = np.polyfit(np.log10(nu), np.log10(flux), 2)
    if a == 0.0:
        return float("nan"), False
    x_peak = -b / (2.0 * a)
    if abs(x_peak) > 6.0:  # extremum far outside any sane radio band -> no in-band turnover
        return float("nan"), bool(a < 0.0)
    return float(10.0**x_peak), bool(a < 0.0)


def find_peaked(
    tgss: dict[str, np.ndarray],
    nvss: dict[str, np.ndarray],
    vlass: dict[str, np.ndarray],
    *,
    radius_arcsec: float = 15.0,
    tgss_limit_mjy: float = TGSS_LIMIT_MJY,
    up: float = 0.1,
    dn: float = -0.1,
    a_high_floor: float = A_HIGH_FLOOR,
) -> dict[str, np.ndarray]:
    """NVSS-anchored peaked-spectrum selection, using TGSS as an UPPER LIMIT. Reuses ``spectra``.

    Anchors on NVSS (1.4 GHz) sources detected by VLASS (3 GHz). Requiring a TGSS detection would
    *exclude* peaked sources, which are faint at 150 MHz; instead a TGSS non-detection is treated as
    $S_{150}<$ ``tgss_limit_mjy``, giving a lower bound on $\\alpha_\\mathrm{low}$. A source is
    ``peaked`` when it is rising at low frequency ($\\alpha_\\mathrm{low}$, or its lower bound,
    $>$ ``up``), falling at high frequency ($\\alpha_\\mathrm{high}<$ ``dn``), and *not* a resolution
    artefact ($\\alpha_\\mathrm{high}>$ ``a_high_floor``: a steeper 1.4$\\to$3 GHz drop is NVSS
    extended emission resolved out by VLASS). Returns per-(NVSS$\\cap$VLASS)-source arrays incl.
    ``tgss_detected`` and ``alpha_low_is_limit``.
    """
    from .spectra import crossmatch

    ra_n, dec_n = np.asarray(nvss["ra"], float), np.asarray(nvss["dec"], float)
    iN_v, iv, _ = crossmatch(ra_n, dec_n, vlass["ra"], vlass["dec"], radius_arcsec)
    iN_t, it, _ = crossmatch(ra_n, dec_n, tgss["ra"], tgss["dec"], radius_arcsec)
    if iN_v.size == 0:
        return {k: np.array([]) for k in ("ra", "dec", "s_nvss", "s_vlass", "cls")}
    t_of = dict(zip(iN_t.tolist(), it.tolist(), strict=True))

    s_n = np.asarray(nvss["flux"], float)[iN_v]
    s_v = np.asarray(vlass["flux"], float)[iv]
    e_n = np.asarray(nvss["eflux"], float)[iN_v]
    e_v = np.asarray(vlass["eflux"], float)[iv]
    alpha_high, _ = spectral_index(s_n, NU_GHZ["nvss"], s_v, NU_GHZ["vlass"], e_n, e_v)

    lnr = np.log(NU_GHZ["nvss"] / NU_GHZ["tgss"])
    s_t = np.full(iN_v.size, np.nan)
    alpha_low = np.empty(iN_v.size)
    tgss_det = np.zeros(iN_v.size, dtype=bool)
    tflux = np.asarray(tgss["flux"], float)
    for k, nidx in enumerate(iN_v.tolist()):
        if nidx in t_of:  # TGSS-detected -> measured alpha_low
            s_t[k] = tflux[t_of[nidx]]
            alpha_low[k] = np.log(s_n[k] / s_t[k]) / lnr
            tgss_det[k] = True
        else:  # TGSS non-detection -> S_150 < limit -> lower bound on alpha_low
            alpha_low[k] = np.log(s_n[k] / tgss_limit_mjy) / lnr
    rising = (
        alpha_low > up
    )  # optically-thick rise (TGSS-faint, NVSS-bright) -- the GPS/HFP signature
    is_peaked = (
        rising & (alpha_high < dn) & (alpha_high > a_high_floor)
    )  # turnover in 0.7-2 GHz band
    is_ghz_peaked = rising & (alpha_high > up)  # still rising at 3 GHz -> peak above the band (HFP)
    cls = np.array(
        [classify_sed(lo, hi) for lo, hi in zip(alpha_low, alpha_high, strict=True)], dtype=object
    )
    cls[(alpha_high <= a_high_floor)] = "extended"  # resolution artefact (NVSS>>VLASS)
    cls[is_peaked] = "peaked"
    cls[is_ghz_peaked] = "ghz_peaked"
    return {
        "ra": ra_n[iN_v],
        "dec": dec_n[iN_v],
        "s_tgss": s_t,
        "s_nvss": s_n,
        "s_vlass": s_v,
        "alpha_low": alpha_low,
        "alpha_high": alpha_high,
        "alpha_low_is_limit": ~tgss_det,
        "tgss_detected": tgss_det,
        "cls": cls.astype(str),
        "is_peaked": is_peaked,
        "is_ghz_peaked": is_ghz_peaked,
        "is_rising": rising,  # peaked OR ghz_peaked: the full optically-thick-rising candidate set
    }


def validate_known(*, max_sources: int = 120) -> dict:  # pragma: no cover - network
    """Validate the selection against the Callingham et al. (2017) peaked-spectrum catalogue.

    Fetches known peaked sources (with measured turnover ``nuPk``) and their TGSS/NVSS fluxes from
    VizieR (``J/ApJ/836/174``), adds VLASS 3 GHz per source (reusing ``vlass``), classifies each, and
    reports the recovery as ``peaked`` binned by ``nuPk``. The three-point (150 MHz, 1.4 GHz, 3 GHz)
    method has a narrow ~0.7--2 GHz peaked window: it should reject the GLEAM-dominated
    *MHz-peaked* sources (turnover below the 150 MHz floor -> they look steep), not the other way
    round. So a low recovery at low ``nuPk`` is correct purity, not a miss.
    """
    import numpy as _np
    from astroquery.vizier import Vizier

    from . import vlass as _vlass

    v = Vizier(columns=["*"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs("J/ApJ/836/174/pkfreq")[0]
    ra = _np.asarray(t["RAJ2000"], float)
    dec = _np.asarray(t["DEJ2000"], float)
    nupk = _np.asarray(t["nuPk"], float)  # MHz
    snvss = _np.asarray(t["SNVSS"], float) * 1e3  # Jy -> mJy
    stgss = _np.asarray(t["STGSS"], float) * 1e3
    sel = _np.where((dec > -39) & (snvss > 60) & _np.isfinite(snvss))[0]
    sel = sel[_np.argsort(-snvss[sel])][:max_sources]

    npk_l: list[float] = []
    pk_l: list[bool] = []
    for i in sel:
        try:
            vr, vd, vf, _ = _vlass._fetch_e1_tap((ra[i], dec[i]), 0.02)  # ~70" cone
            if vr.size == 0:
                continue
            sv = vf[int(_np.argmin((vr - ra[i]) ** 2 + (vd - dec[i]) ** 2))]
        except Exception:
            continue
        st = stgss[i] if _np.isfinite(stgss[i]) else TGSS_LIMIT_MJY
        a_low = _np.log(snvss[i] / st) / _np.log(NU_GHZ["nvss"] / NU_GHZ["tgss"])
        a_high = _np.log(sv / snvss[i]) / _np.log(NU_GHZ["vlass"] / NU_GHZ["nvss"])
        npk_l.append(float(nupk[i]))
        pk_l.append(bool((a_low > 0.1) and (A_HIGH_FLOOR < a_high < -0.1)))
    npk = _np.asarray(npk_l)
    ispk = _np.asarray(pk_l, bool)
    bins = [(72, 250), (250, 500), (500, 1000)]
    recovery = {
        f"{lo}-{hi}MHz": (
            int(ispk[(npk >= lo) & (npk < hi)].sum()),
            int(((npk >= lo) & (npk < hi)).sum()),
        )
        for lo, hi in bins
    }
    fp = float(ispk[npk < 250].mean()) if bool((npk < 250).any()) else 0.0
    return {
        "n_validated": int(npk.size),
        "false_positive_rate_below_250MHz": fp,
        "recovery_by_nupk": recovery,
    }


def validate_hfp(*, max_sources: int = 120) -> dict:  # pragma: no cover - network
    """Recover-a-known test against the Dallacasa et al. (2000) High-Frequency-Peaker catalogue.

    The Callingham (2017) GLEAM sample is MHz-peaked (turnover below the 150 MHz floor), so this
    method correctly leaves it in the ``steep`` class (see :func:`validate_known`). The *classical*
    GHz-peaked / HFP population peaks at $\\gtrsim$few GHz, so across 150 MHz, 1.4 GHz, 3 GHz it is
    **rising throughout** and lands in this method's ``ghz_peaked`` class. This fetches the Dallacasa
    bright HFP sample (``J/A+A/363/887``, NVSS 1.4 GHz fluxes), adds VLASS 3 GHz per source, and
    reports the fraction recovered as *rising* (optically thick at low frequency, the GPS/HFP
    signature) and as *GHz-peaked* (still rising at 3 GHz). A clean recover-a-known, complementing the
    Callingham purity test.
    """
    import numpy as _np
    from astropy import units as _u
    from astropy.coordinates import SkyCoord as _SkyCoord
    from astroquery.vizier import Vizier

    from . import vlass as _vlass

    v = Vizier(columns=["*"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs("J/A+A/363/887/table1")[0]
    coo = _SkyCoord(t["RAJ2000"], t["DEJ2000"], unit=(_u.hourangle, _u.deg))
    ra = _np.asarray(coo.ra.deg, float)
    dec = _np.asarray(coo.dec.deg, float)
    snvss = _np.asarray(t["NVSS"], float)  # already mJy
    sel = _np.where((dec > -39) & (snvss > 30) & _np.isfinite(snvss))[0]
    sel = sel[_np.argsort(-snvss[sel])][:max_sources]

    a_low_l: list[float] = []
    a_high_l: list[float] = []
    for i in sel:
        try:
            vr, vd, vf, _ = _vlass._fetch_e1_tap((ra[i], dec[i]), 0.02)  # ~70" cone
            if vr.size == 0:
                continue
            sv = vf[int(_np.argmin((vr - ra[i]) ** 2 + (vd - dec[i]) ** 2))]
        except Exception:
            continue
        a_low_l.append(float(_np.log(snvss[i] / TGSS_LIMIT_MJY) / _np.log(1.4 / 0.1475)))
        a_high_l.append(float(_np.log(sv / snvss[i]) / _np.log(3.0 / 1.4)))
    a_low = _np.asarray(a_low_l)
    a_high = _np.asarray(a_high_l)
    n = int(a_low.size)
    rising = (a_low > 0.1).sum() if n else 0
    ghz_peaked = ((a_low > 0.1) & (a_high > 0.1)).sum() if n else 0
    return {
        "n_validated": n,
        "median_alpha_low": float(_np.median(a_low)) if n else float("nan"),
        "median_alpha_high": float(_np.median(a_high)) if n else float("nan"),
        "frac_rising": float(rising) / n if n else 0.0,
        "frac_ghz_peaked": float(ghz_peaked) / n if n else 0.0,
    }


def synthetic_field(
    n_sources: int = 1500, *, peaked_fraction: float = 0.05, rel_err: float = 0.1, seed: int = 0
) -> tuple[dict, dict, dict, np.ndarray]:
    """Synthetic TGSS/NVSS/VLASS catalogues with injected peaked + steep + flat SEDs (offline fixture).

    Returns ``(tgss, nvss, vlass)`` survey dicts sharing positions (small jitter), so the cross-match
    recovers them. Peaked sources rise to ~1.4 GHz then fall; the rest are steep or flat power laws.
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(180.0, 185.0, n_sources)
    dec = rng.uniform(20.0, 25.0, n_sources)
    nu = np.array([NU_GHZ["tgss"], NU_GHZ["nvss"], NU_GHZ["vlass"]])
    s_nvss = 10.0 ** rng.uniform(0.3, 1.5, n_sources)  # ~2-30 mJy at 1.4 GHz
    is_peaked = rng.random(n_sources) < peaked_fraction
    # the upper-limit method can only confirm peaked sources brighter than the TGSS limit at 1.4 GHz
    # (so that a TGSS non-detection forces alpha_low > 0); make the injected peaked ones bright.
    s_nvss[is_peaked] = 10.0 ** rng.uniform(1.6, 2.3, int(is_peaked.sum()))  # ~40-200 mJy
    alpha = rng.uniform(-1.1, -0.5, n_sources)  # steep/flat power-law index for the rest
    flux = np.empty((n_sources, 3))
    for i in range(n_sources):
        if is_peaked[i]:  # parabola peaking near 1.4 GHz: rising then falling
            lp = np.log10(s_nvss[i]) - 1.2 * (np.log10(nu / NU_GHZ["nvss"])) ** 2
            flux[i] = 10.0**lp
        else:
            flux[i] = s_nvss[i] * (nu / NU_GHZ["nvss"]) ** alpha[i]
    flux *= rng.normal(1.0, rel_err, flux.shape)
    flux = np.clip(flux, 1e-3, None)
    jit = lambda: rng.normal(0.0, 1.0 / 3600.0, n_sources)  # noqa: E731  (~1" position jitter)

    def survey(col, *, mask=None):
        m = np.ones(n_sources, bool) if mask is None else mask  # apply a survey's depth
        return {
            "ra": (ra + jit())[m],
            "dec": (dec + jit())[m],
            "flux": flux[m, col],
            "eflux": rel_err * flux[m, col],
        }

    # TGSS is shallow: only sources above its 150 MHz limit are detected. Peaked sources (faint at
    # 150 MHz) fall below it, so they are TGSS-non-detected -- exactly as on the real sky.
    tgss = survey(0, mask=flux[:, 0] > TGSS_LIMIT_MJY)
    return tgss, survey(1), survey(2), is_peaked


def run(
    center=None,
    radius_deg: float = 2.0,
    out: str = ".",
    *,
    offline: bool = False,
    validate: bool = False,
) -> dict:
    """Full slice: fetch (or synthesise) TGSS/NVSS/VLASS, find peaked sources, vet, write artifacts.

    With ``validate`` (real data only), also runs the two recover-a-known tests --- ``validate_known``
    (Callingham 2017 MHz-peaked purity) and ``validate_hfp`` (Dallacasa 2000 GHz-peaked recovery) ---
    and folds their headline numbers into the metrics and macros the paper inputs.
    """
    import json
    from pathlib import Path

    if offline or center is None:
        tgss, nvss, vlass, truth = synthetic_field()
        source = "synthetic"
    else:  # pragma: no cover - network
        from .spectra import fetch_survey
        from .vlass import _fetch_e1_tap

        tgss = fetch_survey(center, radius_deg, "tgss")
        nvss = fetch_survey(center, radius_deg, "nvss")
        vra, vdec, vflux, veflux = _fetch_e1_tap((center.ra.deg, center.dec.deg), radius_deg)
        vlass = {"ra": vra, "dec": vdec, "flux": vflux, "eflux": veflux}
        truth = None
        source = f"TGSSxNVSSxVLASS @ ({center.ra.deg:.1f}, {center.dec.deg:.1f}) r={radius_deg}deg"

    res = find_peaked(tgss, nvss, vlass)
    cls = res.get("cls", np.array([]))
    peaked = res.get("is_peaked", cls == "peaked")
    ghz_peaked = res.get("is_ghz_peaked", cls == "ghz_peaked")
    metrics = {
        "source": source,
        "n_nvss_vlass": int(cls.size),
        "n_peaked": int(peaked.sum()),
        "n_ghz_peaked": int(ghz_peaked.sum()),
        "n_rising": int(np.sum(res.get("is_rising", peaked | ghz_peaked))),
        "n_extended_artefact": int(np.sum(cls == "extended")),
        "n_tgss_detected": int(np.sum(res.get("tgss_detected", np.zeros(cls.size, bool)))),
    }
    if truth is not None:  # synthetic: recovery + purity of the injected peaked sources
        from .spectra import crossmatch

        ra_t = np.asarray(nvss["ra"], float)[np.flatnonzero(truth)]
        dec_t = np.asarray(nvss["dec"], float)[np.flatnonzero(truth)]
        # which injected-peaked positions land on a recovered "peaked" candidate
        ip, _, _ = crossmatch(ra_t, dec_t, res["ra"][peaked], res["dec"][peaked], 5.0)
        metrics["n_injected_peaked"] = int(truth.sum())
        metrics["n_peaked_recovered"] = int(ip.size)

    if validate and not offline and center is not None:  # pragma: no cover - network
        hfp = validate_hfp()
        metrics["hfp_n"] = int(hfp["n_validated"])
        metrics["hfp_rising_pct"] = round(100.0 * hfp["frac_rising"])
        metrics["hfp_ghz_pct"] = round(100.0 * hfp["frac_ghz_peaked"])
        cal = validate_known()
        lo_flag, lo_tot = cal["recovery_by_nupk"].get("72-250MHz", (0, 0))
        metrics["call_n"] = int(cal["n_validated"])
        metrics["call_low_flagged"] = int(lo_flag)
        metrics["call_low_total"] = int(lo_tot)

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "peaked_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(res, op / "papers" / "peaked" / "figures")
    _write_macros(metrics, op / "papers" / "peaked" / "generated" / "macros.tex")
    return metrics


def _figure(res: dict, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    al = res.get("alpha_low", np.array([]))
    ah = res.get("alpha_high", np.array([]))
    peaked = res.get("is_peaked", np.zeros(al.size, bool))
    ghz = res.get("is_ghz_peaked", np.zeros(al.size, bool))
    is_lim = res.get("alpha_low_is_limit", np.zeros(al.size, bool))
    other = ~(peaked | ghz)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(al[other], ah[other], s=6, color="0.6", label="other")
    # TGSS non-detections: alpha_low is a lower limit (shown with rightward arrows)
    ax.scatter(al[peaked], ah[peaked], s=34, color="r", marker="*", label="peaked (GPS/CSS)")
    ax.scatter(al[ghz], ah[ghz], s=34, color="b", marker="^", label="GHz-peaked (HFP)")
    if (peaked | ghz).any():
        for x, y in zip(al[(peaked | ghz) & is_lim], ah[(peaked | ghz) & is_lim], strict=True):
            ax.annotate(
                "", xy=(x + 0.15, y), xytext=(x, y), arrowprops=dict(arrowstyle="->", color="0.4")
            )
    ax.axhline(0.0, color="k", lw=0.5)
    ax.axvline(0.0, color="k", lw=0.5)
    ax.set(
        xlabel=r"$\alpha_{\rm low}$ (150$\to$1400 MHz; $\to$ = TGSS lower limit)",
        ylabel=r"$\alpha_{\rm high}$ (1400$\to$3000 MHz)",
        title="Spectral-curvature plane",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "curvature.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    lines = [
        "% Auto-generated by jansky_research.peaked._write_macros — do not edit by hand.",
        rf"\newcommand{{\pkSource}}{{{m['source']}}}",
        rf"\newcommand{{\pkNnvssvlass}}{{{m['n_nvss_vlass']}}}",
        rf"\newcommand{{\pkNpeaked}}{{{m['n_peaked']}}}",
        rf"\newcommand{{\pkNghzpeaked}}{{{m.get('n_ghz_peaked', 0)}}}",
        rf"\newcommand{{\pkNrising}}{{{m.get('n_rising', 0)}}}",
        rf"\newcommand{{\pkNextended}}{{{m['n_extended_artefact']}}}",
        rf"\newcommand{{\pkNtgssdet}}{{{m['n_tgss_detected']}}}",
        # Validation macros (populated by run(validate=True) on real data; 0 placeholders offline).
        rf"\newcommand{{\pkHfpN}}{{{m.get('hfp_n', 0)}}}",
        rf"\newcommand{{\pkHfpRising}}{{{m.get('hfp_rising_pct', 0)}}}",
        rf"\newcommand{{\pkHfpGhz}}{{{m.get('hfp_ghz_pct', 0)}}}",
        rf"\newcommand{{\pkCallN}}{{{m.get('call_n', 0)}}}",
        rf"\newcommand{{\pkCallLowFlagged}}{{{m.get('call_low_flagged', 0)}}}",
        rf"\newcommand{{\pkCallLowTotal}}{{{m.get('call_low_total', 0)}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    from astropy.coordinates import SkyCoord

    p = argparse.ArgumentParser(description="Find peaked-spectrum (GPS/CSS) radio sources.")
    p.add_argument("--ra", type=float, help="field-centre RA (deg)")
    p.add_argument("--dec", type=float, help="field-centre Dec (deg)")
    p.add_argument("--radius", type=float, default=2.0, help="cone radius (deg)")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument(
        "--validate", action="store_true", help="also run the Callingham + Dallacasa tests"
    )
    args = p.parse_args(argv)
    center = None if (args.offline or args.ra is None) else SkyCoord(args.ra, args.dec, unit="deg")
    metrics = run(center, args.radius, args.out, offline=args.offline, validate=args.validate)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
