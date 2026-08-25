"""Radio spectral-index tooling — the ultra-steep-spectrum (USS) source hunt.

Cross-matches two radio continuum surveys at different frequencies, computes the two-point
spectral index alpha (S_nu ~ nu^alpha) with propagated errors, classifies each source, and flags
**ultra-steep-spectrum** sources (alpha < -1.3, the classic high-redshift radio-galaxy selection;
De Breuck et al. 2000, A&AS 143, 303) and **anomalous positive-index** sources (alpha > 0 — a
calibration/variability/resolution flag, *not* a GPS/peaked-spectrum claim, which would need a
third frequency). The matched-catalogue approach follows de Gasperin, Intema & Frail (2018,
MNRAS 474, 5008; arXiv:1711.11367).

The default pairing is TGSS ADR1 (Intema et al. 2017, A&A 598, A78; 147.5 MHz) x NVSS (1.4 GHz) —
a decade-wide lever, both all-sky and openly cone-searchable on VizieR with no authentication.
Everything is pure NumPy + astropy; the survey fetch goes through ``astroquery.vizier``, and a
synthetic two-survey field (:func:`synthetic_field`) lets tests/CI run offline against known spectra.

Caveats the analysis must surface (and the write-up must state):
- TGSS ADR1 carries a position-dependent flux-scale SCATTER (Hurley-Walker 2017, arXiv:1703.06635;
  ~15-17%, approximately zero-mean over a field), which contributes ~0.06-0.07 to alpha per source
  -- scatter, not a systematic inflation (de Gasperin et al. 2018 SS4 bound its effect on their
  indices at ~0.06). The reference comparison and the selection-bias model below are how a raw run
  is scored; selection on the noisy alpha, not the flux scale, is what biases a USS candidate list.
- TGSS (25"), NVSS (45") and VLASS (2.5") differ in resolution; alpha from NVSS x TGSS (both recover
  extended flux) is the primary estimate, VLASS (QL Epoch 1, Gordon et al. 2021) a compact-source
  curvature check.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "NVSS_LIMIT_MJY",
    "SURVEYS",
    "Survey",
    "TGSS_LIMIT_MJY",
    "USS_THRESHOLD",
    "classify",
    "crossmatch",
    "fetch_reference_cone",
    "fetch_survey",
    "find_uss",
    "matched_sensitivity",
    "population_offset",
    "reference_crossmatch",
    "reference_spindex",
    "selection_bias_mc",
    "spectral_index",
    "synthetic_field",
    "uss_confusion",
]

USS_THRESHOLD = -1.3  # alpha below this = ultra-steep-spectrum (high-z radio-galaxy candidate)

# Nominal survey detection floors used for the matched-sensitivity truncation analysis.
# NVSS: ~50% completeness near 2.5 mJy (Condon et al. 1998). TGSS ADR1: 7 sigma at the median
# 3.5 mJy/beam rms (Intema et al. 2017). These are documented, tunable inputs, not fits.
NVSS_LIMIT_MJY = 2.5
TGSS_LIMIT_MJY = 24.5


@dataclass(frozen=True)
class Survey:
    """A radio continuum survey on VizieR, normalised to (ra, dec, flux_mjy, eflux_mjy)."""

    name: str
    vizier: str
    freq_mhz: float
    flux_col: str
    eflux_col: str


SURVEYS: dict[str, Survey] = {
    "tgss": Survey("TGSS ADR1", "J/A+A/598/A78", 147.5, "Stotal", "e_Stotal"),
    "nvss": Survey("NVSS", "VIII/65", 1400.0, "S1.4", "e_S1.4"),
    "vlass": Survey("VLASS QL", "J/ApJS/255/30", 3000.0, "Ftot", "e_Ftot"),
}


def spectral_index(
    s_lo: np.ndarray,
    nu_lo: float,
    s_hi: np.ndarray,
    nu_hi: float,
    e_lo: np.ndarray | None = None,
    e_hi: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Two-point spectral index alpha and its 1-sigma error.

    ``alpha = ln(S_hi / S_lo) / ln(nu_hi / nu_lo)`` for ``S_nu ~ nu^alpha``. The error propagates
    the fractional flux errors: ``sigma_alpha = sqrt((e_hi/S_hi)^2 + (e_lo/S_lo)^2) / |ln(nu_hi/nu_lo)|``.
    """
    s_lo = np.asarray(s_lo, float)
    s_hi = np.asarray(s_hi, float)
    lnr = np.log(nu_hi / nu_lo)
    alpha = np.log(s_hi / s_lo) / lnr
    if e_lo is None or e_hi is None:
        return alpha, np.full_like(alpha, np.nan)
    frac = np.sqrt((np.asarray(e_hi, float) / s_hi) ** 2 + (np.asarray(e_lo, float) / s_lo) ** 2)
    return alpha, frac / abs(lnr)


def classify(alpha: float) -> str:
    """Label a spectral index: uss / steep / flat / inverted (anomalous +ve index, not a GPS claim)."""
    if alpha < USS_THRESHOLD:
        return "uss"
    if alpha < -0.5:
        return "steep"
    if alpha < 0.0:
        return "flat"
    return "inverted"


def crossmatch(
    ra_lo: np.ndarray,
    dec_lo: np.ndarray,
    ra_hi: np.ndarray,
    dec_hi: np.ndarray,
    radius_arcsec: float = 15.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Positional cross-match: for each low-frequency source find the nearest high-freq source.

    Returns ``(idx_lo, idx_hi, sep_arcsec)`` for the matches within ``radius_arcsec`` — the
    low-frequency rows that have a counterpart, the matched high-frequency row indices, and the
    separations. Uses astropy's KD-tree matcher.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    lo = SkyCoord(ra_lo * u.deg, dec_lo * u.deg)
    hi = SkyCoord(ra_hi * u.deg, dec_hi * u.deg)
    idx, sep, _ = lo.match_to_catalog_sky(hi)
    keep = sep.arcsec <= radius_arcsec
    return np.flatnonzero(keep), idx[keep], sep.arcsec[keep]


def fetch_survey(
    center, radius_deg: float, survey: str
) -> dict[str, np.ndarray]:  # pragma: no cover - network
    """Cone-search a survey on VizieR and return normalised arrays (ra, dec, flux_mjy, eflux_mjy)."""
    from astropy import units as u
    from astroquery.vizier import Vizier

    spec = SURVEYS[survey]
    v = Vizier(columns=["*"])
    v.ROW_LIMIT = -1
    res = v.query_region(center, radius=radius_deg * u.deg, catalog=spec.vizier)
    if not res:
        return {
            "ra": np.array([]),
            "dec": np.array([]),
            "flux": np.array([]),
            "eflux": np.array([]),
        }
    t = res[0]
    flux = np.asarray(t[spec.flux_col], float)
    eflux = (
        np.asarray(t[spec.eflux_col], float)
        if spec.eflux_col in t.colnames
        else np.full(len(t), np.nan)
    )
    # Some VizieR catalogues serve RAJ2000/DEJ2000 as decimal degrees (TGSS), others as
    # sexagesimal strings (NVSS gives "11 46 14.18"). Parse both robustly.
    try:
        ra = np.asarray(t["RAJ2000"], float)
        dec = np.asarray(t["DEJ2000"], float)
    except (ValueError, TypeError):
        from astropy.coordinates import SkyCoord

        sc = SkyCoord(t["RAJ2000"], t["DEJ2000"], unit=(u.hourangle, u.deg))
        ra, dec = sc.ra.deg, sc.dec.deg
    out = {"ra": ra, "dec": dec, "flux": flux, "eflux": eflux}
    # TGSS ADR1 carries a per-source local-rms column; expose it so callers can derive a
    # field-local detection limit instead of assuming the survey-wide median (round-6 referee).
    if "Noise" in t.colnames:
        out["noise"] = np.asarray(t["Noise"], float)
    return out


def find_uss(
    low: dict[str, np.ndarray],
    high: dict[str, np.ndarray],
    *,
    low_survey: str = "tgss",
    high_survey: str = "nvss",
    radius_arcsec: float = 15.0,
) -> dict[str, np.ndarray]:
    """Cross-match a low- and high-frequency catalogue and compute per-source spectral indices.

    ``low``/``high`` are dicts with ``ra, dec, flux, eflux`` (mJy). Returns column arrays for the
    matched sources: ``ra, dec, s_lo, s_hi, alpha, e_alpha, sep, cls`` plus a boolean ``is_uss``.
    """
    nu_lo, nu_hi = SURVEYS[low_survey].freq_mhz, SURVEYS[high_survey].freq_mhz
    i_lo, i_hi, sep = crossmatch(low["ra"], low["dec"], high["ra"], high["dec"], radius_arcsec)
    s_lo, s_hi = low["flux"][i_lo], high["flux"][i_hi]
    e_lo, e_hi = low["eflux"][i_lo], high["eflux"][i_hi]
    alpha, e_alpha = spectral_index(s_lo, nu_lo, s_hi, nu_hi, e_lo, e_hi)
    cls = np.array([classify(a) for a in alpha])
    out = {
        "ra": low["ra"][i_lo],
        "dec": low["dec"][i_lo],
        "s_lo": s_lo,
        "s_hi": s_hi,
        "alpha": alpha,
        "e_alpha": e_alpha,
        "sep": sep,
        "cls": cls,
        "is_uss": alpha < USS_THRESHOLD,
    }
    if "alpha_true" in low:
        out["alpha_true"] = low["alpha_true"][i_lo]
    return out


def synthetic_field(
    n: int = 300,
    f_uss: float = 0.05,
    f_inverted: float = 0.05,
    seed: int | None = 0,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Generate a synthetic TGSS-like + NVSS-like field with known injected spectra (offline).

    Most sources are ordinary steep-spectrum (alpha ~ -0.8); a fraction are USS and a fraction
    inverted (alpha ~ +0.4). Half the injected USS sit well below the cut (alpha ~ -1.6) and half
    just below it (alpha ~ -1.35) --- an injection 3 sigma from the threshold exercises nothing at
    the boundary where the selection actually operates. The two catalogues share positions (with
    small astrometric jitter) so the cross-match recovers them; the low dict carries ``alpha_true``
    so the offline leg can score the cut against a known truth. Used by the tests and CI.
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(179.5, 180.5, n)
    dec = rng.uniform(29.5, 30.5, n)
    alpha_true = rng.normal(-0.8, 0.15, n)
    n_uss = int(round(f_uss * n))
    n_inv = int(round(f_inverted * n))
    n_deep = n_uss // 2
    alpha_true[:n_deep] = rng.normal(-1.6, 0.1, n_deep)
    alpha_true[n_deep:n_uss] = rng.normal(-1.35, 0.05, n_uss - n_deep)
    alpha_true[n_uss : n_uss + n_inv] = rng.normal(0.4, 0.1, n_inv)

    nu_lo, nu_hi = SURVEYS["tgss"].freq_mhz, SURVEYS["nvss"].freq_mhz
    s_lo = 10 ** rng.uniform(1.5, 3.0, n)  # 30 mJy .. 1 Jy at 150 MHz
    s_hi = s_lo * (nu_hi / nu_lo) ** alpha_true
    e_lo = 0.1 * s_lo
    e_hi = 0.05 * s_hi
    s_lo = rng.normal(s_lo, e_lo)
    s_hi = rng.normal(s_hi, e_hi)
    jit = 2.0 / 3600.0  # ~2" astrometric jitter
    low = {"ra": ra, "dec": dec, "flux": s_lo, "eflux": e_lo, "alpha_true": alpha_true}
    high = {
        "ra": ra + rng.normal(0, jit, n),
        "dec": dec + rng.normal(0, jit, n),
        "flux": s_hi,
        "eflux": e_hi,
    }
    return low, high


def annotate_known(
    ra: np.ndarray, dec: np.ndarray, radius_arcsec: float = 30.0
) -> list[dict]:  # pragma: no cover - network
    """Query NED at each position; return a known-classification record per source.

    For each candidate this returns ``{"name", "type", "z", "known_hzrg"}``. ``known_hzrg`` is True
    when NED already classifies the object as a galaxy/QSO with a measured redshift — i.e. it is
    *not* a novel candidate. Sources with no NED match (or no redshift) are the interesting ones.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.ipac.ned import Ned

    out = []
    for r, d in zip(ra, dec, strict=True):
        rec = {"name": "", "type": "", "z": float("nan"), "known_hzrg": False}
        try:
            t = Ned.query_region(SkyCoord(r * u.deg, d * u.deg), radius=radius_arcsec * u.arcsec)
            if t and len(t):
                row = t[0]
                try:
                    zval = float(row["Redshift"])
                except (TypeError, ValueError):
                    zval = float("nan")
                rec["name"] = str(row["Object Name"])
                rec["type"] = str(row["Type"])
                rec["z"] = zval
                # Positive galaxy/AGN classifications only. NED type "RadioS" (bare radio detection,
                # the default for an NVSS-only source) and "GClstr" are NOT an HzRG identification;
                # "RadioG" (radio galaxy) is. zval==zval is False only for NaN (no redshift).
                rec["known_hzrg"] = (zval == zval) and str(row["Type"]) in {"G", "QSO", "RadioG"}
        except Exception:  # noqa: BLE001 - NED outages must not crash the analysis
            pass
        out.append(rec)
    return out


def reference_spindex(
    ra: np.ndarray, dec: np.ndarray, radius_arcsec: float = 30.0
) -> list[dict]:  # pragma: no cover - network
    """Authoritative TGSS×NVSS spectral index from de Gasperin et al. (2018), per position.

    Looks each position up in the 1.4M-source catalogue (VizieR ``J/MNRAS/474/5008/spidxcat``) and
    returns ``{"spindex", "e_spindex", "sep"}`` for the nearest entry. NOTE: that catalogue is
    built from the SAME uncorrected TGSS ADR1 fluxes as this tool (its §4 bounds the flux-scale
    effect on its indices at ~0.06), so a raw-minus-reference difference does not measure the TGSS
    flux scale — it measures matching/deblending differences plus, for candidates selected on the
    extreme of a noisy alpha, selection bias (:func:`selection_bias_mc`). Prefer
    :func:`fetch_reference_cone` + :func:`reference_crossmatch` for a committable whole-field
    comparison; this per-position helper remains for spot checks.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    v = Vizier(columns=["*"])
    v.ROW_LIMIT = -1
    out = []
    for r, d in zip(np.atleast_1d(ra), np.atleast_1d(dec), strict=True):
        c = SkyCoord(r * u.deg, d * u.deg)
        rec = {"spindex": float("nan"), "e_spindex": float("nan"), "sep": float("nan")}
        try:
            t = v.query_region(
                c, radius=radius_arcsec * u.arcsec, catalog="J/MNRAS/474/5008/spidxcat"
            )
            if t and len(t[0]):
                tt = t[0]
                m = SkyCoord(tt["RAJ2000"], tt["DEJ2000"], unit="deg")
                sep = c.separation(m).arcsec
                k = int(np.argmin(sep))
                rec = {
                    "spindex": float(tt["SpIndex"][k]),
                    "e_spindex": float(tt["e_SpIndex"][k]),
                    "sep": float(sep[k]),
                }
        except Exception:  # noqa: BLE001 - VizieR outages must not crash the analysis
            pass
        out.append(rec)
    return out


def fetch_reference_cone(
    center, radius_deg: float
) -> dict[str, np.ndarray]:  # pragma: no cover - network
    """One cone query of the de Gasperin et al. (2018) spectral-index catalogue.

    Returns ``ra, dec, spindex, e_spindex, scode`` for every entry in the field. ``scode`` is the
    catalogue's own detection code: ``"S"`` for a two-survey detection, ``"L"`` for a limit (one
    survey undetected; ``e_SpIndex`` is 0 there and the index is not a measurement). One query for
    the whole field replaces per-source lookups, so the population comparison is committable.
    """
    from astropy import units as u
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "SpIndex", "e_SpIndex", "Scode"])
    v.ROW_LIMIT = -1
    res = v.query_region(center, radius=radius_deg * u.deg, catalog="J/MNRAS/474/5008/spidxcat")
    if not res:
        return {k: np.array([]) for k in ("ra", "dec", "spindex", "e_spindex", "scode")}
    t = res[0]
    return {
        "ra": np.asarray(t["RAJ2000"], float),
        "dec": np.asarray(t["DEJ2000"], float),
        "spindex": np.asarray(t["SpIndex"], float),
        "e_spindex": np.asarray(t["e_SpIndex"], float),
        "scode": np.array([str(s) for s in t["Scode"]]),
    }


def reference_crossmatch(
    res: dict[str, np.ndarray], ref: dict[str, np.ndarray], radius_arcsec: float = 15.0
) -> dict[str, np.ndarray]:
    """Nearest reference-catalogue entry for every matched source, aligned to ``res``.

    Returns ``ref_alpha, ref_e, ref_scode, ref_sep`` with ``nan``/``""`` where no reference entry
    lies within ``radius_arcsec``.
    """
    n = res["ra"].size
    out = {
        "ref_alpha": np.full(n, np.nan),
        "ref_e": np.full(n, np.nan),
        "ref_scode": np.array([""] * n, dtype=object),
        "ref_sep": np.full(n, np.nan),
    }
    if ref["ra"].size == 0 or n == 0:
        return out
    i_lo, i_hi, sep = crossmatch(res["ra"], res["dec"], ref["ra"], ref["dec"], radius_arcsec)
    out["ref_alpha"][i_lo] = ref["spindex"][i_hi]
    out["ref_e"][i_lo] = ref["e_spindex"][i_hi]
    out["ref_scode"][i_lo] = ref["scode"][i_hi]
    out["ref_sep"][i_lo] = sep
    return out


def population_offset(
    alpha: np.ndarray,
    e_alpha: np.ndarray,
    ref: dict[str, np.ndarray],
    dec: np.ndarray,
    *,
    dec_split: float = 30.0,
) -> dict:
    """Mean raw-minus-reference index offset over the whole matched population.

    Only genuine reference detections (``scode == "S"``) enter --- a limit row's index is not a
    measurement and its ``e = 0`` would be nonsense in any weight. Reports the mean offset with its
    standard error, the per-source scatter against the combined formal errors, and the same offset
    split at ``dec_split`` (the rescaled-TGSS coverage edge) --- the check that measures whether a
    declination-edge mechanism exists rather than asserting it.
    """
    ok = (
        np.isfinite(alpha)
        & np.isfinite(ref["ref_alpha"])
        & (np.asarray([str(s) for s in ref["ref_scode"]]) == "S")
    )
    d = alpha[ok] - ref["ref_alpha"][ok]
    n = int(ok.sum())
    if n < 2:
        return {"n_pairs": n}
    comb = np.sqrt(np.nanmean(e_alpha[ok] ** 2 + ref["ref_e"][ok] ** 2))
    lo = ok & (dec <= dec_split)
    hi = ok & (dec > dec_split)

    def _mean_se(mask: np.ndarray) -> tuple[float | None, float | None, int]:
        if mask.sum() < 2:
            return None, None, int(mask.sum())
        dd = alpha[mask] - ref["ref_alpha"][mask]
        return (
            round(float(np.mean(dd)), 4),
            round(float(np.std(dd, ddof=1) / np.sqrt(dd.size)), 4),
            int(mask.sum()),
        )

    m_lo, se_lo, n_lo = _mean_se(lo)
    m_hi, se_hi, n_hi = _mean_se(hi)
    return {
        "n_pairs": n,
        "mean_offset": round(float(np.mean(d)), 4),
        "se_offset": round(float(np.std(d, ddof=1) / np.sqrt(n)), 4),
        "scatter": round(float(np.std(d, ddof=1)), 4),
        "expected_scatter": round(float(comb), 4),
        "dec_le_mean": m_lo,
        "dec_le_se": se_lo,
        "dec_le_n": n_lo,
        "dec_gt_mean": m_hi,
        "dec_gt_se": se_hi,
        "dec_gt_n": n_hi,
    }


def selection_bias_mc(
    truth_alpha: np.ndarray,
    e_alpha: np.ndarray,
    *,
    threshold: float = USS_THRESHOLD,
    n_mc: int = 2000,
    seed: int = 0,
) -> dict:
    """The no-free-parameter model of selection on the noisy variable.

    Treat the reference indices as truth, add each source's own measurement noise, select
    ``noisy < threshold``, and record the mean (noisy - truth) of the selected set. That is the
    offset a *perfectly calibrated* survey pair shows for sources chosen on the extreme of the
    noisy measurement --- the number the six-candidate offset must be compared against before any
    of it is attributed to a flux scale. Deterministic (fixed seed).
    """
    ok = np.isfinite(truth_alpha) & np.isfinite(e_alpha) & (e_alpha > 0)
    t = truth_alpha[ok]
    e = e_alpha[ok]
    rng = np.random.default_rng(seed)
    means = []
    n_sel = []
    for _ in range(n_mc):
        noisy = t + rng.normal(0.0, e)
        sel = noisy < threshold
        if sel.any():
            means.append(float(np.mean(noisy[sel] - t[sel])))
            n_sel.append(int(sel.sum()))
    if not means:
        return {"n_realizations": 0}
    return {
        "n_realizations": len(means),
        "predicted_offset": round(float(np.mean(means)), 4),
        "predicted_offset_sd": round(float(np.std(means, ddof=1)), 4),
        "mean_n_selected": round(float(np.mean(n_sel)), 1),
    }


def matched_sensitivity(
    *,
    s_lim_lo_mjy: float = TGSS_LIMIT_MJY,
    s_lim_hi_mjy: float = NVSS_LIMIT_MJY,
    low_survey: str = "tgss",
    high_survey: str = "nvss",
    alpha_probe: float = -1.65,
) -> dict:
    """The truncation a joint two-survey detection imposes on steep spectra.

    A source is in the matched sample only if detected in BOTH surveys. The two floors are equally
    sensitive at ``alpha_match = ln(s_lim_hi/s_lim_lo)/ln(nu_hi/nu_lo)``; for anything steeper, the
    high-frequency survey is the limiting one and the required low-frequency flux grows as
    ``s_lim_hi * (nu_lo/nu_hi)^alpha`` --- so USS completeness is poor *by construction*, not by
    flux-scale error. ``alpha_probe`` reports the required 150 MHz flux at a specimen steep index.
    """
    nu_lo, nu_hi = SURVEYS[low_survey].freq_mhz, SURVEYS[high_survey].freq_mhz
    lnr = np.log(nu_hi / nu_lo)
    alpha_match = float(np.log(s_lim_hi_mjy / s_lim_lo_mjy) / lnr)
    s_lo_needed_at_thr = float(s_lim_hi_mjy * (nu_lo / nu_hi) ** USS_THRESHOLD)
    s_lo_needed_probe = float(s_lim_hi_mjy * (nu_lo / nu_hi) ** alpha_probe)
    return {
        "alpha_match": round(alpha_match, 3),
        "s150_needed_at_threshold_mjy": round(s_lo_needed_at_thr, 1),
        "alpha_probe": alpha_probe,
        "s150_needed_at_probe_mjy": round(s_lo_needed_probe, 1),
    }


def uss_confusion(
    res: dict[str, np.ndarray],
    ref: dict[str, np.ndarray],
    *,
    radius_arcsec: float = 15.0,
    threshold: float = USS_THRESHOLD,
) -> dict:
    """The raw USS cut scored against the reference catalogue's own USS population.

    Purity: of the raw-flagged candidates, how many the reference (non-limit entries) also has
    below the threshold. Completeness: of the reference's USS sources in the field, how many the
    raw cut recovers. Both directions matter --- a cut can look 'contaminated' while ALSO missing
    most of the real population, and this measures each against the same catalogue.
    """
    flagged = res["is_uss"] & np.isfinite(res["alpha"])
    det = ref["scode"] == "S"
    ref_uss = det & (ref["spindex"] < threshold)
    n_flagged = int(flagged.sum())
    n_ref_uss = int(ref_uss.sum())
    tp_purity = 0
    if n_flagged and det.any():
        i_lo, i_hi, _sep = crossmatch(
            res["ra"][flagged], res["dec"][flagged], ref["ra"][det], ref["dec"][det], radius_arcsec
        )
        tp_purity = int((ref["spindex"][det][i_hi] < threshold).sum())
    n_recovered = 0
    if n_ref_uss and n_flagged:
        j_lo, _j_hi, _s = crossmatch(
            ref["ra"][ref_uss],
            ref["dec"][ref_uss],
            res["ra"][flagged],
            res["dec"][flagged],
            radius_arcsec,
        )
        n_recovered = int(j_lo.size)
    return {
        "n_flagged": n_flagged,
        "n_flagged_ref_uss": tp_purity,
        "purity": round(tp_purity / n_flagged, 3) if n_flagged else None,
        "n_ref_uss": n_ref_uss,
        "n_ref_uss_recovered": n_recovered,
        "completeness": round(n_recovered / n_ref_uss, 3) if n_ref_uss else None,
    }


def chance_matches(n_lo: int, n_hi: int, area_deg2: float, radius_arcsec: float = 15.0) -> float:
    """Expected number of chance coincidences in the positional cross-match.

    ``n_lo`` targets each searched in a ``radius_arcsec`` disc against a background surface density
    ``n_hi / area_deg2`` --- the arithmetic the papers previously asserted as '~1-2' without
    computing.
    """
    disc_deg2 = np.pi * (radius_arcsec / 3600.0) ** 2
    return float(n_lo * n_hi * disc_deg2 / area_deg2)


def analyze(res: dict[str, np.ndarray], source: str = "unknown") -> dict:
    """Summarise a :func:`find_uss` result into a JSON-serialisable metrics dict."""
    a = res["alpha"]
    finite = np.isfinite(a)
    uss = res["is_uss"] & finite
    return {
        "source": source,
        "n_matched": int(finite.sum()),
        "n_uss": int(uss.sum()),
        "n_inverted": int(np.sum(res["cls"] == "inverted")),
        "alpha_median": float(np.median(a[finite])) if finite.any() else float("nan"),
        "alpha_min": float(np.min(a[finite])) if finite.any() else float("nan"),
        "uss_threshold": USS_THRESHOLD,
    }


def make_figures(
    res: dict[str, np.ndarray], out_dir, rc: dict[str, np.ndarray] | None = None
) -> list:
    """Spectral-index distribution + alpha-vs-flux scatter + (real runs) raw-vs-reference panel."""
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    a = res["alpha"]
    finite = np.isfinite(a)
    uss = res["is_uss"] & finite
    paths = []

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.hist(a[finite], bins=30, color="0.6")
    ax.axvline(USS_THRESHOLD, color="r", ls="--", label=f"USS ($\\alpha<{USS_THRESHOLD}$)")
    ax.set(
        xlabel=r"spectral index $\alpha_{150}^{1400}$",
        ylabel="sources",
        title="Two-point spectral indices",
    )
    ax.legend()
    p = out / "alpha_hist.pdf"
    fig.tight_layout()
    fig.savefig(p)
    plt.close(fig)
    paths.append(p)

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.scatter(res["s_lo"][finite & ~uss], a[finite & ~uss], s=8, color="0.6", label="sources")
    ax.scatter(res["s_lo"][uss], a[uss], s=24, color="r", marker="*", label="USS candidates")
    ax.axhline(USS_THRESHOLD, color="r", ls="--", lw=0.8)
    ax.set(
        xscale="log",
        xlabel="150 MHz flux (mJy)",
        ylabel=r"$\alpha_{150}^{1400}$",
        title="USS candidates vs flux",
    )
    ax.legend()
    p = out / "alpha_vs_flux.pdf"
    fig.tight_layout()
    fig.savefig(p)
    plt.close(fig)
    paths.append(p)

    if rc is not None:  # pragma: no cover - network run only
        ok = finite & np.isfinite(rc["ref_alpha"])
        det = ok & (np.asarray([str(s) for s in rc["ref_scode"]]) == "S")
        lim = ok & ~det
        fig, ax = plt.subplots(figsize=(5, 4.2))
        ax.scatter(rc["ref_alpha"][det & ~uss], a[det & ~uss], s=8, color="0.6", label="matched")
        ax.scatter(
            rc["ref_alpha"][det & uss],
            a[det & uss],
            s=40,
            color="r",
            marker="*",
            label="raw USS candidates",
        )
        if lim.any():
            ax.scatter(
                rc["ref_alpha"][lim],
                a[lim],
                s=26,
                facecolors="none",
                edgecolors="C1",
                label="reference limit (Scode=L)",
            )
        span = [-2.2, 0.8]
        ax.plot(span, span, color="k", lw=0.8, ls=":")
        ax.axhline(USS_THRESHOLD, color="r", ls="--", lw=0.8)
        ax.axvline(USS_THRESHOLD, color="r", ls="--", lw=0.8)
        ax.set(
            xlim=span,
            ylim=span,
            xlabel=r"reference $\alpha$ (de Gasperin et al. 2018)",
            ylabel=r"raw $\alpha$ (this work)",
            title="Raw vs. reference spectral index",
        )
        ax.legend(fontsize=7, loc="upper left")
        p = out / "alpha_compare.pdf"
        fig.tight_layout()
        fig.savefig(p)
        plt.close(fig)
        paths.append(p)
    return paths


def run(center=None, radius_deg: float = 2.0, out: str = ".", *, offline: bool = False) -> dict:
    """Full slice: fetch (or synthesise) a TGSS x NVSS field, find USS sources, write artifacts.

    Writes ``results/uss_metrics.json``, the figures, ``results/uss_candidates.csv``, and (real
    runs) the full-population reference comparison ``results/uss_reference_check.csv``. A failed
    fetch RAISES --- it must never silently substitute the synthetic field for real data. Returns
    the metrics dict.
    """
    import csv
    from pathlib import Path

    ref: dict[str, np.ndarray] | None = None
    if offline or center is None:
        low, high = synthetic_field()
        source = "synthetic"
    else:  # pragma: no cover - network
        low = fetch_survey(center, radius_deg, "tgss")
        high = fetch_survey(center, radius_deg, "nvss")
        if low["ra"].size == 0 or high["ra"].size == 0:
            raise RuntimeError("TGSS/NVSS fetch returned no sources; refusing a synthetic fallback")
        source = f"tgss-x-nvss @ ({center.ra.deg:.2f}, {center.dec.deg:+.2f}) r={radius_deg:g} deg"
        ref = fetch_reference_cone(center, radius_deg)
    res = find_uss(low, high)
    metrics = analyze(res, source)
    metrics["radius_deg"] = radius_deg
    if center is not None:  # pragma: no cover - network
        metrics["field_ra"] = round(float(center.ra.deg), 2)
        metrics["field_dec"] = round(float(center.dec.deg), 2)
    area = np.pi * radius_deg**2
    metrics["chance_matches_expected"] = round(
        chance_matches(int(low["ra"].size), int(high["ra"].size), area), 2
    )
    metrics["matched_sensitivity"] = matched_sensitivity(
        alpha_probe=round(float(np.nanmin(res["alpha"])), 2)
    )
    ms = metrics["matched_sensitivity"]
    s_need = ms["s150_needed_at_probe_mjy"]
    metrics["frac_sample_above_s150_needed"] = round(
        float(np.mean(res["s_lo"][np.isfinite(res["alpha"])] >= s_need)), 3
    )

    rc: dict[str, np.ndarray] | None = None
    if ref is not None:  # pragma: no cover - network
        rc = reference_crossmatch(res, ref)
        metrics["population_offset"] = population_offset(
            res["alpha"], res["e_alpha"], rc, res["dec"]
        )
        det = np.asarray([str(s) for s in rc["ref_scode"]]) == "S"
        metrics["selection_bias_mc"] = selection_bias_mc(rc["ref_alpha"][det], res["e_alpha"][det])
        metrics["uss_confusion"] = uss_confusion(res, ref)
        flagged = res["is_uss"] & np.isfinite(res["alpha"])
        sel = flagged & np.isfinite(rc["ref_alpha"])
        sel_det = sel & det
        metrics["flagged_offset_all"] = (
            round(float(np.mean(res["alpha"][sel] - rc["ref_alpha"][sel])), 3)
            if sel.any()
            else None
        )
        metrics["flagged_offset_detections"] = (
            round(float(np.mean(res["alpha"][sel_det] - rc["ref_alpha"][sel_det])), 3)
            if sel_det.any()
            else None
        )
        metrics["n_flagged_ref_limit"] = int((sel & ~det).sum())
        _write_reference_check(Path(out) / "results" / "uss_reference_check.csv", res, rc)
    # the fixture's own boundary measurement: score the cut and the selection model against the
    # injected truth. Deterministic, and computed in BOTH modes so a real run emits the Syn
    # macros the Methods section cites (the two-namespace pattern; see CLAUDE.md).
    syn_res = res if "alpha_true" in res else find_uss(*synthetic_field())
    truth_uss = syn_res["alpha_true"] < USS_THRESHOLD
    syn_flagged = syn_res["is_uss"] & np.isfinite(syn_res["alpha"])
    tp = int((syn_flagged & truth_uss).sum())
    metrics["syn_cut_purity"] = round(tp / max(int(syn_flagged.sum()), 1), 3)
    metrics["syn_cut_completeness"] = round(tp / max(int(truth_uss.sum()), 1), 3)
    metrics["syn_selection_bias_mc"] = selection_bias_mc(syn_res["alpha_true"], syn_res["e_alpha"])
    metrics["syn_flagged_offset"] = (
        round(float(np.mean(syn_res["alpha"][syn_flagged] - syn_res["alpha_true"][syn_flagged])), 3)
        if syn_flagged.any()
        else None
    )

    outp = Path(out)
    (outp / "results").mkdir(parents=True, exist_ok=True)
    from .report import _results_are_real, write_results

    # a synthetic run must not clobber the real CSV/figures (the JSON and macros carry their own
    # merge guards; the CSV and PDFs do not)
    json_path = outp / "results" / "uss_metrics.json"
    write_artifacts = True
    if source == "synthetic":
        try:
            import json as _json

            write_artifacts = not (
                json_path.is_file() and _results_are_real(_json.loads(json_path.read_text()))
            )
        except Exception:
            write_artifacts = True

    write_results(metrics, json_path)
    if write_artifacts:
        make_figures(res, outp / "papers" / "spectra" / "figures", rc)
    _write_macros(metrics, outp / "papers" / "spectra" / "generated" / "macros.tex")
    if write_artifacts:
        # the candidate table (USS first), for the write-up + follow-up
        order = np.argsort(res["alpha"])
        with open(outp / "results" / "uss_candidates.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(
                ["ra", "dec", "s150_mjy", "s1400_mjy", "alpha", "e_alpha", "sep_arcsec", "class"]
            )
            for i in order:
                if not np.isfinite(res["alpha"][i]):
                    continue
                w.writerow(
                    [
                        f"{res['ra'][i]:.5f}",
                        f"{res['dec'][i]:.5f}",
                        f"{res['s_lo'][i]:.2f}",
                        f"{res['s_hi'][i]:.2f}",
                        f"{res['alpha'][i]:.2f}",
                        f"{res['e_alpha'][i]:.2f}",
                        f"{res['sep'][i]:.1f}",
                        res["cls"][i],
                    ]
                )
    return metrics


def _write_reference_check(path, res, rc) -> None:  # pragma: no cover - network run only
    """Commit the full raw-vs-reference comparison, one row per matched source."""
    import csv
    from pathlib import Path

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "ra",
                "dec",
                "alpha_raw",
                "e_alpha",
                "ref_alpha",
                "ref_e",
                "ref_scode",
                "ref_sep",
                "flagged_uss",
            ]
        )
        for i in range(res["ra"].size):
            if not np.isfinite(res["alpha"][i]):
                continue
            w.writerow(
                [
                    f"{res['ra'][i]:.5f}",
                    f"{res['dec'][i]:.5f}",
                    f"{res['alpha'][i]:.3f}",
                    f"{res['e_alpha'][i]:.3f}",
                    f"{rc['ref_alpha'][i]:.3f}" if np.isfinite(rc["ref_alpha"][i]) else "",
                    f"{rc['ref_e'][i]:.3f}" if np.isfinite(rc["ref_e"][i]) else "",
                    str(rc["ref_scode"][i]),
                    f"{rc['ref_sep'][i]:.1f}" if np.isfinite(rc["ref_sep"][i]) else "",
                    int(bool(res["is_uss"][i])),
                ]
            )


def _write_macros(m: dict, path) -> None:
    """Emit LaTeX ``\\newcommand`` macros so the paper hard-codes no number."""
    from pathlib import Path

    po = m.get("population_offset") or {}
    sb = m.get("selection_bias_mc") or {}
    cf = m.get("uss_confusion") or {}
    ms = m.get("matched_sensitivity") or {}

    def _g(dic: dict, key: str) -> str:
        val = dic.get(key)
        return "--" if val is None else str(val)

    def _m(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    def _pct(dic: dict, key: str) -> str:
        val = dic.get(key)
        return "--" if val is None else f"{100 * float(val):.0f}"

    lines = [
        "% Auto-generated by jansky_research.spectra._write_macros — do not edit by hand.",
        rf"\newcommand{{\usSource}}{{{m['source']}}}",
        rf"\newcommand{{\usNmatched}}{{{m['n_matched']}}}",
        rf"\newcommand{{\usAlphaMedian}}{{{m['alpha_median']:.2f}}}",
        rf"\newcommand{{\usAlphaMin}}{{{m['alpha_min']:.2f}}}",
        rf"\newcommand{{\usNuss}}{{{m['n_uss']}}}",
        rf"\newcommand{{\usNinverted}}{{{m['n_inverted']}}}",
        rf"\newcommand{{\usThreshold}}{{{m['uss_threshold']:.1f}}}",
        rf"\newcommand{{\usFieldRa}}{{{_m('field_ra')}}}",
        rf"\newcommand{{\usFieldDec}}{{{_m('field_dec')}}}",
        rf"\newcommand{{\usRadius}}{{{_m('radius_deg')}}}",
        rf"\newcommand{{\usChanceMatches}}{{{_m('chance_matches_expected')}}}",
        # the full-population reference comparison (real runs)
        rf"\newcommand{{\usRealNPairs}}{{{_g(po, 'n_pairs')}}}",
        rf"\newcommand{{\usRealDeltaAlpha}}{{{_g(po, 'mean_offset')}}}",
        rf"\newcommand{{\usRealDeltaAlphaSE}}{{{_g(po, 'se_offset')}}}",
        rf"\newcommand{{\usRealScatter}}{{{_g(po, 'scatter')}}}",
        rf"\newcommand{{\usRealScatterExp}}{{{_g(po, 'expected_scatter')}}}",
        rf"\newcommand{{\usRealDecLeMean}}{{{_g(po, 'dec_le_mean')}}}",
        rf"\newcommand{{\usRealDecLeSE}}{{{_g(po, 'dec_le_se')}}}",
        rf"\newcommand{{\usRealDecGtMean}}{{{_g(po, 'dec_gt_mean')}}}",
        rf"\newcommand{{\usRealDecGtSE}}{{{_g(po, 'dec_gt_se')}}}",
        rf"\newcommand{{\usRealFlaggedOffset}}{{{_m('flagged_offset_detections')}}}",
        rf"\newcommand{{\usRealFlaggedOffsetAll}}{{{_m('flagged_offset_all')}}}",
        rf"\newcommand{{\usRealNFlaggedLimit}}{{{_m('n_flagged_ref_limit')}}}",
        rf"\newcommand{{\usRealSelPred}}{{{_g(sb, 'predicted_offset')}}}",
        rf"\newcommand{{\usRealSelPredSD}}{{{_g(sb, 'predicted_offset_sd')}}}",
        # the confusion matrix against the reference's own USS population
        rf"\newcommand{{\usRealNRefUss}}{{{_g(cf, 'n_ref_uss')}}}",
        rf"\newcommand{{\usRealNRefUssRec}}{{{_g(cf, 'n_ref_uss_recovered')}}}",
        rf"\newcommand{{\usRealNFlaggedRefUss}}{{{_g(cf, 'n_flagged_ref_uss')}}}",
        rf"\newcommand{{\usRealPurityPct}}{{{_pct(cf, 'purity')}}}",
        rf"\newcommand{{\usRealCompletenessPct}}{{{_pct(cf, 'completeness')}}}",
        # matched-sensitivity truncation (pure arithmetic from documented survey floors)
        rf"\newcommand{{\usAlphaMatch}}{{{_g(ms, 'alpha_match')}}}",
        rf"\newcommand{{\usSNeededThr}}{{{_g(ms, 's150_needed_at_threshold_mjy')}}}",
        rf"\newcommand{{\usAlphaProbe}}{{{_g(ms, 'alpha_probe')}}}",
        rf"\newcommand{{\usSNeededProbe}}{{{_g(ms, 's150_needed_at_probe_mjy')}}}",
        rf"\newcommand{{\usFracAboveSNeeded}}{{{_m('frac_sample_above_s150_needed')}}}",
        # the offline fixture's own boundary measurement (synthetic namespace)
        rf"\newcommand{{\usSynCutPurity}}{{{_m('syn_cut_purity')}}}",
        rf"\newcommand{{\usSynCutCompleteness}}{{{_m('syn_cut_completeness')}}}",
        rf"\newcommand{{\usSynFlaggedOffset}}{{{_m('syn_flagged_offset')}}}",
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

    from astropy.coordinates import SkyCoord

    p = argparse.ArgumentParser(description="Ultra-steep-spectrum source hunt (TGSS x NVSS).")
    p.add_argument("--ra", type=float, help="field centre RA (deg)")
    p.add_argument("--dec", type=float, help="field centre Dec (deg)")
    p.add_argument("--radius", type=float, default=2.0, help="cone radius (deg)")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    center = None if (args.offline or args.ra is None) else SkyCoord(args.ra, args.dec, unit="deg")
    print(json.dumps(run(center, args.radius, args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
