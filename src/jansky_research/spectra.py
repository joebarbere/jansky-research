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
- TGSS ADR1 has a position-dependent flux-scale systematic (Hurley-Walker 2017, arXiv:1703.06635;
  ~15% typical, up to ~40-50% in places) that can bias alpha by ~0.1-0.2; the conservative -1.3 cut
  limits class flips. The rescaled TGSS-RSADR1 covers only Dec <= +30deg, so fields above that use
  the uncorrected ADR1 scale and may carry a larger, unquantified bias.
- TGSS (25"), NVSS (45") and VLASS (2.5") differ in resolution; alpha from NVSS x TGSS (both recover
  extended flux) is the primary estimate, VLASS (QL Epoch 1, Gordon et al. 2021) a compact-source
  curvature check.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "SURVEYS",
    "Survey",
    "USS_THRESHOLD",
    "classify",
    "crossmatch",
    "fetch_survey",
    "find_uss",
    "reference_spindex",
    "spectral_index",
    "synthetic_field",
]

USS_THRESHOLD = -1.3  # alpha below this = ultra-steep-spectrum (high-z radio-galaxy candidate)


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
    return {"ra": ra, "dec": dec, "flux": flux, "eflux": eflux}


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
    return {
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


def synthetic_field(
    n: int = 300,
    f_uss: float = 0.05,
    f_inverted: float = 0.05,
    seed: int | None = 0,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Generate a synthetic TGSS-like + NVSS-like field with known injected spectra (offline).

    Most sources are ordinary steep-spectrum (alpha ~ -0.8); a fraction are USS (alpha ~ -1.6) and
    a fraction inverted (alpha ~ +0.4). The two catalogues share positions (with small astrometric
    jitter) so the cross-match recovers them. Used by the tests and as the offline fallback.
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(179.5, 180.5, n)
    dec = rng.uniform(29.5, 30.5, n)
    alpha_true = rng.normal(-0.8, 0.15, n)
    n_uss = int(round(f_uss * n))
    n_inv = int(round(f_inverted * n))
    alpha_true[:n_uss] = rng.normal(-1.6, 0.1, n_uss)
    alpha_true[n_uss : n_uss + n_inv] = rng.normal(0.4, 0.1, n_inv)

    nu_lo, nu_hi = SURVEYS["tgss"].freq_mhz, SURVEYS["nvss"].freq_mhz
    s_lo = 10 ** rng.uniform(1.5, 3.0, n)  # 30 mJy .. 1 Jy at 150 MHz
    s_hi = s_lo * (nu_hi / nu_lo) ** alpha_true
    e_lo = 0.1 * s_lo
    e_hi = 0.05 * s_hi
    s_lo = rng.normal(s_lo, e_lo)
    s_hi = rng.normal(s_hi, e_hi)
    jit = 2.0 / 3600.0  # ~2" astrometric jitter
    low = {"ra": ra, "dec": dec, "flux": s_lo, "eflux": e_lo}
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
    """Authoritative TGSS×NVSS spectral index from de Gasperin et al. (2018).

    Looks each position up in the **flux-scale-corrected** 1.4M-source catalogue
    (VizieR ``J/MNRAS/474/5008/spidxcat``) and returns ``{"spindex", "e_spindex", "sep"}`` for the
    nearest entry. Compare against this tool's raw-TGSS ``alpha``: a large negative offset means the
    uncorrected TGSS ADR1 flux scale has *inflated* the index and the USS flag is likely spurious —
    the mandatory validation before any USS/HzRG claim (this is exactly what reduced a 6-candidate
    raw list to 1 confirmed USS source in the NGP test field).
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


def make_figures(res: dict[str, np.ndarray], out_dir) -> list:
    """Spectral-index distribution + alpha-vs-flux scatter highlighting USS sources."""
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
    return paths


def run(center=None, radius_deg: float = 2.0, out: str = ".", *, offline: bool = False) -> dict:
    """Full slice: fetch (or synthesise) a TGSS x NVSS field, find USS sources, write artifacts.

    Writes ``results/uss_metrics.json``, the figures, and a ``results/uss_candidates.csv`` table.
    Returns the metrics dict.
    """
    import csv
    import json
    from pathlib import Path

    if offline or center is None:
        low, high = synthetic_field()
        source = "synthetic"
    else:  # pragma: no cover - network
        try:
            low = fetch_survey(center, radius_deg, "tgss")
            high = fetch_survey(center, radius_deg, "nvss")
            source = f"tgss-x-nvss @ {center}"
        except Exception:
            low, high = synthetic_field()
            source = "synthetic"
    res = find_uss(low, high)
    metrics = analyze(res, source)

    outp = Path(out)
    (outp / "results").mkdir(parents=True, exist_ok=True)
    (outp / "results" / "uss_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    make_figures(res, outp / "paper" / "figures")
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
