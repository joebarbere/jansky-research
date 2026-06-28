"""Southern peaked-spectrum (GPS/CSS) selection via GLEAM-X + RACS multi-band spectral curvature.

The northern :mod:`jansky_research.peaked` slice could only *bound* the low-frequency spectral index
(TGSS is shallow, so 150 MHz is an upper limit). In the south, GLEAM-X DR2 (Ross et al. 2024) measures
each source in **20 in-band sub-bands over 72--231 MHz**, so the low-frequency shape is *measured*.
Adding the three RACS bands (887.5, 1367.5, 1655.5 MHz) gives up to ~23 flux points over a factor of
~23 in frequency --- enough to fit a real log-parabola SED and **measure the turnover frequency**
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


def classify_curved(fit: dict, alpha_lo: float, alpha_hi: float, *, rise_min: float = -0.1) -> str:
    """Classify a southern SED: peaked / uss / steep / flat / inverted (reuses ``peaked.classify_sed``).

    ``peaked`` requires *both* a concave in-band log-parabola turnover **and** an optically-thick
    rising low-frequency side ($\\alpha_\\mathrm{lo}>$ ``rise_min``, the GPS/SSA signature) --- the
    latter rejects ordinary steep sources whose mild low-frequency flattening fits a shallow parabola
    peaking spuriously at the GLEAM band edge. Otherwise ``uss`` when both two-point indices are
    ultra-steep ($<$ ``USS_THRESHOLD``; candidate high-z radio galaxy), else the coarse two-index class.
    """
    from .peaked import classify_sed

    if not (np.isfinite(alpha_lo) and np.isfinite(alpha_hi)):
        return "nan"
    if fit.get("is_peaked") and alpha_lo > rise_min:
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

    nu_pk, cls, is_pk, is_uss, a_lo, a_hi = [], [], [], [], [], []
    for k in range(ig.size):
        gi, rj = int(ig[k]), int(ir[k])
        flux = np.concatenate([gflux[gi], rflux[rj]])
        eflux = np.concatenate([geflux[gi], reflux[rj]])
        fit = fit_log_parabola(nu_all, flux, eflux)
        # measured low-freq (GLEAM-X) and high-freq (RACS) two-point indices
        alpha_lo, _ = spectral_index(
            gflux[gi, 0], GLEAMX_NU_GHZ[0], gflux[gi, -1], GLEAMX_NU_GHZ[-1]
        )
        alpha_hi, _ = spectral_index(rflux[rj, 0], RACS_NU_GHZ[0], rflux[rj, -1], RACS_NU_GHZ[-1])
        c = classify_curved(fit, float(alpha_lo), float(alpha_hi))
        if not compact_arr[gi] and c in ("peaked", "inverted"):
            c = "extended"  # GLEAM-extended -> RACS flux loss fakes a turnover; reject
        nu_pk.append(fit["nu_pk_ghz"])
        cls.append(c)
        is_pk.append(c == "peaked")
        is_uss.append(c == "uss")
        a_lo.append(float(alpha_lo))
        a_hi.append(float(alpha_hi))
    return {
        "ra": gra[ig],
        "dec": gdec[ig],
        "nu_pk_ghz": np.asarray(nu_pk),
        "alpha_lo": np.asarray(a_lo),
        "alpha_hi": np.asarray(a_hi),
        "cls": np.asarray(cls, dtype=object).astype(str),
        "is_peaked": np.asarray(is_pk, bool),
        "is_uss": np.asarray(is_uss, bool),
    }


def synthetic_field(
    n_sources: int = 1200,
    *,
    peaked_fraction: float = 0.05,
    uss_fraction: float = 0.05,
    rel_err: float = 0.08,
    seed: int = 0,
) -> tuple[dict, dict, np.ndarray, np.ndarray]:
    """Synthetic GLEAM-X(20-band)+RACS(3-band) catalogues with injected peaked/uss/steep/flat SEDs."""
    rng = np.random.default_rng(seed)
    ra = rng.uniform(20.0, 25.0, n_sources)
    dec = rng.uniform(-40.0, -35.0, n_sources)
    s_ref = 10.0 ** rng.uniform(1.0, 2.5, n_sources)  # ~10-300 mJy at 200 MHz
    nu_ref = 0.2  # GHz
    is_peaked = rng.random(n_sources) < peaked_fraction
    is_uss = (~is_peaked) & (rng.random(n_sources) < uss_fraction / (1 - peaked_fraction))
    alpha = rng.uniform(-0.9, -0.5, n_sources)  # ordinary steep/flat
    alpha[is_uss] = rng.uniform(-1.6, -1.3, int(is_uss.sum()))  # ultra-steep
    nu_pk = rng.uniform(0.3, 0.7, n_sources)  # injected turnover (GHz) for peaked sources
    curv = rng.uniform(0.6, 1.4, n_sources)

    def sed(nu):
        out = np.empty((n_sources, nu.size))
        for i in range(n_sources):
            if is_peaked[i]:
                out[i] = 10.0 ** (np.log10(s_ref[i]) - curv[i] * (np.log10(nu / nu_pk[i])) ** 2)
            else:
                out[i] = s_ref[i] * (nu / nu_ref) ** alpha[i]
        return out

    gflux = sed(GLEAMX_NU_GHZ) * rng.normal(1.0, rel_err, (n_sources, GLEAMX_NU_GHZ.size))
    rflux = sed(RACS_NU_GHZ) * rng.normal(1.0, rel_err, (n_sources, RACS_NU_GHZ.size))
    gflux = np.clip(gflux, 1e-3, None)
    rflux = np.clip(rflux, 1e-3, None)
    jit = lambda: rng.normal(0.0, 2.0 / 3600.0, n_sources)  # noqa: E731  (~2" jitter)
    gleamx = {"ra": ra + jit(), "dec": dec + jit(), "flux": gflux, "eflux": rel_err * gflux}
    racs = {"ra": ra + jit(), "dec": dec + jit(), "flux": rflux, "eflux": rel_err * rflux}
    return gleamx, racs, is_peaked, is_uss


def fetch_gleamx(
    center, radius_deg: float, *, rel_err: float = 0.1
) -> dict:  # pragma: no cover - network
    """Cone-search GLEAM-X DR2 (VizieR ``VIII/113``); 19 in-band sub-band fluxes, **Jy → mJy**.

    GLEAM-X DR2 stores integrated fluxes in **Jy** (RACS is in mJy) — converted here ×1000 to a common
    mJy scale. DR2 exposes no per-band errors via VizieR, so a nominal ``rel_err`` per band is used for
    the fit weights (an honest limitation). Returns ``ra``/``dec`` and ``flux``/``eflux`` of shape
    ``(n, len(GLEAMX_NU_GHZ))`` in the ``GLEAMX_NU_GHZ`` band order.
    """
    import numpy as _np
    from astropy import units as _u
    from astroquery.vizier import Vizier

    nb = len(_GLEAM_BANDS_MHZ)
    cols = ["RAJ2000", "DEJ2000", "Fintwide", "Fpwide", *[f"Fint{b:03d}" for b in _GLEAM_BANDS_MHZ]]
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
    flux = _np.where(_np.isfinite(flux) & (flux > 0), flux, _np.nan)
    # compactness at GLEAM's ~2' resolution: integrated/peak ~ 1 for point sources. Non-compact (>1.2)
    # sources lose flux at RACS's 15-25" resolution, faking a high-frequency turnover (resolution
    # artefact) -- exclude them, the southern analogue of the peaked slice's alpha_high floor.
    fint_w = _np.asarray(t["Fintwide"], float)
    fp_w = _np.asarray(t["Fpwide"], float)
    with _np.errstate(divide="ignore", invalid="ignore"):
        ratio = fint_w / fp_w
    compact = _np.isfinite(ratio) & (ratio < 1.2)
    return {
        "ra": _np.asarray(t["RAJ2000"], float),
        "dec": _np.asarray(t["DEJ2000"], float),
        "flux": flux,
        "eflux": rel_err * flux,
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


def run(center=None, radius_deg: float = 3.0, out: str = ".", *, offline: bool = True) -> dict:
    """Full slice: synthesise (or fetch) GLEAM-X×RACS, fit curvature, classify, write artifacts."""
    import json
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
    metrics = {
        "source": source,
        "n_matched": int(cls.size),
        "n_peaked": int(res["is_peaked"].sum()),
        "n_uss": int(res["is_uss"].sum()),
        "n_extended": int(np.sum(cls == "extended")),
        "median_nu_pk_mhz": (
            round(1e3 * float(np.nanmedian(res["nu_pk_ghz"][res["is_peaked"]])), 1)
            if res["is_peaked"].any()
            else 0.0
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

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "southern_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(gleamx, racs, res, op / "papers" / "southern" / "figures")
    _write_macros(metrics, op / "papers" / "southern" / "generated" / "macros.tex")
    return metrics


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
    from pathlib import Path

    lines = [
        "% Auto-generated by jansky_research.southern._write_macros — do not edit by hand.",
        rf"\newcommand{{\soSource}}{{{m['source']}}}",
        rf"\newcommand{{\soNmatched}}{{{m['n_matched']}}}",
        rf"\newcommand{{\soNpeaked}}{{{m['n_peaked']}}}",
        rf"\newcommand{{\soNuss}}{{{m['n_uss']}}}",
        rf"\newcommand{{\soNextended}}{{{m.get('n_extended', 0)}}}",
        rf"\newcommand{{\soMedianNupk}}{{{m['median_nu_pk_mhz']}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


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
    args = p.parse_args(argv)
    center = None if (args.offline or args.ra is None) else SkyCoord(args.ra, args.dec, unit="deg")
    print(json.dumps(run(center, args.radius, args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
