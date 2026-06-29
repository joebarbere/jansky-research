"""Multi-decade VLBI flux variability of compact AGN, from the Astrogeo database.

The Astrogeo VLBI image database (Petrov; astrogeo.org) holds decades of dual-band **S/X**
(2.3 / 8.4 GHz) observations of ~21k compact sources. Each source's per-session total flux densities
form a **multi-decade, parsec-scale light curve**, to which the standard transient-survey variability
statistics apply --- the same $\\eta$ (weighted reduced $\\chi^2$) and $V$ (coefficient of variation)
we built and tested for the VLASS three-epoch slice. The dual band additionally gives a per-source
**S/X spectral index**, so a source is characterised as variable *and* by its spectrum.

This module composes the tested helpers --- ``vlass.variability_metrics`` / ``vlass.select_candidates``
and ``spectra.spectral_index`` --- and adds a synthetic offline fixture plus the Astrogeo fetch. Pure
NumPy; the real fetch is network-gated. The honest caveat: a VLBI total flux density depends on the
session's ``(u,v)`` coverage and resolved-out flux, so apparent variability can be structural --- hence
a minimum-epoch gate and a literature recover-a-known validation before any source is called variable.
"""

from __future__ import annotations

import numpy as np

from . import spectra, vlass

__all__ = [
    "NU_S_GHZ",
    "NU_X_GHZ",
    "fetch_astrogeo",
    "lightcurve_metrics",
    "run",
    "select_variable",
    "sx_index",
    "synthetic_lightcurves",
    "variability_floor",
]

NU_S_GHZ = 2.3  # Astrogeo S band
NU_X_GHZ = 8.4  # Astrogeo X band
MIN_EPOCHS = 4  # a light curve needs at least this many finite epochs to be tested

# A curated *validation* set (not a blind survey): well-known, well-observed compact AGN whose
# variability is documented, so the run is a recover-a-known. Most are Doppler-boosted blazars expected
# to vary strongly; the four CSOs (compact symmetric objects) lack a boosted core and serve as steady
# negative controls. Caveat: OQ 208 and 2021+614 are documented as atypically variable for CSOs (Wu et
# al. 2013; Taylor et al. 2000), so the control floor is checked for sensitivity to them in the paper.
# The cleaner steady controls are 0108+388 and NGC 3894. J2000 name -> common name.
VALIDATION_SOURCES: dict[str, str] = {
    "J2202+4216": "BL Lac",
    "J0854+2006": "OJ 287",
    "J2253+1608": "3C 454.3",
    "J1256-0547": "3C 279",
    "J2232+1143": "CTA 102",
    "J0238+1636": "AO 0235+164",
    "J1512-0905": "PKS 1510-089",
    "J1224+2122": "4C 21.35",
    "J0841+7053": "4C 71.07",
    "J0006-0623": "PKS 0003-066",
    "J0433+0521": "3C 120",
    "J0319+4130": "3C 84",
    "J1642+3948": "3C 345",
    "J1229+0203": "3C 273",
    "J1407+2827": "OQ 208 (CSO)",
    "J2022+6136": "2021+614 (CSO)",
    "J0111+3906": "0108+388 (CSO)",
    "J1148+5924": "NGC 3894 (CSO)",
}


def lightcurve_metrics(
    fmat: np.ndarray, emat: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-source variability metrics from a ``(n_sources, n_epochs)`` flux/error matrix.

    Each row is one source's light curve with ``nan`` for sessions in which it was not measured. For
    every row with at least :data:`MIN_EPOCHS` finite points we compute, via the tested
    ``vlass.variability_metrics``, the significance $\\eta$ (weighted reduced $\\chi^2$), the amplitude
    $V$ (coefficient of variation), the $\\chi^2$ p-value, the epoch count, and the mean flux. Rows with
    too few epochs get ``nan`` metrics (and ``n_epochs`` counts the finite points regardless).
    """
    f = np.asarray(fmat, float)
    e = np.asarray(emat, float)
    n = f.shape[0]
    eta = np.full(n, np.nan)
    v = np.full(n, np.nan)
    pval = np.full(n, np.nan)
    nep = np.zeros(n, dtype=int)
    mean = np.full(n, np.nan)
    for i in range(n):
        ok = np.isfinite(f[i]) & np.isfinite(e[i]) & (e[i] > 0)
        nep[i] = int(ok.sum())
        if nep[i] < MIN_EPOCHS:
            continue
        m = vlass.variability_metrics(f[i, ok], e[i, ok])
        eta[i], v[i], pval[i], mean[i] = m.eta, m.v, m.p_value, m.mean_flux
    return eta, v, pval, nep, mean


def sx_index(
    flux_s: np.ndarray, flux_x: np.ndarray, e_s: np.ndarray, e_x: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Mean S/X two-point spectral index per source (reuses ``spectra.spectral_index``).

    ``flux_s`` / ``flux_x`` are ``(n_sources, n_epochs)`` matrices; we average each band over its
    finite epochs (the time-mean flux density) and take the two-point index between :data:`NU_S_GHZ`
    and :data:`NU_X_GHZ`. Sources lacking a finite mean in either band get ``nan``.
    """
    s = np.nanmean(np.where(np.isfinite(flux_s), flux_s, np.nan), axis=1)
    x = np.nanmean(np.where(np.isfinite(flux_x), flux_x, np.nan), axis=1)
    es = np.nanmean(np.where(np.isfinite(e_s), e_s, np.nan), axis=1)
    ex = np.nanmean(np.where(np.isfinite(e_x), e_x, np.nan), axis=1)
    good = np.isfinite(s) & np.isfinite(x) & (s > 0) & (x > 0)
    alpha = np.full(s.shape, np.nan)
    aerr = np.full(s.shape, np.nan)
    if good.any():
        a, ae = spectra.spectral_index(s[good], NU_S_GHZ, x[good], NU_X_GHZ, es[good], ex[good])
        alpha[good] = a
        aerr[good] = ae
    return alpha, aerr


def select_variable(
    eta: np.ndarray,
    v: np.ndarray,
    n_epochs: np.ndarray,
    *,
    min_epochs: int = MIN_EPOCHS,
    sigma: float = 3.0,
) -> tuple[np.ndarray, float, float]:
    """Variable candidates: 2-D log-$\\eta$/log-$V$ outliers (``vlass.select_candidates``) with enough epochs.

    Sources with fewer than ``min_epochs`` finite points are excluded before the cut is computed, so
    short light curves cannot define or pass the threshold. Returns ``(mask, eta_thr, v_thr)`` aligned
    to the input length.
    """
    eta = np.asarray(eta, float)
    v = np.asarray(v, float)
    nep = np.asarray(n_epochs)
    testable = (nep >= min_epochs) & np.isfinite(eta) & np.isfinite(v) & (eta > 0) & (v > 0)
    mask = np.zeros(eta.shape, dtype=bool)
    eta_thr = v_thr = float("nan")
    if testable.sum() >= 2:
        sub, eta_thr, v_thr = vlass.select_candidates(eta[testable], v[testable], sigma=sigma)
        mask[testable] = sub
    return mask, eta_thr, v_thr


def variability_floor(
    v: np.ndarray, n_epochs: np.ndarray, is_control: np.ndarray, *, min_epochs: int = MIN_EPOCHS
) -> tuple[float, np.ndarray]:
    """Empirical amplitude-variability floor set by intrinsically steady control sources.

    For VLBI total flux density, the per-session $V$ (coefficient of variation) of a genuinely steady
    source is *not* zero: it is set by amplitude-calibration scatter and by $(u,v)$-coverage /
    resolved-structure differences between sessions. Compact symmetric objects (CSOs), which lack a
    Doppler-boosted core, are such steady controls. Their median $V$ is therefore the floor below which
    $V$ is consistent with non-variability; testable non-control sources with $V$ **above** the floor
    are the amplitude-selected variables. Returns ``(floor, mask_above)`` aligned to the input length.
    """
    v = np.asarray(v, float)
    nep = np.asarray(n_epochs)
    ctrl = np.asarray(is_control, bool)
    okc = ctrl & (nep >= min_epochs) & np.isfinite(v)
    if not okc.any():
        return float("nan"), np.zeros(v.shape, dtype=bool)
    floor = float(np.median(v[okc]))
    above = (~ctrl) & (nep >= min_epochs) & np.isfinite(v) & (v > floor)
    return floor, above


def synthetic_lightcurves(
    n_sources: int = 400,
    n_epochs: int = 10,
    *,
    frac_variable: float = 0.08,
    var_amp: float = 2.0,
    err_frac: float = 0.07,
    miss_frac: float = 0.25,
    seed: int = 0,
) -> dict:
    """Synthetic dual-band VLBI population: steady sources + an injected variable subset.

    Steady sources have a constant mean flux per band (so $\\eta\\approx1$, $V\\approx$ the measurement
    error); the injected variable fraction gets a single-session flare of relative amplitude ``var_amp``
    (high $\\eta$ and $V$). Each source has a flat-ish S/X index, ``err_frac`` fractional errors, and a
    fraction ``miss_frac`` of sessions randomly missing (``nan``) to mimic uneven VLBI sampling. Returns
    a dict with ``flux_x/err_x/flux_s/err_s`` ``(N, M)`` matrices and the boolean ``is_variable`` truth.
    """
    rng = np.random.default_rng(seed)
    n = n_sources
    mean_x = 10.0 ** rng.uniform(-1.0, 0.5, n)  # ~0.1-3 Jy
    alpha = rng.normal(0.0, 0.2, n)  # flat-spectrum compact AGN
    mean_s = mean_x * (NU_S_GHZ / NU_X_GHZ) ** alpha

    is_variable = rng.random(n) < frac_variable
    flare_epoch = rng.integers(0, n_epochs, n)  # shared across bands: a real flare is broadband

    def _band(mean: np.ndarray, rs: np.random.Generator) -> tuple:
        f = mean[:, None] * (1.0 + rs.normal(0.0, err_frac, (n, n_epochs)))
        # inject a single-epoch flare into the variable subset
        boost = np.zeros((n, n_epochs))
        boost[np.arange(n), flare_epoch] = var_amp
        f = f * (1.0 + np.where(is_variable[:, None], boost, 0.0))
        e = err_frac * mean[:, None] * np.ones((n, n_epochs))
        miss = rs.random((n, n_epochs)) < miss_frac
        # never drop a variable's flare epoch -- an undetected flare is just a steady curve, not a
        # measurement of the injected truth, so the fixture keeps the injected signal observable
        miss[np.arange(n), flare_epoch] = np.where(
            is_variable, False, miss[np.arange(n), flare_epoch]
        )
        f = np.where(miss, np.nan, f)
        e = np.where(miss, np.nan, e)
        return f, e

    fx, ex = _band(mean_x, rng)
    fs, es = _band(mean_s, rng)
    return {
        "flux_x": fx,
        "err_x": ex,
        "flux_s": fs,
        "err_s": es,
        "is_variable": is_variable,
    }


ASTROGEO_BASE = "http://astrogeo.org/images"
# Geodetic/absolute-astrometry VLBI has no per-observation flux error; the dominant uncertainty is
# amplitude calibration. We adopt 5% as a common VLBI starting assumption (it is NOT a value prescribed
# by a specific Astrogeo paper) -- and the CSO control-floor analysis then shows the *effective*
# per-session scatter is several times larger. This is THE assumption the variability rests on, so it
# is a documented, tunable parameter and absolute eta/chi^2 is not trusted as a discriminant.
VLBI_CAL_FRAC = 0.05


def _parse_cfd_tab(text: str) -> tuple[float, float] | None:  # pragma: no cover - network
    """Pull (Fl_int, Fl_noi) in Jy from a one-row Astrogeo ``_cfd.tab`` correlated-flux file."""
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c = line.split()
        if len(c) >= 8:
            return float(c[3]), float(c[7])  # Fl_int (total correlated flux), Fl_noi (image noise)
    return None


def fetch_astrogeo(
    sources: list[str],
    *,
    bands: tuple[str, ...] = ("S", "X"),
    cal_frac: float = VLBI_CAL_FRAC,
    pause: float = 0.15,
) -> dict:  # pragma: no cover - network
    """Per-source, per-epoch VLBI flux histories from Astrogeo (Petrov), keyed by band.

    For each J2000 source name (e.g. ``"J2202+4216"``) we read the source's image directory listing,
    pick out the per-epoch ``_cfd.tab`` correlated-flux files for each requested band, and read each
    one's integrated flux density ``Fl_int`` (Jy). The per-point error is
    ``sqrt((cal_frac*Fl_int)^2 + Fl_noi^2)`` --- a calibration-fraction floor (see :data:`VLBI_CAL_FRAC`)
    in quadrature with the image noise. Returns ``{band: (flux, err)}`` with each an aligned
    ``(n_sources, n_epochs)`` matrix padded with ``nan``. Network-gated; tests use the synthetic fixture.
    """
    import re
    import time

    import requests

    sess = requests.Session()
    sess.headers["User-Agent"] = "jansky-research (amateur radio-astronomy research)"
    per: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {b: [] for b in bands}
    max_ep: dict[str, int] = {b: 0 for b in bands}
    for name in sources:
        try:
            idx = sess.get(f"{ASTROGEO_BASE}/{name}/", timeout=60).text
        except Exception:
            idx = ""
        files = sorted(set(re.findall(rf"{re.escape(name)}_[A-Z]_[0-9_]+[a-z]+_cfd\.tab", idx)))
        for b in bands:
            flux: list[float] = []
            err: list[float] = []
            for fn in (f for f in files if f"_{b}_" in f):
                try:
                    parsed = _parse_cfd_tab(
                        sess.get(f"{ASTROGEO_BASE}/{name}/{fn}", timeout=60).text
                    )
                except Exception:
                    parsed = None
                if parsed is not None:
                    fl, noi = parsed
                    flux.append(fl)
                    err.append(float(np.hypot(cal_frac * fl, noi)))
                if pause:
                    time.sleep(pause)
            per[b].append((np.asarray(flux, float), np.asarray(err, float)))
            max_ep[b] = max(max_ep[b], len(flux))
    n = len(sources)
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for b in bands:
        fmat = np.full((n, max(max_ep[b], 1)), np.nan)
        emat = np.full((n, max(max_ep[b], 1)), np.nan)
        for i, (farr, earr) in enumerate(per[b]):
            fmat[i, : farr.size] = farr
            emat[i, : earr.size] = earr
        out[b] = (fmat, emat)
    return out


def run(
    out: str = ".",
    *,
    offline: bool = True,
    sources: list[str] | None = None,
    controls: list[str] | None = None,
) -> dict:
    """Full slice: variability-rank a (synthetic or fetched) VLBI population and write outputs.

    ``controls`` names a subset of ``sources`` known to be intrinsically steady (e.g. CSOs); their
    median $V$ sets the empirical variability floor (:func:`variability_floor`) above which non-control
    sources are the amplitude-selected variables.
    """
    import json
    from pathlib import Path

    names: list[str] | None
    if offline or sources is None:
        pop = synthetic_lightcurves()
        source = "synthetic"
        truth: np.ndarray | None = pop["is_variable"]
        names = None
    else:  # pragma: no cover - network
        data = fetch_astrogeo(sources)
        fx, ex = data["X"]
        fs, es = data["S"]
        pop = {"flux_x": fx, "err_x": ex, "flux_s": fs, "err_s": es}
        source = f"Astrogeo VLBI ({len(sources)} sources)"
        truth = None
        names = list(sources)

    eta, v, pval, nep, mean = lightcurve_metrics(pop["flux_x"], pop["err_x"])
    alpha, _aerr = sx_index(pop["flux_s"], pop["flux_x"], pop["err_s"], pop["err_x"])
    mask, eta_thr, v_thr = select_variable(eta, v, nep)

    testable = (nep >= MIN_EPOCHS) & np.isfinite(eta)
    n_testable = int(testable.sum())
    # control-floor analysis (the meaningful selector for a calibrator-dominated set; the relative
    # outlier cut above is for blind fields and returns ~nothing when most sources vary)
    ctrl_set = set(controls or [])
    is_control = (
        np.array([nm in ctrl_set for nm in names], dtype=bool)
        if names
        else np.zeros(eta.shape, dtype=bool)
    )
    v_floor, above = variability_floor(v, nep, is_control)
    # the recover-a-known anchor: the single most significant source by eta among testable ones
    top = int(np.argmax(np.where(testable, eta, -np.inf))) if testable.any() else -1
    metrics: dict = {
        "source": source,
        "n_sources": int(pop["flux_x"].shape[0]),
        "n_testable": n_testable,
        "n_candidates": int(mask.sum()),
        "eta_thr": round(float(eta_thr), 3) if np.isfinite(eta_thr) else None,
        "v_thr": round(float(v_thr), 3) if np.isfinite(v_thr) else None,
        "median_alpha_sx": round(float(np.nanmedian(alpha)), 3)
        if np.isfinite(alpha).any()
        else None,
    }
    if np.isfinite(v_floor):  # pragma: no cover - only with named controls (network run)
        metrics["n_controls"] = int(is_control.sum())
        metrics["n_noncontrol"] = int((testable & ~is_control).sum())
        metrics["v_floor"] = round(v_floor, 3)
        metrics["n_above_floor"] = int(above.sum())
        metrics["median_v_control"] = round(float(np.median(v[is_control & testable])), 3)
        var_v = v[above]
        metrics["median_v_variable"] = round(float(np.median(var_v)), 3) if var_v.size else None
    if top >= 0:
        metrics["top_variable"] = {
            "name": names[top] if names else f"row{top}",
            "n_epochs": int(nep[top]),
            "eta": round(float(eta[top]), 1),
            "v": round(float(v[top]), 3),
            "p_value": float(f"{pval[top]:.2e}"),
            "mean_flux_jy": round(float(mean[top]), 3),
            "alpha_sx": round(float(alpha[top]), 2) if np.isfinite(alpha[top]) else None,
        }
    if truth is not None:
        tp = int((mask & truth).sum())
        completeness = tp / int(truth.sum()) if truth.sum() else float("nan")
        purity = tp / int(mask.sum()) if mask.sum() else float("nan")
        metrics["n_injected_variable"] = int(truth.sum())
        metrics["completeness"] = round(completeness, 3)
        metrics["purity"] = round(purity, 3)

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "vlbi_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    if names is not None:  # pragma: no cover - network
        _write_candidates(
            op / "results" / "vlbi_candidates.csv",
            names,
            eta,
            v,
            pval,
            nep,
            mean,
            alpha,
            above,
            is_control,
        )
    _figure(eta, v, above, is_control, v_floor, op / "papers" / "vlbi" / "figures")
    _write_macros(metrics, op / "papers" / "vlbi" / "generated" / "macros.tex")
    return metrics


def _write_candidates(
    path, names, eta, v, pval, nep, mean, alpha, above, is_control
) -> None:  # pragma: no cover
    """Write the full variability-ranked table (most significant first), tagging controls/variables."""
    import csv
    from pathlib import Path

    rows = sorted(
        range(len(names)), key=lambda i: eta[i] if np.isfinite(eta[i]) else -np.inf, reverse=True
    )
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "name",
                "control",
                "above_floor",
                "n_epochs",
                "eta",
                "v",
                "p_value",
                "mean_flux_jy",
                "alpha_sx",
            ]
        )
        for i in rows:
            w.writerow(
                [
                    names[i],
                    int(bool(is_control[i])),
                    int(bool(above[i])),
                    int(nep[i]),
                    f"{eta[i]:.2f}" if np.isfinite(eta[i]) else "",
                    f"{v[i]:.3f}" if np.isfinite(v[i]) else "",
                    f"{pval[i]:.2e}" if np.isfinite(pval[i]) else "",
                    f"{mean[i]:.3f}" if np.isfinite(mean[i]) else "",
                    f"{alpha[i]:.2f}" if np.isfinite(alpha[i]) else "",
                ]
            )


def _figure(eta, v, above, is_control, v_floor, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ok = np.isfinite(eta) & np.isfinite(v) & (eta > 0) & (v > 0)
    ctrl = np.asarray(is_control, bool)
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    rest = ok & ~ctrl & ~above
    ax.scatter(eta[rest], v[rest], s=10, c="0.6", label="below floor")
    ax.scatter(eta[ok & above], v[ok & above], s=26, c="C3", label="variable (above floor)")
    if ctrl.any():  # pragma: no cover - only with named controls (network run)
        ax.scatter(
            eta[ok & ctrl],
            v[ok & ctrl],
            s=46,
            marker="s",
            facecolors="none",
            edgecolors="C0",
            label="steady control (CSO)",
        )
    if np.isfinite(v_floor):  # pragma: no cover - only with named controls (network run)
        ax.axhline(v_floor, ls="--", c="C0", lw=0.9, label=f"floor V={v_floor:.2f}")
    ax.set(
        xscale="log",
        xlabel=r"$\eta$ (significance vs.\ constant)",
        ylabel=r"$V$ (fractional amplitude)",
        title="VLBI variability (X band)",
    )
    ax.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "etav.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    # Always emit the FULL union of macros (placeholder "--" for whichever set is inactive) so the
    # paper compiles identically from offline-synthetic macros (CI) and real macros (make reproduce).
    tv = m.get("top_variable") or {}

    def _tv(key: str) -> str:
        val = tv.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.vlbi._write_macros — do not edit by hand.",
        rf"\newcommand{{\viSource}}{{{m['source']}}}",
        rf"\newcommand{{\viN}}{{{m['n_sources']}}}",
        rf"\newcommand{{\viTestable}}{{{m['n_testable']}}}",
        rf"\newcommand{{\viNcand}}{{{m['n_candidates']}}}",
        rf"\newcommand{{\viEtaThr}}{{{_fmt('eta_thr')}}}",
        rf"\newcommand{{\viVThr}}{{{_fmt('v_thr')}}}",
        rf"\newcommand{{\viMedAlpha}}{{{_fmt('median_alpha_sx')}}}",
        rf"\newcommand{{\viNctrl}}{{{_fmt('n_controls')}}}",
        rf"\newcommand{{\viNnoncontrol}}{{{_fmt('n_noncontrol')}}}",
        rf"\newcommand{{\viFloor}}{{{_fmt('v_floor')}}}",
        rf"\newcommand{{\viNabove}}{{{_fmt('n_above_floor')}}}",
        rf"\newcommand{{\viMedVctrl}}{{{_fmt('median_v_control')}}}",
        rf"\newcommand{{\viMedVvar}}{{{_fmt('median_v_variable')}}}",
        rf"\newcommand{{\viInjected}}{{{_fmt('n_injected_variable')}}}",
        rf"\newcommand{{\viCompleteness}}{{{_fmt('completeness')}}}",
        rf"\newcommand{{\viPurity}}{{{_fmt('purity')}}}",
        rf"\newcommand{{\viTopName}}{{{_tv('name')}}}",
        rf"\newcommand{{\viTopEpochs}}{{{_tv('n_epochs')}}}",
        rf"\newcommand{{\viTopEta}}{{{_tv('eta')}}}",
        rf"\newcommand{{\viTopV}}{{{_tv('v')}}}",
        rf"\newcommand{{\viTopFlux}}{{{_tv('mean_flux_jy')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Multi-decade VLBI flux variability (Astrogeo).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument(
        "--online", action="store_true", help="run on the curated VALIDATION_SOURCES set"
    )
    p.add_argument("--sources", nargs="*", help="explicit Astrogeo J2000 source names")
    args = p.parse_args(argv)
    sources = args.sources or (list(VALIDATION_SOURCES) if args.online else None)
    # the steady controls are the CSOs in the curated set
    controls = (
        [j for j, name in VALIDATION_SOURCES.items() if "CSO" in name] if args.online else None
    )
    metrics = run(args.out, offline=args.offline or not sources, sources=sources, controls=controls)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
