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
]

NU_S_GHZ = 2.3  # Astrogeo S band
NU_X_GHZ = 8.4  # Astrogeo X band
MIN_EPOCHS = 4  # a light curve needs at least this many finite epochs to be tested


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


def fetch_astrogeo(sources: list[str], *, band: str = "X") -> dict:  # pragma: no cover - network
    """Per-source flux-density histories from the Astrogeo VLBI database (S or X band).

    For each source name/coordinate, retrieves the per-session total flux densities from Astrogeo
    (astrogeo.org) and returns aligned ``(n_sources, n_epochs)`` flux/error matrices padded with
    ``nan``. Network-gated; the offline tests use :func:`synthetic_lightcurves` instead.
    """
    import requests

    base = "https://astrogeo.smce.nasa.gov/vlbi_images"
    histories: list[tuple[np.ndarray, np.ndarray]] = []
    max_ep = 0
    for name in sources:
        try:
            url = f"{base}/{name}/{name}_{band.lower()}_flux.txt"
            txt = requests.get(url, timeout=120).text
            rows = [r.split() for r in txt.splitlines() if r and not r.startswith("#")]
            flux = np.array([float(r[1]) for r in rows], float)
            err = np.array([float(r[2]) for r in rows], float)
        except Exception:
            flux, err = np.array([]), np.array([])
        histories.append((flux, err))
        max_ep = max(max_ep, flux.size)
    n = len(sources)
    fmat = np.full((n, max_ep), np.nan)
    emat = np.full((n, max_ep), np.nan)
    for i, (flux, err) in enumerate(histories):
        fmat[i, : flux.size] = flux
        emat[i, : err.size] = err
    return {"flux": fmat, "err": emat}


def run(out: str = ".", *, offline: bool = True, sources: list[str] | None = None) -> dict:
    """Full slice: variability-rank a (synthetic or fetched) VLBI population and write outputs."""
    import json
    from pathlib import Path

    if offline or sources is None:
        pop = synthetic_lightcurves()
        source = "synthetic"
        truth: np.ndarray | None = pop["is_variable"]
    else:  # pragma: no cover - network
        fx = fetch_astrogeo(sources, band="X")
        fs = fetch_astrogeo(sources, band="S")
        pop = {"flux_x": fx["flux"], "err_x": fx["err"], "flux_s": fs["flux"], "err_s": fs["err"]}
        source = f"Astrogeo VLBI ({len(sources)} sources)"
        truth = None

    eta, v, _pval, nep, _mean = lightcurve_metrics(pop["flux_x"], pop["err_x"])
    alpha, _aerr = sx_index(pop["flux_s"], pop["flux_x"], pop["err_s"], pop["err_x"])
    mask, eta_thr, v_thr = select_variable(eta, v, nep)

    n_testable = int(((nep >= MIN_EPOCHS) & np.isfinite(eta)).sum())
    cand_alpha = alpha[mask]
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
        "median_cand_alpha_sx": (
            round(float(np.nanmedian(cand_alpha)), 3) if np.isfinite(cand_alpha).any() else None
        ),
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
    _figure(eta, v, mask, eta_thr, v_thr, op / "papers" / "vlbi" / "figures")
    _write_macros(metrics, op / "papers" / "vlbi" / "generated" / "macros.tex")
    return metrics


def _figure(eta, v, mask, eta_thr, v_thr, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ok = np.isfinite(eta) & np.isfinite(v) & (eta > 0) & (v > 0)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.scatter(eta[ok & ~mask], v[ok & ~mask], s=8, c="0.6", label="steady")
    ax.scatter(eta[ok & mask], v[ok & mask], s=22, c="C3", label="candidate")
    if np.isfinite(eta_thr):
        ax.axvline(eta_thr, ls="--", c="0.4", lw=0.8)
    if np.isfinite(v_thr):
        ax.axhline(v_thr, ls="--", c="0.4", lw=0.8)
    ax.set(
        xscale="log",
        yscale="log",
        xlabel=r"$\eta$ (significance)",
        ylabel=r"$V$ (amplitude)",
        title="VLBI variability (X band)",
    )
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "etav.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
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
        rf"\newcommand{{\viCandAlpha}}{{{_fmt('median_cand_alpha_sx')}}}",
    ]
    if "completeness" in m:
        lines += [
            rf"\newcommand{{\viInjected}}{{{m['n_injected_variable']}}}",
            rf"\newcommand{{\viCompleteness}}{{{m['completeness']}}}",
            rf"\newcommand{{\viPurity}}{{{m['purity']}}}",
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
    p.add_argument("--sources", nargs="*", help="Astrogeo source names for a real run")
    args = p.parse_args(argv)
    metrics = run(args.out, offline=args.offline or not args.sources, sources=args.sources)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
