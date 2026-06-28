"""VLASS multi-epoch radio variability — find variable/transient sources across the VLA Sky Survey.

The VLA Sky Survey (VLASS; Lacy et al. 2020) imaged the sky north of Dec $-40\\degree$ at 2--4 GHz,
2.5" resolution, in **three epochs** (2017--2024). Cross-matching the per-epoch CIRADA Quick-Look
component catalogues (Gordon et al. 2021) by position gives each source a 2--3-point radio light
curve. Published variability work stops at Epochs 1--2; adding Epoch 3 gives a three-point curve.

Two metrics, following the transient-survey convention (de Vries et al. 2004; Scheers 2011;
Swinbank et al. 2015; Rowlinson et al. 2019):

- $V = \\sigma_S / \\bar S$ -- the coefficient of variation (amplitude of variability);
- $\\eta = \\frac{1}{N-1}\\sum_i (S_i - \\bar S_w)^2 / \\sigma_i^2$ -- the weighted reduced chi-square
  against a constant flux (significance of variability), with $\\bar S_w$ the inverse-variance
  weighted mean.

Sources that are outliers in **both** $\\log\\eta$ and $\\log V$ are variable/transient candidates.
Pure NumPy/SciPy + an astropy positional cross-match; a synthetic three-epoch population with an
injected variable subset lets the tests run offline.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "VariabilityMetrics",
    "crossmatch_epochs",
    "debiased_modulation_index",
    "eta_metric",
    "run",
    "select_candidates",
    "synthetic_epochs",
    "v_metric",
    "variability_metrics",
]

VLASS_MATCH_ARCSEC = 2.5  # one VLASS resolution element; the default cross-epoch match radius


def _weighted_mean(flux: np.ndarray, err: np.ndarray) -> float:
    w = 1.0 / np.square(err)
    return float(np.sum(flux * w) / np.sum(w))


def eta_metric(flux: np.ndarray, err: np.ndarray) -> float:
    """Weighted reduced chi-square of a light curve against a constant flux (variability significance).

    $\\eta = \\frac{1}{N-1}\\sum_i (S_i - \\bar S_w)^2/\\sigma_i^2$. A steady source measured within its
    errors gives $\\eta \\approx 1$; a genuine variable gives $\\eta \\gg 1$.
    """
    flux = np.asarray(flux, dtype=float)
    err = np.asarray(err, dtype=float)
    n = flux.size
    if n < 2:
        return 0.0
    sw = _weighted_mean(flux, err)
    return float(np.sum(np.square(flux - sw) / np.square(err)) / (n - 1))


def v_metric(flux: np.ndarray) -> float:
    """Coefficient of variation $V = \\sigma_S/\\bar S$ (amplitude of variability)."""
    flux = np.asarray(flux, dtype=float)
    if flux.size < 2:
        return 0.0
    mean = float(np.mean(flux))
    if mean == 0.0:
        return 0.0
    return float(np.std(flux, ddof=1) / mean)


def debiased_modulation_index(flux: np.ndarray, err: np.ndarray) -> float:
    """Modulation index with the measurement noise removed: $m_d = \\sqrt{\\sigma_S^2 - \\langle\\sigma_i^2\\rangle}/\\bar S$.

    Clipped at zero. Positive $m_d$ means variability in excess of the measurement errors -- the
    honest amplitude, since the raw coefficient of variation is inflated by noise on faint sources.
    """
    flux = np.asarray(flux, dtype=float)
    err = np.asarray(err, dtype=float)
    if flux.size < 2:
        return 0.0
    mean = float(np.mean(flux))
    if mean == 0.0:
        return 0.0
    excess = np.var(flux, ddof=1) - float(np.mean(np.square(err)))
    return float(np.sqrt(max(excess, 0.0)) / mean)


@dataclass(frozen=True)
class VariabilityMetrics:
    """Per-source variability summary."""

    n_epochs: int
    mean_flux: float
    eta: float
    v: float
    m_debiased: float
    chi2: float
    p_value: float  # probability the light curve is consistent with a constant flux


def variability_metrics(flux: np.ndarray, err: np.ndarray) -> VariabilityMetrics:
    """Full variability summary for one source's light curve (flux + error per epoch)."""
    from scipy import stats

    flux = np.asarray(flux, dtype=float)
    err = np.asarray(err, dtype=float)
    n = flux.size
    eta = eta_metric(flux, err)
    chi2 = eta * (n - 1) if n >= 2 else 0.0
    p = float(stats.chi2.sf(chi2, df=n - 1)) if n >= 2 else 1.0
    return VariabilityMetrics(
        n_epochs=n,
        mean_flux=float(np.mean(flux)),
        eta=eta,
        v=v_metric(flux),
        m_debiased=debiased_modulation_index(flux, err),
        chi2=chi2,
        p_value=p,
    )


def crossmatch_epochs(
    ra: list[np.ndarray],
    dec: list[np.ndarray],
    flux: list[np.ndarray],
    err: list[np.ndarray],
    *,
    radius_arcsec: float = VLASS_MATCH_ARCSEC,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Cross-match per-epoch catalogues positionally onto the first epoch.

    Each argument is a list with one array per epoch. The first epoch is the anchor; later epochs
    are matched to it within ``radius_arcsec``. Returns ``(ra0, dec0, flux_matrix, err_matrix)`` where
    the matrices have shape ``(n_anchor, n_epochs)`` and unmatched cells are ``nan``.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    n_epochs = len(ra)
    ra0 = np.asarray(ra[0], dtype=float)
    dec0 = np.asarray(dec[0], dtype=float)
    n = ra0.size
    fmat = np.full((n, n_epochs), np.nan)
    emat = np.full((n, n_epochs), np.nan)
    fmat[:, 0] = np.asarray(flux[0], dtype=float)
    emat[:, 0] = np.asarray(err[0], dtype=float)
    anchor = SkyCoord(ra0 * u.deg, dec0 * u.deg)
    for e in range(1, n_epochs):
        cat = SkyCoord(np.asarray(ra[e]) * u.deg, np.asarray(dec[e]) * u.deg)
        idx, sep, _ = anchor.match_to_catalog_sky(cat)
        ok = sep.arcsec <= radius_arcsec
        fmat[ok, e] = np.asarray(flux[e], dtype=float)[idx[ok]]
        emat[ok, e] = np.asarray(err[e], dtype=float)[idx[ok]]
    return ra0, dec0, fmat, emat


def _clipped_threshold(x: np.ndarray, sigma: float, *, n_iter: int = 3) -> float:
    """Sigma-clipped ``mean + sigma*std`` of ``x`` (robust to the variable-source tail)."""
    x = x[np.isfinite(x)]
    keep = np.ones(x.size, dtype=bool)
    mean = std = 0.0
    for _ in range(n_iter):
        mean = float(np.mean(x[keep]))
        std = float(np.std(x[keep]))
        if std == 0.0:
            break
        keep = np.abs(x - mean) < sigma * std
    return mean + sigma * std


def select_candidates(
    eta: np.ndarray, v: np.ndarray, *, sigma: float = 3.0
) -> tuple[np.ndarray, float, float]:
    """Flag variable candidates as outliers in **both** $\\log\\eta$ and $\\log V$.

    Fits a sigma-clipped Gaussian to $\\log_{10}\\eta$ and $\\log_{10}V$ over the (mostly steady)
    population and flags sources above ``sigma`` in both (Rowlinson et al. 2019). Returns
    ``(mask, eta_threshold, v_threshold)`` with the thresholds in linear units.
    """
    eta = np.asarray(eta, dtype=float)
    v = np.asarray(v, dtype=float)
    good = np.isfinite(eta) & np.isfinite(v) & (eta > 0) & (v > 0)
    le = np.log10(np.where(good, eta, np.nan))
    lv = np.log10(np.where(good, v, np.nan))
    eta_thr = _clipped_threshold(le[good], sigma)
    v_thr = _clipped_threshold(lv[good], sigma)
    mask = good & (le > eta_thr) & (lv > v_thr)
    return mask, float(10.0**eta_thr), float(10.0**v_thr)


def synthetic_epochs(
    n_sources: int = 2000,
    *,
    n_epochs: int = 3,
    var_fraction: float = 0.05,
    rel_err: float = 0.1,
    amp_dex: float = 0.8,
    seed: int | None = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Synthetic multi-epoch catalogue with a known injected variable subset (offline fixture).

    Steady sources have a constant intrinsic flux measured with fractional error ``rel_err`` (so
    $\\eta\\approx1$, $V\\approx$ ``rel_err``); the injected variable fraction takes a strong,
    transient-like multiplicative excursion per epoch (``amp_dex`` in dex), giving $\\eta\\gg1$ and
    $V\\gg$ ``rel_err``. With only three epochs the steady $\\eta$/$V$ estimators are themselves
    noisy, so a conservative $3\\sigma$ selection recovers the clearly variable sources, not the
    marginal ones -- the honest behaviour. Returns ``(ra, dec, flux[n,e], err[n,e], is_variable)``;
    fluxes already include measurement noise.
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(150.0, 160.0, n_sources)
    dec = rng.uniform(20.0, 30.0, n_sources)
    base = 10.0 ** rng.uniform(0.0, 1.5, n_sources)  # ~1-30 mJy, log-uniform
    is_var = rng.random(n_sources) < var_fraction
    truth = np.repeat(base[:, None], n_epochs, axis=1)
    swing = 10.0 ** rng.normal(
        0.0, amp_dex, size=(n_sources, n_epochs)
    )  # per-epoch excursion (dex)
    truth = np.where(is_var[:, None], base[:, None] * swing, truth)
    err = rel_err * truth
    flux = truth + rng.normal(0.0, err)
    flux = np.clip(flux, 1e-3, None)
    return ra, dec, flux, err, is_var


def run(out: str = ".", *, offline: bool = False, sigma: float = 3.0) -> dict:
    """Build the multi-epoch variability catalogue (synthetic offline, or real CIRADA). Writes a figure."""
    import json
    from pathlib import Path

    if offline:
        ra, dec, flux, err, truth = synthetic_epochs()
        source = "synthetic"
    else:  # pragma: no cover - network
        ra, dec, flux, err = _fetch_and_match()
        truth = None
        source = "VLASS CIRADA Quick-Look (Epochs 1-3)"

    # per-source metrics across the available epochs
    detected = np.sum(np.isfinite(flux), axis=1) >= 2
    eta = np.array(
        [eta_metric(f[np.isfinite(f)], e[np.isfinite(f)]) for f, e in zip(flux, err, strict=True)]
    )
    v = np.array([v_metric(f[np.isfinite(f)]) for f in flux])
    eta[~detected] = np.nan
    v[~detected] = np.nan
    mask, eta_thr, v_thr = select_candidates(eta, v, sigma=sigma)

    metrics = {
        "source": source,
        "n_sources": int(ra.size),
        "n_detected_multi_epoch": int(detected.sum()),
        "n_candidates": int(mask.sum()),
        "eta_threshold": eta_thr,
        "v_threshold": v_thr,
        "sigma": sigma,
    }
    if truth is not None:  # synthetic: report recovery of the injected variables
        n_true = int(truth.sum())
        metrics["n_injected_variables"] = n_true
        metrics["recovered_fraction"] = float(np.sum(mask & truth) / n_true) if n_true else 0.0
        metrics["false_positive_fraction"] = (
            float(np.sum(mask & ~truth) / max(int(mask.sum()), 1)) if mask.any() else 0.0
        )

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "vlass_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(eta, v, mask, eta_thr, v_thr, op / "papers" / "vlass" / "figures")
    return metrics


def _fetch_and_match():  # pragma: no cover - network
    """Fetch the three CIRADA VLASS epoch catalogues for a region and cross-match them.

    Placeholder for the real-data path: a CADC/CIRADA TAP cone query per epoch
    (``https://cirada.ca`` / CADC ``ivo://cadc.nrc.ca/vlass``), returning per-epoch
    ``(ra, dec, peak_flux, peak_flux_err)`` arrays passed through :func:`crossmatch_epochs`.
    """
    raise NotImplementedError(
        "Real CIRADA fetch is not wired yet; use run(offline=True) for the synthetic fixture."
    )


def _figure(eta, v, mask, eta_thr, v_thr, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    good = np.isfinite(eta) & np.isfinite(v) & (eta > 0) & (v > 0)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.loglog(eta[good & ~mask], v[good & ~mask], ".", ms=2, color="0.6", label="steady")
    ax.loglog(eta[mask], v[mask], "*", ms=7, color="r", label="candidates")
    ax.axvline(eta_thr, color="r", ls=":", lw=0.8)
    ax.axhline(v_thr, color="r", ls=":", lw=0.8)
    ax.set(
        xlabel=r"$\eta$ (variability significance)",
        ylabel=r"$V$ (variability amplitude)",
        title="VLASS multi-epoch variability",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "eta_v.pdf")
    plt.close(fig)


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Build the VLASS multi-epoch variability catalogue.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true", help="use the synthetic fixture (no network)")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
