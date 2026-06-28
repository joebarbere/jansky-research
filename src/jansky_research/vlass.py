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
    "apply_flux_scale",
    "confirm_candidates",
    "crossmatch_epochs",
    "debiased_modulation_index",
    "eta_metric",
    "fetch_vlass_cutout",
    "fetch_vlass_epoch",
    "forced_photometry",
    "image_lightcurve",
    "injection_recovery",
    "isolated_mask",
    "measure_image_flux",
    "run",
    "select_candidates",
    "synthetic_epochs",
    "v_metric",
    "variability_metrics",
    "vet_candidates",
]

VLASS_MATCH_ARCSEC = 2.5  # one VLASS resolution element; the default cross-epoch match radius

# Real-data access. There is no single TAP for all epochs: Epoch 1 is on VizieR TAP
# (Gordon et al. 2021, ApJS 255, 30); Epochs 2-3 are bulk files on the NRAO server with a
# different schema. See survey/data-source-scan.md and VLASS Memos 13 & 22.
VLASS_TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"
VLASS_E1_TABLE = '"J/ApJS/255/30/comp"'
VLASS_BULK_URLS: dict[int, list[str]] = {
    2: [
        "https://vlass-dl.nrao.edu/vlass/quicklook/catalogs/epoch2/CIRADA_VLASS2QLv2_table1_components.csv.gz"
    ],
    3: [
        "https://vlass-dl.nrao.edu/vlass/quicklook/catalogs/epoch3/QL3.1_components.fits",
        "https://vlass-dl.nrao.edu/vlass/quicklook/catalogs/epoch3/QL3.2_components.fits",
    ],
}
# Multiplicative peak-flux corrections onto the Perley-Butler 2017 scale (VLASS Memos 13/22):
# the Quick-Look images *underestimate* peak flux, epoch-dependently. Epoch 1 is the mean of the
# 1.1 (~15% low) and 1.2 (~8% low) campaigns. WITHOUT this, the epoch-to-epoch offset alone makes
# every source look variable -- the single most important systematic for VLASS variability.
VLASS_PEAK_CORRECTION: dict[int, float] = {1: 1.13, 2: 1.075, 3: 1.031}
# Residual epoch-to-epoch flux-scale scatter for bright point sources after correction (~7%, Memo
# 22). Added in quadrature to the catalogue errors so it does not masquerade as variability.
VLASS_SYS_FRAC = 0.07


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


def apply_flux_scale(
    epoch: int, flux: np.ndarray, err: np.ndarray, *, sys_frac: float = VLASS_SYS_FRAC
) -> tuple[np.ndarray, np.ndarray]:
    """Put one epoch's Quick-Look peak fluxes on the common Perley-Butler scale and add the systematic floor.

    Multiplies ``flux`` and ``err`` by the per-epoch correction (:data:`VLASS_PEAK_CORRECTION`), then
    adds the residual cross-epoch scale scatter ``sys_frac``$\\times S$ in quadrature to the error.
    This is the GATE-2-critical step: without the correction the epochs sit on different flux scales
    and every source looks variable; without the systematic floor $\\eta$ is inflated and steady
    sources cross the threshold.
    """
    flux = np.asarray(flux, dtype=float)
    err = np.asarray(err, dtype=float)
    corr = VLASS_PEAK_CORRECTION.get(epoch, 1.0)
    cflux = flux * corr
    cerr = np.hypot(err * corr, sys_frac * cflux)
    return cflux, cerr


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


def isolated_mask(ra: np.ndarray, dec: np.ndarray, *, radius_arcsec: float = 5.0) -> np.ndarray:
    """True for sources with no catalogue neighbour within ``radius_arcsec`` (deblending-safe).

    Multi-epoch "variability" in component catalogues is frequently spurious: a single slightly
    extended source is deblended into a different number of Gaussian components in each epoch, so a
    secondary component's flux jumps around. Image vetting of the first real candidate showed exactly
    this. Requiring isolation drops the crowded sources whose per-epoch deblending is unreliable.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    ra = np.asarray(ra, dtype=float)
    if ra.size < 2:
        return np.ones(ra.size, dtype=bool)
    sc = SkyCoord(ra * u.deg, np.asarray(dec, float) * u.deg)
    _, sep, _ = sc.match_to_catalog_sky(sc, nthneighbor=2)  # distance to the nearest *other* source
    return np.asarray(sep.arcsec > radius_arcsec)


def injection_recovery(
    flux: np.ndarray,
    err: np.ndarray,
    *,
    factors: tuple[float, ...] = (1.25, 1.5, 2.0, 3.0, 5.0, 10.0),
    sigma: float = 3.0,
    n_per_factor: int = 400,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Selection completeness vs single-epoch flare amplitude, by injection into the real light curves.

    For each flare ``factor``, takes steady (all-epoch-detected) light curves, multiplies one random
    epoch by the factor (scaling that epoch's error with it, preserving the fractional error), and
    measures the fraction that cross the field's own 2-D $(\\eta, V)$ threshold. This is the honest,
    data-driven completeness: it uses the real per-source noise and the real selection cut, so the
    observed variable fraction can be corrected for what the selection actually sees. Returns
    ``(factors, recovered_fraction)``.
    """
    flux = np.asarray(flux, dtype=float)
    err = np.asarray(err, dtype=float)
    steady = np.all(np.isfinite(flux), axis=1) & np.all(np.isfinite(err) & (err > 0), axis=1)
    fs, es = flux[steady], err[steady]
    factors = tuple(float(f) for f in factors)
    if fs.shape[0] < 20:
        return np.asarray(factors), np.zeros(len(factors))
    eta = np.array([eta_metric(f, e) for f, e in zip(fs, es, strict=True)])
    v = np.array([v_metric(f) for f in fs])
    _, eta_thr, v_thr = select_candidates(eta, v, sigma=sigma)
    rng = np.random.default_rng(seed)
    recovered = []
    for fac in factors:
        hits = 0
        for _ in range(n_per_factor):
            j = int(rng.integers(fs.shape[0]))
            f = fs[j].copy()
            e = es[j].copy()
            k = int(rng.integers(f.size))
            f[k] *= fac
            e[k] *= fac  # error scales with flux -> fractional error preserved
            if eta_metric(f, e) > eta_thr and v_metric(f) > v_thr:
                hits += 1
        recovered.append(hits / n_per_factor)
    return np.asarray(factors), np.asarray(recovered)


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


def run(
    out: str = ".",
    *,
    offline: bool = False,
    center: tuple[float, float] | None = None,
    radius_deg: float = 1.0,
    epochs: tuple[int, ...] = (1, 2, 3),
    sigma: float = 3.0,
    isolation_arcsec: float = 5.0,
) -> dict:
    """Build the multi-epoch variability catalogue (synthetic offline, or real VLASS). Writes a figure.

    Real-data run: pass ``center=(ra_deg, dec_deg)`` and a ``radius_deg`` cone; the requested VLASS
    ``epochs`` are fetched (Epoch 1 via VizieR TAP, Epochs 2-3 via the NRAO bulk catalogues), put on a
    common flux scale (:func:`apply_flux_scale`), cross-matched, and run through the variability
    selection. Needs the ``vlass`` extra (``uv sync --extra vlass``) for ``pyvo``.
    """
    import json
    from pathlib import Path

    if offline:
        ra, dec, flux, err, truth = synthetic_epochs()
        source = "synthetic"
    else:  # pragma: no cover - network
        if center is None:
            raise ValueError("real-data run needs center=(ra_deg, dec_deg); or pass offline=True")
        ra, dec, flux, err = _fetch_and_match(center, radius_deg, epochs)
        truth = None
        source = (
            f"VLASS Quick-Look epochs {','.join(map(str, epochs))} @ {center} r={radius_deg}deg"
        )

    # per-source metrics across the available epochs
    detected = np.sum(np.isfinite(flux), axis=1) >= 2
    eta = np.array(
        [eta_metric(f[np.isfinite(f)], e[np.isfinite(f)]) for f, e in zip(flux, err, strict=True)]
    )
    v = np.array([v_metric(f[np.isfinite(f)]) for f in flux])
    # Drop deblending-prone crowded sources: their per-epoch component fluxes are unreliable and are
    # the dominant false-positive in QL multi-epoch variability (confirmed by image vetting).
    isolated = isolated_mask(ra, dec, radius_arcsec=isolation_arcsec)
    usable = detected & isolated
    eta[~usable] = np.nan
    v[~usable] = np.nan
    mask, eta_thr, v_thr = select_candidates(eta, v, sigma=sigma)

    # Census statistics: variable fraction + data-driven completeness vs flare amplitude (so the
    # observed fraction can be corrected for what the selection actually sees).
    n_usable = int(usable.sum())
    factors, recovered = injection_recovery(flux, err, sigma=sigma)
    c50 = float(np.interp(0.5, recovered, factors)) if recovered.max() >= 0.5 else float("nan")
    c90 = float(np.interp(0.9, recovered, factors)) if recovered.max() >= 0.9 else float("nan")

    metrics = {
        "source": source,
        "n_sources": int(ra.size),
        "n_detected_multi_epoch": int(detected.sum()),
        "n_excluded_crowded": int(np.sum(detected & ~isolated)),
        "n_usable": n_usable,
        "n_candidates": int(mask.sum()),
        "variable_fraction": float(mask.sum() / n_usable) if n_usable else 0.0,
        "eta_threshold": eta_thr,
        "v_threshold": v_thr,
        "sigma": sigma,
        "completeness_factors": factors.tolist(),
        "completeness_recovered": recovered.tolist(),
        "completeness_50_factor": c50,
        "completeness_90_factor": c90,
        "completeness_max": float(recovered.max()),
    }
    if truth is not None:  # synthetic: report recovery of the injected variables
        n_true = int(truth.sum())
        metrics["n_injected_variables"] = n_true
        metrics["recovered_fraction"] = float(np.sum(mask & truth) / n_true) if n_true else 0.0
        metrics["false_positive_fraction"] = (
            float(np.sum(mask & ~truth) / max(int(mask.sum()), 1)) if mask.any() else 0.0
        )

    # On a real-data run, automatically image-vet each survivor: forced photometry confirms (or
    # rejects) the catalogue variability against the images, and SIMBAD/NED give any counterpart.
    conf = (
        confirm_candidates(ra[mask], dec[mask], epochs=epochs)
        if (not offline and mask.any())
        else None
    )
    vet = vet_candidates(ra[mask], dec[mask]) if (not offline and mask.any()) else None
    if conf is not None:
        metrics["n_image_confirmed"] = int(sum(c["confirmed"] for c in conf))

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "vlass_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    figs = op / "papers" / "vlass" / "figures"
    _figure(eta, v, mask, eta_thr, v_thr, figs)
    _completeness_figure(factors, recovered, c50, figs)
    _write_candidates(
        op / "results" / "vlass_candidates.csv", ra, dec, flux, eta, v, mask, conf, vet
    )
    return metrics


def _completeness_figure(factors, recovered, c50, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.semilogx(factors, recovered, "o-", lw=1, label="recovered")
    ax.axhline(0.5, color="0.6", ls=":")
    if np.isfinite(c50):
        ax.axvline(c50, color="r", ls="--", label=f"50% at {c50:.1f}x")
    ax.set(
        xlabel="injected single-epoch flare factor",
        ylabel="recovered fraction",
        title="VLASS variability selection completeness",
        ylim=(0, 1.02),
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "completeness.pdf")
    plt.close(fig)


def _write_candidates(path, ra, dec, flux, eta, v, mask, conf, vet) -> None:
    """Write the candidates: catalogue light curve, forced-photometry confirmation, and any counterpart."""
    import csv as _csv

    n_epochs = flux.shape[1]
    key = lambda r: (round(r["ra"], 6), round(r["dec"], 6))  # noqa: E731
    cmap = {key(r): r for r in (conf or [])}
    vmap = {key(r): r for r in (vet or [])}
    idx = np.where(mask)[0]
    idx = idx[np.argsort(-eta[idx])]
    with open(path, "w", newline="") as fh:
        w = _csv.writer(fh)
        w.writerow(
            ["ra", "dec", *[f"flux_e{e + 1}" for e in range(n_epochs)], "eta", "v"]
            + [
                "image_confirmed",
                *[f"forced_e{e + 1}" for e in range(n_epochs)],
                "forced_eta",
                "forced_v",
            ]
            + [
                "simbad_name",
                "simbad_type",
                "simbad_sep",
                "ned_name",
                "ned_type",
                "ned_z",
                "ned_sep",
            ]
        )
        for i in idx:
            k = (round(float(ra[i]), 6), round(float(dec[i]), 6))
            c = cmap.get(k, {})
            m = vmap.get(k, {})
            fphot = c.get("forced_flux") or [None] * n_epochs
            w.writerow(
                [f"{ra[i]:.6f}", f"{dec[i]:.6f}"]
                + [f"{flux[i, e]:.3f}" if np.isfinite(flux[i, e]) else "" for e in range(n_epochs)]
                + [f"{eta[i]:.2f}", f"{v[i]:.3f}"]
                + [c.get("confirmed", "")]
                + [("" if fphot[e] is None else fphot[e]) for e in range(n_epochs)]
                + [c.get("forced_eta", ""), c.get("forced_v", "")]
                + [
                    m.get("simbad_name", ""),
                    m.get("simbad_type", ""),
                    m.get("simbad_sep", ""),
                    m.get("ned_name", ""),
                    m.get("ned_type", ""),
                    m.get("ned_z", ""),
                    m.get("ned_sep", ""),
                ]
            )


def _cached_download(url: str):  # pragma: no cover - network
    """Download ``url`` into the dataset cache (idempotent); return the local path."""
    from . import data as _data

    target = _data.data_dir() / url.rsplit("/", 1)[-1]
    if not target.exists():
        _data._download(url, target)
    return target


def _cone_mask(ra, dec, center, radius_deg):
    """Boolean mask of ``(ra, dec)`` within ``radius_deg`` of ``center`` (great-circle)."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(center[0] * u.deg, center[1] * u.deg)
    pts = SkyCoord(np.asarray(ra, float) * u.deg, np.asarray(dec, float) * u.deg)
    return c.separation(pts).deg <= radius_deg


def _fetch_e1_tap(center, radius_deg):  # pragma: no cover - network
    """Epoch 1 peak fluxes via the VizieR TAP cone search, with the standard quality cuts."""
    import pyvo

    ra0, dec0 = center
    query = f"""
        SELECT "RAJ2000", "DEJ2000", "Fpeak", "e_Fpeak"
        FROM {VLASS_E1_TABLE}
        WHERE 1=CONTAINS(POINT('ICRS', "RAJ2000", "DEJ2000"),
                         CIRCLE('ICRS', {ra0}, {dec0}, {radius_deg}))
          AND "DupFlag" < 2 AND "QualFlag" IN (0, 4) AND "SCode" != 'E'
    """
    t = pyvo.dal.TAPService(VLASS_TAP_URL).search(query).to_table()
    # The VizieR table labels Fpeak "mJy/beam" but the values are actually micro-Jy (verified
    # empirically: the catalogue's flux floor is ~600-800, i.e. the ~0.7 mJy VLASS detection
    # threshold). Epochs 2-3 (NRAO bulk files) are genuine mJy, so scale Epoch 1 to mJy to match.
    return (
        np.asarray(t["RAJ2000"], float),
        np.asarray(t["DEJ2000"], float),
        np.asarray(t["Fpeak"], float) * 1e-3,
        np.asarray(t["e_Fpeak"], float) * 1e-3,
    )


def _fetch_e2_csv(center, radius_deg):  # pragma: no cover - network
    """Epoch 2 peak fluxes: stream-filter the bulk CIRADA CSV (no pandas) to the region + quality cuts."""
    import csv
    import gzip

    path = _cached_download(VLASS_BULK_URLS[2][0])
    ra, dec, fp, efp = [], [], [], []
    dlo, dhi = center[1] - radius_deg - 0.1, center[1] + radius_deg + 0.1
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh):
            try:
                d = float(r["DEC"])
                if not (dlo <= d <= dhi):  # cheap Dec pre-filter before the cone test
                    continue
                # flags are stored as floats ("0.0", "2.0") -> parse via float, not int()
                if int(float(r["Duplicate_flag"])) >= 2 or int(float(r["Quality_flag"])) not in (
                    0,
                    4,
                ):
                    continue
                if r["S_Code"].strip() == "E":
                    continue
                ra.append(float(r["RA"]))
                dec.append(d)
                fp.append(float(r["Peak_flux"]))
                efp.append(float(r["E_Peak_flux"]))
            except (KeyError, ValueError):
                continue
    ra, dec, fp, efp = (np.asarray(x, float) for x in (ra, dec, fp, efp))
    m = _cone_mask(ra, dec, center, radius_deg)
    return ra[m], dec[m], fp[m], efp[m]


def _fetch_e3_fits(center, radius_deg):  # pragma: no cover - network
    """Epoch 3 peak fluxes: read the NRAO QL3.1 + QL3.2 FITS catalogues, region + sidelobe cuts."""
    from astropy.io import fits

    ras, decs, fps, efps = [], [], [], []
    for url in VLASS_BULK_URLS[3]:
        with fits.open(_cached_download(url)) as hd:
            d = hd[1].data
            scode = np.asarray(d["S_Code"]).astype(str)
            keep = (np.asarray(d["Flag"]) == 0) & (scode != "E")
            ras.append(np.asarray(d["RA"], float)[keep])
            decs.append(np.asarray(d["DEC"], float)[keep])
            fps.append(np.asarray(d["Peak_flux"], float)[keep])
            efps.append(np.asarray(d["E_Peak_flux"], float)[keep])
    ra, dec, fp, efp = (np.concatenate(x) for x in (ras, decs, fps, efps))
    m = _cone_mask(ra, dec, center, radius_deg)
    return ra[m], dec[m], fp[m], efp[m]


_EPOCH_FETCHERS = {1: _fetch_e1_tap, 2: _fetch_e2_csv, 3: _fetch_e3_fits}


def fetch_vlass_epoch(epoch: int, center: tuple[float, float], radius_deg: float):
    """Fetch one VLASS epoch's quality-cut peak fluxes in a cone: ``(ra, dec, peak_mJy, e_peak_mJy)``.

    Epoch 1 via VizieR TAP; Epochs 2-3 via the cached NRAO bulk catalogues. Fluxes are the raw
    Quick-Look values (mJy/beam) before the per-epoch scale correction (:func:`apply_flux_scale`).
    """
    if epoch not in _EPOCH_FETCHERS:
        raise ValueError(f"unsupported VLASS epoch {epoch!r}; supported: {sorted(_EPOCH_FETCHERS)}")
    return _EPOCH_FETCHERS[epoch](center, radius_deg)  # pragma: no cover - network


def _fetch_and_match(
    center: tuple[float, float], radius_deg: float, epochs: tuple[int, ...]
):  # pragma: no cover - network
    """Fetch each requested epoch, flux-correct it, and cross-match onto the first epoch."""
    ra_l, dec_l, flux_l, err_l = [], [], [], []
    for e in epochs:
        ra, dec, fp, efp = fetch_vlass_epoch(e, center, radius_deg)
        # An empty epoch means a schema/filter mismatch silently dropped everything (not a real
        # empty sky), which would yield zero candidates with no error -- fail loudly instead.
        if ra.size == 0:
            raise ValueError(f"epoch {e}: no components fetched in cone (schema/filter mismatch?)")
        # Guard against a per-epoch flux-unit error (e.g. micro-Jy vs mJy): a flux-limited VLASS
        # catalogue has a median peak flux of order 1 mJy, so a median far outside [0.01, 1000] mJy
        # means the units are wrong and any cross-epoch comparison would be meaningless.
        med = float(np.nanmedian(fp))
        if not (1e-2 < med < 1e3):
            raise ValueError(f"epoch {e}: implausible median peak flux {med:.3g} mJy (unit error?)")
        fp, efp = apply_flux_scale(e, fp, efp)
        ra_l.append(ra)
        dec_l.append(dec)
        flux_l.append(fp)
        err_l.append(efp)
    return crossmatch_epochs(ra_l, dec_l, flux_l, err_l)


def _simbad_nearest(c, radius_arcsec: float) -> dict:  # pragma: no cover - network
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    try:
        from astroquery.simbad import Simbad

        s = Simbad()
        s.add_votable_fields("otype")
        t = s.query_region(c, radius=radius_arcsec * u.arcsec)
        if t is None or not len(t):
            return {}
        mc = SkyCoord(np.asarray(t["ra"], float) * u.deg, np.asarray(t["dec"], float) * u.deg)
        seps = c.separation(mc).arcsec
        i = int(np.argmin(seps))
        return {
            "simbad_name": str(t["main_id"][i]),
            "simbad_type": str(t["otype"][i]),
            "simbad_sep": round(float(seps[i]), 2),
        }
    except Exception:
        return {}


def _ned_nearest(c, radius_arcsec: float) -> dict:  # pragma: no cover - network
    from astropy import units as u

    try:
        from astroquery.ipac.ned import Ned

        t = Ned.query_region(c, radius=radius_arcsec * u.arcsec)
        if t is None or not len(t):
            return {}
        sep = np.asarray(t["Separation"], float) * 60.0  # NED Separation is arcmin -> arcsec
        i = int(np.argmin(sep))
        z = t["Redshift"][i]
        return {
            "ned_name": str(t["Object Name"][i]),
            "ned_type": str(t["Type"][i]),
            "ned_z": float(z) if z is not None and np.isfinite(float(z)) else float("nan"),
            "ned_sep": round(float(sep[i]), 2),
        }
    except Exception:
        return {}


def vet_candidates(
    ra: np.ndarray, dec: np.ndarray, *, radius_arcsec: float = 5.0
) -> list[dict]:  # pragma: no cover - network
    """Nearest SIMBAD + NED counterpart for each candidate position (the GATE-2 cross-check).

    A known AGN/blazar/quasar/radio-star counterpart supports a real variable; *no* catalogued
    counterpart for a single-epoch brightening is a red flag for a Quick-Look artefact, not evidence
    of a new transient. Returns one dict per candidate (name/type/redshift/separation, blank if none).
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    out = []
    for r, d in zip(np.asarray(ra, float), np.asarray(dec, float), strict=True):
        c = SkyCoord(r * u.deg, d * u.deg)
        rec = {"ra": float(r), "dec": float(d)}
        rec.update(_simbad_nearest(c, radius_arcsec))
        rec.update(_ned_nearest(c, radius_arcsec))
        out.append(rec)
    return out


def fetch_vlass_cutout(
    ra: float, dec: float, epoch: int, *, size_arcmin: float = 1.5
):  # pragma: no cover - network
    """VLASS Quick-Look image cutout for one epoch around ``(ra, dec)``: returns ``(image_mJy, wcs)``.

    Uses the CADC SODA cutout service (collection ``VLASS``). The data is peak brightness in
    mJy/beam. This is the ground truth for vetting catalogue variability candidates.
    """
    import io
    import re

    import requests
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.nddata import Cutout2D
    from astropy.wcs import WCS
    from astroquery.cadc import Cadc

    c = SkyCoord(ra * u.deg, dec * u.deg)
    rad = (size_arcmin / 60.0) * u.deg
    cadc = Cadc()
    res = cadc.query_region(c, radius=rad, collection="VLASS")
    ql = [
        u_
        for u_ in cadc.get_image_list(res, c, rad)
        if ".ql." in u_ and re.search(rf"VLASS{epoch}\.", u_)
    ]
    if not ql:
        raise ValueError(f"no VLASS Quick-Look image for epoch {epoch} at ({ra}, {dec})")
    with fits.open(io.BytesIO(requests.get(ql[0], timeout=180).content)) as hd:
        img = np.squeeze(np.asarray(hd[0].data, float)) * 1e3  # Jy/beam -> mJy/beam
        w = WCS(hd[0].header).celestial
    cut = Cutout2D(img, c, size_arcmin * u.arcmin, wcs=w)
    return cut.data, cut.wcs


def image_lightcurve(
    ra: float, dec: float, *, epochs: tuple[int, ...] = (1, 2, 3)
) -> np.ndarray:  # pragma: no cover - network
    """Peak image brightness (mJy/beam) at ``(ra, dec)`` per epoch — the ground-truth variability check.

    A genuine variable varies in these *image* peaks; a flat image light curve under a varying
    *catalogue* light curve means the catalogue "variability" is a component-extraction artefact
    (deblending, cross-match miss, fit inconsistency), not real. ``nan`` where no image is available.
    """
    out = []
    for e in epochs:
        try:
            data, _ = fetch_vlass_cutout(ra, dec, e, size_arcmin=0.5)
            out.append(float(np.nanmax(data)))
        except Exception:
            out.append(float("nan"))
    return np.asarray(out)


def measure_image_flux(
    image: np.ndarray,
    wcs,
    ra: float,
    dec: float,
    *,
    search_arcsec: float = 4.0,
    rms_annulus_arcsec: tuple[float, float] = (15.0, 45.0),
) -> tuple[float, float, float]:
    """Forced peak photometry at a fixed ``(ra, dec)``: brightest pixel within ``search_arcsec``.

    Returns ``(peak, rms, offset_arcsec)`` in the image's units. The small search box absorbs
    sub-beam astrometric shifts without grabbing a neighbour; ``rms`` (local noise from an annulus)
    is the per-epoch error. Measuring at the *same locked position* in every epoch is forced
    photometry: it is immune to the deblending and cross-match failures that corrupt the catalogue
    light curve, so a flat forced light curve exposes a spurious catalogue "variable".
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.wcs.utils import proj_plane_pixel_scales

    image = np.asarray(image, dtype=float)
    px, py = wcs.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    scale = float(np.mean(proj_plane_pixel_scales(wcs)) * 3600.0)  # arcsec/pixel
    ny, nx = image.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(xx - float(px), yy - float(py)) * scale
    near = (rr <= search_arcsec) & np.isfinite(image)
    if not near.any():
        return float("nan"), float("nan"), float("nan")
    region = np.where(near, image, -np.inf)
    iy, ix = np.unravel_index(int(np.argmax(region)), region.shape)
    peak = float(image[iy, ix])
    offset = float(np.hypot(ix - float(px), iy - float(py)) * scale)
    ann = image[(rr > rms_annulus_arcsec[0]) & (rr < rms_annulus_arcsec[1]) & np.isfinite(image)]
    rms = float(np.std(ann)) if ann.size > 20 else float("nan")
    return peak, rms, offset


def forced_photometry(
    ra: float, dec: float, *, epochs: tuple[int, ...] = (1, 2, 3), search_arcsec: float = 4.0
) -> tuple[np.ndarray, np.ndarray]:  # pragma: no cover - network
    """Per-epoch forced peak flux + error (mJy/beam) at a locked position; nan where no image."""
    flux, err = [], []
    for e in epochs:
        try:
            data, w = fetch_vlass_cutout(ra, dec, e, size_arcmin=1.5)
            f, r, _ = measure_image_flux(data, w, ra, dec, search_arcsec=search_arcsec)
            flux.append(f)
            err.append(r if (np.isfinite(r) and r > 0) else 0.1 * abs(f))
        except Exception:
            flux.append(float("nan"))
            err.append(float("nan"))
    return np.asarray(flux), np.asarray(err)


def confirm_candidates(
    ra: np.ndarray,
    dec: np.ndarray,
    *,
    epochs: tuple[int, ...] = (1, 2, 3),
    p_max: float = 0.01,
    v_min: float = 0.3,
) -> list[dict]:  # pragma: no cover - network
    """Forced-photometry confirmation of variability candidates against the images (auto image-vetting).

    For each candidate, measures the forced light curve, recomputes $\\eta$/$V$ on it, and marks it
    ``confirmed`` only if the *image* light curve is significantly variable (p < ``p_max``) with real
    amplitude ($V$ > ``v_min``). Rejects the catalogue artefacts (flat forced light curve) automatically.
    """
    out = []
    for r, d in zip(np.asarray(ra, float), np.asarray(dec, float), strict=True):
        flux, err = forced_photometry(r, d, epochs=epochs)
        good = np.isfinite(flux) & np.isfinite(err) & (err > 0)
        if good.sum() >= 2:
            m = variability_metrics(flux[good], err[good])
            confirmed = bool(m.p_value < p_max and m.v > v_min)
            rec = {
                "ra": float(r),
                "dec": float(d),
                "forced_flux": [round(float(x), 3) if np.isfinite(x) else None for x in flux],
                "forced_eta": round(m.eta, 2),
                "forced_v": round(m.v, 3),
                "forced_p": m.p_value,
                "confirmed": confirmed,
            }
        else:
            rec = {"ra": float(r), "dec": float(d), "forced_flux": None, "confirmed": False}
        out.append(rec)
    return out


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
    p.add_argument("--ra", type=float, help="cone-centre RA (deg) for a real-data run")
    p.add_argument("--dec", type=float, help="cone-centre Dec (deg) for a real-data run")
    p.add_argument("--radius", type=float, default=1.0, help="cone radius (deg)")
    p.add_argument("--epochs", default="1,2,3", help="comma-separated VLASS epochs, e.g. 1,2,3")
    args = p.parse_args(argv)
    center = None if (args.offline or args.ra is None) else (args.ra, args.dec)
    epochs = tuple(int(e) for e in args.epochs.split(","))
    print(
        json.dumps(
            run(
                args.out,
                offline=args.offline,
                center=center,
                radius_deg=args.radius,
                epochs=epochs,
            ),
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
