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
    "epoch_confound",
    "fetch_astrogeo",
    "floor_diagnostics",
    "lightcurve_metrics",
    "run",
    "select_variable",
    "sx_index",
    "synthetic_lightcurves",
    "v_sampling_scatter",
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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-source variability metrics from a ``(n_sources, n_epochs)`` flux/error matrix.

    Each row is one source's light curve with ``nan`` for sessions in which it was not measured. For
    every row with at least :data:`MIN_EPOCHS` finite points we compute, via the tested
    ``vlass.variability_metrics``, the significance $\\eta$ (weighted reduced $\\chi^2$), the amplitude
    $V$ (coefficient of variation), the $\\chi^2$ p-value, the epoch count, the mean flux, and the
    noise-debiased modulation index $m_d$ (on a genuinely steady source, $m_d>0$ directly measures how
    much scatter the assumed error model fails to account for). Rows with too few epochs get ``nan``
    metrics (and ``n_epochs`` counts the finite points regardless).
    """
    f = np.asarray(fmat, float)
    e = np.asarray(emat, float)
    n = f.shape[0]
    eta = np.full(n, np.nan)
    v = np.full(n, np.nan)
    pval = np.full(n, np.nan)
    nep = np.zeros(n, dtype=int)
    mean = np.full(n, np.nan)
    md = np.full(n, np.nan)
    for i in range(n):
        ok = np.isfinite(f[i]) & np.isfinite(e[i]) & (e[i] > 0)
        nep[i] = int(ok.sum())
        if nep[i] < MIN_EPOCHS:
            continue
        m = vlass.variability_metrics(f[i, ok], e[i, ok])
        eta[i], v[i], pval[i], mean[i], md[i] = m.eta, m.v, m.p_value, m.mean_flux, m.m_debiased
    return eta, v, pval, nep, mean, md


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


def floor_diagnostics(
    v: np.ndarray, n_epochs: np.ndarray, is_control: np.ndarray, *, min_epochs: int = MIN_EPOCHS
) -> dict:
    """Everything the median-of-controls floor implies but a point estimate hides.

    A median threshold flags a *steady* source half the time by construction --- that is its
    selection function, and it must be stated, not discovered by a referee. This computes: the
    floor and the individual control $V$s; how many controls sit above their own floor; the full
    drop-one jackknife (each control removed in turn, with the floor and the resulting above-floor
    count per variant --- dropping only the highest controls is a check that can only ever add
    sources); and the above-floor count at the control *maximum*, the most conservative threshold
    the controls themselves can supply.
    """
    v = np.asarray(v, float)
    nep = np.asarray(n_epochs)
    ctrl = np.asarray(is_control, bool)
    okc = ctrl & (nep >= min_epochs) & np.isfinite(v)
    okn = ~ctrl & (nep >= min_epochs) & np.isfinite(v)
    cvs = np.sort(v[okc])
    floor = float(np.median(cvs))

    def _n_above(f: float) -> int:
        return int((v[okn] > f).sum())

    jack = []
    for i in range(cvs.size):
        f_i = float(np.median(np.delete(cvs, i)))
        jack.append({"floor": round(f_i, 5), "n_above": _n_above(f_i)})
    jk_floors = [j["floor"] for j in jack]
    jk_counts = [j["n_above"] for j in jack]
    return {
        "floor": round(floor, 5),
        "control_vs": [round(float(x), 5) for x in cvs],
        "n_above": _n_above(floor),
        "n_controls_above_floor": int((v[okc] > floor).sum()),
        "jackknife": jack,
        "jk_floor_lo": min(jk_floors) if jk_floors else None,
        "jk_floor_hi": max(jk_floors) if jk_floors else None,
        "jk_n_above_lo": min(jk_counts) if jk_counts else None,
        "jk_n_above_hi": max(jk_counts) if jk_counts else None,
        "n_above_at_control_max": _n_above(float(cvs[-1])) if cvs.size else None,
    }


def epoch_confound(
    v: np.ndarray,
    n_epochs: np.ndarray,
    is_control: np.ndarray,
    *,
    min_epochs: int = MIN_EPOCHS,
    n_perm: int = 10000,
    seed: int = 0,
) -> dict:
    """How strongly $V$ depends on the number of epochs, and an epoch-matched floor.

    A light curve sampled at more sessions explores more of the source's variability (and more of
    the $(u,v)$-coverage scatter), so $V$ grows with epoch count for both variables and controls.
    Reports the Spearman correlation of $V$ with $\\log_{10} n$ over all testable sources (with a
    permutation p-value), the OLS fit over the same sources (the confound's size), and an
    **epoch-matched floor**: the controls' own OLS trend evaluated at each non-control's epoch
    count, with the resulting above-floor count. Four controls cannot pin a two-parameter trend,
    so the epoch-matched count is a sensitivity variant, not a replacement headline.
    """
    from scipy import stats

    v = np.asarray(v, float)
    nep = np.asarray(n_epochs)
    ctrl = np.asarray(is_control, bool)
    ok = (nep >= min_epochs) & np.isfinite(v)
    x = np.log10(nep[ok].astype(float))
    y = v[ok]
    rho = float(stats.spearmanr(x, y).statistic)
    rng = np.random.default_rng(seed)
    k = 0
    for _ in range(n_perm):
        r = float(stats.spearmanr(x, rng.permutation(y)).statistic)
        if abs(r) >= abs(rho):
            k += 1
    p_perm = (k + 1) / (n_perm + 1)
    b_all, a_all = np.polyfit(x, y, 1)
    okc = ok & ctrl
    okn = ok & ~ctrl
    out: dict = {
        "spearman_rho": round(rho, 3),
        "spearman_p_perm": float(f"{p_perm:.2e}"),
        "ols_intercept": round(float(a_all), 3),
        "ols_slope": round(float(b_all), 3),
    }
    if okc.sum() >= 2:
        bc, ac = np.polyfit(np.log10(nep[okc].astype(float)), v[okc], 1)
        floors = ac + bc * np.log10(nep[okn].astype(float))
        out["ctrl_ols_intercept"] = round(float(ac), 3)
        out["ctrl_ols_slope"] = round(float(bc), 3)
        out["n_above_epoch_matched"] = int((v[okn] > floors).sum())
    return out


def v_sampling_scatter(
    fmat: np.ndarray,
    emat: np.ndarray,
    *,
    min_epochs: int = MIN_EPOCHS,
    n_boot: int = 2000,
    seed: int = 0,
) -> np.ndarray:
    """Per-source bootstrap standard error of $V$ (epochs resampled with replacement).

    $V$ from a handful of epochs carries sampling scatter of its own; a source "above the floor"
    by less than that scatter is not securely above it. This is the per-source term the z-score
    $(V-{\\rm floor})/\\sigma_V$ needs. Deterministic (fixed seed).
    """
    f = np.asarray(fmat, float)
    e = np.asarray(emat, float)
    n = f.shape[0]
    rng = np.random.default_rng(seed)
    se = np.full(n, np.nan)
    for i in range(n):
        ok = np.isfinite(f[i]) & np.isfinite(e[i]) & (e[i] > 0)
        if int(ok.sum()) < min_epochs:
            continue
        fi = f[i, ok]
        idx = rng.integers(0, fi.size, (n_boot, fi.size))
        boots = fi[idx]
        means = boots.mean(axis=1)
        stds = boots.std(axis=1, ddof=1)
        vs = np.where(means > 0, stds / means, np.nan)
        se[i] = float(np.nanstd(vs))
    return se


def synthetic_lightcurves(
    n_sources: int = 400,
    n_epochs: int = 10,
    *,
    frac_variable: float = 0.08,
    var_amp: float = 2.0,
    err_frac: float = 0.07,
    miss_frac: float = 0.25,
    n_controls: int = 0,
    seed: int = 0,
) -> dict:
    """Synthetic dual-band VLBI population: steady sources + an injected variable subset.

    Steady sources have a constant mean flux per band (so $\\eta\\approx1$, $V\\approx$ the measurement
    error); the injected variable fraction gets a single-session flare of relative amplitude ``var_amp``
    (high $\\eta$ and $V$). Each source has a flat-ish S/X index, ``err_frac`` fractional errors, and a
    fraction ``miss_frac`` of sessions randomly missing (``nan``) to mimic uneven VLBI sampling.
    ``n_controls`` flags that many of the *steady* sources as known-steady controls (``is_control``),
    so the control-floor selector --- the discriminant the real run actually uses --- can be exercised
    and its completeness/purity measured offline. Returns a dict with ``flux_x/err_x/flux_s/err_s``
    ``(N, M)`` matrices and the boolean ``is_variable`` / ``is_control`` truths.
    """
    rng = np.random.default_rng(seed)
    n = n_sources
    mean_x = 10.0 ** rng.uniform(-1.0, 0.5, n)  # ~0.1-3 Jy
    alpha = rng.normal(0.0, 0.2, n)  # flat-spectrum compact AGN
    mean_s = mean_x * (NU_S_GHZ / NU_X_GHZ) ** alpha

    is_variable = rng.random(n) < frac_variable
    is_control = np.zeros(n, dtype=bool)
    if n_controls:
        is_control[np.where(~is_variable)[0][:n_controls]] = True
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
        "is_control": is_control,
    }


def floor_fixture_metrics(
    *, var_amp: float = 2.0, n_sources: int = 400, n_controls: int = 4, seed: int = 0
) -> dict:
    """Completeness/purity of the control-floor selector on the synthetic population.

    This is the validation the real run's discriminant actually needs: the earlier fixture only
    exercised the relative $\\log\\eta$--$\\log V$ outlier cut, which the real (curated,
    everything-varies) sample cannot use. ``var_amp`` scales the injected flare, so sweeping it
    measures how completeness degrades as the injected amplitude approaches the floor.
    """
    pop = synthetic_lightcurves(
        n_sources=n_sources, var_amp=var_amp, n_controls=n_controls, seed=seed
    )
    _eta, v, _p, nep, _mean, _md = lightcurve_metrics(pop["flux_x"], pop["err_x"])
    floor, above = variability_floor(v, nep, pop["is_control"])
    truth = pop["is_variable"] & (nep >= MIN_EPOCHS)
    tp = int((above & truth).sum())
    return {
        "var_amp": var_amp,
        "floor": round(float(floor), 5),
        "n_injected": int(truth.sum()),
        "n_selected": int(above.sum()),
        "completeness": round(tp / max(int(truth.sum()), 1), 3),
        "purity": round(tp / max(int(above.sum()), 1), 3),
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
    years: dict[str, list[int]] = {b: [] for b in bands}
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
                ym = re.search(rf"_{b}_([12][0-9]{{3}})_", fn)
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
                    if ym:
                        years[b].append(int(ym.group(1)))
                if pause:
                    time.sleep(pause)
            per[b].append((np.asarray(flux, float), np.asarray(err, float)))
            max_ep[b] = max(max_ep[b], len(flux))
    n = len(sources)
    out: dict = {}
    for b in bands:
        fmat = np.full((n, max(max_ep[b], 1)), np.nan)
        emat = np.full((n, max(max_ep[b], 1)), np.nan)
        for i, (farr, earr) in enumerate(per[b]):
            fmat[i, : farr.size] = farr
            emat[i, : earr.size] = earr
        out[b] = (fmat, emat)
    out["epoch_years"] = {b: (min(y), max(y)) for b, y in years.items() if y}
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
    from pathlib import Path

    names: list[str] | None
    epoch_years: tuple[int, int] | None = None
    if offline or sources is None:
        pop = synthetic_lightcurves(n_controls=4)
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
        epoch_years = data.get("epoch_years", {}).get("X")

    eta, v, pval, nep, mean, md = lightcurve_metrics(pop["flux_x"], pop["err_x"])
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
        else np.asarray(pop.get("is_control", np.zeros(eta.shape, dtype=bool)), bool)
    )
    v_floor, above = variability_floor(v, nep, is_control)
    v_se = v_sampling_scatter(pop["flux_x"], pop["err_x"])
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
    if np.isfinite(v_floor):
        fd = floor_diagnostics(v, nep, is_control)
        ec = epoch_confound(v, nep, is_control)
        metrics["n_controls"] = int(is_control.sum())
        metrics["n_noncontrol"] = int((testable & ~is_control).sum())
        metrics["v_floor"] = round(v_floor, 5)
        metrics["n_above_floor"] = int(above.sum())
        metrics["median_v_control"] = round(float(np.median(v[is_control & testable])), 3)
        # the honest amplitude comparison: ALL testable non-controls, not the above-floor subset
        # (the median of sources selected by v > median-of-controls exceeds the control median for
        # any input whatsoever -- a comparison that cannot fail measures nothing)
        nonc = v[testable & ~is_control]
        metrics["median_v_noncontrol"] = round(float(np.median(nonc)), 3) if nonc.size else None
        metrics["floor_diagnostics"] = fd
        metrics["epoch_confound"] = ec
        # z-scores against the floor with each source's own V sampling scatter
        with np.errstate(invalid="ignore"):
            z = (v - v_floor) / v_se
        metrics["n_above_z3"] = int((above & (z > 3.0)).sum())
        if above.any():
            marg = int(np.argmin(np.where(above, z, np.inf)))
            metrics["marginal_source"] = {
                "name": names[marg] if names else f"row{marg}",
                "v": round(float(v[marg]), 5),
                "v_se": round(float(v_se[marg]), 5),
                "z": round(float(z[marg]), 2),
            }
        # on genuinely steady controls the debiased modulation index is a direct measurement of
        # the scatter the assumed error model fails to account for
        mdc = md[is_control & testable]
        metrics["median_md_control"] = round(float(np.median(mdc)), 3) if mdc.size else None
        below = testable & ~is_control & ~above
        if names and below.any():  # pragma: no cover - network run only
            metrics["below_floor"] = [
                {
                    "name": names[i],
                    "v": round(float(v[i]), 5),
                    "n_epochs": int(nep[i]),
                    "p_value": float(f"{pval[i]:.2e}"),
                }
                for i in np.where(below)[0]
            ]
    if epoch_years is not None:  # pragma: no cover - network run only
        metrics["epoch_year_min"], metrics["epoch_year_max"] = epoch_years
    # the control-floor selector's own validation, measured on the synthetic population at the
    # default injected amplitude and down a ladder approaching the floor (deterministic; computed
    # in BOTH modes so a real run emits the Syn macros the Methods section cites)
    metrics["syn_floor_validation"] = floor_fixture_metrics()
    metrics["syn_floor_amp_sweep"] = [
        floor_fixture_metrics(var_amp=a) for a in (0.25, 0.5, 1.0, 2.0)
    ]
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
    from .report import _results_are_real, write_results

    # a synthetic run must not clobber the real figure/table (the metrics JSON and macros have
    # their own merge guards; the binary artifacts do not, so check the provenance ourselves)
    json_path = op / "results" / "vlbi_metrics.json"
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
    if names is not None:  # pragma: no cover - network
        _write_candidates(
            op / "results" / "vlbi_candidates.csv",
            names,
            eta,
            v,
            v_se,
            pval,
            nep,
            mean,
            alpha,
            md,
            above,
            is_control,
        )
        _write_table(
            op / "papers" / "vlbi" / "generated" / "sources.tex",
            names,
            eta,
            v,
            v_se,
            nep,
            mean,
            alpha,
            md,
            above,
            is_control,
        )
    if write_artifacts:
        _figure(eta, v, above, is_control, v_floor, op / "papers" / "vlbi" / "figures")
    _write_macros(metrics, op / "papers" / "vlbi" / "generated" / "macros.tex")
    return metrics


def _write_candidates(
    path, names, eta, v, v_se, pval, nep, mean, alpha, md, above, is_control
) -> None:  # pragma: no cover
    """Write the full variability-ranked table (most significant first), tagging controls/variables.

    $V$ and its sampling error carry five decimals so a tie against the floor is resolvable from the
    committed evidence rather than created by the file's own rounding.
    """
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
                "common_name",
                "control",
                "above_floor",
                "n_epochs",
                "eta",
                "v",
                "v_se",
                "m_debiased",
                "p_value",
                "mean_flux_jy",
                "alpha_sx",
            ]
        )
        for i in rows:
            w.writerow(
                [
                    names[i],
                    VALIDATION_SOURCES.get(names[i], ""),
                    int(bool(is_control[i])),
                    int(bool(above[i])),
                    int(nep[i]),
                    f"{eta[i]:.2f}" if np.isfinite(eta[i]) else "",
                    f"{v[i]:.5f}" if np.isfinite(v[i]) else "",
                    f"{v_se[i]:.5f}" if np.isfinite(v_se[i]) else "",
                    f"{md[i]:.5f}" if np.isfinite(md[i]) else "",
                    f"{pval[i]:.2e}" if np.isfinite(pval[i]) else "",
                    f"{mean[i]:.3f}" if np.isfinite(mean[i]) else "",
                    f"{alpha[i]:.2f}" if np.isfinite(alpha[i]) else "",
                ]
            )


def _write_table(
    path, names, eta, v, v_se, nep, mean, alpha, md, above, is_control
) -> None:  # pragma: no cover - network run only
    """Emit the per-source deluxetable body the paper \\inputs, ranked by V (the paper's statistic)."""
    from pathlib import Path

    rows = sorted(
        range(len(names)), key=lambda i: v[i] if np.isfinite(v[i]) else -np.inf, reverse=True
    )
    lines = [
        "% Auto-generated by jansky_research.vlbi._write_table -- do not edit by hand.",
        r"\begin{deluxetable*}{llrrrrrrc}",
        r"\tablecaption{The curated validation set, ranked by the amplitude $V$. $\sigma_V$ is the"
        r" per-source bootstrap sampling error of $V$; $m_d$ is the noise-debiased modulation index"
        r" under the assumed 5\% error model. \label{tab:sources}}",
        r"\tablehead{\colhead{J2000} & \colhead{Common name} & \colhead{$N_{\rm ep}$} &"
        r" \colhead{$\eta$} & \colhead{$V$} & \colhead{$\sigma_V$} & \colhead{$m_d$} &"
        r" \colhead{$\langle S_X\rangle$ (Jy)} & \colhead{$\alpha_{SX}$}}",
        r"\startdata",
    ]
    for i in rows:
        tag = " (control)" if is_control[i] else ""
        star = r"$^{\dagger}$" if above[i] else ""

        def _f(x, fmt: str) -> str:
            return format(x, fmt) if np.isfinite(x) else r"\nodata"

        lines.append(
            f"{names[i]}{star} & {VALIDATION_SOURCES.get(names[i], '')}{tag} & {int(nep[i])} & "
            f"{_f(eta[i], '.1f')} & {_f(v[i], '.3f')} & {_f(v_se[i], '.3f')} & "
            f"{_f(md[i], '.3f')} & {_f(mean[i], '.2f')} & {_f(alpha[i], '+.2f')} \\\\"
        )
    lines += [
        r"\enddata",
        r"\tablenotetext{\dagger}{Above the control floor (median control $V$).}",
        r"\end{deluxetable*}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


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
    if ctrl.any():
        ax.scatter(
            eta[ok & ctrl],
            v[ok & ctrl],
            s=46,
            marker="s",
            facecolors="none",
            edgecolors="C0",
            label="steady control (CSO)",
        )
    if np.isfinite(v_floor):
        ax.axhline(v_floor, ls="--", c="C0", lw=0.9, label=f"floor V={v_floor:.3f}")
    ax.set(
        xscale="log",
        xlabel=r"$\eta$ (significance vs. constant)",
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

    fd = m.get("floor_diagnostics") or {}
    ec = m.get("epoch_confound") or {}
    marg = m.get("marginal_source") or {}
    syn = m.get("syn_floor_validation") or {}

    def _d(dic: dict, key: str, fmt: str | None = None) -> str:
        val = dic.get(key)
        if val is None:
            return "--"
        return format(val, fmt) if fmt else str(val)

    def _p3(key: str) -> str:
        # a 3-decimal presentation of a 5-decimal metrics value, for prose
        val = m.get(key)
        return "--" if val is None else f"{float(val):.3f}"

    floor = m.get("v_floor")
    eff_pct = "--" if floor is None else f"{100 * float(floor):.0f}"
    ratio = "--" if floor is None else f"{float(floor) / VLBI_CAL_FRAC:.1f}"
    cvs = fd.get("control_vs") or []
    cvs_lo = f"{cvs[0]:.3f}" if cvs else "--"
    cvs_hi = f"{cvs[-1]:.3f}" if cvs else "--"
    below = m.get("below_floor") or [{}]
    below_p = below[0].get("p_value")
    below_p_tex = (
        "--"
        if below_p is None
        else rf"{float(f'{below_p:.1e}'.split('e')[0])}\times10^{{{int(f'{below_p:.1e}'.split('e')[1])}}}"
    )
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
        rf"\newcommand{{\viFloor}}{{{_p3('v_floor')}}}",
        rf"\newcommand{{\viFloorFive}}{{{_fmt('v_floor')}}}",
        rf"\newcommand{{\viNabove}}{{{_fmt('n_above_floor')}}}",
        rf"\newcommand{{\viMedVctrl}}{{{_fmt('median_v_control')}}}",
        rf"\newcommand{{\viMedVnonctrl}}{{{_fmt('median_v_noncontrol')}}}",
        rf"\newcommand{{\viInjected}}{{{_fmt('n_injected_variable')}}}",
        rf"\newcommand{{\viCompleteness}}{{{_fmt('completeness')}}}",
        rf"\newcommand{{\viPurity}}{{{_fmt('purity')}}}",
        rf"\newcommand{{\viTopName}}{{{_tv('name')}}}",
        rf"\newcommand{{\viTopEpochs}}{{{_tv('n_epochs')}}}",
        rf"\newcommand{{\viTopEta}}{{{_tv('eta')}}}",
        rf"\newcommand{{\viTopV}}{{{_tv('v')}}}",
        rf"\newcommand{{\viTopFlux}}{{{_tv('mean_flux_jy')}}}",
        # floor uncertainty / selection function
        rf"\newcommand{{\viCtrlVmin}}{{{cvs_lo}}}",
        rf"\newcommand{{\viCtrlVmax}}{{{cvs_hi}}}",
        rf"\newcommand{{\viNctrlAbove}}{{{_d(fd, 'n_controls_above_floor')}}}",
        rf"\newcommand{{\viFloorJkLo}}{{{_d(fd, 'jk_floor_lo', '.3f') if fd.get('jk_floor_lo') is not None else '--'}}}",
        rf"\newcommand{{\viFloorJkHi}}{{{_d(fd, 'jk_floor_hi', '.3f') if fd.get('jk_floor_hi') is not None else '--'}}}",
        rf"\newcommand{{\viNaboveJkLo}}{{{_d(fd, 'jk_n_above_lo')}}}",
        rf"\newcommand{{\viNaboveJkHi}}{{{_d(fd, 'jk_n_above_hi')}}}",
        rf"\newcommand{{\viNaboveCtrlMax}}{{{_d(fd, 'n_above_at_control_max')}}}",
        # epoch-count confound
        rf"\newcommand{{\viSpearmanRho}}{{{_d(ec, 'spearman_rho')}}}",
        rf"\newcommand{{\viSpearmanP}}{{{_d(ec, 'spearman_p_perm')}}}",
        rf"\newcommand{{\viOlsA}}{{{_d(ec, 'ols_intercept')}}}",
        rf"\newcommand{{\viOlsB}}{{{_d(ec, 'ols_slope')}}}",
        rf"\newcommand{{\viNaboveEpochMatched}}{{{_d(ec, 'n_above_epoch_matched')}}}",
        # sampling scatter / the marginal source / debiased index
        rf"\newcommand{{\viNaboveZthree}}{{{_fmt('n_above_z3')}}}",
        rf"\newcommand{{\viMarginalName}}{{{_d(marg, 'name')}}}",
        rf"\newcommand{{\viMarginalV}}{{{_d(marg, 'v')}}}",
        rf"\newcommand{{\viMarginalZ}}{{{_d(marg, 'z')}}}",
        rf"\newcommand{{\viMedMdCtrl}}{{{_fmt('median_md_control')}}}",
        # derived presentation values (from v_floor; not hand-typed)
        rf"\newcommand{{\viEffScatterPct}}{{{eff_pct}}}",
        rf"\newcommand{{\viScatterRatio}}{{{ratio}}}",
        rf"\newcommand{{\viEpochYearMin}}{{{_fmt('epoch_year_min')}}}",
        rf"\newcommand{{\viEpochYearMax}}{{{_fmt('epoch_year_max')}}}",
        rf"\newcommand{{\viBelowFloorP}}{{{below_p_tex}}}",
        # the control-floor selector's synthetic validation (computed in both modes)
        rf"\newcommand{{\viSynFloorCompleteness}}{{{_d(syn, 'completeness')}}}",
        rf"\newcommand{{\viSynFloorPurity}}{{{_d(syn, 'purity')}}}",
        rf"\newcommand{{\viSynFloorNinj}}{{{_d(syn, 'n_injected')}}}",
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
