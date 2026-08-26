"""Solar type III radio bursts: fit the frequency drift, recover the exciter (beam) speed.

A type III burst is a beam of ~keV electrons climbing the corona along open field; it excites radio
emission at the local plasma frequency :math:`f_p\\propto\\sqrt{n_e}` (or its harmonic), and because
density falls with height the burst drifts fast from high to low frequency. Given a coronal density
model, the frequency drift becomes a **height-versus-time track** whose slope is the beam speed --
classically a sizeable fraction of :math:`c`.

This module fits the drift ridge in an e-Callisto dynamic spectrum (Benz et al. 2009; open FITS, no
auth) and inverts it to an exciter speed, reusing the course's coronal-physics helpers
(``jansky.solar.density_from_plasma_frequency`` / ``newkirk_radius``; Newkirk 1961). Pure NumPy/SciPy
with a synthetic offline fixture built from the same forward model, so a clean burst round-trips. The
honest systematics -- fundamental vs harmonic (a factor of two in density), the Newkirk fold factor,
and projection -- are reported, not hidden.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "C_KMS",
    "background_subtract",
    "detect_burst_ridge",
    "exciter_speed",
    "fetch_ecallisto",
    "find_burst_window",
    "fit_drift_rate",
    "run",
    "synthetic_burst",
    "run_candidates",
]

C_KMS = 299792.458  # speed of light (km/s)

# The recover-a-known: a clean, isolated type III from the Monstein/e-Callisto SGD burst list
# (110914 1150.1-1151.3 III B quality-1, 21-78 MHz), recorded at Birr (BIR). Chosen because its
# height-time track is coherent (converged clipped fit R^2 ~ 0.8), unlike the flare-storm
# intervals. Harmonic emission and quiet (1x) Newkirk are the default knobs; the paper reports
# the full systematic grid. pad_s here matches the committed evidence (the run() default);
# an earlier 5.0 disagreed with --recover's pinned 10.0, the stale-grid incident's root cause.
RECOVER_EVENT = {
    "station": "BIR",
    "date": "20110914",
    "hhmm": "1150",
    "harmonic": 2,
    "fold": 1.0,
    "pad_s": 10.0,
}

# The full candidate list behind the event selection, plus the flare-storm negative control:
# the paper's crux is that R^2 separates one clean drift from the rest, so the rejected
# candidates' numbers must be committed evidence at the SAME parameterization as the headline,
# not a prose table from a superseded run (referee finding).
CANDIDATE_EVENTS = (
    {"station": "BIR", "date": "20110914", "hhmm": "1150", "role": "accepted"},
    {"station": "BIR", "date": "20110705", "hhmm": "1054", "role": "candidate"},
    {"station": "BIR", "date": "20110716", "hhmm": "1310", "role": "candidate"},
    {"station": "BIR", "date": "20110919", "hhmm": "0743", "role": "candidate"},
    {"station": "BIR", "date": "20110809", "hhmm": "0805", "role": "storm control (X6.9 flare)"},
)


def run_candidates(out: str = ".") -> dict:  # pragma: no cover - network
    """Run every candidate (and the storm control) at the headline parameterization; commit it.

    Writes ``results/solarbursts_candidates.json`` with the same fields for each event, so the
    event-selection function -- the paper's stated crux -- is auditable from evidence.
    """
    from pathlib import Path

    from .report import write_results

    rows = []
    for ev in CANDIDATE_EVENTS:
        burst = fetch_ecallisto(ev["station"], ev["date"], ev["hhmm"])
        window = find_burst_window(burst["data"], burst["times"], pad_s=10.0)
        rf, rt = detect_burst_ridge(burst["data"], burst["freqs"], burst["times"], window=window)
        spd = exciter_speed(rf, rt, harmonic=2, fold=1.0)
        spd.pop("_keep", None)
        rows.append(
            {
                "event": f"{ev['station']} {ev['date']} {ev['hhmm']}",
                "role": ev["role"],
                "n_ridge": int(rf.size),
                "n_used": spd["n_used"],
                "r2": round(spd["r2"], 3) if np.isfinite(spd["r2"]) else None,
                "drift_mhz_per_s": round(fit_drift_rate(rf, rt), 3),
                "speed_c": round(spd["speed_c"], 4) if np.isfinite(spd["speed_c"]) else None,
            }
        )
    payload = {
        "source": "e-Callisto candidate selection at the headline parameterization "
        "(pad_s=10, snr=5, converged clip)",
        "candidates": rows,
    }
    write_results(payload, Path(out) / "results" / "solarbursts_candidates.json")
    # Fold the selection numbers into the main metrics and regenerate the macros through the
    # normal writer: _write_macros always emits the candidate names (placeholders when the
    # stage has not run), otherwise the merge guard would drop them on the next run() rebuild
    # -- the same accumulate-names-not-just-values trap the vlass archival stage hit.
    import json as _json

    mpath = Path(out) / "results" / "solarbursts_metrics.json"
    metrics = _json.loads(mpath.read_text()) if mpath.is_file() else {"source": "?"}
    metrics["candidate_selection"] = rows
    write_results(metrics, mpath)
    _write_macros(metrics, Path(out) / "papers" / "solarbursts" / "generated" / "macros.tex")
    return payload


def synthetic_burst(
    *,
    speed_c: float = 0.3,
    r0_rsun: float = 1.6,
    harmonic: int = 2,
    fold: float = 1.0,
    f_lo_mhz: float = 20.0,
    f_hi_mhz: float = 90.0,
    n_freq: int = 200,
    duration_s: float = 9.0,
    n_time: int = 400,
    width_mhz: float = 2.0,
    amp: float = 12.0,
    noise: float = 1.0,
    seed: int = 0,
) -> dict:
    """Synthetic type III dynamic spectrum with an injected exciter of known speed.

    The beam starts at heliocentric radius ``r0_rsun`` and climbs at ``speed_c`` times the speed of
    light; at each instant the Newkirk density at its radius sets the plasma frequency and hence the
    (harmonic) emission frequency, tracing a fast high-to-low drift ridge. A Gaussian of width
    ``width_mhz`` is laid along that ridge on a noisy background. Because this forward model uses the
    *same* ``jansky.solar`` mapping the analysis inverts, a clean burst recovers ``speed_c`` exactly.
    Returns a dict with ``data`` (n_freq x n_time), ``freqs`` (MHz, descending), ``times`` (s), and the
    injected ``truth_speed_c`` / ``harmonic`` / ``fold``.
    """
    from jansky import solar

    rng = np.random.default_rng(seed)
    freqs = np.linspace(f_hi_mhz, f_lo_mhz, n_freq)  # descending
    times = np.linspace(0.0, duration_s, n_time)
    v_rsun_per_s = speed_c * C_KMS / solar.R_SUN_KM
    r_t = r0_rsun + v_rsun_per_s * times
    fp_t = solar.plasma_frequency(solar.newkirk_density(r_t, fold))
    f_ridge = harmonic * fp_t  # emission frequency vs time (MHz), decreasing

    data = rng.normal(0.0, noise, (n_freq, n_time))
    for j, fr in enumerate(f_ridge):
        data[:, j] += amp * np.exp(-0.5 * ((freqs - fr) / width_mhz) ** 2)
    return {
        "data": data,
        "freqs": freqs,
        "times": times,
        "truth_speed_c": speed_c,
        "harmonic": harmonic,
        "fold": fold,
    }


def background_subtract(data: np.ndarray) -> np.ndarray:
    """Remove each frequency channel's baseline (its median over time).

    e-Callisto raw spectra carry a strong, channel-dependent instrumental offset; subtracting the
    per-channel median over the sweep leaves the transient burst above a near-zero background.
    """
    arr = np.asarray(data, float)
    return arr - np.median(arr, axis=1, keepdims=True)


def find_burst_window(data: np.ndarray, times: np.ndarray, *, pad_s: float = 10.0) -> np.ndarray:
    """Boolean time mask around the burst: the band-integrated power peak, +/- ``pad_s`` seconds.

    e-Callisto files are 15 minutes long but a type III lasts only seconds; the brightest
    band-integrated (background-subtracted) sample locates the burst, and we keep a short window
    around it so the ridge detector is not swamped by quiescent background and RFI.
    """
    clean = background_subtract(data)
    times = np.asarray(times, float)
    t0 = times[int(np.argmax(clean.sum(axis=0)))]
    return (times >= t0 - pad_s) & (times <= t0 + pad_s)


def detect_burst_ridge(
    data: np.ndarray,
    freqs: np.ndarray,
    times: np.ndarray,
    *,
    window: np.ndarray | None = None,
    snr_threshold: float = 5.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Find the drift ridge: for each frequency channel, the time of peak intensity.

    A type III sweeps through each channel once, so the peak *time* per channel traces the drift
    ridge densely (one point per channel where the burst is detectable). Background-subtracts; the
    per-channel noise is the robust MAD over the **full** sweep (a seconds-long burst does not dominate
    a 15-minute series); the peak is then sought within ``window`` (a boolean time mask from
    :func:`find_burst_window`, or the whole sweep if ``None``). Channels whose in-window peak exceeds
    ``snr_threshold`` times their noise are kept. Returns ``(ridge_freqs_mhz, ridge_times_s)``.
    """
    clean = background_subtract(data)
    freqs = np.asarray(freqs, float)
    times = np.asarray(times, float)
    noise_ch = 1.4826 * np.median(np.abs(clean - np.median(clean, axis=1, keepdims=True)), axis=1)
    sub = clean[:, window] if window is not None else clean
    sub_times = times[window] if window is not None else times
    t_local = np.argmax(sub, axis=1)  # peak time per channel, within the window
    peak = sub[np.arange(sub.shape[0]), t_local]
    nz = np.where(noise_ch > 0, noise_ch, np.inf)
    keep = peak > snr_threshold * nz
    return freqs[keep], sub_times[t_local[keep]]


def _robust_linfit(
    x: np.ndarray, y: np.ndarray, *, n_iter: int = 3, sigma: float = 3.0, converge: bool = False
) -> tuple[float, float, np.ndarray]:
    """Iterative sigma-clipped linear fit ``y = m x + b``; returns ``(m, b, keep_mask)``.

    Re-fits after rejecting points more than ``sigma`` residual-standard-deviations from the line,
    so a few RFI-corrupted ridge points (or a second overlapping burst) do not bias the slope.

    ``converge=True`` iterates to a fixed point (the mask stops changing) and guarantees the
    returned slope was fitted on exactly the returned mask. A referee traced the solarbursts
    headline moving 0.111--0.147 c purely over the legacy hard-coded ``n_iter=3``, whose returned
    slope and mask additionally came from *different* iterations. The legacy fixed-count mode is
    kept as the default because four other slices' committed evidence was produced with it; each
    should migrate when its own evidence is re-run.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    keep = np.ones(x.size, dtype=bool)
    m = b = float("nan")
    limit = 100 if converge else n_iter
    for _ in range(limit):
        if keep.sum() < 3:
            break
        m, b = np.polyfit(x[keep], y[keep], 1)
        resid = y - (m * x + b)
        s = np.std(resid[keep])
        if s == 0:
            break
        new = np.abs(resid) < sigma * s
        if converge and bool(np.all(new == keep)):
            break
        keep = new
    if converge and keep.sum() >= 3:
        m, b = np.polyfit(x[keep], y[keep], 1)  # slope and mask from the same, final fit
    return float(m), float(b), keep


def fit_drift_rate(ridge_freqs: np.ndarray, ridge_times: np.ndarray) -> float:
    """Representative frequency drift rate ``df/dt`` (MHz/s) from a linear fit to the ridge.

    Type III drift is steep and negative (frequency falls with time). The single linear slope is a
    representative value over the band; the physical exciter speed comes from :func:`exciter_speed`.
    """
    f = np.asarray(ridge_freqs, float)
    t = np.asarray(ridge_times, float)
    if f.size < 2 or np.ptp(t) == 0:
        return float("nan")
    slope, _ = np.polyfit(t, f, 1)
    return float(slope)


def exciter_speed(
    ridge_freqs: np.ndarray,
    ridge_times: np.ndarray,
    *,
    harmonic: int = 2,
    fold: float = 1.0,
    clip_sigma: float = 3.0,
) -> dict:
    """Exciter (beam) speed from the drift ridge, via the Newkirk coronal density model.

    Each ridge frequency is taken as the ``harmonic`` of the local plasma frequency, so
    :math:`f_p = f/\\mathrm{harmonic}` gives the density
    (``jansky.solar.density_from_plasma_frequency``) and hence the heliocentric radius
    (``newkirk_radius`` with the active-region ``fold``). Fitting radius versus time gives the radial
    speed; returned in km/s and in units of :math:`c`, with the radius range covered. Harmonic and
    ``fold`` are the two model knobs the result depends on (see the caveats).
    """
    from jansky import solar

    f = np.asarray(ridge_freqs, float)
    t = np.asarray(ridge_times, float)
    fp = f / harmonic
    ne = solar.density_from_plasma_frequency(fp)
    r = solar.newkirk_radius(ne, fold)  # heliocentric radius (R_sun)
    nan = float("nan")
    if r.size < 3 or np.ptp(t) == 0:
        return {
            "speed_kms": nan,
            "speed_c": nan,
            "r_lo": nan,
            "r_hi": nan,
            "n_points": int(r.size),
            "n_used": 0,
            "r2": nan,
        }
    # converged clipped fit: the slope, mask and R^2 all come from the same, final iteration
    dr_dt, icpt, keep = _robust_linfit(t, r, sigma=clip_sigma, converge=True)
    speed_kms = abs(dr_dt) * solar.R_SUN_KM
    # coefficient of determination on the kept (coherent) ridge: ~1 for one clean burst, low for a storm
    rk, tk = r[keep], t[keep]
    ss_res = float(np.sum((rk - (dr_dt * tk + icpt)) ** 2))
    ss_tot = float(np.sum((rk - np.mean(rk)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else nan
    return {
        "speed_kms": float(speed_kms),
        "speed_c": float(speed_kms / C_KMS),
        "r_lo": float(np.min(rk)),
        "r_hi": float(np.max(rk)),
        "n_points": int(r.size),
        "n_used": int(keep.sum()),
        "r2": float(r2),
        "clip_sigma": float(clip_sigma),
        "_keep": keep,  # the final fit's mask (stripped before any JSON write)
        "_slope": float(dr_dt),
        "_icpt": float(icpt),
    }


def _isolated_channels(ridge_freqs: np.ndarray, *, gap_mhz: float = 10.0) -> np.ndarray:
    """True for ridge channels more than ``gap_mhz`` from their nearest ridge neighbour.

    Both of the committed ridge's quoted band extremes rest on single isolated channels (78.94
    MHz sits 16.5 MHz above the next detection and drifts the wrong way for a type III; 10.0
    MHz is the clipped low outlier), so the drop-isolated sensitivity variant is committed with
    the headline.
    """
    f = np.sort(np.asarray(ridge_freqs, float))
    order = np.argsort(np.asarray(ridge_freqs, float))
    gaps = np.full(f.size, np.inf)
    if f.size >= 2:
        d = np.diff(f)
        gaps[:-1] = d
        gaps2 = np.full(f.size, np.inf)
        gaps2[1:] = d
        gaps = np.minimum(gaps, gaps2)
    iso_sorted = gaps > gap_mhz
    iso = np.zeros(f.size, dtype=bool)
    iso[order] = iso_sorted
    return iso


#: The systematics grid the paper quotes: emission interpretation x Newkirk fold factor. The three
#: points bracket the model dependence -- fundamental vs harmonic is a factor of two in density,
#: and fold 4 is the active-region enhancement.
SPEED_GRID = ((1, 1.0), (2, 1.0), (2, 4.0))


def speed_grid(ridge_freqs: np.ndarray, ridge_times: np.ndarray) -> list[dict]:
    """Exciter speed over the harmonic x fold systematics grid, from one fitted ridge.

    The ridge itself is model-independent; only the frequency-to-radius mapping changes with the
    grid point, so one detection yields the whole grid. The paper's Results paragraph used to
    hand-type these three speeds from a superseded run -- 0.137 c for the harmonic/1x point three
    lines under a macro saying 0.1347 -- because nothing emitted them. Now they are computed from
    the same ridge as the headline number and cannot drift from it.
    """
    out = []
    for harmonic, fold in SPEED_GRID:
        spd = exciter_speed(ridge_freqs, ridge_times, harmonic=harmonic, fold=fold)
        out.append(
            {
                "harmonic": harmonic,
                "fold": fold,
                "speed_c": round(spd["speed_c"], 4) if np.isfinite(spd["speed_c"]) else None,
                "r_lo_rsun": round(spd["r_lo"], 3) if np.isfinite(spd["r_lo"]) else None,
                "r_hi_rsun": round(spd["r_hi"], 3) if np.isfinite(spd["r_hi"]) else None,
                # per-point fit accounting: the grid rows differ in n_used and R^2 too
                "n_used": spd["n_used"],
                "r2": round(spd["r2"], 3) if np.isfinite(spd["r2"]) else None,
            }
        )
    return out


def fetch_ecallisto(
    station: str, date_yyyymmdd: str, hhmm: str, *, filename: str | None = None
) -> dict:  # pragma: no cover - network
    """Fetch + parse one e-Callisto 15-minute dynamic spectrum.

    With ``filename`` given, fetches EXACTLY that file --- the path a caller that already listed
    the day must use. Without it, lists the day-directory and picks the ``station`` file whose
    start time most closely precedes ``hhmm``; note that resolution sends ``HHMM59``-starting
    files to the PREVIOUS file and collapses focus-code siblings onto the first match, so a
    caller iterating a day listing through the ``hhmm`` path can analyse a different file than
    the one it listed (the round-11 ecallisto_pipeline referee measured 64% mismatches on a real
    day). Returns ``data`` (n_freq x n_time), ``freqs`` (MHz), ``times`` (s). No authentication.
    """
    import gzip
    import io
    import re

    import requests
    from astropy.io import fits

    base = "http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto"
    yyyy, mm, dd = date_yyyymmdd[:4], date_yyyymmdd[4:6], date_yyyymmdd[6:8]
    day_url = f"{base}/{yyyy}/{mm}/{dd}/"
    best = filename
    if best is None:
        idx = requests.get(day_url, timeout=60).text
        pat = rf"{re.escape(station)}_{date_yyyymmdd}_([0-9]{{6}})_[0-9]+\.fit\.gz"
        want = int(hhmm) * 100  # HHMM -> HHMM00 seconds-of-day key
        best_dt = None
        for m in re.finditer(pat, idx):
            start = int(m.group(0).split("_")[2])
            if start <= want and (best_dt is None or want - start < best_dt):
                best, best_dt = m.group(0), want - start
        if best is None:
            raise RuntimeError(f"no e-Callisto {station} file near {hhmm} on {date_yyyymmdd}")
    raw = gzip.decompress(requests.get(day_url + best, timeout=120).content)
    with fits.open(io.BytesIO(raw)) as hd:
        data = np.asarray(hd[0].data, float)
        freqs = np.asarray(hd[1].data["FREQUENCY"][0], float)
        times = np.asarray(hd[1].data["TIME"][0], float)
    return {"data": data, "freqs": freqs, "times": times, "file": best}


def run(
    out: str = ".",
    *,
    offline: bool = True,
    station: str | None = None,
    date: str | None = None,
    hhmm: str | None = None,
    harmonic: int = 2,
    fold: float = 1.0,
    pad_s: float = 10.0,
) -> dict:
    """Full slice: fit a type III drift (synthetic or fetched) and report the exciter speed.

    ``pad_s`` sets the half-width of the burst time window (:func:`find_burst_window`); a tight window
    isolates a single type III, while a wide one over a burst storm smears the drift (flagged by a low
    fit ``r2``).
    """
    from pathlib import Path

    if offline or station is None:
        burst = synthetic_burst(harmonic=harmonic, fold=fold)
        source = "synthetic"
        truth: float | None = burst["truth_speed_c"]
    else:  # pragma: no cover - network
        if date is None or hhmm is None:
            raise ValueError("a real run needs --station, --date (YYYYMMDD), and --hhmm")
        burst = fetch_ecallisto(station, date, hhmm)
        source = f"e-Callisto {station} {date} {hhmm}"
        truth = None

    window = find_burst_window(burst["data"], burst["times"], pad_s=pad_s)
    rf, rt = detect_burst_ridge(burst["data"], burst["freqs"], burst["times"], window=window)
    drift = fit_drift_rate(rf, rt)
    spd = exciter_speed(rf, rt, harmonic=harmonic, fold=fold)
    used = spd.pop("_keep", np.ones(rf.size, bool))
    fit_slope, fit_icpt = spd.pop("_slope", float("nan")), spd.pop("_icpt", float("nan"))
    # The band and drift of the USED (clip-surviving) channels — the set the headline speed and
    # R^2 describe. The full-ridge drift is kept separately: quoting the full detection band and
    # the unclipped drift next to "n_used channels, R^2" bound numbers to the wrong point set.
    drift_used = fit_drift_rate(rf[used], rt[used]) if used.sum() >= 2 else float("nan")

    def _spdc(rf_, rt_, **kw) -> float | None:
        v = exciter_speed(rf_, rt_, harmonic=harmonic, fold=fold, **kw)["speed_c"]
        return round(float(v), 4) if np.isfinite(v) else None

    # Analysis-choice sensitivity, committed with the headline: window half-width, clip sigma,
    # and the isolated band-edge channels. The headline is quoted with the min-max spread.
    sens: dict[str, float | None] = {}
    for p in (5.0, 10.0):
        wv = find_burst_window(burst["data"], burst["times"], pad_s=p)
        rfp, rtp = detect_burst_ridge(burst["data"], burst["freqs"], burst["times"], window=wv)
        sens[f"pad_{p:g}s"] = _spdc(rfp, rtp)
    for sg in (2.5, 3.5):
        sens[f"clip_sigma_{sg:g}"] = _spdc(rf, rt, clip_sigma=sg)
    iso = _isolated_channels(rf)
    sens["drop_isolated"] = _spdc(rf[~iso], rt[~iso]) if (~iso).sum() >= 3 else None
    spread = [v for v in sens.values() if v is not None]
    if np.isfinite(spd["speed_c"]):
        spread.append(round(float(spd["speed_c"]), 4))

    metrics: dict = {
        "source": source,
        "n_ridge": int(rf.size),
        "n_used": spd.get("n_used"),
        "r2": round(spd["r2"], 3) if np.isfinite(spd.get("r2", float("nan"))) else None,
        "f_lo_mhz": round(float(np.min(rf)), 2) if rf.size else None,
        "f_hi_mhz": round(float(np.max(rf)), 2) if rf.size else None,
        "fit_f_lo_mhz": round(float(np.min(rf[used])), 2) if used.any() else None,
        "fit_f_hi_mhz": round(float(np.max(rf[used])), 2) if used.any() else None,
        "drift_mhz_per_s": round(drift, 3) if np.isfinite(drift) else None,
        "drift_used_mhz_per_s": round(drift_used, 3) if np.isfinite(drift_used) else None,
        "harmonic": harmonic,
        "fold": fold,
        "pad_s": float(pad_s),
        "snr_threshold": 5.0,
        "clip_sigma": spd.get("clip_sigma"),
        "fit_converged": True,
        "r_lo_rsun": round(spd["r_lo"], 3) if np.isfinite(spd["r_lo"]) else None,
        "r_hi_rsun": round(spd["r_hi"], 3) if np.isfinite(spd["r_hi"]) else None,
        # rounded to the hundreds: six significant figures on a value with a percent-level
        # analysis spread was false precision (referee finding)
        "speed_kms": int(round(spd["speed_kms"], -2)) if np.isfinite(spd["speed_kms"]) else None,
        "speed_c": round(spd["speed_c"], 4) if np.isfinite(spd["speed_c"]) else None,
        "speed_sensitivity": sens,
        "speed_c_min": round(min(spread), 4) if spread else None,
        "speed_c_max": round(max(spread), 4) if spread else None,
    }
    metrics["speed_grid"] = speed_grid(rf, rt)
    if truth is not None:
        metrics["truth_speed_c"] = truth
        if np.isfinite(spd["speed_c"]):
            metrics["recovery_ratio"] = round(spd["speed_c"] / truth, 3) if truth else None

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "solarbursts_metrics.json")
    # The raw fitted ridge, committed alongside the summary: exciter_speed needs it, so without it
    # neither the headline speed nor the grid is recomputable from evidence (the innerrc lesson --
    # a results file omitting the numbers its own headline is computed from).
    import csv as _csv

    with (op / "results" / "solarbursts_ridge.csv").open("w", newline="") as fh:
        fh.write(
            f"# {source}; pad_s={pad_s:g}; snr_threshold=5; used=1 marks channels surviving "
            "the converged clip\n"
        )
        w = _csv.writer(fh)
        w.writerow(["freq_mhz", "time_s", "used"])
        w.writerows(zip(np.round(rf, 4), np.round(rt, 4), used.astype(int), strict=True))
    _figure(
        burst,
        rf,
        rt,
        used,
        fit_slope,
        fit_icpt,
        harmonic,
        fold,
        op / "papers" / "solarbursts" / "figures",
    )
    _write_macros(metrics, op / "papers" / "solarbursts" / "generated" / "macros.tex")
    return metrics


def _figure(burst, rf, rt, used, fit_slope, fit_icpt, harmonic, fold, out_dir) -> None:
    from pathlib import Path

    from jansky import solar

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    clean = background_subtract(burst["data"])
    freqs, times = burst["freqs"], burst["times"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
    ax1.imshow(
        clean,
        origin="upper",
        aspect="auto",
        cmap="inferno",
        extent=[times.min(), times.max(), freqs.min(), freqs.max()],
    )
    ax1.plot(rt, rf, ".", color="cyan", ms=3, label="ridge")
    ax1.set(xlabel="time (s)", ylabel="frequency (MHz)", title="Type III dynamic spectrum")
    ax1.legend(loc="upper right", fontsize=8)
    if rf.size >= 2:
        # Draw the line the caption quotes: the converged clipped fit over the kept points.
        # An earlier figure drew an unweighted fit over all points (including the clipped
        # outliers, unmarked), so the plotted slope was 14% below the reported one and the
        # caption's R^2 belonged to a different line.
        r = solar.newkirk_radius(solar.density_from_plasma_frequency(rf / harmonic), fold)
        ax2.plot(rt[used], r[used], "o", color="C0", ms=3, label="kept")
        if (~used).any():
            ax2.plot(rt[~used], r[~used], "x", color="0.5", ms=5, label="clipped")
        if np.isfinite(fit_slope):
            tk = rt[used]
            tt = np.linspace(tk.min(), tk.max(), 20)
            ax2.plot(tt, fit_slope * tt + fit_icpt, "-", color="C3", lw=1, label="clipped fit")
        ax2.legend(fontsize=7)
        ax2.set(
            xlabel="time (s)", ylabel=r"heliocentric radius ($R_\odot$)", title="Height--time track"
        )
    fig.tight_layout()
    fig.savefig(out / "burst.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.solarbursts._write_macros -- do not edit by hand.",
        rf"\newcommand{{\sbSource}}{{{m['source']}}}",
        rf"\newcommand{{\sbNridge}}{{{m['n_ridge']}}}",
        rf"\newcommand{{\sbNused}}{{{_fmt('n_used')}}}",
        rf"\newcommand{{\sbRsq}}{{{_fmt('r2')}}}",
        rf"\newcommand{{\sbFlo}}{{{_fmt('f_lo_mhz')}}}",
        rf"\newcommand{{\sbFhi}}{{{_fmt('f_hi_mhz')}}}",
        rf"\newcommand{{\sbDrift}}{{{_fmt('drift_mhz_per_s')}}}",
        rf"\newcommand{{\sbHarmonic}}{{{m['harmonic']}}}",
        rf"\newcommand{{\sbFold}}{{{m['fold']}}}",
        rf"\newcommand{{\sbRlo}}{{{_fmt('r_lo_rsun')}}}",
        rf"\newcommand{{\sbRhi}}{{{_fmt('r_hi_rsun')}}}",
        rf"\newcommand{{\sbSpeedKms}}{{{_fmt('speed_kms')}}}",
        rf"\newcommand{{\sbSpeedC}}{{{_fmt('speed_c')}}}",
        rf"\newcommand{{\sbSpeedCMin}}{{{_fmt('speed_c_min')}}}",
        rf"\newcommand{{\sbSpeedCMax}}{{{_fmt('speed_c_max')}}}",
        rf"\newcommand{{\sbFitFlo}}{{{_fmt('fit_f_lo_mhz')}}}",
        rf"\newcommand{{\sbFitFhi}}{{{_fmt('fit_f_hi_mhz')}}}",
        rf"\newcommand{{\sbDriftUsed}}{{{_fmt('drift_used_mhz_per_s')}}}",
        rf"\newcommand{{\sbPad}}{{{_fmt('pad_s')}}}",
        rf"\newcommand{{\sbClipSigma}}{{{_fmt('clip_sigma')}}}",
        # synthetic-leg-only quantities carry the Syn namespace so a real run's '--' cannot be
        # mistaken for a blanked value (they were \sbTruth/\sbRatio, unused and rendering as --)
        rf"\newcommand{{\sbSynTruth}}{{{_fmt('truth_speed_c')}}}",
        rf"\newcommand{{\sbSynRatio}}{{{_fmt('recovery_ratio')}}}",
    ]
    # Candidate-selection macros: always emitted (placeholders until run_candidates has run),
    # so the merge guard never drops them on a plain run() rebuild.
    cand_rows = m.get("candidate_selection") or []
    rejects = [r for r in cand_rows if r.get("role") == "candidate"]
    slots = {
        "CandAccepted": next((r for r in cand_rows if r.get("role") == "accepted"), None),
        "CandOne": rejects[0] if len(rejects) > 0 else None,
        "CandTwo": rejects[1] if len(rejects) > 1 else None,
        "CandThree": rejects[2] if len(rejects) > 2 else None,
        "CandStorm": next((r for r in cand_rows if "storm" in r.get("role", "")), None),
    }
    for name, row in slots.items():
        for suffix, key in (("Rsq", "r2"), ("Drift", "drift_mhz_per_s"), ("SpeedC", "speed_c")):
            val = "--" if row is None or row.get(key) is None else str(row[key])
            lines.append(rf"\newcommand{{\sb{name}{suffix}}}{{{val}}}")
    # The systematics grid, one macro per point, so the Results paragraph cannot hand-type them.
    grid = {(g["harmonic"], g["fold"]): g for g in m.get("speed_grid", [])}
    for (h, f), name in (((1, 1.0), "FundOne"), ((2, 1.0), "HarmOne"), ((2, 4.0), "HarmFour")):
        g = grid.get((h, f), {})
        v = g.get("speed_c")
        lines.append(rf"\newcommand{{\sbGrid{name}}}{{{'--' if v is None else v}}}")
        rl, rh = g.get("r_lo_rsun"), g.get("r_hi_rsun")
        lines.append(rf"\newcommand{{\sbGrid{name}Rhi}}{{{'--' if rh is None else rh}}}")
        lines.append(rf"\newcommand{{\sbGrid{name}Rlo}}{{{'--' if rl is None else rl}}}")
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

    p = argparse.ArgumentParser(
        description="Solar type III burst drift -> exciter speed (e-Callisto)."
    )
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--station", help="e-Callisto station, e.g. BIR")
    p.add_argument("--date", help="YYYYMMDD")
    p.add_argument("--hhmm", help="UT start HHMM")
    p.add_argument("--harmonic", type=int, default=2)
    p.add_argument("--fold", type=float, default=1.0)
    p.add_argument("--pad", type=float, default=10.0, help="burst-window half-width (s)")
    p.add_argument("--recover", action="store_true", help="run the canonical recover-a-known event")
    p.add_argument(
        "--candidates",
        action="store_true",
        help="run all candidate events + the storm control; commit the selection evidence",
    )
    args = p.parse_args(argv)
    if args.candidates:
        print(json.dumps(run_candidates(args.out), indent=2))
        return 0
    if args.recover:  # the canonical event in RECOVER_EVENT, spelled out for the type checker
        metrics = run(
            args.out,
            offline=False,
            station="BIR",
            date="20110914",
            hhmm="1150",
            harmonic=2,
            fold=1.0,
            # 10.0, not 5.0: the committed evidence was produced with the run() default, and
            # pad 5 gives a materially different fit (r2 0.897 vs 0.811, speed 0.1368 vs
            # 0.1347). --recover exists to regenerate the committed result, so it must pin the
            # committed parameterization; the pad sensitivity is recorded in the findings file.
            pad_s=10.0,
        )
    else:
        metrics = run(
            args.out,
            offline=args.offline or not args.station,
            station=args.station,
            date=args.date,
            hhmm=args.hhmm,
            harmonic=args.harmonic,
            fold=args.fold,
            pad_s=args.pad,
        )
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
