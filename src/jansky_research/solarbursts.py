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
]

C_KMS = 299792.458  # speed of light (km/s)


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
    if r.size < 2 or np.ptp(t) == 0:
        return {
            "speed_kms": float("nan"),
            "speed_c": float("nan"),
            "r_lo": float("nan"),
            "r_hi": float("nan"),
            "n_points": int(r.size),
        }
    dr_dt, _ = np.polyfit(t, r, 1)  # R_sun per second
    speed_kms = abs(dr_dt) * solar.R_SUN_KM
    return {
        "speed_kms": float(speed_kms),
        "speed_c": float(speed_kms / C_KMS),
        "r_lo": float(np.min(r)),
        "r_hi": float(np.max(r)),
        "n_points": int(r.size),
    }


def fetch_ecallisto(
    station: str, date_yyyymmdd: str, hhmm: str
) -> dict:  # pragma: no cover - network
    """Fetch + parse one e-Callisto 15-minute dynamic spectrum covering ``hhmm`` on ``date``.

    Lists the public archive day-directory, picks the ``station`` file whose start time most closely
    precedes ``hhmm``, downloads the gzipped FITS, and returns ``data`` (n_freq x n_time), ``freqs``
    (MHz), ``times`` (s). No authentication; the archive is open over HTTP.
    """
    import gzip
    import io
    import re

    import requests
    from astropy.io import fits

    base = "http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto"
    yyyy, mm, dd = date_yyyymmdd[:4], date_yyyymmdd[4:6], date_yyyymmdd[6:8]
    day_url = f"{base}/{yyyy}/{mm}/{dd}/"
    idx = requests.get(day_url, timeout=60).text
    pat = rf"{re.escape(station)}_{date_yyyymmdd}_([0-9]{{6}})_[0-9]+\.fit\.gz"
    want = int(hhmm) * 100  # HHMM -> HHMM00 seconds-of-day key
    best, best_dt = None, None
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
) -> dict:
    """Full slice: fit a type III drift (synthetic or fetched) and report the exciter speed."""
    import json
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

    window = find_burst_window(burst["data"], burst["times"])
    rf, rt = detect_burst_ridge(burst["data"], burst["freqs"], burst["times"], window=window)
    drift = fit_drift_rate(rf, rt)
    spd = exciter_speed(rf, rt, harmonic=harmonic, fold=fold)

    metrics: dict = {
        "source": source,
        "n_ridge": int(rf.size),
        "f_lo_mhz": round(float(np.min(rf)), 2) if rf.size else None,
        "f_hi_mhz": round(float(np.max(rf)), 2) if rf.size else None,
        "drift_mhz_per_s": round(drift, 3) if np.isfinite(drift) else None,
        "harmonic": harmonic,
        "fold": fold,
        "r_lo_rsun": round(spd["r_lo"], 3) if np.isfinite(spd["r_lo"]) else None,
        "r_hi_rsun": round(spd["r_hi"], 3) if np.isfinite(spd["r_hi"]) else None,
        "speed_kms": round(spd["speed_kms"], 1) if np.isfinite(spd["speed_kms"]) else None,
        "speed_c": round(spd["speed_c"], 4) if np.isfinite(spd["speed_c"]) else None,
    }
    if truth is not None:
        metrics["truth_speed_c"] = truth
        if np.isfinite(spd["speed_c"]):
            metrics["recovery_ratio"] = round(spd["speed_c"] / truth, 3) if truth else None

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "solarbursts_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(burst, rf, rt, harmonic, fold, op / "papers" / "solarbursts" / "figures")
    _write_macros(metrics, op / "papers" / "solarbursts" / "generated" / "macros.tex")
    return metrics


def _figure(burst, rf, rt, harmonic, fold, out_dir) -> None:
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
        r = solar.newkirk_radius(solar.density_from_plasma_frequency(rf / harmonic), fold)
        ax2.plot(rt, r, "o", color="C0", ms=3)
        slope, icpt = np.polyfit(rt, r, 1)
        ax2.plot(rt, slope * rt + icpt, "-", color="C3", lw=1)
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
        rf"\newcommand{{\sbFlo}}{{{_fmt('f_lo_mhz')}}}",
        rf"\newcommand{{\sbFhi}}{{{_fmt('f_hi_mhz')}}}",
        rf"\newcommand{{\sbDrift}}{{{_fmt('drift_mhz_per_s')}}}",
        rf"\newcommand{{\sbHarmonic}}{{{m['harmonic']}}}",
        rf"\newcommand{{\sbFold}}{{{m['fold']}}}",
        rf"\newcommand{{\sbRlo}}{{{_fmt('r_lo_rsun')}}}",
        rf"\newcommand{{\sbRhi}}{{{_fmt('r_hi_rsun')}}}",
        rf"\newcommand{{\sbSpeedKms}}{{{_fmt('speed_kms')}}}",
        rf"\newcommand{{\sbSpeedC}}{{{_fmt('speed_c')}}}",
        rf"\newcommand{{\sbTruth}}{{{_fmt('truth_speed_c')}}}",
        rf"\newcommand{{\sbRatio}}{{{_fmt('recovery_ratio')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


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
    args = p.parse_args(argv)
    metrics = run(
        args.out,
        offline=args.offline or not args.station,
        station=args.station,
        date=args.date,
        hhmm=args.hhmm,
        harmonic=args.harmonic,
        fold=args.fold,
    )
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
