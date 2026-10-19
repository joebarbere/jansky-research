"""The within-burst type III drift law from LV614 LOFAR dynamic spectra (plan 101, ``lvdrift``).

Šteinbergs, Kolotkov & Nakariakov (arXiv:2609.26535) fit df/dt = A f^a ACROSS bursts -- one
peak-time drift per burst, binned by peak frequency, flux-weighted -- and get A = 0.78 +/- 0.91,
a = 0.43 +/- 0.31 (Fig. 15). Their 290 public type III arrays (1 s x 195 kHz, 29.3-60.7 MHz)
support a WITHIN-burst measurement: track each burst's arrival time across its channels and fit the
integrated law. This module implements the estimator frozen in ``plans/101-lvdrift-within-burst.md``.

For df/dt = -A f^a the arrival time is t(f) = t0 - g(f; a) / A with
g = (f^(1-a) - f_ref^(1-a)) / (1 - a), which tends to ln(f / f_ref) as a -> 1, so the fit is
smooth through a = 1 (f_ref is the geometric-mean frequency of the kept channels).
"""

from __future__ import annotations

from typing import Any

import numpy as np

__all__ = [
    "SNR_MIN",
    "SMOOTH",
    "MIN_CHANNELS",
    "MIN_FREQ_RATIO",
    "MONO_TOL_S",
    "lofar_freqs",
    "burst_ridge",
    "ridge_quality",
    "fit_index",
    "burst_index",
    "synthetic_burst",
    "shuffle_channels",
    "peak_drift",
    "noise_windows",
    "parse_name",
    "p1_real",
]

# ---- frozen estimator constants (plans/101, 2026-10-18) -----------------------------------------
SNR_MIN = 5.0  # ridge point kept if smoothed peak / sigma_ch >= this
SMOOTH = 3  # boxcar, samples x channels
EDGE = 2  # peak must lie this many samples from the array edges
MIN_CHANNELS = 40
MIN_FREQ_RATIO = 1.5
MONO_TOL_S = 2.0  # type III drift high -> low; arrival times may go backwards by at most this
TIMING_ERR_S = 1.0


def lofar_freqs(sb_min: int = 150, sb_max: int = 311) -> np.ndarray:
    """Channel centre frequencies (MHz) of LOFAR subbands ``sb_min..sb_max`` (200 MHz clock,
    1024 subbands: f = sb x 200 / 1024)."""
    return np.arange(sb_min, sb_max + 1) * 200.0 / 1024.0


def burst_ridge(arr: np.ndarray, freqs: np.ndarray, dt: float = 1.0) -> dict[str, np.ndarray]:
    """Per-channel arrival time (s) of the smoothed peak, kept where peak / sigma >= SNR_MIN.

    Background: per-channel median subtracted; sigma: per-channel MAD x 1.4826 of the raw array;
    smoothing: SMOOTH x SMOOTH boxcar. ``arr`` is (time, channel).
    """
    from scipy.ndimage import uniform_filter

    a = np.asarray(arr, float)
    a = a - np.median(a, axis=0, keepdims=True)
    sig = 1.4826 * np.median(np.abs(a - np.median(a, axis=0, keepdims=True)), axis=0)
    sig = np.where(sig > 0, sig, np.inf)
    sm = uniform_filter(a, size=SMOOTH, mode="nearest")
    # smoothing averages SMOOTH^2 independent samples: the noise of `sm` is sig / SMOOTH
    k = np.argmax(sm, axis=0)
    peak = sm[k, np.arange(a.shape[1])]
    snr = peak / (sig / SMOOTH)
    nt = a.shape[0]
    keep = (snr >= SNR_MIN) & (k >= EDGE) & (k <= nt - 1 - EDGE)
    return {"freqs": np.asarray(freqs, float)[keep], "t_peak": k[keep] * dt, "snr": snr[keep]}


def ridge_quality(freqs: np.ndarray, t_peak: np.ndarray) -> dict[str, Any]:
    """The frozen quality cut: >= MIN_CHANNELS kept, frequency span >= MIN_FREQ_RATIO, and
    arrival times monotonic (later at lower frequency) to within MONO_TOL_S."""
    f = np.asarray(freqs, float)
    t = np.asarray(t_peak, float)
    n = int(f.size)
    if n == 0:
        return {"ok": False, "reason": "no ridge", "n": 0}
    order = np.argsort(-f)  # high -> low frequency
    ts = t[order]
    back = float(np.max(np.maximum.accumulate(ts) - ts)) if n > 1 else 0.0
    ratio = float(f.max() / f.min())
    reasons = []
    if n < MIN_CHANNELS:
        reasons.append(f"channels {n} < {MIN_CHANNELS}")
    if ratio < MIN_FREQ_RATIO:
        reasons.append(f"span {ratio:.2f} < {MIN_FREQ_RATIO}")
    if back > MONO_TOL_S:
        reasons.append(f"non-monotonic by {back:.1f} s")
    return {"ok": not reasons, "reason": "; ".join(reasons), "n": n, "ratio": ratio, "back_s": back}


def _g(f: np.ndarray, a: float, f_ref: float) -> np.ndarray:
    if abs(1.0 - a) < 1e-6:
        return np.log(f / f_ref)
    return (f ** (1.0 - a) - f_ref ** (1.0 - a)) / (1.0 - a)


def fit_index(freqs: np.ndarray, t_peak: np.ndarray) -> dict[str, float]:
    """Least-squares fit of t(f) = t0 - g(f; a) / A (uniform 1 s timing error). Returns a, A
    (MHz s^-1 at f in MHz), t0 and the residual rms; a is bounded to [-3, 5] and a fit on the
    bound is flagged (``a_at_bound``), never reported as a measurement."""
    from scipy.optimize import curve_fit

    f = np.asarray(freqs, float)
    t = np.asarray(t_peak, float)
    f_ref = float(np.exp(np.mean(np.log(f))))

    def model(ff, t0, inv_a, a):
        return t0 - inv_a * _g(ff, a, f_ref)

    best = None
    for a0 in (0.0, 0.5, 1.0, 1.5, 2.0):
        try:
            p, cov = curve_fit(
                model,
                f,
                t,
                p0=(float(np.median(t)), 0.1, a0),
                bounds=([-np.inf, 1e-6, -3.0], [np.inf, np.inf, 5.0]),
                sigma=np.full(f.size, TIMING_ERR_S),
                absolute_sigma=True,
                maxfev=20000,
            )
        except (RuntimeError, ValueError):
            continue
        rms = float(np.sqrt(np.mean((model(f, *p) - t) ** 2)))
        if best is None or rms < best[2]:
            best = (p, cov, rms)
    if best is None:
        return {
            "a": float("nan"),
            "A": float("nan"),
            "t0": float("nan"),
            "rms": float("nan"),
            "a_err": float("nan"),
            "a_at_bound": True,
        }
    p, cov, rms = best
    a = float(p[2])
    a_err = float(np.sqrt(cov[2, 2])) if np.isfinite(cov[2, 2]) and cov[2, 2] >= 0 else float("nan")
    return {
        "a": a,
        "A": float(1.0 / p[1]),
        "t0": float(p[0]),
        "rms": rms,
        "a_err": a_err,
        "a_at_bound": bool(a <= -3.0 + 1e-3 or a >= 5.0 - 1e-3),
        # a degenerate design (e.g. one frequency) leaves a free, and SciPy then reports an error
        # of exactly 0 -- never physical, so zero/undefined errors and a negligible span are flagged
        "unconstrained": bool(not np.isfinite(a_err) or a_err <= 0.0 or np.ptp(np.log(f)) < 1e-3),
    }


def burst_index(arr: np.ndarray, freqs: np.ndarray, dt: float = 1.0) -> dict[str, Any]:
    """The frozen within-burst estimator: ridge -> quality cut -> fit. A burst failing the cut is
    returned with ``ok=False`` and its reason (counted, never dropped silently)."""
    r = burst_ridge(arr, freqs, dt)
    q = ridge_quality(r["freqs"], r["t_peak"])
    out: dict[str, Any] = {"ok": q["ok"], "reason": q["reason"], "n_channels": q["n"]}
    if q["ok"]:
        fit = fit_index(r["freqs"], r["t_peak"])
        out.update(fit)
        if fit["a_at_bound"]:
            out.update(ok=False, reason="fit at the a bound")
        elif fit["unconstrained"]:
            out.update(ok=False, reason="index unconstrained")
    return out


def synthetic_burst(
    freqs: np.ndarray,
    *,
    a: float,
    A: float,
    n_t: int = 120,
    dt: float = 1.0,
    t_hi: float = 30.0,
    width_s: float = 3.0,
    peak: float = 10.0,
    noise: np.ndarray | None = None,
    sigma: float = 1.0,
    seed: int = 0,
) -> np.ndarray:
    """A type III burst obeying df/dt = -A f^a: arrival at the highest frequency at ``t_hi``,
    Gaussian time profile of FWHM ``width_s`` and height ``peak`` (in units of ``sigma``), added to
    ``noise`` (time, channel) or to white noise of ``sigma``."""
    f = np.asarray(freqs, float)
    t_arr = t_hi + (_g(f.max(), a, f.max()) - _g(f, a, f.max())) / A
    t = np.arange(n_t)[:, None] * dt
    s = width_s / 2.3548
    burst = peak * sigma * np.exp(-0.5 * ((t - t_arr[None, :]) / s) ** 2)
    if noise is None:
        noise = np.random.default_rng(seed).normal(0.0, sigma, (n_t, f.size))
    return np.asarray(noise, float)[:n_t] + burst


def shuffle_channels(arr: np.ndarray, seed: int = 0) -> np.ndarray:
    """N1 null: each channel circularly shifted by an independent random number of samples."""
    a = np.asarray(arr, float)
    rng = np.random.default_rng(seed)
    out = np.empty_like(a)
    for j in range(a.shape[1]):
        out[:, j] = np.roll(a[:, j], int(rng.integers(0, a.shape[0])))
    return out


def peak_drift(arr: np.ndarray, freqs: np.ndarray, dt: float = 1.0, window_s: float = 3.0) -> float:
    """R0: the paper's quantity -- the instantaneous drift (MHz s^-1, positive for falling
    frequency) at the time of peak flux, from a linear regression of ridge frequency on time within
    +/- ``window_s`` of that time (their Sect. III.5; the window is not stated there, ours is 3 s)."""
    r = burst_ridge(arr, freqs, dt)
    if r["freqs"].size < 3:
        return float("nan")
    a = np.asarray(arr, float)
    a = a - np.median(a, axis=0, keepdims=True)
    t_pk = float(np.argmax(a.sum(axis=1)) * dt)
    near = np.abs(r["t_peak"] - t_pk) <= window_s
    if near.sum() < 3 or np.ptp(r["t_peak"][near]) == 0:
        return float("nan")
    slope = np.polyfit(r["t_peak"][near], r["freqs"][near], 1)[0]
    return float(-slope)


NOISE_WINDOW_S = 20
NOISE_MARGIN_S = 30  # a window is used only if the array's summed flux peaks at or after this


def noise_windows(arrays: list[np.ndarray], dt: float = 1.0) -> list[np.ndarray]:
    """P1's real-noise canvases: the first ``NOISE_WINDOW_S`` s of every array whose
    background-subtracted summed flux peaks at >= ``NOISE_MARGIN_S`` s (flux timing only)."""
    out = []
    n = int(NOISE_WINDOW_S / dt)
    for a in arrays:
        a = np.asarray(a, float)
        prof = (a - np.median(a, axis=0, keepdims=True)).sum(axis=1)
        if int(np.argmax(prof)) * dt >= NOISE_MARGIN_S and a.shape[0] >= n:
            out.append(a[:n])
    return out


def p1_real(
    windows: list[np.ndarray],
    freqs: np.ndarray,
    *,
    indices: tuple[float, ...] = (0.0, 0.43, 1.0, 1.23, 1.84),
    n_plants: int = 50,
    seed: int = 101,
) -> dict[str, dict[str, float]]:
    """P1 in real noise: plant bursts (S/N 5-50 against each window's median channel sigma, FWHM
    2-6 s, t_hi = 3 s, A = 4/45^a MHz/s) into randomly chosen windows; median recovered a per
    planted index, its bias and the pass flag (|bias| <= 0.15)."""
    rng = np.random.default_rng(seed)
    res = {}
    for a in indices:
        est = []
        fails = 0
        for _ in range(n_plants):
            w = windows[int(rng.integers(len(windows)))]
            w0 = w - np.median(w, axis=0, keepdims=True)
            sig = float(np.median(1.4826 * np.median(np.abs(w0), axis=0)))
            b = synthetic_burst(
                freqs,
                a=a,
                A=4.0 / 45.0**a,
                n_t=w.shape[0],
                t_hi=3.0,
                width_s=float(rng.uniform(2, 6)),
                peak=float(rng.uniform(5, 50)),
                noise=w,
                sigma=sig,
            )
            r = burst_index(b, freqs)
            if r["ok"]:
                est.append(r["a"])
            else:
                fails += 1
        med = float(np.median(est)) if est else float("nan")
        res[f"{a:g}"] = {
            "median": med,
            "bias": med - a,
            "n_ok": len(est),
            "n_failed_cut": fails,
            "passed": bool(abs(med - a) <= 0.15),
        }
    return res


def parse_name(path: str) -> dict[str, Any]:
    """Date, start/stop and channel frequencies from an LV614 file name, e.g.
    ``type3_burst_start_time_2023-04-28_07:02:39_1100_stop_time_1103_subband_min_150_subband_max_311.npy``."""
    import re
    from pathlib import Path

    m = re.search(
        r"type(\d)_burst_start_time_(\d{4}-\d{2}-\d{2})_(\d{2}:\d{2}:\d{2})_(\d{4})_stop_time_(\d{4})"
        r"_subband_min_(\d+)_subband_max_(\d+)",
        Path(path).name,
    )
    if not m:
        raise ValueError(f"unrecognised LV614 file name: {Path(path).name}")
    return {
        "type": int(m.group(1)),
        "date": m.group(2),
        "start": m.group(3),
        "hhmm_start": m.group(4),
        "hhmm_stop": m.group(5),
        "freqs": lofar_freqs(int(m.group(6)), int(m.group(7))),
    }
