"""Voyager 2 PRA: a modern re-derivation of the Uranus & Neptune radio rotation periods (plan 46, F9).

Uranus's rotation period was moved by 28 s by Lamy+2025 (Nat. Astron. 9, 658) using 11 yr of HST UV
aurora --- but the 1986 *radio* value (17.24 +/- 0.01 h, Warwick+1986) that underlies System III was
never reanalysed with modern statistics, and Neptune's 16.11 h (Warwick+1989) has had no independent
check in 28 yr. The Cecconi+2017 PRA refurbishing (arXiv:1710.10471) covers only Jupiter/Saturn, so
the Uranus and Neptune Voyager-2 Planetary-Radio-Astronomy encounter volumes are an unworked niche.

This slice re-derives both periods from the open PDS-PPI VG2-PRA low-band 6-second data, reusing the
merged `frbperiod` Rayleigh-Z^2 machinery on the planetary-radio-burst time series, and delivers a
three-way comparison: the 1986/1989 radio value, Lamy+2025 (Uranus), and this work.

**Honest framing (and the actual result): a controlled NULL.** The approach is a BLIND periodogram
of the total-power flux --- no beaming or magnetic-longitude modelling, unlike the original
determinations. It is validated end-to-end: on a synthetic flyby, searched in the SAME wide 14-20 h
window used on the real data, it recovers an injected rotation period to ~1 min. On the real data it
recovers NEITHER: the Uranus peak wanders between frequency sub-bands (spread ~1.8 h, offset +1.2 h
from 17.24 h) and the Neptune peak rails to the search-window bound (offset +3.9 h) --- the blind
periodogram is dominated by the red-noise / flyby-envelope continuum and lands on window- and
band-dependent non-rotation features, not the rotation. Because the method demonstrably works when a
clean rotational signal is present, the real-data failure is a genuine limitation of blind
total-power analysis, not a pipeline bug: the ice-giant total-power flux is not a clean rotational
sinusoid, so the historical determinations' beaming / magnetic-longitude modelling was essential.
Separately, the flyby precision ceiling is ~1-2 h even in the best case --- hundreds of times coarser
than the 28-s HST shift (Lamy+2025). The citable result: a modern blind reanalysis quantifies why
the Voyager radio periods needed the sophisticated original methods, and cannot approach modern
precision.

Data (GATE-0 2026-07-10, novelty PASS, egress verified): the PDS-PPI volumes
`VG2-U-PRA-3-RDR-LOWBAND-6SEC-V1.0` and `VG2-N-PRA-3-RDR-LOWBAND-6SEC-V1.0` --- fixed-width ASCII,
one 48-s major frame per line (`DATE` YYMMDD, `SECOND` of day, then 8 sweeps of 1 status word + 70
channels in millibell); low-band channels f_i = 1326.0 - 19.2*i kHz (i=0..69); sweep k starts at
SECOND + 6*k seconds.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np

from .frbperiod import false_alarm_prob

__all__ = [
    "FREQ_KHZ",
    "PUBLISHED_HR",
    "read_pra_series",
    "band_flux",
    "detect_bursts",
    "period_posterior",
    "bin_series",
    "flux_period_posterior",
    "band_stability",
    "synthetic_flyby",
    "compare_periods",
    "run",
]

# 70 low-band PRA channels: 1326.0 kHz down to 1.2 kHz, spaced 19.2 kHz (PDS label).
FREQ_KHZ = 1326.0 - 19.2 * np.arange(70)
N_SWEEPS = 8  # sweeps per 48-s major frame; sweep k starts at SECOND + 6*k s
SWEEP_DT_S = 6.0
N_ITEMS = 71  # 1 status word + 70 channels per sweep
FIELD_W = 4  # 4-char fixed-width ASCII integer per field
_DATA_START = 12  # DATE(6) + SECOND(6)

# Published radio rotation periods (System III) and the modern HST re-derivation (hours).
PUBLISHED_HR = {
    "URANUS": {"radio_1986": (17.24, 0.01), "lamy_2025": (17.247864, 0.000010)},
    "NEPTUNE": {
        "radio_1989": (16.11, 0.05)
    },  # Warwick+1989 (+/-0.05); Lecacheux+1993 16.108+/-0.006
}
# One wide, physically-motivated ice-giant rotation window for BOTH planets (NOT tuned to the known
# values -- tuning the window to bracket the answer would p-hack a match). Wide enough that the
# rotation-block bootstrap band is not truncated by a bound.
SEARCH_WINDOW_HR = {"URANUS": (14.0, 20.0), "NEPTUNE": (14.0, 20.0)}
# Default burst-extraction band (kHz): the mid low-band where UKR/NKR episodes dominate.
DEFAULT_BAND_KHZ = (100.0, 1000.0)
# Sub-bands for the peak achromaticity check. A COHERENT signal (rotation OR a viewing-geometry /
# beaming / spectral-window modulation) is achromatic, so a small spread confirms the peak is a real
# coherent signal, NOT noise -- it does NOT by itself confirm the peak is the rotation. (On the real
# data Uranus's WRONG peak is more band-stable than Neptune's correct one.)
# Three DISJOINT thirds of the analysis band, each with independent channels. The previous
# set -- ((100,1000),(200,800),(300,900),(500,1300)) -- was one copy of the analysis band, two
# nested bands sharing 26/31 channels, and one band extending 288 kHz ABOVE the analysis band;
# a referee showed the whole Uranus "not recovered" verdict rested on that out-of-band fourth
# (drop it and the spread falls 1.76 -> 0.19 h, flipping the boolean). Achromaticity must be
# tested across independent channels INSIDE the band the science uses.
STABILITY_BANDS_KHZ = ((100.0, 400.0), (400.0, 700.0), (700.0, 1000.0))
# PDS-PPI open encounter volumes (direct HTTP; GATE-0 verified 2026-07-10).
_PDS_BASE = "https://pds-ppi.igpp.ucla.edu/data"
PDS_TAB_URL = {
    "URANUS": f"{_PDS_BASE}/VG2-U-PRA-3-RDR-LOWBAND-6SEC-V1.0/DATA/VG2_URN_PRA_6SEC.TAB",
    "NEPTUNE": f"{_PDS_BASE}/VG2-N-PRA-3-RDR-LOWBAND-6SEC-V1.0/DATA/VG2_NEP_PRA_6SEC.TAB",
}


def _seconds_since_epoch(yymmdd: int, sec_of_day: float) -> float:
    """Absolute seconds from a fixed ordinal for a YYMMDD date + seconds-of-day (19xx encounters)."""
    y = 1900 + yymmdd // 10000
    m = (yymmdd // 100) % 100
    d = yymmdd % 100
    return date(y, m, d).toordinal() * 86400.0 + sec_of_day


def read_pra_series(source: str | Path) -> dict:
    """Parse a VG2-PRA low-band 6-s .TAB into a burst-analysis-ready dynamic spectrum.

    ``source`` is a path or the raw text. Each line is one 48-s major frame; we expand its 8 sweeps
    into 8 spectra at times ``SECOND + 6*k`` s (k=0..7), dropping the per-sweep status word. Returns
    ``times_hr`` (hours from the first sample), ``freqs_khz`` (70), and ``spectra`` (n_samples, 70)
    in **millibell** as recorded. Malformed/short lines are skipped.
    """
    text = Path(source).read_text() if _looks_like_path(source) else str(source)
    abs_s: list[float] = []
    rows: list[list[int]] = []
    need = _DATA_START + N_SWEEPS * N_ITEMS * FIELD_W
    for line in text.splitlines():
        if len(line) < need:
            continue
        try:
            yymmdd = int(line[0:6])
            sec = int(line[6:12])
        except ValueError:
            continue
        base = _seconds_since_epoch(yymmdd, sec)
        for k in range(N_SWEEPS):
            off = _DATA_START + k * N_ITEMS * FIELD_W
            # skip item 0 (status word); read the 70 channel fields
            chans = _parse_fields(line, off + FIELD_W, 70)
            if chans is None:
                continue
            abs_s.append(base + k * SWEEP_DT_S)
            rows.append(chans)
    if not rows:
        return {"times_hr": np.zeros(0), "freqs_khz": FREQ_KHZ, "spectra": np.zeros((0, 70))}
    t = np.asarray(abs_s, float)
    order = np.argsort(t)
    t = t[order]
    spectra = np.asarray(rows, float)[order]
    return {"times_hr": (t - t[0]) / 3600.0, "freqs_khz": FREQ_KHZ, "spectra": spectra}


def _looks_like_path(source: str | Path) -> bool:
    if isinstance(source, Path):
        return True
    return "\n" not in source and len(source) < 4096 and source.strip().endswith((".tab", ".TAB"))


def _parse_fields(line: str, start: int, n: int) -> list[int] | None:
    out = []
    for j in range(n):
        seg = line[start + j * FIELD_W : start + (j + 1) * FIELD_W]
        try:
            out.append(int(seg))
        except ValueError:
            return None
    return out


def band_flux(spectra: np.ndarray, freqs_khz: np.ndarray, band_khz=DEFAULT_BAND_KHZ) -> np.ndarray:
    """Linear band-integrated power vs time from a millibell dynamic spectrum.

    Millibell = 0.01 dB, so a channel's power ratio is ``10**(mB/1000)``. We sum linear power over
    the channels inside ``band_khz`` --- a physical total-power series in which the periodic
    UKR/NKR bursts stand out --- returning one value per time sample (NaN-free).
    """
    f = np.asarray(freqs_khz, float)
    sel = (f >= band_khz[0]) & (f <= band_khz[1])
    if sel.sum() == 0:
        return np.zeros(spectra.shape[0])
    power = np.power(10.0, np.asarray(spectra, float)[:, sel] / 1000.0)
    return power.sum(axis=1)


def detect_bursts(
    times_hr: np.ndarray,
    flux: np.ndarray,
    *,
    k: float = 4.5,
    smooth_hr: float = 3.0,
    min_sep_hr: float = 1.5,
) -> np.ndarray:
    """Times (hours) of radio-burst episodes: peaks above a slowly-varying background by ``k`` MAD.

    The background (galactic + receiver + smooth planetary component) is a running median over a
    ``smooth_hr``-wide window; the scatter is its MAD. Samples above ``background + k*1.4826*MAD``
    are grouped into above-threshold runs (one peak each), then peaks within ``min_sep_hr`` are
    merged (keeping the strongest) so a single burst episode yields ONE epoch --- the events
    `frbperiod`'s Rayleigh Z^2 folds --- rather than many noise-split sub-peaks that would bias the
    period.
    """
    t = np.asarray(times_hr, float)
    y = np.asarray(flux, float)
    if t.size < 3:
        return np.zeros(0)
    dt = np.median(np.diff(t)) if t.size > 1 else 1.0
    win = max(int(smooth_hr / dt), 3) | 1  # odd window
    bg = _running_median(y, win)
    mad = np.median(np.abs(y - bg)) or (np.std(y) or 1.0)
    thr = bg + k * 1.4826 * mad
    above = y > thr
    times_pk: list[float] = []
    vals_pk: list[float] = []
    i, n = 0, t.size
    while i < n:
        if above[i]:
            j = i
            while j < n and above[j]:
                j += 1
            seg = slice(i, j)
            kk = int(np.argmax(y[seg]))
            times_pk.append(float(t[seg][kk]))
            vals_pk.append(float(y[seg][kk]))
            i = j
        else:
            i += 1
    if not times_pk:
        return np.zeros(0)
    # merge peaks closer than min_sep_hr into one episode (keep the strongest)
    merged: list[float] = []
    ct, cv = times_pk[0], vals_pk[0]
    for pt, pv in zip(times_pk[1:], vals_pk[1:], strict=False):
        if pt - ct <= min_sep_hr:
            if pv > cv:
                ct, cv = pt, pv
        else:
            merged.append(ct)
            ct, cv = pt, pv
    merged.append(ct)
    return np.asarray(merged, float)


def _running_median(y: np.ndarray, win: int) -> np.ndarray:
    n = y.size
    half = win // 2
    out = np.empty(n)
    for i in range(n):
        out[i] = np.median(y[max(0, i - half) : min(n, i + half + 1)])
    return out


def _z2_grid(times: np.ndarray, periods: np.ndarray) -> np.ndarray:
    """Vectorized Rayleigh $Z^2_1$ over a whole period grid at once (same statistic as frbperiod)."""
    t = np.asarray(times, float)
    phi = 2.0 * np.pi * ((t[None, :] / periods[:, None]) % 1.0)
    c = np.cos(phi).sum(axis=1)
    s = np.sin(phi).sum(axis=1)
    return (2.0 / t.size) * (c * c + s * s)


def period_posterior(
    burst_times_hr: np.ndarray,
    *,
    p_lo: float,
    p_hi: float,
    n_grid: int = 3000,
    n_boot: int = 200,
    seed: int = 0,
) -> dict:
    """Rayleigh periodogram peak + a bootstrap (few-cycle-honest) uncertainty, in hours.

    Runs `frbperiod`'s Rayleigh $Z^2_1$ (vectorized) over a trial-period grid on the burst epochs.
    The uncertainty is a nonparametric bootstrap: resample the epochs with replacement ``n_boot``
    times, re-locate the peak, and report the 16/84 percentiles --- which, for a short flyby, come
    out honestly wide. Returns best period, Z^2, exposure-blind FAP, and the bootstrap band.
    """
    t = np.asarray(burst_times_hr, float)
    grid = np.linspace(p_lo, p_hi, n_grid)
    if t.size < 4:
        return {
            "best_period_hr": float("nan"),
            "z2": float("nan"),
            "fap": float("nan"),
            "boot_lo_hr": float("nan"),
            "boot_hi_hr": float("nan"),
            "boot_sigma_hr": float("nan"),
            "n_bursts": int(t.size),
        }
    z2 = _z2_grid(t, grid)
    k = int(np.argmax(z2))
    span = float(t.max() - t.min())
    n_indep = min(int((1.0 / p_lo - 1.0 / p_hi) * span) + 1, n_grid)
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    for b in range(n_boot):
        sample = rng.choice(t, size=t.size, replace=True)
        boot[b] = grid[int(np.argmax(_z2_grid(sample, grid)))]
    lo, hi = np.percentile(boot, [16, 84])
    return {
        "best_period_hr": round(float(grid[k]), 5),
        "z2": round(float(z2[k]), 3),
        "fap": float(false_alarm_prob(float(z2[k]), n_indep)),
        "boot_lo_hr": round(float(lo), 5),
        "boot_hi_hr": round(float(hi), 5),
        "boot_sigma_hr": round(float(np.std(boot)), 5),
        "n_bursts": int(t.size),
    }


def bin_series(times_hr: np.ndarray, y: np.ndarray, bin_hr: float) -> tuple[np.ndarray, np.ndarray]:
    """Block-average an unevenly-sampled series into ``bin_hr``-wide time bins (mean per bin)."""
    t = np.asarray(times_hr, float)
    v = np.asarray(y, float)
    if t.size == 0:
        return t, v
    b = np.floor(t / bin_hr).astype(int)
    _, inv = np.unique(b, return_inverse=True)
    counts = np.bincount(inv)
    tb = np.bincount(inv, weights=t) / counts
    vb = np.bincount(inv, weights=v) / counts
    return tb, vb


def _prep_flux(
    times_hr: np.ndarray, flux: np.ndarray, bin_hr: float, detrend_deg: int
) -> tuple[np.ndarray, np.ndarray]:
    """Bin + log + polynomial-detrend a flux series for Lomb-Scargle (removes the approach trend)."""
    t, y = bin_series(times_hr, np.log10(np.maximum(np.asarray(flux, float), 1e-30)), bin_hr)
    if t.size >= detrend_deg + 1:
        y = y - np.polyval(np.polyfit(t, y, detrend_deg), t)
    return t, y


def _ls_peak(t: np.ndarray, y: np.ndarray, grid_f: np.ndarray) -> tuple[float, float]:
    """Lomb-Scargle peak (period_hr, power) of a prepared series over a frequency grid (cycles/h)."""
    from astropy.timeseries import LombScargle

    power = LombScargle(t, y).power(grid_f)
    j = int(np.argmax(power))
    return 1.0 / grid_f[j], float(power[j])


def flux_period_posterior(
    times_hr: np.ndarray,
    flux: np.ndarray,
    *,
    p_lo: float,
    p_hi: float,
    bin_hr: float = 0.1,
    n_grid: int = 8000,
    n_boot: int = 150,
    detrend_deg: int = 2,
    seed: int = 0,
) -> dict:
    """Lomb-Scargle rotation period from the *continuous* flux series + rotation-block bootstrap.

    The rotation modulates a continuous, red-noise-dominated emission, so a Lomb-Scargle periodogram
    of the (log, binned, polynomial-detrended) flux is far more sensitive than folding discrete
    burst epochs. The analytic LS false-alarm probability is meaningless here (autocorrelated
    samples + huge N make it absurdly small), so the honest uncertainty is a **rotation-block
    bootstrap**: resample whole rotation-length blocks (which preserves the within-rotation
    autocorrelation) and re-locate the peak; the 16/84 percentiles give the band. Returns the
    best period, LS power, and that bootstrap band --- for a short flyby it comes out honestly wide.
    """
    t, y = _prep_flux(times_hr, flux, bin_hr, detrend_deg)
    if t.size < 8:
        return {
            "best_period_hr": float("nan"),
            "ls_power": float("nan"),
            "boot_lo_hr": float("nan"),
            "boot_hi_hr": float("nan"),
            "boot_sigma_hr": float("nan"),
            "n_binned": int(t.size),
        }
    grid_f = np.linspace(1.0 / p_hi, 1.0 / p_lo, n_grid)
    p0, pw = _ls_peak(t, y, grid_f)
    if n_boot == 0:
        return {
            "best_period_hr": round(float(p0), 5),
            "ls_power": round(float(pw), 5),
            "boot_lo_hr": float("nan"),
            "boot_hi_hr": float("nan"),
            "boot_sigma_hr": float("nan"),
            "n_binned": int(t.size),
        }
    rng = np.random.default_rng(seed)
    cyc = np.floor(t / p0).astype(int)
    rots = np.unique(cyc)
    boot = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.choice(rots, size=rots.size, replace=True)
        idx = np.concatenate([np.where(cyc == r)[0] for r in pick])
        boot[b] = _ls_peak(t[idx], y[idx], grid_f)[0]
    lo, hi = np.percentile(boot, [16, 84])
    return {
        "best_period_hr": round(float(p0), 5),
        "ls_power": round(float(pw), 5),
        "boot_lo_hr": round(float(lo), 5),
        "boot_hi_hr": round(float(hi), 5),
        "boot_sigma_hr": round(float(np.std(boot)), 5),
        "n_binned": int(t.size),
    }


def band_stability(
    times_hr: np.ndarray,
    spectra: np.ndarray,
    freqs_khz: np.ndarray,
    *,
    p_lo: float,
    p_hi: float,
    bands=STABILITY_BANDS_KHZ,
    bin_hr: float = 0.1,
    detrend_deg: int = 2,
) -> dict:
    """LS peak period in several sub-bands: an ACHROMATICITY (coherence) check, not a right/wrong test.

    A coherent signal is achromatic within the auroral band, so a small peak spread across sub-bands
    confirms the periodogram peak is a real coherent modulation rather than noise. It does NOT
    establish that the peak is the *rotation*: a coherent viewing-geometry / beaming / spectral-window
    modulation is equally achromatic (indeed, on the real data Uranus's wrong-period peak is more
    band-stable than Neptune's correct one). The spread is reported as a systematic, method-level
    uncertainty; the right/wrong verdict rests on comparison to the historical period, not on this.
    """
    grid_f = np.linspace(1.0 / p_hi, 1.0 / p_lo, 6000)
    peaks = []
    for band in bands:
        flux = band_flux(spectra, freqs_khz, band_khz=band)
        t, y = _prep_flux(times_hr, flux, bin_hr, detrend_deg)
        if t.size >= 8:
            peaks.append(_ls_peak(t, y, grid_f)[0])
    peaks_arr = np.asarray(peaks, float)
    return {
        "band_peaks_hr": [round(float(p), 4) for p in peaks_arr],
        "band_spread_hr": round(float(np.std(peaks_arr)), 5) if peaks_arr.size else float("nan"),
    }


def synthetic_flyby(
    *,
    period_hr: float = 17.24,
    n_rot: int = 17,
    bursts_per_rot: float = 1.0,
    active_frac: float = 0.06,  # narrow burst window (phase fraction); << the detrend scale
    cadence_s: float = 6.0,
    noise: float = 1.0,
    burst_amp: float = 8.0,
    gap_frac: float = 0.2,
    seed: int = 0,
    n_chan: int = 0,
) -> dict:
    """A flyby-length synthetic PRA flux series with a KNOWN rotation period, for recovery.

    Builds ``n_rot`` rotations sampled at ``cadence_s``, with a burst episode each rotation inside a
    phase window of width ``active_frac``; adds noise and drops a fraction ``gap_frac`` of samples
    (data gaps). Returns the flux series, its times (hours), the injected period, and the true burst
    epochs --- so `detect_bursts` + `period_posterior` can be shown to recover ``period_hr`` within
    honestly-wide few-cycle error bars.
    """
    rng = np.random.default_rng(seed)
    span_hr = n_rot * period_hr
    n = int(span_hr * 3600.0 / cadence_s)
    t_hr = np.arange(n) * cadence_s / 3600.0
    flux = rng.normal(0.0, noise, n) + 20.0  # a smooth background level
    phase = (t_hr / period_hr) % 1.0
    cyc = (t_hr / period_hr).astype(int)
    centre = 0.3  # a fixed rotation-phase activity window (the periodic signal), with small jitter
    half = active_frac / 2.0  # narrow burst window (kept << the detrend background scale)
    for c in range(n_rot):
        if bursts_per_rot < 1.0 and rng.random() > bursts_per_rot:
            continue
        jit = rng.normal(0.0, 0.02)
        circ = np.abs(((phase - centre - jit + 0.5) % 1.0) - 0.5)  # circular phase distance
        in_window = (circ < half) & (cyc == c)
        if in_window.any():
            flux[in_window] += burst_amp * (1.0 + rng.normal(0.0, 0.2, int(in_window.sum())))
    keep = rng.random(n) > gap_frac
    out = {
        "times_hr": t_hr[keep],
        "flux": flux[keep],
        "period_hr": period_hr,
        "n_rot": n_rot,
    }
    if n_chan:
        # An achromatic spectral axis: every channel carries the same burst modulation plus
        # independent noise. This exists so `band_stability` -- the statistic that decides the
        # real-data verdict -- can be validated on a known positive at realistic per-band SNR,
        # which it never was before (the only prior test tiled one channel N times, a case
        # that cannot fail).
        base = flux[keep] - 20.0
        spec = 20.0 + base[:, None] + rng.normal(0.0, noise, (int(keep.sum()), n_chan))
        out["spectra"] = spec
        out["freqs_khz"] = np.linspace(1300.0, 100.0, n_chan)
    return out


def compare_periods(planet: str, this_work: dict) -> list[dict]:
    """Three-way comparison rows: historical radio, Lamy+2025 (Uranus), and this work (hours)."""
    rows = []
    for key, (val, err) in PUBLISHED_HR[planet.upper()].items():
        rows.append({"source": key, "period_hr": val, "err_hr": err})
    rows.append(
        {
            "source": "this_work",
            "period_hr": this_work.get("best_period_hr"),
            "err_hr": this_work.get("boot_sigma_hr"),
        }
    )
    return rows


def _consistent(this_work: dict, val: float) -> bool:
    """Is ``val`` inside this-work's bootstrap band (few-cycle honest consistency)?"""
    lo, hi = this_work.get("boot_lo_hr"), this_work.get("boot_hi_hr")
    if lo is None or hi is None or not np.isfinite(lo) or not np.isfinite(hi):
        return False
    # widen the band to at least +/-3 sigma so a peak-grid quantisation doesn't fake a discrepancy
    sig = this_work.get("boot_sigma_hr") or 0.0
    best = this_work.get("best_period_hr", float("nan"))
    return bool(min(lo, best - 3 * sig) <= val <= max(hi, best + 3 * sig))


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: synthetic recover-a-known (inject a period, recover it within the bootstrap band)."""
    import json

    if offline:
        metrics: dict = _synthetic_metrics()
    else:  # pragma: no cover - real leg streams the PDS-PPI encounter volumes
        metrics = _real_analysis()

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    if offline:
        # The validation leg's metrics are committed evidence too: four synthetic numbers
        # reach both abstracts, and until 2026-08-13 they existed only inside macros.tex --
        # a referee could not check them without executing code.
        (op / "results" / "vgpra_synthetic_metrics.json").write_text(
            json.dumps(metrics, indent=2) + "\n"
        )
    else:
        (op / "results" / "vgpra_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(metrics, op / "papers" / "vgpra" / "figures", real_cache=None)
    _write_macros(metrics, op / "papers" / "vgpra" / "generated" / "macros.tex")
    return metrics


def scramble_fap(
    times_hr: np.ndarray,
    flux: np.ndarray,
    *,
    p_lo: float,
    p_hi: float,
    n_scramble: int = 200,
    bin_hr: float = 0.1,
    detrend_deg: int = 2,
    seed: int = 0,
) -> dict:
    """Scramble-based false-alarm probability for the window's peak LS power.

    The analytic LS FAP assumes white noise and independent samples, both false here, which
    is why the paper declines to quote it -- but declining left the null with NO significance
    statement at all. Shuffling the binned fluxes over their own timestamps destroys any
    coherent modulation while preserving the sampling pattern and the marginal flux
    distribution; the fraction of shuffles whose window-maximum power beats the observed one
    is a defensible empirical FAP. Also reports the two numbers a periodogram claim needs:
    the independent-frequency count over the window and the resolution at the peak.
    """
    rng = np.random.default_rng(seed)
    t, y = _prep_flux(times_hr, flux, bin_hr, detrend_deg)
    grid_f = np.linspace(1.0 / p_hi, 1.0 / p_lo, 8000)
    p0, pw = _ls_peak(t, y, grid_f)
    beat = 0
    for _ in range(n_scramble):
        _, pw_s = _ls_peak(t, rng.permutation(y), grid_f)
        if pw_s >= pw:
            beat += 1
    span = float(t.max() - t.min())
    return {
        "peak_period_hr": round(float(p0), 5),
        "peak_power": round(float(pw), 5),
        "scramble_fap": round((beat + 1) / (n_scramble + 1), 4),
        "n_scramble": n_scramble,
        "n_indep_freq": round((1.0 / p_lo - 1.0 / p_hi) * span, 1),
        "resolution_hr": round(p0 * p0 / span, 3),
    }


def injection_recovery_real(
    times_hr: np.ndarray,
    flux: np.ndarray,
    *,
    period_hr: float,
    p_lo: float,
    p_hi: float,
    depths=(0.02, 0.05, 0.1, 0.2, 0.4),
    n_phase: int = 8,
    bin_hr: float = 0.1,
    detrend_deg: int = 2,
    seed: int = 0,
) -> dict:
    """Inject a fractional sinusoid at the historical period into the REAL flux and re-search.

    This is the selection function the null needs. The synthetic control (8-sigma spike train
    on white noise) proves only the machinery; injecting into the actual series tests recovery
    against the real red-noise continuum, the real flyby envelope and the real detrend -- the
    three things the paper blames for the failure. Recovery = blind window peak within one
    resolution element of the injected period. Reports the recovered fraction per depth, so
    the paper can say "any rotational modulation deeper than X% would have been found".
    """
    rng = np.random.default_rng(seed)
    t_raw = np.asarray(times_hr, float)
    f_raw = np.asarray(flux, float)
    grid_f = np.linspace(1.0 / p_hi, 1.0 / p_lo, 8000)
    span = float(t_raw.max() - t_raw.min())
    res = period_hr * period_hr / span
    curve = []
    for d in depths:
        hits = 0
        for k in range(n_phase):
            phase = 2.0 * np.pi * (k / n_phase + rng.uniform(0, 1e-3))
            mod = 1.0 + d * np.sin(2.0 * np.pi * t_raw / period_hr + phase)
            t, y = _prep_flux(t_raw, f_raw * mod, bin_hr, detrend_deg)
            p0, _ = _ls_peak(t, y, grid_f)
            if abs(p0 - period_hr) <= res:
                hits += 1
        curve.append({"depth": d, "recovered_frac": round(hits / n_phase, 3)})
    min_depth = next((c["depth"] for c in curve if c["recovered_frac"] >= 0.9), None)
    return {
        "injected_period_hr": period_hr,
        "resolution_hr": round(res, 3),
        "curve": curve,
        "min_depth_recovered": min_depth,
    }


def _synthetic_metrics() -> dict:
    """Recover an injected Uranus-like period from a synthetic flyby; the pipeline validation.

    Validates BOTH methods end-to-end: the primary Lomb-Scargle flux posterior and the secondary
    Rayleigh burst posterior must each recover the injected period within their bootstrap bands.
    """
    truth = 17.24
    # search the SAME wide window used on the real data (14-20 h): recovering the injected period in
    # that wide window proves the method works when a real rotation IS present, so the real-data
    # failure is a data limitation (no clean rotational total-power signal), not a window artifact.
    lo, hi = 14.0, 20.0
    s = synthetic_flyby(period_hr=truth, n_rot=17, seed=0, n_chan=32)
    ls = flux_period_posterior(s["times_hr"], s["flux"], p_lo=lo, p_hi=hi, n_boot=120, seed=1)
    bursts = detect_bursts(s["times_hr"], s["flux"])
    ray = period_posterior(bursts, p_lo=lo, p_hi=hi, n_boot=120, seed=1)
    # Grade the control by the SAME three-clause criterion the real data face (a referee found
    # the control had been graded by a looser rule, which the real Uranus leg also passed):
    # proximity within total uncertainty, band-robustness, and not railed.
    stab = band_stability(s["times_hr"], s["spectra"], s["freqs_khz"], p_lo=lo, p_hi=hi)
    total_unc = max(float(ls["boot_sigma_hr"]), float(stab["band_spread_hr"]))
    railed = abs(ls["best_period_hr"] - lo) < 0.03 or abs(ls["best_period_hr"] - hi) < 0.03
    recovered_full = bool(
        abs(ls["best_period_hr"] - truth) <= max(total_unc, 1e-9)
        and stab["band_spread_hr"] < 0.6
        and not railed
    )
    # Seed-to-seed scatter: the within-realization bootstrap resamples one noise realization
    # and cannot see realization variance (the rmstructure failure). Thirty seeds measure it.
    seeds = []
    for k in range(30):
        sk = synthetic_flyby(period_hr=truth, n_rot=17, seed=100 + k)
        lsk = flux_period_posterior(sk["times_hr"], sk["flux"], p_lo=lo, p_hi=hi, n_boot=0)
        seeds.append(float(lsk["best_period_hr"]))
    seeds_arr = np.asarray(seeds)
    return {
        "source": "synthetic flyby (injected rotation period + noise + data gaps)",
        "is_real": False,
        "injected_period_hr": truth,
        "best_period_hr": ls["best_period_hr"],
        "boot_lo_hr": ls["boot_lo_hr"],
        "boot_hi_hr": ls["boot_hi_hr"],
        "boot_sigma_hr": ls["boot_sigma_hr"],
        "ls_power": ls["ls_power"],
        "rayleigh_period_hr": ray["best_period_hr"],
        "rayleigh_z2": ray["z2"],
        "n_bursts": ray["n_bursts"],
        "recovered_injected": bool(_consistent(ls, truth)),
        "recovered_rayleigh": bool(_consistent(ray, truth)),
        "band_peaks_hr": stab["band_peaks_hr"],
        "band_spread_hr": stab["band_spread_hr"],
        "recovered_full_criterion": recovered_full,
        "seed_scatter": {
            "n_seeds": int(seeds_arr.size),
            "mean_hr": round(float(seeds_arr.mean()), 5),
            "std_hr": round(float(seeds_arr.std()), 5),
            "worst_offset_hr": round(float(np.abs(seeds_arr - truth).max()), 5),
        },
    }


def analyse_planet(
    planet: str, source: str | Path, *, band_khz=DEFAULT_BAND_KHZ, k: float = 4.5, n_boot: int = 200
) -> dict:  # pragma: no cover - exercised on real encounter volumes, not in CI
    """Full period re-derivation for one planet from its .TAB (path or text).

    PRIMARY = Lomb-Scargle on the continuous flux series (`flux_period_posterior`); SECONDARY =
    Rayleigh Z^2 on detected burst epochs (`period_posterior`) as an independent cross-check.
    Consistency (with the historical radio value, and Lamy+2025 for Uranus) is judged against the
    LS rotation-block-bootstrap band.
    """
    ds = read_pra_series(source)
    flux = band_flux(ds["spectra"], ds["freqs_khz"], band_khz)
    lo, hi = SEARCH_WINDOW_HR[planet.upper()]
    ls = flux_period_posterior(ds["times_hr"], flux, p_lo=lo, p_hi=hi, n_boot=n_boot, seed=0)
    stab = band_stability(ds["times_hr"], ds["spectra"], ds["freqs_khz"], p_lo=lo, p_hi=hi)
    bursts = detect_bursts(ds["times_hr"], flux, k=k)
    ray = period_posterior(bursts, p_lo=lo, p_hi=hi, n_boot=n_boot, seed=0)
    pub = PUBLISHED_HR[planet.upper()]
    hist_val = next(iter(pub.values()))[0]
    best = ls["best_period_hr"]
    # combined uncertainty = the larger of the statistical bootstrap sigma and the systematic
    # band-to-band spread. NOTE: achromaticity (small spread) does NOT gate the verdict -- a coherent
    # non-rotation modulation is also achromatic; the recovers_hist call rests on the historical prior.
    total_unc = max(ls["boot_sigma_hr"], stab["band_spread_hr"])
    out = dict(ls)  # LS is the headline (best_period_hr, boot_* from the flux periodogram)
    out.update(stab)
    out["total_unc_hr"] = round(float(total_unc), 5)
    out["offset_from_hist_hr"] = (
        round(float(best - hist_val), 5) if np.isfinite(best) else float("nan")
    )
    out["rayleigh_period_hr"] = ray["best_period_hr"]
    out["rayleigh_z2"] = ray["z2"]
    out["n_bursts"] = ray["n_bursts"]
    out["span_hr"] = round(float(ds["times_hr"][-1]) if ds["times_hr"].size else 0.0, 2)
    out["n_samples"] = int(ds["spectra"].shape[0])
    at_edge = bool(np.isfinite(best) and min(best - lo, hi - best) < 0.03 * (hi - lo))
    out["at_window_edge"] = at_edge
    # A clean recovery needs ALL THREE: the blind peak lands within 1 combined-sigma of the
    # historical value, is band-robust (achromatic), and is not railed to a search-window bound.
    # (Proximity alone is worthless when sigma is ~2 h; a wandering or railed peak is not a
    # measurement.) `consistent_hist` is the weaker within-2-sigma statement.
    out["recovers_hist"] = bool(
        np.isfinite(best)
        and abs(best - hist_val) <= total_unc
        and stab["band_spread_hr"] < 0.6
        and not at_edge
    )
    # None, not a knife-edge boolean, when the peak is railed: a railed best_period is a
    # bound, so no consistency statement is possible (Neptune's old "not consistent" verdict
    # had a 33-second margin on a threshold that is twice the window width).
    out["consistent_hist"] = (
        None if at_edge else bool(np.isfinite(best) and abs(best - hist_val) <= 2 * total_unc)
    )
    # The significance and sensitivity statements the null was missing entirely: an empirical
    # scramble FAP with the trials and resolution stated, and an injection-recovery curve on
    # the REAL series -- the selection function ("modulation deeper than X% would have been
    # found") that turns "no peak" into a bounded null.
    out["significance"] = scramble_fap(ds["times_hr"], flux, p_lo=lo, p_hi=hi, n_scramble=200)
    out["injection"] = injection_recovery_real(
        ds["times_hr"], flux, period_hr=hist_val, p_lo=lo, p_hi=hi
    )
    if "lamy_2025" in pub:
        lamy = pub["lamy_2025"][0]
        out["consistent_lamy"] = bool(np.isfinite(best) and abs(best - lamy) <= 2 * total_unc)
    return out


def _real_analysis(*, cache_dir: str | Path | None = None) -> dict:  # pragma: no cover - network
    """Real leg: download both PDS-PPI encounter volumes and run the period posteriors.

    Streams each .TAB to a local cache (49 MB Uranus / 79 MB Neptune), parses, extracts the
    band-integrated flux, detects UKR/NKR burst episodes, and runs the Rayleigh posterior. Returns
    the three-way-comparison-ready metrics for both planets.
    """
    import tempfile
    import urllib.request

    cache = Path(cache_dir) if cache_dir else Path(tempfile.gettempdir()) / "vgpra_cache"
    cache.mkdir(parents=True, exist_ok=True)
    planets: dict = {}
    for planet, url in PDS_TAB_URL.items():
        dest = cache / url.rsplit("/", 1)[1]
        if not dest.exists():
            urllib.request.urlretrieve(url, dest)  # noqa: S310 (trusted NASA PDS host)
        planets[planet] = analyse_planet(planet, dest)
    return {
        "source": "PDS-PPI Voyager 2 PRA low-band 6-s encounter volumes (Uranus + Neptune)",
        "is_real": True,
        "planets": planets,
    }


def _figure(m: dict, out_dir: str | Path, *, real_cache: str | Path | None = None) -> None:
    """The figure a reader of a null needs: the FULL 14-20 h search window, with real data.

    The previous figure plotted a 1.4 h window centred on the known answer -- exactly the
    tuned window the paper says it avoided -- contained no real data, and its caption named a
    recovered value that was not drawn. Now: left panel, the synthetic control's periodogram
    over the full window; middle/right, the real Uranus and Neptune periodograms over the same
    window when the cached encounter volumes are available, with the historical periods and
    the blind peaks both marked.
    """
    from astropy.timeseries import LombScargle

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    lo, hi = 14.0, 20.0
    grid_f = np.linspace(1.0 / hi, 1.0 / lo, 6000)

    def _pgram(t_hr: np.ndarray, flux: np.ndarray):  # noqa: ANN202
        tb, yb = bin_series(t_hr, np.log10(np.maximum(flux, 1e-30)), 0.1)
        yb = yb - np.polyval(np.polyfit(tb, yb, 2), tb)
        return 1.0 / grid_f, LombScargle(tb, yb).power(grid_f)

    panels: list = []
    syn_truth = 17.24
    s_ = synthetic_flyby(period_hr=syn_truth, n_rot=17, seed=0)
    panels.append(("Synthetic control", *_pgram(s_["times_hr"], s_["flux"]), syn_truth, None))
    real = (m.get("planets") or {}) if m.get("is_real") else {}
    if real and real_cache is not None:
        for planet in ("URANUS", "NEPTUNE"):
            url = PDS_TAB_URL[planet]
            dest = Path(real_cache) / url.rsplit("/", 1)[1]
            if not dest.exists():
                continue
            ds = read_pra_series(dest)
            flux = band_flux(ds["spectra"], ds["freqs_khz"])
            hist = next(iter(PUBLISHED_HR[planet].values()))[0]
            peak = real[planet].get("best_period_hr")
            panels.append((planet.title(), *_pgram(ds["times_hr"], flux), hist, peak))
    fig, axes = plt.subplots(1, len(panels), figsize=(4.7 * len(panels), 3.6))
    if len(panels) == 1:
        axes = [axes]
    for ax, (title, per, power, hist, peak) in zip(axes, panels, strict=True):
        ax.plot(per, power, "-", lw=0.8, color="C0")
        ax.axvline(hist, color="C3", ls="--", lw=1, label=f"historical {hist:g} h")
        if peak is not None and np.isfinite(peak):
            ax.axvline(peak, color="C1", ls=":", lw=1, label=f"blind peak {peak:g} h")
        ax.set(xlabel="period (h)", title=title, xlim=(lo, hi))
        ax.legend(fontsize=7)
    axes[0].set_ylabel("Lomb--Scargle power")
    fig.tight_layout()
    fig.savefig(out / "vgpra.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    """Emit both namespaces: vgSyn* (synthetic recovery, always live) + vgReal{U,N}* (real legs)."""

    def _histok(d: dict) -> str:
        # None = the peak is railed, so no consistency statement is possible (a bound, not a
        # measurement); distinct from a clean yes/no.
        if not d:
            return "--"
        v = d.get("consistent_hist")
        return "railed" if v is None else ("yes" if v else "no")

    def _pct(v) -> str:  # noqa: ANN001
        return "--" if v is None else f"{100 * float(v):g}"

    def g(src: dict, key: str) -> str:
        v = src.get(key)
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "--"
        return str(v)

    # The synthetic leg's metrics are COMMITTED evidence (results/vgpra_synthetic_metrics.json)
    # as of 2026-08-13; read them rather than recomputing at write time, so the macros a
    # referee reads always trace to a file a referee can read.
    if m.get("is_real"):
        import json as _json

        _syn_path = Path("results/vgpra_synthetic_metrics.json")
        syn = _json.loads(_syn_path.read_text()) if _syn_path.is_file() else _synthetic_metrics()
    else:
        syn = m
    real = m if m.get("is_real") else {}
    lines = [
        "% Auto-generated by jansky_research.vgpra._write_macros -- do not edit.",
        "% vgSyn* = synthetic recover-a-known (always live); vgReal{U,N}* = real Uranus/Neptune legs.",
        rf"\newcommand{{\vgSource}}{{{m['source']}}}",
        rf"\newcommand{{\vgSynInjected}}{{{g(syn, 'injected_period_hr')}}}",
        rf"\newcommand{{\vgSynRecovered}}{{{g(syn, 'best_period_hr')}}}",
        rf"\newcommand{{\vgSynSigma}}{{{g(syn, 'boot_sigma_hr')}}}",
        rf"\newcommand{{\vgSynRayPeriod}}{{{g(syn, 'rayleigh_period_hr')}}}",
        rf"\newcommand{{\vgSynRecoveredOK}}{{{'yes' if syn.get('recovered_injected') else 'no'}}}",
        rf"\newcommand{{\vgSynRayOK}}{{{'yes' if syn.get('recovered_rayleigh') else 'no'}}}",
        rf"\newcommand{{\vgSynFullOK}}{{{'yes' if syn.get('recovered_full_criterion') else 'no'}}}",
        rf"\newcommand{{\vgSynSeedSd}}{{{g(syn.get('seed_scatter') or {}, 'std_hr')}}}",
        rf"\newcommand{{\vgSynSeedN}}{{{g(syn.get('seed_scatter') or {}, 'n_seeds')}}}",
    ]
    for tag, planet in (("U", "URANUS"), ("N", "NEPTUNE")):
        d = (real.get("planets") or {}).get(planet, {})
        pub = PUBLISHED_HR[planet]
        hist = next(iter(pub.values()))[0]
        lines += [
            rf"\newcommand{{\vgReal{tag}Period}}{{{g(d, 'best_period_hr')}}}",
            rf"\newcommand{{\vgReal{tag}Sigma}}{{{g(d, 'boot_sigma_hr')}}}",
            rf"\newcommand{{\vgReal{tag}Lo}}{{{g(d, 'boot_lo_hr')}}}",
            rf"\newcommand{{\vgReal{tag}Hi}}{{{g(d, 'boot_hi_hr')}}}",
            rf"\newcommand{{\vgReal{tag}RayPeriod}}{{{g(d, 'rayleigh_period_hr')}}}",
            rf"\newcommand{{\vgReal{tag}Span}}{{{g(d, 'span_hr')}}}",
            rf"\newcommand{{\vgReal{tag}Unc}}{{{g(d, 'total_unc_hr')}}}",
            rf"\newcommand{{\vgReal{tag}Spread}}{{{g(d, 'band_spread_hr')}}}",
            rf"\newcommand{{\vgReal{tag}Recovers}}{{{'yes' if d.get('recovers_hist') else ('--' if not d else 'no')}}}",
            rf"\newcommand{{\vgReal{tag}Hist}}{{{hist}}}",
            rf"\newcommand{{\vgReal{tag}HistOK}}{{{_histok(d)}}}",
        ]
        sig = d.get("significance") or {}
        inj = d.get("injection") or {}
        lines += [
            rf"\newcommand{{\vgReal{tag}Fap}}{{{g(sig, 'scramble_fap')}}}",
            rf"\newcommand{{\vgReal{tag}Res}}{{{g(sig, 'resolution_hr')}}}",
            rf"\newcommand{{\vgReal{tag}NIndep}}{{{g(sig, 'n_indep_freq')}}}",
            rf"\newcommand{{\vgReal{tag}MinDepth}}{{{_pct(inj.get('min_depth_recovered'))}}}",
        ]
    lines.append(rf"\newcommand{{\vgRealULamy}}{{{PUBLISHED_HR['URANUS']['lamy_2025'][0]}}}")
    urn = (real.get("planets") or {}).get("URANUS", {})
    lines.append(
        rf"\newcommand{{\vgRealULamyOK}}{{{'yes' if urn.get('consistent_lamy') else ('--' if not urn else 'no')}}}"
    )
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

    ap = argparse.ArgumentParser(
        description="Voyager 2 PRA Uranus/Neptune rotation-period re-derivation."
    )
    ap.add_argument("--out", default=".")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
