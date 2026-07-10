"""Metric/decametric type II solar radio burst census on OVRO-LWA (plan 50, fable-ideas F13).

Type II bursts are slow-drifting radio emission from CME-driven coronal shocks --- the radio
signature of a shock climbing through the corona. They differ from the fast type III beam bursts
(the `solarbursts` slice) in three testable ways: a **slow** frequency drift (~0.02--1 MHz/s at
20--90 MHz, one-to-two orders slower than type III's tens--hundreds), a **minutes**-long duration
(vs type III seconds), and a **fundamental+harmonic** band pair at frequency ratio ~2 (often each
lane band-split at ~1.2 by the shock). The published OVRO-LWA real-time burst detector
(arXiv:2603.25446, ApJ 1003) is **type-III-only** (a YOLO net trained on type III); no type II
census exists on this archive (GATE-0 2026-07-09: RSTN's Cycle-24 catalogue is a different
instrument, arXiv:2512.21846 is an N=10 case study). This slice ships a slow-drift + harmonic
detector and a coverage-corrected census cross-matched to the LASCO CME catalogue.

**What ships here (honest scope).** The DETECTOR + its synthetic recover-a-known are the validated
deliverable, in core CI: (1) completeness/purity on a mixed synthetic set spanning easy
(strong+harmonic), single-lane, and near-threshold-SNR type IIs --- reported as a
completeness-vs-SNR curve, NOT a single saturated number; (2) a temporal CME cross-match whose
recovered fast-and-wide fractions ECHO the injected Gopalswamy-biased CME distribution --- i.e. a
wiring check of the match logic, not an independent reproduction (that needs real events). There
is **no real census run yet**, but the data is NOT access-blocked: the OVRO-LWA solar dynamic
spectra are on **AWS Open Data** (bucket ``ovro-lwa-solar``, path ``spec_fits/<YYYY>/<YYYYMMDD>.fits``,
directly downloadable with no login) --- the Cloudflare Turnstile only gates the query *UI*, not
the data (`s3_dspec_url`). The daily files are large (~1.7 GB, a 4D I/V dynamic spectrum), so the
real leg **streams** each day: `stream_dspec` opens the S3 FITS lazily (astropy ``use_fsspec``) and
range-reads ONLY the Stokes-I plane in time-chunks, block-averaging each to ~4 s bins before it is
kept, so a day is processed **entirely in memory with nothing written to disk** and peak memory is
one reduced chunk, never the 1.7 GB file. `sweep_day` then slides a burst-sized window across the
day and `real_census` frees each day before streaming the next (`scripts/typeii_real.py`). The
detector below is the shippable deliverable, ready to run; the coverage-corrected
occurrence-vs-cycle-phase piece is deferred with that census (it needs the full event list).

Data: OVRO-LWA solar dynamic spectra (FITS, ~15--85 MHz, ~0.26 s, ~1.7 GB/day, on AWS Open Data);
LASCO CME v2 (CDAW); GOES flares; SILSO for cycle phase (real leg). Reuse:
`solarbursts` (coronal-density drift model + `background_subtract`) drives both the synthetic type
II (a slow shock emitting fundamental+harmonic) and the type III contaminants; the real leg's
occurrence rate will reuse `ecallisto_census.coverage_corrected_rate` (deferred with the census).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .solarbursts import background_subtract

__all__ = [
    "synthetic_typeii",
    "rfi_mask",
    "track_drift_ridge",
    "harmonic_score",
    "classify_burst",
    "detect_typeii",
    "crossmatch_cme",
    "cme_association_fraction",
    "parse_cme_html",
    "parse_hek_flares",
    "crossmatch_flare",
    "occurrence_vs_phase",
    "purity_diagnostics",
    "parse_lwa_dspec",
    "s3_dspec_url",
    "sweep_day",
    "stream_dspec",
    "run",
]

# type II vs type III discriminators (metric/decametric band), pinned at GATE-0
DRIFT_SLOW_MIN = 0.01  # MHz/s: must actually drift (rejects horizontal narrowband RFI)
DRIFT_SLOW_MAX = 2.0  # MHz/s: type II upper bound; a guard band above the physical ~1, in the
# empty gap below type III's tens-hundreds MHz/s (nothing real lives at 1-2 MHz/s)
DURATION_MIN_S = 90.0  # type II lasts minutes; type III seconds
COHERENCE_MIN = 0.55  # |corr(time, ridge freq)|: a real drift is coherent, noise/RFI scatters
N_RIDGE_MIN = 20  # a genuine ridge spans many time columns; sparse noise peaks do not
HARMONIC_RATIO = 2.0  # fundamental -> second-harmonic lane ratio
FAST_CME_KMS = 900.0  # Gopalswamy+2005 fast-CME threshold for type II association
WIDE_CME_DEG = 60.0  # ... and wide-CME threshold
DATA_DIR = Path("data/typeii")


def synthetic_typeii(
    *,
    shock_speed_kms: float = 1000.0,
    r0_rsun: float = 1.4,
    with_harmonic: bool = True,
    band_split: bool = True,
    f_lo_mhz: float = 20.0,
    f_hi_mhz: float = 88.0,
    n_freq: int = 240,
    duration_s: float = 300.0,
    n_time: int = 600,
    width_mhz: float = 1.5,
    amp: float = 10.0,
    noise: float = 1.0,
    seed: int = 0,
) -> dict:
    """A synthetic OVRO-LWA-like dynamic spectrum with an injected type II shock.

    Reuses the `jansky.solar` coronal-density model that `solarbursts` uses, but driven by a slow
    CME shock (``shock_speed_kms``, ~hundreds--few-thousand km/s) instead of a fast beam --- so the
    fundamental emission ``f_p(t)`` drifts slowly. If ``with_harmonic`` a second lane is laid at
    ``HARMONIC_RATIO`` * f_p; if ``band_split`` each lane splits into two sub-bands at ~1.2. Returns
    the data (n_freq x n_time), descending freqs (MHz), times (s), and the truth drift/speed.
    """
    from jansky import solar

    rng = np.random.default_rng(seed)
    freqs = np.linspace(f_hi_mhz, f_lo_mhz, n_freq)
    times = np.linspace(0.0, duration_s, n_time)
    v_rsun_per_s = shock_speed_kms / solar.R_SUN_KM
    r_t = r0_rsun + v_rsun_per_s * times
    fp_t = solar.plasma_frequency(solar.newkirk_density(r_t, 1.0))  # fundamental, MHz vs time

    data = rng.normal(0.0, noise, (n_freq, n_time))
    lanes = [(1.0, amp)]
    if with_harmonic:
        lanes.append((HARMONIC_RATIO, 0.7 * amp))
    for ratio, lane_amp in lanes:
        centers = [ratio * fp_t]
        sub_amps = [lane_amp]
        if band_split:  # upstream/downstream sub-bands at ~1.2
            centers = [ratio * fp_t * 0.94, ratio * fp_t * 1.14]
            sub_amps = [lane_amp, 0.8 * lane_amp]
        for c_t, a in zip(centers, sub_amps, strict=True):
            for j, fc in enumerate(c_t):
                data[:, j] += a * np.exp(-0.5 * ((freqs - fc) / width_mhz) ** 2)
    # drift rate of the fundamental over the window (MHz/s, negative)
    drift = float((fp_t[-1] - fp_t[0]) / (times[-1] - times[0]))
    return {
        "data": data,
        "freqs": freqs,
        "times": times,
        "truth_drift_mhz_s": drift,
        "truth_shock_kms": shock_speed_kms,
    }


def rfi_mask(data: np.ndarray, *, rfi_factor: float = 3.0) -> np.ndarray:
    """Boolean channel mask (True = keep): drop narrowband-RFI channels by variance outlier.

    Narrowband RFI (persistent or intermittently bright in one channel) inflates that channel's
    time-domain scatter far above the typical channel. A drifting burst, by contrast, passes
    THROUGH each channel only briefly (the ridge crosses it for a fraction of the sweep), so its
    per-channel MAD stays near the quiescent level. We flag channels whose robust scatter (MAD over
    time) exceeds ``rfi_factor`` * the median channel MAD. Returns a per-channel keep mask; the
    masked fraction is the reported RFI-contamination number.
    """
    clean = background_subtract(data)
    mad_ch = 1.4826 * np.median(np.abs(clean - np.median(clean, axis=1, keepdims=True)), axis=1)
    typical = float(np.median(mad_ch[mad_ch > 0])) if np.any(mad_ch > 0) else 1.0
    return mad_ch <= rfi_factor * typical


def track_drift_ridge(
    data: np.ndarray,
    freqs: np.ndarray,
    times: np.ndarray,
    *,
    keep: np.ndarray | None = None,
    snr_threshold: float = 4.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Track the drifting band: per time-column, the peak-intensity frequency (if above noise).

    A slow-drift type II sweeps DOWN in frequency over minutes, so the peak frequency per time bin
    traces the ridge densely in time (the natural axis for a slow drift; contrast `solarbursts`
    which tracks peak-time per channel for a fast type III). RFI-masked channels (``keep`` False)
    are excluded. Returns ``(ridge_times_s, ridge_freqs_mhz)`` for columns with a real peak.
    """
    clean = background_subtract(data)
    freqs = np.asarray(freqs, float)
    times = np.asarray(times, float)
    if keep is not None:
        clean = clean[keep]
        freqs = freqs[keep]
    noise_t = 1.4826 * np.median(np.abs(clean - np.median(clean, axis=0, keepdims=True)), axis=0)
    noise_t = np.where(noise_t > 0, noise_t, np.inf)
    peak_ch = np.argmax(clean, axis=0)
    peak_val = clean[peak_ch, np.arange(clean.shape[1])]
    good = peak_val > snr_threshold * noise_t
    return times[good], freqs[peak_ch[good]]


def harmonic_score(
    data: np.ndarray,
    freqs: np.ndarray,
    ridge_times: np.ndarray,
    ridge_freqs: np.ndarray,
    times: np.ndarray,
    *,
    ratio: float = HARMONIC_RATIO,
    tol: float = 0.12,
) -> float:
    """Fraction of ridge points that have co-emission at ``ratio`` x (or 1/ratio x) their frequency.

    The type II fundamental+harmonic signature: for each ridge point (t, f) check whether the
    spectrum at the same time shows power near ``ratio*f`` (harmonic, if the ridge is the
    fundamental) OR ``f/ratio`` (fundamental, if the ridge is the harmonic), within ``tol`` in
    log-frequency and above the local column noise. Returns the fraction with a partner lane
    (0 = no harmonic structure, ~1 = clear fundamental/harmonic pair) --- a purity discriminator.
    """
    clean = background_subtract(data)
    freqs = np.asarray(freqs, float)
    times = np.asarray(times, float)
    if ridge_times.size == 0:
        return 0.0
    hits = 0
    for t, f in zip(ridge_times, ridge_freqs, strict=True):
        col = clean[:, int(np.argmin(np.abs(times - t)))]
        noise = 1.4826 * np.median(np.abs(col - np.median(col)))
        thr = 3.0 * (noise if noise > 0 else np.inf)
        for target in (ratio * f, f / ratio):
            band = np.abs(np.log(freqs / target)) < tol
            if band.any() and np.max(col[band]) > thr:
                hits += 1
                break
    return hits / ridge_times.size


def _robust_drift(ridge_times: np.ndarray, ridge_freqs: np.ndarray) -> float:
    """Robust (Theil-Sen) drift rate df/dt (MHz/s) from the ridge points; NaN if too few."""
    t = np.asarray(ridge_times, float)
    f = np.asarray(ridge_freqs, float)
    if t.size < 5:
        return float("nan")
    n = t.size
    i, j = np.triu_indices(n, k=1)
    dt = t[j] - t[i]
    ok = np.abs(dt) > 1e-6
    return float(np.median((f[j] - f[i])[ok] / dt[ok]))


def classify_burst(
    drift_mhz_s: float, duration_s: float, harm: float, coherence: float, n_ridge: int
) -> str:
    """Classify a tracked ridge: ``type_II`` (slow + long + coherent), ``type_III`` (fast), ``none``.

    Type II: a genuine SLOW downward drift (``DRIFT_SLOW_MIN`` <= |df/dt| <= ``DRIFT_SLOW_MAX``,
    negative) lasting >= ``DURATION_MIN_S``, whose ridge is COHERENT (|corr(t, f)| >=
    ``COHERENCE_MIN``) and dense (>= ``N_RIDGE_MIN`` points). The coherence + density cuts reject
    the two false positives a bare slope would admit: scattered noise/RFI peaks (low coherence) and
    a horizontal RFI line (near-zero drift). Fast drift is type III. The harmonic score is reported
    for purity but not required (single-lane type IIs exist).
    """
    if not np.isfinite(drift_mhz_s):
        return "none"
    rate = abs(drift_mhz_s)
    if rate > DRIFT_SLOW_MAX and coherence >= COHERENCE_MIN and n_ridge >= N_RIDGE_MIN:
        return "type_III"
    if (
        DRIFT_SLOW_MIN <= rate <= DRIFT_SLOW_MAX
        and drift_mhz_s < 0
        and duration_s >= DURATION_MIN_S
        and coherence >= COHERENCE_MIN
        and n_ridge >= N_RIDGE_MIN
    ):
        return "type_II"
    return "none"


def detect_typeii(
    data: np.ndarray, freqs: np.ndarray, times: np.ndarray, *, snr_threshold: float = 4.0
) -> dict:
    """Classify one dynamic spectrum: RFI-mask, track the drift ridge, measure drift/duration/harm.

    Returns the classification and the measured properties (drift rate, duration, frequency span,
    harmonic score, RFI-masked channel fraction). ``detected`` is True only for ``type_II``.
    """
    keep = rfi_mask(data)
    rt, rf = track_drift_ridge(data, freqs, times, keep=keep, snr_threshold=snr_threshold)
    drift = _robust_drift(rt, rf)
    duration = float(rt.max() - rt.min()) if rt.size else 0.0
    coherence = float(abs(np.corrcoef(rt, rf)[0, 1])) if rt.size >= 3 and np.std(rf) > 0 else 0.0
    harm = harmonic_score(data, freqs, rt, rf, times)
    klass = classify_burst(drift, duration, harm, coherence, int(rt.size))
    return {
        "klass": klass,
        "detected": klass == "type_II",
        "drift_mhz_s": round(drift, 4) if np.isfinite(drift) else None,
        "duration_s": round(duration, 1),
        "coherence": round(coherence, 3),
        "freq_hi_mhz": round(float(rf.max()), 1) if rf.size else None,
        "freq_lo_mhz": round(float(rf.min()), 1) if rf.size else None,
        "harmonic_score": round(harm, 3),
        "rfi_masked_frac": round(1.0 - float(keep.mean()), 3),
        "n_ridge": int(rt.size),
    }


def crossmatch_cme(
    burst_time_hr: float, cme_list: list[dict], *, window_hr: float = 2.0
) -> dict | None:
    """Nearest LASCO CME within +/- ``window_hr`` of the burst (in time), or None.

    A metric type II and its driving CME are near-simultaneous but the catalogue times differ in
    sign: the type II is emitted low in the corona (~1.5-3 Rsun) while CDAW logs the CME's FIRST C2
    APPEARANCE (~2.5-6 Rsun), which trails the burst by ~30-90 min; other conventions put the
    extrapolated onset before it. So we use a symmetric time window (default +/-2 h) and take the
    nearest CME, robust to the offset direction. Returns the CME dict (``speed_kms``, ``width_deg``)
    or None if none is within the window.
    """
    cands = [c for c in cme_list if abs(burst_time_hr - float(c["onset_hr"])) <= window_hr]
    if not cands:
        return None
    return min(cands, key=lambda c: abs(burst_time_hr - float(c["onset_hr"])))


def cme_association_fraction(matched_cmes: list[dict | None]) -> dict:
    """Fast-and-wide-CME fractions of a set of type II-associated CMEs (cf. Gopalswamy+2005).

    Reports the fraction driven by fast (>= ``FAST_CME_KMS``) and wide (>= ``WIDE_CME_DEG``) CMEs
    plus the median speed. On the SYNTHETIC census these fractions echo the injected Gopalswamy-
    biased CME distribution --- so the offline result is a **wiring check** of the temporal
    cross-match (does it pick the right CME among decoys), NOT an independent reproduction of the
    association, which requires the real (Turnstile-blocked) event list.
    """
    m = [c for c in matched_cmes if c is not None]
    if not m:
        return {"n_matched": 0, "frac_fast": None, "frac_wide": None, "median_speed_kms": None}
    speeds = np.array([float(c["speed_kms"]) for c in m])
    widths = np.array([float(c["width_deg"]) for c in m])
    return {
        "n_matched": len(m),
        "frac_fast": round(float((speeds >= FAST_CME_KMS).mean()), 3),
        "frac_wide": round(float((widths >= WIDE_CME_DEG).mean()), 3),
        "frac_fast_and_wide": round(
            float(((speeds >= FAST_CME_KMS) & (widths >= WIDE_CME_DEG)).mean()), 3
        ),
        "median_speed_kms": round(float(np.median(speeds)), 1),
    }


EPOCH_ISO = "2024-01-01T00:00:00"  # common clock: all burst/CME/flare times in hours since here


def _iso_to_hours(iso: str) -> float:
    """ISO-8601 UTC timestamp -> hours since ``EPOCH_ISO`` (the common cross-match clock)."""
    from datetime import datetime

    s = iso.strip().replace("Z", "").split(".")[0]
    fmt = "%Y-%m-%dT%H:%M:%S" if "T" in s else "%Y-%m-%d %H:%M:%S"
    dt = datetime.strptime(s, fmt)
    ep = datetime.strptime(EPOCH_ISO, "%Y-%m-%dT%H:%M:%S")
    return (dt - ep).total_seconds() / 3600.0


def parse_cme_html(html: str) -> list[dict]:
    """Parse one CDAW SOHO/LASCO monthly CME-catalogue HTML page -> [{onset_hr, speed_kms, width_deg}].

    Columns (CDAW UNIVERSAL_ver2): date, time [UT], central PA, angular width [deg], linear speed
    [km/s], ... We strip tags and match ``DATE TIME (CPA|Halo) WIDTH SPEED``; rows with an
    unmeasured speed (``----``) are skipped. Onset time -> hours on the ``EPOCH_ISO`` clock.
    """
    import re

    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
    out = []
    for m in re.finditer(
        r"(\d{4})/(\d{2})/(\d{2}) (\d{2}:\d{2}:\d{2}) (?:\d+|Halo) (\d+) (\d+)", text
    ):
        y, mo, d, tod, width, speed = m.groups()
        out.append(
            {
                "onset_hr": _iso_to_hours(f"{y}-{mo}-{d}T{tod}"),
                "width_deg": float(width),
                "speed_kms": float(speed),
            }
        )
    return out


def parse_hek_flares(payload: dict) -> list[dict]:
    """Parse a HEK flare-search JSON payload -> [{peak_hr, goes_class, goes_flux}].

    HEK returns ``result`` rows with ``event_peaktime`` and ``fl_goescls`` (e.g. ``X1.7``); we keep
    SWPC-sourced flares, convert the class letter+number to a log10 W/m^2 flux for magnitude cuts,
    and the peak time to hours on the ``EPOCH_ISO`` clock.
    """
    band = {"A": -8.0, "B": -7.0, "C": -6.0, "M": -5.0, "X": -4.0}
    out = []
    for r in payload.get("result", []):
        cls = (r.get("fl_goescls") or "").strip()
        pk = r.get("event_peaktime")
        if not cls or not pk or cls[0] not in band:
            continue
        # class letter sets the decade (X=1e-4 W/m^2), the trailing number the mantissa
        try:
            mant = float(cls[1:]) if len(cls) > 1 and float(cls[1:]) > 0 else 1.0
        except ValueError:
            mant = 1.0
        flux = band[cls[0]] + float(np.log10(mant))
        out.append(
            {"peak_hr": _iso_to_hours(pk), "goes_class": cls, "goes_log_flux": round(flux, 2)}
        )
    return out


def crossmatch_flare(
    burst_hr: float, flare_list: list[dict], *, window_hr: float = 1.0
) -> dict | None:
    """Nearest GOES flare whose peak is within [-window, +0.25] hr of the burst (flare leads/accompanies)."""
    cands = [f for f in flare_list if -0.25 <= (burst_hr - float(f["peak_hr"])) <= window_hr]
    if not cands:
        return None
    return min(cands, key=lambda f: abs(burst_hr - float(f["peak_hr"])))


def occurrence_vs_phase(
    burst_hrs: list[float], observing_day_hrs: list[float], sunspot_by_month: dict
) -> dict:
    """Coverage-corrected monthly type II rate vs SILSO sunspot number (reuse `ecallisto_census`).

    Bins detections and OBSERVING days by calendar month (both on the ``EPOCH_ISO`` clock),
    forms rate = detections / observing-days (the `ecallisto_census` coverage correction), and
    correlates it with the SILSO monthly mean sunspot number. With the ~2 yr baseline this is
    indicative only (stated). ``sunspot_by_month`` maps ``"YYYY-MM"`` -> mean SN.
    """
    from datetime import datetime, timedelta

    from .ecallisto_census import census_correlation, coverage_corrected_rate

    ep = datetime.strptime(EPOCH_ISO, "%Y-%m-%dT%H:%M:%S")

    def month_key(hr: float) -> str:
        return (ep + timedelta(hours=hr)).strftime("%Y-%m")

    months = sorted({month_key(h) for h in observing_day_hrs})
    n_events = np.array([sum(month_key(b) == mk for b in burst_hrs) for mk in months], float)
    coverage = np.array(
        [sum(month_key(o) == mk for o in observing_day_hrs) for mk in months], float
    )
    rate = coverage_corrected_rate(n_events, coverage)
    sn = np.array([sunspot_by_month.get(mk, np.nan) for mk in months], float)
    corr = census_correlation(rate, sn)
    return {"n_months": len(months), "months": months, **corr}


def _synthetic_census(seed: int = 0) -> dict:
    """Offline recover-a-known: a mixed event set + a CME list with a known fast-and-wide bias.

    Injects type II bursts (each driven by a CME whose speed/width follow the Gopalswamy bias),
    type III contaminants, and RFI-only spectra; runs the detector on each and cross-matches the
    detected type IIs to the CME list. Returns detector completeness/purity + the recovered
    association fractions, all of which the census must reproduce.
    """
    from .solarbursts import synthetic_burst

    rng = np.random.default_rng(seed)
    events, truth = [], []
    cme_list = []
    n_ii, n_iii, n_rfi = 24, 16, 8
    # a realistic MIX of type II difficulty, so completeness/purity is not a saturated easy number:
    # a third strong+harmonic (easy), a third single-lane (no harmonic), a third weak-SNR (amp~4)
    for i in range(n_ii):
        speed = float(rng.normal(1200, 350))  # fast-and-wide-biased CME (Gopalswamy)
        width = float(np.clip(rng.normal(110, 45), 10, 360))
        onset = float(rng.uniform(0, 2000))  # sparse over a realistic baseline (hours)
        cme_list.append({"onset_hr": onset, "speed_kms": speed, "width_deg": width})
        # mixed difficulty: strong+harmonic (easy), single-lane, near-threshold-SNR
        kind = i % 3
        with_harm = kind != 1
        amp = 2.5 if kind == 2 else 10.0
        s = synthetic_typeii(
            shock_speed_kms=max(400.0, speed),
            with_harmonic=with_harm,
            amp=amp,
            seed=rng.integers(1 << 30).item(),
        )
        s["burst_hr"] = onset + rng.uniform(0.05, 0.4)  # burst follows CME onset
        events.append(s)
        truth.append("type_II")
    for _ in range(n_iii):  # fast type III contaminants
        events.append(synthetic_burst(speed_c=0.3, seed=rng.integers(1 << 30).item()))
        truth.append("type_III")
    for _ in range(n_rfi):  # RFI-only: horizontal narrowband lines
        f = np.linspace(88, 20, 240)
        d = rng.normal(0, 1.0, (240, 600))
        for ch in rng.choice(240, 4, replace=False):
            d[ch] += rng.uniform(6, 12)  # persistent bright channel
        events.append({"data": d, "freqs": f, "times": np.linspace(0, 300, 600)})
        truth.append("rfi")
    # add unrelated (decoy) CMEs so cross-match must discriminate, not just pair by order
    for _ in range(30):
        cme_list.append(
            {
                "onset_hr": float(rng.uniform(0, 2000)),
                "speed_kms": float(rng.normal(500, 200)),
                "width_deg": float(np.clip(rng.normal(40, 20), 5, 360)),
            }
        )
    return {
        "events": events,
        "truth": truth,
        "cme_list": cme_list,
        "n_injected_typeii": n_ii,
    }


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: synthetic census (detector completeness/purity + recovered CME bias). Real: OVRO-LWA."""
    import json

    if offline:
        cen = _synthetic_census()
        results = [detect_typeii(e["data"], e["freqs"], e["times"]) for e in cen["events"]]
        truth = cen["truth"]
        det_ii = [i for i, r in enumerate(results) if r["detected"]]
        true_ii = [i for i, t in enumerate(truth) if t == "type_II"]
        tp = len([i for i in det_ii if truth[i] == "type_II"])
        completeness = tp / max(len(true_ii), 1)
        purity = tp / max(len(det_ii), 1)
        matched = [
            crossmatch_cme(cen["events"][i]["burst_hr"], cen["cme_list"])
            for i in det_ii
            if truth[i] == "type_II"
        ]
        assoc = cme_association_fraction(matched)
        # completeness vs injected SNR (amplitude / noise): the honest performance curve, so the
        # headline number is not a single saturated value from only-easy injections
        curve = {}
        for amp in (2.0, 2.5, 3.0, 4.0):
            hits = sum(
                detect_typeii(
                    *[
                        synthetic_typeii(amp=amp, seed=1000 + k)[x]
                        for x in ("data", "freqs", "times")
                    ]
                )["detected"]
                for k in range(24)
            )
            curve[f"snr{amp:g}"] = round(hits / 24, 3)
        source = "synthetic mixed census (type II + III + RFI; Gopalswamy CME bias injected)"
        metrics: dict = {
            "source": source,
            "is_real": False,
            "n_events": len(results),
            "n_typeii_detected": len(det_ii),
            "completeness": round(completeness, 3),
            "purity": round(purity, 3),
            "completeness_vs_snr": curve,
            "median_rfi_masked_frac": round(
                float(np.median([r["rfi_masked_frac"] for r in results])), 3
            ),
            **{f"assoc_{k}": v for k, v in assoc.items()},
        }
    else:  # pragma: no cover - the streaming real census needs explicit dates + a CME table
        raise SystemExit(
            "The real type II census streams days from AWS Open Data and needs explicit dates + a "
            "LASCO CME table -- run `python scripts/typeii_real.py --dates ... --cme ...` instead of "
            "`--offline`."
        )

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "typeii_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(metrics, op / "papers" / "typeii" / "figures")
    _write_macros(metrics, op / "papers" / "typeii" / "generated" / "macros.tex")
    return metrics


TIME_DOWNSAMPLE = 16  # average 16 native ~0.256 s samples -> ~4 s bins (type II is minutes-long)
CADENCE_S = 0.256  # native OVRO-LWA dynamic-spectrum time resolution
SWEEP_WINDOW_S = 900.0  # slide a 15-min window across a day (type II lasts minutes)
SWEEP_STEP_S = 450.0


def _downsample_time(spec: np.ndarray, factor: int) -> np.ndarray:
    """Block-average a (freq, time) spectrum over ``factor`` adjacent time columns.

    A type II lasts minutes, so averaging to ~4 s bins loses nothing and shrinks a full-day
    spectrum ~16x --- the key to holding a streamed day in memory.
    """
    if factor <= 1:
        return spec
    n = (spec.shape[1] // factor) * factor
    return spec[:, :n].reshape(spec.shape[0], n // factor, factor).mean(axis=2)


def _axes_from_header(hdr: dict, n_freq: int, n_time: int, *, factor: int = 1) -> tuple:
    """Descending freq (MHz) + time (s) axes from a dynamic-spectrum header (post-downsample)."""
    fmin = float(hdr.get("FREQMIN", 0.0134)) * 1e3  # GHz -> MHz
    fmax = float(hdr.get("FREQMAX", 0.0869)) * 1e3
    freqs = np.linspace(fmin, fmax, n_freq)
    times = np.arange(n_time) * CADENCE_S * factor
    if freqs[0] < freqs[-1]:  # descending, as solarbursts expects
        freqs = freqs[::-1]
    return freqs, times


def sweep_day(
    data: np.ndarray,
    freqs: np.ndarray,
    times: np.ndarray,
    *,
    window_s: float = SWEEP_WINDOW_S,
    step_s: float = SWEEP_STEP_S,
) -> list[dict]:
    """Slide a burst-sized window across a full day and run `detect_typeii` in each.

    A day-long dynamic spectrum holds many minutes-long candidate intervals; `detect_typeii` is a
    single-window classifier, so we scan overlapping ``window_s`` windows (step ``step_s``) and
    return the type II detections, each tagged with its window-centre time in SECONDS from the start
    of the day (``t_center_s``); `real_census` adds the day's absolute offset for the cross-match.
    Descending frequency is assumed (as from `stream_dspec`).
    """
    t = np.asarray(times, float)
    if t.size < 2:
        return []
    out = []
    t0 = float(t.min())
    while t0 < t.max():
        m = (t >= t0) & (t < t0 + window_s)
        if m.sum() >= N_RIDGE_MIN:
            r = detect_typeii(data[:, m], freqs, t[m] - t[m][0])
            if r["detected"]:
                r["t_center_s"] = float(t[m].mean())
                out.append(r)
        t0 += step_s
    return out


def parse_lwa_dspec(path: str | Path) -> dict:
    """Parse one OVRO-LWA solar dynamic-spectrum FITS -> (Stokes-I data, freqs MHz, times s).

    Confirmed on a real product (AWS Open Data ``spec_fits/YYYY/YYYYMMDD.fits``, 2026): a 4D
    PRIMARY array with FITS axes (NAXIS1=time, NAXIS2=freq, NAXIS3=1, NAXIS4=stokes), i.e. numpy
    shape ``(stokes, 1, freq, time)`` in Jy; the frequency axis spans header ``FREQMIN``--``FREQMAX``
    (GHz) over NAXIS2 channels and the time axis ``DATE_OBS``--``DATE_END``. We take Stokes I (the
    first plane), build the freq/time axes from the header, and return the spectrum oriented
    (freq, time) with descending frequency (as `solarbursts` expects). Falls back to the older
    table layout (2D spectrum HDU + 1D freq/time HDUs) if the primary array is not 4D.
    """
    from astropy.io import fits

    with fits.open(path) as hdul:
        prim = hdul[0]
        if getattr(prim, "data", None) is not None and prim.data.ndim == 4:
            spec = np.asarray(prim.data[0, 0], float)  # Stokes I -> (freq, time)
            hdr = prim.header
            n_f = spec.shape[0]
            fmin = float(hdr.get("FREQMIN", 0.0134)) * 1e3  # GHz -> MHz
            fmax = float(hdr.get("FREQMAX", 0.0869)) * 1e3
            freqs = np.linspace(fmin, fmax, n_f)
            times = np.linspace(0.0, spec.shape[1] * 0.256, spec.shape[1])
        else:  # older documented table layout: 2D spectrum HDU + 1D freq/time HDUs
            spec = np.asarray(
                next(
                    h.data
                    for h in hdul
                    if getattr(h, "data", None) is not None and h.data.ndim == 2
                ),
                float,
            )
            arrays = [
                np.asarray(h.data, float).ravel()
                for h in hdul
                if getattr(h, "data", None) is not None and h.data.ndim == 1
            ]
            freqs = next(
                (a for a in arrays if a.size in spec.shape), np.linspace(13.4, 86.9, spec.shape[0])
            )
            times = next(
                (a for a in arrays if a.size in spec.shape and a is not freqs),
                np.arange(spec.shape[1]) * 0.256,
            )
            if spec.shape[0] != freqs.size and spec.shape[1] == freqs.size:
                spec = spec.T
            if freqs.max() > 1e3:
                freqs = freqs / 1e6  # Hz -> MHz
    if freqs[0] < freqs[-1]:  # descending frequency, as solarbursts expects
        spec, freqs = spec[::-1], freqs[::-1]
    return {"data": spec, "freqs": freqs, "times": times}


S3_SPEC_FITS = "https://ovro-lwa-solar.s3-us-west-2.amazonaws.com/spec_fits"


def s3_dspec_url(date: str) -> str:
    """Public AWS-Open-Data URL for a daily OVRO-LWA dynamic spectrum (``date`` = YYYY-MM-DD).

    The FITS are on AWS Open Data (bucket ``ovro-lwa-solar``), directly downloadable with no login
    and no bot challenge --- the Cloudflare Turnstile only gates the query *UI*, not the data. Path:
    ``spec_fits/<YYYY>/<YYYYMMDD>.fits`` (one/day, ~1.7 GB, 15-85 MHz, ~0.26 s cadence).
    """
    ymd = date.replace("-", "")
    return f"{S3_SPEC_FITS}/{ymd[:4]}/{ymd}.fits"


def stream_dspec(
    date: str, *, time_downsample: int = TIME_DOWNSAMPLE, chunk: int = 128
) -> dict:  # pragma: no cover - network (astropy use_fsspec range reads)
    """Stream one day's dynamic spectrum from AWS Open Data into memory --- no file on disk.

    Opens the S3 FITS lazily (astropy ``use_fsspec``; needs the ``typeii`` extra: fsspec+aiohttp)
    and reads ONLY the Stokes-I plane, in ``chunk``-channel FREQUENCY blocks via HTTP range
    requests. The plane is stored C-contiguous with time fastest-varying, so a frequency slice
    ``[f0:f1, :]`` is one CONTIGUOUS byte range (fast); a time slice would be strided across every
    channel (many small requests). Each block is time-averaged to ~4 s bins before it is kept, so
    peak memory is one raw frequency block (~150 MB) plus the reduced spectrum, never the ~1.7 GB
    file. Returns the reduced (freq, time) spectrum + axes, ready for `sweep_day`.
    """
    from astropy.io import fits

    url = s3_dspec_url(date)
    with fits.open(url, use_fsspec=True) as hdul:
        hdu = hdul[0]
        hdr = dict(hdu.header)
        n_f = int(hdu.shape[-2])  # FITS axis-2 = frequency
        blocks = []
        for a in range(0, n_f, chunk):
            b = min(a + chunk, n_f)
            block = np.asarray(hdu.section[0, 0, a:b, :], float)  # Stokes I, (freq_chunk, time)
            blocks.append(_downsample_time(block, time_downsample))
            del block
        spec = np.concatenate(blocks, axis=0)
    freqs, times = _axes_from_header(hdr, spec.shape[0], spec.shape[1], factor=time_downsample)
    if freqs[0] < freqs[-1]:
        spec = spec[::-1]
    date_obs = str(hdr.get("DATE_OBS", f"{date}T00:00:00"))  # absolute obs start for cross-match
    return {"data": spec, "freqs": freqs, "times": times, "date_obs": date_obs}


def purity_diagnostics(
    detections: list[dict],
    cme_list: list[dict],
    n_days: int,
    *,
    window_hr: float = 2.0,
    window_s: float = SWEEP_WINDOW_S,
    seed: int = 0,
) -> dict:
    """Is the candidate list real type II, or false positives? The honest purity test.

    Real type II are driven by fast (>=900 km/s) wide CMEs (Gopalswamy+2005), so a real sample's
    matched-CME speeds should sit FAR above the background CME population. Three independent tests:

    (1) **Speed** --- matched-CME median vs the background CME median; and the observed CME-match
    rate vs the CHANCE rate (random times, any CME within ``window_hr``). At solar max the CME
    density makes ``chance`` high, so this test is weak on its own and the SPEED comparison carries
    the argument (`association_is_background_like`).
    (2) **Decorrelation** --- the burst's own drift rate encodes the shock kinematics, so a real
    sample would show ``|drift|`` correlated with the matched CME speed. A near-zero
    ``drift_cme_speed_corr`` means the radio properties know nothing about the CME --- the "matches"
    are coincidence (any flare-coincident blob auto-selects a nearby fast CME).
    (3) **Window saturation** --- real type II are minutes-long transients of VARIED duration; a
    high ``window_saturation_frac`` (ridges filling the whole detection window) is persistent
    RFI/background structure tracked as a slow drift, the classic false-positive signature.
    """
    if not cme_list or not detections:
        return {"detection_rate_per_day": None}
    onsets = np.array(sorted(float(c["onset_hr"]) for c in cme_list))
    bg_speed = np.array([float(c["speed_kms"]) for c in cme_list])
    matched = [crossmatch_cme(r["burst_hr"], cme_list, window_hr=window_hr) for r in detections]
    m_speed = np.array([float(c["speed_kms"]) for c in matched if c is not None])
    rng = np.random.default_rng(seed)
    rand = rng.uniform(onsets.min(), onsets.max(), 5000)
    chance = float(np.mean([np.any(np.abs(onsets - t) <= window_hr) for t in rand]))
    obs_rate = sum(c is not None for c in matched) / len(detections)
    bg_med = float(np.median(bg_speed))
    m_med = float(np.median(m_speed)) if m_speed.size else float("nan")
    # (2) does the burst drift know anything about the matched CME speed? (real type II: yes)
    pairs = [
        (abs(float(r["drift_mhz_s"])), float(c["speed_kms"]))
        for r, c in zip(detections, matched, strict=True)
        if c is not None and r.get("drift_mhz_s") is not None
    ]
    if len(pairs) >= 5:
        dr, sp = np.array(pairs).T
        drift_speed_corr = (
            float(np.corrcoef(dr, sp)[0, 1]) if np.std(dr) > 0 and np.std(sp) > 0 else float("nan")
        )
    else:
        drift_speed_corr = float("nan")
    # (3) fraction of detections whose ridge fills >=98% of the detection window (FP signature)
    durs = [float(r["duration_s"]) for r in detections if r.get("duration_s") is not None]
    sat = float(np.mean([d >= 0.98 * window_s for d in durs])) if durs else float("nan")
    return {
        "detection_rate_per_day": round(len(detections) / max(n_days, 1), 3),
        "bg_cme_median_kms": round(bg_med, 1),
        "bg_cme_frac_fast": round(float((bg_speed >= FAST_CME_KMS).mean()), 3),
        "matched_cme_median_kms": round(m_med, 1) if np.isfinite(m_med) else None,
        "chance_cme_match_rate": round(chance, 3),
        "observed_cme_match_rate": round(obs_rate, 3),
        "drift_cme_speed_corr": round(drift_speed_corr, 3)
        if np.isfinite(drift_speed_corr)
        else None,
        "window_saturation_frac": round(sat, 3) if np.isfinite(sat) else None,
        # matches ~ background population, radio props decorrelated from CME speed, and ridges
        # saturate the window -> the candidate list is false-positive-dominated
        "association_is_background_like": bool(
            np.isfinite(m_med) and m_med < 1.5 * bg_med and obs_rate <= chance
        ),
    }


def fetch_lasco_cme(dates: list[str]) -> list[dict]:  # pragma: no cover - network (CDAW HTML)
    """CDAW SOHO/LASCO CMEs for the months spanning ``dates`` -> [{onset_hr, speed_kms, width_deg}]."""
    import urllib.request

    out: list[dict] = []
    for ym in sorted({d[:7] for d in dates}):  # YYYY-MM
        y, mo = ym.split("-")
        url = f"https://cdaw.gsfc.nasa.gov/CME_list/UNIVERSAL_ver2/{y}_{mo}/univ{y}_{mo}.html"
        try:
            html = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "replace")
            out += parse_cme_html(html)
        except Exception as exc:  # noqa: BLE001
            print(f"  CDAW {ym} unavailable: {exc!r}", flush=True)
    return out


def fetch_goes_flares(dates: list[str]) -> list[dict]:  # pragma: no cover - network (HEK)
    """SWPC GOES flares over the span of ``dates`` from the HEK -> [{peak_hr, goes_class, ...}]."""
    import json
    import urllib.parse
    import urllib.request

    lo, hi = min(dates), max(dates)
    q = {
        "cosec": "2",
        "cmd": "search",
        "type": "column",
        "event_type": "fl",
        "event_starttime": f"{lo}T00:00:00",
        "event_endtime": f"{hi}T23:59:59",
        "event_coordsys": "helioprojective",
        "x1": "-1200",
        "x2": "1200",
        "y1": "-1200",
        "y2": "1200",
        "result_limit": "5000",
        "return": "fl_goescls,event_peaktime,frm_name",
        "param0": "FRM_NAME",
        "op0": "=",
        "value0": "SWPC",
    }
    url = "https://www.lmsal.com/hek/her?" + urllib.parse.urlencode(q)
    try:
        payload = json.loads(
            urllib.request.urlopen(url, timeout=120).read().decode("utf-8", "replace")
        )
        return parse_hek_flares(payload)
    except Exception as exc:  # noqa: BLE001
        print(f"  HEK flares unavailable: {exc!r}", flush=True)
        return []


def fetch_sunspots_by_month() -> dict:  # pragma: no cover - network (SILSO)
    """SILSO monthly mean sunspot number -> {"YYYY-MM": SN} (the cycle-phase axis)."""
    import urllib.request

    from .ecallisto_census import SILSO_URL

    txt = urllib.request.urlopen(SILSO_URL, timeout=60).read().decode("utf-8", "replace")
    out = {}
    for line in txt.splitlines():
        parts = line.replace(";", " ").split()
        if len(parts) >= 4:
            try:
                out[f"{int(parts[0]):04d}-{int(parts[1]):02d}"] = float(parts[3])
            except ValueError:
                continue
    return out


def real_census(dates, cme_list=None, flare_list=None, sunspots=None):  # noqa: ANN001  # pragma: no cover
    """Stream each day, sweep for type II bursts, cross-match LASCO CMEs + GOES flares, occurrence.

    ``dates`` is ``YYYY-MM-DD`` strings. The catalogues are fetched if not supplied (CDAW CMEs, HEK
    GOES flares, SILSO sunspots). Each day is streamed into memory (`stream_dspec`), swept
    (`sweep_day`), and freed before the next --- no ~1.7 GB file touches disk. Detections carry an
    absolute time (day ``DATE_OBS`` + window centre), matched to the nearest preceding CME and
    flare; the coverage-corrected monthly rate is correlated with the sunspot number. Emits the full
    plan product set (event list + CME + GOES association + occurrence vs cycle phase).
    """
    cme_list = fetch_lasco_cme(dates) if cme_list is None else cme_list
    flare_list = fetch_goes_flares(dates) if flare_list is None else flare_list
    sunspots = fetch_sunspots_by_month() if sunspots is None else sunspots

    all_det: list[dict] = []
    observing_day_hrs: list[float] = []
    n_failed = 0
    for date in dates:
        try:
            ds = stream_dspec(date)  # streamed, in memory
            day_hr = _iso_to_hours(ds["date_obs"])
            observing_day_hrs.append(day_hr)
            for r in sweep_day(ds["data"], ds["freqs"], ds["times"]):
                r["burst_hr"] = day_hr + r["t_center_s"] / 3600.0
                all_det.append(r)
            del ds  # free the day before streaming the next
        except Exception as exc:  # noqa: BLE001 - one bad/corrupt day must not abort a long run
            n_failed += 1
            print(f"  {date}: FAILED {type(exc).__name__}: {exc}", flush=True)
            continue
        print(f"  {date}: {len(all_det)} type II so far", flush=True)

    matched_cme = [crossmatch_cme(r["burst_hr"], cme_list) for r in all_det]
    matched_fl = [crossmatch_flare(r["burst_hr"], flare_list) for r in all_det]
    assoc = cme_association_fraction(matched_cme)
    occ = occurrence_vs_phase([r["burst_hr"] for r in all_det], observing_day_hrs, sunspots)
    purity = purity_diagnostics(all_det, cme_list, len(observing_day_hrs))
    event_list = [
        {
            "burst_hr": round(r["burst_hr"], 3),
            "drift_mhz_s": r["drift_mhz_s"],
            "duration_s": r["duration_s"],
            "harmonic_score": r["harmonic_score"],
            "cme": matched_cme[i],
            "goes_class": matched_fl[i]["goes_class"] if matched_fl[i] else None,
        }
        for i, r in enumerate(all_det)
    ]
    n_fl = sum(1 for f in matched_fl if f is not None)
    return {
        "source": f"OVRO-LWA dspec (AWS Open Data, streamed), {len(dates)} days",
        "is_real": True,
        "n_events": len(dates),
        "n_days_processed": len(observing_day_hrs),
        "n_days_failed": n_failed,
        "n_typeii_detected": len(all_det),
        "completeness": None,
        "purity": None,
        "n_flare_associated": n_fl,
        "frac_flare_associated": round(n_fl / max(len(all_det), 1), 3),
        "occ_n_months": occ["n_months"],
        "occ_pearson_r": round(occ["pearson_r"], 3) if np.isfinite(occ["pearson_r"]) else None,
        "occ_spearman_rho": round(occ["spearman_rho"], 3)
        if np.isfinite(occ["spearman_rho"])
        else None,
        "event_list": event_list,
        **purity,
        **{f"assoc_{k}": v for k, v in assoc.items()},
    }


def _figure(m: dict, out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    # a representative synthetic type II dynamic spectrum with its tracked ridge
    s = synthetic_typeii(seed=1)
    keep = rfi_mask(s["data"])
    rt, rf = track_drift_ridge(s["data"], s["freqs"], s["times"], keep=keep)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    ax1.imshow(
        background_subtract(s["data"]),
        aspect="auto",
        origin="lower",
        extent=[s["times"][0], s["times"][-1], s["freqs"][-1], s["freqs"][0]],
        cmap="viridis",
    )
    ax1.plot(rt, rf, ".", color="r", ms=2, label="tracked ridge")
    ax1.set(xlabel="time (s)", ylabel="freq (MHz)", title="Synthetic type II + ridge")
    ax1.legend(fontsize=8)
    labels = ["complete.", "purity", "frac fast", "frac wide"]
    vals = [
        m.get("completeness"),
        m.get("purity"),
        m.get("assoc_frac_fast"),
        m.get("assoc_frac_wide"),
    ]
    vals = [v if isinstance(v, (int, float)) else 0.0 for v in vals]
    ax2.bar(labels, vals, color=["C0", "C0", "C3", "C3"])
    ax2.set(ylim=(0, 1.05), ylabel="fraction", title="Detector + CME association")
    fig.tight_layout()
    fig.savefig(out / "typeii.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    def g(key: str) -> str:
        v = m.get(key)
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    pref = "tiiReal" if m.get("is_real") else "tiiSyn"
    lines = [
        "% Auto-generated by jansky_research.typeii._write_macros -- do not edit.",
        "% Synthetic (tiiSyn*) and real (tiiReal*) namespaces are BOTH emitted; the inactive one",
        "% holds placeholders (an offline rebuild resets tiiReal* to '--').",
        rf"\newcommand{{\tiiSource}}{{{m['source']}}}",
        rf"\newcommand{{\tiiNEvents}}{{{g('n_events')}}}",
        rf"\newcommand{{\tiiCompleteness}}{{{g('completeness')}}}",
        rf"\newcommand{{\tiiPurity}}{{{g('purity')}}}",
    ]
    curve = m.get("completeness_vs_snr") or {}
    # LaTeX control sequences are letters-only -> spell out the SNR values
    for snr, tag in (
        ("snr2", "Snrtwo"),
        ("snr2.5", "Snrtwofive"),
        ("snr3", "Snrthree"),
        ("snr4", "Snrfour"),
    ):
        v = curve.get(snr)
        lines.append(rf"\newcommand{{\tiiComp{tag}}}{{{'--' if v is None else v}}}")
    for ns in ("tiiSyn", "tiiReal"):
        live = ns == pref
        for macro, key in (
            ("NTypeII", "n_typeii_detected"),
            ("FracFast", "assoc_frac_fast"),
            ("FracWide", "assoc_frac_wide"),
            ("FracFastWide", "assoc_frac_fast_and_wide"),
            ("MedSpeed", "assoc_median_speed_kms"),
            ("NMatched", "assoc_n_matched"),
            ("FracFlare", "frac_flare_associated"),
            ("OccMonths", "occ_n_months"),
            ("OccPearson", "occ_pearson_r"),
            ("OccSpearman", "occ_spearman_rho"),
            ("DetRate", "detection_rate_per_day"),
            ("BgCmeMed", "bg_cme_median_kms"),
            ("BgFracFast", "bg_cme_frac_fast"),
            ("MatchCmeMed", "matched_cme_median_kms"),
            ("ChanceMatch", "chance_cme_match_rate"),
            ("ObsMatch", "observed_cme_match_rate"),
            ("DriftSpeedCorr", "drift_cme_speed_corr"),
            ("WindowSat", "window_saturation_frac"),
        ):
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(key) if live else '--'}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="OVRO-LWA type II burst census.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
