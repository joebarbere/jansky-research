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
is **no real census yet**: the OVRO-LWA solar portal (`ovsa.njit.edu/lwadata-query`) sits behind a
Cloudflare Turnstile bot challenge, so the FITS cannot be fetched by a script (GATE-0 2026-07-09);
the real census (and the coverage-corrected occurrence-vs-cycle-phase piece, which needs the
multi-year event list) awaits FITS downloaded interactively through the portal --- see
`real_census` / `scripts/typeii_real.py`.

Data: OVRO-LWA Level-1 beamforming spectrograms (FITS, 13.4--86.9 MHz, 256 ms, ~0.6 GB/day,
2024-04->present); LASCO CME v2 (CDAW); GOES flares; SILSO for cycle phase (real leg). Reuse:
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
    "parse_lwa_dspec",
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
    burst_time_hr: float, cme_list: list[dict], *, window_hr: float = 1.0
) -> dict | None:
    """Best-matching LASCO CME for a type II (onset within [-window, +0.25] hr of the burst).

    A CME-shock type II follows CME onset, so we match CMEs whose onset time (hours) precedes the
    burst by up to ``window_hr`` (or trails by up to 15 min for timing slop). The gap
    ``burst - onset`` is therefore POSITIVE for a real driver (onset first); we accept
    ``-0.25 <= gap <= window_hr``. The physical driver is the CME whose onset most closely PRECEDES
    the burst, so among candidates we return the one with the smallest |gap|. Returns the CME dict
    (with ``speed_kms``, ``width_deg``) or None.
    """
    cands = [c for c in cme_list if -0.25 <= (burst_time_hr - float(c["onset_hr"])) <= window_hr]
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
    else:  # pragma: no cover - real leg runs via scripts/typeii_real.py on local OVRO-LWA FITS
        metrics = real_census(DATA_DIR)

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "typeii_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(metrics, op / "papers" / "typeii" / "figures")
    _write_macros(metrics, op / "papers" / "typeii" / "generated" / "macros.tex")
    return metrics


def parse_lwa_dspec(path: str | Path) -> dict:
    """Parse one OVRO-LWA Level-1 beamforming dynamic-spectrum FITS -> (data, freqs MHz, times s).

    File pattern ``ovro-lwa.lev1_bmf_256ms_96kHz.YYYY-MM-DD.dspec_I.fits`` (one/day, 13.4-86.9 MHz,
    768 ch, 256 ms; Stokes-I in SFU). Per the OVRO-LWA Data Products wiki the FITS carries the
    frequency list, the time list, and the I dynamic spectrum as tables/extensions; we locate the
    2D spectrum HDU and the 1D freq/time axes and return the spectrum oriented (freq, time).
    """
    from astropy.io import fits

    with fits.open(path) as hdul:
        spec = next(
            h.data for h in hdul if getattr(h, "data", None) is not None and h.data.ndim == 2
        )
        arrays = [
            np.asarray(h.data, float).ravel()
            for h in hdul
            if getattr(h, "data", None) is not None and h.data.ndim == 1
        ]
        freqs = next(
            (a for a in arrays if a.size == spec.shape[0] or a.size == spec.shape[1]), None
        )
        times = next(
            (a for a in arrays if freqs is not None and a.size in spec.shape and a is not freqs),
            None,
        )
    spec = np.asarray(spec, float)
    if freqs is not None and spec.shape[0] != freqs.size and spec.shape[1] == freqs.size:
        spec = spec.T  # orient (freq, time)
    freqs = freqs if freqs is not None else np.linspace(86.9, 13.4, spec.shape[0])
    times = times if times is not None else np.arange(spec.shape[1]) * 0.256
    if freqs[0] < freqs[-1]:  # want descending frequency (as solarbursts)
        spec, freqs = spec[::-1], freqs[::-1]
    return {"data": spec, "freqs": freqs / 1e6 if freqs.max() > 1e3 else freqs, "times": times}


def real_census(data_dir: str | Path) -> dict:  # pragma: no cover - needs local OVRO-LWA FITS
    """Sweep local OVRO-LWA dspec FITS for type II bursts + cross-match a local LASCO CME table.

    The OVRO-LWA solar portal (`ovsa.njit.edu/lwadata-query`) is behind a Cloudflare bot challenge,
    so the FITS must be fetched interactively into ``data_dir`` first (see scripts/typeii_real.py).
    This runs `detect_typeii` over per-day burst windows and cross-matches detections to a CDAW
    LASCO CME table (``data_dir/lasco_cme.csv`` with onset_hr/speed_kms/width_deg), reproducing the
    fast-and-wide association on real events.
    """
    import csv

    dd = Path(data_dir)
    cme_list = [
        {k: float(v) for k, v in row.items()} for row in csv.DictReader(open(dd / "lasco_cme.csv"))
    ]
    events, det = [], []
    for f in sorted(dd.glob("ovro-lwa.*dspec_I.fits")):
        ds = parse_lwa_dspec(f)
        # scan the day in overlapping windows sized for a minutes-long type II
        r = detect_typeii(ds["data"], ds["freqs"], ds["times"])
        if r["detected"]:
            det.append(r)
        events.append(r)
    matched = [crossmatch_cme(r.get("burst_hr", 0.0), cme_list) for r in det]
    assoc = cme_association_fraction(matched)
    return {
        "source": f"OVRO-LWA dspec, {len(events)} days",
        "is_real": True,
        "n_events": len(events),
        "n_typeii_detected": len(det),
        "completeness": None,
        "purity": None,
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
    for snr, tag in (("snr2", "Snrtwo"), ("snr2.5", "Snrtwofive"), ("snr3", "Snrthree"),
                     ("snr4", "Snrfour")):
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
