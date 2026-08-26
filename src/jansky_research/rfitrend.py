"""e-Callisto as an accidental 15-year RFI observatory: the megaconstellation trend (plan 54, F17).

The e-Callisto solar-spectrograph network has recorded 45--870 MHz dynamic spectra worldwide since
2002. The only published RFI-trend study on it (Prieto et al. 2020, Sol. Phys. 295, 11 --- with
Bussons Gordo, Rodriguez-Pacheco et al.; an earlier draft mis-attributed this to a "Perez-Torres"
who is not an author) is a two-epoch campaign finding increased interference levels between 2012
and 2019 at the Spanish sites --- and it stops in 2019, the year
Starlink deployment scaled. Nobody has mined the *continuous* archive across the megaconstellation
era. This slice does: a robust, solar-burst-immune occupancy metric over configuration-stable
stations 2012--2026, trended in the Starlink unintended-emission (UEM) band and attributed to the
public Starlink constellation-count time series (GATE-0 2026-07-10, novelty PASS).

**The load-bearing idea (the systematics gate).** e-Callisto data is uncalibrated digital
log-power that drifts with station gain, and solar bursts are broadband transients --- both would
fake or mask a trend. Both are **common-mode**: they move the whole 45--870 MHz spectrum together.
So the metric is a **differential** --- the occupancy level in the satellite UEM band MINUS the
level in the FM-broadcast control band (87.5--108 MHz, a strong terrestrial fixed transmitter that
should not trend) at the SAME station. Gain drift and solar bursts cancel in the difference; only
band-specific satellite RFI survives. A frequency-resolved test (the intrinsic Starlink narrowband
lines at 125/135/150/175 MHz, Di Vruno+2023, vs adjacent clean channels) is the cleanest
Starlink-specific attribution.

**Honest scope.** A real RFI rise almost certainly exists (Perez+2020 documents it to 2019); the
question the archive can answer is whether the *post-2019, UEM-band-specific* rise tracks Starlink
growth. Attribution to Starlink specifically is bounded --- general terrestrial RFI also grows, the
143.05 MHz feature is reflected GRAVES radar (not intrinsic, excluded), and the self-normalizing
line test still cancels only *common-mode* gain, not differential local flank contamination --- so
we report the trend + correlation with honest attribution caveats, not causation, and gate every
trend claim on cross-station coherence.

Data: the open e-Callisto FITS archive (`solarbursts.fetch_ecallisto` + `ecallisto_catalog`);
Starlink counts from the public planet4589 time series; coverage via `ecallisto_census`. Reuse is
near-total.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "UEM_BAND_MHZ",
    "UEM_LINES_MHZ",
    "FM_CONTROL_MHZ",
    "occupancy_metric",
    "band_differential",
    "pick_control_band",
    "available_lines",
    "trend_fit",
    "starlink_count",
    "line_vs_adjacent",
    "synthetic_month_stack",
    "summarize_stations",
    "run",
]

# Starlink unintended-emission window in e-Callisto's 45-870 MHz range: LOFAR (Di Vruno+2023) sees
# broadband UEM across 110-188 MHz and narrowband features INTRINSIC to the satellites at 125, 135,
# 150 and 175 MHz. The 143.05 MHz feature is EXCLUDED -- Di Vruno+2023 attribute it to reflected
# GRAVES space-surveillance-radar signals, not intrinsic Starlink emission.
UEM_BAND_MHZ = (110.0, 188.0)
UEM_LINES_MHZ = (125.0, 135.0, 150.0, 175.0)  # Di Vruno+2023 intrinsic Starlink narrowband lines
GRAVES_MHZ = (
    142.0,
    144.0,
)  # French GRAVES radar (143.05 MHz) reflections -- NOT intrinsic; excluded
FM_CONTROL_MHZ = (88.0, 108.0)  # FM broadcast: the cleanest flat terrestrial control
ECALLISTO_STATIONS = (
    "HUMAIN",
    "ALMATY",
    "GLASGOW",
)  # 45-870 MHz, span 2012-2026, sample the UEM band


def occupancy_metric(data: np.ndarray, freqs: np.ndarray, band_mhz: tuple[float, float]) -> float:
    """Persistent RFI level in a frequency band (median per-channel time-median), burst-immune.

    e-Callisto values are uncalibrated log-power. Per channel we take the time MEDIAN --- the
    persistent occupied level --- which a transient solar burst or satellite pass (a small fraction
    of the 15-min sweep) does not move; the band level is the median of those channel floors. The
    raw value carries the station's gain, removed by `band_differential`. Returns NaN if the band is
    outside coverage.
    """
    f = np.asarray(freqs, float)
    d = np.asarray(data, float)
    sel = (f >= band_mhz[0]) & (f <= band_mhz[1])
    if sel.sum() == 0:
        return float("nan")
    floor = np.median(d[sel], axis=1)  # per-channel persistent level (burst-immune)
    return float(np.median(floor))


# Candidate control bands in preference order (name, lo, hi). FM is the cleanest fixed-transmitter
# control, but MANY e-Callisto stations notch out the RFI-heavy FM band entirely (HUMAIN samples no
# FM channels: its committed control_name is "low"; the 84->112 MHz gap figure seen during
# smoke-testing is unlogged and is not committed evidence) -- so the control is picked per station
# from whichever candidate the instrument actually samples. All lie OUTSIDE the 110-170 UEM window.
_CONTROL_CANDIDATES = (
    ("FM", 88.0, 108.0),  # FM broadcast -- cleanest, but often notched out
    ("low", 55.0, 80.0),  # fixed/aeronautical below the UEM window
    ("high", 300.0, 420.0),  # UHF above the UEM window (some TV-switchover risk)
)
_MIN_CONTROL_CH = 4  # a control band needs at least this many sampled channels to be usable


def pick_control_band(freqs: np.ndarray) -> tuple[str, tuple[float, float]]:
    """Best-sampled clean control band for THIS station's channel grid (name, (lo, hi)).

    e-Callisto instruments sample station-specific sub-bands (RFI-avoidance notches), so a fixed FM
    control is not universally observed. Returns the first preference-ordered candidate the grid
    samples with >= ``_MIN_CONTROL_CH`` channels, or ("none", (nan, nan)) if none qualifies.
    """
    f = np.asarray(freqs, float)
    for name, lo, hi in _CONTROL_CANDIDATES:
        if ((f >= lo) & (f <= hi)).sum() >= _MIN_CONTROL_CH:
            return name, (lo, hi)
    return "none", (float("nan"), float("nan"))


def available_lines(
    freqs: np.ndarray, *, lines_mhz=UEM_LINES_MHZ, half_width: float = 1.0
) -> tuple[float, ...]:
    """Which narrowband UEM lines this grid samples with both a core AND flanking channels.

    A line is usable for `line_vs_adjacent` only if the instrument samples channels within
    +/-``half_width`` MHz of it AND in its flanks; stations that notch a line's region simply drop
    it. Any line falling inside the GRAVES radar band (143.05 MHz, not intrinsic Starlink) is
    excluded outright. Config-stability = this set is constant across epochs.
    """
    f = np.asarray(freqs, float)
    out = []
    for line in lines_mhz:
        if GRAVES_MHZ[0] <= line <= GRAVES_MHZ[1]:  # reflected GRAVES radar, not intrinsic Starlink
            continue
        core = (np.abs(f - line) <= half_width) & (f > 0)
        flank = (np.abs(f - line) > half_width) & (np.abs(f - line) <= 3 * half_width)
        if core.any() and flank.any():
            out.append(line)
    return tuple(out)


def band_differential(
    data: np.ndarray,
    freqs: np.ndarray,
    *,
    uem_band: tuple[float, float] = UEM_BAND_MHZ,
    control_band: tuple[float, float] | None = FM_CONTROL_MHZ,
) -> float:
    """UEM-band level MINUS a clean control-band level (log-power), the gain- and burst-robust metric.

    Both station gain drift and broadband solar bursts move the whole spectrum together
    (common-mode), so they cancel in the UEM-vs-control difference; only band-specific satellite RFI
    survives. ``control_band=None`` picks the best-sampled control for this grid via
    `pick_control_band` (real stations often notch FM). NaN if either band is out of coverage.
    """
    if control_band is None:
        _, control_band = pick_control_band(freqs)
    uem = occupancy_metric(data, freqs, uem_band)
    ctrl = occupancy_metric(data, freqs, control_band)
    return uem - ctrl if np.isfinite(uem) and np.isfinite(ctrl) else float("nan")


def line_vs_adjacent(
    data: np.ndarray, freqs: np.ndarray, *, lines_mhz=UEM_LINES_MHZ, half_width: float = 1.0
) -> float:
    """Excess level at the narrowband Starlink UEM lines over their adjacent (clean) channels.

    The cleanest Starlink-specific attribution: a real UEM line at 125/135/150/175 MHz sits ABOVE the
    band around it, whereas station gain / broadband RFI raise line and neighbourhood together. For
    each line we take (level within +/-``half_width`` MHz) minus (level in the flanking channels)
    and average. An earlier draft carried 137.05 MHz here (a transcription of GRAVES 143.05);
    the current line list is :data:`UEM_LINES_MHZ`, and the low lines are reported, not trusted
    alone. NaN if none of the lines is in coverage.
    """
    f = np.asarray(freqs, float)
    d = np.asarray(data, float)
    excesses = []
    for line in lines_mhz:
        core = (np.abs(f - line) <= half_width) & (f > 0)
        flank = (np.abs(f - line) > half_width) & (np.abs(f - line) <= 3 * half_width)
        if core.any() and flank.any():
            excesses.append(
                float(np.median(np.median(d[core], axis=1)))
                - float(np.median(np.median(d[flank], axis=1)))
            )
    return float(np.mean(excesses)) if excesses else float("nan")


def trend_fit(x: np.ndarray, y: np.ndarray) -> dict:
    """Robust (Theil--Sen) trend of ``y`` vs ``x`` + Kendall-tau significance and the fractional change.

    Theil--Sen is outlier-robust (a bad month does not swing it); Kendall tau gives a
    distribution-free monotonic-trend p-value. Returns the slope (per unit x), its p, and the total
    change over the span. NaN-safe.
    """
    from scipy import stats as _stats

    x = np.asarray(x, float)
    y = np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if x.size < 5 or np.ptp(x) == 0:
        return {
            "slope": float("nan"),
            "p_value": float("nan"),
            "n": int(x.size),
            "total_change": float("nan"),
        }
    slope, intercept, lo, hi = _stats.theilslopes(y, x)
    tau = _stats.kendalltau(x, y)
    return {
        "slope": float(slope),
        "slope_lo": float(lo),
        "slope_hi": float(hi),
        "p_value": float(tau.pvalue),
        "kendall_tau": float(tau.statistic),
        "n": int(x.size),
        "total_change": float(slope * (x.max() - x.min())),
    }


def step_analysis(
    years: np.ndarray, values: np.ndarray, *, n_steps: int = 2, min_seg: int = 6
) -> dict:
    """Is a series a trend or a staircase? Find the largest single-interval jumps and test within.

    The round-10 referee showed HUMAIN's fitted "+0.45/yr rise" is the difference between two
    plateaus separated by single-month discontinuities larger than the whole fitted change --- a
    signature of receiver reconfiguration, not of a growing constellation. This locates the
    ``n_steps`` largest |month-to-month jumps|, splits the series into regimes at them, and reports
    each regime's median and Theil--Sen slope/p: a genuine trend rises WITHIN regimes; a staircase
    does not. Segments shorter than ``min_seg`` months are reported but their fits are NaN.
    """
    yr = np.asarray(years, float)
    v = np.asarray(values, float)
    good = np.isfinite(yr) & np.isfinite(v)
    yr, v = yr[good], v[good]
    order = np.argsort(yr)
    yr, v = yr[order], v[order]
    if yr.size < 2 * min_seg:
        return {"steps": [], "segments": []}
    jumps = np.abs(np.diff(v))
    idx = np.argsort(jumps)[::-1][:n_steps]
    idx = np.sort(idx)
    steps = [
        {
            "at_year": round(float(yr[i + 1]), 3),
            "delta": round(float(v[i + 1] - v[i]), 2),
        }
        for i in idx
    ]
    bounds = [0, *[int(i) + 1 for i in idx], yr.size]
    segments = []
    for a, b in zip(bounds[:-1], bounds[1:], strict=False):
        seg_yr, seg_v = yr[a:b], v[a:b]
        fit = trend_fit(seg_yr, seg_v) if seg_yr.size >= min_seg else {}
        segments.append(
            {
                "start": round(float(seg_yr[0]), 3),
                "end": round(float(seg_yr[-1]), 3),
                "n": int(seg_yr.size),
                "median": round(float(np.median(seg_v)), 2),
                "slope": round(fit["slope"], 3) if fit.get("slope") is not None else None,
                "p_value": round(fit["p_value"], 3) if fit.get("p_value") is not None else None,
            }
        )
    return {"steps": steps, "segments": segments}


def coherence_power(
    per_station: dict,
    *,
    amps: tuple[float, ...] = (0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0),
    n_trials: int = 200,
    seed: int = 0,
) -> list[dict]:
    """Power of the paper's own coherence criterion against an injected COMMON Starlink-shaped rise.

    "No coherent rise" divides by sensitivity (the frblens lesson). For each amplitude ``A`` (total
    log-units at the final Starlink count): detrend each line-sampling station's committed series,
    circularly shift its residuals (preserving the real autocorrelation and step structure), add
    ``A``-scaled Starlink-shaped ramps to BOTH stations, and apply the same rule
    ``summarize_stations`` uses --- both stations significant (Kendall p<0.05), both rising, pooled
    Starlink correlation > 0.3. Reports P(coherent_rise) per amplitude. Deterministic.
    """
    rng = np.random.default_rng(seed)
    series = []
    for d in per_station.values():
        if not d.get("stable_lines"):
            continue
        yr = np.asarray(d.get("years", []), float)
        lx = np.asarray(d.get("line_excess", []), float)
        g = np.isfinite(yr) & np.isfinite(lx)
        yr, lx = yr[g], lx[g]
        if yr.size < 10:
            continue
        fit = trend_fit(yr, lx)
        resid = lx - fit["slope"] * (yr - yr.mean())
        resid = resid - np.median(resid)
        series.append((yr, resid))
    if len(series) < 2:
        return []
    smax = float(starlink_count(2026.5))
    out = []
    for a in amps:
        hits = 0
        for _ in range(n_trials):
            slopes, ps = [], []
            pooled_y: list[float] = []
            pooled_lx: list[float] = []
            for yr, resid in series:
                shift = int(rng.integers(0, resid.size))
                inj = np.roll(resid, shift) + a * starlink_count(yr) / smax
                fit = trend_fit(yr, inj)
                slopes.append(fit["slope"])
                ps.append(fit["p_value"])
                pooled_y.extend(yr)
                pooled_lx.extend(inj)
            py = np.asarray(pooled_y)
            plx = np.asarray(pooled_lx)
            pooled_corr = float(np.corrcoef(starlink_count(py), plx)[0, 1])
            sig_pos = [s for s, p in zip(slopes, ps, strict=True) if p < 0.05 and s > 0]
            sig_neg = [s for s, p in zip(slopes, ps, strict=True) if p < 0.05 and s < 0]
            if len(sig_pos) >= 2 and not sig_neg and pooled_corr > 0.3:
                hits += 1
        out.append({"amp": a, "power": round(hits / n_trials, 3)})
    return out


def equal_n_pooled(per_station: dict, *, n_draws: int = 500, seed: int = 0) -> dict:
    """Does the pooled slope survive equalising the station month counts?

    The abstract called the positive pooled slope "a sample-size artifact of whichever station
    dominates". If that were a COUNT artifact, subsampling the larger station to the smaller's
    month count would kill it; measured, it does not --- the pooled slope is dominated by one
    station's AMPLITUDE, which is a different (and correct) reason to distrust it.
    """
    rng = np.random.default_rng(seed)
    keyed = [
        (np.asarray(d["years"], float), np.asarray(d["line_excess"], float))
        for d in per_station.values()
        if d.get("stable_lines")
    ]
    keyed = [(y[np.isfinite(l_)], l_[np.isfinite(l_)]) for y, l_ in keyed]
    keyed = [(y, l_) for y, l_ in keyed if y.size >= 10]
    if len(keyed) < 2:
        return {}
    n_min = min(y.size for y, _ in keyed)
    slopes = []
    for _ in range(n_draws):
        ys, ls = [], []
        for y, l_ in keyed:
            if y.size > n_min:
                pick = rng.choice(y.size, n_min, replace=False)
                ys.append(y[pick])
                ls.append(l_[pick])
            else:
                ys.append(y)
                ls.append(l_)
        fit = trend_fit(np.concatenate(ys), np.concatenate(ls))
        slopes.append(fit["slope"])
    arr = np.asarray(slopes, float)
    return {
        "n_min": int(n_min),
        "n_draws": n_draws,
        "pooled_slope_median": round(float(np.median(arr)), 4),
        "pooled_slope_p05": round(float(np.percentile(arr, 5)), 4),
        "pooled_slope_p95": round(float(np.percentile(arr, 95)), 4),
    }


def regressor_shape_comparison(years: np.ndarray, values: np.ndarray) -> dict:
    """Correlation of a series with three shapes: the Starlink curve, a linear ramp, a 2019.4 step.

    The only quantitative discriminant available for a single-station candidate: a genuine
    constellation signal should prefer the Starlink shape over a generic drift or a single step.
    """
    yr = np.asarray(years, float)
    v = np.asarray(values, float)
    g = np.isfinite(yr) & np.isfinite(v)
    yr, v = yr[g], v[g]
    if yr.size < 10:
        return {}

    def _c(shape: np.ndarray) -> float:
        return round(float(np.corrcoef(shape, v)[0, 1]), 3)

    return {
        "corr_starlink_shape": _c(starlink_count(yr)),
        "corr_linear_ramp": _c(yr),
        "corr_step_2019": _c((yr >= 2019.4).astype(float)),
    }


def window_matched_trend(per_station: dict, target: str, reference: str) -> dict:
    """Fit ``target``'s series restricted to ``reference``'s retained window (span-confound check)."""
    t = per_station.get(target, {})
    r = per_station.get(reference, {})
    ty = np.asarray(t.get("years", []), float)
    tl = np.asarray(t.get("line_excess", []), float)
    ry = np.asarray(r.get("years", []), float)
    if ty.size < 10 or ry.size < 10:
        return {}
    lo, hi = float(np.nanmin(ry)), float(np.nanmax(ry))
    m = (ty >= lo) & (ty <= hi) & np.isfinite(tl)
    fit = trend_fit(ty[m], tl[m])
    star = starlink_count(ty[m])
    corr = float(np.corrcoef(star, tl[m])[0, 1]) if m.sum() > 3 else float("nan")
    return {
        "window": [round(lo, 3), round(hi, 3)],
        "n": int(m.sum()),
        "slope": round(fit["slope"], 4),
        "p_value": round(fit["p_value"], 5),
        "corr_starlink": round(corr, 3),
    }


# Public Starlink cumulative on-orbit count by year-start (planet4589 / J. McDowell statistics,
# rounded); the attribution regressor -- no Space-Track login needed. NOTE the interpolation is
# nonzero from 2019.1 (linear between the 2019.0 and 2019.5 anchors) and CLAMPS at the last
# anchor for months past 2026.0 -- both properties the paper states rather than papering over.
_STARLINK_BY_YEAR = {
    2015: 0,
    2019.0: 0,  # first operational batch launched 2019.4 -> zero at year-start
    2019.5: 120,
    2020: 1000,
    2021: 1900,
    2022: 3300,
    2023: 5000,
    2024: 6900,
    2025: 8500,
    2026: 10400,
}


def starlink_count(decimal_year) -> np.ndarray:
    """Cumulative on-orbit Starlink count at ``decimal_year`` (scalar or array), interpolated."""
    yrs = np.array(sorted(_STARLINK_BY_YEAR))
    cnt = np.array([_STARLINK_BY_YEAR[y] for y in yrs], float)
    return np.interp(np.asarray(decimal_year, float), yrs, cnt)


def synthetic_month_stack(
    *,
    n_months: int = 168,  # 2012-01 .. 2025-12
    n_freq: int = 500,  # over 45-450 MHz -> ~0.8 MHz channels, fine enough to resolve the UEM lines
    n_time: int = 120,
    uem_rise: float = 3.0,  # total UEM-band level rise (log units) tracking Starlink, post-2019
    line_rise: float = 4.0,  # total narrowband-LINE excess rise tracking Starlink (primary metric)
    gain_sigma: float = 2.0,  # per-month station-gain drift (common-mode -- must cancel)
    burst_frac: float = 0.3,  # fraction of months with a broadband solar burst (must not bias)
    flank_rise: float = 0.0,  # local RFI accruing in the lines' FLANKING channels only (see below)
    seed: int = 0,
) -> dict:
    """Synthetic per-month spectra with a KNOWN UEM trend + gain drift + solar bursts, for recovery.

    The recover-a-known: both `band_differential` (secondary) AND `line_vs_adjacent` (PRIMARY) must
    recover the injected Starlink-shaped rise while staying immune to common-mode gain drift and
    broadband bursts. A broadband UEM-band level rise is injected for the differential; NARROWBAND
    excesses at the UEM lines are injected on top for the line test. The FM control carries no
    injected trend. Returns per-month (data, freqs) plus the decimal-year and truth.

    ``flank_rise`` injects the systematic the two arms above CANNOT probe, because the
    line-vs-adjacent difference cancels anything common to core and flank algebraically: local
    terrestrial RFI accruing in the lines' *flanking* channels only, as a linear-in-time ramp
    (local RFI has no reason to track the Starlink deployment curve). This is the mechanism the
    paper names as decisive for the sign disagreement -- a station whose flanks fill in faster
    than its lines registers a FALLING excess under a rising true signal. Before this arm existed
    the validation could not fail in the regime that limits the study.
    """
    rng = np.random.default_rng(seed)
    freqs = np.linspace(45.0, 450.0, n_freq)
    years = 2012.0 + np.arange(n_months) / 12.0
    star = starlink_count(years)
    star_frac = star / star.max() if star.max() > 0 else star  # 0..1, the injected UEM shape
    months = []
    uem_sel = (freqs >= UEM_BAND_MHZ[0]) & (freqs <= UEM_BAND_MHZ[1])
    line_sel = np.zeros(n_freq, bool)  # narrowband cores at the intrinsic UEM lines
    for line in UEM_LINES_MHZ:
        line_sel |= np.abs(freqs - line) <= 1.0
    for i in range(n_months):
        gain = rng.normal(0.0, gain_sigma)  # common-mode station-gain drift this month
        data = rng.normal(10.0, 0.5, (n_freq, n_time)) + gain
        data[uem_sel] += uem_rise * star_frac[i]  # band-specific UEM level (post-2019)
        data[line_sel] += line_rise * star_frac[i]  # narrowband line excess on top (primary metric)
        if flank_rise:
            # local RFI in the flanks only: linear ramp, NOT Starlink-shaped
            flank_sel = np.zeros(n_freq, bool)
            for line in UEM_LINES_MHZ:
                flank_sel |= (np.abs(freqs - line) > 1.0) & (np.abs(freqs - line) <= 3.0)
            data[flank_sel] += flank_rise * (i / max(n_months - 1, 1))
        if rng.random() < burst_frac:  # a broadband solar burst in a few time columns
            t0 = rng.integers(0, n_time - 10)
            data[:, t0 : t0 + 8] += rng.uniform(8, 20)
        months.append({"data": data, "freqs": freqs})
    return {"months": months, "years": years, "star_frac": star_frac, "uem_rise": uem_rise}


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: synthetic recover-a-known (differential recovers the UEM trend, burst/gain-immune)."""

    if offline:
        metrics: dict = _synthetic_metrics()
    else:  # pragma: no cover - real leg streams the e-Callisto archive
        metrics = _real_trend(out)

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "rfitrend_metrics.json")
    _figure(metrics, op / "papers" / "rfitrend" / "figures")
    _write_macros(metrics, op / "papers" / "rfitrend" / "generated" / "macros.tex")
    return metrics


def _synthetic_metrics() -> dict:
    """Synthetic recover-a-known metrics (differential recovers the injected UEM trend, gain/burst-immune).

    Deterministic and offline. Emitted as the ``rfSyn*`` macros in BOTH the offline and the real
    build so the paper's validation section is always populated.
    """
    s = synthetic_month_stack()
    diff = np.array([band_differential(m["data"], m["freqs"]) for m in s["months"]])
    # a null control: two bands BOTH outside the UEM window -> the differential must be flat
    # (no injected trend), proving the metric doesn't manufacture trends from gain drift
    fm_self = np.array(
        [
            occupancy_metric(m["data"], m["freqs"], FM_CONTROL_MHZ)
            - occupancy_metric(m["data"], m["freqs"], (300.0, 400.0))
            for m in s["months"]
        ]
    )
    # PRIMARY metric: the narrowband line-vs-adjacent excess must ALSO recover the injected trend
    line = np.array([line_vs_adjacent(m["data"], m["freqs"]) for m in s["months"]])
    tr = trend_fit(s["years"], diff)
    line_tr = trend_fit(s["years"], line)
    ctrl_tr = trend_fit(s["years"], fm_self)
    corr = float(np.corrcoef(s["star_frac"], diff)[0, 1])
    line_corr = float(np.corrcoef(s["star_frac"], line)[0, 1])
    # The arm that can fail: same injected line rise, plus local RFI accruing in the flanking
    # channels only. The line-vs-adjacent difference cannot cancel it (it is not common to core
    # and flank), so the recovered slope is BIASED -- with a large enough flank ramp, past zero.
    sf = synthetic_month_stack(flank_rise=6.0)
    flank_line = np.array([line_vs_adjacent(m["data"], m["freqs"]) for m in sf["months"]])
    flank_tr = trend_fit(sf["years"], flank_line)
    # the KNOWNS this recover-a-known previously never stated: the Theil-Sen slope of the
    # injected series itself, so the recovery is scored as an amplitude ratio, not a detection
    truth_line = trend_fit(s["years"], 4.0 * s["star_frac"])["slope"]
    truth_diff = trend_fit(s["years"], 3.0 * s["star_frac"])["slope"]
    # the flank sweep: where does the sign cross, and which flank_rise reproduces ALMATY?
    sweep = [{"flank_rise": 0.0, "slope": round(line_tr["slope"], 4), "corr": round(line_corr, 3)}]
    for fr in (2.0, 3.0, 3.5, 5.0, 6.0):
        st = synthetic_month_stack(flank_rise=fr)
        lx = np.array([line_vs_adjacent(m["data"], m["freqs"]) for m in st["months"]])
        fitf = trend_fit(st["years"], lx)
        cf = float(np.corrcoef(st["star_frac"], lx)[0, 1])
        sweep.append({"flank_rise": fr, "slope": round(fitf["slope"], 4), "corr": round(cf, 3)})
    # burst immunity on the PRIMARY metric, measured, not asserted
    s00 = synthetic_month_stack(burst_frac=0.0)
    s90 = synthetic_month_stack(burst_frac=0.9)
    sl00 = trend_fit(
        s00["years"], np.array([line_vs_adjacent(m["data"], m["freqs"]) for m in s00["months"]])
    )["slope"]
    sl90 = trend_fit(
        s90["years"], np.array([line_vs_adjacent(m["data"], m["freqs"]) for m in s90["months"]])
    )["slope"]
    return {
        "source": "synthetic monthly stack (injected UEM trend + gain drift + solar bursts)",
        "is_real": False,
        "n_months": int(s["years"].size),
        "diff_slope_per_yr": round(tr["slope"], 4),
        "diff_trend_p": round(tr["p_value"], 5),
        "diff_total_change": round(tr["total_change"], 3),
        "line_excess_slope_per_yr": round(line_tr["slope"], 4),
        "line_excess_trend_p": round(line_tr["p_value"], 5),
        "line_corr_with_starlink": round(line_corr, 3),
        "control_slope_per_yr": round(ctrl_tr["slope"], 4),
        "control_trend_p": round(ctrl_tr["p_value"], 3),
        "corr_with_starlink": round(corr, 3),
        "recovered_uem_trend": bool(tr["p_value"] < 0.01 and corr > 0.8),
        "recovered_line_trend": bool(line_tr["p_value"] < 0.01 and line_corr > 0.8),
        "control_flat": bool(ctrl_tr["p_value"] > 0.05),
        "truth_line_slope_per_yr": round(truth_line, 4),
        "truth_diff_slope_per_yr": round(truth_diff, 4),
        "line_recovery_ratio": round(line_tr["slope"] / truth_line, 3),
        "diff_recovery_ratio": round(tr["slope"] / truth_diff, 3),
        "line_slope_burst_none": round(sl00, 4),
        "line_slope_burst_heavy": round(sl90, 4),
        # Flank-contamination arm: the recovered slope under flank RFI, and its bias against the
        # clean arm's slope. A sign flip here is the ALMATY mechanism reproduced end-to-end.
        "flank_line_slope_per_yr": round(flank_tr["slope"], 4),
        "flank_slope_bias_per_yr": round(flank_tr["slope"] - line_tr["slope"], 4),
        "flank_sign_flipped": bool(np.sign(flank_tr["slope"]) != np.sign(line_tr["slope"])),
        "flank_sweep": sweep,
    }


def _decimal_year(yyyymm: str) -> float:
    y, m = int(yyyymm[:4]), int(yyyymm[4:6])
    return y + (m - 0.5) / 12.0


def sample_month_metrics(
    station: str, year: int, month: int, *, days=(15, 10, 20, 5, 25), times=("1000", "1200", "0800")
):  # pragma: no cover - network (one e-Callisto FITS fetch)
    """Fetch one representative spectrum for a station-month and reduce it to the RFI metrics.

    Tries a few days/UTs until the station has data. Returns a dict with `line_excess` (primary,
    narrowband UEM lines over flanks), `differential` (UEM minus the auto-picked control),
    `uem_level`, `control_name`, `lines` (which UEM lines were sampled), or None if no spectrum for
    the month. The metrics are per-channel time-medians -> burst-immune; the differential/line tests
    cancel station gain.
    """
    from .solarbursts import fetch_ecallisto

    for day in days:
        for hhmm in times:
            date = f"{year:04d}{month:02d}{day:02d}"
            try:
                ds = fetch_ecallisto(station, date, hhmm)
            except Exception:
                continue
            data, freqs = ds["data"], np.asarray(ds["freqs"], float)
            cname, cband = pick_control_band(freqs)
            lines = available_lines(freqs)
            return {
                "station": station,
                "yyyymm": f"{year:04d}{month:02d}",
                "decimal_year": _decimal_year(f"{year:04d}{month:02d}"),
                "line_excess": line_vs_adjacent(data, freqs, lines_mhz=lines)
                if lines
                else float("nan"),
                "differential": band_differential(data, freqs, control_band=cband),
                "uem_level": occupancy_metric(data, freqs, UEM_BAND_MHZ),
                "control_name": cname,
                "lines": lines,
                "file": ds.get("file", ""),
            }
    return None


def _real_trend(
    out: str,
    *,
    stations=ECALLISTO_STATIONS,
    start_year: int = 2012,
    end_year: int = 2026,
) -> dict:  # pragma: no cover - network (e-Callisto archive + monthly sampling)
    """Real leg: monthly-sampled narrowband-UEM-line trend per station 2012-2026 + Starlink attribution.

    PRIMARY metric = `line_vs_adjacent` (narrowband UEM lines over adjacent clean channels): it is
    self-normalizing within the UEM band, so it cancels station gain AND survives the per-station
    RFI-avoidance notches that make FM unusable at many stations (HUMAIN notches FM entirely). The
    FM/low/high band `differential` is a station-adaptive cross-check. Config-stability is enforced
    by requiring a station keep the SAME sampled UEM lines across its retained months.
    """
    from collections import Counter

    per_station: dict[str, dict] = {}
    for st in stations:
        rows = []
        for yr in range(start_year, end_year + 1):
            for mo in range(1, 13):
                if yr == end_year and mo > 7:  # today is 2026-07
                    break
                r = sample_month_metrics(st, yr, mo)
                if r is not None:
                    rows.append(r)
        if not rows:
            per_station[st] = {"n_months": 0, "note": "no data"}
            continue
        # config-stability: keep the modal set of sampled UEM lines; drop months that differ
        line_sets = Counter(tuple(r["lines"]) for r in rows)
        stable_lines = max(line_sets, key=lambda k: line_sets[k])
        stable = [r for r in rows if tuple(r["lines"]) == stable_lines]
        yrs = np.array([r["decimal_year"] for r in stable])
        lx = np.array([r["line_excess"] for r in stable])
        df = np.array([r["differential"] for r in stable])
        star = starlink_count(yrs)
        lx_tr = trend_fit(yrs, lx)
        df_tr = trend_fit(yrs, df)
        good = np.isfinite(lx) & np.isfinite(star)
        corr = float(np.corrcoef(star[good], lx[good])[0, 1]) if good.sum() > 3 else float("nan")
        # Perez+2020 reproduction: UEM occupancy 2012-2013 vs 2018-2019 (pre-Starlink rise)
        ulev = np.array([r["uem_level"] for r in stable])
        pre = ulev[(yrs >= 2012) & (yrs < 2014)]
        post19 = ulev[(yrs >= 2018) & (yrs < 2020)]
        perez_ratio = (
            float(np.nanmedian(post19) - np.nanmedian(pre))
            if pre.size and post19.size
            else float("nan")
        )
        per_station[st] = {
            "n_months": len(stable),
            "n_months_raw": len(rows),
            "stable_lines": list(stable_lines),
            "control_name": stable[0]["control_name"],
            "line_excess_slope_per_yr": round(lx_tr["slope"], 4),
            "line_excess_p": round(lx_tr["p_value"], 5),
            "diff_slope_per_yr": round(df_tr["slope"], 4),
            "diff_p": round(df_tr["p_value"], 5),
            "corr_line_excess_starlink": round(corr, 3),
            "perez_2012_2019_change": round(perez_ratio, 3),
            "years": [round(float(y), 3) for y in yrs],
            "line_excess": [round(float(v), 3) for v in lx],
        }

    return summarize_stations(per_station)


def summarize_stations(per_station: dict) -> dict:
    """Pool the per-station line-excess series + compute the CROSS-STATION COHERENCE diagnostic.

    The load-bearing honest test. A real *global* Starlink UEM signal must make every station's
    UEM-line excess rise TOGETHER (the constellation is overhead everywhere). So beyond the pooled
    trend + Starlink correlation, we test whether the stations that significantly trend AGREE IN
    SIGN. If they disagree (one rises, one falls), the trends are per-station/local-RFI systematics,
    NOT a coherent megaconstellation signal -- the pooled slope is then a sample-size artifact of
    whichever station dominates, and no Starlink attribution is warranted. This verdict
    (`coherent_rise`) is emitted so the artifact carries the conclusion, not just the prose.
    """
    # pooled trend over every USABLE (finite line-excess) station-month -- exclude stations that
    # sample no UEM line (their line_excess is all-NaN and must not inflate the reported sample)
    yrs_l: list[float] = []
    lx_l: list[float] = []
    for d in per_station.values():
        if d.get("n_months", 0) >= 5 and d.get("stable_lines"):
            for yr, lx in zip(d.get("years", []), d.get("line_excess", []), strict=False):
                if np.isfinite(lx):
                    yrs_l.append(float(yr))
                    lx_l.append(float(lx))
    all_yrs, all_lx = np.array(yrs_l, float), np.array(lx_l, float)
    pooled = trend_fit(all_yrs, all_lx) if all_yrs.size >= 5 else {"slope": float("nan")}
    pstar = starlink_count(all_yrs) if all_yrs.size else np.array([])
    g = np.isfinite(all_lx) & np.isfinite(pstar) if all_yrs.size else np.array([], bool)
    pooled_corr = (
        float(np.corrcoef(pstar[g], all_lx[g])[0, 1])
        if all_yrs.size and g.sum() > 3
        else float("nan")
    )
    # cross-station coherence: signs of the SIGNIFICANT per-station line-excess slopes
    sig_slopes = [
        d["line_excess_slope_per_yr"]
        for d in per_station.values()
        if d.get("stable_lines")
        and np.isfinite(d.get("line_excess_slope_per_yr", float("nan")))
        and d.get("line_excess_p", 1.0) < 0.05
    ]
    n_pos = sum(1 for s in sig_slopes if s > 0)
    n_neg = sum(1 for s in sig_slopes if s < 0)
    signs_agree = len(sig_slopes) >= 2 and (n_pos == 0 or n_neg == 0)
    # a coherent rise needs >=2 stations that trend, all rising, AND pooled correlation with Starlink
    coherent_rise = bool(
        signs_agree and n_pos >= 2 and np.isfinite(pooled_corr) and pooled_corr > 0.3
    )
    n_lines = sum(1 for d in per_station.values() if d.get("stable_lines"))
    return {
        "source": "e-Callisto archive, monthly-sampled narrowband-UEM-line occupancy 2012-2026",
        "is_real": True,
        "n_stations": sum(1 for d in per_station.values() if d.get("n_months", 0) >= 5),
        "n_stations_with_lines": n_lines,
        "n_months": int(all_yrs.size),  # usable (finite line-excess) station-months only
        "line_excess_slope_per_yr": round(float(pooled["slope"]), 4),
        "line_excess_trend_p": round(float(pooled.get("p_value", float("nan"))), 5),
        "corr_with_starlink": round(pooled_corr, 3),
        "n_significant_stations": len(sig_slopes),
        "n_rising": n_pos,
        "n_falling": n_neg,
        "cross_station_signs_agree": signs_agree,
        "coherent_rise": coherent_rise,
        "per_station": per_station,
    }


def derived_real_analyses(per_station: dict) -> dict:
    """Every derived diagnostic computable OFFLINE from the committed per-station arrays.

    This is what lets the round-10 findings be fixed without re-streaming the archive: step
    analysis and shape comparison per line station, Theil--Sen slope intervals, the coherence
    power curve, the equal-n pooled test, and the window-matched mirror check.
    """
    out: dict = {"stations": {}}
    line_sts = []
    for st, d in per_station.items():
        if not (d.get("stable_lines") and d.get("years")):
            continue
        line_sts.append(st)
        yr = np.asarray(d["years"], float)
        lx = np.asarray(d["line_excess"], float)
        fit = trend_fit(yr, lx)
        out["stations"][st] = {
            "slope_ci95": [round(fit["slope_lo"], 4), round(fit["slope_hi"], 4)],
            "steps": step_analysis(yr, lx, n_steps=3),
            "shapes": regressor_shape_comparison(yr, lx),
        }
    out["coherence_power"] = coherence_power(per_station)
    out["equal_n_pooled"] = equal_n_pooled(per_station)
    if len(line_sts) == 2:
        a, b = sorted(line_sts, key=lambda s: len(per_station[s].get("years", [])), reverse=True)
        out["window_matched"] = {
            "target": a,
            "reference": b,
            **window_matched_trend(per_station, a, b),
        }
    return out


def _figure(m: dict, out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.9))

    # left: synthetic recover-a-known (always)
    s = synthetic_month_stack()
    diff = np.array([band_differential(mm["data"], mm["freqs"]) for mm in s["months"]])
    ax1.plot(s["years"], diff, ".", ms=3, color="C0", label="recovered differential")
    ax1.plot(s["years"], s["star_frac"] * s["uem_rise"], "-", color="C3", label="injected")
    ax1.set(xlabel="year", ylabel="UEM$-$control level", title="Synthetic recover-a-known")
    ax1.legend(fontsize=8)

    # right: real per-station line-excess vs year (the incoherence), else synthetic scatter
    if m.get("is_real") and m.get("per_station"):
        for i, (st, d) in enumerate(m["per_station"].items()):
            if d.get("stable_lines") and d.get("line_excess"):
                yr = np.array(d["years"], float)
                lx = np.array(d["line_excess"], float)
                sl = d.get("line_excess_slope_per_yr")
                ax2.plot(yr, lx, ".", ms=3, color=f"C{i}", alpha=0.5, label=f"{st} ({sl:+.2f}/yr)")
        ax2.axvline(2019.4, ls=":", color="grey", lw=1)  # Starlink deployment start
        ax2.set(
            xlabel="year",
            ylabel="UEM-line excess",
            title=f"Real: signs {'agree' if m.get('cross_station_signs_agree') else 'DISAGREE'}",
        )
        ax2.legend(fontsize=7)
    else:
        ax2.plot(s["star_frac"], diff, ".", ms=3, color="C0")
        ax2.set(
            xlabel="Starlink count (norm.)",
            ylabel="differential",
            title=f"r={m.get('corr_with_starlink')}",
        )
    fig.tight_layout()
    fig.savefig(out / "rfitrend.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    """Emit both namespaces: rfSyn* ALWAYS live (synthetic recovery), rfReal* live only in a real run.

    The synthetic validation numbers are recomputed here (deterministic, offline) so the paper's
    recover-a-known section is populated even in the final real-data build; the real leg fills
    rfReal*.
    """

    def g(src: dict, key: str) -> str:
        v = src.get(key)
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    def gp(src: dict, key: str) -> str:
        """Format a p-value as a LaTeX math body (used inside $...$): never prints a bare '0.0'."""
        v = src.get(key)
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "--"
        v = float(v)
        if v < 1e-5:  # round(p,5) collapsed it to 0.0 -> report as an upper bound
            return r"<10^{-5}"
        if v < 1e-3:
            mant, exp = f"{v:.0e}".split("e")
            return rf"{mant}\times10^{{{int(exp)}}}"
        return f"{v:.3f}"

    syn = _synthetic_metrics()
    real = m if m.get("is_real") else {}
    lines = [
        "% Auto-generated by jansky_research.rfitrend._write_macros -- do not edit.",
        "% rfSyn* = synthetic recover-a-known (always live); rfReal* = real e-Callisto leg",
        "% (live only after a real run; '--' placeholders otherwise). *P macros are math bodies.",
        rf"\newcommand{{\rfSource}}{{{m['source']}}}",
        rf"\newcommand{{\rfSynNMonths}}{{{g(syn, 'n_months')}}}",
        rf"\newcommand{{\rfRealNMonths}}{{{g(real, 'n_months')}}}",
        rf"\newcommand{{\rfRealNStations}}{{{g(real, 'n_stations')}}}",
        rf"\newcommand{{\rfRealNLineStations}}{{{g(real, 'n_stations_with_lines')}}}",
        rf"\newcommand{{\rfRealNRising}}{{{g(real, 'n_rising')}}}",
        rf"\newcommand{{\rfRealNFalling}}{{{g(real, 'n_falling')}}}",
        rf"\newcommand{{\rfRealCoherent}}{{{'yes' if real.get('coherent_rise') else 'no'}}}",
        # syn validates the DIFFERENTIAL; real headlines the LINE-EXCESS (the primary metric)
        rf"\newcommand{{\rfSynSlope}}{{{g(syn, 'diff_slope_per_yr')}}}",
        rf"\newcommand{{\rfSynTrendP}}{{{gp(syn, 'diff_trend_p')}}}",
        rf"\newcommand{{\rfSynTotalChange}}{{{g(syn, 'diff_total_change')}}}",
        rf"\newcommand{{\rfSynCtrlSlope}}{{{g(syn, 'control_slope_per_yr')}}}",
        rf"\newcommand{{\rfSynCorr}}{{{g(syn, 'corr_with_starlink')}}}",
        rf"\newcommand{{\rfSynLineExcessSlope}}{{{g(syn, 'line_excess_slope_per_yr')}}}",
        rf"\newcommand{{\rfSynFlankSlope}}{{{g(syn, 'flank_line_slope_per_yr')}}}",
        rf"\newcommand{{\rfSynFlankBias}}{{{g(syn, 'flank_slope_bias_per_yr')}}}",
        rf"\newcommand{{\rfSynLineCorr}}{{{g(syn, 'line_corr_with_starlink')}}}",
        rf"\newcommand{{\rfRealSlope}}{{{g(real, 'line_excess_slope_per_yr')}}}",
        rf"\newcommand{{\rfRealTrendP}}{{{gp(real, 'line_excess_trend_p')}}}",
        rf"\newcommand{{\rfRealLineExcessSlope}}{{{g(real, 'line_excess_slope_per_yr')}}}",
        rf"\newcommand{{\rfRealCorr}}{{{g(real, 'corr_with_starlink')}}}",
    ]
    # per-station specifics for the results table (letters-only macro names: \rfReal<NAME><field>).
    # ALWAYS emit for the canonical station set so the paper's macros are defined even in the offline
    # build (real per_station is empty offline -> every field is the "--" placeholder).
    per_station = real.get("per_station") or {}
    for st in ECALLISTO_STATIONS:
        d = per_station.get(st, {})
        name = "".join(ch for ch in st if ch.isalpha())
        lines.append(rf"\newcommand{{\rfReal{name}Slope}}{{{g(d, 'line_excess_slope_per_yr')}}}")
        lines.append(rf"\newcommand{{\rfReal{name}P}}{{{gp(d, 'line_excess_p')}}}")
        lines.append(
            rf"\newcommand{{\rfReal{name}StarCorr}}{{{g(d, 'corr_line_excess_starlink')}}}"
        )
        lines.append(rf"\newcommand{{\rfReal{name}Perez}}{{{g(d, 'perez_2012_2019_change')}}}")
        lines.append(rf"\newcommand{{\rfReal{name}NMonths}}{{{g(d, 'n_months')}}}")
        # The pre-screening count, so the config-stability screen's attrition is visible. It is
        # not uniform and it is not incidental: ALMATY loses 28% of its months and ALMATY is the
        # station carrying the sign disagreement the paper's null rests on.
        lines.append(rf"\newcommand{{\rfReal{name}NMonthsRaw}}{{{g(d, 'n_months_raw')}}}")
        sl = d.get("stable_lines") or []
        lines.append(
            # "none" for an empty result, NOT "--". A station with no stable lines is a
            # measurement; "--" is this repo's marker for a value the results JSON did not
            # have, and the same table row uses \nodata for genuinely inapplicable cells. All
            # three rendered identically before, so a reader could not tell "we looked and
            # found none" from "this number is missing".
            rf"\newcommand{{\rfReal{name}Lines}}{{{', '.join(f'{x:g}' for x in sl) or 'none'}}}"
        )
        dstat = ((real.get("derived") or {}).get("stations") or {}).get(st, {})
        ci = dstat.get("slope_ci95") or [None, None]
        lines.append(rf"\newcommand{{\rfReal{name}SlopeLo}}{{{g({'v': ci[0]}, 'v')}}}")
        lines.append(rf"\newcommand{{\rfReal{name}SlopeHi}}{{{g({'v': ci[1]}, 'v')}}}")

    # the derived diagnostics (round-10): steps, shapes, power, equal-n, window mirror
    der = real.get("derived") or {}
    hum = (der.get("stations") or {}).get("HUMAIN", {})
    steps = (hum.get("steps") or {}).get("steps") or [{}, {}]
    segs = (hum.get("steps") or {}).get("segments") or []
    seg_slopes = [(s.get("slope"), s.get("p_value")) for s in segs if s.get("slope") is not None]
    max_seg = max(seg_slopes, key=lambda t: t[0]) if seg_slopes else (None, None)
    shapes = hum.get("shapes") or {}
    power = {row["amp"]: row["power"] for row in der.get("coherence_power") or []}
    amp90 = next(
        (row["amp"] for row in der.get("coherence_power") or [] if row["power"] >= 0.9), None
    )
    eq = der.get("equal_n_pooled") or {}
    win = der.get("window_matched") or {}
    sweep = syn.get("flank_sweep") or []
    cross = next((row["flank_rise"] for row in sweep if row["slope"] < 0), None)
    at5: dict = next((row for row in sweep if row["flank_rise"] == 5.0), {})
    lines += [
        rf"\newcommand{{\rfRealHumStepOneYr}}{{{g(steps[0] if steps else {}, 'at_year')}}}",
        rf"\newcommand{{\rfRealHumStepOneDelta}}{{{g(steps[0] if steps else {}, 'delta')}}}",
        rf"\newcommand{{\rfRealHumStepTwoYr}}{{{g(steps[1] if len(steps) > 1 else {}, 'at_year')}}}",
        rf"\newcommand{{\rfRealHumStepTwoDelta}}{{{g(steps[1] if len(steps) > 1 else {}, 'delta')}}}",
        rf"\newcommand{{\rfRealHumStepThreeYr}}{{{g(steps[2] if len(steps) > 2 else {}, 'at_year')}}}",
        rf"\newcommand{{\rfRealHumStepThreeDelta}}{{{g(steps[2] if len(steps) > 2 else {}, 'delta')}}}",
        rf"\newcommand{{\rfRealHumMaxRegimeSlope}}{{{g({'v': max_seg[0]}, 'v')}}}",
        rf"\newcommand{{\rfRealHumMaxRegimeP}}{{{g({'v': max_seg[1]}, 'v')}}}",
        rf"\newcommand{{\rfRealHumShapeStar}}{{{g(shapes, 'corr_starlink_shape')}}}",
        rf"\newcommand{{\rfRealHumShapeRamp}}{{{g(shapes, 'corr_linear_ramp')}}}",
        rf"\newcommand{{\rfRealHumShapeStep}}{{{g(shapes, 'corr_step_2019')}}}",
        rf"\newcommand{{\rfRealPowerAtSix}}{{{g({'v': power.get(6.0)}, 'v')}}}",
        rf"\newcommand{{\rfRealPowerAtTwenty}}{{{g({'v': power.get(20.0)}, 'v')}}}",
        rf"\newcommand{{\rfRealAmpNinety}}{{{g({'v': amp90}, 'v')}}}",
        rf"\newcommand{{\rfRealEqualNSlope}}{{{g(eq, 'pooled_slope_median')}}}",
        rf"\newcommand{{\rfRealEqualNLo}}{{{g(eq, 'pooled_slope_p05')}}}",
        rf"\newcommand{{\rfRealEqualNHi}}{{{g(eq, 'pooled_slope_p95')}}}",
        rf"\newcommand{{\rfRealWinSlope}}{{{g(win, 'slope')}}}",
        rf"\newcommand{{\rfRealWinCorr}}{{{g(win, 'corr_starlink')}}}",
        rf"\newcommand{{\rfSynTruthLineSlope}}{{{g(syn, 'truth_line_slope_per_yr')}}}",
        rf"\newcommand{{\rfSynLineRecovery}}{{{g(syn, 'line_recovery_ratio')}}}",
        rf"\newcommand{{\rfSynBurstNone}}{{{g(syn, 'line_slope_burst_none')}}}",
        rf"\newcommand{{\rfSynBurstHeavy}}{{{g(syn, 'line_slope_burst_heavy')}}}",
        rf"\newcommand{{\rfSynFlankCross}}{{{g({'v': cross}, 'v')}}}",
        rf"\newcommand{{\rfSynFlankAtFiveSlope}}{{{g(at5, 'slope')}}}",
        rf"\newcommand{{\rfSynFlankAtFiveCorr}}{{{g(at5, 'corr')}}}",
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

    p = argparse.ArgumentParser(description="e-Callisto megaconstellation RFI-trend census.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
