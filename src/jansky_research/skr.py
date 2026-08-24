"""Cassini SKR occurrence census + Saturn-proximity duty-cycle law (plan 60).

Saturn Kilometric Radiation (SKR) is Saturn's dominant auroral radio emission (~100--400 kHz,
cyclotron-maser). This slice ports the merged `junodam` Jovian-DAM occurrence-census pattern to
Saturn: a background+k-sigma detection over the Cassini/RPWS 60-s key-parameter flux, folded
against the Cassini--Saturn range from JPL Horizons, to test whether SKR **occurrence/duty-cycle
rises as Cassini approaches Saturn** (a proximity law, by analogy to the junodam ~180x DAM
result).

**Novelty, honestly scoped (GATE-0, 2026-07-08).** The SKR dual-period record (~10.6/10.8 h
north/south) is ALREADY published end-to-end --- Fischer+2015 (Icarus 254, 72; through early
2013), Gurnett+2016 (2012--2015), Provan+2019 (2016 to end of mission, N~10.79/S~10.68 h). So the
period tracking here is **validation only** (a Lomb-Scargle re-derivation of the ~10.7 h rotation
period, anchoring the pipeline to the literature), NOT a new result. The unclaimed angle is the
**occurrence/duty-cycle-vs-Saturn-distance proximity law**, which no one has run.

**The central caveat, stated from the outset:** a proximity--occurrence trend is confounded with
detection sensitivity (closer = stronger signal = more above-threshold bins) AND with SKR beaming
(sub-spacecraft latitude coverage varies with orbit). This is a *visibility* law unless corrected;
`magnetic_latitude_weight` is the stated-model-dependence attempt to separate intrinsic occurrence
from viewing geometry, and gets its own GATE-2 scrutiny. Like junodam, the headline is framed as
proximity-DOMINATED detection, not intrinsic emission.

Data: PDS-PPI `CO-V/E/J/S/SS-RPWS-4-SUMM-KEY60S-V1.0` (volume CORPWS_9002, whole mission), per-day
`RPWS_KEY__<YYYY><DDD>_<n>.TAB` fixed-length ASCII (ROW_BYTES 1175): a 1-row frequency table (115
channels; the first 73 are the electric antenna, 1 Hz--16 MHz at 0.1-decade spacing) then ~1382
1-minute rows of ELECTRIC_SPECTRAL_DENSITIES (73 ch, V^2/m^2/Hz). No pre-integrated SKR flux ---
we band-integrate the electric channels over the SKR band ourselves. Reuse: the `junodam`
detection/occurrence/Horizons pattern; astropy Lomb-Scargle for the period anchor.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "SKR_BAND_HZ",
    "parse_scet_to_jd",
    "read_key_params",
    "band_integrated_flux",
    "detect_skr",
    "dual_period_ls",
    "proximity_duty_cycle",
    "magnetic_latitude_weight",
    "synthetic_skr",
    "run",
]

SKR_BAND_HZ = (1.0e5, 1.2e6)  # SKR band: ~100 kHz core to ~1.2 MHz (Lamy+2008)
N_ELECTRIC = 73  # electric-antenna channels lead the 115-channel KEY60S frequency grid
DATA_DIR = Path("data/skr")
PDS_BASE = (
    "https://pds-ppi.igpp.ucla.edu/data/CO-V_E_J_S_SS-RPWS-4-SUMM-KEY60S-V1.0/DATA/KEY_PARAMS"
)
# Physically-motivated Saturn-rotation search band for the SKR period anchor: the SKR rotation
# periods are established at ~10.6-10.8 h (Gurnett+2009, Fischer+2015, Provan+2019), so the
# validation searches a ~5% window bracketing them. A broader 10.0-11.5 h search shows a stronger
# ~10.34 h feature (an orbital-sampling harmonic of the 6.5-day proximal orbit, not rotation) --
# reported in the findings, deliberately excluded from the rotation-period anchor.
SKR_ROT_HR = (10.4, 11.0)


def parse_scet_to_jd(scet: str) -> float:
    """SCET string ``YYYY-DDDThh:mm:ss.fff`` (day-of-year) -> Julian Date."""
    date, _, tod = scet.strip().partition("T")
    year, doy = date.split("-")
    h, m, s = tod.split(":")
    # JD of Jan 1 00:00 of the year, plus (doy-1) days plus time-of-day
    y = int(year)
    a = (14 - 1) // 12
    yy = y + 4800 - a
    mm = 1 + 12 * a - 3
    jdn = 1 + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045
    jd_jan1 = jdn - 0.5  # Jan 1 00:00 UT
    frac = (int(doy) - 1) + (int(h) + int(m) / 60.0 + float(s) / 3600.0) / 24.0
    return jd_jan1 + frac


def read_key_params(path: str | Path, *, band_hz: tuple[float, float] = SKR_BAND_HZ) -> dict:
    """Parse one daily KEY60S ``.TAB`` -> per-minute JD + SKR-band electric flux.

    Row 1 is the frequency table; rows 2+ are 1-minute spectral-density rows. Reads the electric
    frequency grid (first ``N_ELECTRIC`` of the 115 channels), then band-integrates each row's
    ELECTRIC_SPECTRAL_DENSITIES over ``band_hz``. Rows with DATA_QUALITY_FLAG != 0 are dropped.
    """
    lines = open(path).read().splitlines()
    freq_all = _read_items(lines[0], 24, 115)
    freqs = freq_all[:N_ELECTRIC]
    jd, flux, dqf = [], [], []
    for ln in lines[1:]:
        if len(ln) < 23 + N_ELECTRIC * 10:  # electric items start at index 23 (byte 24)
            continue
        q = ln[22:23].strip()
        dens = _read_items(ln, 24, N_ELECTRIC)
        jd.append(parse_scet_to_jd(ln[:21]))
        flux.append(band_integrated_flux(freqs, dens, band_hz))
        dqf.append(int(q) if q.isdigit() else 9)
    jd_a, flux_a, dqf_a = np.array(jd), np.array(flux), np.array(dqf)
    good = dqf_a == 0
    # DQF==0 rows whose band flux is still NaN (too few valid channels) are counted: they sit in
    # the duty-cycle denominator as inactive bins, so their number must be visible evidence.
    n_nan = int(np.sum(~np.isfinite(flux_a[good])))
    return {"jd": jd_a[good], "flux": flux_a[good], "freqs": freqs, "n_flux_nan": n_nan}


def _read_items(line: str, start_byte: int, n: int, *, item_bytes: int = 10) -> np.ndarray:
    """Parse ``n`` fixed-width (``item_bytes``) ASCII_REAL items starting at 1-based ``start_byte``."""
    off = start_byte - 1
    out = np.empty(n)
    for i in range(n):
        tok = line[off + i * item_bytes : off + (i + 1) * item_bytes].strip()
        try:
            out[i] = float(tok)
        except ValueError:
            out[i] = np.nan
    return out


def band_integrated_flux(
    freqs: np.ndarray, dens: np.ndarray, band_hz: tuple[float, float] = SKR_BAND_HZ
) -> float:
    r"""Integrate electric spectral density (V^2/m^2/Hz) over ``band_hz`` -> band flux (V^2/m^2).

    :math:`\int S(f)\,df` by the trapezoid rule over the log-spaced channels inside the band. NaN
    channels (fill) are skipped; returns NaN if fewer than two valid channels fall in the band.
    """
    f = np.asarray(freqs, float)
    s = np.asarray(dens, float)
    inb = (f >= band_hz[0]) & (f <= band_hz[1]) & np.isfinite(s)
    if inb.sum() < 2:
        return float("nan")
    return float(np.trapezoid(s[inb], f[inb]))


def detect_skr(flux: np.ndarray, *, k: float = 3.0, baseline_pct: float = 25.0) -> np.ndarray:
    """SKR-active bins: band flux exceeds a robust background + ``k`` sigma (the junodam pattern).

    The background is the ``baseline_pct`` percentile of the (log) flux and sigma is the robust
    MAD scatter below it, so the threshold is set by the quiescent floor, not the SKR-active bins
    themselves. Works in log space (SKR intensity spans decades). NaN bins are inactive.
    """
    f = np.asarray(flux, float)
    good = np.isfinite(f) & (f > 0)
    if good.sum() < 3:
        return np.zeros(f.shape, bool)
    lg = np.log10(f[good])
    bg = np.percentile(lg, baseline_pct)
    mad = np.median(np.abs(lg[lg <= bg] - bg)) if np.any(lg <= bg) else np.std(lg)
    sigma = 1.4826 * mad if mad > 0 else np.std(lg)
    thresh = bg + k * sigma
    out = np.zeros(f.shape, bool)
    out[good] = lg > thresh
    return out


def dual_period_ls(
    jd: np.ndarray,
    flux: np.ndarray,
    *,
    period_band_hr: tuple[float, float] = SKR_ROT_HR,
    n_freq: int = 4000,
    n_peaks: int = 2,
) -> dict:
    """Lomb-Scargle of the SKR intensity series; return the top rotation-band period peak(s).

    The published SKR rotation period (~10.7 h) is the pipeline anchor: this re-derivation
    validates the flux series before any occurrence claim. Returns the strongest ``n_peaks``
    periods (hr) within ``period_band_hr`` and the peak power + astropy false-alarm probability.
    """
    from astropy.timeseries import LombScargle

    f = np.asarray(flux, float)
    t = np.asarray(jd, float) * 24.0  # hours
    good = np.isfinite(f) & (f > 0)
    t, y = t[good], np.log10(f[good])
    y = y - y.mean()
    freq_hr = np.linspace(1.0 / period_band_hr[1], 1.0 / period_band_hr[0], n_freq)
    ls = LombScargle(t, y)
    power = ls.power(freq_hr)
    periods = 1.0 / freq_hr
    # local maxima, strongest first
    ismax = np.r_[False, (power[1:-1] > power[:-2]) & (power[1:-1] > power[2:]), False]
    idx = np.where(ismax)[0]
    idx = idx[np.argsort(power[idx])[::-1]][:n_peaks]
    peaks = [(round(float(periods[i]), 3), round(float(power[i]), 4)) for i in idx]
    best_pow = float(power.max())
    return {
        "peak_periods_hr": [p for p, _ in peaks],
        "peak_powers": [pw for _, pw in peaks],
        "best_period_hr": peaks[0][0] if peaks else float("nan"),
        "best_power": best_pow,
        "fap": float(ls.false_alarm_probability(best_pow)),
        "periods_hr": periods,
        "power": power,
    }


def proximity_duty_cycle(active: np.ndarray, range_rs: np.ndarray, *, n_bins: int = 4) -> dict:
    """SKR-active duty cycle binned by Cassini--Saturn range quartile (the junodam proximity test).

    ``range_rs`` is the Cassini--Saturn distance (Saturn radii) per bin. Returns the duty cycle in
    each of ``n_bins`` equal-count range bins (nearest first) and the near/far ratio --- the
    proximity law, reported as a DETECTION-occurrence trend (visibility-confounded; see
    `magnetic_latitude_weight`), not intrinsic emission.
    """
    a = np.asarray(active, bool)
    r = np.asarray(range_rs, float)
    good = np.isfinite(r)
    a, r = a[good], r[good]
    if a.size < n_bins:
        return {"duty_by_bin": [], "range_edges_rs": [], "near_far_ratio": float("nan")}
    edges = np.quantile(r, np.linspace(0, 1, n_bins + 1))
    duty, centers = [], []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        m = (r >= lo) & (r <= hi) if i == n_bins - 1 else (r >= lo) & (r < hi)
        duty.append(float(a[m].mean()) if m.any() else float("nan"))
        centers.append(float(np.median(r[m])) if m.any() else float("nan"))
    near, far = duty[0], duty[-1]
    ratio = near / far if far and np.isfinite(far) and far > 0 else float("inf")
    return {
        "duty_by_bin": [round(d, 4) for d in duty],
        "range_centers_rs": [round(c, 2) for c in centers],
        "range_edges_rs": [round(float(e), 2) for e in edges],
        "near_far_ratio": round(ratio, 2) if np.isfinite(ratio) else None,
    }


def distance_correct_flux(
    flux: np.ndarray, range_rs: np.ndarray, *, ref_rs: float | None = None
) -> np.ndarray:
    r"""Rescale each bin's SKR band flux to a common reference range: :math:`S\,(r/r_\mathrm{ref})^2`.

    SKR power falls as :math:`1/r^2` with observer distance, so at smaller range the *same*
    intrinsic emission clears a fixed threshold more often --- a pure sensitivity effect. Dividing
    it out (correcting every bin to ``ref_rs``, default the median range) is the null model: if the
    proximity duty-cycle trend is only sensitivity, the corrected occurrence is flat with range.
    Any residual near/far ratio after correction is the intrinsic+beaming part, bounded honestly.
    """
    r = np.asarray(range_rs, float)
    ref = float(np.nanmedian(r)) if ref_rs is None else ref_rs
    return np.asarray(flux, float) * (r / ref) ** 2


def distance_correct_excess(
    flux: np.ndarray,
    range_rs: np.ndarray,
    *,
    ref_rs: float | None = None,
    baseline_pct: float = 25.0,
) -> np.ndarray:
    r"""The 1/r^2 sensitivity null, applied to the emission EXCESS over a range-independent floor.

    ``distance_correct_flux`` rescales the *total* band flux, but the quiescent floor in this band
    is instrumental + galactic background that does not scale with the Cassini--Saturn distance ---
    rescaling it imposes an ~0.86 dex range-dependent offset on the floor itself, which biases the
    corrected far-bin duty cycle up and the near-bin down, i.e. in the direction that manufactures
    a collapse (referee finding, 2026-08-24; the sibling ``junodam`` census corrects the SNR, the
    equivalent choice). This version holds the floor fixed: the floor is the ``baseline_pct``
    percentile of the finite flux, the excess above it is scaled by :math:`(r/r_\mathrm{ref})^2`,
    and the floor is added back, so re-detection sees a range-independent background.
    """
    f = np.asarray(flux, float)
    r = np.asarray(range_rs, float)
    ref = float(np.nanmedian(r)) if ref_rs is None else ref_rs
    good = np.isfinite(f) & (f > 0)
    floor = float(np.percentile(f[good], baseline_pct)) if good.any() else 0.0
    # Bins at or below the floor keep their ORIGINAL values: clipping them onto the floor would
    # pile ~baseline_pct% of bins at one exact value, degenerate the MAD the re-detection
    # threshold is built from, and inflate the threshold through detect_skr's std fallback. The
    # quiescent noise distribution is range-independent, so leaving it untouched is the point.
    return np.where(f > floor, floor + (f - floor) * (r / ref) ** 2, f)


def orbit_segments(jd: np.ndarray, range_rs: np.ndarray, *, min_sep_days: float = 3.0) -> list:
    """Split the series into periapsis passes at the apoapsis maxima of the range series.

    The natural resampling unit of this census is the orbit (~9 proximal passes in 59 days), not
    the 83,382 autocorrelated one-minute bins. Returns a list of ``(start, stop)`` index pairs,
    one per segment between consecutive apoapses (ends included as partial segments).
    """
    t = np.asarray(jd, float)
    r = np.asarray(range_rs, float)
    n = t.size
    if n < 10:
        return [(0, n)]
    # Apoapsis candidates are the top decile of the range series; picked greedily by range with
    # a minimum separation so each pass contributes exactly one boundary. Restricting to the top
    # decile is what stops the greedy pass from eventually accepting periapsis bins once every
    # true apoapsis is taken.
    finite = np.isfinite(r)
    if not finite.any():
        return [(0, n)]
    hi_thresh = np.nanpercentile(r, 90.0)
    cand = np.where(finite & (r >= hi_thresh))[0]
    cand = cand[np.argsort(-r[cand])]
    apo: list[int] = []
    for i in cand:
        if all(abs(t[i] - t[j]) > min_sep_days for j in apo):
            apo.append(int(i))
    apo = sorted(a for a in apo if 0 < a < n - 1)
    bounds = [0] + apo + [n]
    return [(lo, hi) for lo, hi in zip(bounds[:-1], bounds[1:], strict=True) if hi - lo > 60]


def jackknife_near_far(
    jd: np.ndarray,
    flux: np.ndarray,
    range_rs: np.ndarray,
    *,
    mode: str = "common",
    k: float | None = None,
    baseline_pct: float = 25.0,
) -> dict:
    """Leave-one-orbit-out jackknife of the near/far duty-cycle ratio.

    Every sentence in the census is a claim about the *size* of a ratio, and the orbit-to-orbit
    scatter is the variance a within-realization estimate cannot see (the ``rmstructure``
    lesson). Drops one periapsis pass at a time, recomputes the ratio (quartile edges re-derived
    per sample, consistent with the estimator), and returns the full-sample value, the jackknife
    mean, and the jackknife standard error.
    """

    def _active(fx: np.ndarray, rr: np.ndarray) -> np.ndarray:
        if mode == "common":
            return common_sensitivity_active(
                fx, rr, k=6.0 if k is None else k, baseline_pct=baseline_pct
            )
        if mode == "excess":
            return detect_skr(
                distance_correct_excess(fx, rr, baseline_pct=baseline_pct),
                k=3.0 if k is None else k,
                baseline_pct=baseline_pct,
            )
        return detect_skr(fx, k=3.0 if k is None else k, baseline_pct=baseline_pct)

    segs = orbit_segments(jd, range_rs)
    full = proximity_duty_cycle(_active(flux, range_rs), range_rs)
    vals = []
    for lo, hi in segs:
        m = np.ones(np.asarray(jd).size, bool)
        m[lo:hi] = False
        p = proximity_duty_cycle(_active(flux[m], range_rs[m]), range_rs[m])
        v = p["near_far_ratio"]
        if v is not None and np.isfinite(v):
            vals.append(float(v))
    n = len(vals)
    if n < 2:
        return {"full": full["near_far_ratio"], "jk_mean": None, "jk_se": None, "n_orbits": n}
    arr = np.asarray(vals)
    jk_mean = float(arr.mean())
    jk_se = float(np.sqrt((n - 1) / n * np.sum((arr - jk_mean) ** 2)))
    return {
        "full": full["near_far_ratio"],
        "jk_mean": round(jk_mean, 2),
        "jk_se": round(jk_se, 2),
        "n_orbits": n,
    }


def threshold_sweep(
    flux: np.ndarray,
    range_rs: np.ndarray,
    *,
    ks=(4.0, 6.0, 8.0, 10.0),
    baselines=(10.0, 25.0, 50.0),
) -> list[dict]:
    """The common-sensitivity near/far ratio over the detection-rule grid.

    The headline is a ratio of threshold-crossing rates on a steeply falling flux distribution,
    which is where threshold choice bites hardest; a number quoted from one (k, baseline) pair
    carries an unmeasured rule dependence. One row per grid point.
    """
    out = []
    for bp in baselines:
        for k in ks:
            p = proximity_duty_cycle(
                common_sensitivity_active(flux, range_rs, k=k, baseline_pct=bp), range_rs
            )
            out.append(
                {
                    "k": k,
                    "baseline_pct": bp,
                    "corrected_near_far": p["near_far_ratio"],
                }
            )
    return out


def block_permutation_p(
    jd: np.ndarray,
    flux: np.ndarray,
    *,
    n_perm: int = 199,
    block_days: float = 1.0,
    period_band_hr: tuple[float, float] = SKR_ROT_HR,
    seed: int = 0,
) -> dict:
    """A day-block permutation null for the rotation-band LS peak.

    Astropy's analytic false-alarm probability assumes independent samples; 83,382 one-minute
    bins of an emission that persists for hours are massively autocorrelated, so it underflows
    regardless of the truth. Permuting whole day blocks (flux blocks shuffled onto the fixed time
    grid) preserves the within-day autocorrelation and destroys only the multi-day phase
    coherence a rotation period requires. p = (k+1)/(n+1), floored at 1/(n+1).
    """
    rng = np.random.default_rng(seed)
    t = np.asarray(jd, float)
    f = np.asarray(flux, float)
    obs = dual_period_ls(t, f, period_band_hr=period_band_hr)["best_power"]
    day = np.floor(t / block_days).astype(int)
    blocks = [f[day == d] for d in np.unique(day)]
    exceed = 0
    for _ in range(n_perm):
        rng.shuffle(blocks)
        fp = np.concatenate(blocks)
        if fp.size != f.size:  # ragged edge days -- pad/trim defensively
            fp = fp[: f.size] if fp.size > f.size else np.pad(fp, (0, f.size - fp.size))
        if dual_period_ls(t, fp, period_band_hr=period_band_hr)["best_power"] >= obs:
            exceed += 1
    return {"p_perm": round((exceed + 1) / (n_perm + 1), 3), "n_perm": n_perm, "obs_power": obs}


def common_sensitivity_active(
    flux: np.ndarray,
    range_rs: np.ndarray,
    *,
    k: float = 6.0,
    baseline_pct: float = 25.0,
    far_pct: float = 95.0,
) -> np.ndarray:
    r"""SKR-active bins at a COMMON sensitivity: the census trimmed to the far-range threshold.

    The 1/r^2 null has two failure modes this estimator avoids. Rescaling the *total* flux moves
    the range-independent noise floor with the signal; rescaling the *excess* upward at far range
    amplifies far-range noise deviations and manufactures false positives there (both measured on
    controls in ``tests/test_skr.py``). The clean question is dr20radio's common-limit move:
    *which bins would still be detected if every bin were observed from the far reference range?*
    Each bin's excess over the range-independent floor is scaled by :math:`(r/r_\mathrm{far})^2
    \le 1` --- only ever *down*, so noise is never promoted --- and compared against the one
    threshold the RAW series defines. Near bins keep only emission bright enough to be seen from
    ``far_pct``-percentile range; the resulting duty-cycle ratio compares intrinsic occurrence
    above a common effective luminosity limit.

    ``k`` defaults to 6, stiffer than the raw census's 3, because ``detect_skr``'s MAD-below-
    percentile sigma underestimates the true noise width by ~30%: at k=3 a few percent of pure
    noise crosses, uniformly in range for the raw census but range-ASYMMETRICALLY once the real
    signal is censored to its bright tail. k=6 is where the control test's noise crossings reach
    zero while a flat intrinsic rate comes out flat (``tests/test_skr.py``).
    """
    f = np.asarray(flux, float)
    r = np.asarray(range_rs, float)
    good = np.isfinite(f) & (f > 0) & np.isfinite(r)
    if good.sum() < 3:
        return np.zeros(f.shape, bool)
    floor = float(np.percentile(f[good], baseline_pct))
    lg = np.log10(f[good])
    bg = np.percentile(lg, baseline_pct)
    mad = np.median(np.abs(lg[lg <= bg] - bg)) if np.any(lg <= bg) else np.std(lg)
    sigma = 1.4826 * mad if mad > 0 else np.std(lg)
    thresh_excess = 10.0 ** (bg + k * sigma) - floor  # the raw series' excess threshold
    r_far = float(np.nanpercentile(r[good], far_pct))
    excess = np.clip(f - floor, 0.0, None)
    # The census is a CONJUNCTION: the bin must be raw-active (which controls the noise-crossing
    # rate in the log domain, identically at every range) AND its excess must survive the scale
    # to the far reference. Censoring alone would let the raw detector's small noise-crossing
    # rate dominate the far bins once the real signal is trimmed to its bright tail.
    out = detect_skr(f, k=k, baseline_pct=baseline_pct)
    out[good] &= excess[good] * np.clip(r[good] / r_far, None, 1.0) ** 2 > thresh_excess
    return out


def latitude_by_range_bin(
    range_rs: np.ndarray, sub_lat_deg: np.ndarray, *, n_bins: int = 4
) -> dict:
    """Sub-spacecraft |latitude| range spanned by each Cassini-range quartile.

    The latitude confound is BETWEEN range bins (range and latitude are correlated along the
    orbit): if the near and far quartiles sample different latitudes, the latitude dependence of
    SKR visibility is entangled with range and is NOT removed by any within-bin reweighting. This
    reports each bin's median and span of ``|sub_lat|`` so the confound is visible, not asserted.
    """
    r = np.asarray(range_rs, float)
    lat = np.abs(np.asarray(sub_lat_deg, float))
    good = np.isfinite(r) & np.isfinite(lat)
    r, lat = r[good], lat[good]
    if r.size < n_bins:
        return {"abs_lat_median_by_bin": [], "abs_lat_span_deg": float("nan")}
    edges = np.quantile(r, np.linspace(0, 1, n_bins + 1))
    meds = []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        m = (r >= lo) & (r <= hi) if i == n_bins - 1 else (r >= lo) & (r < hi)
        meds.append(round(float(np.median(lat[m])), 1) if m.any() else float("nan"))
    return {
        "abs_lat_median_by_bin": meds,
        "abs_lat_span_deg": round(float(max(meds) - min(meds)), 1),
    }


def magnetic_latitude_weight(
    active: np.ndarray, range_rs: np.ndarray, sub_lat_deg: np.ndarray, *, n_bins: int = 4
) -> dict:
    """Latitude-weighted proximity duty cycle: divide out SKR viewing-geometry visibility.

    SKR is auroral and beamed, brightest when the sub-spacecraft point is at low-to-mid latitude
    on the emitting hemisphere; occurrence-vs-range is therefore confounded with the latitude
    coverage of each range bin. This applies a simple, STATED-model visibility weight
    ``w(lat) = |sin(lat)|`` (a dipole-auroral-beaming proxy: emission favours higher magnetic
    latitude viewing) and reports the weighted duty cycle per range bin. The model dependence is
    the point: if the near/far ratio survives latitude weighting, WITHIN-BIN latitude
    coverage is not the whole story --- but the confound is BETWEEN bins (see
    `latitude_by_range_bin`), which no within-bin reweighting removes, so this metric cannot
    clear proximity of the latitude confound and the paper does not cite it. This is a proxy, not a radiative-transfer model --- caveated.
    """
    a = np.asarray(active, bool)
    r = np.asarray(range_rs, float)
    lat = np.asarray(sub_lat_deg, float)
    good = np.isfinite(r) & np.isfinite(lat)
    a, r, lat = a[good], r[good], lat[good]
    w = np.abs(np.sin(np.radians(lat))) + 1e-3  # visibility weight, floored to avoid /0
    if a.size < n_bins:
        return {"weighted_duty_by_bin": [], "weighted_near_far_ratio": float("nan")}
    edges = np.quantile(r, np.linspace(0, 1, n_bins + 1))
    wduty = []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        m = (r >= lo) & (r <= hi) if i == n_bins - 1 else (r >= lo) & (r < hi)
        if m.any():
            wduty.append(float(np.sum(a[m] * w[m]) / np.sum(w[m])))
        else:
            wduty.append(float("nan"))
    near, far = wduty[0], wduty[-1]
    ratio = near / far if far and np.isfinite(far) and far > 0 else float("inf")
    return {
        "weighted_duty_by_bin": [round(d, 4) for d in wduty],
        "weighted_near_far_ratio": round(ratio, 2) if np.isfinite(ratio) else None,
    }


def synthetic_skr(
    *,
    n_days: float = 8.0,
    cadence_min: float = 1.0,
    period_hr: float = 10.7,
    second_period_hr: float = 10.6,
    near_far_contrast: float = 6.0,
    seed: int = 0,
) -> dict:
    """A synthetic SKR series with a KNOWN dual period + range-dependent occurrence, for recovery.

    Builds a 1-min flux series modulated at two close periods, on a range track that sweeps from
    near to far; the SKR-active probability rises toward periapsis by ``near_far_contrast``. The
    recover-a-known: `dual_period_ls` must find ~``period_hr`` and `proximity_duty_cycle` must
    recover a near/far ratio near ``near_far_contrast``.
    """
    rng = np.random.default_rng(seed)
    n = int(n_days * 24 * 60 / cadence_min)
    jd = 2.456e6 + np.arange(n) * (cadence_min / 60.0 / 24.0)
    t_hr = np.arange(n) * cadence_min / 60.0
    # a range track: sinusoidal periapsis sweep (Rs), 3 -> 60 Rs
    range_rs = 31.5 - 28.5 * np.cos(2 * np.pi * t_hr / (n_days * 24.0))
    sub_lat = 20.0 * np.sin(2 * np.pi * t_hr / (n_days * 24.0 / 3))
    # dual-period intensity modulation
    mod = 0.5 * np.sin(2 * np.pi * t_hr / period_hr) + 0.5 * np.sin(
        2 * np.pi * t_hr / second_period_hr
    )
    base = 10.0 ** (rng.normal(-14.0, 0.15, n))  # quiescent electric spectral floor scale
    # occurrence probability: higher near periapsis (proximity) and near intensity-mod maxima
    p_near = 1.0 / range_rs
    p_active = (p_near / p_near.max()) * 0.5 * (1 + np.tanh(3 * mod))
    p_active *= near_far_contrast * 0.05
    active_true = rng.random(n) < np.clip(p_active, 0, 1)
    flux = base * (1 + active_true * 10 ** rng.uniform(1.5, 3.0, n))
    return {
        "jd": jd,
        "flux": flux,
        "range_rs": range_rs,
        "sub_lat_deg": sub_lat,
        "active_true": active_true,
        "period_hr": period_hr,
        "near_far_contrast": near_far_contrast,
    }


def fetch_geometry(
    jd_start: float, jd_stop: float, *, step: str = "10m"
) -> dict:  # pragma: no cover - network
    """Cassini--Saturn range (Saturn radii) + sub-Cassini latitude from JPL Horizons.

    TARGET=699 (Saturn centre), CENTER=500@-82 (Cassini): quantity 20 -> delta (range),
    quantity 14 -> observer sub-lon/lat on Saturn.
    """
    from astroquery.jplhorizons import Horizons

    obj = Horizons(
        id="699",
        location="500@-82",
        epochs={"start": f"JD{jd_start}", "stop": f"JD{jd_stop}", "step": step},
    )
    eph = obj.ephemerides(quantities="14,20")
    au_per_rs = 60268.0 / 1.495978707e8  # Saturn equatorial radius (km) in AU
    return {
        "jd": np.asarray(eph["datetime_jd"], float),
        "range_rs": np.asarray(eph["delta"], float) / au_per_rs,
        "sub_lat_deg": np.asarray(eph["PDObsLat"], float),
    }


def fetch_rpws_key60s(
    year: int, doy: int, seq: int = 3, *, out_dir: str | Path = DATA_DIR
) -> Path:  # pragma: no cover - network
    """Download one daily KEY60S ``.TAB`` from PDS-PPI into ``out_dir``; returns the local path."""
    import urllib.request

    month_dir = f"T{year}{doy // 100}XX"  # PDS-PPI monthly directory convention for KEY60S
    name = f"RPWS_KEY__{year}{doy:03d}_{seq}.TAB"
    url = f"{PDS_BASE}/{month_dir}/{name}"
    dest = Path(out_dir) / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        urllib.request.urlretrieve(url, dest)
    return dest


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: synthetic recover-a-known. Real: parse local KEY60S days + Horizons geometry."""

    if offline:
        s = synthetic_skr()
        jd, flux, range_rs, sub_lat = s["jd"], s["flux"], s["range_rs"], s["sub_lat_deg"]
        source = "synthetic SKR series (injected dual period + proximity trend)"
        n_flux_nan = 0
        expected_ratio = s["near_far_contrast"]
        expected_period = s["period_hr"]
    else:  # pragma: no cover - data files + network
        files = sorted(DATA_DIR.glob("RPWS_KEY__*.TAB"))
        parts = [read_key_params(f) for f in files]
        jd = np.concatenate([p["jd"] for p in parts])
        flux = np.concatenate([p["flux"] for p in parts])
        order = np.argsort(jd)
        jd, flux = jd[order], flux[order]
        starts = [0] + list(np.where(np.diff(jd) > 2.0)[0] + 1) + [jd.size]
        range_rs = np.full(jd.size, np.nan)
        sub_lat = np.full(jd.size, np.nan)
        for a, b in zip(starts[:-1], starts[1:], strict=True):
            eph = fetch_geometry(float(jd[a]) - 0.02, float(jd[b - 1]) + 0.02)
            range_rs[a:b] = np.interp(jd[a:b], eph["jd"], eph["range_rs"])
            sub_lat[a:b] = np.interp(jd[a:b], eph["jd"], eph["sub_lat_deg"])
        source = f"Cassini/RPWS KEY60S, {len(files)} days"
        n_flux_nan = int(sum(p.get("n_flux_nan", 0) for p in parts))
        expected_ratio = float("nan")
        expected_period = float("nan")

    active = detect_skr(flux)
    ls = dual_period_ls(jd, flux)
    prox = proximity_duty_cycle(active, range_rs)
    latw = magnetic_latitude_weight(active, range_rs, sub_lat)
    latbin = latitude_by_range_bin(range_rs, sub_lat)
    # The 1/r^2 sensitivity null, applied to the EXCESS over a range-independent floor (the
    # 2026-08-24 referee fix: rescaling total flux moved the noise floor with the signal, in the
    # direction that manufactures a collapse). The old total-flux variant is kept as a recorded
    # comparison so the change is auditable.
    # PRIMARY: the common-sensitivity census (only ever scales excesses DOWN; noise is never
    # promoted; k=6 where the control's noise crossings vanish). The two superseded null models
    # are kept as recorded comparisons: the total-flux rescale moves the noise floor with the
    # signal (biased toward a collapse), the excess-up-rescale amplifies far-range noise
    # (biased toward a reversal). Their spread is the method-dependence evidence.
    active_corr = common_sensitivity_active(flux, range_rs)
    prox_corr = proximity_duty_cycle(active_corr, range_rs)
    prox_corr_excess = proximity_duty_cycle(
        detect_skr(distance_correct_excess(flux, range_rs)), range_rs
    )
    prox_corr_flux = proximity_duty_cycle(
        detect_skr(distance_correct_flux(flux, range_rs)), range_rs
    )
    # Orbit-level uncertainty, the detection-rule sweep, the day-block permutation null for the
    # rotation peak, and the broad-band periodogram the excluded ~10.34 h claim needs evidence for.
    jk_raw = jackknife_near_far(jd, flux, range_rs, mode="raw")
    jk_corr = jackknife_near_far(jd, flux, range_rs, mode="common")
    sweep = threshold_sweep(flux, range_rs)
    perm = block_permutation_p(jd, flux)
    ls_broad = dual_period_ls(jd, flux, period_band_hr=(10.0, 11.5))
    # Robustness variants: exclude the near-Saturn bins where far-field 1/r^2 fails (and where
    # ring-plane dust contamination lives), and a per-day background in place of the pooled one.
    far_field = ~(np.isfinite(range_rs) & (range_rs < 5.0))
    prox_ff = proximity_duty_cycle(
        common_sensitivity_active(flux[far_field], range_rs[far_field]), range_rs[far_field]
    )
    day_idx = np.floor(np.asarray(jd, float)).astype(int)
    active_perday = np.zeros(flux.size, bool)
    for d in np.unique(day_idx):
        m = day_idx == d
        active_perday[m] = common_sensitivity_active(flux[m], range_rs[m])
    prox_perday = proximity_duty_cycle(active_perday, range_rs)

    metrics: dict = {
        "source": source,
        "is_real": not offline,
        "n_bins": int(active.size),
        "n_active": int(active.sum()),
        "n_flux_nan": n_flux_nan,
        "duty_cycle_pct": round(100.0 * float(active.mean()), 3),
        "anchor_period_hr": ls["best_period_hr"],
        # deviation of the recovered dominant period from the published late-mission SKR periods
        # (Provan+2019: S 10.68 h, N 10.79 h) -- the meaningful validation, not the search-band
        # tautology. Matched to the nearer of the two.
        "anchor_dev_pct": round(
            100.0 * min(abs(ls["best_period_hr"] - p) / p for p in (10.68, 10.79)), 2
        )
        if np.isfinite(ls["best_period_hr"])
        else None,
        "peak_periods_hr": ls["peak_periods_hr"],
        "period_two_hr": ls["peak_periods_hr"][1] if len(ls["peak_periods_hr"]) > 1 else None,
        "duty_by_range_bin": prox["duty_by_bin"],
        "range_centers_rs": prox["range_centers_rs"],
        "range_near_rs": prox["range_centers_rs"][0] if prox["range_centers_rs"] else None,
        "range_far_rs": prox["range_centers_rs"][-1] if prox["range_centers_rs"] else None,
        "duty_near_pct": round(100 * prox["duty_by_bin"][0], 1) if prox["duty_by_bin"] else None,
        "duty_far_pct": round(100 * prox["duty_by_bin"][-1], 1) if prox["duty_by_bin"] else None,
        "near_far_ratio": prox["near_far_ratio"],
        "range_edges_rs": prox["range_edges_rs"],
        # The primary null: the common-sensitivity census (see common_sensitivity_active).
        "sensitivity_corrected_near_far": prox_corr["near_far_ratio"],
        "duty_by_range_bin_corrected": prox_corr["duty_by_bin"],
        "n_active_corrected": int(active_corr.sum()),
        "common_census_k": 6.0,
        # Superseded null models, kept as method-dependence evidence: the total-flux rescale is
        # biased toward a collapse (it moves the noise floor with the signal), the excess-up
        # rescale toward a reversal (it amplifies far-range noise).
        "sensitivity_corrected_near_far_fluxscaled": prox_corr_flux["near_far_ratio"],
        "sensitivity_corrected_near_far_excess_up": prox_corr_excess["near_far_ratio"],
        # Orbit-level uncertainty (leave-one-periapsis-out), the number every size claim needs.
        "near_far_jackknife": jk_raw,
        "corrected_near_far_jackknife": jk_corr,
        # Detection-rule dependence of the corrected ratio.
        "threshold_sweep": sweep,
        # Day-block permutation significance for the rotation-band peak (replaces the analytic
        # ls_fap, which assumes independent samples and underflows on autocorrelated minutes).
        "rotation_peak_perm": perm,
        # Broad-band periodogram peaks, so the excluded ~10.34 h feature is committed evidence.
        "broad_peak_periods_hr": ls_broad["peak_periods_hr"],
        "broad_peak_powers": ls_broad["peak_powers"],
        # Robustness variants of the corrected ratio.
        "corrected_near_far_farfield_only": prox_ff["near_far_ratio"],
        "corrected_near_far_perday_background": prox_perday["near_far_ratio"],
        # The detection rule and null-model parameters, so the census is reproducible from the
        # committed evidence alone.
        "detect_k": 3.0,
        "detect_baseline_pct": 25.0,
        "band_hz": list(SKR_BAND_HZ),
        "ref_rs": round(float(np.nanmedian(range_rs)), 2),
        "jd_range": [round(float(np.nanmin(jd)), 2), round(float(np.nanmax(jd)), 2)],
        "weighted_near_far_ratio": latw["weighted_near_far_ratio"],
        "abs_lat_median_by_bin": latbin["abs_lat_median_by_bin"],
        "abs_lat_span_deg": latbin["abs_lat_span_deg"],
        "expected_near_far_ratio": round(expected_ratio, 2)
        if np.isfinite(expected_ratio)
        else None,
        "expected_period_hr": round(expected_period, 2) if np.isfinite(expected_period) else None,
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "skr_metrics.json")
    _figure(ls, prox, op / "papers" / "skr" / "figures")
    _write_macros(metrics, op / "papers" / "skr" / "generated" / "macros.tex")
    return metrics


def _figure(ls: dict, prox: dict, out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.9))
    ax1.plot(ls["periods_hr"], ls["power"], color="C0", lw=0.8)
    ax1.axvspan(10.5, 10.9, color="C3", alpha=0.2, label="published ~10.6-10.8 h")
    ax1.set(xlabel="period (hr)", ylabel="LS power", title="SKR rotation-period anchor")
    ax1.legend(fontsize=8)
    duty = prox.get("duty_by_bin") or []
    centers = prox.get("range_centers_rs") or list(range(len(duty)))
    if duty:
        ax2.plot(centers, [100 * d for d in duty], "o-", color="C0")
        ax2.set(
            xlabel="Cassini--Saturn range (Rs)",
            ylabel="SKR-active duty cycle (%)",
            title="Proximity duty-cycle law",
        )
    else:
        ax2.set_axis_off()
    fig.tight_layout()
    fig.savefig(out / "skr.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    def g(key: str) -> str:
        v = m.get(key)
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    pref = "skrReal" if m.get("is_real") else "skrSyn"
    lines = [
        "% Auto-generated by jansky_research.skr._write_macros -- do not edit.",
        "% Synthetic (skrSyn*) and real (skrReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders (an offline rebuild resets skrReal* to '--').",
        rf"\newcommand{{\skrSource}}{{{m['source']}}}",
        rf"\newcommand{{\skrNBins}}{{{m['n_bins']}}}",
        rf"\newcommand{{\skrDutyPct}}{{{g('duty_cycle_pct')}}}",
        rf"\newcommand{{\skrAnchorPeriod}}{{{g('anchor_period_hr')}}}",
    ]
    for ns in ("skrSyn", "skrReal"):
        live = ns == pref
        for macro, key in (
            ("NearFar", "near_far_ratio"),
            ("SensCorrNearFar", "sensitivity_corrected_near_far"),
            ("WeightedNearFar", "weighted_near_far_ratio"),
            ("AbsLatSpan", "abs_lat_span_deg"),
            ("AnchorPeriod", "anchor_period_hr"),
            ("AnchorDev", "anchor_dev_pct"),
            ("PeriodTwo", "period_two_hr"),
            ("RangeNear", "range_near_rs"),
            ("RangeFar", "range_far_rs"),
            ("DutyNear", "duty_near_pct"),
            ("DutyFar", "duty_far_pct"),
            ("DutyPct", "duty_cycle_pct"),
            ("NFluxNan", "n_flux_nan"),
        ):
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(key) if live else '--'}}}")
        # The revision's evidence: jackknife errors, the rule sweep's spread, the permutation
        # significance, the broad-band peaks, and the robustness variants.
        jr = m.get("near_far_jackknife") or {}
        jc = m.get("corrected_near_far_jackknife") or {}
        sw = m.get("threshold_sweep") or []
        pm = m.get("rotation_peak_perm") or {}
        vals = [s["corrected_near_far"] for s in sw if s.get("corrected_near_far") is not None]
        derived = (
            ("NearFarJkSe", jr.get("jk_se") if live else None),
            ("CorrNearFarJkMean", jc.get("jk_mean") if live else None),
            ("CorrNearFarJkSe", jc.get("jk_se") if live else None),
            ("NOrbits", jc.get("n_orbits") if live else None),
            ("SweepMin", round(min(vals), 2) if (live and vals) else None),
            ("SweepMax", round(max(vals), 2) if (live and vals) else None),
            ("PermP", pm.get("p_perm") if live else None),
            (
                "CorrNearFarFluxScaled",
                m.get("sensitivity_corrected_near_far_fluxscaled") if live else None,
            ),
            ("CorrNearFarFarField", m.get("corrected_near_far_farfield_only") if live else None),
            ("CorrNearFarPerDay", m.get("corrected_near_far_perday_background") if live else None),
            ("RangeMin", m.get("range_edges_rs", [None])[0] if live else None),
            ("RangeMax", m.get("range_edges_rs", [None])[-1] if live else None),
            ("BroadPeak", (m.get("broad_peak_periods_hr") or [None])[0] if live else None),
        )
        for macro, v in derived:
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{'--' if v is None else v}}}")
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

    p = argparse.ArgumentParser(description="Cassini SKR occurrence + proximity census.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
