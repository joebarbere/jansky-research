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
    return {"jd": jd_a[good], "flux": flux_a[good], "freqs": freqs}


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
    the point: if the near/far ratio survives latitude weighting, proximity is not purely a
    latitude-coverage artefact. This is a proxy, not a radiative-transfer model --- caveated.
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

    month_dir = f"T{year}2XX"  # PDS-PPI monthly directory convention for KEY60S
    name = f"RPWS_KEY__{year}{doy:03d}_{seq}.TAB"
    url = f"{PDS_BASE}/{month_dir}/{name}"
    dest = Path(out_dir) / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        urllib.request.urlretrieve(url, dest)
    return dest


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: synthetic recover-a-known. Real: parse local KEY60S days + Horizons geometry."""
    import json

    if offline:
        s = synthetic_skr()
        jd, flux, range_rs, sub_lat = s["jd"], s["flux"], s["range_rs"], s["sub_lat_deg"]
        source = "synthetic SKR series (injected dual period + proximity trend)"
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
        expected_ratio = float("nan")
        expected_period = float("nan")

    active = detect_skr(flux)
    ls = dual_period_ls(jd, flux)
    prox = proximity_duty_cycle(active, range_rs)
    latw = magnetic_latitude_weight(active, range_rs, sub_lat)
    latbin = latitude_by_range_bin(range_rs, sub_lat)
    # the 1/r^2 sensitivity null: correct flux to a common range, re-detect, re-bin. If the
    # near/far trend is pure visibility, this ratio -> ~1; any residual is intrinsic+beaming.
    active_corr = detect_skr(distance_correct_flux(flux, range_rs))
    prox_corr = proximity_duty_cycle(active_corr, range_rs)

    metrics: dict = {
        "source": source,
        "is_real": not offline,
        "n_bins": int(active.size),
        "n_active": int(active.sum()),
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
        "ls_fap": round(ls["fap"], 4) if np.isfinite(ls["fap"]) else None,
        "duty_by_range_bin": prox["duty_by_bin"],
        "range_centers_rs": prox["range_centers_rs"],
        "range_near_rs": prox["range_centers_rs"][0] if prox["range_centers_rs"] else None,
        "range_far_rs": prox["range_centers_rs"][-1] if prox["range_centers_rs"] else None,
        "duty_near_pct": round(100 * prox["duty_by_bin"][0], 1) if prox["duty_by_bin"] else None,
        "duty_far_pct": round(100 * prox["duty_by_bin"][-1], 1) if prox["duty_by_bin"] else None,
        "near_far_ratio": prox["near_far_ratio"],
        "sensitivity_corrected_near_far": prox_corr["near_far_ratio"],
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
    (op / "results" / "skr_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
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
            ("Fap", "ls_fap"),
            ("DutyPct", "duty_cycle_pct"),
        ):
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(key) if live else '--'}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


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
