"""Jovian decametric (DAM) occurrence census from Juno/Waves in the (CML, Io-phase) plane (plan 37).

The classic ground-based result (Bigg 1964; Carr, Desch & Alexander 1983; Marques et al. 2017's
26-yr Nancay catalogue) is that Jovian DAM occurrence is organised by the observer's System III
central meridian longitude (CML) and Io's orbital phase into the Io-A/B/C/D source regions. The
public Juno/Waves Estimated Flux Density Dataset (doi:10.25935/6jg4-mk86; daily CDFs, 110 channels
to 40.5 MHz, per-channel ``Background``/``Sigma``) lets an occurrence census be built **from
Juno's vantage** --- the recover-a-known is Io-controlled enhancement; the vantage is the new part.

Two ephemeris subtleties (GATE 0, live-verified): the sub-Juno CML must come from JPL Horizons
(``PDObsLon``; the naive IAU :math:`W_{III}` rotation formula is wrong for Juno by up to ~40 deg),
and Io's phase uses the Lieske (1987) mean longitude
:math:`l_1 = 106.07719 + 203.488955790\\,(JD-2451545)` deg minus the CML. Units are
:math:`V^2 m^{-2} Hz^{-1}` (spectral power density, not W-flux) --- fine for an *occurrence*
census, which needs only detection above the shipped background.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "io_mean_longitude",
    "io_phase",
    "read_waves_cdf",
    "detect_active",
    "occurrence_map",
    "io_region_contrast",
    "synthetic_orbit",
    "run",
]

J2000_JD = 2451545.0
#: Lieske (1987, A&A 176, 146) E2x5 mean longitude of Io (deg, deg/day) --- <0.2 deg vs Horizons.
IO_L0, IO_RATE = 106.07719, 203.488955790
DAM_BAND_MHZ = (3.0, 40.5)
DATA_DIR = Path("data/junodam")
RJ_AU = 71492.0 / 1.495978707e8  # Jupiter equatorial radius in AU (for range display)
CDF_URL = (
    "https://maser.obspm.fr/repository/juno/waves/data/l3a/data/cdf/"
    "{y}/{m:02d}/jno_wav_cdr_lesia_{y}{m:02d}{d:02d}_v01.cdf"
)


def io_mean_longitude(jd: np.ndarray) -> np.ndarray:
    """Io's mean orbital longitude (deg, System III-adjacent frame; Lieske 1987)."""
    return (IO_L0 + IO_RATE * (np.asarray(jd, float) - J2000_JD)) % 360.0


def io_phase(jd: np.ndarray, cml_deg: np.ndarray) -> np.ndarray:
    r"""Conventional Io phase :math:`\Phi_{Io}` --- Io's departure from SUPERIOR conjunction.

    :math:`\Phi_{Io} = \mathrm{CML} + 180^\circ - \Lambda_{Io}` (Bigg 1964 convention as
    codified by Carr et al. 1983): zero when Io is diametrically opposite the observer as seen
    from Jupiter. GATE-2 caught an earlier sign/offset error here (:math:`\Lambda_{Io} -
    \mathrm{CML}` = :math:`180^\circ - \Phi_{Io}`), which displaced every canonical box.
    """
    return (np.asarray(cml_deg, float) + 180.0 - io_mean_longitude(jd)) % 360.0


def fetch_cml_horizons(
    jd_start: float, jd_stop: float, *, step: str = "15m"
) -> dict:  # pragma: no cover - network
    """Sub-Juno System III CML from JPL Horizons (``PDObsLon``; observer ``500@-61``)."""
    from astroquery.jplhorizons import Horizons

    obj = Horizons(
        id="599",
        location="500@-61",
        epochs={"start": f"JD{jd_start}", "stop": f"JD{jd_stop}", "step": step},
    )
    eph = obj.ephemerides(quantities="14,20")
    return {
        "jd": np.asarray(eph["datetime_jd"], float),
        "cml": np.asarray(eph["PDObsLon"], float),
        "delta_au": np.asarray(eph["delta"], float),  # Juno--Jupiter range
    }


def read_waves_cdf(path: str | Path, *, bin_s: int = 15) -> dict:  # pragma: no cover - data file
    """Read one daily Juno/Waves L3a CDF -> DAM-band activity per ``bin_s`` time bin.

    Returns per-bin JD and the fraction of DAM-band channels whose median power in the bin
    exceeds ``Background + 5 Sigma`` (the shipped per-channel statistics).
    """
    import cdflib

    cdf = cdflib.CDF(str(path))
    epoch = cdflib.cdfepoch.to_datetime(cdf.varget("Epoch"))
    jd = np.array([e.astype("datetime64[s]").astype(float) for e in epoch]) / 86400.0 + 2440587.5
    freq_mhz = np.asarray(cdf.varget("Frequency"), float) / 1e3  # kHz -> MHz
    band = (freq_mhz >= DAM_BAND_MHZ[0]) & (freq_mhz <= DAM_BAND_MHZ[1])
    data = np.asarray(cdf.varget("Data"), float)[:, band]
    bg = np.asarray(cdf.varget("Background"), float)[band]
    sig = np.asarray(cdf.varget("Sigma"), float)[band]

    n = (data.shape[0] // bin_s) * bin_s
    d = data[:n].reshape(-1, bin_s, band.sum())
    med = np.median(d, axis=1)
    floor = bg + 5.0 * sig
    active_frac = (med > floor).mean(axis=1)
    # per-bin 90th-percentile channel SNR vs the 5-sigma floor: snr_p90 >= 1 APPROXIMATES
    # active_frac >= 0.1 (>=10% of DAM channels above floor) -- exact except for linear-interpolation
    # ties at the 10% boundary, where p90 is marginally stricter. Its value for the null is that it
    # scales LINEARLY under the 1/r^2 distance correction (percentile(c*r,90) = c*percentile(r,90)),
    # so the corrected detector is a self-consistent p90 rule applied identically at every range.
    snr_p90 = np.percentile(med / floor, 90, axis=1)
    jd_bin = jd[:n].reshape(-1, bin_s).mean(axis=1)
    return {"jd": jd_bin, "active_frac": active_frac, "snr_p90": snr_p90}


def detect_active(active_frac: np.ndarray, *, min_frac: float = 0.1) -> np.ndarray:
    """A time bin is 'DAM active' when >= ``min_frac`` of DAM-band channels exceed background."""
    return np.asarray(active_frac, float) >= min_frac


def sensitivity_corrected_active(
    snr_p90: np.ndarray, dist_au: np.ndarray, *, ref_au: float | None = None
) -> np.ndarray:
    r"""The 1/r^2 sensitivity null: distance-correct each bin's DAM SNR, then re-detect.

    DAM power falls as :math:`1/r^2` with Juno--Jupiter range, so near perijove the *same*
    intrinsic emission clears the background+5$\sigma$ floor more often --- a pure sensitivity
    effect, not intrinsic occurrence. Correcting each bin's 90th-percentile SNR to a reference
    range (:math:`S\to S\,(r/r_\mathrm{ref})^2`, default ``ref_au`` = median range) and
    re-thresholding at :math:`\ge 1` is the null: if the proximity duty-cycle trend is only
    sensitivity, the corrected occurrence is flat with range. Any residual near/far ratio after
    correction bounds the intrinsic$+$beaming part. Mirrors the `skr` slice's null model.
    """
    s = np.asarray(snr_p90, float)
    d = np.asarray(dist_au, float)
    ref = float(np.nanmedian(d)) if ref_au is None else ref_au
    return (s * (d / ref) ** 2) >= 1.0


def occurrence_map(
    cml_deg: np.ndarray,
    io_phase_deg: np.ndarray,
    active: np.ndarray,
    *,
    n_bins: int = 18,
    min_exposure: int = 3,
) -> dict:
    """Occurrence probability in the (CML, Io-phase) plane, with per-cell exposure.

    Cells visited fewer than ``min_exposure`` times are masked NaN (a one-month orbit does not
    cover the plane uniformly --- exposure must be reported, not hidden).
    """
    edges = np.linspace(0.0, 360.0, n_bins + 1)
    exp_map, _, _ = np.histogram2d(cml_deg, io_phase_deg, bins=[edges, edges])
    act_map, _, _ = np.histogram2d(cml_deg[active], io_phase_deg[active], bins=[edges, edges])
    with np.errstate(invalid="ignore", divide="ignore"):
        occ = act_map / exp_map
    occ[exp_map < min_exposure] = np.nan
    return {"occ": occ, "exposure": exp_map, "edges": edges}


#: Canonical Io-controlled regions in (CML, conventional Io phase), after Carr et al. (1983)
#: Fig. 7.32 / Marques et al. (2017) Table 2. CML ranges with c0 > c1 WRAP through 0 (Io-C).
IO_REGIONS = {
    "Io-A": ((200.0, 290.0), (205.0, 260.0)),
    "Io-B": ((90.0, 200.0), (80.0, 110.0)),
    "Io-C": ((300.0, 20.0), (225.0, 260.0)),
    "Io-D": ((0.0, 200.0), (95.0, 130.0)),
}


def _in_box(cml, pha, box):
    """Membership with CML wrap support (c0 > c1 wraps through 0 deg)."""
    (c0, c1), (p0, p1) = box
    in_c = ((cml >= c0) | (cml < c1)) if c0 > c1 else ((cml >= c0) & (cml < c1))
    return in_c & (pha >= p0) & (pha < p1)


def io_region_contrast(m: dict) -> dict:
    """Mean occurrence inside the canonical Io boxes vs outside --- the recover-a-known statistic."""
    edges = m["edges"]
    cen = 0.5 * (edges[:-1] + edges[1:])
    cml_g, pha_g = np.meshgrid(cen, cen, indexing="ij")
    inside = np.zeros_like(cml_g, bool)
    for box in IO_REGIONS.values():
        inside |= _in_box(cml_g, pha_g, box)
    occ = m["occ"]
    good = np.isfinite(occ)
    inside_occ = float(np.nanmean(occ[inside & good])) if (inside & good).any() else float("nan")
    outside_occ = float(np.nanmean(occ[~inside & good])) if (~inside & good).any() else float("nan")
    return {
        "occ_io_regions": inside_occ,
        "occ_elsewhere": outside_occ,
        "contrast": inside_occ / outside_occ if outside_occ > 0 else float("nan"),
        "cells_used": int(good.sum()),
    }


def sensitivity_censored_active(
    snr_p90: np.ndarray, dist_au: np.ndarray, *, far_pct: float = 95.0
) -> np.ndarray:
    r"""The common-sensitivity census: detections censored DOWN to the far-range threshold.

    ``sensitivity_corrected_active`` multiplies far-range SNRs by :math:`(d/d_\mathrm{ref})^2 > 1`
    and re-thresholds, which promotes sub-threshold far-range values --- noise included --- across
    the floor (measured: the far-quartile duty rose 0.07% -> 0.112% under it, with no committed
    false-positive rate; the sibling ``skr`` slice hit the same failure and its controls are in
    ``tests/test_skr.py``). This estimator only ever scales down: a bin counts when it is raw-active
    (``snr_p90 >= 1``) AND its SNR rescaled to the ``far_pct``-percentile range still clears 1, so
    every kept detection would have been visible from the far reference and noise is never promoted.
    """
    s = np.asarray(snr_p90, float)
    d = np.asarray(dist_au, float)
    good = np.isfinite(s) & np.isfinite(d)
    d_far = float(np.nanpercentile(d[good], far_pct)) if good.any() else float("nan")
    out = np.zeros(s.shape, bool)
    out[good] = (s[good] >= 1.0) & (s[good] * (d[good] / d_far) ** 2 >= 1.0)
    return out


def monthly_contrast_test(contrasts: list | np.ndarray) -> dict:
    """The one-sample test the abstract cites, computed and committed rather than asserted.

    A ratio is tested against unity on the log scale (a t-test on raw ratios against 1 is the
    wrong parameterisation), with a TWO-sided sign test alongside. Returns the per-month list,
    log-scale mean/sd/t, the two-sided sign-test p, and the 95% CI on the mean monthly contrast
    (log-scale, exponentiated) --- the number that shows what "does not reject" is worth: an
    interval spanning [0.5, 1.6] does not reject a 1.6x enhancement either.
    """
    c = np.asarray(contrasts, float)
    c = c[np.isfinite(c) & (c > 0)]
    n = int(c.size)
    if n < 2:
        return {"n": n, "contrasts": [round(float(x), 2) for x in c]}
    lg = np.log(c)
    mean, sd = float(lg.mean()), float(lg.std(ddof=1))
    se = sd / np.sqrt(n)
    t_stat = mean / se if se > 0 else float("nan")
    # two-sided sign test against unity: P(X <= k) + P(X >= n-k) under Binomial(n, 1/2)
    from math import comb

    k = int(np.sum(c > 1.0))
    lo = min(k, n - k)
    p_sign = min(1.0, 2.0 * sum(comb(n, i) for i in range(lo + 1)) / 2.0**n)
    # 95% CI on the mean log contrast, using the t quantile for small n
    from scipy import stats as _st

    tq = float(_st.t.ppf(0.975, n - 1))
    ci = (float(np.exp(mean - tq * se)), float(np.exp(mean + tq * se)))
    return {
        "n": n,
        "contrasts": [round(float(x), 2) for x in c],
        "log_mean": round(mean, 3),
        "log_sd": round(sd, 3),
        "t": round(t_stat, 2),
        "p_sign_two_sided": round(p_sign, 3),
        "geo_mean": round(float(np.exp(mean)), 2),
        "ci95": [round(ci[0], 2), round(ci[1], 2)],
    }


def box_shift_scan(
    cml: np.ndarray, pha: np.ndarray, active: np.ndarray, *, shifts=None
) -> list[dict]:
    """Io-box contrast under a rigid shift of the box set in Io phase: the frame-convention probe.

    The one systematic capable of manufacturing this census's null is a residual offset in the
    Io-phase convention (``io_phase`` mixes an eastward inertial mean longitude with a westward
    System III sub-observer longitude, and an earlier convention bug produced exactly a washed-out
    contrast). A synthetic round-trip cannot catch it: injection and recovery share the frame. This
    scan can: shift every canonical box rigidly in phase and recompute the contrast. If the true
    boxes organise the emission, the maximum sits at zero shift; a maximum at a large offset says
    the convention is wrong; a flat scan says the boxes do not organise this data from any offset.
    """
    if shifts is None:
        shifts = np.arange(-180.0, 181.0, 30.0)
    out = []
    for sh in shifts:
        m = occurrence_map(cml, (np.asarray(pha, float) - sh) % 360.0, active)
        con = io_region_contrast(m)
        out.append(
            {
                "shift_deg": float(sh),
                "contrast": round(con["contrast"], 2) if np.isfinite(con["contrast"]) else None,
            }
        )
    return out


def day_block_bootstrap_contrast(
    jd: np.ndarray,
    cml: np.ndarray,
    pha: np.ndarray,
    active: np.ndarray,
    *,
    n_boot: int = 200,
    seed: int = 0,
) -> dict:
    """A day-block bootstrap error for the Io-box contrast.

    The 15-s bins are not independent --- the committed map shows contiguous emission episodes ---
    so a per-bin error would be fiction. Days are the natural block (episodes do not span the
    ~9.9 h rotation many times within one). Resamples whole days with replacement and returns the
    bootstrap SE and the 2.5/97.5 percentiles of the contrast.
    """
    rng = np.random.default_rng(seed)
    day = np.floor(np.asarray(jd, float)).astype(int)
    udays = np.unique(day)
    idx_by_day = {d: np.where(day == d)[0] for d in udays}
    vals = []
    for _ in range(n_boot):
        pick = rng.choice(udays, size=udays.size, replace=True)
        idx = np.concatenate([idx_by_day[d] for d in pick])
        con = io_region_contrast(occurrence_map(cml[idx], pha[idx], active[idx]))
        if np.isfinite(con["contrast"]):
            vals.append(con["contrast"])
    if len(vals) < 10:
        return {"se": None, "ci95": None, "n_boot_ok": len(vals)}
    arr = np.asarray(vals)
    return {
        "se": round(float(arr.std(ddof=1)), 2),
        "ci95": [
            round(float(np.percentile(arr, 2.5)), 2),
            round(float(np.percentile(arr, 97.5)), 2),
        ],
        "n_boot_ok": len(vals),
    }


def episode_stats(active: np.ndarray) -> dict:
    """Run-length statistics of the active flag: the effective-N behind every error bar."""
    a = np.asarray(active, bool).astype(int)
    d = np.diff(np.r_[0, a, 0])
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    lengths = ends - starts
    return {
        "n_episodes": int(lengths.size),
        "median_episode_bins": int(np.median(lengths)) if lengths.size else 0,
        "max_episode_bins": int(lengths.max()) if lengths.size else 0,
    }


def per_region_contrast(m: dict) -> dict:
    """The contrast per canonical Io box, not just the union: averaging an enhanced Io-B against
    boxes on the dark side of the plane's structure can return ~1 while a per-region signal is
    present."""
    edges = m["edges"]
    cen = 0.5 * (edges[:-1] + edges[1:])
    cml_g, pha_g = np.meshgrid(cen, cen, indexing="ij")
    occ = m["occ"]
    good = np.isfinite(occ)
    out = {}
    outside = np.ones_like(cml_g, bool)
    for box in IO_REGIONS.values():
        outside &= ~_in_box(cml_g, pha_g, box)
    base = float(np.nanmean(occ[outside & good])) if (outside & good).any() else float("nan")
    for name, box in IO_REGIONS.items():
        inside = _in_box(cml_g, pha_g, box)
        v = float(np.nanmean(occ[inside & good])) if (inside & good).any() else float("nan")
        out[name] = round(v / base, 2) if (np.isfinite(v) and base > 0) else None
    return out


def synthetic_orbit(
    n_days: float = 28.0,
    *,
    bin_s: int = 15,
    p_in: float = 0.35,
    p_out: float = 0.04,
    seed: int = 0,
) -> dict:
    """Synthetic month: real CML/Io-phase rates, DAM active preferentially in the Io boxes.

    CML advances at Jupiter's System III rate (870.536 deg/day; Archinal et al. 2018) and Io's
    longitude at the Lieske rate, so the (CML, phase) plane fills exactly as it does for a real
    observer. Activity is Bernoulli: ``p_in`` inside the canonical boxes, ``p_out`` outside ---
    the census must recover a contrast ~``p_in/p_out``.
    """
    rng = np.random.default_rng(seed)
    jd = J2000_JD + np.arange(0.0, n_days, bin_s / 86400.0)
    cml = (284.95 + 870.5360000 * (jd - J2000_JD)) % 360.0
    pha = io_phase(jd, cml)
    inside = np.zeros(jd.size, bool)
    for box in IO_REGIONS.values():
        inside |= _in_box(cml, pha, box)
    active = rng.random(jd.size) < np.where(inside, p_in, p_out)
    return {"jd": jd, "cml": cml, "io_phase": pha, "active": active, "p_in": p_in, "p_out": p_out}


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: synthetic-orbit recover-a-known. Real: one month of Juno/Waves CDFs + Horizons."""

    if offline:
        s = synthetic_orbit()
        cml, pha, active = s["cml"], s["io_phase"], s["active"]
        source = "synthetic orbit (canonical Io boxes injected)"
        expected = s["p_in"] / s["p_out"]
        # The recovered-vs-injected curve, including contrasts near unity -- the regime the real
        # measurement (1.12) lives in and the single 8.75-point calibration said nothing about.
        # The estimator's boundary-cell dilution contracts every point toward 1; committing the
        # curve makes the contraction a measured correction instead of an unstated bias.
        curve = []
        for p_in_c, p_out_c in ((0.05, 0.04), (0.06, 0.04), (0.08, 0.04), (0.35, 0.04)):
            sc = synthetic_orbit(p_in=p_in_c, p_out=p_out_c, seed=1)
            cc = io_region_contrast(occurrence_map(sc["cml"], sc["io_phase"], sc["active"]))
            curve.append(
                {
                    "injected": round(p_in_c / p_out_c, 2),
                    "recovered": round(cc["contrast"], 2) if np.isfinite(cc["contrast"]) else None,
                }
            )
        extra_syn = {"recovery_curve": curve}
    else:  # pragma: no cover - data files + network
        files = sorted(DATA_DIR.glob("jno_wav_cdr_lesia_*_v0?.cdf"))
        parts = [read_waves_cdf(f) for f in files]
        months = [f.name.split("_")[4][:6] for f in files]
        jd = np.concatenate([p["jd"] for p in parts])
        af = np.concatenate([p["active_frac"] for p in parts])
        snr_p90 = np.concatenate([p["snr_p90"] for p in parts])
        # fetch Horizons per contiguous data segment (a >2-day gap starts a new segment),
        # so multi-month runs stay inside per-query epoch limits
        order = np.argsort(jd)
        jd, af, snr_p90 = jd[order], af[order], snr_p90[order]
        starts = [0] + list(np.where(np.diff(jd) > 2.0)[0] + 1) + [jd.size]
        cml = np.empty_like(jd)
        dist = np.empty_like(jd)
        for a, b in zip(starts[:-1], starts[1:], strict=True):
            eph = fetch_cml_horizons(float(jd[a]) - 0.02, float(jd[b - 1]) + 0.02)
            cml_unwrap = np.unwrap(np.radians(eph["cml"]))
            cml[a:b] = np.degrees(np.interp(jd[a:b], eph["jd"], cml_unwrap)) % 360.0
            dist[a:b] = np.interp(jd[a:b], eph["jd"], eph["delta_au"])
        pha = io_phase(jd, cml)
        active = detect_active(af)
        active_p90 = snr_p90 >= 1.0  # the p90 surrogate detector (uncorrected baseline)
        # PRIMARY null: the censored (downward-only) census, which cannot promote noise. The
        # upward-rescale variant is kept as a recorded comparison; it promotes sub-threshold
        # far-range values across the floor (its far-quartile duty EXCEEDS the raw one).
        active_cens = sensitivity_censored_active(snr_p90, dist)
        active_corr = sensitivity_corrected_active(snr_p90, dist)  # superseded, recorded
        source = f"Juno/Waves L3a v01+v02, {len(files)} days"
        expected = float("nan")
        # the vantage dimension: Juno--Jupiter range dominates detection (proximity, not clock)
        far = dist > np.median(dist)
        m_far = occurrence_map(cml[far], pha[far], active[far])
        con_far = io_region_contrast(m_far)
        extra = {
            "activity_near_half_pct": round(100 * float(active[~far].mean()), 2),
            "activity_far_half_pct": round(100 * float(active[far].mean()), 2),
            "io_contrast_far_half": round(con_far["contrast"], 2)
            if np.isfinite(con_far["contrast"])
            else None,
        }
        # distance-RESOLVED Io contrast (the paper's scoped test): quartiles of Juno range
        qs = np.quantile(dist, [0.25, 0.5, 0.75])
        qmasks = [
            dist <= qs[0],
            (dist > qs[0]) & (dist <= qs[1]),
            (dist > qs[1]) & (dist <= qs[2]),
            dist > qs[2],
        ]
        for k, msk in enumerate(qmasks, start=1):
            mq = occurrence_map(cml[msk], pha[msk], active[msk])
            cq = io_region_contrast(mq)
            extra[f"io_contrast_q{k}"] = (
                round(cq["contrast"], 2) if np.isfinite(cq["contrast"]) else None
            )
            extra[f"activity_q{k}_pct"] = round(100 * float(active[msk].mean()), 2)
            # like-for-like: the SAME p90 detector raw, then under each null
            extra[f"activity_q{k}_raw_p90_pct"] = round(100 * float(active_p90[msk].mean()), 3)
            extra[f"activity_q{k}_cens_pct"] = round(100 * float(active_cens[msk].mean()), 3)
            extra[f"activity_q{k}_corr_pct"] = round(100 * float(active_corr[msk].mean()), 3)
            extra[f"n_active_q{k}"] = int(active[msk].sum())
            extra[f"cells_used_q{k}"] = cq["cells_used"]
        # raw vs corrected near/far ratios, all on the SAME p90 detector so the correction is
        # isolated from a detector swap. NOTE the two raw baselines (active_frac vs p90) differ
        # by ~70% (196 vs 330) -- they are different detectors and are reported as such.
        near_raw, far_raw = float(active[qmasks[0]].mean()), float(active[qmasks[3]].mean())
        near_p, far_p = float(active_p90[qmasks[0]].mean()), float(active_p90[qmasks[3]].mean())
        near_c, far_c = float(active_corr[qmasks[0]].mean()), float(active_corr[qmasks[3]].mean())
        near_z, far_z = float(active_cens[qmasks[0]].mean()), float(active_cens[qmasks[3]].mean())
        extra["near_far_raw"] = round(near_raw / far_raw, 1) if far_raw > 0 else None
        extra["near_far_raw_p90"] = round(near_p / far_p, 1) if far_p > 0 else None
        extra["near_far_corrected"] = round(near_c / far_c, 2) if far_c > 0 else None
        extra["near_far_censored"] = round(near_z / far_z, 2) if far_z > 0 else None
        extra["range_near_rj"] = round(float(np.median(dist[qmasks[0]])) / RJ_AU, 1)
        extra["range_far_rj"] = round(float(np.median(dist[qmasks[3]])) / RJ_AU, 1)
        # per-month Io contrast: the robustness spread across orbital configurations
        pm = []
        month_rows = []
        versions = [f.name.split("_")[-1].split(".")[0] for f in files]
        for ym in sorted(set(months)):
            sel = np.zeros(jd.size, bool)
            vset = set()
            n_days_m = 0
            for f_m, ver, p in zip(months, versions, parts, strict=True):
                if f_m == ym:
                    sel[np.searchsorted(jd, p["jd"])] = True
                    vset.add(ver)
                    n_days_m += 1
            cm = io_region_contrast(occurrence_map(cml[sel], pha[sel], active[sel]))
            row = {
                "month": ym,
                "n_days": n_days_m,
                "n_bins": int(sel.sum()),
                "n_active": int(active[sel].sum()),
                "contrast": round(cm["contrast"], 2) if np.isfinite(cm["contrast"]) else None,
                "versions": sorted(vset),
            }
            month_rows.append(row)
            if np.isfinite(cm["contrast"]):
                pm.append(cm["contrast"])
        extra["per_month"] = month_rows
        if pm:
            extra["n_months"] = len(pm)
            extra["io_contrast_month_median"] = round(float(np.median(pm)), 2)
            extra["io_contrast_month_min"] = round(float(np.min(pm)), 2)
            extra["io_contrast_month_max"] = round(float(np.max(pm)), 2)
            # The one-sample test the abstract cites, committed rather than asserted, with the
            # CI that shows what "does not reject" is worth.
            extra["monthly_contrast_test"] = monthly_contrast_test(pm)
            # ...and the same statistics with the single v01 month excluded, so "the version mix
            # is benign" is a measurement instead of a hypothesis.
            pm_v02 = [
                r["contrast"]
                for r in month_rows
                if r["contrast"] is not None and r["versions"] == ["v02"]
            ]
            if len(pm_v02) >= 2:
                extra["monthly_contrast_test_v02_only"] = monthly_contrast_test(pm_v02)
        # the frame-convention probe, the day-block error bar, the effective N, per-region
        extra["box_shift_scan"] = box_shift_scan(cml, pha, active)
        extra["contrast_day_bootstrap"] = day_block_bootstrap_contrast(jd, cml, pha, active)
        extra["episodes"] = episode_stats(active)
        extra["per_region_contrast"] = per_region_contrast(occurrence_map(cml, pha, active))

    m = occurrence_map(cml, pha, active)
    con = io_region_contrast(m)
    metrics = {
        "source": source,
        "n_bins_total": int(active.size),
        "n_active": int(active.sum()),
        "occ_io_regions": round(con["occ_io_regions"], 3),
        "occ_elsewhere": round(con["occ_elsewhere"], 3),
        "io_contrast": round(con["contrast"], 2) if np.isfinite(con["contrast"]) else None,
        "expected_contrast": round(expected, 2) if np.isfinite(expected) else None,
        "cells_used": con["cells_used"],
    }
    if offline:
        metrics.update(extra_syn)
    else:  # pragma: no cover - real-leg extras
        metrics.update(extra)
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "junodam_metrics.json")
    _figure(m, cml, pha, active, op / "papers" / "junodam" / "figures")
    _write_macros(metrics, op / "papers" / "junodam" / "generated" / "macros.tex")
    return metrics


def _figure(m: dict, cml, pha, active, out_dir) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0))
    im = ax1.pcolormesh(m["edges"], m["edges"], m["occ"].T, cmap="viridis", shading="auto")
    fig.colorbar(im, ax=ax1, label="occurrence probability")
    for name, ((c0, c1), (p0, p1)) in IO_REGIONS.items():
        spans = [(c0, c1)] if c0 < c1 else [(c0, 360.0), (0.0, c1)]  # wrap-aware (Io-C)
        for a, b in spans:
            ax1.add_patch(
                plt.Rectangle((a, p0), b - a, p1 - p0, fill=False, ec="w", lw=1.0, ls="--")
            )
        ax1.text(spans[0][0] + 3, p1 - 12, name, color="w", fontsize=7)
    ax1.set(xlabel="CML (System III, deg)", ylabel="Io phase (deg)", title="DAM occurrence map")
    ax2.plot(cml[~active][::7], pha[~active][::7], ".", ms=1, color="0.8")
    ax2.plot(cml[active], pha[active], ".", ms=2, color="C3")
    ax2.set(
        xlabel="CML (deg)",
        ylabel="Io phase (deg)",
        title="Active bins (red) vs coverage",
        xlim=(0, 360),
        ylim=(0, 360),
    )
    fig.tight_layout()
    fig.savefig(out / "junodam.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    """Emit the contrast macros in BOTH namespaces; only the active mode's carry values.

    ``io_contrast`` and ``expected_contrast`` mean different things in the two run modes: on
    the real leg the contrast is the measured Io-region enhancement (1.12), on the offline leg
    it is the recovery of an injected one (~7.2 for an injected 8.75). They shared one macro
    name, so an offline rebuild would have written the synthetic recovery into the macro the
    paper uses for the real measurement -- the documented ``\tiiNEvents`` clobber, which
    ``preserve_live_macros`` cannot arbitrate because both runs write real values.
    """

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    _real = not str(m.get("source", "")).lower().startswith("synthetic")
    lines = [
        "% Auto-generated by jansky_research.junodam._write_macros -- do not edit.",
        rf"\newcommand{{\jdSource}}{{{m['source']}}}",
        # The census counts and occupancies are mode-dependent too (the 2026-08-24 recurrence:
        # a synthetic re-run flipped \jdNbins to the 28-day fixture's 161280 and the abstract's
        # first sentence rendered fixture counts for the 210-day census). Both namespaces are
        # emitted every run -- the merge accumulates values, not names.
        rf"\newcommand{{\jdRealNbins}}{{{_fmt('n_bins_total') if _real else '--'}}}",
        rf"\newcommand{{\jdSynNbins}}{{{_fmt('n_bins_total') if not _real else '--'}}}",
        rf"\newcommand{{\jdRealNact}}{{{_fmt('n_active') if _real else '--'}}}",
        rf"\newcommand{{\jdSynNact}}{{{_fmt('n_active') if not _real else '--'}}}",
        rf"\newcommand{{\jdRealOccIo}}{{{_fmt('occ_io_regions') if _real else '--'}}}",
        rf"\newcommand{{\jdSynOccIo}}{{{_fmt('occ_io_regions') if not _real else '--'}}}",
        rf"\newcommand{{\jdRealOccOut}}{{{_fmt('occ_elsewhere') if _real else '--'}}}",
        rf"\newcommand{{\jdSynOccOut}}{{{_fmt('occ_elsewhere') if not _real else '--'}}}",
        rf"\newcommand{{\jdRealCells}}{{{_fmt('cells_used') if _real else '--'}}}",
        rf"\newcommand{{\jdSynCells}}{{{_fmt('cells_used') if not _real else '--'}}}",
        # Mode-dependent: see the docstring. Only the active mode's namespace is filled.
        rf"\newcommand{{\jdRealContrast}}{{{_fmt('io_contrast') if _real else '--'}}}",
        rf"\newcommand{{\jdSynContrast}}{{{_fmt('io_contrast') if not _real else '--'}}}",
        rf"\newcommand{{\jdSynExpContrast}}{{{_fmt('expected_contrast') if not _real else '--'}}}",
        rf"\newcommand{{\jdActNear}}{{{_fmt('activity_near_half_pct')}}}",
        rf"\newcommand{{\jdActFar}}{{{_fmt('activity_far_half_pct')}}}",
        rf"\newcommand{{\jdContrastFar}}{{{_fmt('io_contrast_far_half')}}}",
        rf"\newcommand{{\jdCqA}}{{{_fmt('io_contrast_q1')}}}",
        rf"\newcommand{{\jdCqB}}{{{_fmt('io_contrast_q2')}}}",
        rf"\newcommand{{\jdCqC}}{{{_fmt('io_contrast_q3')}}}",
        rf"\newcommand{{\jdCqD}}{{{_fmt('io_contrast_q4')}}}",
        rf"\newcommand{{\jdAqA}}{{{_fmt('activity_q1_pct')}}}",
        rf"\newcommand{{\jdAqD}}{{{_fmt('activity_q4_pct')}}}",
        rf"\newcommand{{\jdAqAcorr}}{{{_fmt('activity_q1_corr_pct')}}}",
        rf"\newcommand{{\jdAqDcorr}}{{{_fmt('activity_q4_corr_pct')}}}",
        rf"\newcommand{{\jdNearFarRaw}}{{{_fmt('near_far_raw')}}}",
        rf"\newcommand{{\jdNearFarRawPninety}}{{{_fmt('near_far_raw_p90')}}}",
        rf"\newcommand{{\jdNearFarCorr}}{{{_fmt('near_far_corrected')}}}",
        rf"\newcommand{{\jdAqBcorr}}{{{_fmt('activity_q2_corr_pct')}}}",
        rf"\newcommand{{\jdAqCcorr}}{{{_fmt('activity_q3_corr_pct')}}}",
        rf"\newcommand{{\jdRangeNear}}{{{_fmt('range_near_rj')}}}",
        rf"\newcommand{{\jdRangeFar}}{{{_fmt('range_far_rj')}}}",
        rf"\newcommand{{\jdNmonths}}{{{_fmt('n_months')}}}",
        rf"\newcommand{{\jdCmMed}}{{{_fmt('io_contrast_month_median')}}}",
        rf"\newcommand{{\jdCmMin}}{{{_fmt('io_contrast_month_min')}}}",
        rf"\newcommand{{\jdCmMax}}{{{_fmt('io_contrast_month_max')}}}",
    ]
    # Revision evidence (2026-08-24): the committed test, its CI, the day-block error, the
    # censored census, per-region contrasts, the episode count, and the near-unity injection.
    mt = m.get("monthly_contrast_test") or {}
    ci = mt.get("ci95") or [None, None]
    boot = m.get("contrast_day_bootstrap") or {}
    bci = boot.get("ci95") or [None, None]
    reg = m.get("per_region_contrast") or {}
    epi = m.get("episodes") or {}
    curve = {c["injected"]: c["recovered"] for c in m.get("recovery_curve", [])}
    derived = [
        ("jdMonthSignP", mt.get("p_sign_two_sided")),
        ("jdMonthT", mt.get("t")),
        ("jdMonthGeoMean", mt.get("geo_mean")),
        ("jdMonthCIlo", ci[0]),
        ("jdMonthCIhi", ci[1]),
        ("jdContrastSe", boot.get("se")),
        ("jdContrastCIlo", bci[0]),
        ("jdContrastCIhi", bci[1]),
        ("jdNearFarCens", m.get("near_far_censored")),
        ("jdNEpisodes", epi.get("n_episodes")),
        ("jdRegA", reg.get("Io-A")),
        ("jdRegB", reg.get("Io-B")),
        ("jdRegC", reg.get("Io-C")),
        ("jdRegD", reg.get("Io-D")),
        ("jdSynRecNearUnity", curve.get(1.25)),
        ("jdSynRecMid", curve.get(2.0)),
    ]
    for name, v in derived:
        lines.append(rf"\newcommand{{\{name}}}{{{'--' if v is None else v}}}")
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

    p = argparse.ArgumentParser(description="Juno/Waves DAM occurrence census (CML x Io phase).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
