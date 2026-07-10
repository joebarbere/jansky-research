"""Tests for jansky_research.typeii -- OVRO-LWA type II burst detector + census. Offline."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import typeii as t2
from jansky_research.solarbursts import synthetic_burst


def test_synthetic_typeii_drifts_slowly():
    s = t2.synthetic_typeii(shock_speed_kms=1000.0, seed=0)
    # a CME-shock drift is slow: well inside the type II band, orders below type III
    assert t2.DRIFT_SLOW_MIN < abs(s["truth_drift_mhz_s"]) < t2.DRIFT_SLOW_MAX
    assert s["data"].shape == (s["freqs"].size, s["times"].size)


def test_detect_typeii_flags_slow_drift():
    s = t2.synthetic_typeii(seed=1)
    r = t2.detect_typeii(s["data"], s["freqs"], s["times"])
    assert r["klass"] == "type_II" and r["detected"]
    assert r["drift_mhz_s"] < 0 and abs(r["drift_mhz_s"]) < t2.DRIFT_SLOW_MAX
    assert r["duration_s"] >= t2.DURATION_MIN_S
    assert r["coherence"] >= t2.COHERENCE_MIN
    assert r["harmonic_score"] > 0.3  # fundamental+harmonic injected


def test_detect_rejects_fast_type_iii():
    s3 = synthetic_burst(speed_c=0.3, seed=2)  # fast electron beam
    r = t2.detect_typeii(s3["data"], s3["freqs"], s3["times"])
    assert r["klass"] == "type_III" and not r["detected"]


def test_detect_rejects_noise_and_rfi():
    # pure noise + a few constant-offset RFI channels: no coherent slow drift -> not type II
    rng = np.random.default_rng(0)
    f = np.linspace(88, 20, 240)
    d = rng.normal(0, 1.0, (240, 600))
    for ch in rng.choice(240, 4, replace=False):
        d[ch] += 8.0
    r = t2.detect_typeii(d, f, np.linspace(0, 300, 600))
    assert not r["detected"]
    assert r["coherence"] < t2.COHERENCE_MIN  # scattered peaks, no coherent ridge


def test_rfi_mask_flags_persistent_channels():
    rng = np.random.default_rng(3)
    d = rng.normal(0, 1.0, (100, 400))
    # 5 channels persistently bright (time-variable, so background subtraction keeps them)
    rfi_ch = [10, 20, 30, 40, 50]
    for ch in rfi_ch:
        d[ch] += rng.uniform(5, 15, 400)
    keep = t2.rfi_mask(d)
    assert not keep[rfi_ch].any()  # all RFI channels masked
    assert keep.sum() > 90  # clean channels kept


def test_classify_burst_boundaries():
    # slow + long + coherent + dense -> type II
    assert t2.classify_burst(-0.1, 300, 0.7, 0.9, 100) == "type_II"
    # fast + coherent -> type III
    assert t2.classify_burst(-8.0, 10, 0.5, 0.9, 50) == "type_III"
    # slow but incoherent (noise) -> none
    assert t2.classify_burst(-0.1, 300, 0.1, 0.2, 100) == "none"
    # slow but short -> none
    assert t2.classify_burst(-0.1, 30, 0.7, 0.9, 100) == "none"
    # near-zero drift (horizontal RFI) -> none
    assert t2.classify_burst(-0.005, 300, 0.0, 0.9, 100) == "none"
    # positive drift (rising) -> none (type II drifts down)
    assert t2.classify_burst(0.1, 300, 0.7, 0.9, 100) == "none"


def test_harmonic_score_detects_pair():
    s = t2.synthetic_typeii(with_harmonic=True, band_split=False, seed=4)
    keep = t2.rfi_mask(s["data"])
    rt, rf = t2.track_drift_ridge(s["data"], s["freqs"], s["times"], keep=keep)
    with_h = t2.harmonic_score(s["data"], s["freqs"], rt, rf, s["times"])
    s0 = t2.synthetic_typeii(with_harmonic=False, band_split=False, seed=4)
    keep0 = t2.rfi_mask(s0["data"])
    rt0, rf0 = t2.track_drift_ridge(s0["data"], s0["freqs"], s0["times"], keep=keep0)
    no_h = t2.harmonic_score(s0["data"], s0["freqs"], rt0, rf0, s0["times"])
    assert with_h > no_h  # the fundamental+harmonic pair scores higher than a single lane


def test_crossmatch_cme_nearest_preceding():
    cme = [
        {"onset_hr": 10.0, "speed_kms": 1200, "width_deg": 120},  # true driver, 0.2 hr before
        {"onset_hr": 5.0, "speed_kms": 1800, "width_deg": 200},  # faster but far earlier (decoy)
        {"onset_hr": 50.0, "speed_kms": 900, "width_deg": 80},  # unrelated
    ]
    m = t2.crossmatch_cme(10.2, cme)
    assert m is not None and m["onset_hr"] == 10.0  # nearest-preceding, not the fastest
    # a burst with no preceding CME in window -> None
    assert t2.crossmatch_cme(200.0, cme) is None


def test_cme_association_recovers_fast_wide_bias():
    matched = [
        {"speed_kms": 1200, "width_deg": 120},
        {"speed_kms": 1000, "width_deg": 90},
        {"speed_kms": 400, "width_deg": 30},  # one slow/narrow
        None,  # unmatched burst
    ]
    a = t2.cme_association_fraction(matched)
    assert a["n_matched"] == 3
    assert a["frac_fast"] > 0.5 and a["frac_wide"] > 0.5  # fast-and-wide bias
    assert a["median_speed_kms"] == 1000.0


def test_detects_single_lane_and_marginal_snr():
    # a single-lane type II (no harmonic) still classifies type II (harmonic not required)
    s1 = t2.synthetic_typeii(with_harmonic=False, seed=6)
    assert t2.detect_typeii(s1["data"], s1["freqs"], s1["times"])["klass"] == "type_II"
    # the completeness curve is monotone and non-saturated at low SNR (an honest performance floor)
    m = t2.run("/tmp/t2_snrtest", offline=True)
    c = m["completeness_vs_snr"]
    assert c["snr2"] < c["snr3"] <= c["snr4"]  # rises with SNR; weak SNR genuinely fails some
    assert c["snr2"] < 0.9  # near-threshold is NOT perfectly recovered


def test_parse_lwa_dspec_reads_real_4d_format(tmp_path):
    from astropy.io import fits

    # the REAL AWS-Open-Data layout (confirmed on 20240514.fits): 4D primary array with FITS axes
    # (time, freq, 1, stokes) -> numpy (stokes, 1, freq, time), freq from FREQMIN/FREQMAX (GHz)
    n_f, n_t = 30, 50
    arr = np.random.default_rng(0).normal(0, 1, (2, 1, n_f, n_t))  # (stokes, 1, freq, time)
    arr[0, 0, 15, :] += 20.0  # bright Stokes-I channel to check orientation survives
    hdr = fits.Header({"FREQMIN": 0.0150, "FREQMAX": 0.0849, "BUNIT": "jy"})  # GHz
    p = tmp_path / "20240514.fits"
    fits.PrimaryHDU(arr, header=hdr).writeto(p)
    d = t2.parse_lwa_dspec(p)
    assert d["data"].shape == (n_f, n_t)  # Stokes I, (freq, time)
    assert d["freqs"][0] > d["freqs"][-1]  # descending MHz (as solarbursts expects)
    assert 14.5 < d["freqs"].min() < 15.5 and 84.0 < d["freqs"].max() < 86.0  # GHz->MHz
    assert d["times"].size == n_t


def test_parse_lwa_dspec_reads_table_fallback(tmp_path):
    from astropy.io import fits

    # the older documented table layout: 2D spectrum + 1D freq (Hz) + 1D time HDUs
    n_f, n_t = 30, 50
    spec = np.random.default_rng(1).normal(0, 1, (n_f, n_t))
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(spec, name="DSPEC"),
            fits.ImageHDU(np.linspace(13.4e6, 86.9e6, n_f), name="FREQ"),
            fits.ImageHDU(np.arange(n_t) * 0.256, name="TIME"),
        ]
    ).writeto(tmp_path / "d.fits")
    d = t2.parse_lwa_dspec(tmp_path / "d.fits")
    assert d["data"].shape == (n_f, n_t)
    assert d["freqs"][0] > d["freqs"][-1] and 13.0 < d["freqs"].min() < 14.0  # Hz->MHz, descending


def test_s3_dspec_url():
    u = t2.s3_dspec_url("2024-05-14")
    assert u == "https://ovro-lwa-solar.s3-us-west-2.amazonaws.com/spec_fits/2024/20240514.fits"


def test_downsample_time_block_averages():
    spec = np.arange(4 * 20, dtype=float).reshape(4, 20)  # (freq=4, time=20)
    out = t2._downsample_time(spec, 4)
    assert out.shape == (4, 5)  # 20 time cols -> 5 bins of 4
    assert np.allclose(out[0], [1.5, 5.5, 9.5, 13.5, 17.5])  # mean of each 4-block
    assert t2._downsample_time(spec, 1) is spec  # factor 1 is a no-op
    # a ragged length is truncated to a whole number of blocks
    assert t2._downsample_time(np.ones((2, 23)), 4).shape == (2, 5)


def test_sweep_day_finds_windowed_type_ii():
    # build a day-long spectrum (descending freq) with ONE injected type II in a 15-min window
    s = t2.synthetic_typeii(seed=2, duration_s=600.0, n_time=150)  # a burst window
    n_f = s["data"].shape[0]
    rng = np.random.default_rng(0)
    day = rng.normal(0, 1.0, (n_f, 1500))  # ~quiet day, same freq grid
    day[:, 300:450] += s["data"]  # drop the burst in at ~one window
    times = np.arange(1500) * 4.0  # 4 s bins
    dets = t2.sweep_day(day, s["freqs"], times, window_s=900.0, step_s=450.0)
    assert len(dets) >= 1  # the injected type II is found in its window
    assert all(d["klass"] == "type_II" for d in dets)
    # tagged with a window-centre time (seconds from day start) for the absolute cross-match
    assert all("t_center_s" in d for d in dets)
    assert all(1200 <= d["t_center_s"] <= 1800 for d in dets)  # the burst sits at 300-450 * 4 s


def test_sweep_day_empty_on_quiet_day():
    rng = np.random.default_rng(1)
    quiet = rng.normal(0, 1.0, (120, 1500))
    times = np.arange(1500) * 4.0
    assert t2.sweep_day(quiet, np.linspace(85, 15, 120), times) == []


def test_iso_to_hours_common_clock():
    assert t2._iso_to_hours("2024-01-01T00:00:00") == 0.0
    assert t2._iso_to_hours("2024-01-02T06:00:00") == 30.0  # 1 day + 6 h
    assert abs(t2._iso_to_hours("2024-05-14T13:45:00.096") - (134 * 24 + 13.75)) < 0.01


def test_parse_cme_html():
    # a minimal CDAW-format fragment: date time CPA width speed ...
    html = (
        "<td>2024/05/14</td> <td>02:24:05</td> <td>261</td> <td>26</td> <td>1200</td> <td>x</td>"
        "<tr><td>2024/05/14</td><td>10:00:00</td><td>Halo</td><td>360</td><td>1800</td>"
        "<td>2024/05/14</td><td>11:00:00</td><td>90</td><td>40</td><td>----</td>"  # unmeasured -> skipped
    )
    cmes = t2.parse_cme_html(html)
    assert len(cmes) == 2  # the ---- speed row is dropped
    assert cmes[0]["width_deg"] == 26.0 and cmes[0]["speed_kms"] == 1200.0
    assert cmes[1]["width_deg"] == 360.0 and cmes[1]["speed_kms"] == 1800.0  # halo parsed
    assert cmes[0]["onset_hr"] == t2._iso_to_hours("2024-05-14T02:24:05")


def test_parse_hek_flares():
    payload = {
        "result": [
            {"frm_name": "SWPC", "event_peaktime": "2024-05-14T02:09:00", "fl_goescls": "X1.7"},
            {"frm_name": "SWPC", "event_peaktime": "2024-05-14T15:28:00", "fl_goescls": "C9.6"},
            {
                "frm_name": "SWPC",
                "event_peaktime": "2024-05-14T20:00:00",
                "fl_goescls": "",
            },  # dropped
        ]
    }
    fl = t2.parse_hek_flares(payload)
    assert len(fl) == 2 and fl[0]["goes_class"] == "X1.7"
    assert abs(fl[0]["goes_log_flux"] - (-4 + np.log10(1.7))) < 0.01  # X1.7 -> ~-3.77
    assert fl[0]["peak_hr"] == t2._iso_to_hours("2024-05-14T02:09:00")


def test_crossmatch_flare_nearest_preceding():
    flares = [
        {"peak_hr": 100.0, "goes_class": "M2.0"},  # 0.2 h before the burst
        {"peak_hr": 95.0, "goes_class": "X5.0"},  # bigger but far earlier
    ]
    m = t2.crossmatch_flare(100.2, flares)
    assert m is not None and m["goes_class"] == "M2.0"  # nearest, not biggest
    assert t2.crossmatch_flare(200.0, flares) is None  # nothing in window


def test_purity_diagnostics_flags_background_like_matches():
    # a dense background CME population (mostly slow) + detections at random times: the CME matches
    # should look like the background (not fast-and-wide) and no better than chance -> flagged
    rng = np.random.default_rng(0)
    cme = [
        {"onset_hr": float(h), "speed_kms": float(rng.normal(400, 150)), "width_deg": 40.0}
        for h in rng.uniform(0, 5000, 3000)  # ~dense, slow background
    ]
    # detections at random times, with window-filling durations + drift uncorrelated to any CME
    dets = [
        {"burst_hr": float(h), "drift_mhz_s": float(rng.uniform(-0.05, -0.01)), "duration_s": 897.0}
        for h in rng.uniform(0, 5000, 200)
    ]
    d = t2.purity_diagnostics(dets, cme, n_days=200)
    assert d["association_is_background_like"] is True
    assert d["matched_cme_median_kms"] < 1.5 * d["bg_cme_median_kms"]
    assert d["observed_cme_match_rate"] <= d["chance_cme_match_rate"] + 0.05
    assert abs(d["drift_cme_speed_corr"]) < 0.2  # radio drift knows nothing about the CME speed
    assert d["window_saturation_frac"] == 1.0  # every ridge fills the window (FP signature)

    # a REAL type II sample: detections each just after a FAST wide CME -> not background-like
    fast_cme = [
        {"onset_hr": float(h), "speed_kms": 1500.0, "width_deg": 300.0} for h in range(0, 400, 20)
    ]
    real_dets = [{"burst_hr": float(h) + 0.3} for h in range(0, 400, 20)]
    d2 = t2.purity_diagnostics(real_dets, fast_cme, n_days=20)
    assert (
        d2["matched_cme_median_kms"] > 1.5 * d2["bg_cme_median_kms"]
        or not d2["association_is_background_like"]
    )


def test_occurrence_vs_phase_reuses_coverage_correction():
    # detections and observing days across 6 months; rate should track the injected sunspot number
    from datetime import datetime

    ep = datetime.strptime(t2.EPOCH_ISO, "%Y-%m-%dT%H:%M:%S")

    def hr(y, mo, d):
        return (datetime(y, mo, d) - ep).total_seconds() / 3600.0

    months = [(2024, m) for m in range(1, 7)]
    sunspots = {f"2024-{m:02d}": 50.0 + 30.0 * m for _, m in months}  # rising activity
    obs, bursts = [], []
    for _, mo in months:
        for d in range(1, 21):  # 20 observing days/month (flat coverage)
            obs.append(hr(2024, mo, d))
        for k in range(mo):  # more bursts in higher-activity months
            bursts.append(hr(2024, mo, 1 + k))
    occ = t2.occurrence_vs_phase(bursts, obs, sunspots)
    assert occ["n_months"] == 6
    assert occ["pearson_r"] > 0.9  # coverage-corrected rate correlates with the sunspot number


def test_run_offline_completeness_purity_and_bias(tmp_path):
    m = t2.run(str(tmp_path), offline=True)
    assert m["completeness"] >= 0.8 and m["purity"] >= 0.9  # honest mixed-difficulty number
    assert m["n_typeii_detected"] > 0
    # the recovered CME association ECHOES the injected fast-and-wide bias (a wiring check)
    assert m["assoc_frac_fast"] >= 0.6 and m["assoc_frac_wide"] >= 0.6
    saved = json.loads((tmp_path / "results" / "typeii_metrics.json").read_text())
    assert saved["completeness"] == m["completeness"]
    assert (tmp_path / "papers" / "typeii" / "figures" / "typeii.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "typeii" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\tiiRealFracFast}{--}" in macros  # real namespace placeholder offline


def test_write_macros_dual_namespace(tmp_path):
    p = tmp_path / "m.tex"
    t2._write_macros(
        {
            "source": "x",
            "is_real": True,
            "n_events": 10,
            "completeness": 0.9,
            "purity": 0.9,
            "n_typeii_detected": 5,
            "assoc_frac_fast": 0.8,
        },
        p,
    )
    txt = p.read_text()
    assert r"\newcommand{\tiiRealFracFast}{0.8}" in txt
    assert r"\newcommand{\tiiSynFracFast}{--}" in txt
