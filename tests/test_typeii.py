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


def test_parse_lwa_dspec_reads_fits(tmp_path):
    from astropy.io import fits

    # a wiki-format dspec FITS: 2D spectrum + 1D freq (Hz, ascending) + 1D time HDUs
    n_f, n_t = 30, 50
    freqs_hz = np.linspace(13.4e6, 86.9e6, n_f)  # ascending Hz
    times = np.arange(n_t) * 0.256
    spec = np.random.default_rng(0).normal(0, 1, (n_f, n_t))
    spec[15, :] += 20.0  # a bright channel to check orientation survives
    p = tmp_path / "ovro-lwa.lev1_bmf_256ms_96kHz.2024-05-14.dspec_I.fits"
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(spec, name="DSPEC"),
            fits.ImageHDU(freqs_hz, name="FREQ"),
            fits.ImageHDU(times, name="TIME"),
        ]
    ).writeto(p)
    d = t2.parse_lwa_dspec(p)
    assert d["data"].shape == (n_f, n_t)
    assert d["freqs"][0] > d["freqs"][-1]  # descending MHz (as solarbursts expects)
    assert 13.0 < d["freqs"].min() < 14.0 and 86.0 < d["freqs"].max() < 87.0  # Hz->MHz
    assert d["times"].size == n_t


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
