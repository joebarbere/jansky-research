"""Tests for jansky_research.solarbursts -- type III drift -> exciter speed. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import solarbursts


def test_synthetic_burst_shape_and_truth():
    b = solarbursts.synthetic_burst(speed_c=0.25, seed=0)
    assert b["data"].shape == (200, 400)
    assert b["truth_speed_c"] == 0.25
    # frequencies descending, all within the band
    assert b["freqs"][0] > b["freqs"][-1]
    assert 20.0 <= b["freqs"].min() and b["freqs"].max() <= 90.0


def test_background_subtract_zeroes_baseline():
    data = np.tile(np.arange(50.0)[:, None], (1, 30)) + 5.0  # constant per channel
    clean = solarbursts.background_subtract(data)
    assert np.allclose(np.median(clean, axis=1), 0.0)


def test_find_burst_window_localizes():
    rng = np.random.default_rng(0)
    data = rng.normal(0.0, 1.0, (50, 400))
    times = np.linspace(0.0, 100.0, 400)
    data[:, 200:210] += 20.0  # a burst around t = 50 s
    mask = solarbursts.find_burst_window(data, times, pad_s=5.0)
    assert mask.dtype == bool and mask.sum() > 0
    assert 45.0 < times[mask].mean() < 55.0  # window centred on the burst


def test_detect_ridge_drifts_high_to_low():
    b = solarbursts.synthetic_burst(seed=1)
    rf, rt = solarbursts.detect_burst_ridge(b["data"], b["freqs"], b["times"])
    assert rf.size > 50
    # a type III drifts from high to low frequency as time increases
    slope, _ = np.polyfit(rt, rf, 1)
    assert slope < 0


def test_fit_drift_rate_negative():
    b = solarbursts.synthetic_burst(seed=2)
    rf, rt = solarbursts.detect_burst_ridge(b["data"], b["freqs"], b["times"])
    drift = solarbursts.fit_drift_rate(rf, rt)
    assert drift < 0  # MHz/s, frequency falling
    assert np.isfinite(drift)


def test_exciter_speed_recovers_injected_speed():
    """The forward fixture and the inverse share the Newkirk mapping, so a clean burst round-trips."""
    for v in (0.2, 0.3, 0.4):
        b = solarbursts.synthetic_burst(speed_c=v, harmonic=2, seed=3)
        rf, rt = solarbursts.detect_burst_ridge(b["data"], b["freqs"], b["times"])
        spd = solarbursts.exciter_speed(rf, rt, harmonic=2)
        assert abs(spd["speed_c"] - v) / v < 0.1  # within 10%
        assert 1.0 < spd["r_lo"] < spd["r_hi"] < 5.0  # plausible coronal heights


def test_harmonic_assumption_changes_the_speed():
    """Analysing a 2f burst as fundamental (f=fp) changes the inferred density/height -> speed."""
    b = solarbursts.synthetic_burst(speed_c=0.3, harmonic=2, seed=4)
    rf, rt = solarbursts.detect_burst_ridge(b["data"], b["freqs"], b["times"])
    matched = solarbursts.exciter_speed(rf, rt, harmonic=2)["speed_c"]
    mismatched = solarbursts.exciter_speed(rf, rt, harmonic=1)["speed_c"]
    assert abs(matched - 0.3) / 0.3 < 0.1
    assert abs(mismatched - matched) > 0.02  # the factor-2 systematic is real


def test_robust_linfit_rejects_outliers():
    x = np.linspace(0.0, 10.0, 50)
    y = 2.0 * x + 1.0
    y[5] += 50.0  # two gross outliers (e.g. RFI-corrupted ridge points)
    y[30] -= 40.0
    m, b, keep = solarbursts._robust_linfit(x, y)
    assert abs(m - 2.0) < 0.05 and abs(b - 1.0) < 0.2  # slope/intercept recovered
    assert not keep[5] and not keep[30]  # the outliers are rejected


def test_exciter_speed_reports_quality():
    b = solarbursts.synthetic_burst(speed_c=0.3, seed=7)
    rf, rt = solarbursts.detect_burst_ridge(b["data"], b["freqs"], b["times"])
    spd = solarbursts.exciter_speed(rf, rt, harmonic=2)
    assert spd["r2"] > 0.9  # a clean single burst is a tight straight height-time track
    assert 0 < spd["n_used"] <= spd["n_points"]


def test_run_offline(tmp_path):
    m = solarbursts.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_ridge"] > 50
    assert m["drift_mhz_per_s"] < 0
    assert 0.1 < m["speed_c"] < 0.5  # canonical type III exciter speed
    assert 0.85 < m["recovery_ratio"] < 1.15
    assert m["r2"] > 0.9 and m["n_used"] > 50  # clean, coherent ridge
    assert (tmp_path / "results" / "solarbursts_metrics.json").exists()
    assert (tmp_path / "papers" / "solarbursts" / "figures" / "burst.pdf").exists()
    macros = (tmp_path / "papers" / "solarbursts" / "generated" / "macros.tex").read_text()
    assert r"\sbSpeedC" in macros and r"\sbRatio" in macros
