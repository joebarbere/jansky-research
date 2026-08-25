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
    assert r"\sbSpeedC" in macros and r"\sbSynRatio" in macros


def test_speed_grid_from_one_ridge():
    """The grid re-maps one ridge through the model choices; it must bracket the headline point.

    The paper's grid was hand-typed from a superseded run, so its harmonic/1x value (0.137)
    contradicted \\sbSpeedC (0.1347) three lines above it. Computed from the same ridge, the
    middle grid point IS the headline number by construction.
    """
    burst = solarbursts.synthetic_burst()
    window = solarbursts.find_burst_window(burst["data"], burst["times"])
    rf_, rt_ = solarbursts.detect_burst_ridge(
        burst["data"], burst["freqs"], burst["times"], window=window
    )
    grid = solarbursts.speed_grid(rf_, rt_)
    assert [(g["harmonic"], g["fold"]) for g in grid] == [(1, 1.0), (2, 1.0), (2, 4.0)]
    fund, harm, harm4 = (g["speed_c"] for g in grid)
    spd = solarbursts.exciter_speed(rf_, rt_, harmonic=2, fold=1.0)
    assert harm == round(spd["speed_c"], 4)  # the middle point is the headline number
    # fundamental emission puts the plasma level closer in -> slower; fold 4 pushes it out -> faster
    assert fund < harm < harm4


def test_run_offline_emits_grid_macros_and_ridge(tmp_path):
    solarbursts.run(out=str(tmp_path), offline=True)
    macros = (tmp_path / "papers" / "solarbursts" / "generated" / "macros.tex").read_text()
    for name in (r"\sbGridFundOne", r"\sbGridHarmOne", r"\sbGridHarmFour"):
        assert name in macros, name
    ridge = (tmp_path / "results" / "solarbursts_ridge.csv").read_text().splitlines()
    assert ridge[0].startswith("#") and ridge[1] == "freq_mhz,time_s,used" and len(ridge) > 10


def test_robust_linfit_converges_to_a_fixed_point():
    # The referee traced the headline moving 0.111-0.147 c over the legacy hard-coded n_iter=3,
    # with the returned slope and mask from different iterations. converge=True must reach a
    # mask fixed point and return the slope fitted on exactly that mask.
    rng = np.random.default_rng(2)
    x = np.linspace(0.0, 10.0, 80)
    y = 2.0 * x + 1.0 + rng.normal(0.0, 0.3, x.size)
    y[::7] += rng.normal(0.0, 6.0, x[::7].size)  # heavy outlier tail -> multiple clip rounds
    m, b, keep = solarbursts._robust_linfit(x, y, converge=True)
    # fixed point: one more clip round changes nothing
    resid = y - (m * x + b)
    s = np.std(resid[keep])
    assert np.array_equal(np.abs(resid) < 3.0 * s, keep)
    # slope was fitted on exactly the returned mask
    m2, b2 = np.polyfit(x[keep], y[keep], 1)
    assert np.isclose(m, m2) and np.isclose(b, b2)
    assert abs(m - 2.0) < 0.1


def test_isolated_channels_flags_band_edge_singletons():
    rf = np.array([10.0, 25.5, 26.0, 26.5, 27.0, 62.4, 78.9])
    iso = solarbursts._isolated_channels(rf, gap_mhz=10.0)
    assert iso[0] and iso[-1]  # 10 MHz and 78.9 MHz are isolated
    assert not iso[1] and not iso[2]


def test_run_offline_commits_sensitivity_and_provenance(tmp_path):
    m = solarbursts.run(out=str(tmp_path), offline=True)
    assert m["pad_s"] == 10.0 and m["clip_sigma"] == 3.0 and m["fit_converged"]
    assert set(m["speed_sensitivity"]) >= {"pad_5s", "pad_10s", "clip_sigma_2.5", "clip_sigma_3.5"}
    assert m["speed_c_min"] <= m["speed_c"] <= m["speed_c_max"]
    assert m["fit_f_lo_mhz"] >= m["f_lo_mhz"] and m["fit_f_hi_mhz"] <= m["f_hi_mhz"]
    for row in m["speed_grid"]:
        assert "n_used" in row and "r2" in row
    ridge = (tmp_path / "results" / "solarbursts_ridge.csv").read_text().splitlines()
    assert ridge[0].startswith("#") and "pad_s=10" in ridge[0]
    assert ridge[1] == "freq_mhz,time_s,used"
