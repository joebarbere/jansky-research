"""Tests for jansky_research.driftsearch — injection-recovery benchmark. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import driftsearch


def test_completeness_snr_interpolates():
    snrs = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    p = np.array([0.0, 0.2, 0.5, 0.8, 1.0])
    assert abs(driftsearch.completeness_snr(snrs, p, 0.5) - 2.0) < 1e-9
    # never reaches the level -> nan
    assert np.isnan(driftsearch.completeness_snr(snrs, np.full(5, 0.3), 0.5))


def test_false_positive_rate_low_at_threshold():
    # Noise-only best S/N peaks near ~5; threshold 10 should almost never fire.
    fpr = driftsearch.false_positive_rate(n_trials=150, threshold=10.0, seed=0)
    assert fpr < 0.05


def test_injection_recovery_separates_weak_from_strong():
    res = driftsearch.injection_recovery(
        np.array([0.25, 3.0]),
        drift_rates=np.array([0.3]),
        n_trials=15,
        fpr_trials=100,
        seed=0,
    )
    pm = res.p_detect.mean(axis=1)
    assert pm[0] < 0.5  # a very weak tone is rarely recovered
    assert pm[1] > 0.9  # a strong tone is almost always recovered
    assert res.false_positive_rate < 0.05
    assert res.inj_snrs[0] <= res.completeness_snr_50 <= res.inj_snrs[-1]


def test_run_offline_writes_artifacts(tmp_path):
    m = driftsearch.run(out=str(tmp_path), n_trials=10)
    assert 0.0 <= m["false_positive_rate"] < 0.1
    assert np.isfinite(m["completeness_snr_50"])  # the grid spans the transition
    assert (tmp_path / "results" / "drift_metrics.json").exists()
    assert (tmp_path / "papers" / "driftsearch" / "figures" / "drift_recovery.pdf").exists()
    assert (tmp_path / "papers" / "driftsearch" / "generated" / "macros.tex").exists()


def test_locate_carrier_ignores_the_dc_spike():
    """The regression test the DC-spike lesson never had: a huge time-constant spike at n//2
    must not be reported as the carrier; the drifting tone elsewhere must be."""
    rng = np.random.default_rng(0)
    n_time, n_freq = 64, 4096
    wf = rng.normal(10.0, 1.0, (n_time, n_freq))
    wf[:, n_freq // 2] += 1.0e4  # the band-centre DC artifact, time-independent
    tone = 300
    for t in range(n_time):
        wf[t, tone + int(round(0.2 * t))] += 30.0  # a genuine drifting tone
    loc = driftsearch.locate_carrier(wf, dc_halfwidth=16)
    assert loc["dc_channel"] == n_freq // 2 and loc["dc_is_band_centre"]
    assert abs(loc["channel"] - (tone + 6)) < 20  # the tone, not the spike
    # and the search at the located channel recovers it while the DC window screams artifact
    d = driftsearch.measure_drift(wf, loc["channel"], halfwidth=40)
    assert abs(d["chan_per_sample"] - 0.2) < 0.05


def test_measure_drift_recovers_slope():
    rng = np.random.default_rng(1)
    wf = rng.normal(0.0, 1.0, (32, 512))
    for t in range(32):
        wf[t, 200 + int(round(1.5 * t))] += 25.0
    d = driftsearch.measure_drift(wf, 220, halfwidth=80)
    assert abs(d["chan_per_sample"] - 1.5) < 0.1
    assert d["fit_rms_chan"] < 1.0


def test_noise_peak_stats_far_from_threshold():
    s = driftsearch.noise_peak_stats(n_draws=60)
    # the noise-only best S/N clusters near 4, an order of magnitude below threshold 10
    assert 3.5 < s["mean"] < 4.8
    assert s["max"] < 8.0
    assert 0.03 < s["fpr_upper_95_one_sided"] < 0.06  # 1 - 0.05^(1/60)


def test_run_commits_config_matrix_and_off_grid(tmp_path):
    m = driftsearch.run(str(tmp_path), n_trials=10)
    assert m["n_trials"] == 10
    assert m["config"]["n_search_drifts"] == 41 and m["config"]["fpr_trials"] == 400
    # the FULL matrix is committed, per-drift crossings included
    assert len(m["p_detect"]) == len(m["inj_snrs"])
    assert len(m["p_detect"][0]) == len(m["drift_rates"])
    assert len(m["completeness_snr_50_per_drift"]) == len(m["drift_rates"])
    # the off-grid check exists and lands near the on-grid value
    assert abs(m["off_grid_check"]["completeness_snr_50"] - m["completeness_snr_50"]) < 0.3
    macros = (tmp_path / "papers" / "driftsearch" / "generated" / "macros.tex").read_text()
    assert r"\dsNoiseMean" in macros and r"\dsFprBound" in macros
    # the Voyager namespace exists as placeholders and is never fabricated by the benchmark
    assert r"\newcommand{\dsVoySnr}{--}" in macros
    assert m["source"].startswith("synthetic")


def test_real_voyager_recovery_when_file_cached():
    """The assertion that would have caught the 0.92 MHz targeting error: gated on the cache."""
    import pathlib

    import pytest

    path = pathlib.Path("data/Voyager1.single_coarse.fine_res.h5")
    if not path.exists():
        pytest.skip("Voyager file not cached")
    pytest.importorskip("h5py")
    pytest.importorskip("hdf5plugin")
    m = driftsearch.validate_voyager(str(path))
    assert m["recovered"]
    assert m["carrier"]["snr"] > 100  # three orders of magnitude above the noise floor
    assert abs(m["carrier"]["measured_drift_hz_s"] + 0.374) < 0.05
    assert m["dc_spike"]["is_band_centre"]
    # the legacy asserted frequency remains blank sky -- the recorded lesson
    assert m["legacy_asserted"]["snr"] < m["carrier"]["snr"] / 50
