"""Tests for jansky_research.windwaves -- interplanetary type III (Wind/WAVES). No network."""

from __future__ import annotations

import numpy as np

from jansky_research import windwaves


def test_leblanc_density_at_1au():
    # the Leblanc model is normalised to ~7.2 cm^-3 at 1 AU (215 R_sun)
    assert abs(windwaves.leblanc_density(215.0) - 7.2) < 0.5


def test_leblanc_radius_roundtrips():
    for r in (2.0, 10.0, 50.0, 200.0):
        ne = windwaves.leblanc_density(r)
        assert abs(windwaves.leblanc_radius(ne) - r) / r < 0.02


def test_emission_radius_monotonic():
    # higher frequency -> closer to the Sun (denser plasma)
    r_hi = windwaves.emission_radius(10.0, harmonic=2)
    r_lo = windwaves.emission_radius(0.5, harmonic=2)
    assert 1.5 < r_hi < r_lo  # 10 MHz nearer the Sun than 0.5 MHz


def test_beam_speed_recovers_injected():
    """Forward (Leblanc) fixture and inverse share the model, so a clean burst round-trips."""
    for v in (0.1, 0.15, 0.2):
        b = windwaves.synthetic_ip_burst(speed_c=v, harmonic=2, seed=1)
        rf, rt = windwaves.solarbursts.detect_burst_ridge(b["data"], b["freqs"], b["times"])
        spd = windwaves.beam_speed(rf, rt, harmonic=2)
        assert abs(spd["speed_c"] - v) / v < 0.15
        assert spd["r2"] > 0.9 and spd["r_hi"] > spd["r_lo"] > 1.0


def test_synthetic_burst_shape():
    b = windwaves.synthetic_ip_burst(seed=0)
    assert b["data"].shape == (256, 600)
    assert b["freqs"][0] > b["freqs"][-1]  # descending
    assert b["truth_speed_c"] == 0.15


def test_run_offline(tmp_path):
    m = windwaves.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_ridge"] > 30
    assert 0.05 < m["speed_c"] < 0.4  # interplanetary beam speed
    assert 0.85 < m["recovery_ratio"] < 1.15
    assert m["r_hi_rsun"] > m["r_lo_rsun"] > 1.0
    assert (tmp_path / "results" / "windwaves_metrics.json").exists()
    assert (tmp_path / "papers" / "windwaves" / "figures" / "ipburst.pdf").exists()
    macros = (tmp_path / "papers" / "windwaves" / "generated" / "macros.tex").read_text()
    assert r"\wwSynSpeedC" in macros and r"\wwSynRhiAU" in macros
    assert r"\newcommand{\wwRealSpeedC}{--}" in macros


def test_emission_radius_density_scale_degeneracy():
    # f_p^2 ~ n: harmonic=1 at density_scale=4 must equal harmonic=2 at scale=1 exactly —
    # the emission-mode and density-enhancement systematics are one axis, not two.
    f = np.array([0.5, 1.0, 5.0, 13.8])
    a = windwaves.emission_radius(f, harmonic=1, density_scale=4.0)
    b = windwaves.emission_radius(f, harmonic=2, density_scale=1.0)
    assert np.allclose(a, b, rtol=1e-9)


def test_beam_speed_reports_jackknife_and_estimator_bracket():
    b = windwaves.synthetic_ip_burst(speed_c=0.12, n_time=41, duration_s=2400.0, seed=2)
    from jansky_research import solarbursts

    w = solarbursts.find_burst_window(b["data"], b["times"], pad_s=2400.0)
    rf, rt = solarbursts.detect_burst_ridge(b["data"], b["freqs"], b["times"], window=w)
    spd = windwaves.beam_speed(rf, rt)
    assert np.isfinite(spd["speed_c_se"]) and spd["speed_c_se"] > 0
    assert spd["n_time_cols"] >= 4
    # estimator family brackets the truth on a clean constant-speed burst
    assert 0.09 < spd["speed_c"] < 0.15
    assert np.isfinite(spd["speed_c_inverse"]) and np.isfinite(spd["speed_c_points"])
    grid = windwaves.speed_grid(rf, rt)
    assert len(grid) == len(windwaves.SPEED_GRID)
    # the degenerate pair: fundamental x4 density is not in the grid, but harmonic rows must
    # order by density scale (higher scale -> larger radii -> faster)
    h2 = {g["density_scale"]: g["speed_c"] for g in grid if g["harmonic"] == 2}
    assert h2[1.0] < h2[2.0] < h2[4.0]


def test_matched_cadence_fixture_measures_quantisation_bias():
    # At the real one-minute cadence (n_time=41 over 2400 s) the recovered slope is biased low
    # by a few percent relative to the fine-cadence fixture — a bias the old 3 s fixture could
    # not see (round-6 referee). The estimator must stay within 10% of truth even so.
    from jansky_research import solarbursts

    ratios = []
    for n_time, dur in ((41, 2400.0), (600, 1800.0)):
        b = windwaves.synthetic_ip_burst(speed_c=0.15, n_time=n_time, duration_s=dur, seed=3)
        w = solarbursts.find_burst_window(b["data"], b["times"], pad_s=dur)
        rf, rt = solarbursts.detect_burst_ridge(b["data"], b["freqs"], b["times"], window=w)
        spd = windwaves.beam_speed(rf, rt)
        ratios.append(spd["speed_c"] / 0.15)
    assert 0.90 < ratios[0] < 1.05  # coarse cadence: small negative bias tolerated
    assert 0.95 < ratios[1] < 1.05  # fine cadence: unbiased
