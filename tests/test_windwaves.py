"""Tests for jansky_research.windwaves -- interplanetary type III (Wind/WAVES). No network."""

from __future__ import annotations

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
    assert r"\wwSpeedC" in macros and r"\wwRhiAU" in macros
