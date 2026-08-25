"""Tests for jansky_research.triangulate -- two-spacecraft type III triangulation. No network."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import triangulate


def test_direction_unit_known_axes():
    # colatitude 90, azimuth 0 -> +X; azimuth 90 -> +Y; colatitude 0 -> +Z
    assert np.allclose(triangulate.direction_unit(0.0, 90.0), [1, 0, 0], atol=1e-9)
    assert np.allclose(triangulate.direction_unit(90.0, 90.0), [0, 1, 0], atol=1e-9)
    assert np.allclose(triangulate.direction_unit(0.0, 0.0), [0, 0, 1], atol=1e-9)


def test_mean_direction_weighted_and_nan_safe():
    az = np.array([10.0, 10.0, 10.0, np.nan, 50.0])
    col = np.array([90.0, 90.0, 90.0, 90.0, 90.0])
    w = np.array([1.0, 3.0, 2.0, 1.0, np.nan])
    u, n = triangulate.mean_direction(az, col, w)
    assert n == 3  # the NaN-direction and NaN-weight samples are dropped, leaving three
    assert u is not None and np.isclose(np.linalg.norm(u), 1.0)
    # too few samples -> None
    u2, n2 = triangulate.mean_direction(az[:1], col[:1], w[:1])
    assert u2 is None and n2 == 0


def test_triangulate_rays_exact_intersection():
    # two rays that meet exactly at (10, 0, 0): zero miss, both t>0
    p1 = np.array([0.0, -5.0, 0.0])
    u1 = np.array([10.0, 5.0, 0.0])
    u1 = u1 / np.linalg.norm(u1)
    p2 = np.array([0.0, 5.0, 0.0])
    u2 = np.array([10.0, -5.0, 0.0])
    u2 = u2 / np.linalg.norm(u2)
    tri = triangulate.triangulate_rays(p1, u1, p2, u2)
    assert np.allclose(tri["source"], [10, 0, 0], atol=1e-6)
    assert tri["miss"] < 1e-6
    assert tri["t1"] > 0 and tri["t2"] > 0


def test_triangulate_rays_parallel_is_nan():
    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([0.0, 1.0, 0.0])
    u = np.array([1.0, 0.0, 0.0])
    tri = triangulate.triangulate_rays(p1, u, p2, u)
    assert not np.isfinite(tri["source"]).any()


def test_synthetic_event_schema_and_truth():
    ev = triangulate.synthetic_event(seed=1)
    for spec in (ev["spec_a"], ev["spec_b"]):
        assert spec["az"].shape == spec["col"].shape == spec["sfu"].shape
        assert spec["pos"].shape[1] == 3
        assert spec["freqs"][0] > spec["freqs"][-1]  # descending
    assert ev["truth"]["sep_deg"] > 90  # wide STEREO-like baseline


def test_run_offline_recovers_longitude_and_correlates(tmp_path):
    m = triangulate.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_tri"] >= 8
    # the injected source longitude/latitude are recovered within the DF noise
    assert m["lon_err_deg"] < 10.0
    assert m["lat_err_deg"] < 10.0
    # the geometric distance tracks the independent plasma-frequency distance
    assert m["corr_geom_plasma"] > 0.8
    # outputs land where the paper expects them
    assert (tmp_path / "results" / "triangulate_metrics.json").exists()
    assert (tmp_path / "papers" / "triangulate" / "figures" / "triangulate.pdf").exists()
    macros = (tmp_path / "papers" / "triangulate" / "generated" / "macros.tex").read_text()
    assert r"\triCorr" in macros and r"\triLon" in macros


def test_triangulate_track_drops_backward_and_far_misses():
    # a clean event, but force a tiny max_miss so noisy channels are rejected
    ev = triangulate.synthetic_event(seed=2, noise_deg=0.0)
    track = triangulate.triangulate_track(ev["spec_a"], ev["spec_b"], max_miss_rsun=1.0)
    # with zero noise every triangulated channel intersects almost exactly
    assert track["freq_mhz"].size > 0
    assert np.all(track["miss"] < 1.0)


def test_miss_sweep_is_pure_filtering_and_matches_the_analysis_cut():
    """The 60-row of the sweep must equal the headline analysis on the same open track."""
    ev = triangulate.synthetic_event()
    open_track = triangulate.triangulate_track(
        ev["spec_a"], ev["spec_b"], max_miss_rsun=float("inf")
    )
    sweep = triangulate.miss_sweep(open_track)
    assert [s["max_miss_rsun"] for s in sweep] == [15.0, 30.0, 60.0, 100.0]
    # n is monotone in the threshold: loosening a cut can only admit channels
    ns = [s["n"] for s in sweep]
    assert ns == sorted(ns)
    cut60 = triangulate.triangulate_track(ev["spec_a"], ev["spec_b"], max_miss_rsun=60.0)
    row60 = next(s for s in sweep if s["max_miss_rsun"] == 60.0)
    assert row60["n"] == int(cut60["freq_mhz"].size)
    if row60["n"] >= 3:
        rg, rp = cut60["r_geom"], cut60["r_plasma"]
        assert row60["corr_geom_plasma"] == pytest.approx(
            float(np.corrcoef(rg, rp)[0, 1]), abs=1e-3
        )


def test_run_offline_commits_channels_and_sweep(tmp_path):
    triangulate.run(out=str(tmp_path), offline=True)
    rows = (tmp_path / "results" / "triangulate_channels.csv").read_text().splitlines()
    assert rows[0].startswith("freq_mhz,r_geom_rsun,r_plasma_rsun,miss_rsun,lon_deg,lat_deg")
    assert rows[0].endswith("ua_x,ua_y,ua_z,ub_x,ub_y,ub_z,n_a,n_b")  # auditability columns
    assert len(rows) > 3
    macros = (tmp_path / "papers" / "triangulate" / "generated" / "macros.tex").read_text()
    for name in (r"\triSweepFifteenCorr", r"\triSweepSixtyRatio", r"\triSweepHundredN"):
        assert name in macros, name


def test_circular_median_survives_the_branch_cut():
    """The published longitude was a scalar median across +-180: 16 negative values near
    -175 and 22 positive near +175 medianed to the edge of one cluster. The circular median
    must land near 180 for such a sample, and agree with np.median away from the cut."""
    rng = np.random.default_rng(0)
    a = np.concatenate([rng.normal(176.0, 3.0, 22), rng.normal(-176.0, 3.0, 16)])
    a = ((a + 180.0) % 360.0) - 180.0
    med, scatter = triangulate.circular_median_deg(a)
    assert abs(abs(med) - 180.0) < 4.0  # near the cut, not at the 8th percentile
    assert scatter < 10.0
    # away from the wrap it reduces to the ordinary median
    b = rng.normal(35.0, 5.0, 50)
    med_b, _ = triangulate.circular_median_deg(b)
    assert med_b == pytest.approx(float(np.median(b)), abs=1e-9)


def test_run_offline_recovers_a_far_side_longitude():
    """The old fixture longitude (35 deg) was 145 deg from the branch cut, so the wrap bug
    was invisible to every test. This fixture sits ON the cut at the real event's geometry."""
    ev = triangulate.synthetic_event(lon_deg=179.0, sep_deg=82.0, seed=8)
    track = triangulate.triangulate_track(ev["spec_a"], ev["spec_b"])
    m = triangulate._metrics(track, "synthetic", 2, ev["truth"])
    assert m["lon_err_deg"] < 10.0  # fails under the scalar median


def test_constant_df_bias_produces_a_constant_additive_offset():
    """A few-degree constant azimuth bias on one spacecraft — the systematic the noise
    budget cannot represent — must show up as an ADDITIVE radial offset with OLS slope
    near 1, not as a multiplicative density-like scaling."""
    clean = triangulate.synthetic_event(lon_deg=179.0, sep_deg=82.0, noise_deg=3.0, seed=9)
    biased = triangulate.synthetic_event(
        lon_deg=179.0, sep_deg=82.0, noise_deg=3.0, bias_deg_a=3.0, seed=9
    )
    t_clean = triangulate.triangulate_track(
        clean["spec_a"], clean["spec_b"], max_miss_rsun=float("inf")
    )
    t_bias = triangulate.triangulate_track(
        biased["spec_a"], biased["spec_b"], max_miss_rsun=float("inf")
    )
    s_clean = triangulate.additive_vs_multiplicative(t_clean["r_geom"], t_clean["r_plasma"])
    s_bias = triangulate.additive_vs_multiplicative(t_bias["r_geom"], t_bias["r_plasma"])
    assert abs(s_bias["diff_med_rsun"]) > abs(s_clean["diff_med_rsun"]) + 3.0
    assert 0.7 < s_bias["ols_slope"] < 1.4  # additive, not multiplicative
    assert s_bias["rms_additive_rsun"] < s_bias["rms_multiplicative_rsun"]


def test_spacecraft_time_offset_degrades_the_track():
    """A per-file epoch-origin mismatch shifts B's burst window off A's — the failure mode
    the absolute time base in fetch_stereo_df now prevents. The synthetic knob demonstrates
    the sensitivity the old per-file origin was exposed to."""
    aligned = triangulate.synthetic_event(seed=10)
    shifted = triangulate.synthetic_event(t0_offset_b_s=2400.0, seed=10)
    t_ok = triangulate.triangulate_track(aligned["spec_a"], aligned["spec_b"])
    t_bad = triangulate.triangulate_track(shifted["spec_a"], shifted["spec_b"])
    assert t_bad["freq_mhz"].size < t_ok["freq_mhz"].size  # window misses most of B's burst


def test_harmonic_density_grid_shows_the_degeneracy():
    ev = triangulate.synthetic_event(seed=11)
    track = triangulate.triangulate_track(ev["spec_a"], ev["spec_b"])
    grid = triangulate.harmonic_density_grid(track["freq_mhz"], track["r_geom"])
    # f_p^2 prop n_e: harmonic=1 at 4x density is exactly harmonic=2 at 1x
    assert grid["h1_s4"]["ratio_med"] == pytest.approx(grid["h2_s1"]["ratio_med"], rel=1e-6)
    # fundamental emission means a higher plasma frequency, hence a smaller Leblanc
    # radius, hence a LARGER geometric/plasma ratio (the referee's 2.18 -> 3.89 direction)
    assert grid["h1_s1"]["ratio_med"] > grid["h2_s1"]["ratio_med"]


def test_synthetic_run_does_not_clobber_real_artifacts(tmp_path):
    """A bare offline run must not overwrite a real channels CSV/figure while the JSON
    marker keeps saying STEREO — the marker-that-lies failure."""
    import json

    res = tmp_path / "results"
    res.mkdir(parents=True)
    (res / "triangulate_metrics.json").write_text(
        json.dumps({"source": "STEREO-A+B L3 DF 20130515", "n_tri": 38})
    )
    (res / "triangulate_channels.csv").write_text("real,evidence\n")
    triangulate.run(out=str(tmp_path), offline=True)
    assert (res / "triangulate_channels.csv").read_text() == "real,evidence\n"  # untouched
    # and with no real file on disk, the synthetic artifacts ARE written (fresh checkout)
    triangulate.run(out=str(tmp_path / "fresh"), offline=True)
    assert (tmp_path / "fresh" / "results" / "triangulate_channels.csv").exists()
