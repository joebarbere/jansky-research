"""Tests for jansky_research.triangulate -- two-spacecraft type III triangulation. No network."""

from __future__ import annotations

import numpy as np

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
