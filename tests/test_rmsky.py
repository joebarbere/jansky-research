"""Tests for jansky_research.rmsky -- the Galactic Faraday rotation sky. No network."""

from __future__ import annotations

import numpy as np
from jansky import polarization

from jansky_research import rmsky


def test_rm_from_angles_roundtrips_faraday_rotate():
    wl = np.array([0.18, 0.20, 0.21, 0.22])  # metres, no npi wrap over this span
    rm_true, chi0 = 30.0, 0.4
    angles = polarization.faraday_rotate(chi0, rm_true, wl)
    rm = rmsky.rm_from_angles(wl, angles)
    assert abs(rm - rm_true) < 0.5


def test_latitude_profile_declines_from_plane():
    sky = rmsky.synthetic_rm_sky(seed=1)
    prof = rmsky.latitude_profile(sky["rm"], sky["b"])
    med = [p["median_abs_rm"] for p in prof]
    assert med[0] > med[-1]  # plane brighter than poles
    assert all(np.isfinite(med)) and all(p["n"] > 0 for p in prof)


def test_enhancement_ratio_above_unity():
    sky = rmsky.synthetic_rm_sky(seed=2)
    r = rmsky.enhancement_ratio(sky["rm"], sky["b"])
    assert r > 3.0  # strong disk enhancement injected


def test_sign_asymmetry_recovers_antisymmetry():
    sky = rmsky.synthetic_rm_sky(seed=3)
    a = rmsky.sign_asymmetry(sky["rm"], sky["l"], sky["b"])
    assert a["inner_north"] > 0 and a["inner_south"] < 0  # north +, south - in the inner Galaxy
    assert a["inner_north_n"] > 0 and a["inner_south_n"] > 0
    assert a["inner_north_se"] > 0 and np.isfinite(a["inner_south_se"])  # standard errors reported


def test_synthetic_sky_shape_and_truth():
    sky = rmsky.synthetic_rm_sky(n_sources=2500, seed=4)
    assert sky["rm"].shape == (2500,) and sky["truth_disk_amp"] == 60.0
    assert sky["l"].min() >= 0 and sky["l"].max() <= 360
    assert -90 <= sky["b"].min() and sky["b"].max() <= 90


def test_run_offline(tmp_path):
    m = rmsky.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_sources"] > 1000
    assert m["enhancement_ratio"] > 3.0
    assert m["enhancement_ratio_se"] is not None and m["enhancement_ratio_se"] > 0
    assert m["median_abs_rm_plane"] > m["median_abs_rm_pole"]
    assert m["inner_north_rm"] > 0 and m["inner_south_rm"] < 0
    assert len(m["profile"]) == 4
    assert (tmp_path / "results" / "rmsky_metrics.json").exists()
    assert (tmp_path / "papers" / "rmsky" / "figures" / "rmsky.pdf").exists()
    macros = (tmp_path / "papers" / "rmsky" / "generated" / "macros.tex").read_text()
    assert r"\rmRatio" in macros and r"\rmInnerNorth" in macros
