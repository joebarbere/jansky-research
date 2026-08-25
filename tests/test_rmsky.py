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


def test_sign_fraction_recovers_organisation_and_bounds_the_alias():
    sky = rmsky.synthetic_rm_sky(seed=5)
    f = rmsky.sign_fraction(sky["rm"], sky["l"], sky["b"])
    assert f["inner_north"] > 0.5 > f["inner_south"]  # organised sky departs from 0.5
    # the alias-exposure bound: wrapping k sources by +-652.9 moves a region fraction by at
    # most k / n_region (each wrap can flip at most its own source's sign)
    rm2 = sky["rm"].copy()
    big = np.argsort(-np.abs(rm2))[:50]
    rm2[big] -= np.sign(rm2[big]) * rmsky.NPI_ALIAS_RAD_M2
    f2 = rmsky.sign_fraction(rm2, sky["l"], sky["b"])
    for name in ("inner_north", "inner_south", "outer_north", "outer_south"):
        assert abs(f2[name] - f[name]) <= 50.0 / f[f"{name}_n"] + 1e-12


def test_jackknife_exceeds_bootstrap_on_a_correlated_field():
    """The i.i.d. bootstrap is the correct error model on rmsky's own deterministic-signal
    fixture BY CONSTRUCTION, so no test on that fixture can expose the correlated-sky failure.
    This one can: on a Gaussian-blob field with a known coherence scale (the rmstructure
    fixture), the sky-block jackknife must exceed the i.i.d. source bootstrap."""
    from jansky_research.rmstructure import spatial_block_jackknife, synthetic_rm_screen

    s = synthetic_rm_screen(n_sources=2500, coherence_deg=3.0, plane_boost=5.0, seed=2)
    rm, gl, gb = s["rm"], s["gal_l"], s["gal_b"]
    boot = rmsky._ratio_bootstrap_se(rm, gb, plane_deg=5.0, pole_deg=12.0)
    jk = spatial_block_jackknife(
        gl,
        gb,
        lambda keep: rmsky.enhancement_ratio(rm[keep], gb[keep], plane_deg=5.0, pole_deg=12.0),
        block_deg=10.0,
    )
    assert np.isfinite(boot) and np.isfinite(jk["se"])
    assert jk["se"] > boot  # the bootstrap resamples within correlated patches


def test_run_offline(tmp_path):
    m = rmsky.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_sources"] > 1000
    assert m["enhancement_ratio"] > 3.0
    assert m["enhancement_ratio_se"] is not None and m["enhancement_ratio_se"] > 0
    # the headline error is the sky-block jackknife, committed with its block count
    assert m["enhancement_ratio_jk_se"] is not None and m["enhancement_ratio_jk_se"] > 0
    assert m["jk_n_blocks"] >= 5
    assert m["median_abs_rm_plane"] > m["median_abs_rm_pole"]
    assert m["inner_north_rm"] > 0 and m["inner_south_rm"] < 0
    # every region ships its SE, jackknife SE, count, and alias-immune sign fraction
    for name in ("inner_north", "inner_south", "outer_north", "outer_south"):
        assert m[f"{name}_se"] is not None
        assert m[f"{name}_jk_se"] is not None
        assert m[f"{name}_n"] > 0
        assert 0.0 <= m[f"{name}_frac_pos"] <= 1.0
    assert m["inner_north_frac_pos"] > 0.5 > m["inner_south_frac_pos"]
    # robustness variants committed
    assert m["enhancement_ratio_cut300"] is not None
    assert m["enhancement_ratio_alt_bins_5_70"] is not None
    assert len(m["profile"]) == 4
    assert all(p["jk_se"] is not None for p in m["profile"])
    assert (tmp_path / "results" / "rmsky_metrics.json").exists()
    assert (tmp_path / "papers" / "rmsky" / "figures" / "rmsky.pdf").exists()
    macros = (tmp_path / "papers" / "rmsky" / "generated" / "macros.tex").read_text()
    assert r"\rmRatio" in macros and r"\rmInnerNorth" in macros
    assert r"\rmRatioJkErr" in macros and r"\rmInnerNorthFrac" in macros
    assert r"\rmOuterNorthErr" in macros  # previously computed and thrown away
    assert r"\rmTruth" not in macros  # synthetic-only key no longer ships in the macro file


def test_write_catalogue_roundtrips(tmp_path):
    import csv
    import gzip

    sky = {
        "rm": np.array([10.0, -20.5]),
        "l": np.array([1.0, 2.0]),
        "b": np.array([3.0, -4.0]),
        "e_rm": np.array([1.5, np.nan]),
    }
    path = tmp_path / "cat.csv.gz"
    rmsky._write_catalogue(sky, path)
    with gzip.open(path, "rt") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["gal_l_deg", "gal_b_deg", "rm_rad_m2", "e_rm_rad_m2"]
    assert rows[1][2] == "10.0" and rows[2][3] == ""  # nan e_RM -> empty field
