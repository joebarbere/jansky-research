"""Tests for jansky_research.hi — tangent-point rotation curve. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import hi


def test_tangent_point_formula():
    r, v = hi.tangent_point(30.0, 100.0, R0=8.0, V0=200.0)
    assert np.isclose(r, 8.0 * 0.5)  # R0 sin l
    assert np.isclose(v, 100.0 + 200.0 * 0.5)  # v_term + V0 sin l


def test_terminal_velocity_finds_edge():
    vel = np.linspace(-50.0, 300.0, 400)
    spec = np.where(vel < 150.0, 20.0, 0.0)
    spec[vel < 0] = 0.0
    vt = hi.terminal_velocity(vel, spec, threshold_k=2.0)
    assert abs(vt - 150.0) < 2.0  # within a channel of the edge
    assert np.isnan(hi.terminal_velocity(vel, np.zeros_like(vel)))  # no emission -> nan


def test_synthetic_rotation_curve_recovers_flat():
    # With a realistic (5 km/s) edge width the fixed-threshold crossing sits on the wing, a
    # measured ~w*ln(A/T_thr) ~ 13 km/s above the injected edge -- the fixture now exposes the
    # dominant real-data systematic instead of hiding it behind a 0.6 km/s edge.
    longs = np.array([20.0, 40.0, 60.0, 80.0])
    slices = [hi.synthetic_lv_slice(ell, v_flat=230.0, seed=i) for i, ell in enumerate(longs)]
    R, V = hi.rotation_curve(longs, slices)
    assert np.all(np.diff(R) > 0)  # sorted by radius
    bias = V - 230.0
    assert np.all(bias > 5.0) and np.all(bias < 20.0)  # wing overshoot, roughly w*ln(A/thr)
    # the edge estimator recovers the injected curve without the overshoot
    Re, Ve = hi.rotation_curve(longs, slices, estimator="edge")
    assert np.all(np.abs(Ve - 230.0) < 6.0)


def test_terminal_velocity_edge_recovers_half_height():
    vel = np.linspace(-50.0, 300.0, 500)
    spec = 30.0 / (1.0 + np.exp((vel - 150.0) / 5.0))
    v0, w = hi.terminal_velocity_edge(vel, spec)
    assert abs(v0 - 150.0) < 2.0  # half-height point, not the wing crossing
    assert 2.0 < w < 12.0  # erf-equivalent width of a logistic-5 edge is ~8
    vt = hi.terminal_velocity(vel, spec)
    assert vt - v0 > 8.0  # the threshold crossing overshoots on the wing
    assert np.isnan(hi.terminal_velocity_edge(vel, np.zeros_like(vel))[0])  # no emission


def test_run_offline(tmp_path):
    m = hi.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    # threshold estimator carries the measured wing overshoot; the edge estimator does not
    assert abs(m["V_flat_edge_mean_kms"] - 230.0) < 6.0
    assert 5.0 < m["threshold_minus_edge_kms"] < 20.0
    assert m["V_flat_scatter_kms"] < 8.0  # flat
    assert m["V_flat_sem_kms"] < m["V_flat_scatter_kms"]
    # dense by default: 10-80 deg every 1 deg is 71 sightlines, 51 of them beyond the bar cut
    assert m["step_deg"] == 1.0
    assert len(m["longitudes_deg"]) == 71
    assert m["n_flat"] == 51
    # The fixture injects a genuinely flat curve, so here the slope *must* come out consistent
    # with zero -- this is the control for the real run, where it does not.
    assert m["slope_sigma"] < 3.0
    assert m["slope_edge_sigma"] < 3.0
    assert m["keplerian_at_rmax_kms"] < m["V_flat_mean_kms"]  # the non-Keplerian contrast
    assert len(m["threshold_sweep_vflat"]) == 7
    assert (tmp_path / "results" / "rotation_curve.json").exists()
    assert (tmp_path / "papers" / "hi" / "figures" / "rotation_curve.pdf").exists()
    macros = (tmp_path / "papers" / "hi" / "generated" / "macros.tex").read_text()
    # namespaced: an offline run fills the synthetic side only
    assert r"\newcommand{\hiSynVflat}" in macros
    assert r"\newcommand{\hiRealVflat}{--}" in macros
    assert r"\hiSynThrMinusEdge" in macros and r"\hiSynKepler" in macros
    assert r"\hiSource" in macros
    # A non-finite metric must never reach the paper as the literal string "nan".
    assert "nan" not in macros
    # The reference comparison runs offline too, against the fixture's injected curve.
    assert r"\hiSynRefEdgeOffset" in macros and r"\newcommand{\hiRealRefEdgeOffset}{--}" in macros


def test_read_lab_slice_parses_wcs(tmp_path):
    from astropy.io import fits

    nb, nv = 5, 11
    data = np.arange(nb * nv, dtype=float).reshape(1, nb, nv)
    h = fits.PrimaryHDU(data).header
    h["CTYPE1"], h["CRVAL1"], h["CDELT1"], h["CRPIX1"] = (
        "VELO-LSR",
        0.0,
        20000.0,
        1,
    )  # m/s, 20 km/s/ch
    h["CUNIT1"] = "M/S"
    h["CTYPE2"], h["CRVAL2"], h["CDELT2"], h["CRPIX2"] = "GLAT-CAR", -2.0, 1.0, 1
    p = tmp_path / "slice.fits"
    fits.PrimaryHDU(data, h).writeto(p)
    lat, vel, d = hi.read_lab_slice(p)
    assert d.shape == (nb, nv)
    assert np.isclose(vel[0], 0.0) and np.isclose(vel[1], 20.0)  # 20000 m/s -> 20 km/s
    assert np.isclose(lat[0], -2.0) and np.isclose(lat[-1], 2.0)


def test_run_offline_recovers_the_injected_reference_curve(tmp_path):
    """The offline leg exercises the same comparison path a real run uses, with a known answer.

    On a real run the reference is the VGPS curve of McClure-Griffiths & Dickey (2016); here it
    is the curve the fixture injected, so the edge estimator must reproduce it and the threshold
    estimator must sit above it by the wing overshoot. Without this the comparison code would
    only ever be exercised by the network leg.
    """
    m = hi.run(out=str(tmp_path), offline=True)
    edge, thr = m["compare_edge"], m["compare_threshold"]
    assert edge["n"] > 50 and thr["n"] == edge["n"]
    assert abs(edge["mean"]) < 1.0  # the edge fit recovers the injected terminal velocities
    assert 5.0 < thr["mean"] < 25.0  # the 2 K crossing sits out on the wing
    # and the gap closes as the threshold climbs the edge -- the erfc model's prediction
    sweep = m["threshold_sweep_minus_reference"]
    assert sweep["20"] < sweep["10"] < sweep["5"] < sweep["2"]


def test_offline_fixture_exposes_the_width_bias_relation(tmp_path):
    """The threshold overshoot scales with the width of the profile edge.

    The fixture varies its edge width across longitudes precisely so this relation is
    measurable offline; a single-width fixture reports a meaningless correlation and cannot
    catch a regression in the estimator pair.
    """
    m = hi.run(out=str(tmp_path), offline=True)
    assert m["width_bias_corr"] > 0.8
    assert 1.0 < m["width_bias_slope"] < 3.0


def test_compare_terminal_velocities_matches_and_summarises():
    ref_l = np.arange(20.0, 40.0, 0.065)
    ref_v = 100.0 - 2.0 * (ref_l - 20.0)
    lons = np.array([10.0, 25.0, 30.0, 35.0, 60.0])  # 10 and 60 fall outside the reference
    v_term = np.interp(lons, ref_l, ref_v) + 3.0  # a uniform +3 km/s offset
    c = hi.compare_terminal_velocities(lons, v_term, ref_l, ref_v)
    assert c["n"] == 3  # only the covered longitudes are compared, never extrapolated
    assert c["longitudes_deg"] == [25.0, 30.0, 35.0]
    assert abs(c["mean"] - 3.0) < 0.1
    assert c["sd"] < 0.1


def test_compare_terminal_velocities_drops_nan_and_thin_coverage():
    ref_l = np.array([30.0, 30.1, 30.2])  # only 3 samples: below min_points
    ref_v = np.array([100.0, 100.0, 100.0])
    c = hi.compare_terminal_velocities(np.array([30.0]), np.array([103.0]), ref_l, ref_v)
    assert c["n"] == 0 and np.isnan(c["mean"])
    c2 = hi.compare_terminal_velocities(
        np.array([30.0]), np.array([103.0]), ref_l, ref_v, min_points=3
    )
    assert c2["n"] == 1 and np.isclose(c2["mean"], 3.0)
    # a failed estimator (nan) is dropped rather than poisoning the mean
    c3 = hi.compare_terminal_velocities(
        np.array([30.0]), np.array([np.nan]), ref_l, ref_v, min_points=3
    )
    assert c3["n"] == 0


def test_read_mgd2016_parses_vizier_tsv(tmp_path):
    p = tmp_path / "mgd.tsv"
    p.write_text(
        "#comment\n\nGLON\tvLSR\ndeg\tkm/s\n------\t------\n"
        "18.412\t134.17\n18.478\t134.17\n66.968\t 18.73\n"
    )
    ell, vel = hi.read_mgd2016(p)
    assert ell.size == 3 and np.isclose(ell[0], 18.412) and np.isclose(vel[-1], 18.73)
    import pytest

    bad = tmp_path / "empty.tsv"
    bad.write_text("#only a comment\n")
    with pytest.raises(ValueError, match="no data rows"):
        hi.read_mgd2016(bad)


def test_contiguity_counts_runs():
    vel = np.linspace(0.0, 100.0, 101)
    spec = np.zeros_like(vel)
    spec[10:20] = 5.0
    assert hi._contiguity(vel, spec, 2.0) == 1
    spec[40:45] = 5.0  # a detached second run, e.g. a local-emission island
    assert hi._contiguity(vel, spec, 2.0) == 2
    assert hi._contiguity(vel, np.zeros_like(vel), 2.0) == 0
