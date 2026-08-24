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
    assert m["n_flat"] >= 4
    # flatness is quantified: slope consistent with zero within ~2 sigma on the fixture
    assert abs(m["slope_kms_per_kpc"]) < 3.0 * m["slope_se_kms_per_kpc"] + 3.0
    assert m["keplerian_at_rmax_kms"] < m["V_flat_mean_kms"]  # the non-Keplerian contrast
    assert len(m["threshold_sweep_vflat"]) == 4
    assert (tmp_path / "results" / "rotation_curve.json").exists()
    assert (tmp_path / "papers" / "hi" / "figures" / "rotation_curve.pdf").exists()
    macros = (tmp_path / "papers" / "hi" / "generated" / "macros.tex").read_text()
    # namespaced: an offline run fills the synthetic side only
    assert r"\newcommand{\hiSynVflat}" in macros
    assert r"\newcommand{\hiRealVflat}{--}" in macros
    assert r"\hiSynThrMinusEdge" in macros and r"\hiSynKepler" in macros
    assert r"\hiSource" in macros


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
