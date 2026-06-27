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
    longs = np.array([20.0, 40.0, 60.0, 80.0])
    slices = [hi.synthetic_lv_slice(ell, v_flat=230.0, seed=i) for i, ell in enumerate(longs)]
    R, V = hi.rotation_curve(longs, slices)
    assert np.all(np.diff(R) > 0)  # sorted by radius
    assert np.all(np.abs(V - 230.0) < 6.0)  # recovers the injected flat curve


def test_run_offline(tmp_path):
    m = hi.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert abs(m["V_flat_mean_kms"] - 230.0) < 6.0
    assert m["V_flat_std_kms"] < 6.0  # flat
    assert (tmp_path / "results" / "rotation_curve.json").exists()
    assert (tmp_path / "paper" / "figures" / "rotation_curve.pdf").exists()


def test_read_lab_slice_parses_wcs(tmp_path):
    from astropy.io import fits

    nb, nv = 5, 11
    data = np.arange(nb * nv, dtype=float).reshape(1, nb, nv)
    h = fits.PrimaryHDU(data).header
    h["CTYPE1"], h["CRVAL1"], h["CDELT1"], h["CRPIX1"] = "VELO-LSR", 0.0, 20000.0, 1  # m/s, 20 km/s/ch
    h["CUNIT1"] = "M/S"
    h["CTYPE2"], h["CRVAL2"], h["CDELT2"], h["CRPIX2"] = "GLAT-CAR", -2.0, 1.0, 1
    p = tmp_path / "slice.fits"
    fits.PrimaryHDU(data, h).writeto(p)
    lat, vel, d = hi.read_lab_slice(p)
    assert d.shape == (nb, nv)
    assert np.isclose(vel[0], 0.0) and np.isclose(vel[1], 20.0)  # 20000 m/s -> 20 km/s
    assert np.isclose(lat[0], -2.0) and np.isclose(lat[-1], 2.0)
