"""Tests for jansky_research.rmdipole -- the RM dipole/anisotropy test. Offline."""

from __future__ import annotations

import json

import numpy as np
import pytest

from jansky_research import rmdipole as rmd


def test_dipole_fit_exact_on_noiseless_model():
    # y built exactly from m + p.n must return m, p, amp, apex to machine-ish precision
    rng = np.random.default_rng(0)
    ra, dec = rng.uniform(0, 360, 4000), np.degrees(np.arcsin(rng.uniform(-1, 0.7, 4000)))
    n = rmd._unit_vectors(ra, dec)
    p_true = np.array([0.3, -0.2, 0.1])
    y = 2.0 + n @ p_true
    out = rmd.dipole_fit(ra, dec, y)
    assert out["monopole"] == pytest.approx(2.0, abs=1e-9)
    assert out["amp"] == pytest.approx(np.linalg.norm(p_true) / 2.0, abs=1e-9)
    apex = rmd._unit_vectors(np.array([out["apex_ra"]]), np.array([out["apex_dec"]]))[0]
    assert apex @ (p_true / np.linalg.norm(p_true)) > 0.99999


def test_synthetic_dipole_round_trip_amplitude_and_direction():
    syn = rmd.synthetic_dipole_catalogue(n_sources=60000, amp=0.4, seed=2)
    fit = rmd.fit_dipole(syn["ra"], syn["dec"], syn["resid"], syn["rm_err"], n_boot=40)
    assert fit["amp"] == pytest.approx(0.4, abs=3 * max(fit["amp_se"], 0.02))
    sep = np.degrees(
        np.arccos(
            rmd._unit_vectors(np.array([fit["apex_ra"]]), np.array([fit["apex_dec"]]))[0]
            @ rmd._unit_vectors(np.array([syn["true_apex"][0]]), np.array([syn["true_apex"][1]]))[0]
        )
    )
    assert sep < 15.0


def test_power_stat_debiases_measurement_noise():
    # amp=0: the debiased power monopole must recover sigma0^2, not sigma0^2 + noise^2
    syn = rmd.synthetic_dipole_catalogue(n_sources=40000, amp=0.0, sigma0=10.0, noise=6.0, seed=3)
    fit = rmd.fit_dipole(syn["ra"], syn["dec"], syn["resid"], syn["rm_err"], stat="power")
    assert fit["monopole"] == pytest.approx(100.0, rel=0.05)


def test_scramble_null_flat_without_dipole_and_kills_injected_one():
    quiet = rmd.synthetic_dipole_catalogue(n_sources=8000, amp=0.0, seed=4)
    fit0 = rmd.fit_dipole(quiet["ra"], quiet["dec"], quiet["resid"], quiet["rm_err"], n_boot=0)
    null0 = rmd.footprint_scramble_null(
        quiet["ra"], quiet["dec"], fit0["y"], amp_obs=fit0["amp"], n_scramble=60, seed=4
    )
    assert null0["p_value"] > 0.02  # no false detection on the quiet sky
    loud = rmd.synthetic_dipole_catalogue(n_sources=8000, amp=0.6, seed=5)
    fit1 = rmd.fit_dipole(loud["ra"], loud["dec"], loud["resid"], loud["rm_err"], n_boot=0)
    null1 = rmd.footprint_scramble_null(
        loud["ra"], loud["dec"], fit1["y"], amp_obs=fit1["amp"], n_scramble=60, seed=5
    )
    assert null1["p_value"] < 0.05  # the injected dipole must beat the footprint null


def test_extragalactic_residuals_nn_and_latitude_paths():
    n = 400
    rng = np.random.default_rng(6)
    s = {
        "ra": rng.uniform(0, 360, n),
        "dec": rng.uniform(-60, 40, n),
        "gal_b": rng.uniform(-80, 80, n),
        "rm": rng.normal(10.0, 5.0, n),  # 10 rad/m^2 foreground offset
        "rm_err": np.full(n, 1.0),
        "nn_rm_med": np.full(n, 10.0),
        "nn_rm_count": np.full(n, 9),
    }
    out_nn = rmd.extragalactic_residuals(s, b_min=30.0, method="nn")
    assert np.all(np.abs(out_nn["gal_b"]) >= 30.0)
    assert abs(np.mean(out_nn["resid"])) < 1.5  # the 10 rad/m^2 offset is subtracted
    out_lat = rmd.extragalactic_residuals(s, b_min=30.0, method="latitude")
    assert abs(np.median(out_lat["resid"])) < 1.5
    with pytest.raises(ValueError):
        rmd.extragalactic_residuals(s, method="nope")


def test_extragalactic_residuals_drops_thin_nn():
    s = {
        "ra": np.array([1.0, 2.0]),
        "dec": np.array([0.0, 0.0]),
        "gal_b": np.array([50.0, 50.0]),
        "rm": np.array([5.0, 5.0]),
        "rm_err": np.array([1.0, 1.0]),
        "nn_rm_med": np.array([1.0, np.nan]),
        "nn_rm_count": np.array([9, 9]),
    }
    out = rmd.extragalactic_residuals(s, method="nn")
    assert out["resid"].size == 1


def test_fit_dipole_rejects_unknown_stat():
    syn = rmd.synthetic_dipole_catalogue(n_sources=100, seed=7)
    with pytest.raises(ValueError):
        rmd.fit_dipole(syn["ra"], syn["dec"], syn["resid"], syn["rm_err"], stat="median")


def test_noise_stat_fits_the_error_map_not_the_sky():
    # a dipolar sigma map with isotropic residuals: 'noise' finds it, debiased 'power' must not
    rng = np.random.default_rng(9)
    n = 30000
    ra = rng.uniform(0, 360, n)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 0.75, n)))
    mu = rmd._unit_vectors(ra, dec) @ rmd._unit_vectors(np.array([170.0]), np.array([-7.0]))[0]
    err = np.sqrt(4.0 * (1.0 + 0.5 * mu))
    resid = rng.normal(0.0, 8.0, n) + rng.normal(0.0, err)
    noise_fit = rmd.fit_dipole(ra, dec, resid, err, stat="noise", n_boot=0)
    assert noise_fit["amp"] == pytest.approx(0.5, abs=0.05)
    power_fit = rmd.fit_dipole(ra, dec, resid, err, stat="power", n_boot=40)
    assert power_fit["amp"] < 3 * max(power_fit["amp_se"], 0.01)


def test_clip_quantile_tames_an_outlier_driven_dipole():
    # isotropic core + a handful of huge outliers piled at one apex: clipping kills the "dipole"
    rng = np.random.default_rng(10)
    n = 20000
    ra = rng.uniform(0, 360, n)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 0.75, n)))
    resid = rng.normal(0.0, 5.0, n)
    hot = slice(0, 200)
    ra[hot], dec[hot] = 170.0 + rng.normal(0, 5, 200), -7.0 + rng.normal(0, 5, 200)
    resid[hot] = 200.0
    err = np.full(n, 1.0)
    raw = rmd.fit_dipole(ra, dec, resid, err, stat="power", n_boot=0)
    clipped = rmd.fit_dipole(ra, dec, resid, err, stat="power", clip_quantile=0.98, n_boot=0)
    assert clipped["n_clipped"] >= 200
    assert clipped["amp"] < 0.3 * raw["amp"]
    assert clipped["ra_used"].size == clipped["y"].size < n


def test_compare_directions_zero_at_cmb_apex():
    out = rmd.compare_directions(*rmd.CMB_DIPOLE_RA_DEC)
    assert out["sep_cmb_deg"] == pytest.approx(0.0, abs=1e-9)
    anti = rmd.compare_directions(rmd.CMB_DIPOLE_RA_DEC[0] - 180.0, -rmd.CMB_DIPOLE_RA_DEC[1])
    assert anti["sep_cmb_deg"] == pytest.approx(180.0, abs=1e-9)


def test_synthetic_catalogue_respects_footprint_and_real_positions():
    syn = rmd.synthetic_dipole_catalogue(n_sources=3000, dec_max=49.0, seed=8)
    assert syn["dec"].max() <= 49.0
    ra0, dec0 = np.array([10.0, 20.0, 30.0]), np.array([-10.0, 0.0, 10.0])
    on_real = rmd.synthetic_dipole_catalogue(ra_deg=ra0, dec_deg=dec0, seed=8)
    assert np.array_equal(on_real["ra"], ra0) and on_real["resid"].size == 3


def test_run_offline_writes_artifacts(tmp_path):
    m = rmd.run(str(tmp_path), offline=True, n_scramble=40)
    assert m["source"].startswith("synthetic dipole")
    lead, ctrl = m["legs"][0], m["legs"][1]
    # recover-a-known: injected amp within 3 sigma-ish, apex within 15 deg, null control quiet
    assert lead["amp"] == pytest.approx(m["true_amp"], abs=0.06)
    assert m["sep_true_deg"] < 15.0
    assert lead["p_scramble"] < 0.05 < ctrl["p_scramble"]
    saved = json.loads((tmp_path / "results" / "rmdipole_metrics.json").read_text())
    assert saved == m
    assert (tmp_path / "papers" / "rmdipole" / "figures" / "rmdipole.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "rmdipole" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\rmdSynAmp}" in macros and r"\newcommand{\rmdRealAmp}{--}" in macros
    table = (tmp_path / "papers" / "rmdipole" / "generated" / "legs_table.tex").read_text()
    assert table.count(r"\\") == len(m["legs"]) and "power" in table


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    rmd._write_macros({"source": "x", "is_real": True, "legs": [{"amp": None}]}, p)
    txt = p.read_text()
    assert r"\newcommand{\rmdRealAmp}{--}" in txt and r"\newcommand{\rmdSynAmp}{--}" in txt
