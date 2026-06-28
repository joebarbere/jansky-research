"""Tests for jansky_research.southern — GLEAM-X+RACS multi-band curvature. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import southern


def test_fit_log_parabola_recovers_turnover():
    nu = np.concatenate([southern.GLEAMX_NU_GHZ, southern.RACS_NU_GHZ])
    # a log-parabola peaking at 0.5 GHz
    flux = 10.0 ** (2.0 - 1.0 * (np.log10(nu / 0.5)) ** 2)
    fit = southern.fit_log_parabola(nu, flux, 0.05 * flux)
    assert fit["is_peaked"]
    assert abs(fit["nu_pk_ghz"] - 0.5) < 0.05  # measured turnover, not a bound
    assert fit["a"] < 0  # concave
    assert fit["n_points"] == nu.size


def test_fit_log_parabola_power_law_not_peaked():
    nu = np.concatenate([southern.GLEAMX_NU_GHZ, southern.RACS_NU_GHZ])
    flux = (nu / 0.2) ** -0.8  # straight power law
    fit = southern.fit_log_parabola(nu, flux)
    assert not fit["is_peaked"]
    # too few points -> graceful nan/False
    fit2 = southern.fit_log_parabola(nu[:2], flux[:2])
    assert not fit2["is_peaked"] and fit2["n_points"] == 2


def test_classify_curved():
    peaked_fit = {"is_peaked": True}
    assert southern.classify_curved(peaked_fit, -0.5, -0.7) == "peaked"
    flat_fit = {"is_peaked": False}
    assert southern.classify_curved(flat_fit, -1.4, -1.4) == "uss"  # ultra-steep throughout
    assert southern.classify_curved(flat_fit, -0.7, -0.7) == "steep"
    assert southern.classify_curved(flat_fit, np.nan, -0.7) == "nan"


def test_find_peaked_south_recovers_injected():
    from jansky_research.spectra import crossmatch

    gleamx, racs, truth_pk, truth_uss = southern.synthetic_field(n_sources=1500, seed=1)
    res = southern.find_peaked_south(gleamx, racs)
    pk = res["is_peaked"]
    assert pk.sum() >= 1
    # recovery: injected-peaked positions landing on a flagged peaked candidate
    rp = np.flatnonzero(truth_pk)
    i, _, _ = crossmatch(gleamx["ra"][rp], gleamx["dec"][rp], res["ra"][pk], res["dec"][pk], 5.0)
    assert i.size / truth_pk.sum() > 0.6
    # purity
    j, _, _ = crossmatch(res["ra"][pk], res["dec"][pk], gleamx["ra"][rp], gleamx["dec"][rp], 5.0)
    assert j.size / pk.sum() > 0.6
    # USS sources are flagged
    assert res["is_uss"].sum() >= 1


def test_run_offline(tmp_path):
    m = southern.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_matched"] > 1000
    assert m["n_peaked"] >= 1
    assert m["n_peaked_recovered"] >= 1
    assert m["n_peaked_recovered"] <= m["n_injected_peaked"]
    assert m["median_nu_pk_mhz"] > 0  # a measured median turnover frequency
    assert (tmp_path / "results" / "southern_metrics.json").exists()
    assert (tmp_path / "papers" / "southern" / "figures" / "seds.pdf").exists()
    macros = (tmp_path / "papers" / "southern" / "generated" / "macros.tex").read_text()
    assert r"\soNpeaked" in macros and r"\soMedianNupk" in macros
