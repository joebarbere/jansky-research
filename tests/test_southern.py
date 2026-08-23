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
    assert southern.classify_curved(peaked_fit, 0.5, -0.7) == "peaked"  # concave + rising low side
    # a concave fit but a FALLING low side is a steep-source artefact, not a real peak
    assert southern.classify_curved(peaked_fit, -0.8, -0.7) == "steep"
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
    # An offline run fills the SYNTHETIC namespace only. This test previously asserted the
    # un-namespaced names, i.e. it required the behaviour that let an offline rebuild write
    # synthetic counts into the macros the paper uses for the real cone.
    assert r"\soSynNpeaked" in macros and r"\soSynMedianNupk" in macros
    assert r"\newcommand{\soRealNpeaked}{--}" in macros


def test_macros_are_namespaced_and_merged(tmp_path):
    """An offline rebuild must not touch the real namespace.

    `make figures` runs every slice offline in the repo root. Before this, southern's macros
    were shared between modes and unmerged, so one such run replaced 1545 real matches with
    the synthetic field's count and wrote `\\soCallTried{0}` over a real 50 -- the documented
    `\\tiiNEvents` clobber, which ships a wrong number rather than a hole.
    """
    mac = tmp_path / "macros.tex"
    real = {
        "source": "GLEAM-X DR2 x RACS @ (30.0,-30.0) r=3.0deg",
        "n_matched": 1545,
        "n_peaked": 90,
        "n_uss": 59,
        "n_extended": 28,
        "median_nu_pk_mhz": 210.9,
        "call_tried": 50,
        "call_recovered": 38,
        "call_dlog": 0.112,
        "call_within_pct": 92,
    }
    southern._write_macros(real, mac)
    assert r"\newcommand{\soRealNmatched}{1545}" in mac.read_text()

    southern._write_macros(
        {
            "source": "synthetic",
            "n_matched": 40,
            "n_peaked": 6,
            "n_uss": 4,
            "n_extended": 2,
            "median_nu_pk_mhz": 150.0,
        },
        mac,
    )
    text = mac.read_text()
    assert r"\newcommand{\soRealNmatched}{1545}" in text, "offline run clobbered the real count"
    assert r"\newcommand{\soRealCallTried}{50}" in text, "offline run blanked the validation"
    assert r"\newcommand{\soSynNmatched}{40}" in text, "synthetic count lost its own namespace"


def test_callingham_macros_placeholder_rather_than_zero():
    """The validation runs on the real path only; a 0 there is a wrong number, not a hole."""
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        mac = Path(d) / "m.tex"
        southern._write_macros(
            {
                "source": "synthetic",
                "n_matched": 1,
                "n_peaked": 1,
                "n_uss": 1,
                "n_extended": 1,
                "median_nu_pk_mhz": 1.0,
            },
            mac,
        )
        assert r"\newcommand{\soRealCallTried}{--}" in mac.read_text()
