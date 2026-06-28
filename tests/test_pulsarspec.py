"""Tests for jansky_research.pulsarspec — pulsar radio spectral indices. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import pulsarspec


def test_pulsar_alpha_textbook():
    # S ~ nu^-1.8 between 400 and 1400 MHz
    s400 = np.array([10.0])
    s1400 = s400 * (1.4 / 0.4) ** -1.8
    a, _ = pulsarspec.pulsar_alpha(s400, s1400)
    assert np.isclose(a[0], -1.8, atol=1e-6)


def test_is_millisecond():
    p = np.array([0.003, 0.02, 0.5, 1.0])
    assert pulsarspec.is_millisecond(p).tolist() == [True, True, False, False]
    assert pulsarspec.is_millisecond(np.array([0.05]), p_max=0.1).tolist() == [True]


def test_spectral_distribution():
    a = np.array([-1.0, -2.0, -3.0, np.nan])
    s = pulsarspec.spectral_distribution(a)
    assert s["n"] == 3 and np.isclose(s["mean"], -2.0) and np.isclose(s["median"], -2.0)


def test_find_spectra_skips_incomplete():
    psr = {
        "p0": np.array([0.5, 0.002, 1.0]),
        "s400": np.array([10.0, 5.0, np.nan]),  # 3rd has no S400 -> skipped
        "s1400": np.array([2.0, 1.5, 3.0]),
    }
    res = pulsarspec.find_spectra(psr)
    assert res["alpha"].size == 2  # only the two with both fluxes
    assert res["is_msp"].tolist() == [False, True]


def test_run_offline(tmp_path):
    m = pulsarspec.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n"] > 500
    assert -2.2 < m["mean_alpha"] < -1.4  # recovers the injected steep mean
    assert m["n_msp"] >= 1
    assert (tmp_path / "results" / "pulsarspec_metrics.json").exists()
    assert (tmp_path / "papers" / "pulsarspec" / "figures" / "spectra.pdf").exists()
    macros = (tmp_path / "papers" / "pulsarspec" / "generated" / "macros.tex").read_text()
    assert r"\psrMeanAlpha" in macros and r"\psrMeanAlphaMsp" in macros
