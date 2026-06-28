"""Tests for jansky_research.peaked — peaked-spectrum (GPS/CSS) selection. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import peaked


def test_classify_sed():
    assert peaked.classify_sed(0.5, -0.6) == "peaked"  # rising then falling -> turnover in band
    assert peaked.classify_sed(-0.8, -0.8) == "steep"  # falling throughout
    assert peaked.classify_sed(0.4, 0.3) == "inverted"  # still rising at 3 GHz
    assert peaked.classify_sed(-0.2, -0.2) == "flat"
    assert peaked.classify_sed(np.nan, -0.6) == "nan"


def test_two_point_indices_textbook():
    # a flat power law S ~ nu^-0.7 gives alpha_low = alpha_high = -0.7
    nu = peaked.NU_GHZ
    s_t = (nu["tgss"] / nu["nvss"]) ** -0.7
    s_n = 1.0
    s_v = (nu["vlass"] / nu["nvss"]) ** -0.7
    al, ah = peaked.two_point_indices(np.array([s_t]), np.array([s_n]), np.array([s_v]))
    assert np.isclose(al[0], -0.7, atol=1e-6) and np.isclose(ah[0], -0.7, atol=1e-6)


def test_peak_frequency_recovers_turnover():
    # a parabola in log-log peaking at 1 GHz
    nu = np.array([0.15, 1.0, 5.0])
    flux = 10.0 ** (1.0 - 0.8 * (np.log10(nu / 1.0)) ** 2)
    nu_peak, is_peak = peaked.peak_frequency(flux, nu)
    assert is_peak and abs(nu_peak - 1.0) < 0.05
    # a monotonic power law is not a peak
    _, is_peak2 = peaked.peak_frequency(nu**-0.7, nu)
    assert not is_peak2


def test_find_peaked_recovers_injected():
    tgss, nvss, vlass, truth = peaked.synthetic_field(n_sources=2000, peaked_fraction=0.05, seed=1)
    res = peaked.find_peaked(tgss, nvss, vlass)
    assert res["ra"].size > 1500  # most sources triple-matched
    peaked_mask = res["cls"] == "peaked"
    # the peaked class is enriched and not dominated by the steep majority
    assert peaked_mask.sum() >= 1
    frac_peaked = peaked_mask.mean()
    assert frac_peaked < 0.2  # selection is selective, not flagging everything
    # curvature of peaked candidates is negative (concave: rises then falls)
    assert np.median(res["curvature"][peaked_mask]) < 0.0


def test_run_offline(tmp_path):
    m = peaked.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_matched"] > 1000  # most of the synthetic field triple-matches
    assert m["n_peaked"] >= 1
    assert m["n_injected_peaked"] >= 1
    assert (tmp_path / "results" / "peaked_metrics.json").exists()
    assert (tmp_path / "papers" / "peaked" / "figures" / "curvature.pdf").exists()
    macros = (tmp_path / "papers" / "peaked" / "generated" / "macros.tex").read_text()
    assert r"\pkNpeaked" in macros
