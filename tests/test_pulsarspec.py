"""Tests for jansky_research.pulsarspec — pulsar radio spectral indices. No network."""

from __future__ import annotations

import json

import numpy as np
import pytest

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


def test_compare_subsamples_reports_sensitivity():
    """A null needs a resolution, not two means that happen to agree."""
    rng = np.random.default_rng(0)
    a = rng.normal(-1.6, 0.75, 43)
    b = rng.normal(-1.6, 0.75, 430)
    c = pulsarspec.compare_subsamples(a, b)
    assert c["n_a"] == 43 and c["n_b"] == 430
    # SE on the difference is dominated by the small arm: 0.75/sqrt(43) = 0.114
    assert 0.09 < c["se_diff"] < 0.15
    assert c["resolvable"] == pytest.approx(2.0 * c["se_diff"])
    assert c["n_sigma_observed"] == pytest.approx(abs(c["diff"]) / c["se_diff"])
    # The point of the number: a real offset below `resolvable` cannot be claimed either way.
    assert c["resolvable"] > 0.15


def test_compare_subsamples_degenerate_arms_do_not_fabricate_a_limit():
    c = pulsarspec.compare_subsamples([-1.7], [-1.6, -1.8, -1.9])
    assert c["n_a"] == 1
    assert not np.isfinite(c["se_diff"]) and not np.isfinite(c["resolvable"])


def test_offline_recovers_the_injected_msp_offset():
    """The fixture injects MSPs 0.2 flatter; nothing asserted it was ever recovered.

    This is the leg the Methods section claims is validated ("injected steep spectra AND a
    flatter-millisecond sub-population"). Without this test the split was exercised but never
    checked, and the real sample's own 2-sigma resolution is coarser than the injected offset --
    so the test that matters is whether the *fixture*, which has ~1500 sources, recovers it.
    """
    psr = pulsarspec.synthetic_field(n_sources=4000, seed=3)
    res = pulsarspec.find_spectra(psr)
    c = pulsarspec.compare_subsamples(res["alpha"][res["is_msp"]], res["alpha"][~res["is_msp"]])
    # Injected: alpha_msp = mean_alpha + 0.2, i.e. MSPs flatter (less negative).
    assert c["diff"] == pytest.approx(0.2, abs=0.12)
    # ...and at this size the fixture can actually see it, unlike the real catalogue.
    assert c["n_sigma_observed"] > 2.0
    assert c["resolvable"] < 0.2


def test_run_offline_emits_the_sensitivity_macros(tmp_path):
    pulsarspec.run(out=str(tmp_path), offline=True)
    macros = (tmp_path / "papers" / "pulsarspec" / "generated" / "macros.tex").read_text()
    for name in (
        r"\psrAlphaDiff",
        r"\psrAlphaDiffSE",
        r"\psrAlphaDiffSigma",
        r"\psrAlphaResolvable",
        r"\psrNnormal",
        r"\psrNcatalogue",
    ):
        assert name in macros, name
    m = json.loads((tmp_path / "results" / "pulsarspec_metrics.json").read_text())
    assert m["n_normal"] + m["n_msp"] == m["n"]
    assert m["alpha_resolvable_2sigma"] == pytest.approx(2 * m["alpha_diff_se"], abs=0.011)
