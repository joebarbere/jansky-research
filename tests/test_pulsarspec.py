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
    # The fixture injects mean_alpha = -1.8 with a +0.2-flatter 10% MSP arm, so the expected
    # population mean is ~-1.78. A window of a few SE (SE ~ 0.02 at this n) can actually fail
    # on a wrong transform; the old +/-0.4 window was ~28 SE wide and could not.
    assert m["mean_alpha"] == pytest.approx(-1.78, abs=0.06)
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


def test_find_spectra_counts_the_nonpositive_flux_rejections():
    """J0540-6919 tabulates S400 = 0.0; the positivity guard dropped it invisibly."""
    psr = {
        "p0": np.array([0.05, 0.5, 0.7]),
        "s400": np.array([0.0, 10.0, 5.0]),  # first row: tabulated but non-positive
        "s1400": np.array([24.0, 2.0, 1.0]),
    }
    res = pulsarspec.find_spectra(psr)
    assert res["alpha"].size == 2
    assert res["n_rejected_nonpositive_flux"] == 1


def test_completeness_limited_removes_the_joint_detection_bias():
    """Steep sources near the S400 end are lost to the S1400 floor; the cut must steepen the mean.

    Build a scale-free steep population, impose a hard S1400 floor (the joint-detection
    truncation), and check that (a) the truncated sample's mean is biased flat relative to the
    injected mean, and (b) the completeness-limited subsample recovers it.
    """
    rng = np.random.default_rng(1)
    n = 20000
    alpha = rng.normal(-1.8, 0.5, n)
    s400 = 10.0 ** rng.uniform(-1.0, 2.0, n)
    s1400 = s400 * (1.4 / 0.4) ** alpha
    floor = 0.1
    kept = s1400 >= floor  # the truncation a joint catalogue applies
    biased = float(np.mean(alpha[kept]))
    assert biased > -1.8  # flat bias, as the paper asserts
    comp = pulsarspec.completeness_limited(
        alpha[kept], s400[kept], s1400_floor_mjy=floor, alpha_min=-3.0
    )
    assert comp["s400_cut_mjy"] == pytest.approx(floor * (3.5**3.0), rel=0.01)
    assert comp["n"] > 100
    # above the cut the truncation cannot operate for any alpha >= -3, so the mean returns
    # to the injected value within a few SE, and is steeper than the biased mean
    assert comp["mean"] < biased
    assert comp["mean"] == pytest.approx(-1.8, abs=0.05)


def test_msp_cut_sweep_reassigns_membership():
    psr = pulsarspec.synthetic_field(n_sources=3000, seed=5)
    res = pulsarspec.find_spectra(psr)
    sweep = pulsarspec.msp_cut_sweep(res["alpha"], res["period_s"])
    ns = [sweep[k]["n_msp"] for k in ("0.01", "0.02", "0.03", "0.05", "0.1")]
    assert all(a <= b for a, b in zip(ns, ns[1:], strict=False))  # wider cut, more MSPs
    assert sweep["0.03"]["two_se"] is not None


def test_permutation_pvalue_calibrates():
    rng = np.random.default_rng(2)
    a = rng.normal(-1.8, 0.5, 500)
    is_msp = np.zeros(500, bool)
    is_msp[:50] = True  # a random split: no real difference
    p_null = pulsarspec.permutation_pvalue(a, is_msp, n_perm=2000, seed=0)
    assert p_null > 0.05  # no false detection on an exchangeable split
    a2 = a.copy()
    a2[is_msp] += 0.5  # a large injected offset must be detected
    p_alt = pulsarspec.permutation_pvalue(a2, is_msp, n_perm=2000, seed=0)
    assert p_alt < 0.01
    assert np.isnan(pulsarspec.permutation_pvalue([-1.0], [True]))


def test_run_with_names_writes_the_source_table(tmp_path, monkeypatch):
    psr = pulsarspec.synthetic_field(n_sources=400, seed=7)
    psr["name"] = np.asarray([f"J{i:04d}+0000" for i in range(400)], dtype=object)
    psr["catalogue_version"] = "test fixture"
    psr["fetched_utc"] = "2026-08-25T00:00:00Z"
    monkeypatch.setattr(pulsarspec, "fetch_atnf", lambda: psr)
    m = pulsarspec.run(out=str(tmp_path), offline=False)
    rows = (tmp_path / "results" / "pulsarspec_sources.csv").read_text().strip().splitlines()
    assert rows[0].startswith("psrj,")
    assert len(rows) - 1 == m["n"]
    assert m["catalogue_version"] == "test fixture"
    assert m["n_msp_catalogue"] + m["n_normal_catalogue"] <= 400
    macros = (tmp_path / "papers" / "pulsarspec" / "generated" / "macros.tex").read_text()
    for name in (r"\psrCIlo", r"\psrCompleteMean", r"\psrPermP", r"\psrNSfourHundred"):
        assert name in macros, name


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
