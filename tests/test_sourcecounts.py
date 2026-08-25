"""Tests for jansky_research.sourcecounts -- NVSS 1.4 GHz Euclidean source counts. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import sourcecounts as sc


def test_hopkins_reference_matches_published_anchors():
    # the polynomial should reproduce the canonical 1.4 GHz normalised counts at known fluxes
    # (Hopkins et al. 2003): ~42 Jy^1.5/sr at 10 mJy, ~160 at 55 mJy
    assert 35 < sc.hopkins2003_counts(0.010) < 50
    assert 130 < sc.hopkins2003_counts(0.055) < 190
    # vectorised
    out = sc.hopkins2003_counts(np.array([0.01, 0.055]))
    assert out.shape == (2,)


def test_compute_counts_recovers_hopkins_on_synthetic():
    sky = sc.synthetic_sky(area_sr=0.05, seed=1)
    res = sc.compute_counts(sky["fluxes_jy"], sky["area_sr"], s_min_jy=0.0035)
    # drawn from Hopkins, so the measured normalised counts track Hopkins within the Poisson scatter
    assert 0.8 < res["ratio_med"] < 1.25
    assert res["ratio_scatter_dex"] < 0.15
    # the 1.4 GHz counts are sub-Euclidean at these fluxes (slope flatter than -2.5)
    assert -2.5 < res["slope_diff"] < -1.5


def test_compute_counts_small_sample_returns_gracefully():
    res = sc.compute_counts(np.array([1.0, 2.0, 3.0]), 0.01, s_min_jy=0.5)
    assert res["n_sources"] == 3
    assert res["ratio_med"] is None


def test_synthetic_sky_is_flux_limited():
    sky = sc.synthetic_sky(area_sr=0.05, s_min_jy=0.0035, s_max_jy=5.0, seed=0)
    f = sky["fluxes_jy"]
    assert f.min() >= 0.0035 and f.max() <= 5.0
    assert f.size > 100  # a populated sky


def test_run_offline_writes_outputs_and_recovers(tmp_path):
    m = sc.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_sources"] > 500
    assert 0.8 < m["hopkins_ratio_med"] < 1.25
    assert m["slope_diff"] < -1.5
    assert (tmp_path / "results" / "sourcecounts_metrics.json").exists()
    assert (tmp_path / "papers" / "sourcecounts" / "figures" / "sourcecounts.pdf").exists()
    macros = (tmp_path / "papers" / "sourcecounts" / "generated" / "macros.tex").read_text()
    assert r"\scRatio" in macros and r"\scSlope" in macros


def test_compute_counts_closed_form_power_law_normalisation():
    """The one genuinely new piece of code — the solid-angle normalisation — had no coverage:
    generator and comparator shared the same reference AND the same area expression, so a
    2x-wrong normalisation and even a flat reference passed every test. This closed form
    cannot be fooled: deterministic fluxes drawn from dN/dS = k S^-2.5 over a known area must
    return en = k in absolute units. It fails if the area, the S^2.5 factor, or the bin-width
    division is wrong."""
    k = 60.0  # Jy^1.5 / sr
    area = 0.05  # sr
    s_lo, s_hi = 0.004, 2.0
    # integral of k S^-2.5 dS from s to s_hi = (2k/3)(s^-1.5 - s_hi^-1.5); invert the CDF at
    # deterministic quantile midpoints so the counts carry no sampling noise
    n_tot = (2.0 * k / 3.0) * (s_lo**-1.5 - s_hi**-1.5) * area
    n = int(round(n_tot))
    q = (np.arange(n) + 0.5) / n
    inv = (s_lo**-1.5 - q * (s_lo**-1.5 - s_hi**-1.5)) ** (-1.0 / 1.5)
    res = sc.compute_counts(inv, area, s_min_jy=s_lo)
    good = res["good"]
    assert good.sum() >= 5
    ratio = res["en"][good] / k
    assert np.all(np.abs(ratio - 1.0) < 0.05), ratio


def test_fixed_bin_grid_is_independent_of_the_brightest_source():
    """The old geomspace-to-max grid made every statistic a function of one Poisson-random
    object; the fixed 0.2-dex grid anchored at the cut must leave the compared (sub-Jy) bins
    byte-identical when a 10 Jy source is appended."""
    sky = sc.synthetic_sky(seed=3)
    base = sc.compute_counts(sky["fluxes_jy"], sky["area_sr"], s_min_jy=0.0035)
    plus = sc.compute_counts(
        np.append(sky["fluxes_jy"], 10.0), sky["area_sr"], s_min_jy=0.0035
    )
    n = base["good"].sum()
    assert plus["good"].sum() == n  # the same sub-Jy bins are used
    np.testing.assert_allclose(base["en"][base["good"]], plus["en"][plus["good"]])
    assert base["ratio_med"] == plus["ratio_med"]
    assert base["ratio_scatter_dex"] == plus["ratio_scatter_dex"]


def test_merge_components_sums_fluxes():
    ra = np.array([180.0, 180.0 + 50.0 / 3600.0, 180.5, 181.0])  # first two within 60"
    dec = np.array([30.0, 30.0, 30.0, 30.0])
    fx = np.array([0.010, 0.005, 0.020, 0.001])
    merged = sc.merge_components(ra, dec, fx, link_arcsec=60.0)
    assert len(merged) == 3
    assert np.isclose(sorted(merged)[-1], 0.020) and np.isclose(sorted(merged)[1], 0.015)


def test_clustering_variance_scales_sensibly():
    lo = sc.clustering_variance_pct(8.0, amp=1.0e-3, n_pairs=100_000)
    hi = sc.clustering_variance_pct(8.0, amp=1.6e-3, n_pairs=100_000)
    small_field = sc.clustering_variance_pct(2.0, amp=1.0e-3, n_pairs=100_000)
    assert 0.5 < lo < 5.0  # percent-level for a 8-degree cone, as the referee computed
    assert hi > lo  # scales with the amplitude
    assert small_field > lo  # more variance in a smaller field


def test_run_offline_commits_the_budget_and_the_bins(tmp_path):
    m = sc.run(out=str(tmp_path), offline=True)
    assert m["chi2_unity"] is not None and m["chi2_with_ref_resid"] is not None
    assert m["chi2_with_ref_resid"] < m["chi2_unity"]  # the reference error can only help
    assert m["expected_scatter_dex"] >= m["poisson_scatter_dex"]
    assert m["ratio_med_boot_err"] is not None and m["slope_boot_err"] is not None
    assert m["slope_ref_same_bins"] is not None  # the recover-a-known's own comparand
    assert "3.5" in m["cut_sweep_mjy"] and "20" in m["cut_sweep_mjy"]
    assert m["cosmic_variance_pct_hi"] > m["cosmic_variance_pct_lo"]
    assert len(m["bins"]) == m["n_bins_used"] + sum(1 for b in m["bins"] if not b["used"])
    macros = (tmp_path / "papers" / "sourcecounts" / "generated" / "macros.tex").read_text()
    for name in (r"\scRatioErr", r"\scExpectedDex", r"\scChiRef", r"\scCvLo", r"\scSlopeRef"):
        assert name in macros, name
