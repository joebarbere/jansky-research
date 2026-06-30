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
