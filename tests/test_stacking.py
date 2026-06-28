"""Tests for jansky_research.stacking — sub-threshold radio stacking. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import stacking


def test_gaussian_psf():
    psf = stacking.gaussian_psf(21, fwhm_pix=3.0, amp=2.0)
    assert psf.shape == (21, 21)
    assert np.isclose(psf[10, 10], 2.0)  # peak at centre
    assert psf.max() <= 2.0 + 1e-9


def test_median_stack_recovers_and_rejects_outlier():
    rng = np.random.default_rng(0)
    truth = stacking.gaussian_psf(41, 3.0, amp=1.0)
    cube = truth[None, :, :] + rng.normal(0, 0.5, (200, 41, 41))
    stack = stacking.median_stack(cube)
    assert abs(stack[20, 20] - 1.0) < 0.15  # central value recovered
    # a single very bright interloper cutout is rejected by the sigma-clip
    cube[0] += 1000.0
    stack2 = stacking.median_stack(cube)
    assert abs(stack2[20, 20] - 1.0) < 0.2


def test_measure_stacked_flux():
    rng = np.random.default_rng(1)
    img = stacking.gaussian_psf(51, 2.5, amp=1.0) + rng.normal(0, 0.02, (51, 51))
    m = stacking.measure_stacked_flux(img)
    assert abs(m["peak"] - 1.0) < 0.1 and m["rms"] > 0 and m["snr"] > 10


def test_injection_recovery_unbiased_for_clean_model():
    rng = np.random.default_rng(2)
    bg = rng.normal(0, 0.1, (300, 41, 41))  # noise-only background
    cal = stacking.injection_recovery(bg, inject_amp=0.05)
    assert abs(cal["ratio"] - 1.0) < 0.2  # clean Gaussian model -> near-unity recovery


def test_synthetic_population_stack_beats_noise():
    cube = stacking.synthetic_population(n_sources=600, source_flux=0.05, noise=0.12, seed=3)
    assert cube.shape == (600, 51, 51)
    # individually undetected (noise >> source) but the stack recovers the mean at high SNR
    assert np.std(cube[0]) > 1.5 * 0.05
    m = stacking.measure_stacked_flux(stacking.median_stack(cube))
    assert abs(m["peak"] - 0.05) < 0.02 and m["snr"] > 5


def test_stack_in_bins_recovers_trend():
    # a faint and a bright sub-population, tagged by a binning value
    faint = stacking.synthetic_population(400, source_flux=0.04, seed=1)
    bright = stacking.synthetic_population(400, source_flux=0.10, seed=2)
    cube = np.concatenate([faint, bright])
    values = np.concatenate(
        [np.full(400, 20.5), np.full(400, 18.5)]
    )  # mags: faint=20.5, bright=18.5
    bins = stacking.stack_in_bins(cube, values, n_bins=2)
    assert len(bins) == 2
    bybright = sorted(bins, key=lambda b: b["value_med"])  # brightest (low mag) first
    assert bybright[0]["debiased"] > bybright[1]["debiased"]  # bright bin has more radio flux
    assert abs(bybright[0]["debiased"] - 0.10) < 0.03 and abs(bybright[1]["debiased"] - 0.04) < 0.02


def test_run_offline(tmp_path):
    m = stacking.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_stacked"] == 600
    assert m["stacked_snr"] > 5
    assert abs(m["debiased_flux"] - m["injected_truth"]) < 0.02  # recovers the injected mean
    assert 0.7 < m["recovery_ratio"] < 1.3
    assert m["n_bins"] == 3 and len(m["bins"]) == 3  # magnitude-binned trend produced
    assert m["n_zbins"] == 3 and len(m["zbins"]) == 3  # redshift-binned trend produced
    assert all("z_med" in b for b in m["zbins"])
    assert (tmp_path / "results" / "stacking_metrics.json").exists()
    assert (tmp_path / "papers" / "stacking" / "figures" / "stack.pdf").exists()
    macros = (tmp_path / "papers" / "stacking" / "generated" / "macros.tex").read_text()
    assert r"\stDebiased" in macros and r"\stRatio" in macros
    assert r"\stLowzFlux" in macros and r"\stHighzFlux" in macros and r"\stPeakzFlux" in macros
