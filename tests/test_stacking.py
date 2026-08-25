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
    assert abs(m["flux"] - 1.0) < 0.1 and m["rms"] > 0 and m["snr"] > 10
    assert m["peak_searched"] >= m["flux"]  # the searched max is only a labelled diagnostic


def test_injection_recovery_can_fail_now():
    rng = np.random.default_rng(2)
    bg = rng.normal(0, 0.1, (300, 41, 41))  # noise-only background
    cal = stacking.injection_recovery(bg, inject_amp=0.05)
    # near-unity for a well-centred population, but no longer the shift-equivariance identity:
    # the ratio varies across jitter draws and is not exactly 1.0
    assert abs(cal["ratio"] - 1.0) < 0.2
    assert cal["ratio"] != 1.0
    assert cal["ratio_sd"] > 0
    # a badly-centred population loses flux at the forced centre pixel -- the test CAN fail
    lossy = stacking.injection_recovery(bg, inject_amp=0.05, jitter_pix=2.0)
    assert lossy["ratio"] < 0.8


def test_forced_photometry_goes_negative_half_the_time():
    # CLAUDE.md: "Test a photometry routine on pure noise: forced should go negative about half
    # the time." The searched peak fails this by construction; the forced centre passes.
    rng = np.random.default_rng(7)
    forced = []
    searched = []
    for _ in range(60):
        img = rng.normal(0, 1.0, (51, 51))
        m = stacking.measure_stacked_flux(img)
        forced.append(m["flux"])
        searched.append(m["peak_searched"])
    frac_neg = np.mean(np.asarray(forced) < 0)
    assert 0.3 < frac_neg < 0.7
    assert np.mean(np.asarray(searched) > 0) > 0.9  # the searched peak is positive ~always


def test_individually_detected_flags_bright_sources():
    rng = np.random.default_rng(8)
    cube = rng.normal(0, 0.1, (50, 51, 51))
    cube[3] += stacking.gaussian_psf(51, 2.5, amp=2.0)  # a 20-sigma detection
    det = stacking.individually_detected(cube)
    assert det[3]
    assert det.sum() <= 2  # noise cutouts essentially never trip a 5-sigma searched test


def test_skewed_population_stack_is_the_median_not_the_mean():
    cube = stacking.synthetic_population(
        n_sources=800, source_flux=0.05, noise=0.08, flux_scatter_dex=0.6, seed=9
    )
    m = stacking.measure_stacked_flux(stacking.median_stack(cube))
    # log-normal: median 0.05, mean 0.05*10^(0.6^2 ln10 / 2) ~ 2.6x larger. The clipped-median
    # stack tracks the MEDIAN; calling it a "mean" would be wrong by that factor.
    assert abs(m["flux"] - 0.05) < 0.02
    assert m["flux"] < 0.6 * float(
        np.mean(0.05 * 10.0 ** np.random.default_rng(9).normal(0, 0.6, 800))
    )


def test_synthetic_population_stack_beats_noise():
    cube = stacking.synthetic_population(n_sources=600, source_flux=0.05, noise=0.12, seed=3)
    assert cube.shape == (600, 51, 51)
    # individually undetected (noise >> source) but the stack recovers the mean at high SNR
    assert np.std(cube[0]) > 1.5 * 0.05
    m = stacking.measure_stacked_flux(stacking.median_stack(cube))
    assert abs(m["flux"] - 0.05) < 0.02 and m["snr"] > 5


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
    assert bybright[0]["flux"] > bybright[1]["flux"]  # bright bin has more radio flux
    assert abs(bybright[0]["flux"] - 0.10) < 0.03 and abs(bybright[1]["flux"] - 0.04) < 0.02


def test_run_offline(tmp_path):
    m = stacking.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_stacked"] + m["n_detected_excluded"] == 600
    assert m["stacked_snr"] > 5
    assert abs(m["stacked_flux"] - m["injected_truth"]) < 0.02  # recovers the injected median
    inj = m["injection"]
    assert 0.7 < inj["ratio"] < 1.3 and inj["ratio_sd"] >= 0
    # the off-source control stack: identical pipeline on noise-only positions, consistent with 0
    ctl = m["control"]
    assert ctl["n"] == 600
    assert abs(ctl["flux"]) < 3.5 * ctl["rms"]
    assert abs(ctl["flux"]) < 0.5 * m["stacked_flux"]  # the science signal is not a pedestal
    assert m["n_bins"] == 3 and len(m["bins"]) == 3  # magnitude-binned trend produced
    assert m["n_zbins"] == 3 and len(m["zbins"]) == 3  # redshift-binned trend produced
    assert all("z_med" in b for b in m["zbins"])
    assert (tmp_path / "results" / "stacking_metrics.json").exists()
    assert (tmp_path / "papers" / "stacking" / "figures" / "stack.pdf").exists()
    macros = (tmp_path / "papers" / "stacking" / "generated" / "macros.tex").read_text()
    assert r"\stFlux" in macros and r"\stInjRatio" in macros and r"\stCtlFlux" in macros
    assert r"\stLowzFlux" in macros and r"\stHighzFlux" in macros and r"\stPeakzFlux" in macros
    assert r"\stNqueried" in macros and r"\stNdetExcl" in macros


def test_synthetic_run_does_not_clobber_real_artifacts(tmp_path):
    import json

    (tmp_path / "results").mkdir(parents=True)
    (tmp_path / "results" / "stacking_metrics.json").write_text(
        json.dumps({"source": "SDSS DR16Q x VLASS-SE @ (180.0,25.0) r=3 deg", "n_stacked": 236})
    )
    figdir = tmp_path / "papers" / "stacking" / "figures"
    figdir.mkdir(parents=True)
    (figdir / "stack.pdf").write_bytes(b"REAL")
    stacking.run(out=str(tmp_path), offline=True)
    assert (figdir / "stack.pdf").read_bytes() == b"REAL"
    kept = json.loads((tmp_path / "results" / "stacking_metrics.json").read_text())
    assert kept["n_stacked"] == 236
