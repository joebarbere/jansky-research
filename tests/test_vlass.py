"""Tests for jansky_research.vlass — multi-epoch radio variability. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import vlass


def test_eta_v_constant_source():
    # a steady source measured within its errors: eta ~ 1, V small
    flux = np.array([10.0, 10.2, 9.8])
    err = np.array([0.3, 0.3, 0.3])
    assert vlass.eta_metric(flux, err) < 3.0
    assert vlass.v_metric(flux) < 0.05
    assert vlass.debiased_modulation_index(flux, err) >= 0.0


def test_eta_v_variable_source():
    # a large-amplitude variable with small errors: eta >> 1, V large
    flux = np.array([5.0, 20.0, 40.0])
    err = np.array([0.5, 0.5, 0.5])
    assert vlass.eta_metric(flux, err) > 100.0
    assert vlass.v_metric(flux) > 0.5
    m = vlass.variability_metrics(flux, err)
    assert m.p_value < 1e-3 and m.n_epochs == 3  # inconsistent with a constant flux
    assert m.m_debiased > 0.4


def test_metrics_edge_single_epoch():
    assert vlass.eta_metric(np.array([5.0]), np.array([0.5])) == 0.0
    assert vlass.v_metric(np.array([5.0])) == 0.0
    assert vlass.debiased_modulation_index(np.array([5.0]), np.array([0.5])) == 0.0


def test_crossmatch_epochs_positional():
    ra = [
        np.array([150.0, 151.0]),
        np.array([151.00001, 150.00001]),
    ]  # epoch 2 reordered, ~tiny offset
    dec = [np.array([20.0, 21.0]), np.array([21.00001, 20.00001])]
    flux = [np.array([10.0, 5.0]), np.array([6.0, 12.0])]
    err = [np.array([0.3, 0.2]), np.array([0.3, 0.3])]
    ra0, dec0, fmat, emat = vlass.crossmatch_epochs(ra, dec, flux, err)
    assert fmat.shape == (2, 2)
    # source 0 (150,20) should match epoch-2 entry with flux 12.0; source 1 -> 6.0
    assert np.isclose(fmat[0, 0], 10.0) and np.isclose(fmat[0, 1], 12.0)
    assert np.isclose(fmat[1, 0], 5.0) and np.isclose(fmat[1, 1], 6.0)


def test_crossmatch_no_match_is_nan():
    ra = [np.array([150.0]), np.array([200.0])]  # far apart -> no match within 2.5"
    dec = [np.array([20.0]), np.array([-10.0])]
    flux = [np.array([10.0]), np.array([5.0])]
    err = [np.array([0.3]), np.array([0.3])]
    _, _, fmat, _ = vlass.crossmatch_epochs(ra, dec, flux, err)
    assert np.isnan(fmat[0, 1])


def test_select_candidates_separates_injected():
    ra, dec, flux, err, truth = vlass.synthetic_epochs(n_sources=3000, var_fraction=0.05, seed=1)
    eta = np.array([vlass.eta_metric(f, e) for f, e in zip(flux, err, strict=True)])
    v = np.array([vlass.v_metric(f) for f in flux])
    mask, eta_thr, v_thr = vlass.select_candidates(eta, v, sigma=3.0)
    assert eta_thr > 1.0 and v_thr > 0.0
    recovered = np.sum(mask & truth) / max(truth.sum(), 1)
    contamination = np.sum(mask & ~truth) / max(mask.sum(), 1)
    # a conservative 3-sigma cut over only 3 noisy epochs recovers the majority of strong
    # injected variables at very low contamination (it is deliberately incomplete, not impure)
    assert recovered > 0.6
    assert contamination < 0.2


def test_apply_flux_scale_corrects_and_floors():
    f, e = vlass.apply_flux_scale(1, np.array([10.0]), np.array([0.5]))
    assert np.isclose(f[0], 11.3)  # epoch-1 peak-flux correction = 1.13
    assert e[0] > 0.5 * 1.13  # systematic floor added in quadrature
    f2, e2 = vlass.apply_flux_scale(3, np.array([10.0]), np.array([0.4]), sys_frac=0.0)
    assert np.isclose(f2[0], 10.31) and np.isclose(e2[0], 0.4 * 1.031)  # no floor -> just scaling


def test_flux_correction_removes_epoch_scale_offset():
    # a STEADY source whose raw per-epoch fluxes differ ONLY by the Quick-Look scale offset
    true = 10.0
    raw = np.array([true / vlass.VLASS_PEAK_CORRECTION[e] for e in (1, 2, 3)])
    err = 0.01 * raw
    eta_raw = vlass.eta_metric(raw, err)  # uncorrected: the scale offset masquerades as variability
    cf, ce = [], []
    for ep, fr, er in zip((1, 2, 3), raw, err, strict=True):
        a, b = vlass.apply_flux_scale(ep, np.array([fr]), np.array([er]), sys_frac=0.0)
        cf.append(a[0])
        ce.append(b[0])
    eta_corr = vlass.eta_metric(np.array(cf), np.array(ce))
    assert eta_raw > 10.0  # a steady source would be wrongly flagged without the correction
    assert eta_corr < 1e-6  # the correction collapses it back to a constant


def test_isolated_mask_flags_crowded():
    # three well-separated sources + one 4" neighbour pair (deblending-prone)
    ra = np.array([150.0, 150.5, 151.0, 151.0 + 4.0 / 3600.0])
    dec = np.array([20.0, 20.5, 21.0, 21.0])
    iso = vlass.isolated_mask(ra, dec, radius_arcsec=5.0)
    assert iso[0] and iso[1]  # isolated
    assert not iso[2] and not iso[3]  # the 4" pair is flagged crowded
    assert vlass.isolated_mask(np.array([1.0]), np.array([1.0]))[0]  # single source is isolated


def test_measure_image_flux_forced_photometry():
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.crpix = [50, 50]
    w.wcs.cdelt = [-1 / 3600.0, 1 / 3600.0]  # 1 arcsec/pixel
    w.wcs.crval = [150.0, 20.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    yy, xx = np.mgrid[0:100, 0:100]
    img = 5.0 * np.exp(-((xx - 49) ** 2 + (yy - 49) ** 2) / (2 * 1.5**2))  # 5 mJy source at crval
    img = img + np.random.default_rng(0).normal(0.0, 0.1, img.shape)
    peak, rms, off = vlass.measure_image_flux(img, w, 150.0, 20.0, search_arcsec=4.0)
    assert abs(peak - 5.0) < 0.4  # recovers the source peak
    assert off < 2.5  # peak found at the locked position
    assert 0.05 < rms < 0.2  # local noise estimated from the annulus
    # off-source position (no source within the search box) -> only noise, no 5 mJy peak
    p2, _, _ = vlass.measure_image_flux(img, w, 150.0 + 20 / 3600.0, 20.0, search_arcsec=3.0)
    assert p2 < 1.0


def test_run_offline(tmp_path):
    m = vlass.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_candidates"] >= 1
    assert m["recovered_fraction"] > 0.55
    assert m["false_positive_fraction"] < 0.2
    assert (tmp_path / "results" / "vlass_metrics.json").exists()
    assert (tmp_path / "papers" / "vlass" / "figures" / "eta_v.pdf").exists()
