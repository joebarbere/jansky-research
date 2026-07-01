"""Tests for jansky_research.stokesv — Stokes-V coherent-emitter selection. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import stokesv
from jansky_research.stokesv import match_targets_to_radio, measure_circular_pol


def test_fractional_circular_pol():
    frac, err = stokesv.fractional_circular_pol(
        np.array([-3.0, 0.5]), np.array([10.0, 10.0]), np.array([0.2, 0.2]), np.array([0.5, 0.5])
    )
    # |V|/I uses the magnitude, so a negative V gives a positive fraction
    assert np.allclose(frac, [0.3, 0.05])
    # error never divides by V: sqrt((e_v/I)^2 + (frac*e_i/I)^2)
    expect0 = np.sqrt((0.2 / 10) ** 2 + (0.3 * 0.5 / 10) ** 2)
    assert np.isclose(err[0], expect0)


def test_handedness():
    assert stokesv.handedness(2.0) == "RCP"
    assert stokesv.handedness(-2.0) == "LCP"


def test_leakage_floor():
    frac = np.array([0.01, 0.02, 0.03, np.nan])  # median (ignoring nan) = 0.02
    assert np.isclose(stokesv.leakage_floor(frac, n_sigma=7.0), 0.14)
    assert np.isnan(stokesv.leakage_floor(np.array([np.nan])))


def test_select_circular_pol():
    i_flux = np.array([10.0, 10.0, 10.0])
    v_flux = np.array([3.0, 0.05, 4.0])  # frac = 0.3, 0.005, 0.4
    e_i = np.array([0.5, 0.5, 0.5])
    e_v = np.array([0.2, 0.2, 2.0])  # last has low V SNR (4/2 = 2)
    mask, frac = stokesv.select_circular_pol(
        i_flux, v_flux, e_i, e_v, leakage_threshold=0.05, v_snr_min=5.0
    )
    assert mask.tolist() == [True, False, False]  # #1 passes; #2 below floor; #3 fails SNR
    assert np.isclose(frac[0], 0.3)


def test_proper_motion_confirm():
    # a high-PM star: 200 mas/yr over 20 yr = 4" shift. Radio at the propagated position confirms.
    ra, dec = np.array([45.0]), np.array([-30.0])
    pmra, pmdec = np.array([200.0]), np.array([0.0])
    dt = 20.0
    cosd = np.cos(np.radians(-30.0))
    ra_prop = 45.0 + (200.0 * dt / 1000 / 3600) / cosd
    ok, sep = stokesv.proper_motion_confirm(
        np.array([ra_prop]), dec, ra, dec, pmra, pmdec, dt, match_arcsec=2.5
    )
    assert ok[0] and sep[0] < 0.1
    # radio at the (static) catalogue position fails — it does not track the proper motion
    ok2, sep2 = stokesv.proper_motion_confirm(ra, dec, ra, dec, pmra, pmdec, dt, match_arcsec=2.5)
    assert not ok2[0] and sep2[0] > 3.0


def test_classify_emitter():
    assert stokesv.classify_emitter(4.0, 10.0) == "highly_circular"  # 0.4
    assert stokesv.classify_emitter(1.0, 10.0) == "circular"  # 0.1
    assert stokesv.classify_emitter(0.1, 10.0) == "weak"  # 0.01
    assert stokesv.classify_emitter(np.nan, 10.0) == "nan"


def test_synthetic_field_shapes():
    stars, dt = stokesv.synthetic_field(n_stars=300, pol_fraction=0.1, seed=2)
    assert dt > 0
    assert stars["is_emitter"].sum() >= 1
    for k in ("ra", "dec", "pmra", "pmdec", "ra_radio", "dec_radio", "i_flux", "v_flux"):
        assert stars[k].shape == (300,)
    # injected emitters have deep circular polarization; the leakage population does not
    frac = np.abs(stars["v_flux"]) / stars["i_flux"]
    assert frac[stars["is_emitter"]].mean() > 0.15
    assert np.median(frac[~stars["is_emitter"]]) < 0.05


def test_match_targets_to_radio():
    # target 0 has a RACS-I component 3" away; target 1 has none within 15"
    t_ra = np.array([45.0, 50.0])
    t_dec = np.array([-30.0, -30.0])
    r_ra = np.array([45.0 + 3.0 / 3600.0 / np.cos(np.radians(-30.0)), 10.0])
    r_dec = np.array([-30.0, 10.0])
    r_i = np.array([12.0, 99.0])
    r_ei = np.array([0.5, 0.5])
    out = match_targets_to_radio(t_ra, t_dec, r_ra, r_dec, r_i, r_ei, radius_arcsec=15.0)
    assert out["matched"].tolist() == [True, False]
    assert np.isclose(out["i_flux"][0], 12.0) and np.isnan(out["i_flux"][1])
    assert out["sep_arcsec"][0] < 4.0
    # empty radio catalogue -> all unmatched, no crash
    empty = match_targets_to_radio(
        t_ra, t_dec, np.array([]), np.array([]), np.array([]), np.array([])
    )
    assert not empty["matched"].any()


def test_measure_circular_pol():
    from astropy.wcs import WCS

    # a 2.5"/pixel SIN-projection image centred on the target
    ra0, dec0 = 45.0, -30.0
    w = WCS(naxis=2)
    w.wcs.crpix = [50.0, 50.0]
    w.wcs.cdelt = [-2.5 / 3600.0, 2.5 / 3600.0]
    w.wcs.crval = [ra0, dec0]
    w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
    ny = nx = 100
    yy, xx = np.mgrid[0:ny, 0:nx]
    sig = 2.0  # beam sigma in pixels
    gauss = np.exp(-((xx - 49.0) ** 2 + (yy - 49.0) ** 2) / (2 * sig**2))
    image_i = 20.0 * gauss  # 20 mJy point source at the target
    image_v = -0.4 * image_i  # 40% circularly polarized, LCP (negative V)
    rng = np.random.default_rng(0)
    image_i = image_i + rng.normal(0, 0.05, (ny, nx))
    image_v = image_v + rng.normal(0, 0.05, (ny, nx))
    m = measure_circular_pol(image_i, image_v, w, ra0, dec0, search_arcsec=12.0)
    assert abs(m["i_peak"] - 20.0) < 0.5
    assert m["v_peak"] < 0  # recovers the LCP sign
    assert abs(m["frac_pol"] - 0.4) < 0.05  # recovers injected |V|/I
    assert m["offset_arcsec"] < 3.0
    # a target off the image returns all-NaN, no crash
    off = measure_circular_pol(image_i, image_v, w, ra0 + 5.0, dec0, search_arcsec=12.0)
    assert np.isnan(off["frac_pol"])


def test_run_offline(tmp_path):
    m = stokesv.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_injected"] >= 1
    assert m["n_recovered"] / m["n_injected"] > 0.8  # recovers most injected emitters
    assert m["purity"] > 0.8  # candidates are mostly genuine (leakage + PM cuts work)
    # the floor is estimated from bright sources, so it sits at a physical few-% leakage level,
    # not the tens-of-% a noise-dominated faint-source median would give
    assert 0.0 < m["leakage_floor_pct"] < 12.0
    assert m["n_bright_ref"] > 0
    assert (tmp_path / "results" / "stokesv_metrics.json").exists()
    assert (tmp_path / "papers" / "stokesv" / "figures" / "circular_pol.pdf").exists()
    macros = (tmp_path / "papers" / "stokesv" / "generated" / "macros.tex").read_text()
    assert r"\svNcandidates" in macros
    assert r"\svLeakFloorPct" in macros
    # the real forced-photometry macros are present (placeholders offline) so the paper compiles
    assert r"\svIrec" in macros and r"\svFracVcirc" in macros


def test_racs_science_mask_excludes_noisemap():
    from astropy.table import Table

    from jansky_research.stokesv import _racs_science_mask

    t = Table(
        {
            "filename": [
                "image.v.RACS_0141-18.SB1.cont.taylor.0.restored.conv.fits",  # science V
                "noiseMap.image.v.RACS_0141-18.SB1.cont.taylor.0.restored.conv.fits",  # noise map
                "meanMap.image.v.RACS_0141-18.SB1.cont.taylor.0.restored.conv.fits",  # mean map
                "image.i.RACS_0141-18.SB1.cont.taylor.0.restored.conv.fits",  # science I (wrong stokes)
            ]
        }
    )
    assert _racs_science_mask(t, "v").tolist() == [True, False, False, False]


def test_run_real_merges_forced_photometry(monkeypatch, tmp_path):
    # mock the CASDA forced-photometry (network) with representative rows: I recovered, V variability-limited
    rows = [
        {
            "cat_i": 19.2,
            "cat_frac": 0.90,
            "img_i": 10.5,
            "img_v": 0.35,
            "img_frac": 0.03,
            "offset_arcsec": 8,
        },
        {
            "cat_i": 16.7,
            "cat_frac": 0.59,
            "img_i": 9.1,
            "img_v": 4.2,
            "img_frac": 0.46,
            "offset_arcsec": 5,
        },
        {
            "cat_i": 14.2,
            "cat_frac": 0.42,
            "img_i": 11.0,
            "img_v": 0.2,
            "img_frac": 0.02,
            "offset_arcsec": 6,
        },
        {
            "cat_i": 12.4,
            "cat_frac": 0.77,
            "img_i": 7.0,
            "img_v": 2.4,
            "img_frac": 0.34,
            "offset_arcsec": 4,
        },
    ]
    monkeypatch.setattr(stokesv, "forced_photometry_recover", lambda **k: rows)
    m = stokesv.run(out=str(tmp_path), offline=False)
    # merged metrics carry BOTH the synthetic-validation and the real forced-photometry results
    assert m["source"] == "RACS-low DR1 (CASDA)"
    assert m["purity"] > 0.8  # synthetic selection machinery still ran
    assert m["n_measured"] == 4
    assert 0.3 < m["i_recovery_ratio"] < 1.0  # Stokes I recovered at the known positions
    assert m["n_v_circular"] == 2 and m["frac_v_circular"] == 0.5  # variability-limited V
    macros = (tmp_path / "papers" / "stokesv" / "generated" / "macros.tex").read_text()
    assert r"\svNmeasured}{4}" in macros and r"\svIrec" in macros
