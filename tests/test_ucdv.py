"""Tests for jansky_research.ucdv -- the plan-59 controls and census helpers. Offline."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import ucdv
from jansky_research.stokesv import measure_circular_pol


def test_assign_set_is_balanced_and_data_blind():
    # Gaia DR2 ids in the Reyle table are all even, so parity cannot split them; the hash can
    even_ids = np.arange(10**15, 10**15 + 2 * 4000, 2, dtype=np.int64)
    s = ucdv.assign_set(even_ids)
    assert set(s) == {"A", "B"}
    assert abs(np.mean(s == "B") - 0.5) < 0.03
    assert list(ucdv.assign_set(even_ids[:5])) == list(s[:5])  # deterministic


def test_null_positions_geometry():
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    ra, dec = np.array([10.0, 200.0]), np.array([-30.0, 45.0])
    rn, dn = ucdv.null_positions(ra, dec)
    assert rn.shape == (2 * ucdv.NULL_PER_TARGET,)
    t = SkyCoord(np.repeat(ra, 5) * u.deg, np.repeat(dec, 5) * u.deg)
    sep = t.separation(SkyCoord(rn * u.deg, dn * u.deg)).arcmin
    assert np.allclose(sep, ucdv.NULL_OFFSET_ARCMIN, atol=1e-6)
    # neighbouring nulls of one target are ~141" apart: several 15" beams, i.e. independent
    first = SkyCoord(rn[:5] * u.deg, dn[:5] * u.deg)
    gaps = first[:-1].separation(first[1:]).arcsec
    assert np.all(gaps > 130.0)
    # deterministic under the frozen seed
    rn2, _ = ucdv.null_positions(ra, dec)
    assert np.array_equal(rn, rn2)


def test_inject_point_source_peak_and_orientation():
    img_i, img_v, w, beam = ucdv.synthetic_cutout(rms_mjy=1e-9)
    n = img_v.shape[0]
    ra_c, dec_c = (float(x) for x in w.pixel_to_world_values(n / 2.0, n / 2.0))
    out = ucdv.inject_point_source(img_v, w, ra_c, dec_c, 7.0, ucdv.Beam(30.0, 10.0, 0.0))
    assert out.max() == pytest.approx(7.0, rel=1e-3)
    assert np.allclose(img_v, ucdv.synthetic_cutout(rms_mjy=1e-9)[1])  # input not mutated
    # PA 0: the major axis runs north-south, so the profile is wider along y than x
    iy, ix = np.unravel_index(np.argmax(out), out.shape)
    assert out[iy + 4, ix] > out[iy, ix + 4]
    rot = ucdv.inject_point_source(img_v, w, ra_c, dec_c, 7.0, ucdv.Beam(30.0, 10.0, 90.0))
    assert rot[iy, ix + 4] > rot[iy + 4, ix]


def test_forced_measure_is_genuinely_forced():
    img_i, img_v, w, _ = ucdv.synthetic_cutout()
    n = img_i.shape[0]
    ra_c, dec_c = (float(x) for x in w.pixel_to_world_values(n / 2.0, n / 2.0))
    m = ucdv.forced_measure(img_i, img_v, w, ra_c, dec_c)
    assert m["forced"] is True
    assert m["offset_arcsec"] < 2.5  # within one 2.5" pixel of the position


def test_null_check_passes_gaussian_fails_searching_and_bad_sigma():
    rng = np.random.default_rng(0)
    z = rng.normal(size=400)
    assert ucdv.null_check(z, np.ones(400))["passed"]
    # a searching measurement is positive essentially always -> sign test fails
    pos = ucdv.null_check(np.abs(z), np.ones(400))
    assert not pos["sign_ok"] and not pos["passed"]
    # an RMS underestimated by 1.5x inflates std(V/sigma) out of range
    assert not ucdv.null_check(z, np.full(400, 1 / 1.5))["std_ok"]
    # a planted 5-sigma excess fails the excess test
    z2 = z.copy()
    z2[:6] = 8.0
    assert not ucdv.null_check(z2, np.ones(400))["excess_ok"]
    assert ucdv.null_check(np.array([np.nan]), np.array([1.0])) == {"n": 0, "passed": False}


def test_search_mode_fails_the_null_forced_mode_passes():
    """The stokesv_discovery lesson, as a regression test on real photometry code."""
    v_s, s_s, v_f, s_f = [], [], [], []
    for k in range(40):
        img_i, img_v, w, _ = ucdv.synthetic_cutout(seed=500 + k)
        n = img_i.shape[0]
        ra_c, dec_c = (float(x) for x in w.pixel_to_world_values(n / 2.0, n / 2.0))
        for rn, dn in zip(
            *ucdv.null_positions(np.array([ra_c]), np.array([dec_c]), seed=k), strict=True
        ):
            srch = measure_circular_pol(img_i, img_v, w, float(rn), float(dn), search_arcsec=12.0)
            frc = ucdv.forced_measure(img_i, img_v, w, float(rn), float(dn))
            v_s.append(srch["i_peak"])
            s_s.append(srch["i_rms"])
            v_f.append(frc["i_peak"])
            s_f.append(frc["i_rms"])
    assert not ucdv.null_check(np.array(v_s), np.array(s_s))["sign_ok"]
    assert ucdv.null_check(np.array(v_f), np.array(s_f))["sign_ok"]


def test_injection_check_rules():
    lev = np.repeat([5.0, 7.0, 10.0], 10)
    inj = lev.copy()
    meas = inj.copy()
    r = ucdv.injection_check(lev, inj, meas, np.ones_like(lev))
    assert r["passed"] and r["completeness"][7.0] == 1.0
    meas2 = inj.copy()
    meas2[(lev == 7.0)] = 1.0  # nothing recovered at 7 sigma
    assert not ucdv.injection_check(lev, inj, meas2, np.ones_like(lev))["completeness_ok"]
    meas3 = inj * 0.8  # biased low at >= 10 sigma
    assert not ucdv.injection_check(lev, inj, meas3, np.ones_like(lev))["flux_ok"]


def test_recover_known_check_counts_sign_flips():
    r = ucdv.recover_known_check(np.array([1.0, 2.1, -3.0, -1.0]), np.array([1.0, 2.0, -3.0, 1.0]))
    assert r["n"] == 4 and r["n_sign_flip"] == 1
    assert r["passed"]
    assert not ucdv.recover_known_check(np.array([0.5, 1.0]), np.array([1.0, 2.0]))["passed"]


def test_detection_rule_each_criterion():
    v = np.array([6.0, 4.0, 6.0, 6.0, -6.0])
    sv = np.ones(5)
    i = np.array([1.0, 1.0, 600.0, 1.0, 1.0])  # third: I detected, |V|/I = 1% (leakage)
    si = np.ones(5)
    bright = np.array([False, False, False, True, False])
    r = ucdv.apply_detection_rule(v, sv, i, si, 0.05, bright)
    assert list(r["candidate"]) == [True, False, False, False, True]
    assert not r["v_significant"][1]
    assert not r["clears_leakage"][2]
    assert not r["clear_of_bright"][3]


def test_limits():
    assert np.allclose(ucdv.upper_limit(np.array([-1.0, 0.5]), np.array([0.3, 0.3])), [1.9, 1.4])
    # 1 mJy at 10 pc: 4 pi (3.0857e19 cm)^2 1e-26 = 1.1965e14 erg/s/Hz
    assert ucdv.luminosity_limit(1.0, 10.0) == pytest.approx(1.1965e14, rel=1e-3)


def test_two_epoch_compare_scale_and_spectrum():
    # a 20% drop at high S/N: significant on noise alone, not with the 10% scale term
    r = ucdv.two_epoch_compare(np.array([10.0]), np.array([0.1]), np.array([8.0]), np.array([0.1]))
    assert not r["differs"][0]
    # a big change within what |alpha| <= 3 allows over 887->943 MHz (factor 1.20)
    r = ucdv.two_epoch_compare(
        np.array([10.0]), np.array([0.01]), np.array([11.9]), np.array([0.01]), scale=0.0
    )
    assert r["differs"][0] and not r["beyond_spectral"][0]
    # a sign flip is always beyond spectral
    r = ucdv.two_epoch_compare(np.array([5.0]), np.array([0.2]), np.array([-5.0]), np.array([0.2]))
    assert r["differs"][0] and r["beyond_spectral"][0]
    # appearance from nothing (ratio far outside the spectral band)
    r = ucdv.two_epoch_compare(np.array([0.1]), np.array([0.2]), np.array([6.0]), np.array([0.2]))
    assert r["beyond_spectral"][0]


def test_scope_gate():
    d = np.arange(1.0, 101.0)  # 100 targets at 1..100 pc
    full = ucdv.scope_gate(60.0, 0.0, d, budget_hours=10.0)  # 100 x 2 x 60 s = 3.33 h
    assert full["full_sample"] and full["n"] == 100
    cut = ucdv.scope_gate(60.0, 0.5, d, budget_hours=2.0)  # 2 x 60 / 0.5 = 240 s/target -> 30 fit
    assert not cut["full_sample"] and cut["n"] == 30 and cut["d_max_pc"] == 31.0
    assert cut["hours"] <= 2.0
    assert ucdv.scope_gate(1e6, 0.0, d, budget_hours=1.0)["n"] == 0
    with pytest.raises(ValueError):
        ucdv.scope_gate(60.0, 1.0, d)


def test_validate_synthetic_passes_frozen_controls():
    r = ucdv.validate_synthetic(n_cutouts=80)
    assert r["null"]["n"] == 400
    assert r["null"]["passed"]
    assert r["injection"]["passed"]
    assert r["injection"]["completeness"][3.0] < 0.5  # 3 sigma is below a 5-sigma threshold


def test_parse_polygon_and_in_footprint():
    v = ucdv.parse_polygon("POLYGON ICRS 85.0 38.0 74.2 38.0 73.6 45.4 85.7 45.4")
    assert v.shape == (4, 2)
    inside = ucdv.in_footprint(np.array([80.0, 60.0, 80.0]), np.array([41.0, 41.0, 50.0]), v)
    assert list(inside) == [True, False, False]
    # opposite winding gives the same answer
    assert list(ucdv.in_footprint(np.array([80.0]), np.array([41.0]), v[::-1])) == [True]
    # a tile straddling RA = 0, where a planar test in RA would fail
    w = ucdv.parse_polygon("POLYGON ICRS 355.0 -5.0 5.0 -5.0 5.0 5.0 355.0 5.0")
    assert list(ucdv.in_footprint(np.array([0.5, 359.5, 180.0]), np.zeros(3), w)) == [
        True,
        True,
        False,
    ]
    with pytest.raises(ValueError):
        ucdv.parse_polygon("CIRCLE ICRS 1 2 3")


def test_pick_tile_prefers_most_central_containing_tile():
    edge = "POLYGON ICRS 9.0 -1.0 11.0 -1.0 11.0 1.0 9.0 1.0"  # target near its west edge
    centre = "POLYGON ICRS 8.0 -1.0 10.0 -1.0 10.0 1.0 8.0 1.0"  # target near its centre
    away = "POLYGON ICRS 20.0 -1.0 22.0 -1.0 22.0 1.0 20.0 1.0"  # does not contain it
    assert ucdv.pick_tile(9.1, 0.0, [edge, centre, away]) == 1
    assert ucdv.pick_tile(9.1, 0.0, [away]) == -1


def test_stokes_of_reads_the_stokes_axis():
    from astropy.io import fits

    h = fits.Header({"NAXIS": 4, "CTYPE3": "FREQ", "CTYPE4": "STOKES", "CRVAL4": 4.0})
    assert ucdv._stokes_of(h) == "v"
    h["CRVAL4"] = 1.0
    assert ucdv._stokes_of(h) == "i"
    with pytest.raises(RuntimeError):
        ucdv._stokes_of(fits.Header({"NAXIS": 2, "CTYPE1": "RA---SIN"}))


def test_load_cutout_roundtrip_and_controls_on_synthetic(tmp_path):
    cut = []
    tg = []
    for k in range(80):
        i, v, w, b = ucdv.synthetic_cutout(seed=900 + k)
        p = tmp_path / f"c{k}.npz"
        np.savez_compressed(
            p,
            i=i,
            v=v,
            wcs_header=w.to_header_string(),
            beam=np.array([b.bmaj_arcsec, b.bmin_arcsec, b.bpa_deg]),
            mjd=59700.0,
            obs_id="X",
        )
        c = ucdv.load_cutout(p)
        n = i.shape[0]
        ra, dec = (float(x) for x in c["wcs"].pixel_to_world_values(n / 2.0, n / 2.0))
        cut.append(c)
        tg.append((ra, dec, 0.0, 0.0))
    assert cut[0]["beam"] == ucdv.Beam(15.0, 12.0, 30.0)
    r = ucdv.controls_c1_c2(cut, tg)
    assert r["c1_null_v"]["n"] == 400 and r["c1_null_v"]["passed"]
    assert r["c2_injection"]["passed"]
    assert r["n_injections_per_level"] == 480 and r["n_target_pixel_per_level"] == 80


def test_bright_neighbour_guard():
    img_i, _, w, beam = ucdv.synthetic_cutout()
    n = img_i.shape[0]
    ra_c, dec_c = (float(x) for x in w.pixel_to_world_values(n / 2.0, n / 2.0))
    assert not ucdv.bright_neighbour(img_i, w, ra_c, dec_c)
    # a 200 mJy source 2' away trips the 3' guard; 4' away it does not
    for off_arcmin, expect in ((2.0, True), (4.0, False)):
        ra_b = ra_c + off_arcmin / 60.0 / np.cos(np.radians(dec_c))
        bright = ucdv.inject_point_source(img_i, w, ra_b, dec_c, 200.0, beam)
        assert ucdv.bright_neighbour(bright, w, ra_c, dec_c) is expect
