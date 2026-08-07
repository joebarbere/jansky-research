"""Tests for jansky_research.dr20radio — DR20 BHM radio census. No network."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import dr20radio


def test_select_quasars_filters_and_flags_cartons():
    cls = np.array(["QSO", "GALAXY", "QSO ", "STAR", "QSO"])
    zw = np.array([0, 0, 0, 0, 5])
    carton = np.array(
        ["bhm_spiders_agn", "x", "openfibertargets_bhm_racsradio_boss", "y", "bhm_csc"]
    )
    quasar, radio = dr20radio.select_quasars(cls, zw, carton)
    assert quasar.tolist() == [True, False, True, False, False]  # ZWARNING!=0 excluded
    assert radio.tolist() == [False, False, True, False, False]


def test_crossmatch_known_geometry():
    ra_q = np.array([180.0, 181.0])
    dec_q = np.array([10.0, 10.0])
    # radio source 1" from quasar 0; nothing near quasar 1
    ra_r = np.array([180.0 + 1.0 / 3600.0 / np.cos(np.deg2rad(10.0)), 200.0])
    dec_r = np.array([10.0, -5.0])
    m, sep = dr20radio.crossmatch(ra_q, dec_q, ra_r, dec_r, radius_arcsec=2.5)
    assert m.tolist() == [True, False]
    assert sep[0] == pytest.approx(1.0, abs=0.05)


def test_wilson_interval_properties():
    p, lo, hi = dr20radio.wilson_interval(10, 100)
    assert p == pytest.approx(0.1)
    assert lo < p < hi
    assert dr20radio.wilson_interval(0, 0)[0] != dr20radio.wilson_interval(0, 0)[0]  # nan


def test_synthetic_round_trip_recovers_fraction_and_carton_circularity():
    s = dr20radio.synthetic_survey(seed=1)
    quasar, radio_carton = dr20radio.select_quasars(s["cls"], s["zwarning"], s["firstcarton"])
    assert quasar.all()
    census = ~radio_carton
    m, _ = dr20radio.crossmatch(
        s["ra_q"][census],
        s["dec_q"][census],
        s["ra_r"],
        s["dec_r"],
        radius_arcsec=s["radius_arcsec"],
    )
    fm = dr20radio.false_match_rate(
        s["ra_q"][census],
        s["dec_q"][census],
        s["ra_r"],
        s["dec_r"],
        radius_arcsec=s["radius_arcsec"],
        n_trials=5,
        seed=2,
    )
    corrected = float(np.mean(m)) - fm["rate"]
    assert corrected == pytest.approx(s["true_fraction"], abs=0.03)
    # the radio-carton subset matches ~always (counterpart by construction) — the
    # circularity that must stay OUT of the census fraction
    mc, _ = dr20radio.crossmatch(
        s["ra_q"][radio_carton],
        s["dec_q"][radio_carton],
        s["ra_r"],
        s["dec_r"],
        radius_arcsec=s["radius_arcsec"],
    )
    assert float(np.mean(mc)) > 0.88  # 1" Rayleigh jitter puts ~4% beyond 2.5"
    assert float(np.mean(mc)) - float(np.mean(m)) > 0.5  # excluding them matters


def test_false_match_rate_scales_with_density():
    s = dr20radio.synthetic_survey(seed=3, n_radio=2000, counterpart_fraction=0.0)
    fm_lo = dr20radio.false_match_rate(
        s["ra_q"], s["dec_q"], s["ra_r"], s["dec_r"], n_trials=4, seed=4, radius_arcsec=30.0
    )
    s2 = dr20radio.synthetic_survey(seed=3, n_radio=20000, counterpart_fraction=0.0)
    fm_hi = dr20radio.false_match_rate(
        s2["ra_q"], s2["dec_q"], s2["ra_r"], s2["dec_r"], n_trials=4, seed=4, radius_arcsec=30.0
    )  # 30" radius gives enough chance matches for the scaling comparison to have power
    assert fm_hi["rate"] > fm_lo["rate"]  # denser radio catalog -> more chance matches
    assert fm_hi["rate"] == pytest.approx(10 * fm_lo["rate"], rel=0.6)


def test_detection_fraction_binning():
    z = np.array([0.5, 0.6, 1.5, 1.6, 3.0])
    matched = np.array([True, False, True, True, False])
    out = dr20radio.detection_fraction(z, matched, bins=np.array([0.0, 1.0, 2.0, 4.0]))
    assert out["n"] == [2, 2, 1]
    assert out["k"] == [1, 2, 0]
    assert out["frac"][1] == pytest.approx(1.0)


def test_read_spall_quasars_on_synthetic_fits(tmp_path):
    from astropy.io import fits

    n = 6
    cols = [
        fits.Column(name="RACAT", format="D", array=np.linspace(10, 15, n)),
        fits.Column(name="DECCAT", format="D", array=np.linspace(-50, 50, n)),
        fits.Column(name="Z", format="D", array=np.linspace(0.5, 3.0, n)),
        fits.Column(name="ZWARNING", format="J", array=np.array([0, 0, 0, 4, 0, 0])),
        fits.Column(
            name="CLASS",
            format="10A",
            array=np.array(["QSO", "QSO", "GALAXY", "QSO", "QSO", "QSO"]),
        ),
        fits.Column(
            name="FIRSTCARTON",
            format="40A",
            array=np.array(
                ["bhm_spiders", "openfibertargets_bhm_racsradio_boss", "x", "y", "z", "w"]
            ),
        ),
        fits.Column(
            name="OBS", format="3A", array=np.array(["APO", "LCO", "APO", "APO", "LCO", "APO"])
        ),
    ]
    path = tmp_path / "spall.fits"
    fits.BinTableHDU.from_columns(cols).writeto(path)
    q = dr20radio.read_spall_quasars(str(path))
    assert q["n_total_rows"] == 6
    assert q["ra"].size == 4  # 6 - GALAXY - ZWARNING!=0
    assert q["radio_carton"].sum() == 1
    assert set(q["obs"]) <= {"APO", "LCO"}


def test_two_survey_synthetic_fixes_the_carton_blind_spot():
    s = dr20radio.synthetic_two_surveys(seed=5, fade_fraction=0.35)
    _, radio_carton = dr20radio.select_quasars(s["cls"], s["zwarning"], s["firstcarton"])
    # vs the SELECTING survey: ~100% (counterpart by construction)
    m_sel, _ = dr20radio.crossmatch(
        s["ra_q"][radio_carton],
        s["dec_q"][radio_carton],
        s["ra_r"],
        s["dec_r"],
        radius_arcsec=s["radius_arcsec"],
    )
    assert float(np.mean(m_sel)) > 0.88
    # vs the other-frequency survey: ~fade_fraction (the increment-1 blind spot, now modeled)
    m_oth, _ = dr20radio.crossmatch(
        s["ra_q"][radio_carton],
        s["dec_q"][radio_carton],
        s["ra_r2"],
        s["dec_r2"],
        radius_arcsec=s["radius_arcsec"],
    )
    assert float(np.mean(m_oth)) == pytest.approx(s["fade_fraction"], abs=0.12)
    assert float(np.mean(m_sel)) - float(np.mean(m_oth)) > 0.3


def test_parse_racs_csv():
    text = "ra,dec,peak_flux\n3.14,-41.5,2.5\nbad,row\n10.0,-50.0,1.1\n"
    ra, dec, flux = dr20radio.parse_racs_csv(text)
    assert ra.tolist() == [3.14, 10.0]
    assert dec.tolist() == [-41.5, -50.0]
    assert flux.tolist() == [2.5, 1.1]
