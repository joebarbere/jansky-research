"""Tests for jansky_research.junodam -- DAM occurrence census. Offline."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import junodam as jdm


def test_io_mean_longitude_rate_and_range():
    jd = jdm.J2000_JD + np.array([0.0, 1.0])
    lon = jdm.io_mean_longitude(jd)
    assert 0 <= lon[0] < 360 and abs(lon[0] - 106.07719) < 1e-6
    # one day advances by the Lieske rate mod 360
    assert abs(((lon[1] - lon[0]) % 360.0) - (203.488955790 % 360.0)) < 1e-6


def test_io_phase_wraps():
    jd = np.array([jdm.J2000_JD])
    assert 0 <= jdm.io_phase(jd, np.array([350.0]))[0] < 360


def test_detect_active_threshold():
    af = np.array([0.0, 0.05, 0.1, 0.5])
    assert jdm.detect_active(af).tolist() == [False, False, True, True]


def test_occurrence_map_masks_low_exposure():
    cml = np.array([10.0] * 5 + [200.0])
    pha = np.array([10.0] * 5 + [200.0])
    act = np.array([True, True, False, False, False, True])
    m = jdm.occurrence_map(cml, pha, act, n_bins=18, min_exposure=3)
    assert np.isfinite(m["occ"][0, 0]) and abs(m["occ"][0, 0] - 0.4) < 1e-9
    assert np.isnan(m["occ"][10, 10])  # single visit -> masked


def test_synthetic_orbit_recovers_io_contrast():
    s = jdm.synthetic_orbit(seed=1)
    m = jdm.occurrence_map(s["cml"], s["io_phase"], s["active"])
    con = jdm.io_region_contrast(m)
    expected = s["p_in"] / s["p_out"]
    assert con["contrast"] > 4.0  # strong recovered Io-region enhancement
    assert con["contrast"] < 1.5 * expected  # and not inflated
    assert con["cells_used"] > 200  # a month covers most of the plane


def test_run_offline_writes_artifacts(tmp_path):
    m = jdm.run(str(tmp_path), offline=True)
    assert m["io_contrast"] and m["io_contrast"] > 4.0
    saved = json.loads((tmp_path / "results" / "junodam_metrics.json").read_text())
    assert saved == m
    assert (tmp_path / "papers" / "junodam" / "figures" / "junodam.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "junodam" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\jdContrast}" in macros


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    jdm._write_macros({"source": "x", "io_contrast": None}, p)
    assert r"\newcommand{\jdContrast}{--}" in p.read_text()
