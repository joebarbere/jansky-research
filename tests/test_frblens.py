"""Tests for jansky_research.frblens -- the lensed-repeater delay search. Offline."""

from __future__ import annotations

import json

import numpy as np
import pytest

from jansky_research import frblens as fl
from jansky_research.frbwait import SIDEREAL_DAY, synthetic_repeater_set


def test_all_pairs_delays_excludes_subday():
    t = 59000.0 + np.array([0.0, 0.01, 2.0, 5.0])
    d = fl.all_pairs_delays(t)
    assert d.min() >= fl.MIN_DELAY_DAYS
    assert d.size == 5  # 6 pairs minus the 0.01-day one


def test_match_count_finds_disjoint_lensed_pairs():
    # 3 intrinsic bursts, each with an image exactly 2 sidereal days later
    base = 59000.0 + np.array([0.0, 7.0 * SIDEREAL_DAY, 20.0 * SIDEREAL_DAY])
    delay = 2.0 * SIDEREAL_DAY
    t = np.sort(np.concatenate([base, base + delay]))
    m = fl.match_count_at_delay(t, delay, tol_days=300.0 / 86400.0)
    assert m == 3


def test_match_count_respects_dm_tolerance():
    t = 59000.0 + np.array([0.0, 2.0])
    dm = np.array([500.0, 510.0])
    assert fl.match_count_at_delay(t, 2.0, tol_days=0.01, dm=dm, dm_tol=1.0) == 0
    assert fl.match_count_at_delay(t, 2.0, tol_days=0.01, dm=dm, dm_tol=20.0) == 1


def test_recurring_delay_search_recovers_injection():
    base = synthetic_repeater_set(k=0.7, duty=1.0, mean_wait=12.0, seed=0)
    t0 = base["mjd"]
    inj = fl.inject_lensed_train(
        t0, np.full(t0.size, 10.0), delay=3.0 * SIDEREAL_DAY, mag_ratio=1.0
    )
    assert inj["detectable"] and inj["n_images_detected"] >= 2
    found = fl.recurring_delay_search(inj["mjd"])
    assert abs(found["best_delay"] - 3.0 * SIDEREAL_DAY) < 300.0 / 86400.0
    assert found["m_max"] >= inj["n_images_detected"]


def test_scramble_fap_separates_lensed_from_control():
    base = synthetic_repeater_set(k=0.7, duty=1.0, mean_wait=12.0, seed=0)
    t0 = base["mjd"]
    inj = fl.inject_lensed_train(
        t0, np.full(t0.size, 10.0), delay=3.0 * SIDEREAL_DAY, mag_ratio=1.0
    )
    found = fl.recurring_delay_search(inj["mjd"])
    fap = fl.scramble_fap(inj["mjd"], found["m_max"], n_scramble=40)
    assert fap["p_value"] < 0.05
    found0 = fl.recurring_delay_search(t0)
    fap0 = fl.scramble_fap(t0, found0["m_max"], n_scramble=40, seed=1)
    assert fap0["p_value"] > 0.05


def test_phase_scramble_preserves_days_and_intrinsic_clustering_is_null():
    # the false-positive mode caught on first real-catalogue contact: an intrinsically
    # CLUSTERED train (consecutive-day activity epoch, random phases in the transit window)
    # must be null under the phase scramble even though a day scramble would flag it
    rng = np.random.default_rng(7)
    days = np.repeat(np.arange(0, 40), 3)  # a dense 40-day epoch, 3 bursts/day
    phases = 0.3 + rng.uniform(0.0, 15.0 / (24 * 60), days.size)
    t = np.sort(days * SIDEREAL_DAY + phases + 59000.0)
    scr = fl.phase_scramble(t, rng)
    assert np.array_equal(
        np.floor(np.sort(t) / SIDEREAL_DAY), np.floor(scr / SIDEREAL_DAY)
    )  # day assignments preserved
    found = fl.recurring_delay_search(t)
    fap = fl.scramble_fap(t, found["m_max"], n_scramble=40, seed=8)
    assert fap["p_value"] > 0.05  # clustering alone cannot beat the phase-permutation null


def test_inject_transit_geometry_kills_offgrid_delays():
    base = synthetic_repeater_set(k=0.7, duty=1.0, mean_wait=12.0, seed=0)
    t0 = base["mjd"]
    f0 = np.full(t0.size, 10.0)
    # half a sidereal day off: the trailing images transit outside the beam window
    inj = fl.inject_lensed_train(t0, f0, delay=2.5 * SIDEREAL_DAY, mag_ratio=1.0)
    assert inj["n_images_detected"] == 0
    # too faint: magnification ratio below the empirical fluence floor
    faint = fl.inject_lensed_train(t0, f0, delay=3.0 * SIDEREAL_DAY, mag_ratio=0.01)
    assert faint["n_images_detected"] == 0


def test_sensitivity_map_lights_up_ongrid_cells_only():
    base = synthetic_repeater_set(k=0.7, duty=1.0, mean_wait=12.0, seed=0)
    t0 = base["mjd"]
    f0 = np.full(t0.size, 10.0)
    sens = fl.sensitivity_map(
        t0,
        f0,
        delays=np.array([2.0 * SIDEREAL_DAY, 45.7]),
        mag_ratios=np.array([1.0]),
        n_scramble=20,
    )
    assert sens["sensitive"][0, 0]  # on-comb, equal magnification: detected
    assert not sens["sensitive"][1, 0]  # off-comb: images never transit -> dark cell


def test_barycentric_offsets_have_annual_roemer_amplitude():
    fn = fl.barycentric_offset_fn(50.0, 20.0)
    mjd = 59000.0 + np.arange(0.0, 365.0, 5.0)
    off_s = fn(mjd) * 86400.0
    assert np.max(np.abs(off_s)) < 600.0  # bounded by the light travel time to 1 au
    assert np.ptp(off_s) > 300.0  # annual Roemer swing for a low-ecliptic-latitude source


def test_barycentric_delay_fixed_in_bary_frame_drifts_in_topo():
    # the GATE-2 catch: a barycentre-fixed delay is NOT fixed in topocentric TOAs
    base = synthetic_repeater_set(k=0.7, duty=1.0, mean_wait=12.0, seed=0)
    t0 = base["mjd"]
    f0 = np.full(t0.size, 10.0)
    fn = fl.barycentric_offset_fn(50.0, 20.0)
    delay = 20.0 * SIDEREAL_DAY
    inj = fl.inject_lensed_train(t0, f0, delay=delay, mag_ratio=1.0, bary_fn=fn)
    assert inj["n_images_detected"] >= 5  # Roemer wobble stays inside the transit window
    # barycentric search recovers the fixed delay ...
    tb = inj["mjd"] + inj["bary_offsets"]
    found_b = fl.recurring_delay_search(tb)
    assert abs(found_b["best_delay"] - delay) < fl.TOL_DAYS
    assert found_b["m_max"] >= inj["n_images_detected"] - 1
    fap = fl.scramble_fap(
        inj["mjd"], found_b["m_max"], bary_offsets=inj["bary_offsets"], n_scramble=40
    )
    assert fap["p_value"] < 0.05
    # ... while the same search on raw topocentric times misses most pairs at 5-s tolerance
    found_t = fl.recurring_delay_search(inj["mjd"])
    assert found_t["m_max"] < found_b["m_max"]


def test_injection_dm_scatter_costs_the_dm_cut():
    base = synthetic_repeater_set(k=0.7, duty=1.0, mean_wait=12.0, seed=0)
    t0 = base["mjd"]
    f0 = np.full(t0.size, 10.0)
    dm0 = np.full(t0.size, 500.0)
    inj = fl.inject_lensed_train(
        t0, f0, delay=3.0 * SIDEREAL_DAY, mag_ratio=1.0, dm=dm0, dm_sigma=3.0, seed=3
    )
    assert inj["dm"].size == inj["mjd"].size
    # a too-tight DM cut (the plan's 1 pc/cc) rejects most injected image pairs
    tight = fl.recurring_delay_search(inj["mjd"], dm=inj["dm"], dm_tol=1.0)
    wide = fl.recurring_delay_search(inj["mjd"], dm=inj["dm"], dm_tol=3.0 * np.sqrt(2.0) * 3.0)
    assert wide["m_max"] > tight["m_max"]
    assert wide["m_max"] >= inj["n_images_detected"] - 1


def test_lensed_fraction_limit_values():
    assert fl.lensed_fraction_limit(10, 0) == pytest.approx(-np.log(0.05) / 10)
    assert fl.lensed_fraction_limit(1, 0) == pytest.approx(2.9957, abs=1e-3)
    # 1 detection raises the limit above the 0-detection value
    assert fl.lensed_fraction_limit(10, 1) > fl.lensed_fraction_limit(10, 0)


def test_run_offline_writes_artifacts(tmp_path):
    m = fl.run(str(tmp_path), offline=True, n_scramble=40)
    by = {r["name"]: r for r in m["rows"]}
    assert by["SYN-LENSED"]["p_value"] < 0.05 < by["SYN-CONTROL"]["p_value"]
    assert m["recovered_delay_err_s"] < 300.0
    assert m["n_detections"] == 1 and m["n_searched"] == 1
    saved = json.loads((tmp_path / "results" / "frblens_metrics.json").read_text())
    assert saved["n_detections"] == 1
    assert (tmp_path / "papers" / "frblens" / "figures" / "frblens.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "frblens" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\flSynNDet}" in macros and r"\newcommand{\flRealNDet}{--}" in macros


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    fl._write_macros({"source": "x", "is_real": True, "n_searched": None}, p)
    txt = p.read_text()
    assert r"\newcommand{\flRealNSearched}{--}" in txt
    assert r"\newcommand{\flSynNSearched}{--}" in txt
