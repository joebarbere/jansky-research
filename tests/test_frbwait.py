"""Tests for jansky_research.frbwait -- the uniform Cat-2 repeater census. Offline."""

from __future__ import annotations

import json

import numpy as np
import pytest

from jansky_research import frbwait as fw


def test_sidereal_scramble_preserves_phase_and_span():
    rng = np.random.default_rng(0)
    mjd = np.sort(59000.0 + rng.uniform(0, 400, 60))
    scr = fw.sidereal_scramble(mjd, rng)
    assert scr.size == mjd.size
    # sidereal phases are preserved as a multiset
    ph = np.sort(mjd % fw.SIDEREAL_DAY)
    ph_s = np.sort(scr % fw.SIDEREAL_DAY)
    assert np.allclose(ph, ph_s, atol=1e-9)
    # scrambled days stay within the original day span
    assert scr.min() >= mjd.min() - fw.SIDEREAL_DAY
    assert scr.max() <= mjd.max() + fw.SIDEREAL_DAY


def test_periodogram_finds_injected_period_and_scramble_kills_it():
    syn = fw.synthetic_repeater_set(period=16.35, duty=0.3, seed=2)
    pg = fw.scramble_fap_periodogram(syn["mjd"], n_scramble=60, seed=2)
    # the recovered peak is the injected period (or a 1-cycle alias of the grid resolution)
    assert abs(pg["best_period"] - 16.35) < 0.35
    assert pg["p_scramble"] < 0.05
    # a Poisson control with no activity window must not be significant
    ctrl = fw.synthetic_repeater_set(k=1.0, duty=1.0, seed=3)
    pg0 = fw.scramble_fap_periodogram(ctrl["mjd"], n_scramble=60, seed=3)
    assert pg0["p_scramble"] > 0.05


def test_activity_window_recovers_duty_cycle():
    syn = fw.synthetic_repeater_set(period=16.35, duty=0.3, seed=4)
    win = fw.activity_window(syn["mjd"], 16.35, containment=1.0)
    # full containment window ~ injected duty (transit quantisation broadens it slightly)
    assert 0.2 < win["duty_cycle"] < 0.45
    assert win["n_in_arc"] == syn["mjd"].size


def test_census_recovers_k_and_flags_poisson_control():
    trains = {
        "SYN-CLUSTERED": fw._train_from_mjd(fw.synthetic_repeater_set(k=0.4, seed=5)["mjd"]),
        "SYN-POISSON": fw._train_from_mjd(
            fw.synthetic_repeater_set(k=1.0, duty=1.0, seed=6)["mjd"]
        ),
    }
    rows = fw.census(trains, n_scramble=40, seed=5)
    by = {r["name"]: r for r in rows}
    assert by["SYN-CLUSTERED"]["weibull_k"] < 0.8
    assert by["SYN-CLUSTERED"]["clustered"]
    # the Poisson control's k CI must include 1 (not clustered)
    assert not by["SYN-POISSON"]["clustered"]


def test_census_respects_completeness_cut():
    few = fw._train_from_mjd(59000.0 + np.array([1.0, 30.0, 200.0, 300.0, 500.0]))
    rows = fw.census({"SYN-FEW": few}, min_bursts_stats=10, n_scramble=10)
    assert rows[0]["n_bursts"] == 5
    assert "weibull_k" not in rows[0]  # below the cut: rates only, honestly
    assert np.isfinite(rows[0]["rate_per_hr"])


def test_load_catalog2_dedupes_and_masks(tmp_path):
    csv_text = (
        "tns_name,previous_name,repeater_name,event_id,sub_num,ra,ra_err,dec,dec_err,"
        "ra_dec_notes,gl,gb,exp_up,exp_up_err,exp_low,exp_low_err,exp_notes,bonsai_snr,"
        "bonsai_dm,low_ft_68,up_ft_68,low_ft_95,up_ft_95,snr_fitb,dm_fitb,dm_fitb_err,"
        "dm_exc_ne2001,dm_exc_ymw16,bc_width,scat_time,scat_time_err,flux,flux_err,fluence,"
        "fluence_err,fluence_notes,fluence_win_extended,mjd_400,mjd_400_err,mjd_inf,"
        "mjd_inf_err,width_fitb,width_fitb_err,sp_idx,sp_idx_err,sp_run,sp_run_err,high_freq,"
        "low_freq,peak_freq,chi_sq,dof,flag_frac,notes_fitb,intrachan_flag,excluded_flag,"
        "sidelobe_flag,citizen_science_flag,catalog1_flag,catalog1_param_flag\n"
        "FRB1,,FRBR1,ev1,0,10,0.1,20,0.1,,100,50,30,1,20,1,,12,500,0,0,0,0,15,500.5,0.2,"
        "400,390,5,-9999,0,1,0.1,5,1,,0,59000.5,1e-05,59000.5,1e-05,2,0.1,0,0,0,0,800,400,"
        "600,1,1,0,,0,0,0,0,0,0\n"
        "FRB1,,FRBR1,ev1,1,10,0.1,20,0.1,,100,50,30,1,20,1,,12,500,0,0,0,0,15,500.6,0.2,"
        "400,390,5,-9999,0,1,0.1,5,1,,0,59000.6,1e-05,59000.6,1e-05,2,0.1,0,0,0,0,800,400,"
        "600,1,1,0,,0,0,0,0,0,0\n"
        "FRB2,,-9999,ev2,0,11,0.1,21,0.1,,101,51,31,1,21,1,,13,600,0,0,0,0,16,600.5,0.2,"
        "500,490,5,0.001,0,1,0.1,6,1,,0,59001.5,1e-05,59001.5,1e-05,2,0.1,0,0,0,0,800,400,"
        "600,1,1,0,,0,0,0,0,0,0\n"
    )
    p = tmp_path / "cat2.csv"
    p.write_text(csv_text)
    cat = fw.load_catalog2(p)
    assert cat["mjd"].size == 2  # ev1 sub-bursts deduped
    assert np.isnan(cat["scat_time"][0])  # -9999 -> NaN
    assert cat["repeater_name"][1] == "-9999"


def test_repeater_trains_groups_and_filters():
    cat = {
        "repeater_name": np.array(["A", "A", "A", "-9999", "B"]),
        "mjd": np.array([59002.0, 59001.0, 59003.0, 59000.0, 59005.0]),
        "dm": np.ones(5),
        "dm_err": np.ones(5),
        "fluence": np.ones(5),
        "width": np.ones(5),
        "exp_up_hr": np.full(5, 40.0),
    }
    trains = fw.repeater_trains(cat, min_bursts=3)
    assert list(trains) == ["A"]
    assert np.all(np.diff(trains["A"]["mjd"]) > 0)  # sorted


def test_run_offline_writes_artifacts(tmp_path):
    m = fw.run(str(tmp_path), offline=True, n_scramble=40)
    assert m["source"].startswith("synthetic")
    by = {r["name"]: r for r in m["rows"]}
    inj = by["SYN-INJECTED"]
    assert abs(inj["best_period"] - m["true_period"]) < 0.35
    assert inj["p_scramble"] < 0.05
    assert inj["clustered"] and abs(inj["weibull_k"] - m["true_k"]) < 0.25
    assert "duty_cycle" in inj and 0.2 < inj["duty_cycle"] < 0.45
    assert by["SYN-POISSON"]["p_scramble"] > 0.05
    saved = json.loads((tmp_path / "results" / "frbwait_metrics.json").read_text())
    assert saved["n_periodic_p01"] == m["n_periodic_p01"]
    assert (tmp_path / "papers" / "frbwait" / "figures" / "frbwait.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "frbwait" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\fwSynMedianK}" in macros and r"\newcommand{\fwRealMedianK}{--}" in macros
    table = (tmp_path / "papers" / "frbwait" / "generated" / "census_table.tex").read_text()
    assert "SYN-INJECTED" in table


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    fw._write_macros({"source": "x", "is_real": True, "anchor_period": None}, p)
    txt = p.read_text()
    assert r"\newcommand{\fwRealAnchorPeriod}{--}" in txt
    assert r"\newcommand{\fwSynAnchorPeriod}{--}" in txt


def test_default_period_grid_bounds():
    g = fw.default_period_grid(1800.0)
    assert g.min() == pytest.approx(1.5, rel=1e-6)
    assert g.max() == pytest.approx(600.0, rel=1e-6)
    assert np.all(np.diff(g) > 0)
