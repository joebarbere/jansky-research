"""Tests for jansky_research.frbperiod — Rayleigh periodogram recovers injected periods. No network."""

from __future__ import annotations

import numpy as np

from jansky_research import frbperiod


def test_rayleigh_peaks_at_true_period():
    t = frbperiod.synthetic_periodic_arrivals(
        period=16.35, n=80, active_frac=0.12, span=400, seed=1
    )
    at_true = frbperiod.rayleigh_z2(t, 16.35)
    off = frbperiod.rayleigh_z2(t, 9.1)  # an unrelated period
    assert at_true > 4 * off
    assert at_true > 20  # strongly concentrated


def test_period_search_recovers_injected_period():
    t = frbperiod.synthetic_periodic_arrivals(
        period=16.35, n=80, active_frac=0.12, span=400, seed=2
    )
    grid = np.linspace(5.0, 40.0, 7000)
    res = frbperiod.period_search(t, grid)
    assert abs(res.best_period - 16.35) < 0.2  # recovered within grid resolution
    assert res.z2.shape == grid.shape
    assert res.best_z2 == res.z2.max()
    assert res.fap < 1e-3  # a clear signal has a tiny false-alarm probability


def test_fap_bounds_and_monotonic():
    assert 0.0 <= frbperiod.false_alarm_prob(50.0, 1000) <= 1.0
    # higher peak -> lower false-alarm probability
    assert frbperiod.false_alarm_prob(40.0, 1000) < frbperiod.false_alarm_prob(10.0, 1000)
    # random arrival times give an unremarkable peak (FAP not vanishingly small)
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0, 400, 25))
    res = frbperiod.period_search(t, np.linspace(5, 40, 3000))
    assert res.fap > 1e-4


def test_synthetic_shape():
    t = frbperiod.synthetic_periodic_arrivals(n=30, seed=0)
    assert t.size == 30
    assert np.all(np.diff(t) >= 0)  # sorted


def test_search_repeaters_skips_sparse():
    mjd = np.array([1.0, 2.0, 3.0, 100.0, 101.0])
    names = np.array(["A", "A", "A", "B", "B"])
    rows, _ = frbperiod.search_repeaters(mjd, names, min_bursts=3)
    by = {r["name"]: r for r in rows}
    assert by["A"]["searched"] and by["A"]["n"] == 3
    assert not by["B"]["searched"] and by["B"]["best_period"] is None


def test_run_offline_recovers_periodic_repeater(tmp_path):
    m = frbperiod.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    # the injected 16.35-day repeater is detected; the random one is not
    assert m["n_significant"] == 1
    det = m["detections"][0]
    assert det["name"] == "SYN-PER"
    assert abs(det["period_days"] - 16.35) < 0.2
    assert (tmp_path / "results" / "period_metrics.json").exists()
    # offline writes the per-source CSV to results/ (git-ignored), never the tracked survey/ showcase
    assert (tmp_path / "results" / "period_results_synthetic.csv").exists()
    assert not (tmp_path / "survey" / "period_results.csv").exists()
    assert (tmp_path / "papers" / "frbperiod" / "figures" / "periodogram.pdf").exists()
    assert (tmp_path / "papers" / "frbperiod" / "generated" / "macros.tex").exists()


def test_collapse_transits_merges_same_transit_bursts():
    t = np.array([100.0, 100.003, 100.01, 101.2, 105.0, 105.049])
    e = frbperiod.collapse_transits(t)
    assert e.size == 3
    assert abs(e[0] - np.mean([100.0, 100.003, 100.01])) < 1e-9
    # collapsing duplicate phases DEFLATES an inflated Z2 (the referee's finding, in miniature)
    rng = np.random.default_rng(0)
    base = np.sort(rng.uniform(0, 300, 12))
    dup = np.sort(np.concatenate([base, base[:5] + 0.003]))
    grid = np.linspace(2, 100, 2000)
    z_dup = frbperiod._z2_grid(dup, grid).max()
    z_col = frbperiod._z2_grid(frbperiod.collapse_transits(dup), grid).max()
    assert z_col < z_dup


def test_measured_duty_cycle():
    # bursts confined to a narrow phase window of a 16-day cycle
    d = frbperiod.measured_duty_cycle(np.array([0.1, 16.2, 32.05, 48.15]), 16.0)
    assert d < 0.3


def test_mc_null_and_summary_probabilities():
    grid = np.linspace(2, 100, 3000)
    null = frbperiod.mc_null(14, 260.0, grid, n_trials=150)
    s = frbperiod.null_summary(null, observed_z2=23.9, target_period=16.35, target_tol=0.15)
    assert s["p_exceed"] < 0.1  # a Z2 of 23.9 is rare under the uniform null
    assert s["p_coincidence"] < 0.05  # and the peak almost never lands on 16.35 by luck
    # a clustered donor null has a HEAVIER max-Z2 tail than the uniform one
    donor_iv = np.abs(np.random.default_rng(3).exponential(2.0, 40)) + 0.2
    clu = frbperiod.mc_null(14, 260.0, grid, donor_intervals=donor_iv, n_trials=150, seed=2)
    assert np.median(clu["max_z2"]) > np.median(null["max_z2"])


def test_detection_power_declines_with_duty_cycle():
    grid = np.linspace(2, 100, 3000)
    p = frbperiod.detection_power(11, 316.0, grid, duty_cycles=(0.05, 0.3), n_inj=60)
    tab = {row["duty_cycle"]: row["power"] for row in p}
    assert tab[0.05] > 0.9  # a tight window is easy
    assert tab[0.3] < tab[0.05]  # a wide window is hard -- the sensitivity statement is real


def test_period_bootstrap_gives_an_uncertainty():
    t = frbperiod.synthetic_periodic_arrivals(16.35, n=20, active_frac=0.2, span=300, seed=4)
    grid = np.linspace(2, 100, 4000)
    b = frbperiod.period_bootstrap(t, grid, n_boot=100)
    assert b["period_sd"] > 0.005  # larger than any "grid resolution" claim
    assert b["frac_near_peak"] > 0.7


def test_run_offline_emits_null_and_power_blocks(tmp_path):
    m = frbperiod.run(out=str(tmp_path), offline=True, mc_trials=100, n_inj=30)
    assert m["mc_nulls"]["uniform"]["p_coincidence"] is not None
    assert "clustered" in m["mc_nulls"]
    assert m["null_source_power"]["name"] == "SYN-RND"
    assert m["duty_cycle_at_peak"] < 0.5
    assert m["period_bootstrap"]["period_sd"] > 0
    # all searched sources are in the JSON, nulls included
    assert len(m["searched"]) == m["n_searchable"]
    macros = (tmp_path / "papers" / "frbperiod" / "generated" / "macros.tex").read_text()
    assert r"\fpCluPcoinc" in macros and r"\fpNullZtwo" in macros and r"\fpNepochs" in macros
    assert r"\fpNunsearchable" in macros and r"\fpPmin" in macros


def test_synthetic_run_does_not_clobber_real_figure(tmp_path):
    import json

    (tmp_path / "results").mkdir(parents=True)
    (tmp_path / "results" / "period_metrics.json").write_text(
        json.dumps({"source": "CHIME/FRB Catalog 1", "n_sources": 18})
    )
    figdir = tmp_path / "papers" / "frbperiod" / "figures"
    figdir.mkdir(parents=True)
    (figdir / "periodogram.pdf").write_bytes(b"REAL")
    frbperiod.run(out=str(tmp_path), offline=True, mc_trials=50, n_inj=20)
    assert (figdir / "periodogram.pdf").read_bytes() == b"REAL"
    kept = json.loads((tmp_path / "results" / "period_metrics.json").read_text())
    assert kept["n_sources"] == 18
