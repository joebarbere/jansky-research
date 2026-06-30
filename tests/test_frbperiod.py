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
    assert (tmp_path / "results" / "period_results.csv").exists()
    assert not (tmp_path / "survey" / "period_results.csv").exists()
    assert (tmp_path / "papers" / "frbperiod" / "figures" / "periodogram.pdf").exists()
    assert (tmp_path / "papers" / "frbperiod" / "generated" / "macros.tex").exists()
