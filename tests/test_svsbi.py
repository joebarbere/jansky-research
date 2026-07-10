"""Tests for jansky_research.svsbi -- the physics/forward-model leg (pure NumPy, core CI).

The NPE/SBC leg needs the `sbi` extra and is exercised from the ROCm venv, not here.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from jansky_research import svsbi as sv


def _parent(n=60, seed=0):
    return sv._synthetic_parent(n=n, seed=seed)


def test_prior_bounds_ordered():
    lo, hi = sv.prior_bounds()
    assert np.all(hi > lo) and lo.shape == (3,)


def test_luminosity_draw_respects_slope_and_break():
    rng = np.random.default_rng(0)
    steep = sv._sample_luminosity(20000, 3.0, 13.5, rng)
    flat = sv._sample_luminosity(20000, 1.5, 13.5, rng)
    # a steeper slope -> more weight at low L -> lower median luminosity
    assert np.median(steep) < np.median(flat)
    # the break caps the bright end near 10^13.5
    assert np.percentile(flat, 99) < 10 ** (13.5 + 1.0)


def test_draw_population_beaming_fraction():
    rng = np.random.default_rng(1)
    d = np.full(5000, 10.0)
    pop_hi = sv.draw_population(np.array([2.0, 13.5, 0.5]), d, rng)
    pop_lo = sv.draw_population(np.array([2.0, 13.5, 0.05]), d, rng)
    assert 0.45 < pop_hi["beams"].mean() < 0.55
    assert pop_lo["beams"].mean() < 0.1
    # non-beaming stars have zero observed flux
    assert np.all(pop_hi["flux_mjy"][~pop_hi["beams"]] == 0.0)


def test_forward_model_more_beaming_more_detections():
    parent = _parent()
    rng = np.random.default_rng(2)

    def ndet(fb):
        return np.mean(
            [
                sv.forward_model(np.array([2.0, 14.0, fb]), parent, rng)["detected"].sum()
                for _ in range(120)
            ]
        )

    assert ndet(0.5) > ndet(0.05)


def test_forward_model_brighter_lf_more_detections():
    parent = _parent()
    rng = np.random.default_rng(3)

    def ndet(loglb):
        return np.mean(
            [
                sv.forward_model(np.array([2.0, loglb, 0.5]), parent, rng)["detected"].sum()
                for _ in range(120)
            ]
        )

    assert ndet(14.5) > ndet(12.5)


def test_summary_stats_shape_and_empty():
    empty = sv.summary_stats({"detected": np.zeros(10, bool), "v_obs": np.zeros(10)})
    assert empty.shape == (4,) and empty[0] == 0
    det = np.array([True, False, True])
    s = sv.summary_stats({"detected": det, "v_obs": np.array([2.0, 0.1, 5.0])})
    assert s[0] == 2 and s[1] == pytest.approx(np.log10(5.0), abs=1e-4)


def test_simulate_deterministic_in_seed():
    parent = _parent()
    a = sv.simulate(np.array([2.0, 13.5, 0.3]), parent, seed=7)
    b = sv.simulate(np.array([2.0, 13.5, 0.3]), parent, seed=7)
    assert np.array_equal(a, b)


def test_observed_summary_detects_bright_v():
    parent = {
        "v_best_obs": np.array([5.0, 0.1, 0.2]),
        "v_rms": np.array([0.2, 0.2, 0.2]),
        "leakage_floor_mjy": np.array([0.1, 0.1, 0.1]),
    }
    s = sv.observed_summary(parent)
    assert s[0] == 1  # only the 5 mJy source clears 5*rms and the floor


def test_parent_from_census_parses_and_dedupes(tmp_path):
    # a self-contained census fixture (the real merged CSV is git-ignored, so absent in a fresh
    # CI checkout) exercising the parser + the GJ-65-style unresolved-binary dedup
    import csv as _csv

    rows = [
        # GJ 65: an unresolved binary (two names, byte-identical bright photometry) that IS the
        # detection -- deduped to one physical source (dedup fires only for bright v_best > 1)
        {
            "name": "GJ65A",
            "gaia_id": "4",
            "i_mjy": "3.0",
            "e_i": "0.2",
            "v_mjy": "5.0",
            "e_v": "0.2",
        },
        {
            "name": "GJ65A",
            "gaia_id": "4",
            "i_mjy": "3.0",
            "e_i": "0.2",
            "v_mjy": "0.1",
            "e_v": "0.2",
        },
        {
            "name": "GJ65B",
            "gaia_id": "5",
            "i_mjy": "3.0",
            "e_i": "0.2",
            "v_mjy": "5.0",
            "e_v": "0.2",
        },
        {
            "name": "GJ65B",
            "gaia_id": "5",
            "i_mjy": "3.0",
            "e_i": "0.2",
            "v_mjy": "0.1",
            "e_v": "0.2",
        },
        # two quiet targets (no detection), 2 epochs each
        {"name": "Q1", "gaia_id": "2", "i_mjy": "1.0", "e_i": "0.2", "v_mjy": "0.1", "e_v": "0.2"},
        {"name": "Q1", "gaia_id": "2", "i_mjy": "1.0", "e_i": "0.2", "v_mjy": "0.15", "e_v": "0.2"},
        {"name": "Q2", "gaia_id": "3", "i_mjy": "1.0", "e_i": "0.2", "v_mjy": "0.05", "e_v": "0.2"},
        {"name": "Q2", "gaia_id": "3", "i_mjy": "1.0", "e_i": "0.2", "v_mjy": "0.2", "e_v": "0.2"},
    ]
    p = tmp_path / "census.csv"
    with open(p, "w", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=["name", "gaia_id", "i_mjy", "e_i", "v_mjy", "e_v"])
        w.writeheader()
        w.writerows(rows)
    parent = sv.parent_from_census(p, fetch_distances=False)
    assert parent["v_rms"].size == 3  # GJ65A/GJ65B deduped -> 3 physical targets (GJ65, Q1, Q2)
    assert np.all(parent["v_rms"] > 0)
    assert parent["v_best_obs"].size == parent["v_rms"].size
    assert sv.observed_summary(parent)[0] == 1  # only the deduped GJ 65 clears 5*rms and the floor


def test_run_offline_writes_artifacts(tmp_path):
    m = sv.run(str(tmp_path), offline=True)
    assert m["source"].startswith("synthetic")
    assert m["beaming_monotonic"] and m["luminosity_monotonic"]
    assert m["ndet_high_beaming"] > m["ndet_low_beaming"]
    saved = json.loads((tmp_path / "results" / "svsbi_metrics.json").read_text())
    assert saved["n_targets"] == m["n_targets"]
    assert (tmp_path / "papers" / "svsbi" / "figures" / "svsbi.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "svsbi" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\svbSynNTargets}" not in macros  # NTargets is namespace-free
    assert r"\newcommand{\svbRealFbeam}{--}" in macros


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    sv._write_macros(
        {"source": "x", "is_real": True, "n_targets": 60, "posterior_median": {"f_beam": None}}, p
    )
    txt = p.read_text()
    assert r"\newcommand{\svbRealFbeam}{--}" in txt and r"\newcommand{\svbSynFbeam}{--}" in txt


def test_ks_uniform_zero_for_uniform_ranks():
    ranks = np.arange(150)  # perfectly uniform
    assert sv._ks_uniform(ranks, 150) < 0.02
