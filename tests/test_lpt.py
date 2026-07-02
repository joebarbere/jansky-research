"""Tests for jansky_research.lpt -- the LPT population catalogue. Offline (vendored CSV)."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import lpt


def test_load_sample_shapes_and_flags():
    s = lpt.load_sample()
    assert s["period_s"].size == 13
    assert (s["period_s"] > 60).all()  # every period > 1 minute (unit sanity)
    assert s["pdot_is_measurement"].sum() == 2  # CHIME J0630+25 + CHIME/ILT J1634+44
    assert s["is_wd_binary"].sum() == 6
    # the one clear spin-up is negative
    neg = s["pdot"][s["pdot_is_measurement"]] < 0
    assert neg.sum() == 1


def test_population_table_death_line_headline():
    s = lpt.load_sample()
    pop = lpt.population_table(s)
    assert pop["n_lpt"] == 13
    # every Pdot-constrained object sits below the pulsar death line -- the class puzzle
    assert pop["n_below_death_line"] == pop["n_pdot_constrained"] == 9
    assert pop["period_min_min"] == 7.0 and pop["period_max_hr"] > 6


def test_period_split_stat_honest_at_small_n():
    s = lpt.load_sample()
    out = lpt.period_split_stat(s["period_s"], s["is_wd_binary"])
    assert out["delta_log_median"] > 0  # WD binaries do sit at longer periods...
    assert out["p_perm"] > 0.05  # ...but NOT significantly at N=13 (the honest result)


def test_split_stat_round_trips_injected_split():
    p, wd = lpt.synthetic_lpt_population(seed=0)
    out = lpt.period_split_stat(p, wd)
    assert out["p_perm"] < 0.05  # a REAL split registers -> the real non-detection is informative


def test_split_stat_degenerate():
    out = lpt.period_split_stat(np.array([100.0, 200.0]), np.array([True, False]))
    assert np.isnan(out["p_perm"])


def test_run_writes_artifacts(tmp_path):
    m = lpt.run(str(tmp_path), offline=True)
    assert m["n_lpt"] == 13 and m["n_pdot_measurements"] == 2
    saved = json.loads((tmp_path / "results" / "lpt_metrics.json").read_text())
    assert saved == m
    assert (tmp_path / "papers" / "lpt" / "figures" / "lpt_ppdot.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "lpt" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\lptNdeath}{9}" in macros


def test_write_macros_placeholder(tmp_path):
    p = tmp_path / "m.tex"
    lpt._write_macros({"n_lpt": None}, p)
    assert r"\newcommand{\lptN}{--}" in p.read_text()
