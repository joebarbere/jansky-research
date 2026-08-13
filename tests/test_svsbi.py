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
    # NTargets IS namespaced, as of 2026-08-12. This assertion previously required the
    # opposite and so locked the defect in: un-namespaced, an offline rebuild wrote the
    # synthetic parent size (400) into the macro the abstract uses for the real census (38).
    assert r"\newcommand{\svbSynNTargets}{400}" in macros
    assert r"\newcommand{\svbNTargets}" not in macros
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


def test_source_and_ntargets_are_namespaced():
    """These two were emitted once, outside the mode loop, from whichever leg ran.

    An offline rebuild therefore wrote the SYNTHETIC parent size (400 stars) into the same
    macro the abstract uses for the real census (38 M dwarfs) -- a wrong number, not a blank,
    so neither the '--' placeholder nor the arXiv assembler's check could catch it. That is
    the `\tiiNEvents` incident (768 real observing days -> 48 synthetic events) exactly.
    """
    import json
    from pathlib import Path

    macros = Path("papers/svsbi/generated/macros.tex")
    if not macros.exists():  # pragma: no cover - absent in a bare checkout
        pytest.skip("generated macros not present")
    text = macros.read_text()
    # the un-namespaced names must not exist at all
    assert r"\newcommand{\svbNTargets}" not in text
    assert r"\newcommand{\svbSource}" not in text
    assert r"\newcommand{\svbRealNTargets}" in text
    assert r"\newcommand{\svbSynNTargets}" in text
    # and the real one must match the committed census
    m = json.loads(Path("results/svsbi_metrics.json").read_text())
    assert rf"\newcommand{{\svbRealNTargets}}{{{m['n_targets']}}}" in text


def test_prior_sensitivity_is_committed_and_flags_the_break_as_prior_driven():
    """The experiment that decides what this paper may claim.

    A posterior median that follows its prior wall is reporting the box, not the data. The
    committed run must show that for log_Lbreak -- if it ever stops showing it, the paper's
    retraction of the "pins" claim needs revisiting, in either direction.
    """
    import json
    from pathlib import Path

    m = json.loads(Path("results/svsbi_metrics.json").read_text())
    ps = m.get("prior_sensitivity")
    assert ps, "prior_sensitivity is quoted in the paper and must be committed"
    pub, wide = ps["published_box"], ps["wide_box"]
    # the wide box must actually be wider, or the test is vacuous
    assert wide["prior_high"][1] > pub["prior_high"][1]
    assert wide["prior_low"][0] < pub["prior_low"][0]
    # the break must be flagged prior-driven, and by a margin over the seed noise
    assert ps["prior_driven"]["log_Lbreak"] is True
    shift = abs(ps["median_shift_wide_minus_published"]["log_Lbreak"])
    seed_sd = max(pub["median_seed_sd"]["log_Lbreak"], wide["median_seed_sd"]["log_Lbreak"])
    assert shift > 3 * seed_sd, f"shift {shift} not clear of seed noise {seed_sd}"
    # every parameter must have a per-seed record, or the seed scatter is not evidenced
    for box in (pub, wide):
        for name in ("lf_slope", "log_Lbreak", "f_beam"):
            assert len(box["median_per_seed"][name]) == len(ps["seeds"])


def test_prior_box_is_committed_so_width_ratios_are_checkable():
    """The paper quotes posterior/prior width ratios; without the bounds in the evidence file
    a reader has to open the source to check them."""
    import json
    from pathlib import Path

    m = json.loads(Path("results/svsbi_metrics.json").read_text())
    assert m.get("prior_low") and m.get("prior_high")
    assert m["theta_names"] == list(sv.THETA_NAMES)
    lo, hi = sv.prior_bounds()
    assert m["prior_low"] == [float(v) for v in lo]
    assert m["prior_high"] == [float(v) for v in hi]
