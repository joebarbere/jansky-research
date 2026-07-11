"""Tests for jansky_research.glitchpop -- JBO glitch waiting-time census. Offline."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from jansky_research import glitchpop as gp

FIXTURE = Path(__file__).parent / "data" / "jbo_gtable_sample.html"


def test_parse_glitch_table_real_fixture():
    rows = gp.parse_glitch_table(FIXTURE.read_text())
    assert len(rows) == 25  # the vendored real rows
    r0 = rows[0]
    assert r0["jname"] == "J0007+7303" and r0["mjd"] == 54953.0 and r0["size"] == 554.0
    # a 'New' provenance flag is picked up
    assert any(r["is_new"] for r in rows)
    # missing ('X') fields become NaN, not a crash
    assert any(np.isnan(r["size"]) or np.isnan(r["dnudot"]) for r in rows)


def test_group_by_pulsar_sorts_by_mjd():
    glitches = [
        {"jname": "JA", "mjd": 55000.0, "size": 5.0, "dnudot": 1.0, "refs": "", "is_new": False},
        {"jname": "JA", "mjd": 54000.0, "size": 3.0, "dnudot": 1.0, "refs": "", "is_new": False},
        {"jname": "JB", "mjd": 56000.0, "size": 9.0, "dnudot": 1.0, "refs": "", "is_new": False},
    ]
    by = gp.group_by_pulsar(glitches)
    assert set(by) == {"JA", "JB"}
    assert by["JA"]["n"] == 2 and by["JA"]["mjd"][0] == 54000.0  # sorted ascending


def test_waiting_times_years():
    w = gp.waiting_times(np.array([50000.0, 50000.0 + 365.25, 50000.0 + 3 * 365.25]))
    assert np.allclose(w, [1.0, 2.0])


def test_waiting_time_fit_classifies_three_processes():
    assert (
        gp.classify_pulsar(gp.synthetic_glitch_series(kind="exponential", n=30, seed=3))["klass"]
        == "exponential"
    )
    qp = gp.waiting_time_fit(
        gp.waiting_times(gp.synthetic_glitch_series(kind="quasi_periodic", n=30, seed=3)["mjd"])
    )
    assert qp["klass"] == "quasi_periodic" and qp["cv"] < 0.5 and qp["p_regular"] < 0.05
    # clustered is detectable only with monitoring-gap excision OFF (its long gaps look like gaps)
    cl = gp.waiting_time_fit(
        gp.waiting_times(gp.synthetic_glitch_series(kind="clustered", n=30, seed=3)["mjd"]),
        gap_factor=1e9,
    )
    assert cl["klass"] == "clustered" and cl["cv"] > 1 and cl["p_clustered"] < 0.05


def test_monitoring_gap_excision_recovers_regularity():
    # a regular (quasi-periodic) series with one huge monitoring gap injected: the raw CV is wrecked,
    # but gap excision recovers the quasi-periodic class (the J0537 scenario)
    s = gp.synthetic_glitch_series(kind="quasi_periodic", n=25, mean_yr=0.3, seed=1)
    mjd = s["mjd"].copy()
    mjd[13:] += 2000.0  # inject a ~5.5-yr monitoring gap partway through
    w = gp.waiting_times(mjd)
    naive = gp.waiting_time_fit(w, gap_factor=1e9)  # no excision -> gap wrecks it
    fixed = gp.waiting_time_fit(w)  # default excision -> recovered
    assert fixed["n_gaps_excised"] >= 1
    assert fixed["klass"] == "quasi_periodic"
    assert naive["klass"] != "quasi_periodic"  # the naive (no-excision) fit is misled by the gap


def test_waiting_time_fit_insufficient():
    r = gp.waiting_time_fit(np.array([1.0, 2.0]), min_waits=4)
    assert r["klass"] == "insufficient" and r["n_waits"] == 2


def test_classify_pulsar_below_cut():
    d = gp.synthetic_glitch_series(kind="exponential", n=3, seed=0)
    assert gp.classify_pulsar(d, min_glitches=5)["klass"] == "insufficient"


def test_population_census_sorts_and_drops():
    by = {
        "J-big": gp.synthetic_glitch_series(kind="quasi_periodic", n=25, seed=1),
        "J-mid": gp.synthetic_glitch_series(kind="exponential", n=8, seed=2),
        "J-small": gp.synthetic_glitch_series(kind="exponential", n=3, seed=3),
    }
    rows = gp.population_census(by, min_glitches=5)
    names = [r["jname"] for r in rows]
    assert "J-small" not in names  # below the cut
    assert rows[0]["n"] >= rows[-1]["n"]  # sorted by glitch count
    assert next(r for r in rows if r["jname"] == "J-big")["klass"] == "quasi_periodic"


def test_classification_delta_newly_and_flip():
    # a pulsar with only 3 pre-2019 glitches but 6 total -> newly qualified
    pre = 50000.0 + np.arange(3) * 400.0  # all before split
    post = gp.BASU_END_MJD + np.arange(3) * 400.0
    mjd = np.concatenate([pre, post])
    by = {
        "J-new": {"mjd": mjd, "size": np.full(6, 5.0), "n": 6},
        # a pulsar with >=5 both epochs, quasi-periodic throughout -> stable, no flip
        "J-stable": gp.synthetic_glitch_series(
            kind="quasi_periodic", n=14, mean_yr=1.0, start_mjd=54000.0, seed=1
        ),
    }
    for d in by.values():
        d["mjd"] = np.asarray(d["mjd"], float)
        d["size"] = np.asarray(d.get("size", np.full(d["mjd"].size, 5.0)), float)
    delta = gp.classification_delta(by, min_glitches=5)
    assert delta["n_newly_classifiable"] >= 1
    assert any(x["jname"] == "J-new" for x in delta["newly_classifiable"])
    assert delta["n_stable_sample"] >= 1


def test_population_significance_flags_aggregate_excess():
    # 10 quasi-periodic among 33 vs ~1.65 expected -> highly significant aggregate excess
    rows = [{"klass": "quasi_periodic"}] * 10 + [{"klass": "exponential"}] * 23
    s = gp.population_significance(rows)
    assert abs(s["expected_false_qp"] - 33 * 0.05) < 1e-6
    assert s["qp_excess_significant"] and s["qp_binomial_p"] < 1e-3
    # a fraction at the chance level -> not significant
    chance = [{"klass": "quasi_periodic"}] * 1 + [{"klass": "exponential"}] * 19
    assert not gp.population_significance(chance)["qp_excess_significant"]


def test_inject_recover_accuracy_and_low_fp():
    rec = gp.inject_recover(counts=(6, 20, 40), n_each=25, seed=0)
    assert rec["exponential_accuracy_vs_count"]["n40"] >= 0.85
    assert rec["quasiperiodic_completeness_vs_count"]["n40"] >= 0.85
    assert rec["exponential_false_positive_rate"] < 0.2


def test_run_offline_recovers_and_writes(tmp_path):
    m = gp.run(str(tmp_path), offline=True)
    assert m["recovered"]  # all three processes classified correctly
    assert m["exponential_false_positive_rate"] < 0.2
    saved = json.loads((tmp_path / "results" / "glitchpop_metrics.json").read_text())
    assert saved["source"] == m["source"]
    assert (tmp_path / "papers" / "glitchpop" / "figures" / "glitchpop.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "glitchpop" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\gpSynRecovered}{yes}" in macros
    assert r"\newcommand{\gpRealNFit}{--}" in macros  # real placeholder offline


def test_write_macros_real_namespace(tmp_path):
    real = {
        "source": "JBO",
        "is_real": True,
        "n_glitches": 727,
        "n_pulsars": 222,
        "n_qualified_full": 33,
        "n_exponential": 24,
        "n_quasiperiodic": 4,
        "n_clustered": 5,
        "n_newly_classifiable": 6,
        "n_stable_sample": 27,
        "n_flipped": 3,
        "n_magnetars_dropped": 3,
        "expected_false_qp": 1.55,
        "qp_binomial_p": 1e-6,
        "known_quasiperiodic_ok": "yes",
    }
    p = tmp_path / "m.tex"
    gp._write_macros(real, p)
    txt = p.read_text()
    assert r"\newcommand{\gpRealNFlipped}{3}" in txt
    assert r"\newcommand{\gpRealNQp}{4}" in txt
    assert r"\newcommand{\gpRealNMagnetars}{3}" in txt
    assert r"\newcommand{\gpRealKnownQpOK}{yes}" in txt
    assert r"\newcommand{\gpSynRecovered}{yes}" in txt  # synthetic still live in a real build
