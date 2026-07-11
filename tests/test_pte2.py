"""Tests for jansky_research.pte2 -- PTE-II single-pulse heavy-tail census. Offline."""

from __future__ import annotations

import json
import sqlite3

import numpy as np

from jansky_research import pte2


def test_fit_lognormal_recovers_params():
    rng = np.random.default_rng(0)
    x = np.exp(rng.normal(2.0, 0.5, 5000))
    mu, sigma = pte2.fit_lognormal(x)
    assert abs(mu - 2.0) < 0.05 and abs(sigma - 0.5) < 0.05


def test_vuong_prefers_powerlaw_for_powerlaw_data():
    rng = np.random.default_rng(1)
    # a pure power-law tail (index 2.5) above xmin=10
    u = rng.uniform(0, 1, 4000)
    x = 10.0 * (1 - u) ** (-1.0 / (2.5 - 1.0))
    mu, sigma = pte2.fit_lognormal(x)
    v_pl, _ = pte2.vuong_powerlaw_vs_lognormal(x, 10.0, 2.5, mu, sigma)
    # and a pure log-normal tail
    y = np.exp(rng.normal(2.5, 0.5, 4000))
    yt = y[y >= 10.0]
    muy, sy = pte2.fit_lognormal(yt)
    v_ln, _ = pte2.vuong_powerlaw_vs_lognormal(yt, 10.0, 3.0, muy, sy)
    assert v_pl > v_ln  # power-law data favours power-law more than log-normal data does


def test_fit_energy_tail_classifies_heavy_and_lognormal():
    heavy = pte2.fit_energy_tail(pte2.synthetic_pulses(kind="heavy", n=800, seed=1))
    assert heavy["heavy_tailed"] and heavy["preferred"] == "power_law"
    assert heavy["n_giant"] >= 3 and heavy["p_excess"] < 0.05 and np.isfinite(heavy["gamma"])
    ln = pte2.fit_energy_tail(pte2.synthetic_pulses(kind="lognormal", n=800, seed=2))
    assert not ln["heavy_tailed"]


def test_fit_energy_tail_insufficient_below_cut():
    r = pte2.fit_energy_tail(np.exp(np.random.default_rng(0).normal(2, 0.5, 20)), min_pulses=50)
    assert r["preferred"] == "insufficient" and not r["heavy_tailed"] and r["n"] == 20


def test_synthetic_pulses_heavy_has_a_bigger_tail():
    ln = pte2.synthetic_pulses(kind="lognormal", n=3000, seed=4)
    hv = pte2.synthetic_pulses(kind="heavy", n=3000, seed=4)
    assert hv.max() > 3 * ln.max()  # injected giants dominate the extreme
    assert ln.min() >= 6.0 and hv.min() >= 6.0  # detection floor applied


def test_inject_recover_completeness_and_low_fp():
    rec = pte2.inject_recover_tail(counts=(100, 400, 800), n_each=25, seed=0)
    cov = rec["completeness_vs_count"]
    assert cov["n800"] >= 0.85  # high count -> the method recovers injected heavy tails
    assert cov["n800"] >= cov["n100"]  # honestly count-dependent
    assert rec["false_positive_rate"] < 0.15  # pure log-normals rarely called heavy


def test_census_sorts_and_drops_insufficient():
    per = {
        "J-heavy": {"sn": pte2.synthetic_pulses(kind="heavy", n=700, seed=1), "p0": 0.005},
        "J-ln": {"sn": pte2.synthetic_pulses(kind="lognormal", n=700, seed=2)},
        "J-small": {"sn": pte2.synthetic_pulses(kind="lognormal", n=30, seed=3)},
    }
    rows = pte2.census(per, min_pulses=50)
    names = [r["jname"] for r in rows]
    assert "J-small" not in names  # below the count cut
    assert rows[0]["excess"] >= rows[-1]["excess"]  # sorted by heavy-tail strength
    assert any(r["jname"] == "J-heavy" and r["heavy_tailed"] for r in rows)


def test_count_confound_detects_power_bias():
    # heavy sources deliberately given MORE pulses -> the confound test must flag count_limited
    rng = np.random.default_rng(0)
    rows = []
    for i in range(40):
        heavy = i >= 28
        rows.append(
            {
                "n": int(rng.integers(400, 900) if heavy else rng.integers(60, 300)),
                "heavy_tailed": heavy,
                "excess": 1.0 if heavy else 0.0,
                "gamma": 10.0 if heavy else float("nan"),
            }
        )
    c = pte2.count_confound(rows)
    assert c["n_heavy"] == 12 and c["n_fit"] == 40
    assert c["median_n_heavy"] > c["median_n_nonheavy"]
    assert c["heavy_frac_highcount"] > c["heavy_frac_lowcount"]
    assert c["count_limited"] and c["count_mw_p"] < 0.05
    assert abs(c["median_gamma_heavy"] - 10.0) < 1e-6
    # no bias -> not count_limited
    flat = [
        {"n": 200, "heavy_tailed": i % 2 == 0, "excess": 0.0, "gamma": float("nan")}
        for i in range(20)
    ]
    assert not pte2.count_confound(flat).get("count_limited", False)


def test_tail_vs_edot_correlation_keys():
    # build a census where heavy-tailed sources have higher Edot
    rng = np.random.default_rng(0)
    rows, edot = [], {}
    for i in range(20):
        heavy = i >= 12
        rows.append(
            {
                "jname": f"J{i}",
                "excess": (1.0 if heavy else 0.0) + rng.normal(0, 0.1),
                "heavy_tailed": heavy,
            }
        )
        edot[f"J{i}"] = 10 ** (34.0 + (2.0 if heavy else 0.0) + rng.normal(0, 0.2))
    out = pte2.tail_vs_edot(rows, edot)
    assert out["n_matched"] == 20
    assert out["spearman_excess_logedot"] > 0.3  # excess tracks Edot by construction
    assert "mannwhitney_p" in out and out["logedot_heavy_median"] > out["logedot_nonheavy_median"]


def test_load_pulse_sn_from_synthetic_sqlite(tmp_path):
    # a minimal DB mirroring the documented PTE-II schema
    db = tmp_path / "mini.db"
    con = sqlite3.connect(db)
    con.executescript(
        """
        CREATE TABLE pulsar (pulsarID INTEGER PRIMARY KEY, jname TEXT, p0 REAL, s1400 REAL);
        CREATE TABLE file (pfLinkID INTEGER PRIMARY KEY, pulsarID INTEGER, timeStartMJD REAL);
        CREATE TABLE fileSegment (segID INTEGER PRIMARY KEY, pfLinkID INTEGER, snr_max REAL);
        """
    )
    con.execute("INSERT INTO pulsar VALUES (1,'J0001+0001',0.5,2.0)")
    con.execute("INSERT INTO pulsar VALUES (2,'J0002+0002',0.005,5.0)")
    con.execute("INSERT INTO file VALUES (10,1,50000.0)")
    con.execute("INSERT INTO file VALUES (11,2,50001.0)")
    for i, sn in enumerate([7.0, 9.0, 30.0]):
        con.execute("INSERT INTO fileSegment VALUES (?,10,?)", (100 + i, sn))
    con.execute("INSERT INTO fileSegment VALUES (200,11,8.0)")
    con.commit()
    con.close()
    out = pte2.load_pulse_sn(db)
    assert set(out) == {"J0001+0001", "J0002+0002"}
    assert out["J0001+0001"]["n"] == 3 and out["J0001+0001"]["p0"] == 0.5
    assert sorted(out["J0001+0001"]["sn"].tolist()) == [7.0, 9.0, 30.0]


def test_run_offline_recovers_and_writes(tmp_path):
    m = pte2.run(str(tmp_path), offline=True)
    assert m["recovered"]  # heavy classified heavy, log-normal not
    assert m["false_positive_rate"] < 0.2
    saved = json.loads((tmp_path / "results" / "pte2_metrics.json").read_text())
    assert saved["source"] == m["source"]
    assert (tmp_path / "papers" / "pte2" / "figures" / "pte2.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "pte2" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\ptSynRecovered}{yes}" in macros
    assert r"\newcommand{\ptRealNFit}{--}" in macros  # real namespace placeholder offline


def test_write_macros_real_namespace(tmp_path):
    real = {
        "source": "PTE-II",
        "is_real": True,
        "n_fit": 210,
        "n_heavy": 34,
        "heavy_fraction": 0.162,
        "n_matched": 205,
        "spearman_excess_logedot": 0.28,
        "spearman_p": 3e-5,
        "logedot_heavy_median": 34.9,
        "logedot_nonheavy_median": 33.1,
        "mannwhitney_p": 1e-4,
    }
    p = tmp_path / "m.tex"
    pte2._write_macros(real, p)
    txt = p.read_text()
    assert r"\newcommand{\ptRealNHeavy}{34}" in txt
    assert r"\newcommand{\ptRealSpearman}{0.28}" in txt
    assert r"\newcommand{\ptSynRecovered}{yes}" in txt  # synthetic still live in a real build
