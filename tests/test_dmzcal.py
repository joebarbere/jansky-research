"""Tests for jansky_research.dmzcal -- the p(z|DM) coverage audit (plan 99). Offline."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import dmzcal as dz

SMALL = {"nz": 120, "zmax": 3.0, "dm_step": 5.0, "dm_max": 12000.0}


@pytest.fixture(scope="module")
def table() -> dz.LikelihoodTable:
    return dz.build_table(dz.HOFFMANN_EMIN25, **SMALL)


def test_mean_dm_cosmic_matches_macquart_scale():
    # ~ 900-1000 pc/cc per unit z at z ~ 1 (Macquart+2020 fig. 3); monotone; zero at z = 0
    m = dz.mean_dm_cosmic(np.array([0.0, 0.5, 1.0, 2.0]))
    assert m[0] == pytest.approx(0.0, abs=1e-9)
    assert np.all(np.diff(m) > 0)
    assert 850 < m[2] < 1050


def test_solve_c0_gives_unit_mean():
    for sigma in (0.2, 0.5, 1.5):
        d = np.geomspace(1e-3, 200, 6000)
        c0 = dz.solve_c0(sigma, d)
        pdf = dz._pdelta_unnorm(d, sigma, c0)
        assert np.trapezoid(pdf * d, d) / np.trapezoid(pdf, d) == pytest.approx(1.0, abs=1e-6)


def test_solve_c0_extends_bracket_at_large_sigma():
    d = np.geomspace(1e-3, 200, 6000)
    c0 = dz.solve_c0(40.0, d)
    assert c0 < -100
    pdf = dz._pdelta_unnorm(d, 40.0, c0)
    assert np.trapezoid(pdf * d, d) / np.trapezoid(pdf, d) == pytest.approx(1.0, abs=1e-6)


def test_bin_probs_sum_and_host_redshifting():
    edges = np.arange(0, 6001, 2.0)
    pc = dz.cosmic_bin_probs(0.5, edges)
    assert 0.99 < pc.sum() <= 1.0 + 1e-12
    mids = edges[:-1] + 1
    h0, h1 = dz.host_bin_probs(0.0, edges), dz.host_bin_probs(1.0, edges)
    assert h0.sum() == pytest.approx(1.0, abs=1e-6)
    # observed host DM halves at z = 1
    assert (h1 * mids).sum() / (h0 * mids).sum() == pytest.approx(0.5, rel=0.02)


def test_table_mean_tracks_cosmic_plus_host(table):
    mids = table.edges[:-1] + 0.5 * table.dm_step
    p = dz.HOFFMANN_EMIN25
    host_mean = 10 ** (p.lmean + 0.5 * np.log(10) * p.lsigma**2)
    for zz in (0.1, 1.0):
        i = int(np.argmin(abs(table.z - zz)))
        want = dz.mean_dm_cosmic(np.array([table.z[i]]))[0] + host_mean / (1 + table.z[i])
        assert (table.prob[i] * mids).sum() == pytest.approx(want, rel=0.03)


def test_posterior_normalised_and_moves_with_dm(table):
    post = dz.posterior(table, np.array([300.0, 1500.0]))
    assert post.sum(axis=0) == pytest.approx([1.0, 1.0])
    zmean = (post * table.z[:, None]).sum(axis=0)
    assert zmean[1] > zmean[0] + 0.5
    assert dz.posterior(table, 300.0).shape == table.z.shape


def test_interval_contains_median_and_widens_with_q(table):
    lo68, hi68 = dz.central_interval(table, 800.0, 0.68)
    lo95, hi95 = dz.central_interval(table, 800.0, 0.95)
    assert lo95 < lo68 < hi68 < hi95
    p = dz.pit(table, np.array([800.0, 800.0]), np.array([lo68, hi68]))
    assert p == pytest.approx([0.16, 0.84], abs=0.01)


def test_null_pits_are_uniform(table):
    rng = np.random.default_rng(0)
    z, dm = dz.simulate_sample(table, 3000, rng)
    pits = dz.pit(table, dm, z)
    assert dz.coverage(pits, 0.68) == pytest.approx(0.68, abs=0.03)
    assert dz.coverage(pits, 0.95) == pytest.approx(0.95, abs=0.015)


def test_rule_flags_overconfidence_and_bias():
    rule = dz.CalibrationRule(bands={0.68: (0.55, 0.80), 0.95: (0.88, 1.0)})
    rng = np.random.default_rng(1)
    assert rule.evaluate(rng.random(200))["verdict"] == "CALIBRATED"
    over = np.concatenate([rng.random(100) * 0.05, 0.95 + rng.random(100) * 0.05])
    assert rule.evaluate(over)["verdict"] == "OVERCONFIDENT"
    under = 0.5 + (rng.random(200) - 0.5) * 0.5
    assert rule.evaluate(under)["verdict"] == "UNDERCONFIDENT"
    biased = np.clip(rng.random(200) * 0.75 + 0.2, 0, 1)  # in-band coverage, shifted
    out = rule.evaluate(biased)
    assert not out["passed"]


def test_null_rule_bands_bracket_nominal():
    rng = np.random.default_rng(2)
    rule = dz.null_rule(rng.random((400, 36)))
    lo, hi = rule.bands[0.68]
    assert lo < 0.68 < hi


def test_run_controls_small(monkeypatch):
    out = dz.run_controls(n=20, m=30, seed=5, table_kw=SMALL)
    assert set(out) >= {"rule", "C0", "C1"}
    assert 0.0 <= out["C0"]["false_fail_rate"] <= 1.0
    # F x 2 is the easy case: power must be well above the false-fail rate
    assert out["C1"]["F_x2"]["power"] > out["C0"]["false_fail_rate"]


def test_run_controls_widens_when_gate_fails():
    out = dz.run_controls(n=20, m=30, seed=6, max_false_fail=-1.0, table_kw=SMALL)
    assert out["rule"]["widened"] is True
    assert out["rule"]["band_quantile"] == 0.99


def test_dedupe_merges_only_matching_dms():
    js = [
        {"_file": "FRB20210320", "DM": {"value": 384.8}},
        {"_file": "FRB20210320C", "DM": {"value": 384.6}, "z": 0.28},
        {"_file": "FRB20240101A", "DM": 300.0},
        {"_file": "FRB20240101B", "DM": 900.0},
    ]
    kept, merges = dz.dedupe(js)
    names = sorted(j["_file"] for j in kept)
    assert names == ["FRB20210320C", "FRB20240101A", "FRB20240101B"]
    assert merges == [{"date": "20210320", "kept": "FRB20210320C", "dropped": ["FRB20210320"]}]


def test_dedupe_rejects_nameless():
    with pytest.raises(ValueError):
        dz.dedupe([{"_file": "FRBX"}])


@pytest.mark.parametrize(
    ("tel", "name", "side"),
    [
        ("DSA", "FRB20220207C", "production"),
        ("DSA", "FRB20220307B", "production"),
        ("DSA", "FRB20230101A", "certification"),
        ("ASKAP", "FRB20231231A", "production"),
        ("ASKAP", "FRB20240101A", "certification"),
        ("CHIME", "FRB20231204A", "certification"),
        ("MeerKAT", "FRB20240304B", "certification"),
        ("VLA", "FRB20121102A", "no_model"),
        ("", "FRB20181119A", "no_model"),
    ],
)
def test_provenance_side(tel, name, side):
    s, why = dz.provenance_side(tel, name, ["20220207C", "20220307B"], ["20220207C"], "20231231")
    assert s == side
    if name == "FRB20220307B":
        assert "DM only" in why


def test_solve_c0_gives_up_eventually():
    with pytest.raises(ValueError):
        dz.solve_c0(1e5)


def test_tail_fractions_one_sided():
    pits = np.concatenate([np.full(44, 0.05), np.full(56, 0.5)])
    t = dz.tail_fractions(pits)
    assert t["below"] == pytest.approx(0.44)
    assert t["above"] == 0.0
    assert t["expected_each"] == pytest.approx(0.16)


def test_dm_pit_uniform_on_null(table):
    rng = np.random.default_rng(7)
    z, dm = dz.simulate_sample(table, 3000, rng)
    p = dz.dm_pit(table, dm, z)
    assert dz.coverage(p, 0.68) == pytest.approx(0.68, abs=0.03)
    assert np.median(p) == pytest.approx(0.5, abs=0.03)


def test_dm_pit_detects_host_excess(table):
    rng = np.random.default_rng(8)
    z, dm = dz.simulate_sample(table, 500, rng)
    shifted = dz.dm_pit(table, dm + 300.0, z)
    assert np.median(shifted) > 0.7


RULE = dz.CalibrationRule(bands={0.68: (0.55, 0.80), 0.95: (0.88, 1.0)})


def test_rule_names_biased_when_coverage_in_band():
    # 70% in [0.16, 0.84] (60% of it high), 95% in [0.025, 0.975]; median ~0.62
    pits = np.concatenate(
        [
            np.linspace(0.56, 0.83, 120),
            np.linspace(0.17, 0.5, 20),
            np.linspace(0.85, 0.97, 50),
            np.linspace(0.03, 0.15, 0),
            np.full(10, 0.99),
        ]
    )
    out = RULE.evaluate(pits)
    assert not out["passed"]
    assert out["verdict"] == "BIASED"


def test_rule_names_ks_fail_for_clumped_symmetric_pits():
    pits = np.concatenate(
        [
            np.full(140, 0.2),
            np.full(140, 0.8),
            np.linspace(0.03, 0.15, 54),
            np.linspace(0.85, 0.97, 54),
            np.full(6, 0.01),
            np.full(6, 0.99),
        ]
    )
    out = RULE.evaluate(pits)
    assert out["verdict"] == "KS_FAIL"


def test_rule_precedence_overconfident_before_biased():
    # one-sided low PITs: both coverage and location fail; the label is OVERCONFIDENT
    pits = np.concatenate([np.full(60, 0.01), np.linspace(0.2, 0.8, 40)])
    out = RULE.evaluate(pits)
    assert out["verdict"] == "OVERCONFIDENT"
    assert out["median_pit"] < 0.4
