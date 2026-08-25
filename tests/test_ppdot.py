"""Tests for jansky_research.ppdot -- the pulsar P-Pdot diagram. No network."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import ppdot


def test_magnetic_field_and_age_crab():
    # Crab pulsar: P = 0.0334 s, Pdot = 4.2e-13 -> B ~ 3.8e12 G, tau ~ 1260 yr
    b = ppdot.magnetic_field(0.0334, 4.2e-13)
    tau = ppdot.characteristic_age(0.0334, 4.2e-13)
    assert 3.0e12 < b < 4.5e12
    assert 1100 < tau < 1400


def test_spindown_luminosity_positive_scaling():
    e1 = ppdot.spindown_luminosity(0.1, 1e-15)
    e2 = ppdot.spindown_luminosity(0.1, 2e-15)
    assert e2 > e1 > 0  # Edot ~ Pdot


def test_death_line_orders_alive_dead():
    # a normal pulsar sits well above the death line; a long-period low-Pdot one below
    assert 1e-15 > ppdot.death_line(0.5)  # alive
    assert 1e-16 < ppdot.death_line(8.0)  # dead (death-line Pdot exceeds the source's)


def test_classify_three_populations():
    # normal, millisecond, magnetar archetypes
    p = np.array([0.5, 0.005, 5.0])
    pd = np.array([1.6e-15, 3e-20, 1e-11])
    cls = ppdot.classify(p, pd)
    assert cls[0] == "normal" and cls[1] == "msp" and cls[2] == "magnetar"


def test_synthetic_population_separates_in_field():
    pop = ppdot.synthetic_population(n_each=400, seed=1)
    stats = ppdot.population_stats(pop["period_s"], pop["pdot"])
    # the three classes separate cleanly in median surface field
    assert stats["msp"]["median_log_b"] < stats["normal"]["median_log_b"]
    assert stats["normal"]["median_log_b"] < stats["magnetar"]["median_log_b"]
    assert stats["msp"]["median_log_b"] < 10.0  # recycled, low field
    assert stats["magnetar"]["median_log_b"] > 13.0


def test_run_offline(tmp_path):
    m = ppdot.run(out=str(tmp_path), offline=True)
    assert m["source"] == "synthetic"
    assert m["n_pulsars"] == 1200
    assert m["classify_accuracy"] > 0.9  # the injected classes are recovered
    assert m["median_log_b_msp"] < m["median_log_b_normal"] < m["median_log_b_magnetar"]
    assert m["frac_above_death"] > 0.9  # the injected populations are radio-loud
    assert (tmp_path / "results" / "ppdot_metrics.json").exists()
    assert (tmp_path / "papers" / "ppdot" / "figures" / "ppdot.pdf").exists()
    macros = (tmp_path / "papers" / "ppdot" / "generated" / "macros.tex").read_text()
    assert r"\ppLogBmsp" in macros and r"\ppFracAlive" in macros


def test_named_pulsar_derived_reads_the_row_not_a_constant():
    """The Crab check must run on the fetched table, or it tests only the algebra.

    Before this, `fetch_atnf_ppdot` requested PSRJ and discarded it, so no named pulsar could be
    located and the paper's Crab numbers were typed from the literature. A check that cannot see
    the catalogue cannot catch a units, parsing or cut error on the real path.
    """
    pop = {
        "name": np.array(["J0534+2200", "J0537-6910", "J1939+2134"], dtype=object),
        "period_s": np.array([0.0333924, 0.016, 0.00156]),
        "pdot": np.array([4.21e-13, 5.2e-14, 1.05e-19]),
    }
    crab = ppdot.named_pulsar_derived(pop, ppdot.CRAB_PSRJ)
    assert crab["found"]
    assert crab["b_gauss"] == pytest.approx(3.8e12, rel=0.02)
    assert crab["age_yr"] == pytest.approx(1257, rel=0.02)
    # the leading J is optional on either side of the match
    assert ppdot.named_pulsar_derived(pop, "0534+2200")["found"]


def test_death_line_sweep_tracks_the_constant():
    pop = ppdot.synthetic_population(n_each=400, seed=0)
    sweep = ppdot.death_line_sweep(pop["period_s"], pop["pdot"])
    fracs = [sweep[k] for k in ("0.05", "0.1", "0.2", "0.4", "1.0")]
    # a harsher death line (larger constant) can only kill more pulsars
    assert all(a >= b for a, b in zip(fracs, fracs[1:], strict=False))
    # the paper's constant reproduces population_stats exactly
    stats = ppdot.population_stats(pop["period_s"], pop["pdot"])
    assert sweep["0.2"] == pytest.approx(stats["frac_above_death"], abs=1e-4)


def test_magnetar_threshold_sweep_median_tracks_its_own_cut():
    pop = ppdot.synthetic_population(n_each=400, seed=0)
    sweep = ppdot.magnetar_threshold_sweep(pop["period_s"], pop["pdot"])
    meds, ns = [], []
    for t in (3.0e12, 5.0e12, 1.0e13, 2.0e13, 3.0e13):
        entry = sweep[str(t)]
        # the median of a B-selected class is above its own threshold by construction
        assert entry["median_log_b"] > np.log10(t)
        meds.append(entry["median_log_b"])
        ns.append(entry["n"])
    # and it rises with the cut while the class shrinks -- the referee's point
    assert all(a <= b for a, b in zip(meds, meds[1:], strict=False))
    assert all(a >= b for a, b in zip(ns, ns[1:], strict=False))


def test_msp_period_sweep_is_stable_where_the_paper_says_so():
    pop = ppdot.synthetic_population(n_each=400, seed=0)
    sweep = ppdot.msp_period_sweep(pop["period_s"], pop["pdot"])
    m20, m50 = sweep["0.02"], sweep["0.05"]
    assert m20["n_msp"] <= m50["n_msp"]  # a wider cut can only admit more
    # the fixture's clouds are far from the cut, so the medians barely move
    assert abs(m20["median_log_b_msp"] - m50["median_log_b_msp"]) < 0.3
    assert abs(m20["median_log_b_normal"] - m50["median_log_b_normal"]) < 0.3


def test_named_pulsar_derived_ambiguous_name_raises():
    pop = {
        "name": np.array(["J0534+2200", "0534+2200"], dtype=object),  # same after J-strip
        "period_s": np.array([0.0334, 0.5]),
        "pdot": np.array([4.2e-13, 1e-15]),
    }
    with pytest.raises(ValueError, match="ambiguous"):
        ppdot.named_pulsar_derived(pop, "J0534+2200")


def _real_shaped_pop():
    """A synthetic population dressed as a real fetch: names, S1400, anchors, and known
    discards, so the real-leg paths (anchors, flux cuts, CSV, discard counts) are testable."""
    base = ppdot.synthetic_population(n_each=100, seed=3)
    n = base["period_s"].size
    rng = np.random.default_rng(0)
    # literature (P, Pdot) for the four anchors, inside their gate windows
    anchors = {
        "J0534+2200": (0.0333924, 4.21e-13),
        "J1939+2134": (0.00155780, 1.051e-19),
        "J1550-5418": (2.06983, 2.318e-11),
        "J2144-3933": (8.50983, 4.96e-16),
    }
    # three rows exercising the discard breakdown: null, negative, and exactly-zero Pdot
    extra_p = np.array([a[0] for a in anchors.values()] + [1.0, 1.0, 1.0])
    extra_pd = np.array([a[1] for a in anchors.values()] + [np.nan, -1e-15, 0.0])
    p = np.concatenate([base["period_s"], extra_p])
    pd = np.concatenate([base["pdot"], extra_pd])
    names = np.asarray(
        [f"J{i:04d}+0000" for i in range(n)] + list(anchors) + ["JN1", "JN2", "JN3"],
        dtype=object,
    )
    s1400 = rng.uniform(0.05, 20.0, p.size)
    s1400[rng.random(p.size) < 0.3] = np.nan
    return {
        "period_s": p,
        "pdot": pd,
        "name": names,
        "s1400_mjy": s1400,
        "n_catalogue": int(p.size),
        "catalogue_version": "test fixture",
        "fetched_utc": "2026-08-25T00:00:00Z",
        "catalogue_max_p0_s": float(np.nanmax(p)),
    }


def test_run_real_shaped_writes_anchors_sweeps_and_catalogue(tmp_path, monkeypatch):
    monkeypatch.setattr(ppdot, "fetch_atnf_ppdot", _real_shaped_pop)
    m = ppdot.run(out=str(tmp_path), offline=False)
    # all four anchors read from the table and inside their literature windows
    assert set(m["anchors"]) == set(ppdot.ANCHOR_PSRJS)
    assert m["anchors"]["J1939+2134"]["role"] == "msp"
    assert m["anchors"]["J0534+2200"]["b_gauss"] == pytest.approx(3.8e12, rel=0.03)
    # the injected discards are counted, not assumed
    assert m["n_pdot_null"] == 1 and m["n_pdot_negative"] == 1 and m["n_pdot_zero"] == 1
    assert m["n_p0_missing"] == 0
    assert (
        m["n_pulsars"]
        + m["n_pdot_null"]
        + m["n_pdot_negative"]
        + m["n_pdot_zero"]
        + m["n_p0_missing"]
        == m["n_catalogue"]
    )
    # sweeps and flux-cut robustness are committed
    assert "0.2" in m["death_line_sweep_b12_over_p2"]
    assert "flux_cut_medians" in m and m["flux_cut_max_excursion_dex"] >= 0.0
    # the per-pulsar table exists, with one row per analysed pulsar
    csv_path = tmp_path / "results" / "ppdot_pulsars.csv"
    rows = csv_path.read_text().strip().splitlines()
    assert rows[0].startswith("psrj,")
    assert len(rows) - 1 == m["n_pulsars"]
    assert any(",magnetar," in r for r in rows[1:])  # class labels not truncated
    macros = (tmp_path / "papers" / "ppdot" / "generated" / "macros.tex").read_text()
    assert r"\ppMspAnchorB" in macros and r"\ppFracAlivePctLo" in macros
    assert r"\ppAccuracy" not in macros  # renamed away: offline-only, nothing cites it


def test_run_real_shaped_anchor_gate_catches_a_units_regression(tmp_path, monkeypatch):
    def bad_fetch():
        pop = _real_shaped_pop()
        pop["pdot"] = pop["pdot"] * 1.0e6  # a wrong-units fetch
        return pop

    monkeypatch.setattr(ppdot, "fetch_atnf_ppdot", bad_fetch)
    with pytest.raises(RuntimeError, match="outside literature window"):
        ppdot.run(out=str(tmp_path), offline=False)


def test_named_pulsar_derived_absent_or_unusable_rows():
    pop = {
        "name": np.array(["J0534+2200", "J1234+5678"], dtype=object),
        "period_s": np.array([0.0333924, 0.5]),
        "pdot": np.array(
            [4.21e-13, -1.0e-15]
        ),  # negative Pdot: cluster acceleration, not spin-down
    }
    assert not ppdot.named_pulsar_derived(pop, "J9999+9999")["found"]
    assert not ppdot.named_pulsar_derived(pop, "J1234+5678")["found"]
    # a fetch that drops the name column must not silently claim a match
    assert not ppdot.named_pulsar_derived({"period_s": [1.0], "pdot": [1e-15]}, "J0534+2200")[
        "found"
    ]
