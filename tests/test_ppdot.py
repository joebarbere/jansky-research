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
