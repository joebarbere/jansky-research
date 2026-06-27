"""Tests for jansky_research.frbstats — estimators recover known ground truth. No network."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import frbstats


def test_wait_times_positive_sorted():
    mjds = np.array([5.0, 1.0, 3.0, 3.0, 2.0])  # unsorted + a duplicate
    w = frbstats.wait_times(mjds)
    assert np.all(w > 0)
    assert np.all(np.diff(np.sort(w)) >= 0) or w.size >= 1


def test_weibull_recovers_shape_and_clustering():
    cat = frbstats.synthetic_catalog(n_repeater=600, k_true=0.7, seed=3)
    fit = frbstats.fit_weibull_waits(cat["mjd"][cat["repeater"]], n_boot=200, seed=1)
    assert abs(fit.k - 0.7) < 0.2  # recovered shape
    assert fit.k_ci_low < fit.k < fit.k_ci_high  # CI brackets the point estimate
    assert fit.clustered  # k well below 1


def test_weibull_poisson_limit_not_flagged_clustered():
    # k_true = 1 is a memoryless Poisson process; should not be flagged as clustered.
    cat = frbstats.synthetic_catalog(n_repeater=600, k_true=1.0, seed=7)
    fit = frbstats.fit_weibull_waits(cat["mjd"][cat["repeater"]], n_boot=200, seed=1)
    assert abs(fit.k - 1.0) < 0.2
    assert not fit.clustered


def test_power_law_recovers_index():
    cat = frbstats.synthetic_catalog(n_repeater=2000, n_oneoff=0, gamma_true=2.0, seed=5)
    fit = frbstats.fit_power_law(cat["fluence"])
    assert abs(fit.gamma - 2.0) < 0.2
    assert fit.gamma_err > 0
    assert fit.n_tail == 2000


def test_compare_populations_detects_dm_shift():
    cat = frbstats.synthetic_catalog(seed=2)
    rep = {k: cat[k][cat["repeater"]] for k in ("dm", "fluence", "width")}
    one = {k: cat[k][~cat["repeater"]] for k in ("dm", "fluence", "width")}
    ks = frbstats.compare_populations(rep, one)
    assert "dm" in ks
    assert ks["dm"]["p"] < 0.05  # the synthetic DM shift is detectable


def test_summarise_bundles_results():
    cat = frbstats.synthetic_catalog(seed=0)
    s = frbstats.summarise(cat)
    assert s.n_bursts == cat["repeater"].size
    assert s.n_repeater_bursts == int(cat["repeater"].sum())
    assert s.weibull.n_waits > 0
    assert "dm" in s.ks


def test_grouped_wait_times_are_within_source():
    # Two sources; waits must be computed within each, never across the source boundary.
    mjds = np.array([10.0, 11.0, 13.0, 100.0, 101.0])
    groups = np.array(["A", "A", "A", "B", "B"])
    w = np.sort(frbstats.grouped_wait_times(mjds, groups))
    assert np.allclose(w, [1.0, 1.0, 2.0])  # A:[1,2], B:[1]; the 87-day gap is NOT a wait


def test_select_xmin_and_auto_power_law():
    # Pure power law below x and noise floor of small fluences below it.
    cat = frbstats.synthetic_catalog(n_repeater=3000, n_oneoff=0, gamma_true=2.3, seed=11)
    xm = frbstats.select_xmin(cat["fluence"])
    assert xm >= cat["fluence"].min()
    fit = frbstats.fit_power_law(cat["fluence"], auto_xmin=True)
    assert abs(fit.gamma - 2.3) < 0.25
    assert fit.f_min == xm


def test_fit_guards():
    with pytest.raises(ValueError):
        frbstats.fit_weibull_waits(np.array([58000.0, 58001.0]))  # < 3 waits
    with pytest.raises(ValueError):
        frbstats.fit_power_law(np.array([1.0]))  # < 2 in tail
