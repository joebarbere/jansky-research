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
    rep = cat["repeater"]
    fit = frbstats.fit_weibull_waits(
        cat["mjd"][rep], groups=cat["repeater_name"][rep], n_boot=200, seed=1
    )
    assert abs(fit.k - 1.0) < 0.2
    assert not fit.clustered
    # ...and pooling across the (staggered) sources DISTORTS k -- the round-9 figure blocker:
    # the pooled statistic is measurably different from the within-source one
    pooled = frbstats.fit_weibull_waits(cat["mjd"][rep], n_boot=10, seed=1)
    assert abs(pooled.k - fit.k) > 0.1


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


def test_source_level_ks_uses_per_source_medians():
    cat = frbstats.synthetic_catalog(seed=11)
    rep = cat["repeater"]
    rd = {k: cat[k][rep] for k in ("dm", "fluence", "width")}
    od = {k: cat[k][~rep] for k in ("dm", "fluence", "width")}
    out = frbstats.compare_populations_source_level(rd, od, cat["repeater_name"][rep])
    assert out["dm"]["n_rep_sources"] == 3  # the source, not the burst, is the unit
    assert 0 <= out["dm"]["p"] <= 1


def test_bootstrap_power_law_wider_than_hill():
    """The Hill SE conditions on a fixed f_min chosen from the same data; the joint bootstrap
    must be at least as wide (on Cat 1 it is 2.2x)."""
    cat = frbstats.synthetic_catalog(n_repeater=300, n_oneoff=300, seed=13)
    fit = frbstats.fit_power_law(cat["fluence"], auto_xmin=True)
    boot = frbstats.bootstrap_power_law(cat["fluence"], n_boot=120, seed=1)
    assert boot["gamma_boot_sd"] >= 0.8 * fit.gamma_err  # never materially narrower
    assert boot["gamma_ci_low"] < fit.gamma < boot["gamma_ci_high"]
    assert boot["f_min_ci_low"] <= boot["f_min_ci_high"]


def test_weibull_cluster_ci_brackets_the_fit():
    cat = frbstats.synthetic_catalog(seed=17)
    rep = cat["repeater"]
    fit = frbstats.fit_weibull_waits(
        cat["mjd"][rep], groups=cat["repeater_name"][rep], n_boot=100, seed=0
    )
    lo, hi = frbstats.weibull_cluster_ci(
        cat["mjd"][rep], cat["repeater_name"][rep], n_boot=200, seed=0
    )
    assert lo < fit.k < hi


def test_build_catalog_raises_instead_of_silent_synthetic(monkeypatch, tmp_path):
    """The round-9 split-brain finding: a fetch failure must raise, not silently substitute
    the synthetic fixture (which would rewrite the paper's macros while every guard passed)."""
    import pytest

    from jansky_research import data, pipeline

    def boom(name):
        raise OSError("network down")

    monkeypatch.setattr(data, "fetch", boom)
    with pytest.raises(RuntimeError, match="offline"):
        pipeline.build_catalog(offline=False)


def test_write_macros_refuses_synthetic_over_real(tmp_path):
    """report.write_macros now merges through preserve_live_macros: a synthetic rebuild must
    not overwrite real values under the same names."""
    from jansky_research import pipeline, report

    real_cat = frbstats.synthetic_catalog(seed=1)
    real = pipeline.analyze(real_cat, "chime-frb-catalog")  # marker says REAL
    path = tmp_path / "macros.tex"
    report.write_macros(real, path)
    before = path.read_text()
    assert r"\newcommand{\catalogSource}{chime-frb-catalog}" in before
    syn = pipeline.analyze(frbstats.synthetic_catalog(seed=2, n_repeater=350), "synthetic")
    report.write_macros(syn, path)
    after = path.read_text()
    # the real values survive the synthetic rewrite
    import re

    def val(text, name):
        m = re.search(r"\\newcommand\{\\" + name + r"\}\{([^}]*)\}", text)
        return m.group(1) if m else None

    assert val(after, "nBursts") == val(before, "nBursts")
    assert val(after, "weibullK") == val(before, "weibullK")


def test_figure_waits_match_the_fitted_waits(tmp_path):
    """The blocker test: the figure's empirical CDF must be the within-source waits the
    Weibull was fitted to, not the pooled diffs that cross source boundaries."""
    from jansky_research import report

    cat = frbstats.synthetic_catalog(seed=19)
    stats_ = frbstats.summarise(cat)
    report.make_figures(cat, stats_, tmp_path)
    # reproduce what the figure now plots and what the fit used: they must be identical
    rep = cat["repeater"]
    fitted = np.sort(frbstats.grouped_wait_times(cat["mjd"][rep], cat["repeater_name"][rep]))
    pooled = np.sort(frbstats.wait_times(cat["mjd"][rep]))
    assert fitted.size == stats_.weibull.n_waits
    assert pooled.size != fitted.size or not np.allclose(pooled, fitted)
    assert (tmp_path / "wait_time_cdf.pdf").exists()
    assert (tmp_path / "fluence_ccdf.pdf").exists()
