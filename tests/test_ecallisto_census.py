"""Tests for jansky_research.ecallisto_census -- the type III occurrence census. No network."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import ecallisto_census as census


def test_coverage_corrected_rate_normalises_and_flags_zero():
    n = np.array([10.0, 0.0, 6.0])
    c = np.array([5.0, 2.0, 0.0])
    rate = census.coverage_corrected_rate(n, c)
    assert rate[0] == 2.0  # 10 / 5
    assert rate[1] == 0.0  # 0 / 2
    assert np.isnan(rate[2])  # zero coverage -> NaN, not a divide-by-zero


def test_census_correlation_recovers_a_linear_trend():
    s = np.arange(20.0)
    rate = 0.05 * s + 1.0  # exactly linear, slope 0.05
    corr = census.census_correlation(rate, s)
    assert corr["n_periods"] == 20
    assert corr["pearson_r"] > 0.999
    assert corr["spearman_rho"] > 0.999
    assert abs(corr["slope"] - 0.05) < 1e-9


def test_census_correlation_degenerate_returns_nan():
    # too few finite points -> NaN metrics, no crash
    corr = census.census_correlation(np.array([1.0, np.nan]), np.array([1.0, 2.0]))
    assert corr["n_periods"] == 1
    assert np.isnan(corr["pearson_r"])
    # zero variance in the rate -> NaN too
    flat = census.census_correlation(np.ones(10), np.arange(10.0))
    assert np.isnan(flat["pearson_r"])


def test_synthetic_sunspots_has_a_cycle_and_is_nonnegative():
    ss = census.synthetic_sunspots(n_months=180, seed=0)
    assert ss.shape == (180,)
    assert (ss >= 0).all()
    assert ss.max() > 50.0  # a real activity swing, not flat noise


def test_synthetic_census_round_trips_the_injected_slope():
    # the whole point: rate = N/C recovers the injected proportionality to the sunspot number
    sunspot = census.synthetic_sunspots(n_months=180, seed=0)
    n_events, coverage = census.synthetic_census(sunspot, k=0.03, seed=1)
    assert (coverage >= 2.0).all()
    rate = census.coverage_corrected_rate(n_events, coverage)
    corr = census.census_correlation(rate, sunspot)
    assert corr["pearson_r"] > 0.9  # strong recovered correlation
    assert abs(corr["slope"] - 0.03) < 0.01  # recovers k within noise


def test_parse_silso_reads_semicolon_columns_and_drops_missing():
    text = "\n".join(
        [
            "2020;01;2020.042;120.5;3.4;30;1",
            "2020;02;2020.125;-1.0;0.0;0;0",  # -1 sunspot = missing, dropped
            "2020;02;NaNyear;bad;x;y;z",  # non-numeric cols -> ValueError, skipped
            "garbage line",
            "2020;03;2020.208;95.0;2.1;28;1",
        ]
    )
    out = census.parse_silso(text)
    assert out["sunspot"].tolist() == [120.5, 95.0]
    assert np.allclose(out["decimal_year"], [2020.042, 2020.208])


def test_run_offline_writes_artifacts_and_recovers_correlation(tmp_path):
    metrics = census.run(str(tmp_path), offline=True)
    assert metrics["source"] == "synthetic"
    assert metrics["n_periods"] == 180
    assert metrics["n_events_total"] > 0
    assert metrics["pearson_r"] > 0.9
    assert abs(metrics["slope"] - 0.03) < 0.01

    saved = json.loads((tmp_path / "results" / "ecallisto_census_metrics.json").read_text())
    assert saved == metrics
    fig = tmp_path / "papers" / "ecallisto_census" / "figures" / "census.pdf"
    assert fig.exists() and fig.stat().st_size > 0
    macros = (tmp_path / "papers" / "ecallisto_census" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\ecsPearson}" in macros
    assert r"\newcommand{\ecsSource}{synthetic}" in macros


def test_write_macros_placeholders_for_missing_keys(tmp_path):
    # the offline/real macro union: a None metric renders as the LaTeX-safe placeholder
    path = tmp_path / "macros.tex"
    census._write_macros(
        {"source": "e-Callisto x SILSO (0 days)", "pearson_r": None, "slope": None}, path
    )
    text = path.read_text()
    assert r"\newcommand{\ecsPearson}{--}" in text
    assert r"\newcommand{\ecsSlope}{--}" in text
