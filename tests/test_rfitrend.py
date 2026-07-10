"""Tests for jansky_research.rfitrend -- e-Callisto megaconstellation RFI trend. Offline."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import rfitrend as rf


def test_occupancy_metric_band_select_and_burst_immune():
    freqs = np.linspace(45, 870, 200)
    rng = np.random.default_rng(0)
    data = rng.normal(10.0, 0.5, (200, 120))
    # raise the UEM band's persistent level
    data[(freqs >= 110) & (freqs <= 170)] += 5.0
    base = rf.occupancy_metric(data, freqs, rf.UEM_BAND_MHZ)
    assert base > rf.occupancy_metric(data, freqs, rf.FM_CONTROL_MHZ)  # UEM band is elevated
    # a transient burst (few time columns) must NOT move the median-based level
    burst = data.copy()
    burst[:, 10:18] += 30.0
    assert abs(rf.occupancy_metric(burst, freqs, rf.UEM_BAND_MHZ) - base) < 0.3
    # out-of-coverage band -> NaN
    assert np.isnan(rf.occupancy_metric(data, freqs, (900.0, 1000.0)))


def test_band_differential_cancels_gain_and_bursts():
    freqs = np.linspace(45, 870, 200)
    rng = np.random.default_rng(1)
    base = rng.normal(10.0, 0.5, (200, 120))
    d0 = rf.band_differential(base, freqs)
    # add a common-mode gain offset (whole spectrum up) -> differential unchanged
    d_gain = rf.band_differential(base + 7.0, freqs)
    assert abs(d_gain - d0) < 0.1
    # add a broadband burst (all channels, few columns) -> differential unchanged
    b = base.copy()
    b[:, 20:28] += 25.0
    assert abs(rf.band_differential(b, freqs) - d0) < 0.3


def test_line_vs_adjacent_detects_narrowband_excess():
    freqs = np.linspace(45, 870, 400)
    rng = np.random.default_rng(2)
    data = rng.normal(10.0, 0.3, (400, 120))
    for line in rf.UEM_LINES_MHZ:  # inject a narrowband excess at each UEM line
        data[np.abs(freqs - line) <= 0.8] += 6.0
    excess = rf.line_vs_adjacent(data, freqs)
    assert excess > 3.0  # lines sit well above their flanks
    # no injected lines -> excess ~ 0
    assert abs(rf.line_vs_adjacent(rng.normal(10, 0.3, (400, 120)), freqs)) < 1.0


def test_pick_control_band_prefers_fm_but_adapts_to_notches():
    full = np.linspace(45, 870, 300)  # samples FM -> prefers it
    assert rf.pick_control_band(full)[0] == "FM"
    # a HUMAIN-like grid that notches out FM (gap 84->112) must fall back to a sampled band
    notched = np.concatenate([np.linspace(45, 84, 60), np.linspace(112, 437, 140)])
    name, band = rf.pick_control_band(notched)
    assert name in ("low", "high") and ((notched >= band[0]) & (notched <= band[1])).sum() >= 4
    # no clean band sampled at all -> "none"
    assert rf.pick_control_band(np.linspace(120, 160, 40))[0] == "none"


def test_available_lines_drops_notched_lines():
    # a grid that skips the 137 MHz region (like HUMAIN) keeps 150 + 175, drops 137
    grid = np.concatenate([np.linspace(112, 123, 30), np.linspace(141, 200, 120)])
    lines = rf.available_lines(grid)
    assert 137.05 not in lines and 150.0 in lines and 175.0 in lines
    # a full grid keeps all three
    assert set(rf.available_lines(np.linspace(45, 870, 800))) == set(rf.UEM_LINES_MHZ)


def test_band_differential_auto_control_on_notched_grid():
    # FM-notched grid: control_band=None must auto-pick a sampled control and return finite
    grid = np.concatenate([np.linspace(45, 84, 80), np.linspace(112, 437, 200)])
    rng = np.random.default_rng(7)
    data = rng.normal(10.0, 0.5, (grid.size, 100))
    data[(grid >= 110) & (grid <= 170)] += 4.0
    d = rf.band_differential(data, grid, control_band=None)
    assert np.isfinite(d) and d > 0  # UEM elevated over the auto-picked control


def test_trend_fit_robust_and_significant():
    x = np.linspace(2012, 2026, 60)
    rng = np.random.default_rng(3)
    y = 0.2 * (x - 2012) + rng.normal(0, 0.3, 60)
    y[10] += 50.0  # an outlier month -- Theil-Sen must shrug it off
    tr = rf.trend_fit(x, y)
    assert abs(tr["slope"] - 0.2) < 0.05  # robust slope recovers the injected 0.2/yr
    assert tr["p_value"] < 1e-6  # clearly significant
    # a flat series -> not significant
    assert rf.trend_fit(x, rng.normal(0, 1, 60))["p_value"] > 0.05


def test_starlink_count_monotone_and_zero_pre_2019():
    assert rf.starlink_count(2015.0) == 0.0
    assert rf.starlink_count(2018.5) == 0.0
    c = rf.starlink_count(np.array([2019.0, 2022.0, 2026.0]))
    assert c[0] < c[1] < c[2] and c[2] > 8000  # grows through the megaconstellation era


def test_summarize_stations_coherence_verdict():
    yrs = list(np.linspace(2019, 2026, 40))
    rise = list(np.linspace(0, 5, 40))
    fall = list(np.linspace(5, 0, 40))
    # INCOHERENT: two significant stations that DISAGREE in sign (the real HUMAIN vs ALMATY case)
    incoh = {
        "A": {
            "n_months": 40,
            "stable_lines": [150.0],
            "line_excess_slope_per_yr": 0.7,
            "line_excess_p": 1e-4,
            "years": yrs,
            "line_excess": rise,
        },
        "B": {
            "n_months": 40,
            "stable_lines": [137.05],
            "line_excess_slope_per_yr": -0.7,
            "line_excess_p": 1e-4,
            "years": yrs,
            "line_excess": fall,
        },
        "C": {"n_months": 30, "stable_lines": [], "years": [], "line_excess": []},  # no lines
    }
    r = rf.summarize_stations(incoh)
    assert r["n_significant_stations"] == 2 and r["n_rising"] == 1 and r["n_falling"] == 1
    assert not r["cross_station_signs_agree"] and not r["coherent_rise"]
    assert r["n_stations_with_lines"] == 2
    # COHERENT: two significant stations that BOTH rise with the Starlink count
    coh = {
        "A": {
            "n_months": 40,
            "stable_lines": [150.0],
            "line_excess_slope_per_yr": 0.7,
            "line_excess_p": 1e-4,
            "years": yrs,
            "line_excess": rise,
        },
        "B": {
            "n_months": 40,
            "stable_lines": [175.0],
            "line_excess_slope_per_yr": 0.6,
            "line_excess_p": 1e-3,
            "years": yrs,
            "line_excess": rise,
        },
    }
    r2 = rf.summarize_stations(coh)
    assert r2["cross_station_signs_agree"] and r2["n_rising"] == 2
    assert r2["coherent_rise"] and r2["corr_with_starlink"] > 0.3


def test_run_offline_recovers_trend_and_control_flat(tmp_path):
    m = rf.run(str(tmp_path), offline=True)
    assert m["recovered_uem_trend"]  # differential recovers the injected UEM trend
    assert m["diff_trend_p"] < 0.01 and m["corr_with_starlink"] > 0.9
    # the PRIMARY metric (narrowband line-vs-adjacent excess) also recovers the injected trend
    assert m["recovered_line_trend"] and m["line_corr_with_starlink"] > 0.9
    assert m["control_flat"]  # the null control (two clean bands) does NOT trend
    saved = json.loads((tmp_path / "results" / "rfitrend_metrics.json").read_text())
    assert saved["n_months"] == m["n_months"]
    assert (tmp_path / "papers" / "rfitrend" / "figures" / "rfitrend.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "rfitrend" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\rfRealSlope}{--}" in macros  # real namespace placeholder offline


def test_synthetic_burst_immunity():
    # the recovered trend is the same whether or not months carry broadband solar bursts
    s0 = rf.synthetic_month_stack(burst_frac=0.0, seed=5)
    s1 = rf.synthetic_month_stack(burst_frac=0.9, seed=5)
    d0 = np.array([rf.band_differential(x["data"], x["freqs"]) for x in s0["months"]])
    d1 = np.array([rf.band_differential(x["data"], x["freqs"]) for x in s1["months"]])
    assert (
        abs(rf.trend_fit(s0["years"], d0)["slope"] - rf.trend_fit(s1["years"], d1)["slope"]) < 0.03
    )


def test_write_macros_dual_namespace(tmp_path):
    p = tmp_path / "m.tex"
    rf._write_macros(
        {
            "source": "x",
            "is_real": True,
            "n_months": 210,
            "n_stations": 3,
            "line_excess_slope_per_yr": 0.09,
            "line_excess_trend_p": 2e-6,
        },
        p,
    )
    txt = p.read_text()
    # the real headline slope IS the line-excess (primary metric)
    assert r"\newcommand{\rfRealSlope}{0.09}" in txt
    assert r"\newcommand{\rfRealLineExcessSlope}{0.09}" in txt
    assert r"\newcommand{\rfRealNStations}{3}" in txt
    assert r"\newcommand{\rfRealTrendP}{<10^{-5}}" in txt  # p<1e-5 -> upper-bound math body
    # rfSyn* is ALWAYS live (recomputed synthetic recovery), not a placeholder
    assert r"\newcommand{\rfSynSlope}{--}" not in txt
    assert r"\newcommand{\rfSynNMonths}{168}" in txt


def test_write_macros_offline_leaves_real_placeholders(tmp_path):
    p = tmp_path / "m.tex"
    rf._write_macros(rf._synthetic_metrics(), p)  # is_real=False
    txt = p.read_text()
    assert r"\newcommand{\rfRealSlope}{--}" in txt  # no real run -> placeholder
    assert r"\newcommand{\rfSynSlope}{--}" not in txt  # synthetic always live
