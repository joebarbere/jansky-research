"""Tests for jansky_research.rfitrend -- e-Callisto megaconstellation RFI trend. Offline."""

from __future__ import annotations

import json

import numpy as np
import pytest

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
    # a grid that skips the low-line region keeps 150 + 175, drops the low lines
    grid = np.concatenate([np.linspace(112, 123, 30), np.linspace(141, 200, 120)])
    lines = rf.available_lines(grid)
    assert 125.0 not in lines and 135.0 not in lines  # the sparse grid genuinely drops them
    assert 150.0 in lines and 175.0 in lines
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


def test_flank_contamination_biases_the_line_slope_and_can_flip_its_sign():
    """The validation arm that can fail in the regime that limits the study.

    The clean arms inject only what the line-vs-adjacent difference cancels algebraically
    (common-mode gain, broadband bursts), so they could not fail. Local RFI accruing in the
    flanking channels only is the systematic the paper names as decisive for the HUMAIN/ALMATY
    sign disagreement, and the difference cannot cancel it: the recovered slope must be biased
    low, and with a large enough flank ramp must come out FALLING under a rising injected line.
    """
    clean = rf.synthetic_month_stack()
    dirty = rf.synthetic_month_stack(flank_rise=6.0)
    line_c = np.array([rf.line_vs_adjacent(m["data"], m["freqs"]) for m in clean["months"]])
    line_d = np.array([rf.line_vs_adjacent(m["data"], m["freqs"]) for m in dirty["months"]])
    tr_c = rf.trend_fit(clean["years"], line_c)
    tr_d = rf.trend_fit(dirty["years"], line_d)
    assert tr_c["slope"] > 0  # the clean arm recovers the injected rise
    assert tr_d["slope"] < tr_c["slope"]  # flank RFI biases the slope low...
    assert tr_d["slope"] < 0  # ...and here flips its sign: the ALMATY mechanism end-to-end


def test_synthetic_metrics_report_the_flank_bias():
    m = rf._synthetic_metrics()
    assert m["flank_sign_flipped"] is True
    assert m["flank_slope_bias_per_yr"] < 0
    # the bias is the difference of the two committed slopes, not an independent number
    assert m["flank_slope_bias_per_yr"] == pytest.approx(
        m["flank_line_slope_per_yr"] - m["line_excess_slope_per_yr"], abs=1e-3
    )


def test_trend_fit_carries_its_interval():
    rng = np.random.default_rng(0)
    x = np.linspace(2013, 2026, 100)
    y = 0.5 * x + rng.normal(0, 1.0, 100)
    fit = rf.trend_fit(x, y)
    assert fit["slope_lo"] < fit["slope"] < fit["slope_hi"]
    assert fit["slope_lo"] < 0.5 < fit["slope_hi"]


def test_step_analysis_tells_a_staircase_from_a_trend():
    yr = 2013 + np.arange(160) / 12.0
    stair = np.where(yr < 2018.9, 0.0, -12.0) + np.where(yr < 2022.2, 0.0, 15.0)
    rng = np.random.default_rng(1)
    stair = stair + rng.normal(0, 0.5, yr.size)
    d = rf.step_analysis(yr, stair)
    ats = sorted(s["at_year"] for s in d["steps"])
    assert abs(ats[0] - 2018.9) < 0.15 and abs(ats[1] - 2022.2) < 0.15
    # within regimes there is no significant rise
    for seg in d["segments"]:
        if seg["p_value"] is not None:
            assert seg["p_value"] > 0.01 or seg["slope"] <= 0.2
    # a genuine trend, by contrast, rises within its (arbitrary) segments
    trend = 1.0 * (yr - yr.min()) + rng.normal(0, 0.3, yr.size)
    dt = rf.step_analysis(yr, trend)
    rising = [s for s in dt["segments"] if s["slope"] is not None and s["slope"] > 0.5]
    assert len(rising) >= 2


def _fake_per_station(slope_a=0.0, slope_b=0.0, n=150, seed=0):
    rng = np.random.default_rng(seed)
    yr = 2013 + np.arange(n) / 12.0
    out = {}
    for name, sl in (("A", slope_a), ("B", slope_b)):
        lx = sl * (yr - yr.min()) + rng.normal(0, 1.0, n)
        out[name] = {
            "n_months": n,
            "stable_lines": [150.0, 175.0],
            "years": yr.tolist(),
            "line_excess": lx.tolist(),
            "line_excess_slope_per_yr": sl,
            "line_excess_p": 0.5,
        }
    return out


def test_coherence_power_rises_with_injected_amplitude():
    ps = _fake_per_station()
    curve = rf.coherence_power(ps, amps=(0.0, 20.0), n_trials=40)
    p0 = next(r["power"] for r in curve if r["amp"] == 0.0)
    p20 = next(r["power"] for r in curve if r["amp"] == 20.0)
    assert p0 < 0.2  # near the criterion's false-positive rate on null series
    assert p20 > 0.8  # a huge common signal is detected


def test_equal_n_pooled_reports_quantiles():
    ps = _fake_per_station(slope_a=1.0, slope_b=1.0)
    ps["A"]["years"] = ps["A"]["years"][:100]
    ps["A"]["line_excess"] = ps["A"]["line_excess"][:100]
    eq = rf.equal_n_pooled(ps, n_draws=50)
    assert eq["n_min"] == 100
    assert eq["pooled_slope_p05"] <= eq["pooled_slope_median"] <= eq["pooled_slope_p95"]
    assert eq["pooled_slope_median"] > 0.5  # a genuine common trend survives equal-n


def test_regressor_shape_comparison_prefers_the_generating_shape():
    yr = 2013 + np.arange(160) / 12.0
    star_shaped = rf.starlink_count(yr) / rf.starlink_count(2026.5)
    rng = np.random.default_rng(2)
    v = 5.0 * star_shaped + rng.normal(0, 0.3, yr.size)
    sh = rf.regressor_shape_comparison(yr, v)
    assert sh["corr_starlink_shape"] > sh["corr_linear_ramp"]
    assert sh["corr_starlink_shape"] > 0.9


def test_window_matched_trend_restricts_the_span():
    ps = _fake_per_station(slope_a=1.0, slope_b=0.0)
    # B starts later: truncate its first 30 months
    ps["B"]["years"] = ps["B"]["years"][30:]
    ps["B"]["line_excess"] = ps["B"]["line_excess"][30:]
    w = rf.window_matched_trend(ps, "A", "B")
    assert w["n"] == 120
    assert w["window"][0] > 2015.0
    assert abs(w["slope"] - 1.0) < 0.3  # the trend survives the window restriction


def test_synthetic_metrics_truths_and_sweep():
    m = rf._synthetic_metrics()
    # the recover-a-known now states the known and scores an amplitude ratio
    assert abs(m["line_recovery_ratio"] - 1.0) < 0.1
    assert abs(m["diff_recovery_ratio"] - 1.0) < 0.15
    # burst immunity measured on the primary metric
    assert abs(m["line_slope_burst_none"] - m["line_slope_burst_heavy"]) < 0.03
    # the flank sweep crosses zero between 3 and 6, and flank_rise=5 is in it
    sweep = {row["flank_rise"]: row["slope"] for row in m["flank_sweep"]}
    assert sweep[0.0] > 0 and sweep[6.0] < 0
    assert any(r["flank_rise"] == 5.0 for r in m["flank_sweep"])


def test_derived_real_analyses_from_committed_shapes():
    ps = _fake_per_station(slope_a=0.5, slope_b=-0.1)
    d = rf.derived_real_analyses(ps)
    assert set(d["stations"]) == {"A", "B"}
    assert "slope_ci95" in d["stations"]["A"]
    assert d["coherence_power"] and d["equal_n_pooled"]
    assert d["window_matched"]["n"] > 0
