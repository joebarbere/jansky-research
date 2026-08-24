"""Tests for jansky_research.skr -- Cassini SKR occurrence + proximity census. Offline."""

from __future__ import annotations

import json

import numpy as np

from jansky_research import skr


def test_scet_to_jd_known_dates():
    # 2013 day-of-year 293 = 2013-10-20; JD at 00:00 UT = 2456585.5
    assert skr.parse_scet_to_jd("2013-293T00:00:00.000") == 2456585.5
    assert abs(skr.parse_scet_to_jd("2013-293T12:00:00.000") - 2456586.0) < 1e-9
    assert abs(skr.parse_scet_to_jd("2004-001T00:00:00.000") - 2453005.5) < 1e-6


def test_band_integrated_flux_selects_band():
    # log-spaced electric grid; only channels in [1e5, 1.2e6] should contribute
    freqs = np.array([1e3, 1e4, 1e5, 3e5, 1e6, 1e7])
    dens = np.array([9.0, 9.0, 1.0, 1.0, 1.0, 9.0])  # in-band all 1.0
    val = skr.band_integrated_flux(freqs, dens, (1e5, 1.2e6))
    # trapezoid of constant 1.0 over [1e5, 1e6] ~= 9e5
    assert abs(val - 9e5) / 9e5 < 0.01
    # fewer than 2 in-band channels -> NaN
    assert np.isnan(skr.band_integrated_flux(np.array([1e3, 1e7]), np.array([1.0, 1.0])))


def _synthetic_tab(tmp_path, n_rows=20):
    """A format-faithful KEY60S .TAB: 1 frequency row + n_rows spectral-density rows."""
    # 115-channel grid: first 73 electric (1 Hz..1e7.2), then 42 magnetic (restart at 1 Hz)
    elec = [10 ** (k * 0.1) for k in range(73)]
    mag = [10 ** (k * 0.1) for k in range(42)]
    freqs = elec + mag

    def items(vals):
        return "".join(f"{v:10.3e}" for v in vals)

    lines = [f"{'2013-293T00:00:00.000':<21} 0" + items(freqs)]
    for i in range(n_rows):
        scet = f"2013-293T{i // 60:02d}:{i % 60:02d}:00.000"
        # inject strong SKR-band (electric ch 50-60) power in even rows
        dens = [1e-15] * 73
        if i % 2 == 0:
            for ch in range(50, 61):
                dens[ch] = 1e-12
        dqf = "0"
        lines.append(f"{scet:<21} {dqf}" + items(dens))
    p = tmp_path / "RPWS_KEY__2013293_3.TAB"
    p.write_text("\n".join(lines) + "\n")
    return p


def test_read_key_params_parses_real_format(tmp_path):
    d = skr.read_key_params(_synthetic_tab(tmp_path, n_rows=20))
    assert d["jd"].size == 20  # all DQF=0
    assert d["freqs"].size == 73  # electric channels only
    assert d["jd"][0] == 2456585.5
    # even rows have injected SKR power -> higher band flux than odd rows
    even = d["flux"][::2]
    odd = d["flux"][1::2]
    assert np.nanmedian(even) > 100 * np.nanmedian(odd)


def test_read_key_params_drops_bad_quality(tmp_path):
    p = _synthetic_tab(tmp_path, n_rows=10)
    lines = p.read_text().splitlines()
    lines[3] = lines[3][:22] + "9" + lines[3][23:]  # mark row 3 DQF=9
    p.write_text("\n".join(lines) + "\n")
    d = skr.read_key_params(p)
    assert d["jd"].size == 9  # one dropped


def test_detect_skr_thresholds_on_floor():
    rng = np.random.default_rng(0)
    flux = 10.0 ** rng.normal(-14.0, 0.1, 500)
    flux[::10] *= 1e3  # 50 strong bins
    active = skr.detect_skr(flux, k=3.0)
    assert 30 < active.sum() < 80  # the injected strong bins, not the floor
    assert not active[1] and active[0]  # bin 0 boosted, bin 1 floor


def test_dual_period_ls_recovers_injected_period():
    s = skr.synthetic_skr(seed=1)
    ls = skr.dual_period_ls(s["jd"], s["flux"])
    # the dominant rotation-band period should be near the injected 10.6-10.7 h
    assert 10.4 <= ls["best_period_hr"] <= 10.9
    assert ls["best_power"] > 0


def test_proximity_duty_cycle_recovers_near_far_trend():
    s = skr.synthetic_skr(seed=2, near_far_contrast=6.0)
    active = s["active_true"]
    prox = skr.proximity_duty_cycle(active, s["range_rs"])
    assert len(prox["duty_by_bin"]) == 4
    # nearest bin duty cycle exceeds farthest (proximity trend injected)
    assert prox["duty_by_bin"][0] > prox["duty_by_bin"][-1]
    assert prox["near_far_ratio"] > 1.5


def test_distance_correction_flattens_pure_sensitivity_trend():
    # a range-INDEPENDENT active population seen through 1/r^2 sensitivity: raw occurrence rises
    # near-in (the same emission clears the threshold more often); distance-correction removes it,
    # so the null-corrected near/far trend collapses toward flat
    rng = np.random.default_rng(4)
    n = 6000
    range_rs = np.linspace(5, 25, n)
    rng.shuffle(range_rs)
    active_true = rng.random(n) < 0.25  # 25% intrinsically active, range-independent
    base = 10.0 ** rng.normal(-14.0, 0.1, n)
    flux = base * (1 + active_true * 300.0) / range_rs**2  # observed through 1/r^2
    raw = skr.proximity_duty_cycle(skr.detect_skr(flux), range_rs)
    corr = skr.proximity_duty_cycle(
        skr.detect_skr(skr.distance_correct_flux(flux, range_rs)), range_rs
    )
    assert raw["near_far_ratio"] > 1.5  # sensitivity manufactures a proximity trend
    # the corrected trend is much closer to flat (1.0) than the raw trend
    assert abs(corr["near_far_ratio"] - 1.0) < abs(raw["near_far_ratio"] - 1.0)


def test_latitude_by_range_bin_reports_span():
    s = skr.synthetic_skr(seed=5)
    lb = skr.latitude_by_range_bin(s["range_rs"], s["sub_lat_deg"])
    assert len(lb["abs_lat_median_by_bin"]) == 4
    assert lb["abs_lat_span_deg"] >= 0


def test_magnetic_latitude_weight_returns_ratio():
    s = skr.synthetic_skr(seed=3)
    latw = skr.magnetic_latitude_weight(s["active_true"], s["range_rs"], s["sub_lat_deg"])
    assert len(latw["weighted_duty_by_bin"]) == 4
    assert latw["weighted_near_far_ratio"] is None or latw["weighted_near_far_ratio"] > 0


def test_run_offline_recovers_and_writes(tmp_path):
    m = skr.run(str(tmp_path), offline=True)
    assert m["source"].startswith("synthetic")
    assert m["anchor_dev_pct"] < 2.0  # ~10.7 h recovered within 2% of the published periods
    assert m["near_far_ratio"] > 1.5  # proximity trend recovered
    assert "sensitivity_corrected_near_far" in m and "abs_lat_span_deg" in m
    assert m["n_active"] > 0
    saved = json.loads((tmp_path / "results" / "skr_metrics.json").read_text())
    assert saved["n_bins"] == m["n_bins"]
    assert (tmp_path / "papers" / "skr" / "figures" / "skr.pdf").stat().st_size > 0
    macros = (tmp_path / "papers" / "skr" / "generated" / "macros.tex").read_text()
    assert r"\newcommand{\skrRealNearFar}{--}" in macros  # real namespace placeholder offline


def test_write_macros_dual_namespace(tmp_path):
    p = tmp_path / "m.tex"
    skr._write_macros(
        {
            "source": "x",
            "is_real": True,
            "n_bins": 100,
            "duty_cycle_pct": 3.0,
            "anchor_period_hr": 10.7,
            "near_far_ratio": 5.0,
            "weighted_near_far_ratio": 3.0,
            "ls_fap": 0.001,
        },
        p,
    )
    txt = p.read_text()
    assert r"\newcommand{\skrRealNearFar}{5.0}" in txt
    assert r"\newcommand{\skrSynNearFar}{--}" in txt


def test_excess_correction_flattens_pure_sensitivity_with_range_independent_floor():
    """The control the old test lacked: the floor does NOT scale with range.

    The old null (`distance_correct_flux`) rescaled the total flux, floor included; its only
    test built `flux = base/r**2` so the correction exactly inverted the injection and could
    fail only by making things worse. Here the floor is range-independent (as the real
    instrumental+galactic floor is), the signal alone falls as 1/r^2, and the intrinsic
    emission rate is flat with range -- so a correct null must return the corrected ratio to
    ~1, and the OLD flux-scaled null must NOT (it drags the floor with the signal).
    """
    rng = np.random.default_rng(7)
    n = 8 * 1440
    t_hr = np.arange(n) / 60.0
    range_rs = 31.5 - 28.5 * np.cos(2 * np.pi * t_hr / (8 * 24.0))
    jd = 2.456e6 + np.arange(n) / 1440.0
    floor = 10.0 ** rng.normal(-14.0, 0.1, n)  # range-INDEPENDENT floor
    active_true = rng.random(n) < 0.25  # intrinsic rate flat with range
    # Per-event luminosities span the detection threshold across the range sweep: at periapsis
    # most events are visible, at apoapsis only the bright tail -- the realistic case where a
    # 1/r^2 sensitivity trend actually exists for a threshold detector.
    amp = 10.0 ** rng.normal(-11.8, 0.7, n)
    flux = floor + amp * active_true / range_rs**2

    raw = skr.proximity_duty_cycle(skr.detect_skr(flux), range_rs)["near_far_ratio"]
    com = skr.common_sensitivity_active(flux, range_rs)
    corr = skr.proximity_duty_cycle(com, range_rs)["near_far_ratio"]
    old = skr.proximity_duty_cycle(
        skr.detect_skr(skr.distance_correct_flux(flux, range_rs)), range_rs
    )["near_far_ratio"]
    up = skr.proximity_duty_cycle(
        skr.detect_skr(skr.distance_correct_excess(flux, range_rs)), range_rs
    )["near_far_ratio"]
    assert raw > 2.0  # the sensitivity trend is there to remove
    # the common census flattens a flat intrinsic rate, with no noise promoted
    assert abs(corr - 1.0) < 0.4
    assert (com & ~active_true).sum() == 0
    # and both superseded null models are measurably worse on the same control: the total-flux
    # rescale collapses (or empties a bin entirely); the excess-up rescale reverses the trend
    assert old is None or abs(corr - 1.0) < abs(old - 1.0)
    assert up is None or abs(corr - 1.0) < abs(up - 1.0)
    del jd


def test_common_census_recovers_a_real_proximity_trend():
    """The estimator must not flatten everything: an injected intrinsic trend survives.

    The injected rate is 0.25*(1 + 2*(1/r)/(1/3)); between the NEAR-BIN and FAR-BIN medians
    (~7.8 and ~55 R_S on this sweep) that is an intrinsic bin contrast of ~1.6, which is what
    the census should return -- the endpoint contrast (2.7 at r=3 vs 60) is not the estimand.
    """
    rng = np.random.default_rng(7)
    n = 32 * 1440
    t_hr = np.arange(n) / 60.0
    range_rs = 31.5 - 28.5 * np.cos(2 * np.pi * t_hr / (8 * 24.0))
    floor = 10.0 ** rng.normal(-14.0, 0.1, n)
    p_near = 1.0 / range_rs
    active_true = rng.random(n) < 0.25 * (1 + 2 * p_near / p_near.max())
    flux = floor + 10.0 ** rng.normal(-11.8, 0.7, n) * active_true / range_rs**2
    corr = skr.proximity_duty_cycle(skr.common_sensitivity_active(flux, range_rs), range_rs)[
        "near_far_ratio"
    ]
    assert corr > 1.3  # the intrinsic bin contrast (~1.6) is recovered, not flattened


def test_orbit_segments_and_jackknife_report_orbit_level_scatter():
    s = skr.synthetic_skr(n_days=16.0)
    segs = skr.orbit_segments(s["jd"], s["range_rs"])
    assert 2 <= len(segs) <= 5  # ~2 apoapses in 16 days of an 8-day sweep -> 2-3 segments
    jk = skr.jackknife_near_far(s["jd"], s["flux"], s["range_rs"], mode="raw")
    assert jk["n_orbits"] == len(segs)
    assert jk["jk_se"] is not None and jk["jk_se"] >= 0.0


def test_threshold_sweep_covers_the_rule_grid():
    s = skr.synthetic_skr()
    sw = skr.threshold_sweep(s["flux"], s["range_rs"], ks=(4.0, 6.0), baselines=(25.0,))
    assert [(r["k"], r["baseline_pct"]) for r in sw] == [(4.0, 25.0), (6.0, 25.0)]
    assert all(r["corrected_near_far"] is not None for r in sw)


def test_block_permutation_p_flags_a_real_rotation_and_passes_noise():
    s = skr.synthetic_skr(n_days=8.0)
    real = skr.block_permutation_p(s["jd"], s["flux"], n_perm=19, seed=1)
    assert real["p_perm"] <= 0.10  # injected dual period survives day-block permutation
    rng = np.random.default_rng(2)
    noise = 10.0 ** rng.normal(-14.0, 0.15, s["flux"].size)
    null = skr.block_permutation_p(s["jd"], noise, n_perm=19, seed=1)
    assert null["p_perm"] > 0.10  # pure noise does not
