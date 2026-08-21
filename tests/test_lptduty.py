"""Offline tests for the LPT duty-cycle constraint (plan 90)."""

from __future__ import annotations

import math

import pytest

from jansky_research import lptduty as ld

HEADER = "name,epoch,obs_id,band,epoch_mjd,duration_s,i_mjy,e_i,v_mjy,e_v,offset_arcsec,note\n"


def _row(name="SRC", obs="A1", mjd=59000.0, dur=730.0, v=0.1, ev=0.2, note=""):
    return f"{name},low,{obs},low,{mjd},{dur},0.5,0.2,{v},{ev},1.0,{note}\n"


def test_load_epochs_drops_rows_without_a_measurement(tmp_path):
    p = tmp_path / "epochs.csv"
    p.write_text(
        HEADER
        + _row(obs="ok")
        + _row(obs="unreleased", v="", ev="")  # archive never released it
        + _row(obs="zero_err", ev="0")  # a zero error is not a measurement
    )
    rows = ld.load_epochs(p)
    assert [r.obs_id for r in rows] == ["ok"]


def test_epoch_efficiency_is_monotonic_and_bounded():
    # A pulse far above the threshold is always seen; far below, never.
    assert ld.epoch_efficiency(0.2, 100.0) == pytest.approx(1.0)
    assert ld.epoch_efficiency(0.2, 0.0) == pytest.approx(0.0, abs=1e-6)
    # Exactly at threshold the measurement scatters either side: efficiency is 1/2.
    assert ld.epoch_efficiency(0.2, 5.0 * 0.2) == pytest.approx(0.5)
    # Deeper epochs (smaller sigma) are more efficient for a fixed pulse.
    assert ld.epoch_efficiency(0.1, 0.6) > ld.epoch_efficiency(0.5, 0.6)
    assert ld.epoch_efficiency(0.0, 1.0) == 0.0


def test_denominator_is_efficiency_weighted_not_epoch_count(tmp_path):
    """The frblens lesson: shallow epochs must not pad the exposure."""
    p = tmp_path / "e.csv"
    # one deep epoch (sigma 0.1) plus nine useless ones (sigma 50) -- no detections
    body = _row(obs="deep", v=0.0, ev=0.1) + "".join(
        _row(obs=f"shallow{i}", v=0.0, ev=50.0) for i in range(9)
    )
    p.write_text(HEADER + body)
    rows = ld.load_epochs(p)
    c = ld.duty_constraint(rows, pulse_mjy=1.0)
    assert c.n_epochs == 10
    # only the deep epoch could have seen a 1 mJy pulse
    assert c.effective_epochs == pytest.approx(1.0, abs=1e-3)
    assert c.n_detections == 0
    assert c.p_point is None
    # limit uses ~1 effective epoch, NOT 10 -- an epoch-count denominator would give 0.30
    assert c.p_upper_95 == pytest.approx(ld.POISSON_ZERO_95, rel=1e-2)


def test_no_sensitive_epoch_gives_an_infinite_limit_not_a_confident_one(tmp_path):
    p = tmp_path / "e.csv"
    p.write_text(HEADER + "".join(_row(obs=f"s{i}", v=0.0, ev=100.0) for i in range(20)))
    rows = ld.load_epochs(p)
    c = ld.duty_constraint(rows, pulse_mjy=0.5)
    assert c.effective_epochs == pytest.approx(0.0, abs=1e-9)
    assert math.isinf(c.p_upper_95)
    assert c.p_point is None


def test_detections_counted_by_threshold_and_point_estimate(tmp_path):
    p = tmp_path / "e.csv"
    body = (
        _row(obs="det", v=2.0, ev=0.2)  # 10 sigma
        + _row(obs="neg_det", v=-1.5, ev=0.2)  # -7.5 sigma: |V| counts, V can be either sign
        + "".join(_row(obs=f"n{i}", v=0.05, ev=0.2) for i in range(8))
    )
    p.write_text(HEADER + body)
    rows = ld.load_epochs(p)
    c = ld.duty_constraint(rows, pulse_mjy=2.0)
    assert c.n_detections == 2
    assert c.effective_epochs == pytest.approx(10.0, rel=1e-3)
    assert c.p_point == pytest.approx(0.2, rel=1e-3)
    assert c.p_upper_95 > c.p_point


def test_limit_weakens_as_the_assumed_pulse_gets_fainter(tmp_path):
    """Every limit is a function of the assumed pulse flux; state it, don't hide it."""
    p = tmp_path / "e.csv"
    p.write_text(HEADER + "".join(_row(obs=f"e{i}", v=0.0, ev=0.3) for i in range(30)))
    rows = ld.load_epochs(p)
    bright = ld.duty_constraint(rows, pulse_mjy=5.0)
    faint = ld.duty_constraint(rows, pulse_mjy=1.4)
    assert bright.p_upper_95 < faint.p_upper_95
    assert bright.effective_epochs > faint.effective_epochs


def test_rows_must_belong_to_one_source(tmp_path):
    p = tmp_path / "e.csv"
    p.write_text(HEADER + _row(name="A") + _row(name="B"))
    rows = ld.load_epochs(p)
    with pytest.raises(ValueError, match="more than one source"):
        ld.duty_constraint(rows, pulse_mjy=1.0)
    with pytest.raises(ValueError, match="no epochs"):
        ld.duty_constraint([], pulse_mjy=1.0)


def test_read_periods_from_the_vendored_catalogue():
    periods = ld.read_periods("data/lpt_sample.csv")
    assert periods["ASKAP J183950.5-075635"] == pytest.approx(23221.74)
    assert all(v > 0 for v in periods.values())


def _rows_at(mjds, name="SRC", ev=0.2):
    return [
        ld.EpochRow(name=name, obs_id=f"o{i}", epoch_mjd=m, duration_s=730.0, v_mjy=0.0, e_v=ev)
        for i, m in enumerate(mjds)
    ]


def test_phase_sampling_flags_a_cadence_locked_to_the_period():
    """The failure mode this exists to catch: cadence an exact multiple of the period."""
    period = 3600.0  # 1 h
    step_days = 10 * period / 86400.0  # every 10 whole periods -> identical phase every time
    rows = _rows_at([59000.0 + i * step_days for i in range(60)], ev=0.2)
    ps = ld.phase_sampling(rows, period)
    assert ps.rayleigh_z > 50  # massively concentrated
    assert ps.rayleigh_p < 1e-6
    assert ps.kuiper_v > 0.5


def test_phase_sampling_passes_for_incommensurate_sampling():
    period = 2656.2554  # ASKAP J1832-0911
    # irrational-ish spacing in units of the period -> phases spread
    rows = _rows_at([59000.0 + i * 13.719 for i in range(120)])
    ps = ld.phase_sampling(rows, period)
    assert ps.rayleigh_p > 0.01
    assert ps.kuiper_v < 0.4


def test_phase_sampling_is_invariant_to_the_arbitrary_zero_point():
    """No ephemeris is needed to test uniformity: a phase shift cannot change clustering."""
    period = 4186.3285
    base = [59000.0 + i * 7.3 for i in range(80)]
    a = ld.phase_sampling(_rows_at(base), period)
    shifted = [m + 0.25 * period / 86400.0 for m in base]
    b = ld.phase_sampling(_rows_at(shifted), period)
    assert a.rayleigh_z == pytest.approx(b.rayleigh_z, rel=1e-9)
    assert a.kuiper_v == pytest.approx(b.kuiper_v, rel=1e-9)


def test_required_period_precision_scales_as_p_squared_over_baseline():
    rows = _rows_at([59000.0, 60460.0])  # 1460 d baseline
    ps = ld.phase_sampling(rows, 3600.0)
    assert ps.baseline_days == pytest.approx(1460.0)
    assert ps.required_period_precision_s == pytest.approx(0.1 * 3600.0**2 / (1460.0 * 86400.0))
    # a longer period needs a *less* stringent absolute precision for the same baseline
    ps2 = ld.phase_sampling(rows, 7200.0)
    assert ps2.required_period_precision_s > ps.required_period_precision_s


def test_phase_sampling_rejects_bad_input():
    rows = _rows_at([59000.0, 59001.0])
    with pytest.raises(ValueError, match="period must be positive"):
        ld.phase_sampling(rows, 0.0)
    with pytest.raises(ValueError, match="at least two epochs"):
        ld.phase_sampling(rows[:1], 3600.0)


def _eph(P=1000.0, sigma=1e-6, pepoch=59000.0, width=100.0):
    return ld.Ephemeris(
        name="SRC",
        period_s=P,
        sigma_period_s=sigma,
        pepoch_mjd=pepoch,
        pulse_width_s=width,
        reference="test",
    )


def test_phase_uncertainty_grows_with_time_and_is_inf_when_unreported():
    e = _eph(P=1000.0, sigma=1e-5, pepoch=59000.0)
    near = e.phase_uncertainty_at(59001.0)
    far = e.phase_uncertainty_at(59100.0)
    assert far == pytest.approx(100 * near)
    assert far == pytest.approx(99.0 * 86400.0 * 1e-5 / 1000.0**2 + near, rel=1e-6)
    unknown = ld.Ephemeris("SRC", 1000.0, None, 59000.0, 100.0, "no sigma published")
    assert math.isinf(unknown.phase_uncertainty_at(59001.0))


def test_phase_resolved_recovers_f_active_when_the_source_is_always_on():
    """Snapshots placed on the pulse, all detected -> f_active ~ 1."""
    e = _eph(P=1000.0, sigma=1e-9, width=200.0)
    # place each snapshot at whole periods after PEPOCH -> phase 0 == the pulse
    rows = [
        ld.EpochRow("SRC", f"o{i}", 59000.0 + i * 1000.0 / 86400.0, 100.0, 5.0, 0.2)
        for i in range(20)
    ]
    pr = ld.phase_resolved_activity(rows, e, pulse_mjy=5.0)
    assert pr.usable
    assert pr.n_on_window == 20
    assert pr.n_detections_in_window == 20
    assert pr.f_active_point == pytest.approx(1.0, rel=0.05)


def test_phase_resolved_gives_an_upper_limit_when_on_window_snapshots_see_nothing():
    e = _eph(P=1000.0, sigma=1e-9, width=200.0)
    rows = [
        ld.EpochRow("SRC", f"o{i}", 59000.0 + i * 1000.0 / 86400.0, 100.0, 0.0, 0.2)
        for i in range(30)
    ]
    pr = ld.phase_resolved_activity(rows, e, pulse_mjy=5.0)
    assert pr.n_on_window == 30
    assert pr.n_detections_in_window == 0
    assert pr.f_active_point is None
    assert pr.f_active_upper_95 == pytest.approx(ld.POISSON_ZERO_95 / 30, rel=1e-3)


def test_detections_outside_the_window_are_reported_not_swallowed():
    """A pulse at the wrong phase means the ephemeris is wrong; that must be visible."""
    e = _eph(P=1000.0, sigma=1e-9, width=50.0)
    # snapshots half a period out of phase, yet bright
    rows = [
        ld.EpochRow("SRC", f"o{i}", 59000.0 + (i + 0.5) * 1000.0 / 86400.0, 60.0, 9.0, 0.2)
        for i in range(10)
    ]
    pr = ld.phase_resolved_activity(rows, e, pulse_mjy=9.0)
    assert pr.n_on_window == 0
    assert pr.n_detections_outside == 10


def test_unusable_when_phase_has_drifted_beyond_tolerance():
    # 1e-3 s period error accumulates fast on a 1000 s period over ~5 years
    e = _eph(P=1000.0, sigma=1e-3, pepoch=59000.0)
    rows = [ld.EpochRow("SRC", "o", 61000.0, 730.0, 0.0, 0.2)] * 2
    pr = ld.phase_resolved_activity(list(rows), e, pulse_mjy=5.0)
    assert not pr.usable
    assert pr.max_phase_uncertainty > 0.1


def test_window_fraction_includes_the_snapshot_length():
    """The window is (w + T)/P, not w/P: a snapshot covers phase while it integrates."""
    e = _eph(P=1000.0, sigma=1e-9, width=100.0)
    rows = [
        ld.EpochRow("SRC", "o", 59000.0, 400.0, 0.0, 0.2),
        ld.EpochRow("SRC", "p", 59001.0, 400.0, 0.0, 0.2),
    ]
    pr = ld.phase_resolved_activity(rows, e, pulse_mjy=5.0)
    assert pr.window_fraction == pytest.approx((100.0 + 400.0) / 1000.0)


def test_poisson_upper_limits_match_the_exact_cdf_solution():
    """2.996 + k is wrong for k>0 and this pins the right values (see poisson_upper_95)."""
    expected = {0: 2.9957, 1: 4.7439, 2: 6.2958, 3: 7.7537}
    for k, want in expected.items():
        assert ld.poisson_upper_95(k) == pytest.approx(want, abs=1e-3)
    # the naive formula is low by ~19% at k=1 -- guard against a regression to it
    assert ld.poisson_upper_95(1) > ld.POISSON_ZERO_95 + 1
    with pytest.raises(ValueError, match="non-negative"):
        ld.poisson_upper_95(-1)


def test_leakage_veto_rejects_a_significant_but_leakage_scale_v(tmp_path):
    """lptv's criterion is significance AND |V| > 0.006|I|; this module claimed both."""
    p = tmp_path / "e.csv"
    header = "name,epoch,obs_id,band,epoch_mjd,duration_s,i_mjy,e_i,v_mjy,e_v,offset_arcsec,note\n"
    p.write_text(
        header
        + "S,low,leak,low,59000.0,730.0,5000.0,1.0,20.0,0.2,1.0,\n"  # 100 sigma but |V| < 0.6% of I
        + "S,low,real,low,59001.0,730.0,50.0,1.0,20.0,0.2,1.0,\n"  # same V, modest I -> genuine
    )
    rows = ld.load_epochs(p)
    assert [r.i_mjy for r in rows] == [5000.0, 50.0]
    c = ld.duty_constraint(rows, pulse_mjy=20.0)
    assert c.n_detections == 1  # the leakage-scale one is vetoed
