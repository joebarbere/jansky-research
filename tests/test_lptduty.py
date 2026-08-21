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
