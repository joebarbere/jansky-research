"""Offline tests for LPT linear polarization (plan 92)."""

from __future__ import annotations

import math

import pytest

from jansky_research import lptlin as ll


def test_debias_removes_the_ricean_floor():
    """L = sqrt(Q^2+U^2) is positive even for an unpolarized source; that must not survive."""
    # a pure-noise position: L ~ sigma, so the debiased value must be zero, not "small"
    assert ll.debias_linear(1.0, 1.0) == 0.0
    assert ll.debias_linear(0.8, 1.0) == 0.0
    # a real signal keeps most of its amplitude
    assert ll.debias_linear(10.0, 1.0) == pytest.approx(math.sqrt(99.0))
    with pytest.raises(ValueError, match="sigma_qu must be positive"):
        ll.debias_linear(1.0, 0.0)


def test_unpolarized_noise_does_not_manufacture_polarization():
    """The failure this exists to prevent: L/I positive for every faint source."""
    p = ll.linear_fraction(q_mjy=0.2, u_mjy=-0.15, i_mjy=100.0, sigma_qu_mjy=0.3)
    assert p.l_raw_mjy > 0  # the raw quantity is positive, as it always is
    assert p.l_debiased_mjy == 0.0  # and the debiased one is not
    assert p.linear_fraction == 0.0
    assert not p.detected


def test_leakage_veto_rejects_a_significant_but_instrumental_linear_signal():
    """Same two-part criterion lptv applies to V: significance AND above the leakage floor."""
    # L is 20 sigma but only 0.4% of a very bright I -> below the 0.6% floor
    p = ll.linear_fraction(q_mjy=6.0, u_mjy=0.0, i_mjy=1500.0, sigma_qu_mjy=0.3)
    assert p.significance > 15
    assert not p.above_leakage
    assert not p.detected
    # same L against a fainter I is a genuine detection
    q = ll.linear_fraction(q_mjy=6.0, u_mjy=0.0, i_mjy=100.0, sigma_qu_mjy=0.3)
    assert q.above_leakage and q.detected


def test_fraction_survives_pulse_dilution_while_the_pulse_is_above_the_noise():
    """The property that makes this work where plan 91 failed -- and its limit.

    A pulse on for a fraction f of the synthesis is diluted to f*I, f*Q, f*U alike, so the
    raw ratio is preserved exactly. Debiasing breaks the exact cancellation, because the
    noise does not dilute with the signal: the recovered fraction is biased slightly DOWN
    as the pulse approaches the noise floor. Down is the safe direction -- it cannot invent
    polarization.
    """
    full = ll.linear_fraction(q_mjy=30.0, u_mjy=40.0, i_mjy=200.0, sigma_qu_mjy=0.5)
    ratios = {}
    for f in (0.5, 0.2, 0.05):
        d = ll.linear_fraction(q_mjy=30.0 * f, u_mjy=40.0 * f, i_mjy=200.0 * f, sigma_qu_mjy=0.5)
        ratios[f] = d.linear_fraction / full.linear_fraction
        # the angle is completely dilution-invariant: no debiasing enters it
        assert d.angle_deg == pytest.approx(full.angle_deg, abs=1e-9)
    assert ratios[0.5] == pytest.approx(1.0, abs=0.001)
    assert ratios[0.2] == pytest.approx(1.0, abs=0.005)
    assert 0.95 < ratios[0.05] < 1.0  # measurably low, and low rather than high
    # the deviation grows monotonically as the pulse is diluted further
    assert ratios[0.5] > ratios[0.2] > ratios[0.05]


def test_polarization_angle_wraps_to_the_conventional_range():
    assert ll.polarization_angle(1.0, 0.0) == pytest.approx(0.0)
    assert ll.polarization_angle(0.0, 1.0) == pytest.approx(45.0)
    assert ll.polarization_angle(-1.0, 0.0) == pytest.approx(90.0)
    # the EVPA is a 180-degree quantity: Q,U and -Q,-U are the same angle
    assert ll.polarization_angle(3.0, 4.0) == pytest.approx(
        ll.polarization_angle(-3.0, -4.0) - 90.0
    ) or ll.polarization_angle(3.0, 4.0) == pytest.approx(
        (ll.polarization_angle(-3.0, -4.0) + 90.0) % 180.0
    )
    assert 0.0 <= ll.polarization_angle(-2.0, -5.0) < 180.0


def test_rejects_undefined_inputs():
    with pytest.raises(ValueError, match="Stokes I is zero"):
        ll.linear_fraction(1.0, 1.0, 0.0, 0.3)
    with pytest.raises(ValueError, match="sigma_qu must be positive"):
        ll.linear_fraction(1.0, 1.0, 10.0, 0.0)


def test_total_polarization_is_a_physical_bound_check():
    """Q^2+U^2+V^2 <= I^2: a total fraction above 1 means the measurement is wrong."""
    # the real ASKAP J1832-0911 numbers: 10.9% linear, 5.6% circular
    frac = ll.total_polarization(l_mjy=27.212, v_mjy=14.04, i_mjy=249.968)
    assert frac == pytest.approx(0.1224, abs=0.001)
    assert frac < 1.0
    # an unphysical combination is detectable by the same call
    assert ll.total_polarization(l_mjy=90.0, v_mjy=60.0, i_mjy=100.0) > 1.0
    with pytest.raises(ValueError, match="Stokes I is zero"):
        ll.total_polarization(1.0, 1.0, 0.0)
