"""Offline tests for the Taylor-term in-band spectral index (plan 91)."""

from __future__ import annotations

import math

import pytest

from jansky_research import lptspec as ls


def test_taylor1_is_about_eleven_times_noisier_for_racs():
    """The short lever arm is the whole reason a 5-sigma detection is not enough."""
    ratio = ls.taylor1_noise_ratio()
    assert ratio == pytest.approx(math.sqrt(12.0) * ls.RACS_NU0_HZ / ls.RACS_BANDWIDTH_HZ)
    assert 10.0 < ratio < 11.5
    # a wider fractional band would help: doubling the bandwidth halves the penalty
    assert ls.taylor1_noise_ratio(bandwidth_hz=2 * ls.RACS_BANDWIDTH_HZ) == pytest.approx(ratio / 2)


def test_usable_threshold_matches_the_uncertainty_formula():
    for target in (0.1, 0.3, 0.5):
        # the idealised threshold is where the idealised uncertainty hits the target...
        snr_ideal = ls.usable_snr_threshold(target, realistic=False)
        assert ls.alpha_uncertainty(snr_ideal, 1.0, 0.0) == pytest.approx(target, rel=1e-9)
        # ...and the default (realistic) one is where the penalised uncertainty does
        snr_real = ls.usable_snr_threshold(target)
        assert ls.realistic_sigma_alpha(snr_real, 1.0, 0.0) == pytest.approx(target, rel=1e-9)
    with pytest.raises(ValueError, match="must be positive"):
        ls.usable_snr_threshold(0.0)


def test_alpha_uncertainty_scales_inversely_with_snr():
    a = ls.alpha_uncertainty(100.0, 1.0, -70.0)
    b = ls.alpha_uncertainty(200.0, 1.0, -140.0)  # same alpha, twice the S/N
    assert b == pytest.approx(a / 2, rel=1e-9)
    with pytest.raises(ValueError, match="positive"):
        ls.alpha_uncertainty(0.0, 1.0)


def test_taylor_alpha_and_zero_denominator():
    assert ls.taylor_alpha(10.0, -7.0) == pytest.approx(-0.7)
    with pytest.raises(ValueError, match="undefined"):
        ls.taylor_alpha(0.0, 1.0)


def test_injection_recovery_scatter_matches_the_analytic_prediction():
    """At high S/N the linearised uncertainty should match the simulated scatter."""
    snr, alpha = 150.0, -0.7
    r = ls.injection_recovery(snr, alpha, n_trials=40000, seed=3)
    predicted = ls.alpha_uncertainty(snr, 1.0, alpha * snr)
    assert r.alpha_std == pytest.approx(predicted, rel=0.05)
    assert abs(r.bias) < 0.01
    assert r.frac_within_0p3 > 0.99


def test_injection_recovery_detects_bias_when_the_denominator_gets_noisy():
    """A bias estimator that never finds bias proves nothing; this shows it can."""
    strong = ls.injection_recovery(300.0, -0.7, n_trials=60000, seed=4)
    weak = ls.injection_recovery(8.0, -0.7, n_trials=60000, seed=4)
    assert abs(strong.bias) < 0.002
    assert abs(weak.bias) > 5 * abs(strong.bias)
    assert weak.alpha_std > 10 * strong.alpha_std


def test_recovery_degrades_monotonically_as_snr_falls():
    prev = None
    for snr in (300.0, 100.0, 50.0, 25.0):
        r = ls.injection_recovery(snr, -0.7, n_trials=20000, seed=5)
        if prev is not None:
            assert r.alpha_std > prev
        prev = r.alpha_std


def test_measure_flags_usability_at_the_real_pulse_snrs():
    """The four measured VAST pulses: only the brightest two clear a 0.3 bar."""
    for snr, expected in ((262.2, True), (57.3, True), (36.1, True), (21.0, False)):
        m = ls.measure(t0=snr, sigma0=1.0, t1=-0.7 * snr)
        assert m.snr_taylor0 == pytest.approx(snr)
        assert m.usable is expected, (snr, m.sigma_alpha)


def test_injection_recovery_rejects_bad_input():
    with pytest.raises(ValueError, match="snr must be positive"):
        ls.injection_recovery(0.0, -0.7)


def test_realistic_sigma_alpha_reproduces_the_published_mt_mfs_behaviour():
    """Rashid+2024: MT-MFS in-band alpha errors are >~0.2 for SNR <~100."""
    assert ls.realistic_sigma_alpha(100.0, 1.0, 0.0) == pytest.approx(0.213, abs=0.01)
    # the idealised floor is optimistic by exactly the penalty
    assert ls.realistic_sigma_alpha(100.0, 1.0, 0.0) == pytest.approx(
        ls.MT_MFS_PENALTY * ls.alpha_uncertainty(100.0, 1.0, 0.0)
    )


def test_realistic_threshold_is_stricter_than_the_idealised_one():
    ideal = ls.usable_snr_threshold(0.3, realistic=False)
    real = ls.usable_snr_threshold(0.3)
    assert real == pytest.approx(ls.MT_MFS_PENALTY * ideal)
    assert 34 < ideal < 38 and 69 < real < 73


def test_the_penalty_changes_which_pulses_qualify():
    """It is the difference between five usable pulses and three -- pin it."""
    snrs = [262.2, 213.9, 108.8, 57.3, 36.1, 21.0, 15.0]
    ideal_ok = sum(ls.alpha_uncertainty(s, 1.0, 0.0) <= 0.3 for s in snrs)
    real_ok = sum(ls.realistic_sigma_alpha(s, 1.0, 0.0) <= 0.3 for s in snrs)
    assert ideal_ok == 5
    assert real_ok == 3


class _FakeCol(list):
    pass


class _FakeTable:
    """Minimal stand-in for a CASDA ObsCore table (only ``filename`` is read)."""

    def __init__(self, names):
        self._d = {"filename": _FakeCol(names)}

    def __getitem__(self, k):
        return self._d[k]


def test_taylor_mask_selects_the_right_sbid_and_term():
    """The SBID filter is the whole point: the same position has many epochs."""
    t = _FakeTable(
        [
            "image.i.VAST_1824-06.SB60804.cont.taylor.0.restored.conv.fits",  # want
            "image.i.VAST_1824-06.SB60804.cont.taylor.1.restored.conv.fits",  # other term
            "image.i.VAST_1824-06.SB99999.cont.taylor.0.restored.conv.fits",  # other epoch
            "image.v.VAST_1824-06.SB60804.cont.taylor.0.restored.conv.fits",  # other stokes
            "noiseMap.image.i.VAST_1824-06.SB60804.cont.taylor.0.restored.conv.fits",  # noise
            "image.i.VAST_1824-06.SB60804.cont.taylor.0.fits",  # not restored/conv
        ]
    )
    m = ls.taylor_science_mask(t, term=0, sbid="60804")
    assert m.tolist() == [True, False, False, False, False, False]
    m1 = ls.taylor_science_mask(t, term=1, sbid="60804")
    assert m1.tolist() == [False, True, False, False, False, False]
    assert not ls.taylor_science_mask(t, term=0, sbid="12345").any()
