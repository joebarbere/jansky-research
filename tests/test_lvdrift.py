"""Tests for jansky_research.lvdrift -- the frozen within-burst drift estimator. Offline."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import lvdrift as lv

F = lv.lofar_freqs()


def test_lofar_freqs_cover_the_lv614_band():
    assert F.size == 162
    assert F[0] == pytest.approx(29.297, abs=1e-3) and F[-1] == pytest.approx(60.742, abs=1e-3)


@pytest.mark.parametrize("a", [0.0, 0.43, 1.0, 1.23, 1.84])
def test_planted_index_is_recovered(a):
    A = 4.0 / 45**a  # ~4 MHz/s near 45 MHz, the paper's type III mean scale
    est = [lv.burst_index(lv.synthetic_burst(F, a=a, A=A, peak=20, seed=s), F) for s in range(6)]
    ok = [r["a"] for r in est if r["ok"]]
    assert len(ok) >= 5
    assert np.median(ok) == pytest.approx(a, abs=0.15)


def test_fit_is_smooth_through_a_equals_one():
    t = 10.0 - lv._g(F, 1.0, float(np.exp(np.mean(np.log(F))))) / 0.1
    r = lv.fit_index(F, t)
    assert r["a"] == pytest.approx(1.0, abs=1e-3) and r["A"] == pytest.approx(0.1, rel=1e-3)


def test_pure_noise_and_shuffled_bursts_fail_the_quality_cut():
    noise = np.random.default_rng(3).normal(size=(120, F.size))
    assert not lv.burst_index(noise, F)["ok"]
    burst = lv.synthetic_burst(F, a=1.23, A=4.0 / 45**1.23, peak=20, seed=4)
    assert lv.burst_index(burst, F)["ok"]
    assert not lv.burst_index(lv.shuffle_channels(burst, seed=1), F)["ok"]


def test_quality_reasons_are_reported():
    q = lv.ridge_quality(np.array([50.0, 40.0]), np.array([0.0, 5.0]))
    assert not q["ok"] and "channels" in q["reason"]
    q = lv.ridge_quality(F, np.linspace(0, 10, F.size))  # time rising WITH frequency: wrong way
    assert not q["ok"] and "non-monotonic" in q["reason"]
    assert lv.ridge_quality(np.array([]), np.array([]))["reason"] == "no ridge"


def test_peak_drift_matches_the_local_slope():
    a, A = 1.0, 4.0 / 45.0
    burst = lv.synthetic_burst(F, a=a, A=A, peak=30, width_s=2.0, seed=7)
    d = lv.peak_drift(burst, F)
    # the burst's peak flux is near the top of the band at t_hi; the local drift there is A f^a
    assert 2.0 < d < 8.0
    assert np.isnan(lv.peak_drift(np.zeros((50, F.size)), F))


def test_fit_failure_is_flagged_not_reported():
    r = lv.fit_index(np.array([50.0, 50.0, 50.0]), np.array([1.0, 2.0, 3.0]))
    assert r["a_at_bound"] or r["unconstrained"]
