"""Tests for jansky_research.oamspin -- the frozen twisted-light estimator. Offline."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import oamspin as om

PIX = 0.48  # uas per pixel: 120 uas field on 250 pixels


def _crescent(asym=0.6, phi0=1.0):
    return om.crescent_model(om.N_GRID, PIX, 21.0, 8.0, asym, phi0)


def test_identical_pair_gives_zero_phase_and_q_one():
    img = _crescent()
    assert not np.any(om.tie_phase(img, img, 1.0, PIX))
    assert om.q_estimator(img, img, PIX, 1.0)["q"] == pytest.approx(1.0, abs=1e-12)


def test_real_brightness_alone_has_symmetric_spectrum():
    img = _crescent(asym=0.9)  # strongly lopsided, but real: |c_+1| = |c_-1|
    spec = om.spiral_spectrum(np.sqrt(img).astype(complex), PIX, om.ring_radius(img, PIX))
    assert om.q_from_spectrum(spec) == pytest.approx(1.0, abs=1e-9)
    assert spec[1] > 0.01 * spec[0]  # the asymmetry does put power at |m| = 1


def test_planted_helical_modes_are_recovered_with_their_sign():
    for ell, expect_big in ((+1, True), (-1, False)):
        f = om.laguerre_gauss_field(om.N_GRID, ell)
        rp = om.ring_radius(np.abs(f) ** 2, 1.0)
        q = om.q_from_spectrum(om.spiral_spectrum(f, 1.0, rp))
        assert (q > 100.0) if expect_big else (q < 0.01)


def _forward_didz(i0, p, h):
    """Spectral forward model of SM Eq. 11a: dI/dz = -(1/2pi) div(I grad P)."""
    n = i0.shape[0]
    k = 2 * np.pi * np.fft.fftfreq(n, d=h)
    kx, ky = np.meshgrid(k, k)
    pk = np.fft.fft2(p)
    gx = np.real(np.fft.ifft2(1j * kx * pk))
    gy = np.real(np.fft.ifft2(1j * ky * pk))
    div = np.real(np.fft.ifft2(1j * kx * np.fft.fft2(i0 * gx) + 1j * ky * np.fft.fft2(i0 * gy)))
    return -div / (2 * np.pi)


def _p_true(n):
    y, x = np.mgrid[0:n, 0:n] / n
    return 0.3 * np.sin(2 * np.pi * x) + 0.2 * np.cos(4 * np.pi * y)


def test_tie_recovers_a_known_phase_exactly_for_uniform_intensity():
    n, dz = 128, 1e-3
    i0 = np.ones((n, n))
    pt = _p_true(n)
    d = _forward_didz(i0, pt, 1.0 / n)
    p = om.tie_phase(i0 - 0.5 * dz * d, i0 + 0.5 * dz * d, dz, 1.0 / n)
    assert np.max(np.abs((p - p.mean()) - (pt - pt.mean()))) < 1e-6


def test_tie_is_approximate_for_nonuniform_intensity():
    # Teague's FFT solution drops the curl part of I grad P; with I varying by a factor 3 the
    # recovered phase is close but not exact -- the frozen method's stated approximation.
    n, dz = 128, 1e-3
    y, x = np.mgrid[0:n, 0:n] / n
    i0 = 1.0 + 0.5 * np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y)
    pt = _p_true(n)
    d = _forward_didz(i0, pt, 1.0 / n)
    p = om.tie_phase(i0 - 0.5 * dz * d, i0 + 0.5 * dz * d, dz, 1.0 / n)
    err = np.std((p - p.mean()) - (pt - pt.mean())) / np.std(pt)
    assert 1e-4 < err < 0.15


def test_q_estimator_scales_with_dz_for_a_fixed_difference():
    # the estimator's phase scales as 1/dz: the same image pair at a smaller dz gives a larger phase
    a = _crescent(phi0=1.0)
    b = _crescent(phi0=1.05)  # a slightly rotated crescent: a nonzero day-to-day difference
    big = om.q_estimator(a, b, PIX, 1.0)
    small = om.q_estimator(a, b, PIX, 10.0)
    assert big["phase_rms"] == pytest.approx(10.0 * small["phase_rms"], rel=1e-6)


def test_ring_radius_and_crescent():
    img = _crescent()
    assert om.ring_radius(img, PIX) == pytest.approx(21.0, abs=PIX)
    assert np.all(img >= 0)
