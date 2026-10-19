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


def test_ring_center_finds_an_offset_crescent_not_its_centroid():
    from scipy.ndimage import shift

    img = shift(_crescent(asym=0.8, phi0=0.0), (6.0, -9.0), order=1)
    cy, cx = om.ring_center(img, search_px=14, step=0.5)
    c = (om.N_GRID - 1) / 2.0
    assert abs(cy - (c + 6.0)) <= 1.0 and abs(cx - (c - 9.0)) <= 1.0
    # the centroid is pulled toward the bright side (phi0 = 0 -> +x), the ring centre is not
    n = om.N_GRID
    y, x = np.mgrid[0:n, 0:n]
    assert (img * x).sum() / img.sum() > (c - 9.0) + 2.0


def test_center_on_ring_uses_one_common_shift():
    from scipy.ndimage import shift

    a = shift(_crescent(phi0=1.0), (5.0, 4.0), order=1)
    b = shift(_crescent(phi0=1.1), (5.0, 4.0), order=1)
    ca, cb = om.center_on_ring(a, b)
    assert om.ring_radius(0.5 * (ca + cb), PIX) == pytest.approx(21.0, abs=1.0)
    # identical inputs stay identical after centring -> q is still exactly 1
    assert om.q_estimator(a, a, PIX, 1.0)["q"] == pytest.approx(1.0, abs=1e-12)


def test_patch_pipeline_changes_exactly_one_line():
    src = "a = 1\nttype     = 'nfft'   # comment\nb = 2\n"
    out = om.patch_pipeline(src)
    assert out.count("\n") == src.count("\n")
    assert "ttype = 'fast'" in out and "nfft'" not in out.split("\n")[1]
    with pytest.raises(ValueError):
        om.patch_pipeline("a = 1\n")


def test_units_and_dz_are_physical():
    # 1 uas at 16.8 Mpc ~ 2.5e9 km; the Earth's line-of-sight motion in a day is << that
    assert om.uas_to_km(1.0) == pytest.approx(2.513e9, rel=1e-3)
    dz = om.earth_dz_one_day(57848.0)  # 2017-04-05
    assert 1e4 < dz < 2.6e6  # at most the Earth's orbital speed x 1 day (~2.57e6 km)


def test_crescent_params_recovers_the_model():
    img = _crescent(asym=0.5, phi0=0.7)
    p = om.crescent_params(img, PIX)
    assert p["radius"] == pytest.approx(21.0, abs=PIX)
    assert p["asym"] == pytest.approx(0.5, abs=0.05)
    assert p["phi0"] == pytest.approx(0.7, abs=0.05)


def test_q_sweep_keys_and_identity():
    img = _crescent()
    s = om.q_sweep(img, img, PIX, 1e5)
    assert set(s) == {f"q_x{f:g}" for f in om.DZ_SWEEP}
    assert all(v == pytest.approx(1.0, abs=1e-12) for v in s.values())


def test_tie_is_blind_to_a_phase_vortex_on_a_symmetric_ring():
    # a twisted wavefront on an azimuthally symmetric ring changes no intensity: the TIE sees nothing
    r = om.vortex_injection(om.N_GRID, 0.512, 19.97, 15.87, 0.0, 0.0, 0.3)
    assert r["q_true"] > 1.5
    assert r["q_recovered"] == pytest.approx(1.0, abs=1e-9) and r["phase_rms"] == 0.0


def test_tie_misrecovers_a_vortex_on_a_crescent():
    r = om.vortex_injection(om.N_GRID, 0.512, 19.97, 15.87, 0.535, -2.15, 1.0)
    assert r["q_true"] > 50.0  # an l = 1 vortex: power almost all in m = +1
    assert r["q_recovered"] < 1.2  # recovered field barely asymmetric


def test_paper_metrics_from_committed_results(tmp_path):
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "results"
    m = om.paper_metrics(str(root))
    assert m["syn_n1_n"] == 800 and m["syn_pzero_sym_rec"] == 1.0
    assert m["syn_n1_lo_pone"] <= 1.393 <= m["syn_n1_hi_pone"]
    out = tmp_path / "macros.tex"
    om.write_paper_macros(str(out), str(root))
    tex = out.read_text()
    assert r"\oamspinSynNoneN}{800}" in tex and "--" not in tex.split("Source")[1]
    assert om.paper_metrics(str(tmp_path))["syn_pzero_sym_rec"] == 1.0  # P0 needs no files


def test_paper_figure_renders(tmp_path):
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "results"
    out = tmp_path / "f.pdf"
    om.paper_figure(str(root), str(out))
    assert out.stat().st_size > 1000
