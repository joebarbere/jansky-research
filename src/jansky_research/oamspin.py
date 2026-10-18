"""Can the "twisted light" estimator see black-hole spin in EHT 2017 data? (plan 100, ``oamspin``).

Tamburini, Thide & Della Valle (2020, MNRAS 492, L22; arXiv:1904.07923) recover a phase P from
two intensity images one day apart with the Transport-of-Intensity Equation (their SM Eq. 11a),
build the field sqrt(I) e^{iP}, Fourier transform it in azimuth on circles of pixels ("spiral
spectrum") and read spin from q = |c_{+1}| / |c_{-1}|. This module implements that estimator once,
as the plan freezes it, so that it can be run on non-spinning nulls, simultaneous-band pairs and
rescaled Delta-z -- a test of the estimator, not a reproduction of the authors' code (which is not
public).

Two properties the tests pin down: a real brightness alone gives |c_{+1}| = |c_{-1}| (so q = 1 when
P = 0), and identical input images give dI/dz = 0, hence P = 0 and q = 1 exactly.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "N_GRID",
    "TIKHONOV_FRAC",
    "ANNULUS_HALF_UAS",
    "N_AZ",
    "tie_phase",
    "ring_radius",
    "spiral_spectrum",
    "q_from_spectrum",
    "q_estimator",
    "laguerre_gauss_field",
    "crescent_model",
]

# ---- frozen estimator constants (plans/100, 2026-10-17) -----------------------------------------
N_GRID = 250  # the paper's 250 x 250 intensity matrices (SM Sect. 4.2)
TIKHONOV_FRAC = 1e-3  # floor on I in the TIE division, as a fraction of max(I)
ANNULUS_HALF_UAS = 5.0  # circles within r_pk +/- this are summed
N_AZ = 360  # azimuthal samples per circle


def tie_phase(i1: np.ndarray, i2: np.ndarray, dz: float, pixel: float) -> np.ndarray:
    r"""Phase P from the Transport-of-Intensity Equation, :math:`\partial_z I = -\frac{1}{2\pi}
    \nabla\cdot(I\nabla P)` (Tamburini et al. 2020, SM Eq. 11a), by the standard FFT method.

    With :math:`I\nabla P = \nabla\psi` (Teague's auxiliary potential), :math:`\nabla^2\psi =
    -2\pi\,\partial_z I`, then :math:`\nabla P = \nabla\psi / I` and :math:`\nabla^2 P =
    \nabla\cdot(\nabla\psi/I)`. Both Poisson problems are solved spectrally (periodic boundary,
    zero-mean). ``I`` is the mean of the two images with a Tikhonov floor of ``TIKHONOV_FRAC``
    x max(I). ``dz`` and ``pixel`` share a length unit. Identical images give P = 0 exactly.

    The FFT solution assumes :math:`I\nabla P` is curl-free (exact for uniform I); for a
    lopsided ring the curl part is dropped, the method's standard approximation (tests bound it).
    """
    i1 = np.asarray(i1, float)
    i2 = np.asarray(i2, float)
    didz = (i2 - i1) / dz
    if not np.any(didz):
        return np.zeros_like(i1)
    ii = 0.5 * (i1 + i2)
    ii = np.maximum(ii, TIKHONOV_FRAC * float(np.max(ii)))
    ny, nx = ii.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=pixel)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=pixel)
    kxx, kyy = np.meshgrid(kx, ky)
    k2 = kxx**2 + kyy**2
    k2[0, 0] = 1.0

    def inv_lap(f: np.ndarray) -> np.ndarray:
        fk = np.fft.fft2(f) / (-k2)
        fk[0, 0] = 0.0
        return np.real(np.fft.ifft2(fk))

    def grad(f: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        fk = np.fft.fft2(f)
        return np.real(np.fft.ifft2(1j * kxx * fk)), np.real(np.fft.ifft2(1j * kyy * fk))

    psi = inv_lap(-2.0 * np.pi * didz)
    gx, gy = grad(psi)
    gx, gy = gx / ii, gy / ii
    div = np.real(np.fft.ifft2(1j * kxx * np.fft.fft2(gx) + 1j * kyy * np.fft.fft2(gy)))
    return inv_lap(div)


def ring_radius(image: np.ndarray, pixel: float) -> float:
    """Radius (same unit as ``pixel``) of the peak of the azimuthally averaged brightness."""
    img = np.asarray(image, float)
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    r = np.hypot(xx - cx, yy - cy)
    rb = np.round(r).astype(int)
    prof = np.bincount(rb.ravel(), img.ravel()) / np.maximum(np.bincount(rb.ravel()), 1)
    return float(np.argmax(prof[: min(nx, ny) // 2])) * pixel


def spiral_spectrum(
    field: np.ndarray,
    pixel: float,
    r_peak: float,
    *,
    half_width: float = ANNULUS_HALF_UAS,
    n_az: int = N_AZ,
    m_max: int = 10,
) -> dict[int, float]:
    """Azimuthal Fourier amplitudes |c_m| of a complex field, summed over circles in an annulus.

    Circles are every pixel-step radius in ``r_peak +/- half_width`` about the image centre;
    each is sampled at ``n_az`` azimuths by bilinear interpolation (real and imaginary parts
    separately). Returns ``{m: sum_r |c_m(r)|}`` for ``|m| <= m_max``.
    """
    from scipy.ndimage import map_coordinates

    f = np.asarray(field, complex)
    ny, nx = f.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    phi = 2.0 * np.pi * np.arange(n_az) / n_az
    radii = np.arange(max(r_peak - half_width, pixel), r_peak + half_width + 1e-9, pixel)
    acc = np.zeros(2 * m_max + 1)
    for r in radii:
        rp = r / pixel
        # x to the right, y up: azimuth measured counter-clockwise from +x in pixel space
        coords = np.vstack([cy + rp * np.sin(phi), cx + rp * np.cos(phi)])
        re = map_coordinates(f.real, coords, order=1, mode="nearest")
        im = map_coordinates(f.imag, coords, order=1, mode="nearest")
        c = np.fft.fft(re + 1j * im) / n_az  # c_m for e^{+i m phi}
        ms = np.fft.fftfreq(n_az, d=1.0 / n_az).astype(int)
        for m, cm in zip(ms, c, strict=True):
            if abs(m) <= m_max:
                acc[m + m_max] += abs(cm)
    return {m: float(acc[m + m_max]) for m in range(-m_max, m_max + 1)}


def q_from_spectrum(spec: dict[int, float]) -> float:
    """The paper's asymmetry parameter: q = |c_{+1}| / |c_{-1}| (SM p. 8; Table 1 caption)."""
    return float(spec[1] / spec[-1]) if spec[-1] > 0 else float("inf")


def q_estimator(i1: np.ndarray, i2: np.ndarray, pixel: float, dz: float) -> dict[str, float]:
    """The full frozen estimator: TIE phase, field, annulus spiral spectrum, q."""
    i1 = np.asarray(i1, float) / float(np.sum(i1))
    i2 = np.asarray(i2, float) / float(np.sum(i2))
    p = tie_phase(i1, i2, dz, pixel)
    ii = 0.5 * (i1 + i2)
    u = np.sqrt(ii) * np.exp(1j * p)
    r_pk = ring_radius(ii, pixel)
    spec = spiral_spectrum(u, pixel, r_pk)
    return {"q": q_from_spectrum(spec), "r_peak": r_pk, "phase_rms": float(np.std(p))}


def laguerre_gauss_field(n: int, ell: int, *, w: float = 0.25) -> np.ndarray:
    """Test field: a ring-shaped amplitude r^|ell| e^{-r^2/w^2} with helical phase e^{i ell phi}.

    Coordinates are normalised to [-1, 1] across ``n`` pixels; azimuth counter-clockwise from +x,
    matching :func:`spiral_spectrum`, so a pure mode puts all its power at ``m = ell``.
    """
    y, x = np.mgrid[-1.0 : 1.0 : n * 1j, -1.0 : 1.0 : n * 1j]
    r = np.hypot(x, y)
    return (r / w) ** abs(ell) * np.exp(-((r / w) ** 2)) * np.exp(1j * ell * np.arctan2(y, x))


def crescent_model(
    n: int, pixel: float, radius: float, width: float, asym: float, phi0: float
) -> np.ndarray:
    """N1's non-spinning crescent: a Gaussian ring of ``radius``/``width`` (FWHM) with brightness
    1 + ``asym`` cos(phi - ``phi0``). Real and non-negative -- no phase, no spin."""
    y, x = (np.mgrid[0:n, 0:n] - (n - 1) / 2.0) * pixel
    r = np.hypot(x, y)
    phi = np.arctan2(y, x)
    ring = np.exp(-4.0 * np.log(2.0) * ((r - radius) / width) ** 2)
    return ring * (1.0 + asym * np.cos(phi - phi0))
