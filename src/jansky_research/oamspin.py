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
    "ring_center",
    "center_on_ring",
    "spiral_spectrum",
    "q_from_spectrum",
    "q_estimator",
    "laguerre_gauss_field",
    "crescent_model",
    "crescent_params",
    "vortex_injection",
    "q_sweep",
    "paper_metrics",
    "paper_figure",
    "write_paper_macros",
    "earth_dz_one_day",
    "patch_pipeline",
    "uas_to_km",
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


def ring_center(
    image: np.ndarray, *, search_px: int = 12, step: float = 0.5, n_bins: int = 36
) -> tuple[float, float]:
    """Pixel centre (y, x) about which the ring's radius is most nearly the same in every
    direction: minimises the scatter of the per-azimuth peak radius over ``n_bins`` sectors.

    Brightness cannot bias this (a lopsided crescent has the same radius on its faint side), unlike
    the image centroid or the peak of the azimuthal average, both of which the bright side pulls.
    RML images are only located up to a shift (closure phases carry no absolute position), so the
    centre must come from the ring's shape. Positions only; no flux threshold.
    """
    from scipy.ndimage import map_coordinates

    img = np.asarray(image, float)
    ny, nx = img.shape
    c0y, c0x = (ny - 1) / 2.0, (nx - 1) / 2.0
    radii = np.arange(2.0, min(nx, ny) / 2.0 - search_px - 1.0, 0.5)
    phis = (np.arange(n_bins) + 0.5) * 2.0 * np.pi / n_bins
    rr, pp = np.meshgrid(radii, phis)
    offs = np.arange(-search_px, search_px + 1e-9, step)
    best, best_c = np.inf, (c0y, c0x)
    for dy in offs:
        for dx in offs:
            ys = c0y + dy + rr * np.sin(pp)
            xs = c0x + dx + rr * np.cos(pp)
            prof = map_coordinates(img, [ys.ravel(), xs.ravel()], order=1).reshape(rr.shape)
            r_pk = radii[np.argmax(prof, axis=1)]
            score = float(np.std(r_pk))
            if score < best:
                best, best_c = score, (c0y + dy, c0x + dx)
    return best_c


def center_on_ring(*images: np.ndarray) -> list[np.ndarray]:
    """Shift every image by ONE common offset so the ring centre of their mean sits on the grid
    centre. A common shift keeps the day-to-day difference that drives the TIE unchanged."""
    from scipy.ndimage import shift as _shift

    mean = np.mean([np.asarray(i, float) for i in images], axis=0)
    cy, cx = ring_center(mean)
    ny, nx = mean.shape
    dy, dx = (ny - 1) / 2.0 - cy, (nx - 1) / 2.0 - cx
    return [np.clip(_shift(np.asarray(i, float), (dy, dx), order=1), 0.0, None) for i in images]


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
    i1, i2 = center_on_ring(i1, i2)
    i1 = i1 / float(np.sum(i1))
    i2 = i2 / float(np.sum(i2))
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


def vortex_injection(
    n: int, pixel: float, radius: float, width: float, asym: float, phi0: float, alpha: float
) -> dict[str, float]:
    """Positive control P0 (post hoc, 2026-10-18): inject a phase vortex P = alpha * phi into a
    crescent through the forward TIE (SM Eq. 11a) and see what the inverse TIE recovers.

    For P = alpha*phi, div(I grad P) = (alpha / r^2) dI/dphi, so the two-plane intensity change is
    proportional to the crescent's azimuthal asymmetry -- zero for a symmetric ring. Returns q of
    the true vortex field and q of the TIE-recovered field. The same dz is used forward and back,
    so no unit convention enters.
    """
    img = crescent_model(n, pixel, radius, width, asym, phi0)
    norm = float(img.sum())
    img = img / norm
    y, x = (np.mgrid[0:n, 0:n] - (n - 1) / 2.0) * pixel
    r = np.hypot(x, y)
    phi = np.arctan2(y, x)
    ring = np.exp(-4.0 * np.log(2.0) * ((r - radius) / width) ** 2) / norm
    didz = (
        -(alpha / (2.0 * np.pi)) * ring * (-asym * np.sin(phi - phi0)) / np.maximum(r, pixel) ** 2
    )
    dz = 1e-3
    p = tie_phase(img - 0.5 * dz * didz, img + 0.5 * dz * didz, dz, pixel)
    r_pk = ring_radius(img, pixel)
    q_true = q_from_spectrum(spiral_spectrum(np.sqrt(img) * np.exp(1j * alpha * phi), pixel, r_pk))
    q_rec = q_from_spectrum(spiral_spectrum(np.sqrt(img) * np.exp(1j * p), pixel, r_pk))
    return {"q_true": q_true, "q_recovered": q_rec, "phase_rms": float(np.std(p))}


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


# ---- EHT data, frozen imaging recipe, simulations (ehtim; the `eht` extra; not unit-tested) ------
EHT2017_URL = (
    "https://data.cyverse.org/dav-anon/iplant/projects/commons_repo/curated/"
    "EHTC_FirstM87Results_Apr2019/EHTC_FirstM87Results_Apr2019_uvfits.tgz"
)
FOV_UAS = 120.0  # the simulation grid (0.48 uas px); radii are physical, so this is not the EHT FOV
TOTAL_FLUX_JY = 0.6
SYS_AMP = 0.02  # Paper III abstract / Sect. 8.4.5: non-closing amplitude systematic
SYS_PHASE_DEG = 1.0  # ... and phase
GAIN_AMP = 0.10  # Paper III Sect. 6.1: a lower limit, used as the null's per-scan gain scatter


def fetch_eht2017(dest="data/oamspin") -> list[str]:  # pragma: no cover - network
    """Download and unpack the EHT 2017 M87* Stokes-I uvfits; return their paths."""
    import tarfile
    from pathlib import Path

    from . import data as _data

    d = Path(dest)
    d.mkdir(parents=True, exist_ok=True)
    tgz = d / EHT2017_URL.rsplit("/", 1)[-1]
    if not tgz.exists():
        _data._download(EHT2017_URL, tgz)
    out = d / "EHTC_FirstM87Results_Apr2019"
    if not out.exists():
        with tarfile.open(tgz) as tf:
            tf.extractall(d, filter="data")
    return sorted(str(p) for p in (out / "uvfits").glob("*.uvfits"))


def uvfits_path(day: int, band: str, root: str = "data/oamspin") -> str:
    """Path of the released file for day-of-year ``day`` (95, 96, 100, 101) and band hi/lo."""
    from pathlib import Path

    return str(
        Path(root)
        / "EHTC_FirstM87Results_Apr2019"
        / "uvfits"
        / f"SR1_M87_2017_{day:03d}_{band}_hops_netcal_StokesI.uvfits"
    )


def pixel_uas() -> float:
    """Pixel size of the 250 x 250 analysis grid in uas."""
    return FOV_UAS / N_GRID


# ---- Amendment A1: the EHT Collaboration's released imaging pipeline --------------------------
EHT_PIPELINE_REPO = "eventhorizontelescope/2019-D01-02"
EHT_PIPELINE_COMMIT = "80d230e76c08edd31548e7cfe61dc7cf94300780"
M87_DISTANCE_MPC = 16.8  # Tamburini et al. 2020, p. 3
PC_KM = 3.0856775814913673e13


def patch_pipeline(source: str) -> str:
    """The one change A1 makes to the EHT script: ``ttype = 'nfft'`` -> ``'fast'``. Exactly one
    line may change, so an upstream edit cannot slip through unnoticed."""
    lines = source.splitlines(keepends=True)
    hits = [k for k, ln in enumerate(lines) if ln.replace(" ", "").startswith("ttype='nfft'")]
    if len(hits) != 1:
        raise ValueError(f"expected one `ttype = 'nfft'` line, found {len(hits)}")
    lines[hits[0]] = "ttype = 'fast'  # oamspin A1: NFFT library not installed\n"
    return "".join(lines)


def fetch_eht_pipeline(dest: str = "data/oamspin/eht2019") -> str:  # pragma: no cover - network
    """Download the pinned pipeline (GPL-3.0; not vendored) and write the A1-patched copy."""
    import urllib.request
    from pathlib import Path

    d = Path(dest)
    d.mkdir(parents=True, exist_ok=True)
    out = d / "pipeline_fast.py"
    if not out.exists():
        url = (
            f"https://raw.githubusercontent.com/{EHT_PIPELINE_REPO}/{EHT_PIPELINE_COMMIT}/"
            "eht-imaging/eht-imaging_pipeline.py"
        )
        with urllib.request.urlopen(url, timeout=60) as r:
            src = r.read().decode()
        (d / "eht-imaging_pipeline.py").write_text(src)
        out.write_text(patch_pipeline(src))
    return str(out)


def run_eht_pipeline(
    infile: str, outfile: str, infile2: str | None = None
) -> None:  # pragma: no cover - subprocess + ehtim
    """Image one uvfits (or two bands together) with the patched EHT pipeline."""
    import subprocess
    import sys

    cmd = [sys.executable, "-I", fetch_eht_pipeline(), "-i", infile, "-o", outfile]
    if infile2:
        cmd += ["-i2", infile2]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def load_fits_image(path: str) -> tuple[np.ndarray, float]:  # pragma: no cover - file I/O
    """FITS image -> (unit-flux array resampled to N_GRID x N_GRID, pixel size in uas)."""
    from astropy.io import fits
    from scipy.ndimage import zoom

    with fits.open(path) as hd:
        arr = np.squeeze(np.asarray(hd[0].data, float))
        cdelt = abs(float(hd[0].header["CDELT1"])) * 3.6e9  # deg -> uas
    n0 = arr.shape[0]
    out = np.clip(zoom(arr, N_GRID / n0, order=1), 0.0, None)
    return out / out.sum(), cdelt * n0 / N_GRID


def uas_to_km(theta_uas: float, distance_mpc: float = M87_DISTANCE_MPC) -> float:
    """Transverse length at the source subtended by ``theta_uas`` (the fiducial x unit)."""
    return theta_uas * 4.84813681109536e-12 * distance_mpc * 1e6 * PC_KM


def earth_dz_one_day(mjd: float) -> float:
    """Fiducial dz (km): change in the Earth's barycentric distance along the M87 line of sight
    over one day starting at ``mjd`` (astropy built-in ephemeris)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord, get_body_barycentric
    from astropy.time import Time

    m87 = SkyCoord(187.70593075 * u.deg, 12.39112331 * u.deg)  # M87, ICRS (SIMBAD; J2000)
    n = m87.cartesian.xyz.value
    t0, t1 = Time(mjd, format="mjd"), Time(mjd + 1.0, format="mjd")
    e0 = get_body_barycentric("earth", t0).xyz.to(u.km).value
    e1 = get_body_barycentric("earth", t1).xyz.to(u.km).value
    return float(abs(np.dot(e1 - e0, n)))


def crescent_params(image: np.ndarray, pixel: float) -> dict[str, float]:
    """Radius, asymmetry A and orientation phi0 of a ring image (for the N1 null model): centre
    on the ring, then the m = 1 Fourier term of the brightness on the peak circle."""
    (img,) = center_on_ring(image)
    r_pk = ring_radius(img, pixel)
    from scipy.ndimage import map_coordinates

    n = img.shape[0]
    c = (n - 1) / 2.0
    phi = 2.0 * np.pi * np.arange(N_AZ) / N_AZ
    ring = map_coordinates(img, [c + r_pk / pixel * np.sin(phi), c + r_pk / pixel * np.cos(phi)])
    c1 = np.mean(ring * np.exp(-1j * phi))
    c0 = np.mean(ring)
    # for 1 + A cos(phi - phi0), c1 = (A/2) e^{-i phi0}
    return {"radius": r_pk, "asym": float(2 * abs(c1) / c0), "phi0": float(-np.angle(c1))}


# x fiducial dz. 0.5/2/10 are the plan's N3. In physical units the TIE phase is ~1e16 rad (fully
# wrapped), so the post hoc decades 1e12-1e18 reach the P ~ 1 rad regime where q is defined.
DZ_SWEEP = (0.5, 1.0, 2.0, 10.0, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18)
PAIRS = ((95, 96), (100, 101))
BANDS = ("lo", "hi")


def q_sweep(i1: np.ndarray, i2: np.ndarray, pixel_uas: float, dz_km: float) -> dict[str, float]:
    """q at the fiducial dz and at every ``DZ_SWEEP`` multiple (x in km at M87, z in km)."""
    px = uas_to_km(pixel_uas)
    return {f"q_x{f:g}": q_estimator(i1, i2, px, dz_km * f)["q"] for f in DZ_SWEEP}


def simulate_uvfits(
    template_uvfits: str,
    model: np.ndarray,
    pixel_uas: float,
    out: str,
    seed: int,
    noise: bool,
    gains: bool = True,
) -> None:  # pragma: no cover - ehtim
    """Observe ``model`` (unit-flux array scaled to TOTAL_FLUX_JY) on the exact coverage of a real
    uvfits. ``noise``: thermal noise, 10% per-scan amplitude gains, random station phases, and the
    2% / 1 deg non-closing systematics of Paper III; without it, a noiseless, calibrated copy."""
    import ehtim as eh

    obs = eh.obsdata.load_uvfits(template_uvfits)
    # ehtim's image rows run north-down relative to our array (round trip measured phi0 -> -phi0),
    # so the model is flipped vertically to put its bright side where the real image's is
    im = eh.image.Image(
        np.flipud(np.asarray(model, float)) * TOTAL_FLUX_JY,
        pixel_uas * eh.RADPERUAS,
        obs.ra,
        obs.dec,
        rf=obs.rf,
        source=obs.source,
        mjd=obs.mjd,
    )
    sim = im.observe_same(
        obs,
        ttype="fast",
        add_th_noise=noise,
        ampcal=not (noise and gains),
        phasecal=not (noise and gains),
        gainp=GAIN_AMP,
        gain_offset=0.0,
        seed=seed,
    )
    if noise:
        rng = np.random.default_rng(seed)
        n = len(sim.data)
        g = (1.0 + SYS_AMP * rng.normal(size=n)) * np.exp(
            1j * np.radians(SYS_PHASE_DEG) * rng.normal(size=n)
        )
        sim.data["vis"] = sim.data["vis"] * g
    sim.save_uvfits(out)


def n1_runner(
    *,
    n_real: int = 200,
    out_jsonl: str = "data/oamspin/n1.jsonl",
    work: str = "data/oamspin/n1",
    gains: bool = True,
) -> int:  # pragma: no cover - ehtim + subprocess, ~1 day
    """N1: the non-spinning crescent (params from results/oamspin_crescent.json) observed on both
    days of each pair and band with independent noise, imaged by the EHT pipeline; resumable."""
    import json
    from pathlib import Path

    cp = json.loads(Path("results/oamspin_crescent.json").read_text())
    pix = FOV_UAS / N_GRID
    model = crescent_model(N_GRID, pix, cp["radius"], cp["width"], cp["asym"], cp["phi0"])
    model = model / model.sum()
    done = set()
    outp = Path(out_jsonl)
    if outp.exists():
        done = {
            (r["pair"], r["band"], r["k"]) for r in map(json.loads, outp.read_text().splitlines())
        }
    Path(work).mkdir(parents=True, exist_ok=True)
    n = 0
    with outp.open("a") as fh:
        for k in range(n_real):
            for a, b in PAIRS:
                for band in BANDS:
                    key = (f"{a}-{b}", band, k)
                    if key in done:
                        continue
                    imgs = []
                    for j, day in enumerate((a, b)):
                        uv = str(Path(work) / f"sim_{a}{b}_{band}_{k}_{j}.uvfits")
                        fi = uv.replace(".uvfits", ".fits")
                        seed = 100_000 * k + 1000 * a + 10 * j + (0 if band == "lo" else 1)
                        simulate_uvfits(uvfits_path(day, band), model, pix, uv, seed, True, gains)
                        run_eht_pipeline(uv, fi)
                        imgs.append(load_fits_image(fi))
                        Path(uv).unlink()
                    (i1, p1), (i2, _) = imgs
                    rec = {
                        "pair": key[0],
                        "band": band,
                        "k": k,
                        **q_sweep(i1, i2, p1, cp["dz_km"][key[0]]),
                    }
                    fh.write(json.dumps(rec) + "\n")
                    fh.flush()
                    n += 1
    return n


DAYS = (95, 96, 100, 101)
EHT3C279_URL = (
    "https://data.cyverse.org/dav-anon/iplant/projects/commons_repo/curated/"
    "EHTC_First3C279Results_May2020/EHTC_First3C279Results_May2020_uvfits.tgz"
)


def fetch_3c279(dest: str = "data/oamspin") -> str:  # pragma: no cover - network
    """Download and unpack the EHT 2017 3C 279 Stokes-I uvfits; return the uvfits directory."""
    import tarfile
    from pathlib import Path

    from . import data as _data

    d = Path(dest)
    tgz = d / EHT3C279_URL.rsplit("/", 1)[-1]
    if not tgz.exists():
        _data._download(EHT3C279_URL, tgz)
    out = d / "EHTC_First3C279Results_May2020"
    if not out.exists():
        with tarfile.open(tgz) as tf:
            tf.extractall(d, filter="data")
    return str(out / "uvfits")


def combined_runner(out_json: str = "results/oamspin_combined.json") -> dict:  # pragma: no cover
    """R1 on band-combined images (the pipeline's own ``-i2``: lo and hi imaged together per day),
    closer to the released fiducial images than single-band imaging."""
    import json
    import time
    from pathlib import Path

    cp = json.loads(Path("results/oamspin_crescent.json").read_text())
    imgdir = Path("data/oamspin/images")
    imgdir.mkdir(parents=True, exist_ok=True)
    imgs = {}
    for day in DAYS:
        f = imgdir / f"{day:03d}_lohi.fits"
        if not f.exists():
            run_eht_pipeline(uvfits_path(day, "lo"), str(f), uvfits_path(day, "hi"))
        imgs[day] = load_fits_image(str(f))
    rows = []
    for a, b in PAIRS:
        (i1, p1), (i2, _) = imgs[a], imgs[b]
        rows.append(
            {"pair": f"{a}-{b}", "band": "lo+hi", **q_sweep(i1, i2, p1, cp["dz_km"][f"{a}-{b}"])}
        )
    out = {
        "source": "real: EHT 2017 M87*, each day imaged with lo+hi combined (-i2) by the EHT "
        f"pipeline {EHT_PIPELINE_REPO}@{EHT_PIPELINE_COMMIT[:7]} (ttype fast), frozen estimator",
        "date": time.strftime("%Y-%m-%d"),
        "R1_combined": rows,
    }
    Path(out_json).write_text(json.dumps(out, indent=1))
    return out


def real_runner(out_json: str = "results/oamspin_real.json") -> dict:  # pragma: no cover - ehtim
    """R1 (real day pairs), N2 (simultaneous lo/hi of each day) and N3 (dz rescaling, inside the
    sweep), all through the A1 pipeline and the frozen estimator. Images are cached under
    ``data/oamspin/images``."""
    import json
    import time
    from pathlib import Path

    cp = json.loads(Path("results/oamspin_crescent.json").read_text())
    imgdir = Path("data/oamspin/images")
    imgdir.mkdir(parents=True, exist_ok=True)
    imgs = {}
    for day in DAYS:
        for band in BANDS:
            f = imgdir / f"{day:03d}_{band}.fits"
            if not f.exists():
                run_eht_pipeline(uvfits_path(day, band), str(f))
            imgs[(day, band)] = load_fits_image(str(f))
    r1, n2 = [], []
    for a, b in PAIRS:
        dz = cp["dz_km"][f"{a}-{b}"]
        for band in BANDS:
            (i1, p1), (i2, _) = imgs[(a, band)], imgs[(b, band)]
            r1.append({"pair": f"{a}-{b}", "band": band, **q_sweep(i1, i2, p1, dz)})
        for day in (a, b):
            (lo, p1), (hi, _) = imgs[(day, "lo")], imgs[(day, "hi")]
            n2.append({"day": day, "pair_dz": f"{a}-{b}", **q_sweep(lo, hi, p1, dz)})
    # N4: 3C 279 through the same chain (the M87 pipeline's fixed zbl/prior/FOV are not tuned for
    # a bright jet -- run as frozen, reported with that caveat; q has no spin meaning here)
    n4: list[dict[str, object]] = []
    uvd = Path(fetch_3c279())
    for a, b in PAIRS:
        ii = []
        for day in (a, b):
            f = imgdir / f"3c279_{day:03d}_lo.fits"
            if not f.exists():
                try:
                    run_eht_pipeline(
                        str(uvd / f"SR1_3C279_2017_{day:03d}_lo_hops_netcal_StokesI.uvfits"), str(f)
                    )
                except Exception as exc:
                    n4.append({"pair": f"{a}-{b}", "error": f"{type(exc).__name__}: {exc}"[:300]})
                    break
            ii.append(load_fits_image(str(f)))
        if len(ii) == 2:
            (i1, p1), (i2, _) = ii
            n4.append(
                {"pair": f"{a}-{b}", "band": "lo", **q_sweep(i1, i2, p1, cp["dz_km"][f"{a}-{b}"])}
            )
    out = {
        "source": "real: EHT 2017 M87* uvfits (CyVerse doi:10.25739/g85n-f134), EHT eht-imaging "
        f"pipeline {EHT_PIPELINE_REPO}@{EHT_PIPELINE_COMMIT[:7]} (ttype fast), frozen estimator",
        "date": time.strftime("%Y-%m-%d"),
        "paper_q": {"q1": 1.417, "q2": 1.369, "mean": 1.393},
        "R1_day_pairs": r1,
        "N2_simultaneous_bands": n2,
        "N4_3c279": n4,
    }
    Path(out_json).write_text(json.dumps(out, indent=1))
    return out


def paper_metrics(results_dir: str = "results") -> dict:
    """Flat metrics for the note, read from the committed results files only."""
    import json
    from pathlib import Path

    rd = Path(results_dir)
    m: dict = {}

    def load(name: str):
        p = rd / f"oamspin_{name}.json"
        return json.loads(p.read_text()) if p.exists() else None

    cr = load("crescent")
    if cr:
        m.update(
            real_ring_radius=cr["radius"],
            real_ring_asym=cr["asym"],
            real_dz_km=round(cr["dz_km"]["95-96"]),
        )
    ff = load("falsefailure")
    if ff:
        m.update(syn_ff_q=round(ff["q_x1"], 2), syn_ff_q_pone=round(ff["q_x1e+16"], 2))
    n1 = load("n1")
    if n1:
        s = n1["summary"]
        rows = n1["realisations"]
        m.update(
            syn_n1_n=n1["n"],
            syn_n1_lo_fid=round(s["q_x1"]["p2.5"], 2),
            syn_n1_hi_fid=round(s["q_x1"]["p97.5"], 2),
            syn_n1_lo_pone=round(s["q_x1e+16"]["p2.5"], 2),
            syn_n1_hi_pone=round(s["q_x1e+16"]["p97.5"], 2),
            syn_n1_med_pone=round(s["q_x1e+16"]["median"], 2),
            syn_n1_exceed_pone=round(
                100 * float(np.mean([r["q_x1e+16"] >= 1.393 for r in rows])), 1
            ),
        )
    cal = load("n1cal")
    if cal:
        s = cal["summary"]["q_x1e+16"]
        m.update(
            syn_cal_n=cal["n"],
            syn_cal_lo_pone=round(s["p2.5"], 2),
            syn_cal_hi_pone=round(s["p97.5"], 2),
        )
    real = load("real")
    if real:
        r1 = {(x["pair"], x["band"]): x for x in real["R1_day_pairs"]}
        m.update(
            real_r1_a_lo=round(r1[("95-96", "lo")]["q_x1"], 2),
            real_r1_a_hi=round(r1[("95-96", "hi")]["q_x1"], 2),
            real_r1_b_lo=round(r1[("100-101", "lo")]["q_x1"], 2),
            real_r1_b_hi=round(r1[("100-101", "hi")]["q_x1"], 2),
        )
    comb = load("combined")
    if comb:
        c = {x["pair"]: x for x in comb["R1_combined"]}
        m.update(
            real_comb_a_pone=round(c["95-96"]["q_x1e+16"], 2),
            real_comb_b_pone=round(c["100-101"]["q_x1e+16"], 2),
        )
    p0s = vortex_injection(N_GRID, 0.512, 19.97, 15.87, 0.0, 0.0, 0.3)
    p0c = vortex_injection(N_GRID, 0.512, 19.97, 15.87, 0.535, -2.15, 1.0)
    m.update(
        syn_pzero_sym_true=round(p0s["q_true"], 2),
        syn_pzero_sym_rec=round(p0s["q_recovered"], 2),
        syn_pzero_cres_true=round(p0c["q_true"]),
        syn_pzero_cres_rec=round(p0c["q_recovered"], 2),
    )
    return m


_MACROS = {
    "RealRingRadius": "real_ring_radius",
    "RealDz": "real_dz_km",
    "SynFfQ": "syn_ff_q",
    "SynNoneN": "syn_n1_n",
    "SynNoneLoFid": "syn_n1_lo_fid",
    "SynNoneHiFid": "syn_n1_hi_fid",
    "SynNoneLoPone": "syn_n1_lo_pone",
    "SynNoneHiPone": "syn_n1_hi_pone",
    "SynNoneMedPone": "syn_n1_med_pone",
    "SynNoneExceedPone": "syn_n1_exceed_pone",
    "SynCalN": "syn_cal_n",
    "SynCalLoPone": "syn_cal_lo_pone",
    "SynCalHiPone": "syn_cal_hi_pone",
    "RealROneALo": "real_r1_a_lo",
    "RealROneAHi": "real_r1_a_hi",
    "RealROneBLo": "real_r1_b_lo",
    "RealROneBHi": "real_r1_b_hi",
    "RealCombAPone": "real_comb_a_pone",
    "RealCombBPone": "real_comb_b_pone",
    "SynPzeroSymTrue": "syn_pzero_sym_true",
    "SynPzeroSymRec": "syn_pzero_sym_rec",
    "SynPzeroCresTrue": "syn_pzero_cres_true",
    "SynPzeroCresRec": "syn_pzero_cres_rec",
}


def write_paper_macros(
    path: str = "papers/oamspin/generated/macros.tex", results_dir: str = "results"
) -> None:
    """Macros for the note from the committed results (merged through preserve_live_macros)."""
    from pathlib import Path

    from .report import preserve_live_macros

    m = paper_metrics(results_dir)
    lines = [
        "% Auto-generated by jansky_research.oamspin.write_paper_macros -- do not edit.",
        r"\newcommand{\oamspinSource}{real EHT 2017 data + simulations on its coverage}",
    ]
    for name, key in _MACROS.items():
        v = m.get(key)
        lines.append(rf"\newcommand{{\oamspin{name}}}{{{'--' if v is None else v}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def paper_figure(
    results_dir: str = "results", out: str = "papers/oamspin/figures/oamspin_q.pdf"
) -> None:
    """The note's one exhibit: q vs Delta-z for the non-spinning null (median, 95% band; with and
    without station gains), the real day pairs, and the paper's 1.393."""
    import json
    from pathlib import Path

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rd = Path(results_dir)
    n1 = json.loads((rd / "oamspin_n1.json").read_text())["summary"]
    cal = json.loads((rd / "oamspin_n1cal.json").read_text())["summary"]
    real = json.loads((rd / "oamspin_real.json").read_text())["R1_day_pairs"]
    comb = json.loads((rd / "oamspin_combined.json").read_text())["R1_combined"]
    keys = [f"q_x{f:g}" for f in DZ_SWEEP]
    xs = np.array(DZ_SWEEP)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    # error bars at the sampled dz only (no interpolation across the unsampled 10 - 1e12 gap)
    for summ, col, lab, dx in (
        (n1, "0.45", "null, 10% gains", 0.82),
        (cal, "tab:blue", "null, no gains", 1.22),
    ):
        med = np.array([summ[k]["median"] for k in keys])
        lo = np.array([summ[k]["p2.5"] for k in keys])
        hi = np.array([summ[k]["p97.5"] for k in keys])
        ax.errorbar(
            xs * dx,
            med,
            yerr=[med - lo, hi - med],
            fmt="_",
            color=col,
            capsize=2,
            lw=1,
            label=f"{lab}: median, 95%",
        )
    for x in real:
        ax.plot(xs, [x[k] for k in keys], "o", ms=3, color="k", alpha=0.6)
    for x in comb:
        ax.plot(xs, [x[k] for k in keys], "s", ms=3.5, mfc="none", color="tab:red")
    ax.plot([], [], "o", ms=3, color="k", label="M87* day pairs, single band")
    ax.plot(
        [], [], "s", ms=3.5, mfc="none", color="tab:red", label="M87* day pairs, bands combined"
    )
    ax.axhline(1.393, color="tab:orange", ls="--", lw=1, label="Tamburini et al. (2020)")
    ax.axvspan(3e13, 3e16, color="tab:green", alpha=0.08)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta z$ / fiducial")
    ax.set_ylabel(r"$q = |c_{+1}|/|c_{-1}|$")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
