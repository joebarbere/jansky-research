"""Stokes-V (circular-polarization) selection of coherent radio emitters in RACS.

Circularly polarized radio emission is a near-unambiguous flag of a *coherent* process
(electron-cyclotron maser, coherent plasma emission). Extragalactic AGN are at most a fraction of a
percent circularly polarized, whereas flaring M dwarfs, ultracool/brown dwarfs, magnetic
chemically-peculiar stars, RS\\,CVn binaries and pulsars reach tens of percent. So a high
$|V|/I$ in a wide-field survey is a clean coherent-emitter finder (Pritchard et al. 2021).

The dominant false positive is instrumental: off-axis **Stokes-I$\\to$V leakage** puts spurious
circular polarization at the $\\sim$1\\% level on bright sources, rising to $\\sim$10\\% near beam
edges. We defeat it with a per-region **leakage floor** estimated from the (assumed unpolarised) field
population --- a candidate must clear $n_\\sigma\\times\\mathrm{median}(|V/I|)$, not merely the image
noise (the RACS convention, $n_\\sigma=7$). Genuine stellar association is then confirmed by **proper
motion**: a real radio star sits at the catalogue position *propagated to the radio epoch*, not at the
optical-epoch position --- a discriminant a chance background match fails.

This composes :mod:`jansky_research.spectra` (cross-match, two-point index) and the
:mod:`jansky_research.vlass` forced-photometry/vetting pattern rather than reimplementing them. Pure
NumPy + a synthetic offline fixture for tests; the real run does forced photometry of a curated
late-type-star / ultracool-dwarf target list in CASDA RACS Stokes-V cutouts.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "LEAKAGE_NSIGMA",
    "RACS_NU_GHZ",
    "classify_emitter",
    "fetch_racs_cutout",
    "fetch_racs_i",
    "fetch_radio_star_measurements",
    "fetch_radio_stars",
    "forced_photometry_recover",
    "fractional_circular_pol",
    "handedness",
    "leakage_floor",
    "match_targets_to_radio",
    "measure_circular_pol",
    "proper_motion_confirm",
    "run",
    "select_circular_pol",
    "synthetic_field",
    "validate_srsc",
]

# RACS band reference frequencies (GHz): RACS-low, RACS-mid, RACS-high.
RACS_NU_GHZ = {"low": 0.8875, "mid": 1.3675, "high": 1.6555}
# Candidate threshold in units of the median field |V/I| (the RACS-low2 / Pritchard convention).
LEAKAGE_NSIGMA = 7.0
# VizieR catalogue IDs. RACS-low DR1 Stokes-I (Hale+2021, split into the extragalactic |b| cut and the
# galactic region) and the Sydney Radio Star Catalogue (Driessen+2024) for the recover-a-known.
VIZIER_RACS_LOW_I = ("J/other/PASA/38.58/galcut", "J/other/PASA/38.58/galreg")
VIZIER_SRSC = "J/other/PASA/41.84/stars"


def fractional_circular_pol(
    v_flux: np.ndarray,
    i_flux: np.ndarray,
    e_v: np.ndarray | None = None,
    e_i: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Fractional circular polarization $|V|/I$ and its 1-$\sigma$ error.

    $f=|V|/I$. The error propagates both flux errors without dividing by the (possibly zero or
    sign-flipping) $V$: $\sigma_f=\sqrt{(e_V/I)^2+(f\,e_I/I)^2}$.
    """
    v = np.asarray(v_flux, float)
    i = np.asarray(i_flux, float)
    frac = np.abs(v) / i
    if e_v is None or e_i is None:
        return frac, np.full_like(frac, np.nan)
    e_v = np.asarray(e_v, float)
    e_i = np.asarray(e_i, float)
    sigma = np.sqrt((e_v / i) ** 2 + (frac * e_i / i) ** 2)
    return frac, sigma


def handedness(v_flux: float) -> str:
    """Sense of circular polarization from the sign of $V$ (``RCP`` if $V>0$, else ``LCP``).

    Note: the absolute $V$ sign convention varies between ASKAP/RACS pipeline versions and epochs, so
    handedness should be cross-checked against the relevant survey paper before physical interpretation.
    """
    return "RCP" if v_flux > 0 else "LCP"


def leakage_floor(frac_pol: np.ndarray, *, n_sigma: float = LEAKAGE_NSIGMA) -> float:
    r"""Per-region leakage threshold from the field $|V/I|$ distribution.

    The field is dominated by unpolarised sources whose apparent circular polarization is
    instrumental Stokes-I$\to$V leakage, so $\mathrm{median}(|V/I|)$ estimates the local leakage
    level (robust to the rare genuine emitters). The credible candidate threshold is
    ``n_sigma`` $\times$ that median --- a candidate must clear *leakage*, not just image noise. Feed
    it the $|V/I|$ of **bright** ($I/\sigma_I\gg1$) sources only: on faint sources $|V/I|\approx
    \sigma_V/I$ is noise-, not leakage-, dominated and would inflate the floor (handle that noise with
    a separate $V$-SNR cut, see :func:`select_circular_pol`).
    """
    f = np.asarray(frac_pol, float)
    f = f[np.isfinite(f)]
    if f.size == 0:
        return float("nan")
    return float(n_sigma * np.median(f))


def select_circular_pol(
    i_flux: np.ndarray,
    v_flux: np.ndarray,
    e_i: np.ndarray,
    e_v: np.ndarray,
    *,
    leakage_threshold: float,
    v_snr_min: float = 5.0,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Select circular-polarization candidates above the leakage floor *and* the $V$ SNR cut.

    A source is a candidate when $|V|/I>$ ``leakage_threshold`` (clears instrumental leakage) **and**
    $|V|/\sigma_V\ge$ ``v_snr_min`` (a real $V$ detection, not noise). Forced photometry of a known
    target can use a looser ``v_snr_min`` (3--4) than a blind 5$\sigma$ extraction. Returns
    ``(mask, frac_pol)``.
    """
    frac, _ = fractional_circular_pol(v_flux, i_flux, e_v, e_i)
    v_snr = np.abs(np.asarray(v_flux, float)) / np.asarray(e_v, float)
    mask = (frac > leakage_threshold) & (v_snr >= v_snr_min)
    return mask, frac


def proper_motion_confirm(
    ra_radio: np.ndarray,
    dec_radio: np.ndarray,
    ra_cat: np.ndarray,
    dec_cat: np.ndarray,
    pmra: np.ndarray,
    pmdec: np.ndarray,
    dt_yr: float,
    *,
    match_arcsec: float = 2.5,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Confirm a stellar association by proper motion.

    Propagates each star's catalogue position to the radio epoch --- ``pmra`` is $\mu_{\alpha*}=
    \mu_\alpha\cos\delta$ (mas/yr), so the RA *coordinate* shift is $\mu_{\alpha*}\,\Delta t/\cos\delta$
    --- and compares it to the measured radio centroid. A genuine radio star matches the *propagated*
    position; a chance background match (or a static source) does not, for an appreciable proper
    motion. Returns ``(confirmed_mask, sep_arcsec)`` where ``sep`` is to the propagated position.
    """
    ra_cat = np.asarray(ra_cat, float)
    dec_cat = np.asarray(dec_cat, float)
    cosd = np.cos(np.radians(dec_cat))
    ra_prop = ra_cat + (np.asarray(pmra, float) * dt_yr / 1000.0 / 3600.0) / cosd
    dec_prop = dec_cat + np.asarray(pmdec, float) * dt_yr / 1000.0 / 3600.0
    dra = (np.asarray(ra_radio, float) - ra_prop) * cosd
    ddec = np.asarray(dec_radio, float) - dec_prop
    sep = np.sqrt(dra**2 + ddec**2) * 3600.0
    return sep <= match_arcsec, sep


def match_targets_to_radio(
    target_ra: np.ndarray,
    target_dec: np.ndarray,
    radio_ra: np.ndarray,
    radio_dec: np.ndarray,
    radio_i: np.ndarray,
    radio_ei: np.ndarray,
    *,
    radius_arcsec: float = 15.0,
) -> dict[str, np.ndarray]:
    """Attach each target star's nearest radio Stokes-I component (reuses ``spectra.crossmatch``).

    Returns per-target arrays: ``matched`` (a radio component within ``radius_arcsec``), the matched
    ``i_flux``/``e_i`` (NaN where unmatched), and ``sep_arcsec``. A target with no match is a radio
    non-detection in this catalogue --- where forced photometry on the images then sets an upper limit.
    """
    from .spectra import crossmatch

    n = np.asarray(target_ra).size
    out = {
        "matched": np.zeros(n, dtype=bool),
        "i_flux": np.full(n, np.nan),
        "e_i": np.full(n, np.nan),
        "sep_arcsec": np.full(n, np.nan),
    }
    if np.asarray(radio_ra).size == 0:
        return out
    it, ir, sep = crossmatch(target_ra, target_dec, radio_ra, radio_dec, radius_arcsec)
    out["matched"][it] = True
    out["i_flux"][it] = np.asarray(radio_i, float)[ir]
    out["e_i"][it] = np.asarray(radio_ei, float)[ir]
    out["sep_arcsec"][it] = sep
    return out


def fetch_radio_stars() -> dict[str, np.ndarray]:  # pragma: no cover - network
    """Fetch the Sydney Radio Star Catalogue (Driessen+2024) from VizieR (the recover-a-known target list).

    Returns positions, Gaia proper motions, the detection ``method`` (which flags Stokes-V detections),
    and the SIMBAD identifier per star.
    """
    import numpy as _np
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "pmRA", "pmDE", "Method", "Survey", "Simbad"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs(VIZIER_SRSC)[0]
    return {
        "ra": _np.asarray(t["RAJ2000"], float),
        "dec": _np.asarray(t["DEJ2000"], float),
        "pmra": _np.asarray(t["pmRA"], float),
        "pmdec": _np.asarray(t["pmDE"], float),
        "method": _np.asarray([str(x) for x in t["Method"]], dtype=object),
        "survey": _np.asarray([str(x) for x in t["Survey"]], dtype=object),
        "simbad": _np.asarray([str(x) for x in t["Simbad"]], dtype=object),
    }


def fetch_racs_i(center, radius_deg: float) -> dict[str, np.ndarray]:  # pragma: no cover - network
    """Cone-search RACS-low DR1 Stokes-I (Hale+2021) on VizieR; returns ra/dec/i_flux/e_i/noise (mJy)."""
    import numpy as _np
    from astropy import units as _u
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "Fpk", "e_Fpk", "Noise"])
    v.ROW_LIMIT = -1
    res = v.query_region(center, radius=radius_deg * _u.deg, catalog=list(VIZIER_RACS_LOW_I))
    if not res:
        return {k: _np.array([]) for k in ("ra", "dec", "i_flux", "e_i", "noise")}
    ra, dec, fpk, efpk, noise = [], [], [], [], []
    for t in res:
        ra.append(_np.asarray(t["RAJ2000"], float))
        dec.append(_np.asarray(t["DEJ2000"], float))
        fpk.append(_np.asarray(t["Fpk"], float))
        efpk.append(_np.asarray(t["e_Fpk"], float))
        noise.append(_np.asarray(t["Noise"], float))
    return {
        "ra": _np.concatenate(ra),
        "dec": _np.concatenate(dec),
        "i_flux": _np.concatenate(fpk),
        "e_i": _np.concatenate(efpk),
        "noise": _np.concatenate(noise),
    }


def measure_circular_pol(
    image_i: np.ndarray,
    image_v: np.ndarray,
    wcs,
    ra: float,
    dec: float,
    *,
    search_arcsec: float = 12.0,
    rms_annulus_arcsec: tuple[float, float] = (30.0, 90.0),
) -> dict[str, float]:
    r"""Forced Stokes-I & V photometry at a locked ``(ra, dec)`` (mirrors ``vlass.measure_image_flux``).

    Finds the Stokes-I peak within ``search_arcsec`` of the target, then reads Stokes V **at that same
    pixel** (V can be either sign) --- the physically correct forced measurement for a point-like
    coherent emitter, where the circular-polarization peak coincides with the total-intensity peak.
    Local RMS for each Stokes comes from an annulus. Returns ``i_peak``, ``v_peak`` (signed),
    ``i_rms``, ``v_rms``, ``frac_pol`` ($|V|/I$), and ``offset_arcsec`` (I-peak vs target). Measuring
    at the locked position reaches *below* the blind extraction threshold and turns a non-detection
    into an honest upper limit, rather than missing the source.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.wcs.utils import proj_plane_pixel_scales

    image_i = np.asarray(image_i, dtype=float)
    image_v = np.asarray(image_v, dtype=float)
    px, py = wcs.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    scale = float(np.mean(proj_plane_pixel_scales(wcs)) * 3600.0)  # arcsec/pixel
    ny, nx = image_i.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(xx - float(px), yy - float(py)) * scale
    nan = float("nan")
    near = (rr <= search_arcsec) & np.isfinite(image_i)
    if not near.any():
        return dict.fromkeys(
            ("i_peak", "v_peak", "i_rms", "v_rms", "frac_pol", "offset_arcsec"), nan
        )
    region = np.where(near, image_i, -np.inf)
    iy, ix = np.unravel_index(int(np.argmax(region)), region.shape)
    i_peak = float(image_i[iy, ix])
    v_peak = float(image_v[iy, ix])  # Stokes V at the I-peak pixel (signed)
    offset = float(np.hypot(ix - float(px), iy - float(py)) * scale)
    ann = (rr > rms_annulus_arcsec[0]) & (rr < rms_annulus_arcsec[1])
    i_ann = image_i[ann & np.isfinite(image_i)]
    v_ann = image_v[ann & np.isfinite(image_v)]
    return {
        "i_peak": i_peak,
        "v_peak": v_peak,
        "i_rms": float(np.std(i_ann)) if i_ann.size > 20 else nan,
        "v_rms": float(np.std(v_ann)) if v_ann.size > 20 else nan,
        "frac_pol": abs(v_peak) / i_peak if i_peak > 0 else nan,
        "offset_arcsec": offset,
    }


def fetch_radio_star_measurements() -> dict[str, np.ndarray]:  # pragma: no cover - network
    """Fetch the SRSC per-detection radio table (Driessen+2024): Stokes I & V peak fluxes per survey.

    Returns ra/dec, ``freq`` (MHz), ``survey``, ``i_flux``/``e_i`` (SpeakI), ``v_flux``/``e_v``
    (SpeakV), and ``rms_v`` (localrmsV) --- the real coherent-emitter measurements used by
    :func:`validate_srsc`.
    """
    import numpy as _np
    from astroquery.vizier import Vizier

    v = Vizier(columns=["*"])
    v.ROW_LIMIT = -1
    cats = v.get_catalogs("J/other/PASA/41.84")
    t = [c for c in cats if c.meta["name"].endswith("/radio")][0]
    return {
        "ra": _np.asarray(t["RAJ2000"], float),
        "dec": _np.asarray(t["DEJ2000"], float),
        "freq": _np.asarray(t["Freq"], float),
        "survey": _np.asarray([str(x) for x in t["Survey"]], dtype=object),
        "i_flux": _np.asarray(t["SpeakI"], float),
        "e_i": _np.asarray(t["e_SpeakI"], float),
        "v_flux": _np.asarray(t["SpeakV"], float),
        "e_v": _np.asarray(t["e_SpeakV"], float),
        "rms_v": _np.asarray(t["localrmsV"], float),
    }


def validate_srsc(*, survey_prefix: str = "RACS") -> dict:  # pragma: no cover - network
    """Recover-a-known: do the SRSC's known V-detected radio stars classify as circular emitters?

    Pulls the SRSC radio measurements (:func:`fetch_radio_star_measurements`), keeps detections from
    surveys whose name starts with ``survey_prefix`` (default the three RACS bands) that have both a
    Stokes I and a Stokes V peak flux, and runs :func:`fractional_circular_pol` /
    :func:`classify_emitter` on them. Coherent emitters should be strongly circularly polarized, so
    the fraction classified ``circular``/``highly_circular`` is the recovery. (This validates the
    selection logic on real data; the leakage-floor + forced-photometry path is exercised separately.)
    """
    import numpy as _np

    m = fetch_radio_star_measurements()
    sel = (
        _np.asarray([str(s).startswith(survey_prefix) for s in m["survey"]])
        & _np.isfinite(m["i_flux"])
        & (m["i_flux"] > 0)
        & _np.isfinite(m["v_flux"])
    )
    frac, _ = fractional_circular_pol(m["v_flux"][sel], m["i_flux"][sel])
    cls = [
        classify_emitter(vv, ii) for vv, ii in zip(m["v_flux"][sel], m["i_flux"][sel], strict=True)
    ]
    n = int(sel.sum())
    circ = sum(c in ("circular", "highly_circular") for c in cls)
    return {
        "n_detections": n,
        "median_frac_pol": float(_np.median(frac)) if n else float("nan"),
        "frac_circular": circ / n if n else 0.0,
        "n_highly_circular": sum(c == "highly_circular" for c in cls),
    }


def classify_emitter(
    v_flux: float, i_flux: float, *, strong: float = 0.3, weak: float = 0.06
) -> str:
    """Coarse circular-polarization class from $|V|/I$.

    ``highly_circular`` ($\\ge$ ``strong``: deep ECME-like polarization, a strong coherent-emitter
    candidate); ``circular`` (``weak``--``strong``); ``weak`` (below ``weak``, leakage-dominated
    regime); ``nan`` for non-finite input.
    """
    if not (np.isfinite(v_flux) and np.isfinite(i_flux)) or i_flux <= 0:
        return "nan"
    frac = abs(v_flux) / i_flux
    if frac >= strong:
        return "highly_circular"
    if frac >= weak:
        return "circular"
    return "weak"


def synthetic_field(
    n_stars: int = 500,
    *,
    pol_fraction: float = 0.05,
    leakage_scale: float = 0.008,
    dt_yr: float = 20.0,
    seed: int = 0,
) -> tuple[dict[str, np.ndarray], float]:
    """Synthetic forced-photometry Stokes-V target sample (offline fixture).

    Models a curated late-type-star target list: each star has a forced $(I,V)$ measurement at the
    radio epoch. A ``pol_fraction`` minority are genuine coherent emitters with deep circular
    polarization ($|V|/I\\sim0.2$--$0.8$); the rest show only instrumental leakage
    ($|V|/I\\sim$ ``leakage_scale``). Emitters have a real radio centroid at the *proper-motion-
    propagated* position; non-emitters' forced centroid scatters about the catalogue position.
    Returns ``(stars, dt_yr)`` where ``stars`` has ra/dec (catalogue epoch), pmra/pmdec, the radio
    centroid ra_radio/dec_radio, fluxes, and the ``is_emitter`` truth mask.
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(40.0, 60.0, n_stars)
    dec = rng.uniform(-40.0, -20.0, n_stars)
    # nearby late-type stars have appreciable proper motions (tens to hundreds of mas/yr)
    pmra = rng.normal(0.0, 150.0, n_stars)
    pmdec = rng.normal(0.0, 150.0, n_stars)
    i_flux = 10.0 ** rng.uniform(0.7, 2.5, n_stars)  # ~5-300 mJy Stokes I
    e_i = 0.05 * i_flux + 0.25  # mJy (~5% + a noise floor)
    e_v = 0.25 * np.ones(n_stars)  # ~0.25 mJy V noise per RACS image

    is_emitter = rng.random(n_stars) < pol_fraction
    sign = rng.choice([-1.0, 1.0], n_stars)
    # leakage component: a half-normal |V/I| centred near leakage_scale for every star
    leak_frac = np.abs(rng.normal(0.0, leakage_scale, n_stars))
    v_flux = sign * leak_frac * i_flux + rng.normal(0.0, e_v)
    # genuine emitters: deep circular polarization replaces the leakage value
    frac_em = rng.uniform(0.2, 0.8, n_stars)
    v_flux[is_emitter] = sign[is_emitter] * frac_em[is_emitter] * i_flux[is_emitter]

    cosd = np.cos(np.radians(dec))
    ra_prop = ra + (pmra * dt_yr / 1000.0 / 3600.0) / cosd
    dec_prop = dec + pmdec * dt_yr / 1000.0 / 3600.0
    jit = lambda s: rng.normal(0.0, s / 3600.0, n_stars)  # noqa: E731
    # emitters sit at the PM-propagated position; non-emitter forced centroids scatter near catalogue
    ra_radio = np.where(is_emitter, ra_prop + jit(0.5), ra + jit(1.0))
    dec_radio = np.where(is_emitter, dec_prop + jit(0.5), dec + jit(1.0))

    stars = {
        "ra": ra,
        "dec": dec,
        "pmra": pmra,
        "pmdec": pmdec,
        "ra_radio": ra_radio,
        "dec_radio": dec_radio,
        "i_flux": i_flux,
        "v_flux": v_flux,
        "e_i": e_i,
        "e_v": e_v,
        "is_emitter": is_emitter,
    }
    return stars, dt_yr


def _casda_session(username: str, pw_path: str):  # pragma: no cover - network
    """A logged-in CASDA session (OPAL). Password comes from ``pw_path`` via a ``getpass`` shim."""
    import getpass
    import pathlib

    from astroquery.casda import Casda

    pw = pathlib.Path(pw_path).expanduser().read_text().strip()
    orig = getpass.getpass
    getpass.getpass = lambda *a, **k: pw  # Casda.login() has no password kwarg
    try:
        casda = Casda()
        casda.login(username=username)
    finally:
        getpass.getpass = orig
    return casda


def _racs_science_mask(table, stokes: str):
    """Boolean mask selecting the RACS **science** restored image for ``stokes`` (i/v).

    The filename must start ``image.<stokes>.`` and be a restored, convolved image --- this excludes
    the ``noiseMap`` / ``meanMap`` products (which also carry ``.i.``/``restored`` and otherwise read as
    a flat ~0.2 mJy noise field).
    """
    import numpy as _np

    fn = _np.array([str(x) for x in table["filename"]])
    return _np.array(
        [f.startswith(f"image.{stokes}.") and "restored" in f and "conv" in f for f in fn]
    )


def fetch_racs_cutout(
    ra: float,
    dec: float,
    *,
    stokes: str = "v",
    radius_deg: float = 0.03,
    casda=None,
    username: str | None = None,
    pw_path: str = "~/.casda_pw",
    retries: int = 3,
):  # pragma: no cover - network
    """Stage and read a RACS (low) Stokes-``i``/``v`` cutout from CASDA → ``(image_mJy, wcs, casda)``.

    Logs in (OPAL), queries the ObsCore images at ``(ra, dec)``, picks the RACS science restored image
    for ``stokes`` (:func:`_racs_science_mask`), stages a SODA cutout, downloads the FITS, and returns
    the image in mJy/beam with its celestial WCS plus the (re)used CASDA session. Retries with a fresh
    login on failure (CASDA intermittently returns HTTP 401 on the datalink step). ``username`` falls
    back to ``$CASDA_USERNAME``. Returns ``None`` if no image or all retries fail.
    """
    import os
    import tempfile

    import astropy.units as _u
    import numpy as _np
    import requests
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.wcs import WCS
    from astroquery.casda import Casda

    username = username or os.environ.get("CASDA_USERNAME")
    if not username:
        raise RuntimeError("set CASDA_USERNAME (OPAL email) for the real Stokes-V cutout fetch")
    coord = SkyCoord(ra * _u.deg, dec * _u.deg)
    for _ in range(retries):
        try:
            if casda is None:
                casda = _casda_session(username, pw_path)
            table = Casda.query_region(coord, radius=0.1 * _u.deg)
            mask = _racs_science_mask(table, stokes)
            if not mask.any():
                return None
            urls = casda.cutout(table[mask][:1], coordinates=coord, radius=radius_deg * _u.deg)
            furl = next(u for u in urls if u.endswith(".fits"))
            raw = requests.get(furl, timeout=200).content
            with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as fh:
                fh.write(raw)
                path = fh.name
            with fits.open(path) as hd:
                data = _np.squeeze(_np.asarray(hd[0].data, float))
                wcs = WCS(hd[0].header).celestial
            os.unlink(path)
            return data * 1000.0, wcs, casda  # Jy/beam -> mJy/beam
        except Exception:
            casda = None  # force a fresh login on retry (handles the intermittent 401)
    return None


def forced_photometry_recover(
    *, max_targets: int = 15, min_i_mjy: float = 3.0, search_arcsec: float = 12.0, username=None
) -> list[dict]:  # pragma: no cover - network
    """Forced photometry of catalogued RACS-LOW Stokes-V emitters in real RACS-low DR1 cutouts.

    Selects the brightest-``I`` SRSC RACS-LOW V-detections (``I>min_i_mjy``, RACS-low sky), stages each
    one's Stokes-I and V science cutouts (:func:`fetch_racs_cutout`), and forced-measures $|V|/I$ at the
    known star position (:func:`measure_circular_pol`). Returns one row per successfully-measured target
    with the catalogue and image $(I, V, |V|/I)$ and the I-peak offset. Coherent stellar emission is
    transient, so a single RACS-low DR1 epoch recovers V only for the subset caught in a polarised state
    --- the honest, variability-limited result.
    """
    import numpy as _np

    m = fetch_radio_star_measurements()
    survey = _np.array([str(s) for s in m["survey"]])
    sel = (
        (survey == "RACS-LOW")
        & _np.isfinite(m["i_flux"])
        & (m["i_flux"] > min_i_mjy)
        & _np.isfinite(m["v_flux"])
        & (m["dec"] < 40.0)
        & (m["dec"] > -85.0)
    )
    idx = _np.where(sel)[0]
    idx = idx[_np.argsort(-m["i_flux"][idx])][:max_targets]
    casda = None
    rows: list[dict] = []
    for i in idx:
        ra, dec = float(m["ra"][i]), float(m["dec"][i])
        gi = fetch_racs_cutout(ra, dec, stokes="i", casda=casda, username=username)
        if gi is None:
            continue
        casda = gi[2]
        gv = fetch_racs_cutout(ra, dec, stokes="v", casda=casda, username=username)
        if gv is None:
            continue
        casda = gv[2]
        meas = measure_circular_pol(gi[0], gv[0], gi[1], ra, dec, search_arcsec=search_arcsec)
        cat_i = float(m["i_flux"][i])
        rows.append(
            {
                "cat_i": cat_i,
                "cat_frac": abs(float(m["v_flux"][i])) / cat_i,
                "img_i": meas["i_peak"],
                "img_v": meas["v_peak"],
                "img_frac": meas["frac_pol"],
                "offset_arcsec": meas["offset_arcsec"],
            }
        )
    return rows


def run(
    out: str = ".", *, offline: bool = True, v_snr_min: float = 5.0, i_snr_ref: float = 10.0
) -> dict:
    """Full slice: forced-photometry V over a target star list, select coherent emitters, write artifacts.

    Offline uses :func:`synthetic_field`; the real path (forced photometry of a curated target list in
    CASDA RACS Stokes-V cutouts) is wired separately. The leakage floor is estimated from the
    **bright** sources only ($I/\\sigma_I>$ ``i_snr_ref``), where $|V/I|$ reflects instrumental leakage
    rather than $V$ noise; the noise on faint targets is handled separately by the ``v_snr_min`` cut.
    Candidates clear the leakage floor *and* the $V$-SNR cut, and are confirmed by proper motion.
    """
    import json
    from pathlib import Path

    op = Path(out)
    # always run the synthetic validation of the selection machinery (fast, no network); the real run
    # additionally does the CASDA forced photometry and merges its macros, so a single `reproduce`
    # build carries both the synthetic-validation and the real forced-photometry numbers.
    metrics = _run_offline(op, v_snr_min=v_snr_min, i_snr_ref=i_snr_ref)
    if not offline:
        metrics = {**metrics, **_run_real(op)}

    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "stokesv_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _write_macros(metrics, op / "papers" / "stokesv" / "generated" / "macros.tex")
    return metrics


def _run_offline(op, *, v_snr_min: float, i_snr_ref: float) -> dict:
    """Synthetic selection path: validate the leakage-floor + V-SNR + PM machinery and its purity."""
    stars, dt_yr = synthetic_field()
    frac, _ = fractional_circular_pol(stars["v_flux"], stars["i_flux"], stars["e_v"], stars["e_i"])
    bright = (stars["i_flux"] / stars["e_i"]) > i_snr_ref  # leakage ~ |V/I| where V noise is small
    threshold = leakage_floor(frac[bright] if bright.any() else frac)
    mask, _ = select_circular_pol(
        stars["i_flux"],
        stars["v_flux"],
        stars["e_i"],
        stars["e_v"],
        leakage_threshold=threshold,
        v_snr_min=v_snr_min,
    )
    pm_ok, _sep = proper_motion_confirm(
        stars["ra_radio"],
        stars["dec_radio"],
        stars["ra"],
        stars["dec"],
        stars["pmra"],
        stars["pmdec"],
        dt_yr,
    )
    candidates = mask & pm_ok
    truth = stars["is_emitter"]
    n_cand = int(candidates.sum())
    _figure(stars, frac, threshold, candidates, op / "papers" / "stokesv" / "figures")
    return {
        "source": "synthetic",
        "n_targets": int(truth.size),
        "n_bright_ref": int(bright.sum()),
        "leakage_floor_pct": round(100.0 * threshold, 3),
        "n_above_floor": int(mask.sum()),
        "n_candidates": n_cand,
        "n_injected": int(truth.sum()),
        "n_recovered": int((candidates & truth).sum()),
        "n_pm_rejected": int((mask & ~pm_ok).sum()),
        "purity": round(float((candidates & truth).sum()) / n_cand, 3) if n_cand else 0.0,
    }


def _run_real(op) -> dict:
    """Real path: forced photometry of RACS-LOW emitters in single-epoch RACS-low DR1 CASDA cutouts."""
    import numpy as _np

    rows = forced_photometry_recover()
    n = len(rows)
    img_i = _np.array([r["img_i"] for r in rows])
    cat_i = _np.array([r["cat_i"] for r in rows])
    img_frac = _np.array([r["img_frac"] for r in rows])
    good_i = _np.isfinite(img_i) & (img_i > 0) & _np.isfinite(cat_i) & (cat_i > 0)
    i_ratio = float(_np.median(img_i[good_i] / cat_i[good_i])) if good_i.any() else float("nan")
    cls = [classify_emitter(r["img_v"], r["img_i"]) for r in rows]
    n_circ = sum(c in ("circular", "highly_circular") for c in cls)
    _real_figure(rows, cls, op / "papers" / "stokesv" / "figures")
    return {
        "source": "RACS-low DR1 (CASDA)",  # overrides the synthetic source in the merged real run
        "n_measured": n,
        "i_recovery_ratio": round(i_ratio, 2),
        "n_v_circular": int(n_circ),
        "frac_v_circular": round(n_circ / n, 3) if n else 0.0,
        "median_img_frac_pct": round(100.0 * float(_np.median(img_frac[_np.isfinite(img_frac)])), 2)
        if _np.isfinite(img_frac).any()
        else None,
    }


def _figure(
    stars: dict, frac: np.ndarray, threshold: float, candidates: np.ndarray, out_dir
) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    i = stars["i_flux"]
    other = ~candidates
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(i[other], frac[other], s=6, color="0.6", label="targets (leakage)")
    ax.scatter(i[candidates], frac[candidates], s=40, color="r", marker="*", label="V candidates")
    ax.axhline(
        threshold, color="b", ls="--", lw=0.8, label=rf"leakage floor ({threshold * 100:.1f}%)"
    )
    ax.set(
        xscale="log",
        yscale="log",
        xlabel=r"Stokes $I$ (mJy)",
        ylabel=r"$|V|/I$",
        title="Circular-polarization selection",
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "circular_pol.pdf")
    plt.close(fig)


def _real_figure(rows: list[dict], cls: list[str], out_dir) -> None:
    from pathlib import Path

    import numpy as _np

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cat_i = _np.array([r["cat_i"] for r in rows])
    img_i = _np.array([r["img_i"] for r in rows])
    img_frac = _np.array([r["img_frac"] for r in rows])
    circ = _np.array([c in ("circular", "highly_circular") for c in cls])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 4.0))
    # Left: forced-photometry I recovered at the known positions (validates the CASDA image pipeline)
    lim = [0.5, max(cat_i.max(), img_i.max(), 1) * 1.3]
    ax1.plot(lim, lim, "k--", lw=0.8, label="1:1")
    ax1.plot(cat_i, img_i, "o", color="C0", ms=5)
    ax1.set(
        xscale="log",
        yscale="log",
        xlabel="catalogue $I$ (mJy)",
        ylabel="forced image $I$ (mJy)",
        title="Stokes $I$ recovered",
        xlim=lim,
        ylim=lim,
    )
    ax1.legend(fontsize=8)
    # Right: image |V|/I per target -- V present only for the subset caught in a polarised state
    order = _np.argsort(-img_frac)
    x = _np.arange(len(rows))
    ax2.bar(x, img_frac[order] * 100, color=["C3" if circ[order][k] else "0.6" for k in x])
    ax2.axhline(6.0, color="b", ls="--", lw=0.8, label="circular threshold (6%)")
    ax2.set(
        xlabel="target (sorted)",
        ylabel=r"image $|V|/I$ (%)",
        title="Single-epoch $V$ (variability-limited)",
    )
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "circular_pol.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _f(key: str, default: str = "--") -> str:
        v = m.get(key)
        return default if v is None else str(v)

    lines = [
        "% Auto-generated by jansky_research.stokesv._write_macros — do not edit by hand.",
        rf"\newcommand{{\svSource}}{{{m['source']}}}",
        rf"\newcommand{{\svNtargets}}{{{_f('n_targets')}}}",
        # synthetic selection-machinery macros
        rf"\newcommand{{\svLeakFloorPct}}{{{_f('leakage_floor_pct')}}}",
        rf"\newcommand{{\svNcandidates}}{{{_f('n_candidates')}}}",
        rf"\newcommand{{\svNinjected}}{{{_f('n_injected')}}}",
        rf"\newcommand{{\svNrecovered}}{{{_f('n_recovered')}}}",
        rf"\newcommand{{\svNpmrejected}}{{{_f('n_pm_rejected')}}}",
        rf"\newcommand{{\svPurity}}{{{_f('purity')}}}",
        # real forced-photometry macros (RACS-low DR1)
        rf"\newcommand{{\svNmeasured}}{{{_f('n_measured')}}}",
        rf"\newcommand{{\svIrec}}{{{_f('i_recovery_ratio')}}}",
        rf"\newcommand{{\svNVcirc}}{{{_f('n_v_circular')}}}",
        rf"\newcommand{{\svFracVcirc}}{{{_f('frac_v_circular')}}}",
        rf"\newcommand{{\svMedImgFrac}}{{{_f('median_img_frac_pct')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Select Stokes-V coherent radio emitters (RACS).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true", help="synthetic run (no network/CASDA)")
    p.add_argument("--v-snr-min", type=float, default=5.0)
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline, v_snr_min=args.v_snr_min), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
