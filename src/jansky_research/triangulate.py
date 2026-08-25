"""3D triangulation of a solar type III radio source with two spacecraft (STEREO-A + STEREO-B).

The ``swaves`` / ``windwaves`` slices infer a type III beam's heliocentric distance from its emission
frequency, *through a density model* (the observed frequency is the harmonic of the local plasma
frequency; a density model maps that to a radius). This slice gets the distance a completely different
way --- **geometrically**. STEREO/WAVES provides a Level-3 *direction-finding* product (Cecconi et al. 2008; Krupar et al.
2012): at each time and frequency, the goniopolarimetry gives the **direction of arrival** of the
radio emission --- the direction *toward the source* from the spacecraft --- as an azimuth and
colatitude in the heliocentric HEEQ frame, together with each spacecraft's HEEQ position. (The CDF
labels these the "wave-vector" angles; the L3 product stores the arrival/source direction, which is
what we triangulate, and the ``t>0`` forward gate would catch a sign flip.) Two spacecraft, each giving a line of sight to the source, locate it in 3D as the
least-squares intersection of the two rays --- no density model needed.

That makes the two distances **independent**: the geometric radius from triangulation versus the
plasma-frequency radius from the Leblanc model. Comparing them cross-validates both, and the geometry
adds what the drift method cannot --- the source's heliographic **longitude and latitude**, locating
the beam in 3D.

The honest catch is direction-finding noise: a single type III has an apparent source size of tens of
degrees, so per-sample directions scatter by ~10 deg, which over a ~1 AU baseline is tens of R_sun of
positional error. The tooling therefore (i) **intensity-weighted vector-averages** the direction over
a drift-tracking window (scalar angle averaging is wrong near the azimuth wrap), (ii) keeps only
channels with enough good samples and a forward (both-rays-in-front) intersection, and (iii) reports
the per-channel **miss distance** as the consistency diagnostic. Pure NumPy with a synthetic offline
fixture; the real fetch reads the STEREO L3 CDFs (needs the ``windwaves`` extra, ``cdflib``) and is
network-gated.
"""

from __future__ import annotations

import numpy as np

from . import windwaves

__all__ = [
    "RSUN_KM",
    "additive_vs_multiplicative",
    "circular_median_deg",
    "direction_unit",
    "fetch_stereo_df",
    "harmonic_density_grid",
    "mean_direction",
    "noise_bias_calibration",
    "run",
    "synthetic_event",
    "triangulate_rays",
    "triangulate_track",
]

RSUN_KM = 695700.0


def direction_unit(azimuth_deg: np.ndarray, colatitude_deg: np.ndarray) -> np.ndarray:
    """Unit vector(s) from HEEQ wave-vector azimuth and colatitude (degrees).

    Spherical-to-Cartesian with colatitude :math:`\\theta` from :math:`+Z` (the solar rotation axis)
    and azimuth :math:`\\phi` from :math:`+X` (the Sun--Earth meridian) in the XY plane:
    :math:`(\\sin\\theta\\cos\\phi,\\ \\sin\\theta\\sin\\phi,\\ \\cos\\theta)`. Returns an array with a
    trailing length-3 axis.
    """
    th = np.radians(np.asarray(colatitude_deg, float))
    ph = np.radians(np.asarray(azimuth_deg, float))
    return np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], axis=-1)


def mean_direction(
    azimuth_deg: np.ndarray, colatitude_deg: np.ndarray, weight: np.ndarray
) -> tuple[np.ndarray | None, int]:
    """Intensity-weighted mean direction over a set of samples → (unit vector, n_used).

    Each sample's (azimuth, colatitude) is converted to a unit vector and the vectors are averaged with
    ``weight`` (e.g. flux), then renormalised --- the correct way to average directions, immune to the
    azimuth wrap that breaks scalar angle averaging. NaNs (fill values, dropped channels) are ignored.
    Returns ``(None, 0)`` if fewer than three good samples remain or the vectors cancel.
    """
    az = np.asarray(azimuth_deg, float)
    col = np.asarray(colatitude_deg, float)
    w = np.asarray(weight, float)
    good = np.isfinite(az) & np.isfinite(col) & np.isfinite(w) & (w > 0)
    if int(good.sum()) < 3:
        return None, 0
    v = direction_unit(az[good], col[good])
    m = np.sum(v * w[good][:, None], axis=0)
    n = float(np.linalg.norm(m))
    if n == 0.0:
        return None, 0
    return m / n, int(good.sum())


def triangulate_rays(p1: np.ndarray, u1: np.ndarray, p2: np.ndarray, u2: np.ndarray) -> dict:
    """Least-squares intersection of two rays (spacecraft position + direction-to-source).

    Ray *i* is :math:`p_i + t_i\\,u_i`. The pair of parameters minimising the distance between the two
    lines has the closed form below; the source estimate is the midpoint of the shortest segment and
    the **miss distance** is that segment's length (zero for perfectly intersecting rays). ``t1``/``t2``
    are the signed distances along each ray --- both must be positive for the source to lie *in front*
    of both spacecraft. Returns ``source`` (km, HEEQ), ``miss`` (km), ``t1``, ``t2``; ``source`` is NaN
    for (near-)parallel rays.
    """
    p1 = np.asarray(p1, float)
    p2 = np.asarray(p2, float)
    u1 = np.asarray(u1, float)
    u2 = np.asarray(u2, float)
    w0 = p1 - p2
    b = float(u1 @ u2)
    d = float(u1 @ w0)
    e = float(u2 @ w0)
    den = 1.0 - b * b
    nan3 = np.full(3, np.nan)
    if abs(den) < 1e-9:
        return {"source": nan3, "miss": float("nan"), "t1": float("nan"), "t2": float("nan")}
    t1 = (b * e - d) / den
    t2 = (e - b * d) / den
    pa = p1 + t1 * u1
    pb = p2 + t2 * u2
    return {
        "source": 0.5 * (pa + pb),
        "miss": float(np.linalg.norm(pa - pb)),
        "t1": float(t1),
        "t2": float(t2),
    }


def _window(times: np.ndarray, t_center: float, half_s: float) -> np.ndarray:
    return np.where(np.abs(np.asarray(times, float) - t_center) < half_s)[0]


def _burst_center(spec: dict, *, f_lo: float = 0.5, f_hi: float = 10.0) -> float:
    """Time (s) of peak band-integrated flux --- a crude type III burst finder."""
    f = spec["freqs"]
    band = (f >= f_lo) & (f <= f_hi)
    if not band.any():
        band = np.ones_like(f, bool)
    tot = np.nansum(spec["sfu"][:, band], axis=1)
    return float(spec["times"][int(np.nanargmax(tot))])


def triangulate_track(
    spec_a: dict,
    spec_b: dict,
    *,
    t_center: float | None = None,
    half_s: float = 900.0,
    harmonic: int = 2,
    max_miss_rsun: float = 60.0,
    min_samples: int = 5,
) -> dict:
    """Triangulate the source per frequency over a burst window from two spacecraft spectra.

    Each spectrum dict has ``freqs`` (MHz), ``times`` (s), ``az``/``col`` (time × freq, deg, HEEQ),
    ``sfu`` (time × freq flux), and ``pos`` (time × 3, km, HEEQ). For every shared frequency the
    direction is intensity-weighted vector-averaged (:func:`mean_direction`) over the window on each
    spacecraft, the two rays are triangulated (:func:`triangulate_rays`), and channels are **kept** only
    when both rays point forward (``t1,t2 > 0``), the miss distance is below ``max_miss_rsun``, and each
    spacecraft contributed at least ``min_samples`` good samples. Returns per-kept-channel arrays:
    ``freq_mhz``, ``r_geom`` (R_sun, geometric), ``r_plasma`` (R_sun, Leblanc at this harmonic),
    ``miss`` (R_sun), ``lon``/``lat`` (deg, HEEQ), plus the source XYZ and the mean spacecraft
    positions.
    """
    fa = np.asarray(spec_a["freqs"], float)
    fb = np.asarray(spec_b["freqs"], float)
    if t_center is None:
        t_center = _burst_center(spec_a)
    wa = _window(spec_a["times"], t_center, half_s)
    wb = _window(spec_b["times"], t_center, half_s)
    pa = np.asarray(spec_a["pos"], float)[wa].mean(axis=0)
    pb = np.asarray(spec_b["pos"], float)[wb].mean(axis=0)

    freqs, rg, miss, lon, lat, src = [], [], [], [], [], []
    uas, ubs, nas, nbs = [], [], [], []
    for jf, f in enumerate(fa):
        kb = int(np.argmin(np.abs(fb - f)))
        if abs(fb[kb] - f) > 1e-6:
            continue
        ua, na = mean_direction(spec_a["az"][wa, jf], spec_a["col"][wa, jf], spec_a["sfu"][wa, jf])
        ub, nb = mean_direction(spec_b["az"][wb, kb], spec_b["col"][wb, kb], spec_b["sfu"][wb, kb])
        if ua is None or ub is None or na < min_samples or nb < min_samples:
            continue
        tri = triangulate_rays(pa, ua, pb, ub)
        if not (tri["t1"] > 0 and tri["t2"] > 0):
            continue
        m_rsun = tri["miss"] / RSUN_KM
        if not np.isfinite(m_rsun) or m_rsun > max_miss_rsun:
            continue
        s = tri["source"]
        r = float(np.linalg.norm(s)) / RSUN_KM
        freqs.append(float(f))
        rg.append(r)
        miss.append(m_rsun)
        lon.append(float(np.degrees(np.arctan2(s[1], s[0]))))
        lat.append(float(np.degrees(np.arcsin(s[2] / np.linalg.norm(s)))))
        src.append(s)
        uas.append(ua)
        ubs.append(ub)
        nas.append(na)
        nbs.append(nb)

    farr = np.asarray(freqs, float)
    order = np.argsort(farr)
    r_plasma = (
        windwaves.emission_radius(farr[order], harmonic=harmonic) if farr.size else np.array([])
    )
    return {
        "freq_mhz": farr[order],
        "r_geom": np.asarray(rg)[order],
        "r_plasma": np.asarray(r_plasma, float),
        "miss": np.asarray(miss)[order],
        "lon": np.asarray(lon)[order],
        "lat": np.asarray(lat)[order],
        "source_xyz": (np.asarray(src)[order] if len(src) else np.zeros((0, 3))),
        "u_a": (np.asarray(uas)[order] if len(uas) else np.zeros((0, 3))),
        "u_b": (np.asarray(ubs)[order] if len(ubs) else np.zeros((0, 3))),
        "n_a": np.asarray(nas, int)[order] if len(nas) else np.zeros(0, int),
        "n_b": np.asarray(nbs, int)[order] if len(nbs) else np.zeros(0, int),
        "pos_a": pa,
        "pos_b": pb,
    }


def circular_median_deg(angles_deg: np.ndarray) -> tuple[float, float]:
    """Wrap-safe median (and MAD-based scatter) of angles in degrees.

    A scalar ``np.median`` on angles straddling the +-180 branch cut returns a value at the
    edge of one cluster, not a centre -- the exact failure that put the published longitude
    ~10 deg off (168.9 vs ~179). Centre on the circular mean, take the median of the wrapped
    deviations, and add it back. Returns ``(median, scatter)``, median in (-180, 180].
    """
    a = np.radians(np.asarray(angles_deg, float))
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan"), float("nan")
    centre = np.arctan2(np.mean(np.sin(a)), np.mean(np.cos(a)))
    dev = np.degrees(np.angle(np.exp(1j * (a - centre))))
    med = np.degrees(centre) + float(np.median(dev))
    med = ((med + 180.0) % 360.0) - 180.0
    if med == -180.0:
        med = 180.0
    scatter = float(1.4826 * np.median(np.abs(dev - np.median(dev))))
    return med, scatter


def additive_vs_multiplicative(r_geom: np.ndarray, r_plasma: np.ndarray) -> dict:
    """Which one-parameter model describes the geometric-vs-plasma discrepancy: a constant
    additive offset or a constant ratio?

    The published 'factor 2.18' was a median over a ratio that runs monotonically with
    frequency (1.3 to 3.7); the committed channels are far better described by a constant
    additive offset of order the ray-miss distance -- the signature of a constant
    direction-finding bias, not of a density enhancement (which would multiply, i.e. give an
    OLS slope well above 1). Returns both models' parameters and rms residuals, plus the OLS
    slope/intercept of r_geom on r_plasma.
    """
    rg = np.asarray(r_geom, float)
    rp = np.asarray(r_plasma, float)
    ok = np.isfinite(rg) & np.isfinite(rp)
    rg, rp = rg[ok], rp[ok]
    if rg.size < 3:
        return {}
    diff = rg - rp
    ratio = rg / rp
    med_diff = float(np.median(diff))
    med_ratio = float(np.median(ratio))
    slope, intercept = np.polyfit(rp, rg, 1)
    rms_add = float(np.sqrt(np.mean((rg - (rp + med_diff)) ** 2)))
    rms_mul = float(np.sqrt(np.mean((rg - med_ratio * rp) ** 2)))
    return {
        "diff_med_rsun": round(med_diff, 1),
        "diff_mean_rsun": round(float(np.mean(diff)), 1),
        "diff_std_rsun": round(float(np.std(diff, ddof=1)), 1),
        "ratio_med": round(med_ratio, 2),
        "ratio_min": round(float(ratio.min()), 2),
        "ratio_max": round(float(ratio.max()), 2),
        "ols_slope": round(float(slope), 3),
        "ols_intercept_rsun": round(float(intercept), 2),
        "rms_additive_rsun": round(rms_add, 2),
        "rms_multiplicative_rsun": round(rms_mul, 2),
    }


def _block_jackknife_se(values, stat_fn, *, block: int = 6) -> float:
    """Leave-one-contiguous-block-out jackknife SE for channel statistics: adjacent
    frequency channels share a burst window and a source, so single-channel resampling
    understates the error (the rmstructure lesson, along-track)."""
    n = int(np.asarray(values[0]).size) if isinstance(values, tuple) else int(values.size)
    arrs = values if isinstance(values, tuple) else (values,)
    n_blocks = max(int(np.ceil(n / block)), 2)
    vals = []
    for b in range(n_blocks):
        keep = np.ones(n, bool)
        keep[b * block : (b + 1) * block] = False
        if keep.sum() < 3:
            continue
        vals.append(stat_fn(*[a[keep] for a in arrs]))
    v = np.asarray([x for x in vals if np.isfinite(x)], float)
    k = v.size
    if k < 2:
        return float("nan")
    return float(np.sqrt((k - 1) / k * np.sum((v - v.mean()) ** 2)))


def harmonic_density_grid(
    freq_mhz: np.ndarray,
    r_geom: np.ndarray,
    *,
    grid: tuple = ((1, 1.0), (2, 1.0), (1, 4.0), (2, 4.0)),
) -> dict:
    """The (emission harmonic x density scale) systematic grid for the geometric comparison.

    f_p^2 is proportional to n_e, so harmonic=1 at 4x density is EXACTLY harmonic=2 at 1x --
    the same degeneracy the sibling slices commit. Each cell maps the triangulated channels
    through the Leblanc inverse under that assumption and reports the additive offset and
    median ratio. Keys are ``h<harmonic>_s<scale>``.
    """
    f = np.asarray(freq_mhz, float)
    rg = np.asarray(r_geom, float)
    out = {}
    for h, s in grid:
        rp = windwaves.emission_radius(f, harmonic=h, density_scale=s)
        stats = additive_vs_multiplicative(rg, rp)
        out[f"h{h}_s{s:g}"] = {
            "diff_med_rsun": stats.get("diff_med_rsun"),
            "ratio_med": stats.get("ratio_med"),
            "ols_slope": stats.get("ols_slope"),
        }
    return out


def noise_bias_calibration(
    *,
    sep_deg: float,
    lon_deg: float = 179.0,
    noise_degs: tuple = (9.0, 18.0, 25.0),
    seed: int = 5,
) -> dict:
    """Measured (not hand-typed) triangulation bias of the pipeline at the REAL baseline.

    Synthetic events at the real spacecraft separation and near-far-side longitude, with the
    source placed at the unmodified Leblanc radii, quantify how much outward bias pure
    direction-finding scatter produces at each noise level. Keys are the noise levels.
    """
    out = {}
    for nd in noise_degs:
        ev = synthetic_event(
            lon_deg=lon_deg, sep_deg=sep_deg, noise_deg=nd, f_hi_mhz=2.0, seed=seed
        )
        track = triangulate_track(ev["spec_a"], ev["spec_b"], max_miss_rsun=float("inf"))
        stats = additive_vs_multiplicative(track["r_geom"], track["r_plasma"])
        out[f"{nd:g}"] = {
            "ratio_med": stats.get("ratio_med"),
            "diff_med_rsun": stats.get("diff_med_rsun"),
            "n": int(track["freq_mhz"].size),
        }
    return out


#: The miss-distance thresholds the robustness sweep reports (R_sun). 60 is the analysis cut.
MISS_SWEEP_RSUN = (15.0, 30.0, 60.0, 100.0)


def miss_sweep(track: dict, thresholds=MISS_SWEEP_RSUN) -> list[dict]:
    """The shape cross-check and distance ratio as a function of the miss-distance cut.

    Pure filtering: ``triangulate_track`` already returns the per-channel miss distances, so the
    cut is re-applied to the arrays rather than re-triangulating. Run it on a track built with the
    cut OPEN (``max_miss_rsun=inf``); each row then reports ``corr_geom_plasma``,
    ``ratio_geom_plasma`` and ``n`` for channels with miss below that threshold. The paper's
    Method states a 60 R_sun cut that exceeds the smallest inferred distance, so which side of
    the claim survives a tighter cut is exactly what a referee will ask.
    """
    f = np.asarray(track["freq_mhz"], float)
    rg = np.asarray(track["r_geom"], float)
    rp = np.asarray(track["r_plasma"], float)
    miss = np.asarray(track["miss"], float)
    out = []
    for thr in thresholds:
        k = miss <= thr
        n = int(k.sum())
        corr = (
            float(np.corrcoef(rg[k], rp[k])[0, 1])
            if n >= 3 and np.ptp(rg[k]) > 0 and np.ptp(rp[k]) > 0
            else None
        )
        ratio = float(np.median(rg[k] / rp[k])) if n and np.all(rp[k] > 0) else None
        out.append(
            {
                "max_miss_rsun": float(thr),
                "n": n,
                "corr_geom_plasma": round(corr, 3) if corr is not None else None,
                "ratio_geom_plasma": round(ratio, 2) if ratio is not None else None,
                "f_lo_mhz": round(float(f[k].min()), 4) if n else None,
                "f_hi_mhz": round(float(f[k].max()), 3) if n else None,
            }
        )
    return out


def _baseline_separation_deg(pa: np.ndarray, pb: np.ndarray) -> float:
    """Angular separation (deg) of the two spacecraft as seen from the Sun."""
    ca = pa / np.linalg.norm(pa)
    cb = pb / np.linalg.norm(pb)
    return float(np.degrees(np.arccos(np.clip(ca @ cb, -1.0, 1.0))))


def synthetic_event(
    *,
    lon_deg: float = 35.0,
    lat_deg: float = 5.0,
    r0_rsun: float = 2.0,
    speed_c: float = 0.2,
    harmonic: int = 2,
    f_lo_mhz: float = 0.125,
    f_hi_mhz: float = 2.5,
    n_freq: int = 40,
    n_time: int = 60,
    cadence_s: float = 60.0,
    noise_deg: float = 9.0,
    sep_deg: float = 135.0,
    bias_deg_a: float = 0.0,
    t0_offset_b_s: float = 0.0,
    seed: int = 0,
) -> dict:
    """Synthetic two-spacecraft direction-finding event for a radially outflowing type III.

    A beam climbs radially from ``r0_rsun`` at ``speed_c`` × c along heliographic (``lon_deg``,
    ``lat_deg``); the Leblanc density sets the (harmonic) emission frequency at each radius. Two
    spacecraft sit at 1 AU, separated by ``sep_deg`` in longitude and straddling the source, each
    "observing" the true direction to the beam with ``noise_deg`` of Gaussian angular scatter (mimicking
    the wide apparent source size). ``bias_deg_a`` adds a CONSTANT azimuth bias to spacecraft A
    (a direction-finding calibration error the noise budget cannot represent) and
    ``t0_offset_b_s`` shifts spacecraft B's time axis (a per-file epoch origin mismatch) --- the
    two systematic failure modes the round-8 referee found untestable. Returns
    ``spec_a``/``spec_b`` dicts in the same schema as :func:`fetch_stereo_df`, plus ``truth``
    (the injected longitude/latitude and the radius--frequency mapping).
    """
    rng = np.random.default_rng(seed)
    from jansky import solar

    # source unit direction (heliographic lon/lat) and the radial track in frequency
    lam, phi = np.radians(lon_deg), np.radians(lat_deg)
    s_hat = np.array([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)])
    freqs = np.logspace(np.log10(f_hi_mhz), np.log10(f_lo_mhz), n_freq)  # descending
    r_emit = windwaves.emission_radius(freqs, harmonic=harmonic)  # R_sun per frequency

    au = windwaves.R_AU_RSUN * RSUN_KM  # 1 AU in km
    half = np.radians(sep_deg) / 2.0
    pos_a = au * np.array([np.cos(lam + half), np.sin(lam + half), 0.0])
    pos_b = au * np.array([np.cos(lam - half), np.sin(lam - half), 0.0])
    times = np.arange(n_time, dtype=float) * cadence_s

    t_mid = 0.5 * times[-1] if n_time > 1 else 0.0
    sigma_t = 0.22 * times[-1] if n_time > 1 else 1.0

    def _spec(pos: np.ndarray) -> dict:
        az = np.full((n_time, n_freq), np.nan)
        col = np.full((n_time, n_freq), np.nan)
        sfu = np.zeros((n_time, n_freq))
        # each channel emits over a broad central burst, so the window holds enough samples to beat
        # down the per-sample direction noise (as a real, well-observed type III does)
        prof = np.exp(-0.5 * ((times - t_mid) / sigma_t) ** 2)
        on = prof > 0.2
        for jf in range(n_freq):
            src = r_emit[jf] * RSUN_KM * s_hat
            d_true = src - pos
            d_true = d_true / np.linalg.norm(d_true)
            colat0 = np.degrees(np.arccos(np.clip(d_true[2], -1, 1)))
            az0 = np.degrees(np.arctan2(d_true[1], d_true[0]))
            az[on, jf] = az0 + rng.normal(0, noise_deg, int(on.sum()))
            col[on, jf] = colat0 + rng.normal(0, noise_deg, int(on.sum()))
            sfu[on, jf] = 100.0 * prof[on]
        return {
            "freqs": freqs,
            "times": times,
            "az": az,
            "col": col,
            "sfu": sfu,
            "pos": np.tile(pos, (n_time, 1)),
        }

    _ = solar  # density model used via windwaves.emission_radius
    spec_a = _spec(pos_a)
    spec_b = _spec(pos_b)
    if bias_deg_a:
        spec_a["az"] = spec_a["az"] + bias_deg_a
    if t0_offset_b_s:
        spec_b["times"] = spec_b["times"] + t0_offset_b_s
    return {
        "spec_a": spec_a,
        "spec_b": spec_b,
        "truth": {
            "lon_deg": lon_deg,
            "lat_deg": lat_deg,
            "freqs": freqs,
            "r_emit": r_emit,
            "sep_deg": _baseline_separation_deg(pos_a, pos_b),
        },
    }


def fetch_stereo_df(
    date_yyyymmdd: str, *, spacecraft: str = "a"
) -> dict:  # pragma: no cover - network
    """Fetch a STEREO/WAVES Level-3 HFR direction-finding spectrum for one day from SPDF.

    ``spacecraft`` is ``"a"`` (ahead) or ``"b"`` (behind). Reads the HEEQ wave-vector azimuth/colatitude,
    the per-sample flux (SFU), and the spacecraft HEEQ position from the L3 CDF and returns them in the
    schema :func:`triangulate_track` consumes. Needs the ``windwaves`` extra (``cdflib``).
    """
    import re
    import tempfile

    import cdflib
    import requests

    side = "ahead" if spacecraft.lower() == "a" else "behind"
    sc = "sta" if spacecraft.lower() == "a" else "stb"
    yyyy = date_yyyymmdd[:4]
    base = (
        f"https://spdf.gsfc.nasa.gov/pub/data/stereo/{side}/l3/waves/hfr-direction-finding/{yyyy}/"
    )
    idx = requests.get(base, timeout=60).text
    m = re.findall(rf"{sc}_l3_wav_hfr_{date_yyyymmdd}_v[0-9]+\.cdf", idx)
    if not m:
        raise RuntimeError(f"no STEREO-{spacecraft.upper()} L3 DF file for {date_yyyymmdd}")
    raw = requests.get(base + m[0], timeout=300).content
    with tempfile.NamedTemporaryFile(suffix=".cdf") as fh:
        fh.write(raw)
        fh.flush()
        c = cdflib.CDF(fh.name)
        freqs = np.asarray(c.varget("FREQUENCY"), float) / 1e6  # Hz -> MHz
        ep = cdflib.cdfepoch.to_datetime(c.varget("Epoch"))
        # ABSOLUTE time base: seconds since the file date's UTC midnight, NOT since the
        # file's own first record. The daily L3 CDFs are sparse (only records with DF
        # solutions), so ep[0] need not agree between spacecraft, and a per-file origin
        # would silently misalign the two burst windows -- exactly the constant few-degree
        # direction bias this slice is trying to measure.
        day0 = np.datetime64(
            f"{date_yyyymmdd[:4]}-{date_yyyymmdd[4:6]}-{date_yyyymmdd[6:8]}T00:00:00"
        )
        times = (ep - day0) / np.timedelta64(1, "s")
        az = np.asarray(c.varget("WAVE_AZIMUTH_HEEQ"), float)
        col = np.asarray(c.varget("WAVE_COLATITUDE_HEEQ"), float)
        sfu = np.asarray(c.varget("PSD_SFU"), float)
        pos = np.asarray(c.varget("SC_POS_HEEQ"), float)
    for a in (az, col, sfu):
        a[a < -1e30] = np.nan  # CDF fill value
    return {
        "freqs": freqs,
        "times": np.asarray(times, float),
        "az": az,
        "col": col,
        "sfu": sfu,
        "pos": pos,
        "epoch0_utc": str(ep[0]),
        "t0_offset_s": round(float((ep[0] - day0) / np.timedelta64(1, "s")), 1),
    }


def _metrics(track: dict, source: str, harmonic: int, truth: dict | None) -> dict:
    f = track["freq_mhz"]
    rg = track["r_geom"]
    rp = track["r_plasma"]
    n = int(f.size)
    sep = _baseline_separation_deg(track["pos_a"], track["pos_b"])
    corr = (
        float(np.corrcoef(rg, rp)[0, 1]) if n >= 3 and np.ptp(rg) > 0 and np.ptp(rp) > 0 else None
    )
    ratio = float(np.median(rg / rp)) if n and np.all(rp > 0) else None
    lon_med, lon_scatter = circular_median_deg(track["lon"]) if n else (float("nan"), float("nan"))
    m: dict = {
        "source": source,
        "n_tri": n,
        "harmonic": harmonic,
        "sep_deg": round(sep, 1),
        "pos_a_km_heeq": [round(float(v), 1) for v in track["pos_a"]],
        "pos_b_km_heeq": [round(float(v), 1) for v in track["pos_b"]],
        "f_lo_mhz": round(float(f.min()), 4) if n else None,
        "f_hi_mhz": round(float(f.max()), 3) if n else None,
        "r_lo_rsun": round(float(rg.min()), 1) if n else None,
        "r_hi_rsun": round(float(rg.max()), 1) if n else None,
        "r_lo_au": round(float(rg.min()) / windwaves.R_AU_RSUN, 3) if n else None,
        "r_hi_au": round(float(rg.max()) / windwaves.R_AU_RSUN, 3) if n else None,
        "miss_med_rsun": round(float(np.median(track["miss"])), 1) if n else None,
        # circular statistics: a scalar median across the +-180 wrap put the published
        # longitude at the 8th percentile of its own sample
        "lon_med_deg": round(lon_med, 1) if np.isfinite(lon_med) else None,
        "lon_scatter_deg": round(lon_scatter, 1) if np.isfinite(lon_scatter) else None,
        "lat_med_deg": round(float(np.median(track["lat"])), 1) if n else None,
        "lat_scatter_deg": (
            round(float(1.4826 * np.median(np.abs(track["lat"] - np.median(track["lat"])))), 1)
            if n
            else None
        ),
        "corr_geom_plasma": round(corr, 3) if corr is not None else None,
        "ratio_geom_plasma": round(ratio, 2) if ratio is not None else None,
    }
    if n >= 6:
        # the honest comparison: additive vs multiplicative, log-log shape, block-jackknifed
        m.update(additive_vs_multiplicative(rg, rp))
        lf, lg, lp = np.log10(f), np.log10(rg), np.log10(rp)
        m["loglog_corr"] = round(float(np.corrcoef(lg, lp)[0, 1]), 3)
        m["loglog_slope_geom"] = round(float(np.polyfit(lf, lg, 1)[0]), 3)
        m["loglog_slope_plasma"] = round(float(np.polyfit(lf, lp, 1)[0]), 3)
        m["loglog_slope_geom_vs_plasma"] = round(float(np.polyfit(lp, lg, 1)[0]), 3)
        m["corr_jk_se"] = round(
            _block_jackknife_se(
                (rg, rp), lambda a, b: float(np.corrcoef(a, b)[0, 1]) if a.size >= 3 else np.nan
            ),
            3,
        )
        m["diff_med_jk_se_rsun"] = round(
            _block_jackknife_se((rg, rp), lambda a, b: float(np.median(a - b))), 2
        )
        m["ratio_med_jk_se"] = round(
            _block_jackknife_se((rg, rp), lambda a, b: float(np.median(a / b))), 3
        )
        # the implied constant pointing bias: the additive offset over the mean lever arm
        lever_a = np.linalg.norm(track["source_xyz"] - track["pos_a"], axis=1) / RSUN_KM
        lever_b = np.linalg.norm(track["source_xyz"] - track["pos_b"], axis=1) / RSUN_KM
        lever = float(np.median(0.5 * (lever_a + lever_b)))
        m["lever_med_rsun"] = round(lever, 1)
        m["implied_df_bias_deg"] = round(float(np.degrees(m["diff_med_rsun"] / lever)), 1)
        m["harmonic_density_grid"] = harmonic_density_grid(f, rg)
    if truth is not None:
        m["truth_lon_deg"] = truth["lon_deg"]
        m["truth_lat_deg"] = truth["lat_deg"]
        if n:
            lon_err = abs(((m["lon_med_deg"] - truth["lon_deg"] + 180.0) % 360.0) - 180.0)
            m["lon_err_deg"] = round(lon_err, 1)
            m["lat_err_deg"] = round(abs(m["lat_med_deg"] - truth["lat_deg"]), 1)
    return m


def run(
    out: str = ".",
    *,
    offline: bool = True,
    date: str | None = None,
    harmonic: int = 2,
    half_s: float = 900.0,
    max_miss_rsun: float = 60.0,
) -> dict:
    """Full slice: triangulate a type III in 3D from two spacecraft and cross-check the distance."""
    from pathlib import Path

    if offline or date is None:
        ev = synthetic_event(harmonic=harmonic)
        spec_a, spec_b = ev["spec_a"], ev["spec_b"]
        source = "synthetic"
        truth: dict | None = ev["truth"]
        t_center: float | None = None
    else:  # pragma: no cover - network
        spec_a = fetch_stereo_df(date, spacecraft="a")
        spec_b = fetch_stereo_df(date, spacecraft="b")
        source = f"STEREO-A+B L3 DF {date}"
        truth = None
        t_center = None

    track = triangulate_track(
        spec_a,
        spec_b,
        t_center=t_center,
        half_s=half_s,
        harmonic=harmonic,
        max_miss_rsun=max_miss_rsun,
    )
    metrics = _metrics(track, source, harmonic, truth)
    for spec, tag in ((spec_a, "a"), (spec_b, "b")):
        if "epoch0_utc" in spec:  # pragma: no cover - real leg only
            metrics[f"epoch0_utc_{tag}"] = spec["epoch0_utc"]
            metrics[f"t0_offset_s_{tag}"] = spec["t0_offset_s"]
    # The pipeline's own triangulation bias, measured at the real baseline instead of
    # hand-typed in prose: synthetic events at this separation, source AT the Leblanc radii.
    metrics["noise_bias_calibration"] = noise_bias_calibration(sep_deg=metrics["sep_deg"])
    # The sweep runs on a track with the miss cut OPEN, so every threshold filters the same
    # channel set; the analysis cut (60) is one row of it rather than a separate universe.
    open_track = triangulate_track(
        spec_a,
        spec_b,
        t_center=t_center,
        half_s=half_s,
        harmonic=harmonic,
        max_miss_rsun=float("inf"),
    )
    metrics["miss_sweep"] = miss_sweep(open_track)

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    # CSV and figure are committed evidence with no merge machinery of their own: refuse to
    # let a synthetic run overwrite them when the results JSON on disk is real (previously a
    # bare `python -m jansky_research.triangulate` replaced the real 38-channel CSV with
    # synthetic channels while the JSON kept saying "STEREO-A+B" -- a marker that lies).
    json_path = op / "results" / "triangulate_metrics.json"
    write_artifacts = source != "synthetic"
    if not write_artifacts:
        try:
            import json as _json

            from .report import _results_are_real

            write_artifacts = not (
                json_path.is_file() and _results_are_real(_json.loads(json_path.read_text()))
            )
        except Exception:
            write_artifacts = True

    write_results(metrics, json_path)
    if write_artifacts:
        # Per-channel arrays with the cut open, so the 60 R_sun choice is auditable and any
        # future sweep is offline -- plus the mean directions and sample counts, so a reader
        # can reproduce every r_geom from committed evidence (with pos_a/pos_b in the JSON).
        import csv as _csv

        with (op / "results" / "triangulate_channels.csv").open("w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(
                [
                    "freq_mhz",
                    "r_geom_rsun",
                    "r_plasma_rsun",
                    "miss_rsun",
                    "lon_deg",
                    "lat_deg",
                    "ua_x",
                    "ua_y",
                    "ua_z",
                    "ub_x",
                    "ub_y",
                    "ub_z",
                    "n_a",
                    "n_b",
                ]
            )
            for i in range(open_track["freq_mhz"].size):
                row = [
                    round(float(open_track[key][i]), 4)
                    for key in ("freq_mhz", "r_geom", "r_plasma", "miss", "lon", "lat")
                ]
                row += [round(float(v), 6) for v in open_track["u_a"][i]]
                row += [round(float(v), 6) for v in open_track["u_b"][i]]
                row += [int(open_track["n_a"][i]), int(open_track["n_b"][i])]
                w.writerow(row)
        _figure(track, op / "papers" / "triangulate" / "figures")
    _write_macros(metrics, op / "papers" / "triangulate" / "generated" / "macros.tex")
    return metrics


def _figure(track: dict, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    f, rg, rp, miss = track["freq_mhz"], track["r_geom"], track["r_plasma"], track["miss"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.8))

    # Left: independent distance estimators vs frequency
    if f.size:
        ax1.errorbar(f, rg, yerr=miss, fmt="o", color="C0", ms=4, capsize=2, label="geometric")
        ax1.plot(f, rp, "-", color="C3", lw=1.5, label="plasma-frequency (Leblanc)")
    ax1.set(
        xscale="log",
        yscale="log",
        xlabel="frequency (MHz)",
        ylabel=r"heliocentric distance ($R_\odot$)",
        title="Geometric vs plasma distance",
    )
    ax1.legend(fontsize=8)

    # Right: top-down HEEQ ecliptic geometry
    ax2.plot(0, 0, "*", color="orange", ms=14, label="Sun")
    for p, name, c in ((track["pos_a"], "A", "C0"), (track["pos_b"], "B", "C2")):
        ax2.plot(p[0] / RSUN_KM, p[1] / RSUN_KM, "s", color=c, ms=6)
        ax2.annotate(f"STEREO-{name}", (p[0] / RSUN_KM, p[1] / RSUN_KM), fontsize=8)
    if track["source_xyz"].shape[0]:
        s = track["source_xyz"]
        ax2.plot(s[:, 0] / RSUN_KM, s[:, 1] / RSUN_KM, "o", color="C3", ms=3, label="source track")
        for p, c in ((track["pos_a"], "C0"), (track["pos_b"], "C2")):
            ax2.plot(
                [p[0] / RSUN_KM, s[-1, 0] / RSUN_KM],
                [p[1] / RSUN_KM, s[-1, 1] / RSUN_KM],
                "-",
                color=c,
                lw=0.6,
                alpha=0.6,
            )
    ax2.set(
        xlabel=r"$X_{\rm HEEQ}$ ($R_\odot$)", ylabel=r"$Y_{\rm HEEQ}$ ($R_\odot$)", title="Geometry"
    )
    ax2.set_aspect("equal", "datalim")
    ax2.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(out / "triangulate.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.triangulate._write_macros -- do not edit by hand.",
        rf"\newcommand{{\triSource}}{{{m['source']}}}",
        rf"\newcommand{{\triNtri}}{{{_fmt('n_tri')}}}",
        rf"\newcommand{{\triHarmonic}}{{{m['harmonic']}}}",
        rf"\newcommand{{\triSep}}{{{_fmt('sep_deg')}}}",
        rf"\newcommand{{\triFlo}}{{{_fmt('f_lo_mhz')}}}",
        rf"\newcommand{{\triFhi}}{{{_fmt('f_hi_mhz')}}}",
        rf"\newcommand{{\triRlo}}{{{_fmt('r_lo_rsun')}}}",
        rf"\newcommand{{\triRhi}}{{{_fmt('r_hi_rsun')}}}",
        rf"\newcommand{{\triRloAU}}{{{_fmt('r_lo_au')}}}",
        rf"\newcommand{{\triRhiAU}}{{{_fmt('r_hi_au')}}}",
        rf"\newcommand{{\triMiss}}{{{_fmt('miss_med_rsun')}}}",
        rf"\newcommand{{\triLon}}{{{_fmt('lon_med_deg')}}}",
        rf"\newcommand{{\triLonScatter}}{{{_fmt('lon_scatter_deg')}}}",
        rf"\newcommand{{\triLat}}{{{_fmt('lat_med_deg')}}}",
        rf"\newcommand{{\triLatScatter}}{{{_fmt('lat_scatter_deg')}}}",
        rf"\newcommand{{\triCorr}}{{{_fmt('corr_geom_plasma')}}}",
        rf"\newcommand{{\triCorrJkSe}}{{{_fmt('corr_jk_se')}}}",
        rf"\newcommand{{\triCorrLog}}{{{_fmt('loglog_corr')}}}",
        rf"\newcommand{{\triRatio}}{{{_fmt('ratio_geom_plasma')}}}",
        rf"\newcommand{{\triRatioJkSe}}{{{_fmt('ratio_med_jk_se')}}}",
        rf"\newcommand{{\triRatioMin}}{{{_fmt('ratio_min')}}}",
        rf"\newcommand{{\triRatioMax}}{{{_fmt('ratio_max')}}}",
        # the additive description the committed channels actually support
        rf"\newcommand{{\triDiffMed}}{{{_fmt('diff_med_rsun')}}}",
        rf"\newcommand{{\triDiffMean}}{{{_fmt('diff_mean_rsun')}}}",
        rf"\newcommand{{\triDiffStd}}{{{_fmt('diff_std_rsun')}}}",
        rf"\newcommand{{\triDiffJkSe}}{{{_fmt('diff_med_jk_se_rsun')}}}",
        rf"\newcommand{{\triOlsSlope}}{{{_fmt('ols_slope')}}}",
        rf"\newcommand{{\triOlsIntercept}}{{{_fmt('ols_intercept_rsun')}}}",
        rf"\newcommand{{\triRmsAdd}}{{{_fmt('rms_additive_rsun')}}}",
        rf"\newcommand{{\triRmsMul}}{{{_fmt('rms_multiplicative_rsun')}}}",
        rf"\newcommand{{\triSlopeGeom}}{{{_fmt('loglog_slope_geom')}}}",
        rf"\newcommand{{\triSlopePlasma}}{{{_fmt('loglog_slope_plasma')}}}",
        rf"\newcommand{{\triSlopeGvsP}}{{{_fmt('loglog_slope_geom_vs_plasma')}}}",
        rf"\newcommand{{\triLever}}{{{_fmt('lever_med_rsun')}}}",
        rf"\newcommand{{\triBiasDeg}}{{{_fmt('implied_df_bias_deg')}}}",
    ]
    # the (harmonic x density) grid and the measured noise-bias calibration
    grid = m.get("harmonic_density_grid") or {}
    for key, name in (("h1_s1", "HOne"), ("h2_s1", "HTwo"), ("h1_s4", "HOneSFour")):
        cell = grid.get(key) or {}
        lines.append(
            rf"\newcommand{{\triGrid{name}Ratio}}"
            rf"{{{'--' if cell.get('ratio_med') is None else cell['ratio_med']}}}"
        )
        lines.append(
            rf"\newcommand{{\triGrid{name}Diff}}"
            rf"{{{'--' if cell.get('diff_med_rsun') is None else cell['diff_med_rsun']}}}"
        )
    calib = m.get("noise_bias_calibration") or {}
    for nd, name in (("9", "Nine"), ("18", "Eighteen"), ("25", "TwentyFive")):
        cell = calib.get(nd) or {}
        lines.append(
            rf"\newcommand{{\triCalib{name}Ratio}}"
            rf"{{{'--' if cell.get('ratio_med') is None else cell['ratio_med']}}}"
        )
        lines.append(
            rf"\newcommand{{\triCalib{name}Diff}}"
            rf"{{{'--' if cell.get('diff_med_rsun') is None else cell['diff_med_rsun']}}}"
        )
    # The miss-cut sweep, one macro triple per threshold, so the paper's robustness sentence is
    # regenerable rather than hand-typed.
    sweep = {s["max_miss_rsun"]: s for s in m.get("miss_sweep", [])}
    for thr, name in ((15.0, "Fifteen"), (30.0, "Thirty"), (60.0, "Sixty"), (100.0, "Hundred")):
        s = sweep.get(thr, {})
        for key, suffix in (
            ("corr_geom_plasma", "Corr"),
            ("ratio_geom_plasma", "Ratio"),
            ("n", "N"),
        ):
            v = s.get(key)
            lines.append(rf"\newcommand{{\triSweep{name}{suffix}}}{{{'--' if v is None else v}}}")
    lines += [
        rf"\newcommand{{\triTruthLon}}{{{_fmt('truth_lon_deg')}}}",
        rf"\newcommand{{\triTruthLat}}{{{_fmt('truth_lat_deg')}}}",
        rf"\newcommand{{\triLonErr}}{{{_fmt('lon_err_deg')}}}",
        rf"\newcommand{{\triLatErr}}{{{_fmt('lat_err_deg')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and
    # would otherwise blank the other mode's macros with '--'. `make figures`
    # runs every slice offline in the repo root, so without this an offline
    # rebuild silently empties this paper. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(
        description="3D triangulation of a type III source (STEREO-A+B L3 direction-finding)."
    )
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--date", help="YYYYMMDD")
    p.add_argument("--harmonic", type=int, default=2)
    p.add_argument("--half", type=float, default=900.0, help="burst-window half-width (s)")
    p.add_argument("--max-miss", type=float, default=60.0, help="max miss distance (R_sun)")
    args = p.parse_args(argv)
    metrics = run(
        args.out,
        offline=args.offline or not args.date,
        date=args.date,
        harmonic=args.harmonic,
        half_s=args.half,
        max_miss_rsun=args.max_miss,
    )
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
