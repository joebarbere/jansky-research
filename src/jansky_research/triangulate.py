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
    "direction_unit",
    "fetch_stereo_df",
    "mean_direction",
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


def triangulate_rays(
    p1: np.ndarray, u1: np.ndarray, p2: np.ndarray, u2: np.ndarray
) -> dict:
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

    farr = np.asarray(freqs, float)
    order = np.argsort(farr)
    r_plasma = (
        windwaves.emission_radius(farr[order], harmonic=harmonic)
        if farr.size
        else np.array([])
    )
    return {
        "freq_mhz": farr[order],
        "r_geom": np.asarray(rg)[order],
        "r_plasma": np.asarray(r_plasma, float),
        "miss": np.asarray(miss)[order],
        "lon": np.asarray(lon)[order],
        "lat": np.asarray(lat)[order],
        "source_xyz": (np.asarray(src)[order] if len(src) else np.zeros((0, 3))),
        "pos_a": pa,
        "pos_b": pb,
    }


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
    seed: int = 0,
) -> dict:
    """Synthetic two-spacecraft direction-finding event for a radially outflowing type III.

    A beam climbs radially from ``r0_rsun`` at ``speed_c`` × c along heliographic (``lon_deg``,
    ``lat_deg``); the Leblanc density sets the (harmonic) emission frequency at each radius. Two
    spacecraft sit at 1 AU, separated by ``sep_deg`` in longitude and straddling the source, each
    "observing" the true direction to the beam with ``noise_deg`` of Gaussian angular scatter (mimicking
    the wide apparent source size). Returns ``spec_a``/``spec_b`` dicts in the same schema as
    :func:`fetch_stereo_df`, plus ``truth`` (the injected longitude/latitude and the radius--frequency
    mapping).
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
    return {
        "spec_a": _spec(pos_a),
        "spec_b": _spec(pos_b),
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
        f"https://spdf.gsfc.nasa.gov/pub/data/stereo/{side}/l3/waves/"
        f"hfr-direction-finding/{yyyy}/"
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
        times = (ep - ep[0]) / np.timedelta64(1, "s")
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
    m: dict = {
        "source": source,
        "n_tri": n,
        "harmonic": harmonic,
        "sep_deg": round(sep, 1),
        "f_lo_mhz": round(float(f.min()), 4) if n else None,
        "f_hi_mhz": round(float(f.max()), 3) if n else None,
        "r_lo_rsun": round(float(rg.min()), 1) if n else None,
        "r_hi_rsun": round(float(rg.max()), 1) if n else None,
        "r_lo_au": round(float(rg.min()) / windwaves.R_AU_RSUN, 3) if n else None,
        "r_hi_au": round(float(rg.max()) / windwaves.R_AU_RSUN, 3) if n else None,
        "miss_med_rsun": round(float(np.median(track["miss"])), 1) if n else None,
        "lon_med_deg": round(float(np.median(track["lon"])), 1) if n else None,
        "lat_med_deg": round(float(np.median(track["lat"])), 1) if n else None,
        "corr_geom_plasma": round(corr, 3) if corr is not None else None,
        "ratio_geom_plasma": round(ratio, 2) if ratio is not None else None,
    }
    if truth is not None:
        m["truth_lon_deg"] = truth["lon_deg"]
        m["truth_lat_deg"] = truth["lat_deg"]
        if n:
            m["lon_err_deg"] = round(abs(m["lon_med_deg"] - truth["lon_deg"]), 1)
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
    import json
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

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "triangulate_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
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
    ax2.set(xlabel=r"$X_{\rm HEEQ}$ ($R_\odot$)", ylabel=r"$Y_{\rm HEEQ}$ ($R_\odot$)", title="Geometry")
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
        rf"\newcommand{{\triLat}}{{{_fmt('lat_med_deg')}}}",
        rf"\newcommand{{\triCorr}}{{{_fmt('corr_geom_plasma')}}}",
        rf"\newcommand{{\triRatio}}{{{_fmt('ratio_geom_plasma')}}}",
        rf"\newcommand{{\triTruthLon}}{{{_fmt('truth_lon_deg')}}}",
        rf"\newcommand{{\triTruthLat}}{{{_fmt('truth_lat_deg')}}}",
        rf"\newcommand{{\triLonErr}}{{{_fmt('lon_err_deg')}}}",
        rf"\newcommand{{\triLatErr}}{{{_fmt('lat_err_deg')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


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
