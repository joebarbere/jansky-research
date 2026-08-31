"""Milky Way HI rotation curve via the tangent-point method.

For an inner-Galaxy sightline ($0 < \\ell < 90°$), neutral hydrogen at the *tangent point* —
galactocentric radius $R = R_0\\sin\\ell$ — moves fastest along the line of sight, producing the
**terminal velocity** $v_\\mathrm{term}$, the high-velocity edge of the HI 21 cm profile. The
circular rotation speed there is $V(R) = v_\\mathrm{term} + V_0\\sin\\ell$. Sweeping longitude
traces the rotation curve $V(R)$ — which comes out *flat*, the textbook signature of dark matter.

This module reads the Leiden/Argentine/Bonn (LAB) HI survey $(b, v)$ slices (one per longitude;
Kalberla et al. 2005), extracts the terminal velocity, and builds the curve. Pure NumPy + astropy;
a synthetic $(\\ell, v)$ slice with a known injected curve lets the tests run offline.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "R0_KPC",
    "V0_KMS",
    "compare_terminal_velocities",
    "fetch_lab_longitude",
    "fetch_mgd2016",
    "read_lab_slice",
    "rotation_curve",
    "run",
    "synthetic_lv_slice",
    "synthetic_reference_curve",
    "tangent_point",
    "terminal_velocity",
    "terminal_velocity_edge",
]

R0_KPC = 8.15  # Sun's galactocentric radius (Reid et al. 2019)
V0_KMS = 236.0  # circular rotation speed at the Sun (Reid et al. 2019)


def read_lab_slice(path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a LAB $(b, v)$ FITS slice; return ``(lat_deg, vel_kms, data[lat, vel])``."""
    from astropy.io import fits

    with fits.open(path) as hd:
        h = hd[0].header
        d = np.asarray(hd[0].data, dtype=float).squeeze()  # (lat, vel)
    nb, nv = d.shape
    vel = (np.arange(nv) + 1 - h["CRPIX1"]) * h["CDELT1"] + h["CRVAL1"]
    if str(h.get("CUNIT1", "M/S")).strip().upper() in ("M/S", "M S-1", ""):
        vel = vel / 1000.0  # -> km/s (LAB stores VELO-LSR in m/s)
    if np.nanmax(np.abs(vel)) < 50.0:  # guard against a misread velocity unit
        raise ValueError(f"implausible velocity axis (max |v| = {np.nanmax(np.abs(vel)):.2g} km/s)")
    lat = (np.arange(nb) + 1 - h["CRPIX2"]) * h["CDELT2"] + h["CRVAL2"]
    return lat, vel, d


def terminal_velocity(
    vel_kms: np.ndarray, spectrum: np.ndarray, *, threshold_k: float = 2.0
) -> float:
    """Terminal velocity: the most positive LSR velocity with $T_B$ above ``threshold_k`` (inner Galaxy).

    A fixed brightness-temperature threshold is the simple, standard estimator. Because the
    profile edge has a finite width, the threshold crossing sits out on the wing and
    **overestimates** the terminal velocity; :func:`terminal_velocity_edge` fits the edge itself,
    and :func:`run` measures the offset between the two on the same data rather than citing one.
    """
    vel_kms = np.asarray(vel_kms, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    above = np.isfinite(spectrum) & (spectrum > threshold_k)
    if not above.any():
        return float("nan")
    return float(np.max(vel_kms[above]))


def terminal_velocity_edge(
    vel_kms: np.ndarray, spectrum: np.ndarray, *, threshold_k: float = 2.0
) -> tuple[float, float]:
    """Terminal velocity from an error-function fit to the profile's high-velocity edge.

    Fits $T(v) = \\frac{A}{2}\\,\\mathrm{erfc}\\!\\big((v-v_0)/\\sqrt{2}w\\big) + c$ over a window
    around the threshold crossing and returns ``(v0, w)`` — the edge's half-height velocity and
    its Gaussian width. The half-height point is interior to the edge, so it does not carry the
    wing overshoot of the threshold crossing; the same construction underlies the edge fits of
    McClure-Griffiths & Dickey (2016). Returns ``(nan, nan)`` if no crossing exists or the fit
    fails.
    """
    from scipy.optimize import curve_fit
    from scipy.special import erfc

    vel_kms = np.asarray(vel_kms, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    good = np.isfinite(spectrum)
    vel_kms, spectrum = vel_kms[good], spectrum[good]
    order = np.argsort(vel_kms)
    vel_kms, spectrum = vel_kms[order], spectrum[order]
    v_thr = terminal_velocity(vel_kms, spectrum, threshold_k=threshold_k)
    if not np.isfinite(v_thr):
        return float("nan"), float("nan")
    m = (vel_kms > v_thr - 50.0) & (vel_kms < v_thr + 40.0)
    x, y = vel_kms[m], spectrum[m]
    if x.size < 8:
        return float("nan"), float("nan")

    def f(v, amp, v0, w, c):
        return amp / 2.0 * erfc((v - v0) / (np.sqrt(2.0) * abs(w) + 1e-3)) + c

    p0 = [max(float(np.interp(v_thr - 30.0, x, y)), 5.0), v_thr - 10.0, 8.0, 0.0]
    try:
        p, _ = curve_fit(f, x, y, p0=p0, maxfev=40000)
    except Exception:
        return float("nan"), float("nan")
    return float(p[1]), float(abs(p[2]))


def tangent_point(
    l_deg: float, v_term: float, *, R0: float = R0_KPC, V0: float = V0_KMS
) -> tuple[float, float]:
    """Tangent-point $(R, V)$: $R = R_0\\sin\\ell$ (kpc), $V = v_\\mathrm{term} + V_0\\sin\\ell$ (km/s)."""
    s = np.sin(np.radians(l_deg))
    return float(R0 * s), float(v_term + V0 * s)


def rotation_curve(
    longitudes: np.ndarray, slices, *, threshold_k: float = 2.0, estimator: str = "threshold"
) -> tuple[np.ndarray, np.ndarray]:
    """Build the rotation curve from per-longitude $(b, v)$ slices (uses each $b=0$ spectrum).

    ``slices`` is an iterable of ``(lat_deg, vel_kms, data)`` aligned with ``longitudes``.
    ``estimator`` selects the terminal-velocity estimator: ``"threshold"`` (fixed 2 K crossing)
    or ``"edge"`` (error-function edge fit). Returns ``(R_kpc, V_kms)`` sorted by radius.
    """
    rad, vel = [], []
    for ell, (lat, v, d) in zip(longitudes, slices, strict=True):
        spec = d[int(np.argmin(np.abs(lat)))]  # b = 0
        if estimator == "edge":
            vt, _ = terminal_velocity_edge(v, spec, threshold_k=threshold_k)
        else:
            vt = terminal_velocity(v, spec, threshold_k=threshold_k)
        r, vv = tangent_point(ell, vt)
        rad.append(r)
        vel.append(vv)
    order = np.argsort(rad)
    return np.asarray(rad)[order], np.asarray(vel)[order]


def synthetic_lv_slice(
    l_deg: float,
    *,
    v_flat: float = 230.0,
    R0: float = R0_KPC,
    V0: float = V0_KMS,
    n_lat: int = 21,
    noise_k: float = 0.5,
    edge_width_kms: float = 5.0,
    seed: int | None = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Synthetic LAB-like $(b, v)$ slice with a known *flat* rotation curve (offline fixture).

    Injects HI emission out to the terminal velocity implied by a flat curve $V=$ ``v_flat`` —
    $v_\\mathrm{term} = v_\\mathrm{flat} - V_0\\sin\\ell$ — with a high-velocity edge of realistic
    width ``edge_width_kms`` (real LAB edges are several km/s wide) plus noise. With a finite
    edge width the fixture *measures* the threshold estimator's wing overshoot instead of hiding
    it: :func:`terminal_velocity` lands $\\sim w\\ln(A/T_\\mathrm{thr})$ beyond the half-height
    point, while :func:`terminal_velocity_edge` recovers the injected edge location. An earlier
    fixture used a 0.6 km/s edge, for which threshold $\\approx$ half-height by construction and
    the dominant real-data systematic was invisible offline.
    """
    rng = np.random.default_rng(seed)
    s = np.sin(np.radians(l_deg))
    v_term = v_flat - V0 * s
    vel = np.linspace(-50.0, 320.0, 400)
    lat = np.linspace(-5.0, 5.0, n_lat)
    # bright HI extending to low/negative velocity with a falling edge at the terminal
    # velocity (which is ~0 at high longitudes, where R_tan -> R0 and V(R_tan) -> V0).
    profile = 30.0 / (1.0 + np.exp((vel - v_term) / edge_width_kms))
    data = profile[None, :] * np.exp(-0.5 * (lat[:, None] / 3.0) ** 2)  # peak at b=0
    data = data + rng.normal(0.0, noise_k, size=data.shape)
    return lat, vel, data


def compare_terminal_velocities(
    longitudes: np.ndarray,
    v_term: np.ndarray,
    ref_l_deg: np.ndarray,
    ref_v_kms: np.ndarray,
    *,
    half_width_deg: float = 0.5,
    min_points: int = 5,
) -> dict:
    """Match a terminal-velocity curve against a densely sampled external reference curve.

    The reference (McClure-Griffiths & Dickey 2016, from the VGPS) is sampled every $0.065\\arcdeg$,
    far finer than the LAB beam, so each of our longitudes is compared against the *mean*
    reference value within ``half_width_deg`` — one LAB half-beam. Longitudes with fewer than
    ``min_points`` reference samples in that window (outside the reference's coverage) are
    dropped rather than extrapolated.

    The comparison is made on the terminal velocity itself, not on $V(R)$: $v_\\mathrm{term}$ is a
    directly observed LSR quantity, so no $R_0$/$V_0$ choice enters and the two surveys are
    compared on what each actually measured.

    Returns the per-longitude match plus ``mean``/``sd``/``sem``/``n`` of the difference.
    """
    longitudes = np.asarray(longitudes, dtype=float)
    v_term = np.asarray(v_term, dtype=float)
    ref_l_deg = np.asarray(ref_l_deg, dtype=float)
    ref_v_kms = np.asarray(ref_v_kms, dtype=float)
    matched_l, matched_ref, matched_v = [], [], []
    for ell, vt in zip(longitudes, v_term, strict=True):
        m = np.abs(ref_l_deg - ell) <= half_width_deg
        if int(m.sum()) < min_points or not np.isfinite(vt):
            continue
        matched_l.append(float(ell))
        matched_ref.append(float(np.mean(ref_v_kms[m])))
        matched_v.append(float(vt))
    d = np.asarray(matched_v) - np.asarray(matched_ref)
    n = int(d.size)
    sd = float(np.std(d, ddof=1)) if n > 1 else float("nan")
    return {
        "longitudes_deg": matched_l,
        "reference_v_term_kms": matched_ref,
        "v_term_kms": matched_v,
        "delta_kms": d.tolist(),
        "mean": float(np.mean(d)) if n else float("nan"),
        "median": float(np.median(d)) if n else float("nan"),
        "sd": sd,
        "sem": sd / np.sqrt(n) if n > 1 else float("nan"),
        "n": n,
        "half_width_deg": half_width_deg,
    }


def synthetic_reference_curve(longitudes: np.ndarray, *, v_flat: float = 230.0) -> tuple:
    """Densely sampled *known* terminal-velocity curve for the offline fixture.

    Mirrors what :func:`fetch_mgd2016` supplies on a real run, so the offline leg exercises the
    same comparison code path with a reference whose answer is known exactly: the injected
    curve. The edge estimator should match it, the threshold estimator should sit above it by
    the wing overshoot.
    """
    lo, hi_ = float(np.min(longitudes)) + 1.0, float(np.max(longitudes)) - 1.0
    ell = np.arange(lo, hi_ + 1e-9, 0.065)
    return ell, v_flat - V0_KMS * np.sin(np.radians(ell))


def fetch_mgd2016() -> tuple:  # pragma: no cover - network
    """McClure-Griffiths & Dickey (2016) Table 1: VGPS HI terminal velocities (cached).

    748 rows over $18.4\\arcdeg < \\ell < 67.0\\arcdeg$ at $0.065\\arcdeg$ spacing, from VizieR
    ``J/ApJ/831/124/table1``. Their $v_\\mathrm{LSR}$ is the *fitted* terminal velocity: a sum of
    two error functions seeded from a 20 K threshold crossing, on continuum-masked, latitude-
    averaged spectra. Returns ``(l_deg, v_term_kms)``.
    """
    from . import data as _data

    target = _data.data_dir() / "mgd2016_table1.tsv"
    if not target.exists():
        _data._download(
            "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/831/124/table1"
            "&-out=GLON,vLSR&-out.max=5000",
            target,
        )
    return read_mgd2016(target)


def read_mgd2016(path) -> tuple:
    """Parse a cached VizieR TSV of MG&D 2016 Table 1 into ``(l_deg, v_term_kms)``."""
    from pathlib import Path

    ell, vel = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        try:
            ell.append(float(parts[0]))
            vel.append(float(parts[1]))
        except (ValueError, IndexError):
            continue  # header, unit and rule rows
    if not ell:
        raise ValueError(f"no data rows parsed from {path}")
    return np.asarray(ell), np.asarray(vel)


def fetch_lab_longitude(l_deg: float):  # pragma: no cover - network
    """Download the LAB $(b, v)$ slice at integer/half-degree longitude ``l_deg`` from VizieR (cached)."""
    from . import data as _data

    code = int(round(l_deg * 10))
    url = f"https://vizier.cfa.harvard.edu/ftp/cats/VIII/76/bvmaps/L{code:04d}.fits.gz"
    target = _data.data_dir() / f"lab_L{code:04d}.fits.gz"
    if not target.exists():
        _data._download(url, target)
    return target


def _flat_stats(R: np.ndarray, V: np.ndarray, rmin: float) -> dict:
    """Mean, sample scatter (ddof=1), SEM and count of the flat sample $R>$ ``rmin``."""
    flat = V[R > rmin]
    n = int(flat.size)
    sd = float(np.std(flat, ddof=1)) if n > 1 else float("nan")
    return {
        "mean": float(np.mean(flat)),
        "scatter": sd,
        "sem": sd / np.sqrt(n) if n > 1 else float("nan"),
        "n": n,
    }


def _slope_fit(R: np.ndarray, V: np.ndarray, rmin: float) -> tuple[float, float]:
    """Least-squares slope of $V(R)$ over the flat sample, with its standard error (km/s/kpc)."""
    m = R > rmin
    if m.sum() < 3:
        return float("nan"), float("nan")
    A = np.vstack([R[m], np.ones(int(m.sum()))]).T
    coef, *_ = np.linalg.lstsq(A, V[m], rcond=None)
    resid = V[m] - A @ coef
    se = float(np.sqrt(np.sum(resid**2) / (m.sum() - 2) / np.sum((R[m] - R[m].mean()) ** 2)))
    return float(coef[0]), se


def _contiguity(vel_kms: np.ndarray, spectrum: np.ndarray, threshold_k: float) -> int:
    """Number of separate velocity runs above ``threshold_k`` (1 = a single contiguous edge)."""
    above = np.isfinite(spectrum) & (spectrum > threshold_k)
    if not above.any():
        return 0
    return int(np.sum(np.diff(above.astype(int)) == 1)) + int(above[0])


def run(
    out: str = ".",
    *,
    offline: bool = False,
    threshold_k: float = 2.0,
    step_deg: float = 1.0,
) -> dict:
    """Build the inner-Galaxy rotation curve (real LAB longitudes, or synthetic offline).

    Longitudes are sampled from $10\\arcdeg$ to $80\\arcdeg$ every ``step_deg``; at the default
    $1\\arcdeg$ that is 71 sightlines. The earlier 8-longitude sampling was too sparse to resolve
    the slope of $V(R)$, and reported it as consistent with zero when it is not.

    Both terminal-velocity estimators run on the same spectra — the fixed-threshold crossing and
    the error-function edge fit — so the estimator systematic is a *measured* per-run quantity,
    not a citation. On a real run the curve is additionally compared against the independently
    measured VGPS terminal velocities of McClure-Griffiths & Dickey (2016); offline, the same
    comparison runs against the fixture's injected curve. Sensitivity variants (threshold sweep,
    drop-a-point) are committed with the headline.
    """
    from pathlib import Path

    longitudes = np.arange(10.0, 80.0 + 1e-9, step_deg)
    if offline:
        # Real LAB edges span ~2-12 km/s and the threshold overshoot scales with that width,
        # so a single-width fixture cannot exercise the width-bias relation measured below.
        slices = [
            synthetic_lv_slice(ell, seed=i, edge_width_kms=3.0 + 5.0 * (i % 6) / 5.0)
            for i, ell in enumerate(longitudes)
        ]
        source = "synthetic"
        ref_l, ref_v = synthetic_reference_curve(longitudes)
        reference = "synthetic injected curve"
    else:  # pragma: no cover - network
        slices = [read_lab_slice(fetch_lab_longitude(ell)) for ell in longitudes]
        source = "LAB (Kalberla et al. 2005)"
        ref_l, ref_v = fetch_mgd2016()
        reference = "VGPS (McClure-Griffiths & Dickey 2016, Table 1)"
    rmin = 4.0  # the bar dominates non-circular motions at R < ~4 kpc; exclude it
    R, V = rotation_curve(longitudes, slices, threshold_k=threshold_k)
    Re, Ve = rotation_curve(longitudes, slices, threshold_k=threshold_k, estimator="edge")
    thr_stats = _flat_stats(R, V, rmin)
    edge_stats = _flat_stats(Re, Ve, rmin)
    slope, slope_se = _slope_fit(R, V, rmin)
    slope_e, slope_e_se = _slope_fit(Re, Ve, rmin)
    # Keplerian contrast: a point-mass decline normalised at the innermost flat point.
    m = R > rmin
    kepler = float(V[m][0] * np.sqrt(R[m][0] / R[m][-1])) if m.sum() >= 2 else float("nan")

    # Per-longitude terminal velocities and edge widths, in longitude order.
    v_thr, v_edge, widths, runs = [], [], [], []
    for lat, vv, d in slices:
        spec = d[int(np.argmin(np.abs(lat)))]
        v_thr.append(terminal_velocity(vv, spec, threshold_k=threshold_k))
        ve, w = terminal_velocity_edge(vv, spec, threshold_k=threshold_k)
        v_edge.append(ve)
        widths.append(w)
        runs.append(_contiguity(vv, spec, threshold_k))
    v_thr_a, v_edge_a, widths_a = (np.asarray(x) for x in (v_thr, v_edge, widths))

    # Cross-survey validation: both estimators against the same external reference curve.
    cmp_thr = compare_terminal_velocities(longitudes, v_thr_a, ref_l, ref_v)
    cmp_edge = compare_terminal_velocities(longitudes, v_edge_a, ref_l, ref_v)
    # The reference's own V(R) slope, on its native sampling -- the comparand for ours.
    ref_R = R0_KPC * np.sin(np.radians(ref_l))
    ref_V = ref_v + V0_KMS * np.sin(np.radians(ref_l))
    ref_slope, ref_slope_se = _slope_fit(ref_R, ref_V, rmin)
    ref_slope_n = int(np.sum(ref_R > rmin))

    # Does the threshold estimator's disagreement with the reference close as the threshold
    # rises towards the 20 K the reference itself was seeded from? If the gap is the wing
    # overshoot it must, and that is a prediction the data can refuse.
    sweep, sweep_ref = {}, {}
    for thr in (1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 40.0):
        Rs, Vs = rotation_curve(longitudes, slices, threshold_k=thr)
        sweep[f"{thr:g}"] = round(_flat_stats(Rs, Vs, rmin)["mean"], 1)
        vt_s = np.asarray(
            [
                terminal_velocity(vv, d[int(np.argmin(np.abs(lat)))], threshold_k=thr)
                for (lat, vv, d) in slices
            ]
        )
        sweep_ref[f"{thr:g}"] = round(
            compare_terminal_velocities(longitudes, vt_s, ref_l, ref_v)["mean"], 2
        )

    # Does the answer depend on how the reference is averaged onto our longitudes? The note
    # quotes this range, so it has to be in the committed evidence rather than in a notebook.
    window_sweep = {}
    for half in (0.25, 0.5, 1.0):
        c = compare_terminal_velocities(longitudes, v_edge_a, ref_l, ref_v, half_width_deg=half)
        d = np.asarray(c["delta_kms"])
        window_sweep[f"{half:g}"] = {
            "mean": round(float(np.mean(d)), 3) if d.size else float("nan"),
            "median": round(float(np.median(d)), 3) if d.size else float("nan"),
            "n": c["n"],
        }
    _wv = [v[k] for v in window_sweep.values() for k in ("mean", "median") if np.isfinite(v[k])]

    # Is the threshold bias set by the width of the profile edge, as the erfc model says?
    bias = v_thr_a - v_edge_a
    ok = np.isfinite(bias) & np.isfinite(widths_a)
    if int(ok.sum()) > 2:
        width_corr = float(np.corrcoef(widths_a[ok], bias[ok])[0, 1])
        width_slope = float(np.polyfit(widths_a[ok], bias[ok], 1)[0])
    else:
        width_corr = width_slope = float("nan")

    drop_inner = float(np.mean(np.delete(V[m], 0))) if m.sum() > 1 else float("nan")
    drop_outer_two = float(np.mean(V[m][:-2])) if m.sum() > 2 else float("nan")

    metrics = {
        "source": source,
        "reference": reference,
        "step_deg": step_deg,
        "longitudes_deg": longitudes.tolist(),
        "R_kpc": R.tolist(),
        "V_kms": V.tolist(),
        "V_kms_edge": Ve.tolist(),
        "v_term_threshold_kms": v_thr_a.tolist(),
        "v_term_edge_kms": v_edge_a.tolist(),
        "edge_width_kms": widths_a.tolist(),
        "n_above_threshold_runs": runs,
        "V_flat_mean_kms": thr_stats["mean"],
        "V_flat_scatter_kms": thr_stats["scatter"],
        "V_flat_sem_kms": thr_stats["sem"],
        "n_flat": thr_stats["n"],
        "V_flat_edge_mean_kms": edge_stats["mean"],
        "V_flat_edge_scatter_kms": edge_stats["scatter"],
        "V_flat_edge_sem_kms": edge_stats["sem"],
        # the estimator systematic, measured on the same spectra
        "threshold_minus_edge_kms": thr_stats["mean"] - edge_stats["mean"],
        "slope_kms_per_kpc": slope,
        "slope_se_kms_per_kpc": slope_se,
        "slope_sigma": abs(slope / slope_se) if slope_se else float("nan"),
        "slope_edge_kms_per_kpc": slope_e,
        "slope_edge_se_kms_per_kpc": slope_e_se,
        "slope_edge_sigma": abs(slope_e / slope_e_se) if slope_e_se else float("nan"),
        "reference_slope_kms_per_kpc": ref_slope,
        "reference_slope_se_kms_per_kpc": ref_slope_se,
        "reference_slope_n": ref_slope_n,
        "compare_threshold": cmp_thr,
        "compare_edge": cmp_edge,
        "keplerian_at_rmax_kms": kepler,
        "threshold_sweep_vflat": sweep,
        "threshold_sweep_minus_reference": sweep_ref,
        "compare_edge_window_sweep": window_sweep,
        "compare_edge_window_min": min(_wv) if _wv else float("nan"),
        "compare_edge_window_max": max(_wv) if _wv else float("nan"),
        "width_bias_corr": width_corr,
        "width_bias_slope": width_slope,
        "n_noncontiguous": int(np.sum(np.asarray(runs) > 1)),
        "n_edge_fit_failed": int(np.sum(~np.isfinite(v_edge_a))),
        "vflat_drop_innermost": drop_inner,
        "vflat_drop_outermost_two": drop_outer_two,
        "flat_radius_min_kpc": rmin,
        "threshold_k": threshold_k,
        "R0_kpc": R0_KPC,
        "V0_kms": V0_KMS,
    }
    op = Path(out)
    paper = op / "papers" / "hi"
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "rotation_curve.json")
    _figure(R, V, Re, Ve, metrics, paper / "figures")
    _figure_comparison(metrics, paper / "figures")
    _write_macros(metrics, paper / "generated" / "macros.tex")
    if not offline:  # pragma: no cover - the real table is committed evidence
        _write_table(longitudes, R, V, Ve, rmin, paper / "generated" / "curve_table.tex")
    return metrics


def _write_table(longitudes, R, V, Ve, rmin: float, path, *, every_deg: float = 5.0) -> None:
    """Emit the per-longitude (l, R, V_threshold, V_edge) rows the paper's table \\input{}s.

    The curve is now sampled every degree, which is 71 rows — too many for a two-column
    table and no more informative than a regular subsample, so the printed table steps by
    ``every_deg``. The full per-longitude curve is in the committed results JSON, which is
    the evidence; this is the reader's excerpt.
    """
    from pathlib import Path

    order = np.argsort(np.asarray(R))
    lon = np.asarray(sorted(longitudes, key=lambda x: R0_KPC * np.sin(np.radians(x))))
    lines = []
    for i in order:
        if round(float(lon[i])) % int(every_deg):
            continue
        note = "" if R[i] > rmin else r" (bar; excluded)"
        lines.append(rf"{lon[i]:.0f} & {R[i]:.2f} & {V[i]:.1f} & {Ve[i]:.1f}{note} \\")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _write_macros(m: dict, path) -> None:
    """Emit namespaced LaTeX macros so the paper hard-codes no number.

    The synthetic and real runs give different values for every headline quantity, so the macro
    names are namespaced (``\\hiSyn*`` / ``\\hiReal*``) — the un-namespaced ``\\hiVflat`` was one
    of the mode-dependent names the 2026-08 audit flagged (230 offline vs 257 real), protected
    only by the merge guard's provenance rule until now.
    """
    from pathlib import Path

    def num(value: float, fmt: str = ".2f") -> str:
        """Format a metric, or the repo's placeholder if it is not finite.

        A macro must never carry the literal ``nan``: it renders in the PDF and reads as a
        value. ``--`` is the placeholder the arXiv assembler already blocks on, so a
        non-finite metric that the paper actually cites fails loudly at packaging time
        instead of shipping.
        """
        return "--" if not np.isfinite(value) else format(value, fmt)

    real = not str(m.get("source", "")).lower().startswith("synthetic")
    ns, other = ("hiReal", "hiSyn") if real else ("hiSyn", "hiReal")
    excess = 100.0 * (m["V_flat_mean_kms"] - m["V0_kms"]) / m["V0_kms"]
    excess_edge = 100.0 * (m["V_flat_edge_mean_kms"] - m["V0_kms"]) / m["V0_kms"]
    sweep = m.get("threshold_sweep_vflat", {})
    sweep_ref = m.get("threshold_sweep_minus_reference", {})
    ct = m.get("compare_threshold", {})
    ce = m.get("compare_edge", {})
    values = (
        ("Vflat", num(m["V_flat_mean_kms"], ".0f")),
        ("VflatScatter", num(m["V_flat_scatter_kms"], ".0f")),
        ("VflatSem", num(m["V_flat_sem_kms"], ".1f")),
        ("Nflat", f"{m['n_flat']}"),
        ("VflatEdge", num(m["V_flat_edge_mean_kms"], ".0f")),
        ("VflatEdgeScatter", num(m["V_flat_edge_scatter_kms"], ".0f")),
        ("VflatEdgeSem", num(m["V_flat_edge_sem_kms"], ".1f")),
        ("ThrMinusEdge", num(m["threshold_minus_edge_kms"], ".1f")),
        ("ExcessPct", num(excess, ".0f")),
        ("ExcessPctEdge", num(excess_edge, ".1f")),
        ("Slope", num(m["slope_kms_per_kpc"])),
        ("SlopeErr", num(m["slope_se_kms_per_kpc"])),
        ("SlopeSigma", num(m["slope_sigma"], ".1f")),
        ("SlopeEdge", num(m["slope_edge_kms_per_kpc"])),
        ("SlopeEdgeErr", num(m["slope_edge_se_kms_per_kpc"])),
        ("SlopeEdgeSigma", num(m["slope_edge_sigma"], ".1f")),
        ("RefSlope", num(m["reference_slope_kms_per_kpc"])),
        ("RefSlopeErr", num(m["reference_slope_se_kms_per_kpc"])),
        ("RefSlopeN", f"{m['reference_slope_n']}"),
        ("RefN", f"{ce['n']}"),
        ("RefEdgeOffset", num(ce.get("mean", float("nan")))),
        ("RefEdgeOffsetSd", num(ce.get("sd", float("nan")))),
        ("RefEdgeOffsetSem", num(ce.get("sem", float("nan")))),
        ("RefThrOffset", num(ct.get("mean", float("nan")))),
        ("RefThrOffsetSd", num(ct.get("sd", float("nan")))),
        ("RefThrOffsetSem", num(ct.get("sem", float("nan")))),
        ("SweepRefLo", num(sweep_ref.get("1.5", float("nan")), ".1f")),
        ("SweepRefTwenty", num(sweep_ref.get("20", float("nan")), ".1f")),
        ("SweepRefHi", num(sweep_ref.get("40", float("nan")), ".1f")),
        ("RefWindowMin", num(m.get("compare_edge_window_min", float("nan")))),
        ("RefWindowMax", num(m.get("compare_edge_window_max", float("nan")))),
        ("WidthCorr", num(m["width_bias_corr"])),
        ("WidthSlope", num(m["width_bias_slope"])),
        ("NNoncontig", f"{m['n_noncontiguous']}"),
        ("NEdgeFail", f"{m['n_edge_fit_failed']}"),
        ("Kepler", num(m["keplerian_at_rmax_kms"], ".0f")),
        ("Rmin", f"{m['flat_radius_min_kpc']:.0f}"),
        ("Rmax", num(max(m["R_kpc"]), ".1f")),
        ("SweepLoThr", num(sweep.get("1.5", float("nan")), ".0f")),
        ("SweepHiThr", num(sweep.get("5", float("nan")), ".0f")),
        ("DropInner", num(m["vflat_drop_innermost"], ".0f")),
        ("DropOuterTwo", num(m["vflat_drop_outermost_two"], ".0f")),
        ("Reference", str(m.get("reference", "--")).replace("&", r"\&")),
    )
    lines = [
        "% Auto-generated by jansky_research.hi._write_macros — do not edit by hand.",
        "% Mode-dependent values are namespaced (hiSyn*/hiReal*); the inactive namespace holds",
        "% placeholders so offline CI and real runs never collide.",
        rf"\newcommand{{\hiSource}}{{{m['source']}}}",
        rf"\newcommand{{\hiRzero}}{{{m['R0_kpc']:.2f}}}",
        rf"\newcommand{{\hiVzero}}{{{m['V0_kms']:.0f}}}",
        rf"\newcommand{{\hiNlong}}{{{len(m['longitudes_deg'])}}}",
        rf"\newcommand{{\hiStep}}{{{m.get('step_deg', 1.0):g}}}",
    ]
    for suffix, value in values:
        lines.append(rf"\newcommand{{\{ns}{suffix}}}{{{value}}}")
        lines.append(rf"\newcommand{{\{other}{suffix}}}{{--}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and
    # would otherwise blank the other mode's macros with '--'. `make figures`
    # runs every slice offline in the repo root, so without this an offline
    # rebuild silently empties this paper. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _figure(R, V, Re, Ve, m: dict, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rmin = m["flat_radius_min_kpc"]
    bar = R <= rmin
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(R[~bar], V[~bar], "o", color="C0", ms=3, label="threshold (2 K)")
    ax.plot(R[bar], V[bar], "o", mfc="none", color="C0", ms=3, label="bar region (excluded)")
    ax.plot(Re, Ve, "s", color="C2", ms=3, label="edge fit")
    # The external reference curve the edge fit is validated against, on the same axes.
    ref = m.get("compare_edge", {})
    if ref.get("n"):
        rl = np.asarray(ref["longitudes_deg"])
        ax.plot(
            R0_KPC * np.sin(np.radians(rl)),
            np.asarray(ref["reference_v_term_kms"]) + V0_KMS * np.sin(np.radians(rl)),
            "-",
            color="C3",
            lw=1.2,
            alpha=0.8,
            label="reference curve",
        )
    flat = ~bar
    # Keplerian decline normalised at the innermost flat point, for the non-Keplerian contrast
    rk = np.linspace(R[flat][0], R.max(), 60)
    ax.plot(
        rk, V[flat][0] * np.sqrt(R[flat][0] / rk), color="0.4", ls="-.", lw=1, label="Keplerian"
    )
    ax.axhline(V0_KMS, color="0.6", ls=":", label=f"$V_0={V0_KMS:.0f}$ km/s")
    ax.axhline(
        m["V_flat_mean_kms"],
        color="r",
        ls="--",
        xmin=0.35,
        label=f"flat mean {m['V_flat_mean_kms']:.0f} km/s",
    )
    ax.set(
        xlabel="galactocentric radius $R$ (kpc)",
        ylabel="rotation speed $V$ (km/s)",
        title="Inner Milky Way rotation curve (tangent point)",
        ylim=(0, 300),
    )
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "rotation_curve.pdf")
    plt.close(fig)


def _figure_comparison(m: dict, out_dir) -> None:
    """The note's single exhibit: both estimators against the reference, and the threshold sweep.

    Left, the terminal velocities themselves (not $V(R)$, so no $R_0$/$V_0$ choice enters);
    right, the mean offset from the reference as a function of the brightness-temperature
    threshold, which is the evidence that the gap is the wing overshoot and not a survey
    difference: it closes where the reference's own estimator is seeded.
    """
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ct, ce = m.get("compare_threshold", {}), m.get("compare_edge", {})
    if not ce.get("n"):
        return
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.1, 3.0))

    lo = np.asarray(ce["longitudes_deg"])
    a.plot(lo, ce["reference_v_term_kms"], "-", color="C3", lw=2.0, label="VGPS reference")
    a.plot(
        np.asarray(ct["longitudes_deg"]),
        ct["v_term_kms"],
        "o",
        color="C0",
        ms=3.5,
        label="LAB, threshold (2 K)",
    )
    # open squares, so the reference line stays visible underneath where the two agree
    a.plot(lo, ce["v_term_kms"], "s", color="C2", ms=4, mfc="none", mew=1.1, label="LAB, edge fit")
    a.set(xlabel=r"Galactic longitude $\ell$ (deg)", ylabel=r"$v_\mathrm{term}$ (km s$^{-1}$)")
    a.legend(fontsize=7)

    sweep = m.get("threshold_sweep_minus_reference", {})
    thr = sorted(float(k) for k in sweep)
    val = [sweep[f"{t:g}"] for t in thr]
    ok = [(t, v) for t, v in zip(thr, val, strict=True) if np.isfinite(v)]
    b.plot(
        [t for t, _ in ok], [v for _, v in ok], "o-", color="C0", ms=4, label="threshold estimator"
    )
    b.axhline(ce["mean"], color="C2", ls="--", lw=1.4, label="edge fit")
    b.axhline(0.0, color="0.6", lw=0.8)
    b.set(
        xscale="log",
        xlabel=r"brightness-temperature threshold (K)",
        ylabel=r"mean offset from reference (km s$^{-1}$)",
    )
    # the default log locator crowds these labels into unreadable overlap
    b.set_xticks([t for t, _ in ok])
    b.set_xticklabels([f"{t:g}" for t, _ in ok], fontsize=7)
    b.minorticks_off()
    b.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "vgps_comparison.pdf")
    plt.close(fig)


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Build the Milky Way HI rotation curve.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true", help="use the synthetic fixture (no network)")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
