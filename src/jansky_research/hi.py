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
    "fetch_lab_longitude",
    "read_lab_slice",
    "rotation_curve",
    "run",
    "synthetic_lv_slice",
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


def run(out: str = ".", *, offline: bool = False, threshold_k: float = 2.0) -> dict:
    """Build the inner-Galaxy rotation curve (real LAB longitudes, or synthetic offline). Writes a figure.

    Both terminal-velocity estimators run on the same spectra — the fixed-threshold crossing and
    the error-function edge fit — so the estimator systematic is a *measured* per-run quantity,
    not a citation. Sensitivity variants (threshold sweep, drop-a-point) are committed with the
    headline.
    """
    from pathlib import Path

    longitudes = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0])
    if offline:
        slices = [synthetic_lv_slice(ell, seed=i) for i, ell in enumerate(longitudes)]
        source = "synthetic"
    else:  # pragma: no cover - network
        slices = [read_lab_slice(fetch_lab_longitude(ell)) for ell in longitudes]
        source = "LAB (Kalberla et al. 2005)"
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
    # Sensitivity variants, committed with the headline so the flat level's robustness to the
    # hand-chosen analysis parameters is evidence, not assertion.
    sweep = {}
    for thr in (1.5, 2.0, 3.0, 5.0):
        Rs, Vs = rotation_curve(longitudes, slices, threshold_k=thr)
        sweep[f"{thr:g}"] = round(_flat_stats(Rs, Vs, rmin)["mean"], 1)
    drop_inner = float(np.mean(np.delete(V[m], 0))) if m.sum() > 1 else float("nan")
    drop_outer_two = float(np.mean(V[m][:-2])) if m.sum() > 2 else float("nan")

    metrics = {
        "source": source,
        "longitudes_deg": longitudes.tolist(),
        "R_kpc": R.tolist(),
        "V_kms": V.tolist(),
        "V_kms_edge": Ve.tolist(),
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
        "slope_edge_kms_per_kpc": slope_e,
        "slope_edge_se_kms_per_kpc": slope_e_se,
        "keplerian_at_rmax_kms": kepler,
        "threshold_sweep_vflat": sweep,
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
    _write_macros(metrics, paper / "generated" / "macros.tex")
    if not offline:  # pragma: no cover - the real table is committed evidence
        _write_table(longitudes, R, V, Ve, rmin, paper / "generated" / "curve_table.tex")
    return metrics


def _write_table(longitudes, R, V, Ve, rmin: float, path) -> None:
    """Emit the per-longitude (l, R, V_threshold, V_edge) rows the paper's table \\input{}s."""
    from pathlib import Path

    order = np.argsort(np.asarray(R))
    lon = np.asarray(sorted(longitudes, key=lambda x: R0_KPC * np.sin(np.radians(x))))
    lines = []
    for i in order:
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

    real = not str(m.get("source", "")).lower().startswith("synthetic")
    ns, other = ("hiReal", "hiSyn") if real else ("hiSyn", "hiReal")
    excess = 100.0 * (m["V_flat_mean_kms"] - m["V0_kms"]) / m["V0_kms"]
    excess_edge = 100.0 * (m["V_flat_edge_mean_kms"] - m["V0_kms"]) / m["V0_kms"]
    sweep = m.get("threshold_sweep_vflat", {})
    values = (
        ("Vflat", f"{m['V_flat_mean_kms']:.0f}"),
        ("VflatScatter", f"{m['V_flat_scatter_kms']:.0f}"),
        ("VflatSem", f"{m['V_flat_sem_kms']:.1f}"),
        ("Nflat", f"{m['n_flat']}"),
        ("VflatEdge", f"{m['V_flat_edge_mean_kms']:.0f}"),
        ("VflatEdgeScatter", f"{m['V_flat_edge_scatter_kms']:.0f}"),
        ("VflatEdgeSem", f"{m['V_flat_edge_sem_kms']:.1f}"),
        ("ThrMinusEdge", f"{m['threshold_minus_edge_kms']:.1f}"),
        ("ExcessPct", f"{excess:.0f}"),
        ("ExcessPctEdge", f"{excess_edge:.1f}"),
        ("Slope", f"{m['slope_kms_per_kpc']:.1f}"),
        ("SlopeErr", f"{m['slope_se_kms_per_kpc']:.1f}"),
        ("Kepler", f"{m['keplerian_at_rmax_kms']:.0f}"),
        ("Rmin", f"{m['flat_radius_min_kpc']:.0f}"),
        ("Rmax", f"{max(m['R_kpc']):.1f}"),
        ("SweepLoThr", f"{sweep.get('1.5', float('nan')):.0f}"),
        ("SweepHiThr", f"{sweep.get('5', float('nan')):.0f}"),
        ("DropInner", f"{m['vflat_drop_innermost']:.0f}"),
        ("DropOuterTwo", f"{m['vflat_drop_outermost_two']:.0f}"),
    )
    lines = [
        "% Auto-generated by jansky_research.hi._write_macros — do not edit by hand.",
        "% Mode-dependent values are namespaced (hiSyn*/hiReal*); the inactive namespace holds",
        "% placeholders so offline CI and real runs never collide.",
        rf"\newcommand{{\hiSource}}{{{m['source']}}}",
        rf"\newcommand{{\hiRzero}}{{{m['R0_kpc']:.2f}}}",
        rf"\newcommand{{\hiVzero}}{{{m['V0_kms']:.0f}}}",
        rf"\newcommand{{\hiNlong}}{{{len(m['longitudes_deg'])}}}",
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
    ax.plot(R[~bar], V[~bar], "o", color="C0", label="threshold (2 K)")
    ax.plot(R[bar], V[bar], "o", mfc="none", color="C0", label="bar region (excluded)")
    ax.plot(Re, Ve, "s", color="C2", ms=4, label="edge fit")
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
