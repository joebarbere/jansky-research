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

    A fixed brightness-temperature threshold is the simple, standard estimator; relative to
    inflection-point spectral fitting it tends to **overestimate** the terminal velocity (McClure-
    Griffiths & Dickey 2016 find threshold crossings ~7 km/s higher), biasing $V(R)$ high — a known
    systematic the write-up notes.
    """
    vel_kms = np.asarray(vel_kms, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    above = spectrum > threshold_k
    if not above.any():
        return float("nan")
    return float(np.max(vel_kms[above]))


def tangent_point(
    l_deg: float, v_term: float, *, R0: float = R0_KPC, V0: float = V0_KMS
) -> tuple[float, float]:
    """Tangent-point $(R, V)$: $R = R_0\\sin\\ell$ (kpc), $V = v_\\mathrm{term} + V_0\\sin\\ell$ (km/s)."""
    s = np.sin(np.radians(l_deg))
    return float(R0 * s), float(v_term + V0 * s)


def rotation_curve(
    longitudes: np.ndarray, slices, *, threshold_k: float = 2.0
) -> tuple[np.ndarray, np.ndarray]:
    """Build the rotation curve from per-longitude $(b, v)$ slices (uses each $b=0$ spectrum).

    ``slices`` is an iterable of ``(lat_deg, vel_kms, data)`` aligned with ``longitudes``. Returns
    ``(R_kpc, V_kms)`` sorted by radius.
    """
    rad, vel = [], []
    for ell, (lat, v, d) in zip(longitudes, slices, strict=True):
        spec = d[int(np.argmin(np.abs(lat)))]  # b = 0
        r, vv = tangent_point(ell, terminal_velocity(v, spec, threshold_k=threshold_k))
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
    seed: int | None = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Synthetic LAB-like $(b, v)$ slice with a known *flat* rotation curve (offline fixture).

    Injects HI emission out to the terminal velocity implied by a flat curve $V=$ ``v_flat`` —
    $v_\\mathrm{term} = v_\\mathrm{flat} - V_0\\sin\\ell$ — with a sharp high-velocity edge plus
    noise, so :func:`terminal_velocity` and :func:`tangent_point` recover ``v_flat``.
    """
    rng = np.random.default_rng(seed)
    s = np.sin(np.radians(l_deg))
    v_term = v_flat - V0 * s
    vel = np.linspace(-50.0, 320.0, 400)
    lat = np.linspace(-5.0, 5.0, n_lat)
    # bright HI extending to low/negative velocity with a steep falling edge at the terminal
    # velocity (which is ~0 at high longitudes, where R_tan -> R0 and V(R_tan) -> V0).
    profile = 30.0 / (1.0 + np.exp((vel - v_term) / 0.6))
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


def run(out: str = ".", *, offline: bool = False, threshold_k: float = 2.0) -> dict:
    """Build the inner-Galaxy rotation curve (real LAB longitudes, or synthetic offline). Writes a figure."""
    import json
    from pathlib import Path

    longitudes = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0])
    if offline:
        slices = [synthetic_lv_slice(ell, seed=i) for i, ell in enumerate(longitudes)]
        source = "synthetic"
    else:  # pragma: no cover - network
        slices = [read_lab_slice(fetch_lab_longitude(ell)) for ell in longitudes]
        source = "LAB (Kalberla et al. 2005)"
    R, V = rotation_curve(longitudes, slices, threshold_k=threshold_k)
    flat = V[R > 4.0]  # the bar dominates non-circular motions at R < ~4 kpc; exclude it
    metrics = {
        "source": source,
        "longitudes_deg": longitudes.tolist(),
        "R_kpc": R.tolist(),
        "V_kms": V.tolist(),
        "V_flat_mean_kms": float(np.mean(flat)),
        "V_flat_std_kms": float(np.std(flat)),
        "flat_radius_min_kpc": 4.0,
        "R0_kpc": R0_KPC,
        "V0_kms": V0_KMS,
    }
    op = Path(out)
    paper = op / "papers" / "hi"
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "rotation_curve.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(R, V, metrics["V_flat_mean_kms"], paper / "figures")
    _write_macros(metrics, paper / "generated" / "macros.tex")
    return metrics


def _write_macros(m: dict, path) -> None:
    """Emit LaTeX ``\\newcommand`` macros so the paper hard-codes no number (offline + real share this)."""
    from pathlib import Path

    excess = 100.0 * (m["V_flat_mean_kms"] - m["V0_kms"]) / m["V0_kms"]
    lines = [
        "% Auto-generated by jansky_research.hi._write_macros — do not edit by hand.",
        rf"\newcommand{{\hiSource}}{{{m['source']}}}",
        rf"\newcommand{{\hiVflat}}{{{m['V_flat_mean_kms']:.0f}}}",
        rf"\newcommand{{\hiVflatErr}}{{{m['V_flat_std_kms']:.0f}}}",
        rf"\newcommand{{\hiRmin}}{{{m['flat_radius_min_kpc']:.0f}}}",
        rf"\newcommand{{\hiRmax}}{{{max(m['R_kpc']):.1f}}}",
        rf"\newcommand{{\hiRzero}}{{{m['R0_kpc']:.2f}}}",
        rf"\newcommand{{\hiVzero}}{{{m['V0_kms']:.0f}}}",
        rf"\newcommand{{\hiNlong}}{{{len(m['longitudes_deg'])}}}",
        rf"\newcommand{{\hiExcessPct}}{{{excess:.0f}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _figure(R, V, v_flat, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(R, V, "o-", lw=1)
    ax.axhline(V0_KMS, color="0.6", ls=":", label=f"$V_0={V0_KMS:.0f}$ km/s")
    ax.axhline(v_flat, color="r", ls="--", label=f"mean $V={v_flat:.0f}$ km/s")
    ax.set(
        xlabel="galactocentric radius $R$ (kpc)",
        ylabel="rotation speed $V$ (km/s)",
        title="Inner Milky Way rotation curve (tangent point)",
        ylim=(0, 300),
    )
    ax.legend()
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
