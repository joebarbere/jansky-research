"""The Galactic Faraday rotation sky from extragalactic rotation measures (Taylor+2009).

The polarisation angle of a radio source rotates as :math:`\\chi(\\lambda)=\\chi_0+\\mathrm{RM}\\lambda^2`
on its way through the magnetised interstellar medium, with
:math:`\\mathrm{RM}=0.81\\int n_e B_\\parallel\\,\\mathrm{d}l` (rad m⁻²). The RMs of tens of thousands of
*extragalactic* sources therefore map the **Galactic Faraday rotation sky** (tracing the line-of-sight
integral :math:`\\int n_e B_\\parallel\\,\\mathrm{d}l`; isolating :math:`B_\\parallel` would need an
electron-density model, not applied here). Two large-scale signatures are textbook (Taylor, Stil &
Sunstrum 2009): :math:`|\\mathrm{RM}|` is enhanced toward the Galactic **plane** (the disk path length
grows as :math:`\\csc|b|`), and the RM sky is **sign-organised** with a quadrupole-like antisymmetry
(the disk/halo field).

This module reproduces both from the Taylor+2009 NVSS RM catalogue (VizieR ``J/ApJ/702/1230``, public,
no auth), reusing ``jansky.polarization`` for the underlying :math:`\\lambda^2` measurement. Pure
NumPy with a synthetic offline fixture; the real fetch is network-gated. The catalogue's known limits
(the two-band :math:`n\\pi` ambiguity, intrinsic source RM) are reported, not hidden.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "enhancement_ratio",
    "fetch_taylor2009",
    "latitude_profile",
    "rm_from_angles",
    "run",
    "sign_asymmetry",
    "synthetic_rm_sky",
]


def rm_from_angles(wavelengths_m: np.ndarray, angles_rad: np.ndarray) -> float:
    """Recover an RM (rad m⁻²) from polarisation angle versus wavelength.

    Thin wrapper over ``jansky.polarization.rotation_measure_fit`` (fits
    :math:`\\chi=\\chi_0+\\mathrm{RM}\\lambda^2`) — the foundational step by which every catalogue RM is
    measured. Exposed (and tested) here so the slice is anchored to the polarisation helpers.
    """
    from jansky import polarization

    rm, _angle0 = polarization.rotation_measure_fit(
        np.asarray(wavelengths_m, float), np.asarray(angles_rad, float)
    )
    return float(rm)


def latitude_profile(
    rm: np.ndarray,
    gal_b_deg: np.ndarray,
    *,
    edges: tuple[float, ...] = (0.0, 10.0, 30.0, 60.0, 90.0),
) -> list[dict]:
    """Median :math:`|\\mathrm{RM}|` in bins of Galactic latitude :math:`|b|` (the disk enhancement).

    Returns one dict per :math:`|b|` bin with the bin range, the median :math:`|\\mathrm{RM}|`
    (rad m⁻²), and the source count. Sightlines near the plane traverse more magneto-ionic disk, so
    the median rises toward :math:`|b|=0`.
    """
    rm = np.asarray(rm, float)
    absb = np.abs(np.asarray(gal_b_deg, float))
    out: list[dict] = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        m = (absb >= lo) & (absb < hi) & np.isfinite(rm)
        out.append(
            {
                "b_lo": float(lo),
                "b_hi": float(hi),
                "median_abs_rm": float(np.median(np.abs(rm[m]))) if m.any() else float("nan"),
                "n": int(m.sum()),
            }
        )
    return out


def enhancement_ratio(
    rm: np.ndarray, gal_b_deg: np.ndarray, *, plane_deg: float = 10.0, pole_deg: float = 60.0
) -> float:
    """Median :math:`|\\mathrm{RM}|` for :math:`|b|<\\mathrm{plane}` over that for :math:`|b|>\\mathrm{pole}`.

    A single number for the plane enhancement; ~5 for the Galactic disk in Taylor+2009.
    """
    rm = np.asarray(rm, float)
    absb = np.abs(np.asarray(gal_b_deg, float))
    near = (absb < plane_deg) & np.isfinite(rm)
    far = (absb > pole_deg) & np.isfinite(rm)
    if not near.any() or not far.any():
        return float("nan")
    hi = np.median(np.abs(rm[near]))
    lo = np.median(np.abs(rm[far]))
    return float(hi / lo) if lo > 0 else float("nan")


def sign_asymmetry(rm: np.ndarray, gal_l_deg: np.ndarray, gal_b_deg: np.ndarray) -> dict:
    """Mean RM in the four (north/south × inner/outer-Galaxy) regions — a coarse antisymmetry probe.

    Inner Galaxy is :math:`l<90°` or :math:`l>270°`. The true large-scale structure is a *quadrupole*:
    the :math:`l<90°` and :math:`l>270°` halves carry opposite sign at a given :math:`b`, so this mask
    **conflates** them and the recovered means are partial cancellations — a coarse net-sign indicator,
    not a harmonic decomposition. Returns each region's mean RM (rad m⁻²), its standard error, and count.
    """
    rm = np.asarray(rm, float)
    gl = np.asarray(gal_l_deg, float)
    gb = np.asarray(gal_b_deg, float)
    inner = (gl < 90.0) | (gl > 270.0)
    out: dict = {}
    for name, mask in [
        ("inner_north", inner & (gb > 0)),
        ("inner_south", inner & (gb < 0)),
        ("outer_north", (~inner) & (gb > 0)),
        ("outer_south", (~inner) & (gb < 0)),
    ]:
        m = mask & np.isfinite(rm)
        n = int(m.sum())
        out[name] = float(np.mean(rm[m])) if n else float("nan")
        out[f"{name}_se"] = float(np.std(rm[m]) / np.sqrt(n)) if n > 1 else float("nan")
        out[f"{name}_n"] = n
    return out


def _ratio_bootstrap_se(
    rm: np.ndarray,
    gal_b_deg: np.ndarray,
    *,
    plane_deg: float = 10.0,
    pole_deg: float = 60.0,
    n_boot: int = 500,
    seed: int = 0,
) -> float:
    """Bootstrap standard error on the plane/pole :func:`enhancement_ratio` (median-based)."""
    rm = np.asarray(rm, float)
    absb = np.abs(np.asarray(gal_b_deg, float))
    near = np.abs(rm[(absb < plane_deg) & np.isfinite(rm)])
    far = np.abs(rm[(absb > pole_deg) & np.isfinite(rm)])
    if near.size < 2 or far.size < 2:
        return float("nan")
    rng = np.random.default_rng(seed)
    ratios = []
    for _ in range(n_boot):
        lo = np.median(rng.choice(far, far.size))
        if lo > 0:
            ratios.append(np.median(rng.choice(near, near.size)) / lo)
    return float(np.std(ratios)) if ratios else float("nan")


def synthetic_rm_sky(
    n_sources: int = 4000,
    *,
    disk_amp: float = 60.0,
    outer_frac: float = 0.35,
    extragal_sigma: float = 12.0,
    b_floor_deg: float = 3.0,
    seed: int = 0,
) -> dict:
    """Synthetic RM sky: a :math:`\\csc|b|` disk, a sign-organised field, and extragalactic scatter.

    Sources are uniform on the sphere. The Galactic contribution is
    :math:`A\\,\\mathrm{sign}(b)\\,/\\sin(\\max(|b|,b_\\mathrm{floor}))`, with the amplitude larger in
    the inner Galaxy (full ``disk_amp``) than the outer (``outer_frac``×), reproducing the plane
    enhancement *and* the north-positive/south-negative antisymmetry; extragalactic intrinsic RM is
    added as Gaussian scatter. Returns ``rm`` (rad m⁻²), Galactic ``l``/``b`` (deg), and the truth.
    """
    rng = np.random.default_rng(seed)
    gl = rng.uniform(0.0, 360.0, n_sources)
    gb = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, n_sources)))  # uniform on the sphere
    inner = (gl < 90.0) | (gl > 270.0)
    amp = np.where(inner, disk_amp, outer_frac * disk_amp)
    sinb = np.sin(np.radians(np.maximum(np.abs(gb), b_floor_deg)))
    rm_gal = amp * np.sign(gb) / sinb
    rm = rm_gal + rng.normal(0.0, extragal_sigma, n_sources)
    return {"rm": rm, "l": gl, "b": gb, "truth_disk_amp": disk_amp}


def fetch_taylor2009(max_sources: int = 0) -> dict:  # pragma: no cover - network
    """Fetch the Taylor+2009 NVSS RM catalogue from VizieR; return rm and Galactic coordinates.

    Queries ``J/ApJ/702/1230/catalog`` (RAJ2000, DEJ2000, RM), converts to Galactic :math:`l,b`, and
    returns ``rm`` (rad m⁻²), ``l``, ``b`` (deg). ``max_sources=0`` fetches all 37 543 sources.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "RM"])
    v.ROW_LIMIT = -1 if max_sources <= 0 else max_sources
    t = v.query_constraints(catalog="J/ApJ/702/1230/catalog")[0]
    c = SkyCoord(ra=t["RAJ2000"], dec=t["DEJ2000"], unit=(u.hourangle, u.deg))
    return {
        "rm": np.asarray(t["RM"], float),
        "l": np.asarray(c.galactic.l.deg, float),
        "b": np.asarray(c.galactic.b.deg, float),
    }


def run(out: str = ".", *, offline: bool = True, max_sources: int = 0) -> dict:
    """Full slice: measure the Galactic RM-sky signatures (synthetic or Taylor+2009) and write outputs."""
    import json
    from pathlib import Path

    if offline:
        sky = synthetic_rm_sky()
        source = "synthetic"
        truth: float | None = sky["truth_disk_amp"]
    else:  # pragma: no cover - network
        sky = fetch_taylor2009(max_sources=max_sources)
        source = "Taylor+2009 NVSS RM catalogue"
        truth = None

    rm, gl, gb = sky["rm"], sky["l"], sky["b"]
    prof = latitude_profile(rm, gb)
    ratio = enhancement_ratio(rm, gb)
    ratio_se = _ratio_bootstrap_se(rm, gb)
    asym = sign_asymmetry(rm, gl, gb)

    metrics: dict = {
        "source": source,
        "n_sources": int(np.isfinite(rm).sum()),
        "enhancement_ratio": round(ratio, 2) if np.isfinite(ratio) else None,
        "enhancement_ratio_se": round(ratio_se, 2) if np.isfinite(ratio_se) else None,
        "median_abs_rm_plane": round(prof[0]["median_abs_rm"], 1),
        "median_abs_rm_pole": round(prof[-1]["median_abs_rm"], 1),
        "inner_north_rm": round(asym["inner_north"], 1),
        "inner_north_se": round(asym["inner_north_se"], 1),
        "inner_south_rm": round(asym["inner_south"], 1),
        "inner_south_se": round(asym["inner_south_se"], 1),
        "inner_north_n": asym["inner_north_n"],
        "inner_south_n": asym["inner_south_n"],
        "outer_north_rm": round(asym["outer_north"], 1),
        "outer_south_rm": round(asym["outer_south"], 1),
        "profile": [
            {
                "b_lo": p["b_lo"],
                "b_hi": p["b_hi"],
                "median_abs_rm": round(p["median_abs_rm"], 1),
                "n": p["n"],
            }
            for p in prof
        ],
    }
    if truth is not None:
        metrics["truth_disk_amp"] = truth

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "rmsky_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(rm, gl, gb, prof, op / "papers" / "rmsky" / "figures")
    _write_macros(metrics, op / "papers" / "rmsky" / "generated" / "macros.tex")
    return metrics


def _figure(rm, gl, gb, prof, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    # wrap longitude to [-180,180] radians for the Aitoff projection (l increasing to the left)
    lrad = np.radians(-(np.where(gl > 180.0, gl - 360.0, gl)))
    brad = np.radians(gb)
    clip = np.clip(rm, -150.0, 150.0)
    fig = plt.figure(figsize=(10, 3.8))
    ax1 = fig.add_subplot(1, 2, 1, projection="aitoff")
    sc = ax1.scatter(lrad, brad, c=clip, cmap="coolwarm", s=2, vmin=-150, vmax=150)
    fig.colorbar(sc, ax=ax1, shrink=0.6, label="RM (rad m$^{-2}$)")
    ax1.set_title("Galactic RM sky", pad=12)
    ax1.grid(True)
    ax2 = fig.add_subplot(1, 2, 2)
    bc = [0.5 * (p["b_lo"] + p["b_hi"]) for p in prof]
    med = [p["median_abs_rm"] for p in prof]
    ax2.plot(bc, med, "o-", color="C3")
    ax2.set(
        xlabel=r"$|b|$ (deg)",
        ylabel=r"median $|\mathrm{RM}|$ (rad m$^{-2}$)",
        title="Plane enhancement",
    )
    fig.tight_layout()
    fig.savefig(out / "rmsky.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.rmsky._write_macros -- do not edit by hand.",
        rf"\newcommand{{\rmSource}}{{{m['source']}}}",
        rf"\newcommand{{\rmN}}{{{m['n_sources']}}}",
        rf"\newcommand{{\rmRatio}}{{{_fmt('enhancement_ratio')}}}",
        rf"\newcommand{{\rmRatioErr}}{{{_fmt('enhancement_ratio_se')}}}",
        rf"\newcommand{{\rmPlane}}{{{m['median_abs_rm_plane']}}}",
        rf"\newcommand{{\rmPole}}{{{m['median_abs_rm_pole']}}}",
        rf"\newcommand{{\rmInnerNorth}}{{{m['inner_north_rm']}}}",
        rf"\newcommand{{\rmInnerNorthErr}}{{{_fmt('inner_north_se')}}}",
        rf"\newcommand{{\rmInnerSouth}}{{{m['inner_south_rm']}}}",
        rf"\newcommand{{\rmInnerSouthErr}}{{{_fmt('inner_south_se')}}}",
        rf"\newcommand{{\rmInnerNorthN}}{{{m['inner_north_n']}}}",
        rf"\newcommand{{\rmInnerSouthN}}{{{m['inner_south_n']}}}",
        rf"\newcommand{{\rmOuterNorth}}{{{m['outer_north_rm']}}}",
        rf"\newcommand{{\rmOuterSouth}}{{{m['outer_south_rm']}}}",
        rf"\newcommand{{\rmTruth}}{{{_fmt('truth_disk_amp')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="The Galactic Faraday rotation sky (Taylor+2009).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--max-sources", type=int, default=0, help="0 = all 37,543 sources")
    args = p.parse_args(argv)
    metrics = run(args.out, offline=args.offline, max_sources=args.max_sources)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
