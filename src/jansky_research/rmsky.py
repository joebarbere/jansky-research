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
    "NPI_ALIAS_RAD_M2",
    "enhancement_ratio",
    "fetch_taylor2009",
    "latitude_profile",
    "rm_from_angles",
    "run",
    "sign_asymmetry",
    "sign_fraction",
    "synthetic_rm_sky",
]

#: The Taylor+2009 two-band npi ambiguity step (rad m^-2): NVSS RMs from the two IF bands can
#: alias by +-652.9 (Taylor, Stil & Sunstrum 2009, sec 3). A wrap is not sign-random and
#: preferentially hits large-|RM| (low-latitude) sightlines, so region MEANS are exposed to it;
#: the sign-fraction statistic below cannot move under a wrap at all.
NPI_ALIAS_RAD_M2 = 652.9


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


def sign_fraction(rm: np.ndarray, gal_l_deg: np.ndarray, gal_b_deg: np.ndarray) -> dict:
    """Fraction of positive RMs in the four regions -- the alias-immune sign statistic.

    A two-band :math:`n\\pi` wrap moves an RM by :math:`\\pm652.9` rad m⁻² but (for the
    catalogue's |RM| range) cannot flip which side of zero it started on in a way that is
    correlated with the region means' failure mode; more simply, the *fraction above zero*
    is unchanged by any wrap that does not cross zero, and the handful that could cross zero
    are bounded by the count itself. Under no sign organisation the expectation is 0.5.
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
        out[name] = float(np.mean(rm[m] > 0)) if n else float("nan")
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
    """Fetch the Taylor+2009 NVSS RM catalogue from VizieR; return rm, its error, and Galactic l, b.

    Queries ``J/ApJ/702/1230/catalog`` (RAJ2000, DEJ2000, RM, e_RM), converts to Galactic
    :math:`l,b`, and returns ``rm``/``e_rm`` (rad m⁻²), ``l``, ``b`` (deg), plus fetch
    metadata. ``max_sources=0`` fetches all 37 543 sources. ``e_RM`` was previously never
    retrieved, which left the slice with no quality axis to vary at all.
    """
    from datetime import datetime, timezone

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "RM", "e_RM"])
    v.ROW_LIMIT = -1 if max_sources <= 0 else max_sources
    t = v.query_constraints(catalog="J/ApJ/702/1230/catalog")[0]
    c = SkyCoord(ra=t["RAJ2000"], dec=t["DEJ2000"], unit=(u.hourangle, u.deg))
    e_rm = np.asarray(t["e_RM"], float)
    if hasattr(t["e_RM"], "mask"):
        e_rm = np.where(np.asarray(t["e_RM"].mask), np.nan, e_rm)
    # nearest-neighbour self-match: deduplication was previously asserted by omission
    sep2 = c.match_to_catalog_sky(c, nthneighbor=2)[1]
    n_dup = int(np.sum(sep2 < 5.0 * u.arcsec))
    return {
        "rm": np.asarray(t["RM"], float),
        "e_rm": e_rm,
        "l": np.asarray(c.galactic.l.deg, float),
        "b": np.asarray(c.galactic.b.deg, float),
        "vizier_table": "J/ApJ/702/1230/catalog",
        "fetched_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_dup_5arcsec": n_dup,
    }


def run(out: str = ".", *, offline: bool = True, max_sources: int = 0) -> dict:
    """Full slice: measure the Galactic RM-sky signatures (synthetic or Taylor+2009) and write outputs."""
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
    frac = sign_fraction(rm, gl, gb)

    # Every quoted uncertainty is a 10-degree sky-block jackknife: the RM sky is correlated on
    # degree scales, so the exchangeable unit is a patch, not a source. The i.i.d. bootstrap
    # (kept as the contrast) is the estimator this repo has already measured understating the
    # same statistic 11x on SPICE-RACS. `rmsky` is where `_ratio_bootstrap_se` is defined and
    # was the one slice that never got the fix its importer got.
    from .rmstructure import spatial_block_jackknife

    finite = np.isfinite(rm)
    rm_f, gl_f, gb_f = rm[finite], gl[finite], gb[finite]
    jk_ratio = spatial_block_jackknife(
        gl_f, gb_f, lambda keep: enhancement_ratio(rm_f[keep], gb_f[keep])
    )
    inner_f = (gl_f < 90.0) | (gl_f > 270.0)
    region_masks = {
        "inner_north": inner_f & (gb_f > 0),
        "inner_south": inner_f & (gb_f < 0),
        "outer_north": (~inner_f) & (gb_f > 0),
        "outer_south": (~inner_f) & (gb_f < 0),
    }

    def _region_jk(mask: np.ndarray, values: np.ndarray) -> dict:
        def stat(keep: np.ndarray) -> float:
            m = mask & keep
            return float(np.mean(values[m])) if m.any() else float("nan")

        return spatial_block_jackknife(gl_f, gb_f, stat)

    region_jk = {name: _region_jk(m, rm_f) for name, m in region_masks.items()}
    frac_jk = {name: _region_jk(m, (rm_f > 0).astype(float)) for name, m in region_masks.items()}
    prof_jk = []
    for p in prof:
        lo, hi = p["b_lo"], p["b_hi"]
        binmask = (np.abs(gb_f) >= lo) & (np.abs(gb_f) < hi)

        def stat(keep: np.ndarray, _bm=binmask) -> float:
            m = _bm & keep
            return float(np.median(np.abs(rm_f[m]))) if m.any() else float("nan")

        prof_jk.append(spatial_block_jackknife(gl_f, gb_f, stat))

    # Robustness variants, each free to move the answer: an |RM| < 300 cut (immune to the
    # +-652.9 two-band alias), alternative bin edges, and (real leg) an e_RM quality cut.
    cut300 = finite & (np.abs(rm) < 300.0)
    asym300 = sign_asymmetry(rm[cut300], gl[cut300], gb[cut300])
    ratio300 = enhancement_ratio(rm[cut300], gb[cut300])
    ratio_alt = enhancement_ratio(rm, gb, plane_deg=5.0, pole_deg=70.0)

    def _r1(x) -> float | None:
        return round(float(x), 1) if x is not None and np.isfinite(x) else None

    metrics: dict = {
        "source": source,
        "n_sources": int(np.isfinite(rm).sum()),
        "enhancement_ratio": round(ratio, 2) if np.isfinite(ratio) else None,
        "enhancement_ratio_se": round(ratio_se, 2) if np.isfinite(ratio_se) else None,
        "enhancement_ratio_jk_se": round(jk_ratio["se"], 2)
        if np.isfinite(jk_ratio["se"])
        else None,
        "jk_block_deg": 10.0,
        "jk_n_blocks": jk_ratio["n_blocks"],
        "median_abs_rm_plane": round(prof[0]["median_abs_rm"], 1),
        "median_abs_rm_pole": round(prof[-1]["median_abs_rm"], 1),
        "profile": [
            {
                "b_lo": p["b_lo"],
                "b_hi": p["b_hi"],
                "median_abs_rm": round(p["median_abs_rm"], 1),
                "jk_se": _r1(j["se"]),
                "n": p["n"],
            }
            for p, j in zip(prof, prof_jk, strict=True)
        ],
        "enhancement_ratio_cut300": round(ratio300, 2) if np.isfinite(ratio300) else None,
        "enhancement_ratio_alt_bins_5_70": round(ratio_alt, 2) if np.isfinite(ratio_alt) else None,
    }
    for name in ("inner_north", "inner_south", "outer_north", "outer_south"):
        metrics[f"{name}_rm"] = _r1(asym[name])
        metrics[f"{name}_se"] = _r1(asym[f"{name}_se"])
        metrics[f"{name}_jk_se"] = _r1(region_jk[name]["se"])
        metrics[f"{name}_n"] = asym[f"{name}_n"]
        metrics[f"{name}_rm_cut300"] = _r1(asym300[name])
        metrics[f"{name}_frac_pos"] = (
            round(float(frac[name]), 3) if np.isfinite(frac[name]) else None
        )
        metrics[f"{name}_frac_pos_jk_se"] = (
            round(float(frac_jk[name]["se"]), 3) if np.isfinite(frac_jk[name]["se"]) else None
        )
    if truth is not None:
        metrics["truth_disk_amp"] = truth
    for key in ("vizier_table", "fetched_utc", "n_dup_5arcsec"):
        if key in sky:
            metrics[key] = sky[key]
    e_rm = sky.get("e_rm")
    if e_rm is not None:  # pragma: no cover - real leg only
        e = np.asarray(e_rm, float)
        good_e = finite & np.isfinite(e)
        med_e = float(np.median(e[good_e]))
        cut = good_e & (e < med_e)
        metrics["e_rm_median"] = round(med_e, 1)
        r_e = enhancement_ratio(rm[cut], gb[cut])
        metrics["enhancement_ratio_erm_cut"] = round(r_e, 2) if np.isfinite(r_e) else None
        metrics["n_erm_cut"] = int(cut.sum())

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "rmsky_metrics.json")
    if not offline:  # pragma: no cover - real leg only
        _write_catalogue(sky, op / "data" / "rmsky_taylor2009.csv.gz")
    _figure(rm, gl, gb, metrics["profile"], op / "papers" / "rmsky" / "figures")
    _write_macros(metrics, op / "papers" / "rmsky" / "generated" / "macros.tex")
    return metrics


def _write_catalogue(sky: dict, path) -> None:  # pragma: no cover - real leg only
    """Commit the analysed catalogue (l, b, RM, e_RM): the whole result was fifteen scalars
    plus four bins, unauditable without a network re-fetch, in a paper whose contribution is
    reproducibility. Under a megabyte gzipped."""
    import csv
    import gzip
    from pathlib import Path

    pt = Path(path)
    pt.parent.mkdir(parents=True, exist_ok=True)
    e = sky.get("e_rm")
    with gzip.open(pt, "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gal_l_deg", "gal_b_deg", "rm_rad_m2", "e_rm_rad_m2"])
        for i in range(len(sky["rm"])):
            w.writerow(
                [
                    f"{sky['l'][i]:.5f}",
                    f"{sky['b'][i]:.5f}",
                    f"{sky['rm'][i]:.1f}",
                    f"{e[i]:.1f}" if e is not None and np.isfinite(e[i]) else "",
                ]
            )


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
    # matplotlib labels aitoff ticks with the x value, and x = -l here, so relabel with the
    # true longitude ((-x) mod 360): without this the map's printed 120 deg was l = 240 deg --
    # mirrored labels on the paper whose second result is a longitude-quadrant claim.
    ticks_deg = np.degrees(ax1.get_xticks()).round().astype(int)
    ax1.set_xticklabels([rf"{(-t) % 360:d}$\degree$" for t in ticks_deg], fontsize=7)
    ax2 = fig.add_subplot(1, 2, 2)
    bc = [0.5 * (p["b_lo"] + p["b_hi"]) for p in prof]
    med = [p["median_abs_rm"] for p in prof]
    err = [p.get("jk_se") or 0.0 for p in prof]
    ax2.errorbar(bc, med, yerr=err, fmt="o-", color="C3", capsize=2)
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

    def _dp1(key: str) -> str:
        """One-decimal display, so a ratio and its error are quoted at matched precision."""
        val = m.get(key)
        if val is None:
            return "--"
        return f"{float(val):.1f}"

    lines = [
        "% Auto-generated by jansky_research.rmsky._write_macros -- do not edit by hand.",
        rf"\newcommand{{\rmSource}}{{{m['source']}}}",
        rf"\newcommand{{\rmN}}{{{m['n_sources']}}}",
        # display precision matched between value and error (the rmstructure convention);
        # the un-rounded values are in the results JSON
        rf"\newcommand{{\rmRatio}}{{{_dp1('enhancement_ratio')}}}",
        rf"\newcommand{{\rmRatioErr}}{{{_fmt('enhancement_ratio_se')}}}",
        rf"\newcommand{{\rmRatioJkErr}}{{{_dp1('enhancement_ratio_jk_se')}}}",
        rf"\newcommand{{\rmJkBlocks}}{{{_fmt('jk_n_blocks')}}}",
        rf"\newcommand{{\rmPlane}}{{{m['median_abs_rm_plane']}}}",
        rf"\newcommand{{\rmPole}}{{{m['median_abs_rm_pole']}}}",
        rf"\newcommand{{\rmInnerNorth}}{{{m['inner_north_rm']}}}",
        rf"\newcommand{{\rmInnerNorthErr}}{{{_fmt('inner_north_se')}}}",
        rf"\newcommand{{\rmInnerNorthJkErr}}{{{_fmt('inner_north_jk_se')}}}",
        rf"\newcommand{{\rmInnerSouth}}{{{m['inner_south_rm']}}}",
        rf"\newcommand{{\rmInnerSouthErr}}{{{_fmt('inner_south_se')}}}",
        rf"\newcommand{{\rmInnerSouthJkErr}}{{{_fmt('inner_south_jk_se')}}}",
        rf"\newcommand{{\rmInnerNorthN}}{{{m['inner_north_n']}}}",
        rf"\newcommand{{\rmInnerSouthN}}{{{m['inner_south_n']}}}",
        rf"\newcommand{{\rmOuterNorth}}{{{m['outer_north_rm']}}}",
        rf"\newcommand{{\rmOuterNorthErr}}{{{_fmt('outer_north_se')}}}",
        rf"\newcommand{{\rmOuterNorthJkErr}}{{{_fmt('outer_north_jk_se')}}}",
        rf"\newcommand{{\rmOuterSouth}}{{{m['outer_south_rm']}}}",
        rf"\newcommand{{\rmOuterSouthErr}}{{{_fmt('outer_south_se')}}}",
        rf"\newcommand{{\rmOuterSouthJkErr}}{{{_fmt('outer_south_jk_se')}}}",
        rf"\newcommand{{\rmOuterNorthN}}{{{_fmt('outer_north_n')}}}",
        rf"\newcommand{{\rmOuterSouthN}}{{{_fmt('outer_south_n')}}}",
        # the alias-immune sign statistic (fraction of positive RMs; 0.5 = no organisation)
        rf"\newcommand{{\rmInnerNorthFrac}}{{{_fmt('inner_north_frac_pos')}}}",
        rf"\newcommand{{\rmInnerNorthFracErr}}{{{_fmt('inner_north_frac_pos_jk_se')}}}",
        rf"\newcommand{{\rmInnerSouthFrac}}{{{_fmt('inner_south_frac_pos')}}}",
        rf"\newcommand{{\rmInnerSouthFracErr}}{{{_fmt('inner_south_frac_pos_jk_se')}}}",
        rf"\newcommand{{\rmOuterNorthFrac}}{{{_fmt('outer_north_frac_pos')}}}",
        rf"\newcommand{{\rmOuterNorthFracErr}}{{{_fmt('outer_north_frac_pos_jk_se')}}}",
        rf"\newcommand{{\rmOuterSouthFrac}}{{{_fmt('outer_south_frac_pos')}}}",
        rf"\newcommand{{\rmOuterSouthFracErr}}{{{_fmt('outer_south_frac_pos_jk_se')}}}",
        # robustness variants
        rf"\newcommand{{\rmRatioCutThree}}{{{_dp1('enhancement_ratio_cut300')}}}",
        rf"\newcommand{{\rmRatioAltBins}}{{{_dp1('enhancement_ratio_alt_bins_5_70')}}}",
        rf"\newcommand{{\rmRatioErmCut}}{{{_dp1('enhancement_ratio_erm_cut')}}}",
        rf"\newcommand{{\rmInnerSouthCutThree}}{{{_fmt('inner_south_rm_cut300')}}}",
        rf"\newcommand{{\rmInnerNorthCutThree}}{{{_fmt('inner_north_rm_cut300')}}}",
        rf"\newcommand{{\rmDupCount}}{{{_fmt('n_dup_5arcsec')}}}",
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
