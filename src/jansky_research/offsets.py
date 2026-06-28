"""Radio--optical position offsets of AGN: ICRF3 (VLBI) vs Gaia DR3 (optical).

The VLBI radio position of an AGN marks the synchrotron self-absorbed jet base; the Gaia optical
position can be pulled toward optical jet/host structure. The two therefore do not coincide, and the
*normalised* offset $X=\\sqrt{(\\Delta\\alpha^*/\\sigma_\\alpha)^2+(\\Delta\\delta/\\sigma_\\delta)^2}$
shows a heavy tail far beyond the Rayleigh expectation for pure Gaussian astrometric noise --- a
well-established result (Mignard et al. 2016; Petrov & Kovalev 2017; Kovalev et al. 2017; Lindegren
et al. 2018). This module **reproduces** that excess tail with a small tested tool and builds a
reproducible offset catalogue.

Catalogue-only: ICRF3 (Charlot et al. 2020) cross-matched to Gaia DR3. Composes
:mod:`jansky_research.spectra` (``crossmatch``) and the CDS XMatch service. Pure NumPy + a synthetic
offline fixture for tests.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "fetch_icrf3_gaia",
    "normalised_offset",
    "offset_statistics",
    "radio_optical_offset",
    "run",
    "synthetic_field",
]

DEG_TO_MAS = 3.6e6  # degrees -> milliarcseconds


def radio_optical_offset(
    ra_r: np.ndarray, dec_r: np.ndarray, ra_o: np.ndarray, dec_o: np.ndarray
) -> dict[str, np.ndarray]:
    r"""Radio$\to$optical offset components, total separation (mas), and position angle (deg E of N).

    $\Delta\alpha^*=(\alpha_o-\alpha_r)\cos\delta_r$ and $\Delta\delta=\delta_o-\delta_r$ (both in mas);
    the separation is $\sqrt{\Delta\alpha^{*2}+\Delta\delta^2}$ and the position angle
    $\mathrm{PA}=\mathrm{atan2}(\Delta\alpha^*,\Delta\delta)$ runs from the radio position toward the
    optical, measured East of North.
    """
    ra_r = np.asarray(ra_r, float)
    dec_r = np.asarray(dec_r, float)
    cosd = np.cos(np.radians(dec_r))
    dra = (np.asarray(ra_o, float) - ra_r) * cosd * DEG_TO_MAS
    ddec = (np.asarray(dec_o, float) - dec_r) * DEG_TO_MAS
    offset = np.hypot(dra, ddec)
    pa = np.degrees(np.arctan2(dra, ddec)) % 360.0
    return {"dra_mas": dra, "ddec_mas": ddec, "offset_mas": offset, "pa_deg": pa}


def normalised_offset(
    dra_mas: np.ndarray, ddec_mas: np.ndarray, sig_a_mas: np.ndarray, sig_d_mas: np.ndarray
) -> np.ndarray:
    r"""Significance of an offset: $X=\sqrt{(\Delta\alpha^*/\sigma_\alpha)^2+(\Delta\delta/\sigma_\delta)^2}$.

    The per-axis errors combine the radio and optical formal uncertainties in quadrature
    ($\sigma^2=\sigma_\mathrm{radio}^2+\sigma_\mathrm{Gaia}^2$, supplied by the caller). For pure
    Gaussian astrometric noise $X$ follows a 2-D Rayleigh, so $P(X>x)=e^{-x^2/2}$; a heavier tail is
    real structure.
    """
    a = np.asarray(dra_mas, float) / np.asarray(sig_a_mas, float)
    d = np.asarray(ddec_mas, float) / np.asarray(sig_d_mas, float)
    return np.sqrt(a**2 + d**2)


def offset_statistics(
    x_norm: np.ndarray, offset_mas: np.ndarray | None = None, *, x_cut: float = 3.0
) -> dict:
    r"""Summarise the offset population: the $X>$ ``x_cut`` fraction vs its Rayleigh expectation.

    Returns the count, the median raw offset (mas, if given), the fraction with $X>$ ``x_cut``, the
    Rayleigh expectation $e^{-x_\mathrm{cut}^2/2}$, and their ratio (the *excess* --- the reproduced
    result: AGN show many more significant offsets than Gaussian errors allow).
    """
    x = np.asarray(x_norm, float)
    x = x[np.isfinite(x)]
    n = int(x.size)
    frac = float((x > x_cut).mean()) if n else 0.0
    rayleigh = float(np.exp(-(x_cut**2) / 2.0))
    med = (
        float(np.nanmedian(np.asarray(offset_mas, float)))
        if offset_mas is not None and np.asarray(offset_mas).size
        else float("nan")
    )
    return {
        "n": n,
        "median_offset_mas": med,
        "x_cut": x_cut,
        "frac_x_gt_cut": frac,
        "rayleigh_expectation": rayleigh,
        "excess_ratio": frac / rayleigh if rayleigh > 0 else float("nan"),
    }


def synthetic_field(
    n_sources: int = 4000,
    *,
    structured_fraction: float = 0.15,
    sigma_mas: float = 0.4,
    struct_scale_mas: float = 3.0,
    seed: int = 0,
) -> tuple[dict, dict, np.ndarray]:
    """Synthetic ICRF3-like radio + Gaia-like optical positions with an injected structural-offset tail.

    Most sources have a pure-Gaussian radio$\\to$optical offset (noise at ``sigma_mas`` per axis); a
    ``structured_fraction`` minority carry an additional real offset (exponential magnitude, random
    direction) standing in for optical jet/host structure. Returns ``(radio, optical, is_structured)``
    with ra/dec (deg) and per-axis errors ``e_a``/``e_d`` (mas, on $\\alpha\\cos\\delta$ and $\\delta$).
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(0.0, 360.0, n_sources)
    dec = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, n_sources)))  # uniform on the sphere
    cosd = np.cos(np.radians(dec))
    e_radio = rng.uniform(0.05, 0.5, n_sources)
    e_gaia = rng.uniform(0.05, 0.8, n_sources)
    sig_a = np.hypot(e_radio, e_gaia)
    sig_d = np.hypot(e_radio, e_gaia)
    # noise offset (mas) per axis
    dra = rng.normal(0.0, sig_a)
    ddec = rng.normal(0.0, sig_d)
    is_struct = rng.random(n_sources) < structured_fraction
    mag = rng.exponential(struct_scale_mas, n_sources)
    ang = rng.uniform(0.0, 2 * np.pi, n_sources)
    dra = np.where(is_struct, dra + mag * np.sin(ang), dra)
    ddec = np.where(is_struct, ddec + mag * np.cos(ang), ddec)
    radio = {"ra": ra, "dec": dec, "e_a": e_radio, "e_d": e_radio}
    optical = {
        "ra": ra + (dra / DEG_TO_MAS) / cosd,
        "dec": dec + ddec / DEG_TO_MAS,
        "e_a": e_gaia,
        "e_d": e_gaia,
    }
    return radio, optical, is_struct


def fetch_icrf3_gaia(*, max_arcsec: float = 0.5) -> tuple[dict, dict]:  # pragma: no cover - network
    """ICRF3 S/X (VizieR ``J/A+A/644/A159/table10``) cross-matched to Gaia DR3 via CDS X-Match.

    Returns ``(radio, optical)`` dicts with ra/dec (deg) and per-axis errors ``e_a``/``e_d`` (mas, on
    $\\alpha\\cos\\delta$ and $\\delta$). ICRF3 stores the RA error in **time-seconds** (so
    $\\sigma_{\\alpha^*}=e_\\mathrm{RA}\\times15000\\cos\\delta$ mas) and the Dec error in **arcsec**
    ($\\times1000$ mas); Gaia ``e_RA_ICRS``/``e_DE_ICRS`` are already mas. Keeps the nearest Gaia match
    per ICRF3 source within ``max_arcsec``.
    """
    import numpy as _np
    from astropy import units as _u
    from astropy.coordinates import SkyCoord
    from astropy.table import Table
    from astroquery.vizier import Vizier
    from astroquery.xmatch import XMatch

    v = Vizier(columns=["ICRF", "RAICRS", "DEICRS", "e_RAICRS", "e_DEICRS"])
    v.ROW_LIMIT = -1
    icrf = v.get_catalogs("J/A+A/644/A159/table10")[0]
    # RAICRS/DEICRS are sexagesimal (RA in hours) -> parse to decimal degrees.
    coo = SkyCoord(icrf["RAICRS"], icrf["DEICRS"], unit=(_u.hourangle, _u.deg))
    ra = _np.asarray(coo.ra.deg, float)
    dec = _np.asarray(coo.dec.deg, float)
    keep = _np.isfinite(ra) & _np.isfinite(_np.asarray(icrf["e_RAICRS"], float))
    t1 = Table(
        {
            "icrf_id": _np.arange(int(keep.sum())),
            "RAdeg": ra[keep],
            "DEdeg": dec[keep],
            "e_ra_s": _np.asarray(icrf["e_RAICRS"], float)[keep],
            "e_de_as": _np.asarray(icrf["e_DEICRS"], float)[keep],
        }
    )
    xm = XMatch.query(
        cat1=t1,
        cat2="vizier:I/355/gaiadr3",
        max_distance=max_arcsec * _u.arcsec,
        colRA1="RAdeg",
        colDec1="DEdeg",
    )
    # keep the nearest Gaia match per ICRF3 source
    xm.sort("angDist")
    _, first = _np.unique(_np.asarray(xm["icrf_id"]), return_index=True)
    xm = xm[first]
    decr = _np.asarray(xm["DEdeg"], float)
    cosd = _np.cos(_np.radians(decr))
    radio = {
        "ra": _np.asarray(xm["RAdeg"], float),
        "dec": decr,
        "e_a": _np.asarray(xm["e_ra_s"], float) * 15000.0 * cosd,
        "e_d": _np.asarray(xm["e_de_as"], float) * 1000.0,
    }
    optical = {  # Gaia: RAdeg2/DEdeg2 positions; e_RAdeg already on alpha*cos(dec), mas
        "ra": _np.asarray(xm["RAdeg2"], float),
        "dec": _np.asarray(xm["DEdeg2"], float),
        "e_a": _np.asarray(xm["e_RAdeg"], float),
        "e_d": _np.asarray(xm["e_DEdeg"], float),
    }
    return radio, optical


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Full slice: synthesise (or fetch) ICRF3×Gaia, compute offsets, reproduce the excess tail."""
    import json
    from pathlib import Path

    if offline:
        radio, optical, truth = synthetic_field()
        source = "synthetic"
    else:  # pragma: no cover - network
        radio, optical = fetch_icrf3_gaia()
        truth = None
        source = "ICRF3 x Gaia DR3"

    off = radio_optical_offset(radio["ra"], radio["dec"], optical["ra"], optical["dec"])
    sig_a = np.hypot(radio["e_a"], optical["e_a"])
    sig_d = np.hypot(radio["e_d"], optical["e_d"])
    x = normalised_offset(off["dra_mas"], off["ddec_mas"], sig_a, sig_d)
    stats = offset_statistics(x, off["offset_mas"])
    metrics = {
        "source": source,
        "n": stats["n"],
        "median_offset_mas": round(stats["median_offset_mas"], 3),
        "frac_x_gt3_pct": round(100.0 * stats["frac_x_gt_cut"], 2),
        "rayleigh_pct": round(100.0 * stats["rayleigh_expectation"], 2),
        "excess_ratio": round(stats["excess_ratio"], 1),
    }
    if (
        truth is not None
    ):  # synthetic: the excess tail should sit on the injected structured sources
        metrics["n_structured"] = int(truth.sum())
        metrics["frac_struct_in_tail"] = (
            round(float(truth[x > 3.0].mean()), 3) if (x > 3.0).any() else 0.0
        )

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "offsets_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(x, op / "papers" / "offsets" / "figures")
    _write_macros(metrics, op / "papers" / "offsets" / "generated" / "macros.tex")
    return metrics


def _figure(x_norm: np.ndarray, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    x = np.asarray(x_norm, float)
    x = x[np.isfinite(x)]
    fig, ax = plt.subplots(figsize=(5, 4))
    bins = np.linspace(0, 8, 40)
    ax.hist(x, bins=bins, density=True, color="0.6", label="ICRF3$\\times$Gaia")
    xs = 0.5 * (bins[:-1] + bins[1:])
    ax.plot(xs, xs * np.exp(-(xs**2) / 2.0), "r-", lw=1.2, label="Rayleigh (pure noise)")
    ax.axvline(3.0, color="k", ls=":", lw=0.8, label="$X=3$")
    ax.set(
        xlabel="normalised offset $X$", ylabel="density", title="Radio--optical offset significance"
    )
    ax.legend(fontsize=8)
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(out / "xnorm.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    lines = [
        "% Auto-generated by jansky_research.offsets._write_macros — do not edit by hand.",
        rf"\newcommand{{\offSource}}{{{m['source']}}}",
        rf"\newcommand{{\offN}}{{{m['n']}}}",
        rf"\newcommand{{\offMedian}}{{{m['median_offset_mas']}}}",
        rf"\newcommand{{\offFracTail}}{{{m['frac_x_gt3_pct']}}}",
        rf"\newcommand{{\offRayleigh}}{{{m['rayleigh_pct']}}}",
        rf"\newcommand{{\offExcess}}{{{m['excess_ratio']}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Radio-optical offsets of AGN (ICRF3 x Gaia DR3).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
