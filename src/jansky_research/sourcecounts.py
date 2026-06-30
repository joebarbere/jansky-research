"""Euclidean-normalised radio source counts from NVSS, against the canonical 1.4 GHz reference.

The differential source count :math:`\\mathrm{d}N/\\mathrm{d}S` --- how many radio sources there are per
unit flux per unit sky --- is one of the oldest cosmological tests in radio astronomy: a static
Euclidean universe gives :math:`\\mathrm{d}N/\\mathrm{d}S\\propto S^{-5/2}`, so plotting
:math:`S^{5/2}\\,\\mathrm{d}N/\\mathrm{d}S` (the *Euclidean-normalised* count) flattens that slope out and
any real structure --- the bright-end steepening, the sub-mJy upturn from star-forming galaxies ---
stands out. This slice builds that count from a public NVSS region and compares it to the canonical
1.4 GHz counts (the \\citet{hopkins2003} polynomial fit), as a reproducible recover-a-known.

It reuses the ``jansky.sourcecounts`` helpers (``differential_counts``, ``euclidean_normalised_counts``,
``count_slope``, ``integral_counts``) wholesale; the new code is the NVSS region fetch, the solid-angle
normalisation, and the published reference curve. Pure NumPy/astropy with a synthetic offline fixture;
the real fetch (VizieR NVSS) is network-gated.
"""

from __future__ import annotations

import numpy as np
from jansky import sourcecounts as jsc

__all__ = [
    "compute_counts",
    "fetch_nvss_region",
    "hopkins2003_counts",
    "run",
    "synthetic_sky",
]

# Hopkins et al. (2003, AJ 125, 465) 6th-order polynomial fit to the 1.4 GHz Euclidean-normalised
# differential source counts: log10(S^2.5 dN/dS / [Jy^1.5 sr^-1]) = sum a_i (log10(S/mJy))^i.
HOPKINS2003_COEFFS = (0.859, 0.508, 0.376, -0.049, -0.121, 0.057, -0.008)
HOPKINS2003_SMAX_JY = 1.0  # the published fit is valid only to 1 Jy (0.05--1000 mJy)


def hopkins2003_counts(s_jy: np.ndarray) -> np.ndarray:
    """Canonical 1.4 GHz Euclidean-normalised differential counts (Hopkins et al. 2003).

    Evaluates the published 6th-order polynomial in :math:`x=\\log_{10}(S/\\mathrm{mJy})` and returns
    :math:`S^{5/2}\\,\\mathrm{d}N/\\mathrm{d}S` in Jy\\ :sup:`1.5`\\ sr\\ :sup:`-1`. The fit is valid over
    **0.05--1000 mJy** (50 µJy to 1 Jy) per Hopkins et al. 2003; outside that range the polynomial is an
    extrapolation and not meaningful.
    """
    x = np.log10(np.asarray(s_jy, float) * 1000.0)  # S in mJy
    logc = np.polyval(list(reversed(HOPKINS2003_COEFFS)), x)
    return np.asarray(10.0**logc, float)


def compute_counts(
    fluxes_jy: np.ndarray,
    area_sr: float,
    *,
    s_min_jy: float,
    n_bins: int = 12,
) -> dict:
    """Euclidean-normalised differential counts of a flux-limited sample, vs the Hopkins reference.

    Cuts the sample at ``s_min_jy`` (the completeness limit), bins it logarithmically out to the
    brightest source, and uses ``jansky.sourcecounts`` to form the differential count, divide by the
    survey solid angle ``area_sr`` to get :math:`\\mathrm{d}N/\\mathrm{d}S` per steradian, and
    Euclidean-normalise it. The differential log--log slope and the ratio to :func:`hopkins2003_counts`
    (median and dex scatter) are computed over bins with at least five sources **and below 1 Jy** (the
    published Hopkins validity limit), so neither is contaminated by the Poisson-starved, extrapolated
    bright end. Returns the bin centres (Jy), the normalised counts and Poisson errors
    (Jy\\ :sup:`1.5`\\ sr\\ :sup:`-1`), the per-bin source counts, the slope, and the Hopkins ratio.
    """
    s = np.asarray(fluxes_jy, float)
    s = s[np.isfinite(s) & (s > s_min_jy)]
    nan = float("nan")
    if s.size < 10:
        return {"n_sources": int(s.size), "centres": np.array([]), "ratio_med": None}
    bins = np.geomspace(s_min_jy, s.max() * 1.001, n_bins + 1)
    centres, dn_ds, dn_ds_err = jsc.differential_counts(s, bins)
    per_bin, _ = np.histogram(s, bins)
    dn_ds_sr = dn_ds / area_sr
    err_sr = dn_ds_err / area_sr
    en = jsc.euclidean_normalised_counts(centres, dn_ds_sr)
    en_err = jsc.euclidean_normalised_counts(centres, err_sr)
    ref = hopkins2003_counts(centres)
    # compare/fit only where bins are populated AND within the Hopkins fit's validity (< 1 Jy)
    good = (per_bin >= 5) & (centres < HOPKINS2003_SMAX_JY)
    slope = jsc.count_slope(centres[good], dn_ds_sr[good]) if good.sum() >= 2 else nan
    ratio = en[good] / ref[good]
    ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
    return {
        "n_sources": int(s.size),
        "centres": centres,
        "en": en,
        "en_err": en_err,
        "ref": ref,
        "per_bin": per_bin,
        "s_min_jy": float(s_min_jy),
        "s_max_jy": float(s.max()),
        "n_bins_used": int(good.sum()),
        "slope_diff": float(slope),
        "ratio_med": float(np.median(ratio)) if ratio.size else None,
        "ratio_scatter_dex": float(np.std(np.log10(ratio))) if ratio.size >= 2 else None,
    }


def synthetic_sky(
    *,
    area_sr: float = 0.05,
    s_min_jy: float = 0.0035,
    s_max_jy: float = 5.0,
    seed: int = 0,
) -> dict:
    """Synthetic flux-limited sky drawn from the Hopkins 2003 differential counts (rejection sampling).

    Draws source fluxes whose differential count follows :math:`\\mathrm{d}N/\\mathrm{d}S =
    S^{-5/2}\\,[S^{5/2}\\mathrm{d}N/\\mathrm{d}S]_{\\rm Hopkins}` over ``[s_min_jy, s_max_jy]`` and an area
    ``area_sr``, so the pipeline run on it recovers the Hopkins curve (median ratio ≈ 1). Returns the
    fluxes (Jy) and the area (sr).
    """
    rng = np.random.default_rng(seed)
    grid = np.geomspace(s_min_jy, s_max_jy, 4000)
    dn_ds = hopkins2003_counts(grid) * grid**-2.5  # per Jy per sr
    # expected total number over the area: integrate dN/dS dS * area
    n_exp = float(np.trapezoid(dn_ds, grid) * area_sr)
    n_draw = int(rng.poisson(n_exp))
    # inverse-CDF sampling on the (un-normalised) dN/dS
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (dn_ds[1:] + dn_ds[:-1]) * np.diff(grid))])
    cdf /= cdf[-1]
    u = rng.uniform(0.0, 1.0, n_draw)
    fluxes = np.interp(u, cdf, grid)
    return {"fluxes_jy": fluxes, "area_sr": float(area_sr)}


def fetch_nvss_region(
    ra_deg: float, dec_deg: float, radius_deg: float
) -> dict:  # pragma: no cover - network
    """Fetch NVSS sources in a cone from VizieR (Condon et al. 1998, VIII/65) → fluxes (Jy) + area (sr).

    Returns every NVSS source within ``radius_deg`` of (``ra_deg``, ``dec_deg``) with its 1.4 GHz
    integrated flux density (``S1.4``, mJy → Jy) and the cone solid angle
    :math:`2\\pi(1-\\cos\\theta)`. Needs network access (astroquery/VizieR).
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    v = Vizier(columns=["RAJ2000", "DEJ2000", "S1.4"], row_limit=-1)
    center = SkyCoord(ra_deg, dec_deg, unit="deg")
    tab = v.query_region(center, radius=radius_deg * u.deg, catalog="VIII/65/nvss")[0]
    s_jy = np.asarray(tab["S1.4"], float) / 1000.0  # mJy -> Jy
    area_sr = 2.0 * np.pi * (1.0 - np.cos(np.radians(radius_deg)))
    return {"fluxes_jy": s_jy[np.isfinite(s_jy)], "area_sr": float(area_sr)}


def run(
    out: str = ".",
    *,
    offline: bool = True,
    ra: float = 180.0,
    dec: float = 30.0,
    radius_deg: float = 5.0,
    s_min_mjy: float = 3.5,
    n_bins: int = 12,
) -> dict:
    """Full slice: build the NVSS Euclidean-normalised source counts and compare to Hopkins 2003."""
    import json
    from pathlib import Path

    if offline:
        sky = synthetic_sky(area_sr=2.0 * np.pi * (1.0 - np.cos(np.radians(radius_deg))))
        source = "synthetic"
    else:  # pragma: no cover - network
        sky = fetch_nvss_region(ra, dec, radius_deg)
        source = f"NVSS cone ({ra:.1f}, {dec:+.1f}) r={radius_deg:.1f} deg"

    res = compute_counts(
        sky["fluxes_jy"], sky["area_sr"], s_min_jy=s_min_mjy / 1000.0, n_bins=n_bins
    )
    slope = res.get("slope_diff")
    scatter = res.get("ratio_scatter_dex")
    metrics: dict = {
        "source": source,
        "n_sources": res["n_sources"],
        "area_sr": round(sky["area_sr"], 4),
        "s_min_mjy": round(s_min_mjy, 2),
        "s_max_jy": round(res["s_max_jy"], 2) if "s_max_jy" in res else None,
        "n_bins_used": res.get("n_bins_used"),
        "slope_diff": round(slope, 2) if slope is not None and np.isfinite(slope) else None,
        "hopkins_ratio_med": round(res["ratio_med"], 3)
        if res.get("ratio_med") is not None
        else None,
        "hopkins_scatter_dex": round(scatter, 3) if scatter is not None else None,
    }

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "sourcecounts_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(res, op / "papers" / "sourcecounts" / "figures")
    _write_macros(metrics, op / "papers" / "sourcecounts" / "generated" / "macros.tex")
    return metrics


def _figure(res: dict, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    if len(res.get("centres", [])):
        c = res["centres"]
        ax.errorbar(
            c, res["en"], yerr=res["en_err"], fmt="o", color="C0", ms=4, capsize=2, label="NVSS"
        )
        sm = np.geomspace(c.min(), c.max(), 100)
        ax.plot(sm, hopkins2003_counts(sm), "-", color="C3", lw=1.5, label="Hopkins et al. 2003")
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="flux density (Jy)",
        ylabel=r"$S^{5/2}\,\mathrm{d}N/\mathrm{d}S$ (Jy$^{1.5}$ sr$^{-1}$)",
        title="1.4 GHz Euclidean-normalised counts",
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "sourcecounts.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.sourcecounts._write_macros -- do not edit by hand.",
        rf"\newcommand{{\scSource}}{{{m['source']}}}",
        rf"\newcommand{{\scNsrc}}{{{_fmt('n_sources')}}}",
        rf"\newcommand{{\scArea}}{{{_fmt('area_sr')}}}",
        rf"\newcommand{{\scSmin}}{{{_fmt('s_min_mjy')}}}",
        rf"\newcommand{{\scSmax}}{{{_fmt('s_max_jy')}}}",
        rf"\newcommand{{\scNbins}}{{{_fmt('n_bins_used')}}}",
        rf"\newcommand{{\scSlope}}{{{_fmt('slope_diff')}}}",
        rf"\newcommand{{\scRatio}}{{{_fmt('hopkins_ratio_med')}}}",
        rf"\newcommand{{\scScatter}}{{{_fmt('hopkins_scatter_dex')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="NVSS 1.4 GHz Euclidean-normalised source counts.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--ra", type=float, default=180.0)
    p.add_argument("--dec", type=float, default=30.0)
    p.add_argument("--radius", type=float, default=5.0)
    p.add_argument("--s-min", type=float, default=3.5, help="completeness cut (mJy)")
    p.add_argument("--n-bins", type=int, default=12)
    args = p.parse_args(argv)
    metrics = run(
        args.out,
        offline=args.offline,
        ra=args.ra,
        dec=args.dec,
        radius_deg=args.radius,
        s_min_mjy=args.s_min,
        n_bins=args.n_bins,
    )
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
