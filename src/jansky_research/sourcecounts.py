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
    "clustering_variance_pct",
    "compute_counts",
    "fetch_nvss_region",
    "hopkins2003_counts",
    "merge_components",
    "run",
    "synthetic_sky",
]

#: The Hopkins et al. (2003) fit's own quality: "The residuals from the sixth order fit above
#: have an rms of about 0.04 in the logarithm of the normalised counts" (their sec. 4). The
#: reference curve is not exact, and any agreement budget must carry this term.
HOPKINS2003_RESID_DEX = 0.04
#: NVSS angular correlation function amplitude range at 3.5 mJy, w(theta) ~ A theta_deg^-0.8
#: (Blake & Wall 2002); the range brackets the published fits. Used only for the
#: cosmic-variance term on the count NORMALISATION (clustering shifts all flux bins together).
NVSS_WTHETA_AMP_RANGE = (1.0e-3, 1.6e-3)

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
    bin_dex: float = 0.2,
) -> dict:
    """Euclidean-normalised differential counts of a flux-limited sample, vs the Hopkins reference.

    Cuts the sample at ``s_min_jy`` (the completeness cut), bins it on a FIXED logarithmic grid of
    ``bin_dex`` anchored at the cut --- previously the edges were ``geomspace`` to the brightest
    source, so every reported statistic was a function of one Poisson-random object (dropping the
    two brightest sources, both excluded from the comparison anyway, swung the quoted scatter by a
    factor of two). It uses ``jansky.sourcecounts`` to form the differential count, divides by the
    survey solid angle ``area_sr``, and Euclidean-normalises. The slope and the Hopkins ratio
    (median and dex scatter) are computed over bins with at least five sources **and below 1 Jy**
    (the published Hopkins validity limit). The per-bin Poisson errors are not merely returned:
    ``chi2_unity``/``chi2_with_ref_resid`` weight the comparison by them, with the second folding
    in the Hopkins fit's own quoted 0.04 dex residual --- the agreement is a statistic, not an
    impression. Returns the bin arrays, counts, errors, reference, and summary statistics.
    """
    s = np.asarray(fluxes_jy, float)
    s = s[np.isfinite(s) & (s > s_min_jy)]
    nan = float("nan")
    if s.size < 10:
        return {"n_sources": int(s.size), "centres": np.array([]), "ratio_med": None}
    n_edges = int(np.ceil((np.log10(s.max()) - np.log10(s_min_jy)) / bin_dex)) + 1
    bins = 10.0 ** (np.log10(s_min_jy) + bin_dex * np.arange(n_edges + 1))
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
    with np.errstate(divide="ignore", invalid="ignore"):
        z_unity = (en[good] - ref[good]) / en_err[good]
        ref_err = ref[good] * (10.0**HOPKINS2003_RESID_DEX - 1.0)
        z_ref = (en[good] - ref[good]) / np.hypot(en_err[good], ref_err)
        # the mean per-bin Poisson scatter in dex: the term the observed scatter should carry
        poisson_dex = float(np.mean(en_err[good] / en[good] / np.log(10.0)))
    n_good = int(good.sum())
    return {
        "n_sources": int(s.size),
        "centres": centres,
        "en": en,
        "en_err": en_err,
        "ref": ref,
        "per_bin": per_bin,
        "good": good,
        "s_min_jy": float(s_min_jy),
        "s_max_jy": float(s.max()),
        "bin_dex": float(bin_dex),
        "n_bins_used": n_good,
        "slope_diff": float(slope),
        "ratio_med": float(np.median(ratio)) if ratio.size else None,
        "ratio_scatter_dex": float(np.std(np.log10(ratio))) if ratio.size >= 2 else None,
        "chi2_unity": float(np.sum(z_unity[np.isfinite(z_unity)] ** 2)) if n_good else None,
        "chi2_with_ref_resid": float(np.sum(z_ref[np.isfinite(z_ref)] ** 2)) if n_good else None,
        "poisson_scatter_dex": round(poisson_dex, 3) if n_good else None,
        "expected_scatter_dex": (
            round(float(np.hypot(poisson_dex, HOPKINS2003_RESID_DEX)), 3) if n_good else None
        ),
    }


def merge_components(
    ra_deg: np.ndarray, dec_deg: np.ndarray, fluxes_jy: np.ndarray, *, link_arcsec: float
) -> np.ndarray:
    """Friends-of-friends merge of catalogue components into sources (fluxes summed).

    An NVSS row is a component, not a source: the 45'' beam splits extended sources, and the
    pair-count excess at 60--90'' shows ~2.4% of rows are extra components of an already-counted
    source. A blind FoF over-merges (it links genuinely distinct close neighbours too), so the
    merged counts bound the systematic from the conservative side. Returns the merged flux list.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(np.asarray(ra_deg, float) * u.deg, np.asarray(dec_deg, float) * u.deg)
    i1, i2, _, _ = c.search_around_sky(c, link_arcsec * u.arcsec)
    n = len(c)
    parent = np.arange(n)

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = int(parent[a])
        return a

    for a, b in zip(np.asarray(i1), np.asarray(i2), strict=True):
        if a == b:
            continue
        ra_, rb = find(int(a)), find(int(b))
        if ra_ != rb:
            parent[rb] = ra_
    roots = np.asarray([find(int(i)) for i in range(n)])
    f = np.asarray(fluxes_jy, float)
    sums = np.zeros(n)
    np.add.at(sums, roots, f)
    return sums[np.unique(roots)]


def clustering_variance_pct(
    radius_deg: float,
    *,
    amp: float,
    gamma: float = 0.8,
    n_pairs: int = 400_000,
    theta_min_deg: float = 0.0125,
    seed: int = 0,
) -> float:
    r"""Fractional cosmic-variance on the count normalisation for a circular field (percent).

    Clustering shifts every flux bin together, so its home is the count's NORMALISATION (the
    median ratio), not the bin-to-bin scatter: :math:`\sigma^2_{\rm cv} = \bar w`, the mean of
    the angular correlation function :math:`w(\theta)=A\,\theta_{\rm deg}^{-\gamma}` over pairs
    of positions in the field, Monte-Carlo integrated (``theta_min_deg`` caps the divergent
    small-separation samples at the survey beam scale).
    """
    rng = np.random.default_rng(seed)
    r1 = radius_deg * np.sqrt(rng.uniform(0, 1, n_pairs))
    a1 = rng.uniform(0, 2 * np.pi, n_pairs)
    r2 = radius_deg * np.sqrt(rng.uniform(0, 1, n_pairs))
    a2 = rng.uniform(0, 2 * np.pi, n_pairs)
    dx = r1 * np.cos(a1) - r2 * np.cos(a2)
    dy = r1 * np.sin(a1) - r2 * np.sin(a2)
    theta = np.maximum(np.hypot(dx, dy), theta_min_deg)
    return round(100.0 * float(np.sqrt(np.mean(amp * theta**-gamma))), 2)


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
    coo = SkyCoord(tab["RAJ2000"], tab["DEJ2000"], unit=(u.hourangle, u.deg))
    ok = np.isfinite(s_jy)
    area_sr = 2.0 * np.pi * (1.0 - np.cos(np.radians(radius_deg)))
    return {
        "fluxes_jy": s_jy[ok],
        "ra_deg": np.asarray(coo.ra.deg, float)[ok],
        "dec_deg": np.asarray(coo.dec.deg, float)[ok],
        "center_ra_deg": float(ra_deg),
        "center_dec_deg": float(dec_deg),
        "radius_deg": float(radius_deg),
        "area_sr": float(area_sr),
    }


def run(
    out: str = ".",
    *,
    offline: bool = True,
    ra: float = 180.0,
    dec: float = 30.0,
    radius_deg: float = 5.0,
    s_min_mjy: float = 3.5,
    bin_dex: float = 0.2,
) -> dict:
    """Full slice: build the NVSS Euclidean-normalised source counts and compare to Hopkins 2003."""
    from pathlib import Path

    if offline:
        sky = synthetic_sky(area_sr=2.0 * np.pi * (1.0 - np.cos(np.radians(radius_deg))))
        source = "synthetic"
    else:  # pragma: no cover - network
        sky = fetch_nvss_region(ra, dec, radius_deg)
        source = f"NVSS cone ({ra:.1f}, {dec:+.1f}) r={radius_deg:.1f} deg"

    res = compute_counts(
        sky["fluxes_jy"], sky["area_sr"], s_min_jy=s_min_mjy / 1000.0, bin_dex=bin_dex
    )
    slope = res.get("slope_diff")
    scatter = res.get("ratio_scatter_dex")
    metrics: dict = {
        "source": source,
        "n_sources": res["n_sources"],
        "area_sr": round(sky["area_sr"], 4),
        "s_min_mjy": round(s_min_mjy, 2),
        "s_max_jy": round(res["s_max_jy"], 2) if "s_max_jy" in res else None,
        "bin_dex": bin_dex,
        "n_bins_used": res.get("n_bins_used"),
        "slope_diff": round(slope, 2) if slope is not None and np.isfinite(slope) else None,
        "hopkins_ratio_med": round(res["ratio_med"], 3)
        if res.get("ratio_med") is not None
        else None,
        "hopkins_scatter_dex": round(scatter, 3) if scatter is not None else None,
        "chi2_unity": round(res["chi2_unity"], 1) if res.get("chi2_unity") is not None else None,
        "chi2_with_ref_resid": (
            round(res["chi2_with_ref_resid"], 1)
            if res.get("chi2_with_ref_resid") is not None
            else None
        ),
        "poisson_scatter_dex": res.get("poisson_scatter_dex"),
        "hopkins_resid_dex": HOPKINS2003_RESID_DEX,
        "expected_scatter_dex": res.get("expected_scatter_dex"),
    }
    # the per-bin table: a census committed as scalars is not auditable (the innerrc lesson)
    if len(res.get("centres", [])):
        metrics["bins"] = [
            {
                "centre_jy": round(float(res["centres"][i]), 5),
                "n": int(res["per_bin"][i]),
                "en": round(float(res["en"][i]), 1),
                "en_err": round(float(res["en_err"][i]), 1),
                "ref": round(float(res["ref"][i]), 1),
                "ratio": (
                    round(float(res["en"][i] / res["ref"][i]), 3) if res["ref"][i] > 0 else None
                ),
                "used": bool(res["good"][i]),
            }
            for i in range(len(res["centres"]))
        ]
    # Poisson bootstrap on the committed bin counts: the errors the headline never carried
    rng = np.random.default_rng(0)
    if res.get("ratio_med") is not None:
        g = res["good"]
        n_g = res["per_bin"][g].astype(float)
        en_per_src = res["en"][g] / np.maximum(n_g, 1)
        meds, slopes = [], []
        lc = np.log10(res["centres"][g])
        for _ in range(2000):
            n_b = rng.poisson(n_g)
            en_b = en_per_src * n_b
            ok = (n_b >= 5) & (en_b > 0)
            if ok.sum() >= 2:
                meds.append(np.median(en_b[ok] / res["ref"][g][ok]))
                dn_b = en_b[ok] * res["centres"][g][ok] ** -2.5
                slopes.append(np.polyfit(lc[ok], np.log10(dn_b), 1)[0])
        metrics["ratio_med_boot_err"] = round(float(np.std(meds)), 3) if meds else None
        metrics["slope_boot_err"] = round(float(np.std(slopes)), 3) if slopes else None
    # the recover-a-known's own comparand: Hopkins through the identical estimator/bins
    if res.get("n_bins_used"):
        g = res["good"]
        ref_dnds = res["ref"][g] * res["centres"][g] ** -2.5
        metrics["slope_ref_same_bins"] = round(
            float(np.polyfit(np.log10(res["centres"][g]), np.log10(ref_dnds), 1)[0]), 2
        )
    # completeness-cut sweep: the paper may claim stability only after running it
    sweep = {}
    for cut in (3.5, 5.0, 7.0, 10.0, 20.0):
        r = compute_counts(sky["fluxes_jy"], sky["area_sr"], s_min_jy=cut / 1000.0, bin_dex=bin_dex)
        sweep[f"{cut:g}"] = {
            "n": r.get("n_sources"),
            "ratio_med": round(r["ratio_med"], 3) if r.get("ratio_med") is not None else None,
            "slope": (
                round(r["slope_diff"], 2)
                if r.get("slope_diff") is not None and np.isfinite(r.get("slope_diff", np.nan))
                else None
            ),
        }
    metrics["cut_sweep_mjy"] = sweep
    # cosmic variance on the NORMALISATION (clustering moves all bins together): the term the
    # old prose mis-assigned to the bin-to-bin scatter
    metrics["cosmic_variance_pct_lo"] = clustering_variance_pct(
        radius_deg, amp=NVSS_WTHETA_AMP_RANGE[0]
    )
    metrics["cosmic_variance_pct_hi"] = clustering_variance_pct(
        radius_deg, amp=NVSS_WTHETA_AMP_RANGE[1]
    )
    metrics["poisson_norm_pct"] = round(100.0 / np.sqrt(max(res["n_sources"], 1)), 2)

    if not offline and sky.get("ra_deg") is not None:  # pragma: no cover - real leg only
        # near-threshold completeness profile in fine bins: is the cut above the roll-off?
        fine_edges = np.array([2.1, 2.5, 3.0, 3.5, 4.2, 5.0, 6.0, 7.2]) / 1000.0
        s_all = np.asarray(sky["fluxes_jy"], float)
        prof = []
        for lo, hi in zip(fine_edges[:-1], fine_edges[1:], strict=False):
            msk = (s_all >= lo) & (s_all < hi)
            n_i = int(msk.sum())
            cen = float(np.sqrt(lo * hi))
            dnds = n_i / (hi - lo) / sky["area_sr"]
            ratio_i = (cen**2.5 * dnds) / float(hopkins2003_counts(np.array([cen]))[0])
            prof.append(
                {
                    "lo_mjy": round(1000 * lo, 1),
                    "hi_mjy": round(1000 * hi, 1),
                    "n": n_i,
                    "ratio": round(ratio_i, 3),
                    "err": round(ratio_i / np.sqrt(max(n_i, 1)), 3),
                }
            )
        metrics["threshold_profile"] = prof
        # component-vs-source systematic: friends-of-friends merging sweep
        fof = {}
        for link in (60.0, 80.0, 100.0):
            merged = merge_components(sky["ra_deg"], sky["dec_deg"], s_all, link_arcsec=link)
            r = compute_counts(merged, sky["area_sr"], s_min_jy=s_min_mjy / 1000.0, bin_dex=bin_dex)
            fof[f"{link:g}"] = {
                "n_merged": int(len(s_all) - len(merged)),
                "n_above_cut": r.get("n_sources"),
                "ratio_med": round(r["ratio_med"], 3) if r.get("ratio_med") is not None else None,
            }
        metrics["fof_sweep_arcsec"] = fof
        # surface-density uniformity across four equal-area annuli of the cone
        import astropy.units as u
        from astropy.coordinates import SkyCoord

        centre = SkyCoord(sky["center_ra_deg"] * u.deg, sky["center_dec_deg"] * u.deg)
        seps = centre.separation(SkyCoord(sky["ra_deg"] * u.deg, sky["dec_deg"] * u.deg)).deg
        above = s_all > s_min_mjy / 1000.0
        edges = sky["radius_deg"] * np.sqrt(np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
        dens = []
        for lo, hi in zip(edges[:-1], edges[1:], strict=False):
            in_ann = above & (seps >= lo) & (seps < hi)
            ann_sr = 2 * np.pi * (np.cos(np.radians(lo)) - np.cos(np.radians(hi)))
            dens.append(round(float(in_ann.sum() / ann_sr), 0))
        metrics["annulus_density_per_sr"] = dens

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "sourcecounts_metrics.json")
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
        good = res.get("good", np.ones(c.size, bool))
        ax.errorbar(
            c[good],
            res["en"][good],
            yerr=res["en_err"][good],
            fmt="o",
            color="C0",
            ms=4,
            capsize=2,
            label="NVSS (compared bins)",
        )
        if (~good).any():
            ax.errorbar(
                c[~good],
                res["en"][~good],
                yerr=res["en_err"][~good],
                fmt="o",
                mfc="none",
                color="C0",
                ms=4,
                capsize=2,
                label="excluded ($N<5$ or $>1$ Jy)",
            )
        # solid only inside the published validity; dashed = extrapolation, never assessed
        sm = np.geomspace(c.min(), min(c.max(), HOPKINS2003_SMAX_JY), 100)
        ax.plot(sm, hopkins2003_counts(sm), "-", color="C3", lw=1.5, label="Hopkins et al. 2003")
        if c.max() > HOPKINS2003_SMAX_JY:
            sx = np.geomspace(HOPKINS2003_SMAX_JY, c.max(), 40)
            ax.plot(sx, hopkins2003_counts(sx), "--", color="C3", lw=1.0, alpha=0.6)
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
        # the uncertainties and the computed agreement budget (round-8): the headline carries
        # errors, the scatter has an expectation, and the agreement is a chi-square
        rf"\newcommand{{\scRatioErr}}{{{_fmt('ratio_med_boot_err')}}}",
        rf"\newcommand{{\scSlopeErr}}{{{_fmt('slope_boot_err')}}}",
        rf"\newcommand{{\scSlopeRef}}{{{_fmt('slope_ref_same_bins')}}}",
        rf"\newcommand{{\scPoissonDex}}{{{_fmt('poisson_scatter_dex')}}}",
        rf"\newcommand{{\scExpectedDex}}{{{_fmt('expected_scatter_dex')}}}",
        rf"\newcommand{{\scChiUnity}}{{{_fmt('chi2_unity')}}}",
        rf"\newcommand{{\scChiRef}}{{{_fmt('chi2_with_ref_resid')}}}",
        rf"\newcommand{{\scCvLo}}{{{_fmt('cosmic_variance_pct_lo')}}}",
        rf"\newcommand{{\scCvHi}}{{{_fmt('cosmic_variance_pct_hi')}}}",
        rf"\newcommand{{\scPoissonNormPct}}{{{_fmt('poisson_norm_pct')}}}",
    ]
    # real-leg systematics (placeholders offline): threshold profile edges, FoF, annuli
    fof = m.get("fof_sweep_arcsec") or {}
    for link, name in (("60", "Sixty"), ("80", "Eighty"), ("100", "Hundred")):
        cell = fof.get(link) or {}
        v = cell.get("ratio_med")
        lines.append(rf"\newcommand{{\scFof{name}Ratio}}{{{'--' if v is None else v}}}")
    prof = m.get("threshold_profile") or []
    first_used = next((p_ for p_ in prof if p_["lo_mjy"] >= 3.5), None)
    lines.append(
        rf"\newcommand{{\scFirstBinResid}}"
        rf"{{{'--' if first_used is None else first_used['ratio']}}}"
    )
    dens = m.get("annulus_density_per_sr") or []
    if dens:
        spread = 100.0 * (max(dens) - min(dens)) / (2.0 * (sum(dens) / len(dens)))
        lines.append(rf"\newcommand{{\scAnnulusSpreadPct}}{{{round(spread, 1)}}}")
    else:
        lines.append(r"\newcommand{\scAnnulusSpreadPct}{--}")
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

    p = argparse.ArgumentParser(description="NVSS 1.4 GHz Euclidean-normalised source counts.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--ra", type=float, default=180.0)
    p.add_argument("--dec", type=float, default=30.0)
    p.add_argument("--radius", type=float, default=5.0)
    p.add_argument("--s-min", type=float, default=3.5, help="completeness cut (mJy)")
    p.add_argument("--bin-dex", type=float, default=0.2)
    args = p.parse_args(argv)
    metrics = run(
        args.out,
        offline=args.offline,
        ra=args.ra,
        dec=args.dec,
        radius_deg=args.radius,
        s_min_mjy=args.s_min,
        bin_dex=args.bin_dex,
    )
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
