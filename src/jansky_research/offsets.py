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
    "alignment_stats",
    "fetch_icrf3_gaia",
    "fetch_mojave_jets",
    "jet_axis_angles",
    "match_jets",
    "normalised_offset",
    "offset_statistics",
    "radio_optical_offset",
    "run",
    "synthetic_alignment",
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


def jet_axis_angles(
    offset_pa_deg: np.ndarray, jet_pa_deg: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    r"""Angles between a radio→optical offset PA and the parsec-scale jet PA.

    Returns ``(downstream_deg, axis_deg)``: the *downstream* angle wrapped to $[0,180]$ (0 = the offset
    points along the jet toward the approaching/downstream side, 180 = anti-jet/upstream), and the
    *jet-axis* angle $\min(\Delta,180-\Delta)\in[0,90]$ (0 = aligned with the jet axis, 90 =
    perpendicular). Both PAs are degrees East of North.
    """
    d = np.abs(
        ((np.asarray(offset_pa_deg, float) - np.asarray(jet_pa_deg, float) + 180.0) % 360.0) - 180.0
    )
    return d, np.minimum(d, 180.0 - d)


def alignment_stats(
    offset_pa_deg: np.ndarray, jet_pa_deg: np.ndarray, x_norm: np.ndarray, *, x_cut: float = 2.0
) -> dict:
    r"""Test whether radio→optical offsets align with the jet (the Kovalev/Petrov/Plavin result).

    The **full matched sample** is the primary test: a Kolmogorov--Smirnov comparison of the jet-axis
    angle (:func:`jet_axis_angles`, folded to $[0,90]$) against the uniform distribution expected if
    offsets were randomly oriented (median $45°$, fraction within $30°$ of the axis $=1/3$). The
    fraction of offsets pointing *downstream* (within $45°$ of the jet direction, random $=1/4$) isolates
    the downstream component. As a *qualitative* consistency check the same statistics are reported for
    the significant ($X>$ ``x_cut``) subset, which is expected to align more tightly (weak offsets are
    astrometric noise with random PA). Returns the counts and these statistics.
    """
    from scipy import stats as _stats

    opa = np.asarray(offset_pa_deg, float)
    jpa = np.asarray(jet_pa_deg, float)
    x = np.asarray(x_norm, float)
    good = np.isfinite(opa) & np.isfinite(jpa) & np.isfinite(x)
    opa, jpa, x = opa[good], jpa[good], x[good]
    down, axis = jet_axis_angles(opa, jpa)
    nan = float("nan")
    ks_p = float(_stats.kstest(axis / 90.0, "uniform").pvalue) if axis.size >= 5 else nan
    sig = x > x_cut
    return {
        "n_jet": int(axis.size),
        "median_axis_deg": float(np.median(axis)) if axis.size else nan,
        "frac_axis_lt30": float(np.mean(axis < 30.0)) if axis.size else nan,
        "frac_down_lt45": float(np.mean(down < 45.0)) if axis.size else nan,
        "ks_p": ks_p,
        "x_cut": x_cut,
        "n_jet_signif": int(sig.sum()),
        "median_axis_signif_deg": float(np.median(axis[sig])) if sig.any() else nan,
        "frac_axis_signif": float(np.mean(axis[sig] < 30.0)) if sig.any() else nan,
        "_axis_deg": axis,  # for the figure
    }


def synthetic_alignment(
    *,
    n: int = 420,
    aligned_fraction: float = 0.6,
    downstream_fraction: float = 0.8,
    jet_scatter_deg: float = 18.0,
    seed: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Synthetic matched offset/jet sample with an injected jet-aligned, mostly-downstream population.

    A fraction ``aligned_fraction`` of sources are *structural*: their offset points along the jet (a
    ``downstream_fraction`` majority downstream, the rest anti-jet) with ``jet_scatter_deg`` of scatter,
    and they carry a significant ``X``; the remainder are astrometric noise with a random offset PA and
    low ``X``. Returns ``(offset_pa_deg, jet_pa_deg, x_norm)`` so :func:`alignment_stats` recovers the
    injected alignment offline (no network).
    """
    rng = np.random.default_rng(seed)
    jet_pa = rng.uniform(0.0, 360.0, n)
    aligned = rng.random(n) < aligned_fraction
    downstream = rng.random(n) < downstream_fraction
    along = np.where(downstream, jet_pa, jet_pa + 180.0) + rng.normal(0.0, jet_scatter_deg, n)
    offset_pa = np.where(aligned, along % 360.0, rng.uniform(0.0, 360.0, n))
    x = np.where(aligned, rng.uniform(2.5, 8.0, n), rng.uniform(0.3, 2.0, n))
    return offset_pa, jet_pa, x


def fetch_mojave_jets() -> dict:  # pragma: no cover - network
    """Per-source mean innermost jet position angle from MOJAVE XVIII (Lister et al. 2021).

    VizieR ``J/ApJ/923/30/mojave18``: ``PA`` is the flux-weighted innermost jet PA measured from the
    core toward the approaching (downstream) side, deg East of North; ``delPA`` is its range (jet
    wobble). Returns sky positions (deg), ``jet_pa`` (deg), and ``delpa`` (deg).
    """
    from astroquery.vizier import Vizier

    v = Vizier(columns=["_RA", "_DE", "PA", "delPA"], row_limit=-1)
    t = v.get_catalogs("J/ApJ/923/30")["J/ApJ/923/30/mojave18"]
    return {
        "ra": np.asarray(t["_RA"], float),
        "dec": np.asarray(t["_DE"], float),
        "jet_pa": np.asarray(t["PA"], float),
        "delpa": np.asarray(t["delPA"], float),
    }


def match_jets(
    radio_ra_deg: np.ndarray, radio_dec_deg: np.ndarray, jets: dict, *, max_arcsec: float = 1.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Positionally match radio sources to the jet catalogue → ``(mask, jet_pa, delpa)``.

    Nearest-neighbour match within ``max_arcsec`` (same AGN). ``mask`` selects matched radio sources;
    ``jet_pa``/``delpa`` are aligned to the radio array order. Pure astropy (no network).
    """
    import astropy.units as _u
    from astropy.coordinates import SkyCoord

    rc = SkyCoord(
        np.asarray(radio_ra_deg, float) * _u.deg, np.asarray(radio_dec_deg, float) * _u.deg
    )
    jc = SkyCoord(jets["ra"] * _u.deg, jets["dec"] * _u.deg)
    idx, sep, _ = rc.match_to_catalog_sky(jc)
    mask = sep.arcsec < max_arcsec
    return mask, np.asarray(jets["jet_pa"], float)[idx], np.asarray(jets["delpa"], float)[idx]


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

    # jet-alignment test: does the offset DIRECTION point along the parsec-scale jet?
    if offline:
        a_off_pa, a_jet_pa, a_x = synthetic_alignment()
        align = alignment_stats(a_off_pa, a_jet_pa, a_x)
        jet_source = "synthetic"
    else:  # pragma: no cover - network
        jets = fetch_mojave_jets()
        mask, jet_pa, _delpa = match_jets(radio["ra"], radio["dec"], jets)
        align = alignment_stats(off["pa_deg"][mask], jet_pa[mask], x[mask])
        jet_source = "MOJAVE XVIII"

    metrics = {
        "source": source,
        "n": stats["n"],
        "median_offset_mas": round(stats["median_offset_mas"], 3),
        "frac_x_gt3_pct": round(100.0 * stats["frac_x_gt_cut"], 2),
        "rayleigh_pct": round(100.0 * stats["rayleigh_expectation"], 2),
        "excess_ratio": round(stats["excess_ratio"], 1),
        "jet_source": jet_source,
        "n_jet": align["n_jet"],
        "median_axis_deg": round(align["median_axis_deg"], 1),
        "frac_axis_lt30": round(align["frac_axis_lt30"], 3),
        "frac_down_lt45": round(align["frac_down_lt45"], 3),
        "ks_p": align["ks_p"],
        "n_jet_signif": align["n_jet_signif"],
        "median_axis_signif_deg": round(align["median_axis_signif_deg"], 1),
        "frac_axis_signif": round(align["frac_axis_signif"], 3),
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
    _figure(x, align["_axis_deg"], op / "papers" / "offsets" / "figures")
    _write_macros(metrics, op / "papers" / "offsets" / "generated" / "macros.tex")
    return metrics


def _figure(x_norm: np.ndarray, axis_deg: np.ndarray, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.0))

    # Left: the offset-magnitude excess over Rayleigh (the existing result)
    x = np.asarray(x_norm, float)
    x = x[np.isfinite(x)]
    bins = np.linspace(0, 8, 40)
    ax1.hist(x, bins=bins, density=True, color="0.6", label="ICRF3$\\times$Gaia")
    xs = 0.5 * (bins[:-1] + bins[1:])
    ax1.plot(xs, xs * np.exp(-(xs**2) / 2.0), "r-", lw=1.2, label="Rayleigh (pure noise)")
    ax1.axvline(3.0, color="k", ls=":", lw=0.8, label="$X=3$")
    ax1.set(xlabel="normalised offset $X$", ylabel="density", title="Offset significance")
    ax1.legend(fontsize=8)
    ax1.set_yscale("log")

    # Right: the offset DIRECTION vs the jet axis (the new result) -- a peak at 0 = aligned
    a = np.asarray(axis_deg, float)
    a = a[np.isfinite(a)]
    abins = np.linspace(0, 90, 19)
    ax2.hist(a, bins=abins, density=True, color="C0", label="offset vs jet")
    ax2.axhline(1.0 / 90.0, color="r", ls="-", lw=1.2, label="uniform (random)")
    ax2.set(
        xlabel="jet-axis angle (deg)",
        ylabel="density",
        title="Offset direction vs parsec-scale jet",
    )
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "xnorm.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _texp(p: float) -> str:
        """Format a (tiny) p-value as a LaTeX exponent, e.g. ``3\\times10^{-22}``."""
        if not np.isfinite(p) or p <= 0:
            return "--"
        e = int(np.floor(np.log10(p)))
        mant = p / 10.0**e
        return rf"{mant:.0f}\times10^{{{e}}}" if e < -2 else f"{p:.3f}"

    lines = [
        "% Auto-generated by jansky_research.offsets._write_macros — do not edit by hand.",
        rf"\newcommand{{\offSource}}{{{m['source']}}}",
        rf"\newcommand{{\offN}}{{{m['n']}}}",
        rf"\newcommand{{\offMedian}}{{{m['median_offset_mas']}}}",
        rf"\newcommand{{\offFracTail}}{{{m['frac_x_gt3_pct']}}}",
        rf"\newcommand{{\offRayleigh}}{{{m['rayleigh_pct']}}}",
        rf"\newcommand{{\offExcess}}{{{m['excess_ratio']}}}",
        rf"\newcommand{{\offJetSource}}{{{m['jet_source']}}}",
        rf"\newcommand{{\offJetN}}{{{m['n_jet']}}}",
        rf"\newcommand{{\offJetMedAxis}}{{{m['median_axis_deg']}}}",
        rf"\newcommand{{\offJetFracAxis}}{{{round(100.0 * m['frac_axis_lt30'])}}}",
        rf"\newcommand{{\offJetFracDown}}{{{round(100.0 * m['frac_down_lt45'])}}}",
        rf"\newcommand{{\offJetKsP}}{{{_texp(m['ks_p'])}}}",
        rf"\newcommand{{\offJetNsig}}{{{m['n_jet_signif']}}}",
        rf"\newcommand{{\offJetMedAxisSig}}{{{m['median_axis_signif_deg']}}}",
        rf"\newcommand{{\offJetFracAxisSig}}{{{round(100.0 * m['frac_axis_signif'])}}}",
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
