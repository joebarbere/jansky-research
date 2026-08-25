"""Pulsar radio spectral indices from the ATNF Pulsar Catalogue.

Pulsars have steep radio spectra: $S_\\nu\\propto\\nu^\\alpha$ with a mean two-frequency index
$\\alpha\\approx-1.6$ to $-1.8$ (Maron et al. 2000; Bates et al. 2013; Jankowski et al. 2018), far
steeper than the $\\sim-0.7$ of ordinary synchrotron sources. Using the ATNF Pulsar Catalogue's
tabulated 400 and 1400 MHz flux densities (Manchester et al. 2005) we reproduce that distribution and
compare **millisecond** ($P<30$ ms) with **normal** pulsars.

Catalogue-only and maximal-reuse: composes :func:`jansky_research.spectra.spectral_index` for the
400$\\to$1400 MHz two-point index. Pure NumPy + a synthetic offline fixture for tests.
"""

from __future__ import annotations

import numpy as np

from .spectra import spectral_index

__all__ = [
    "MSP_PERIOD_MAX_S",
    "NU_GHZ",
    "completeness_limited",
    "fetch_atnf",
    "find_spectra",
    "is_millisecond",
    "msp_cut_sweep",
    "permutation_pvalue",
    "pulsar_alpha",
    "run",
    "spectral_distribution",
    "synthetic_field",
]

# ATNF tabulated flux frequencies (GHz).
NU_GHZ = {"s400": 0.4, "s1400": 1.4}
# Millisecond-pulsar period cut (s) — the standard P < 30 ms boundary.
MSP_PERIOD_MAX_S = 0.03


def pulsar_alpha(
    s400: np.ndarray,
    s1400: np.ndarray,
    e400: np.ndarray | None = None,
    e1400: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Two-point pulsar spectral index $\alpha^{1400}_{400}$ and its error (reuses ``spectra``)."""
    return spectral_index(s400, NU_GHZ["s400"], s1400, NU_GHZ["s1400"], e400, e1400)


def is_millisecond(period_s: np.ndarray, *, p_max: float = MSP_PERIOD_MAX_S) -> np.ndarray:
    """Boolean MSP classifier: spin period below ``p_max`` (default 30 ms)."""
    return np.asarray(period_s, float) < p_max


def spectral_distribution(alpha: np.ndarray) -> dict:
    """Summarise an $\\alpha$ distribution: ``n``, ``mean``, ``median``, ``std``, ``se_mean``.

    ``std`` uses ddof=1 (a sample dispersion), matching :func:`compare_subsamples`, so the paper
    does not present two dispersion conventions side by side as the same quantity. ``se_mean`` is
    the statistical standard error on the mean -- the paper's headline previously carried a
    scatter but no uncertainty at all.
    """
    a = np.asarray(alpha, float)
    a = a[np.isfinite(a)]
    n = int(a.size)
    nan = float("nan")
    std = float(np.std(a, ddof=1)) if n > 1 else nan
    return {
        "n": n,
        "mean": float(np.mean(a)) if n else nan,
        "median": float(np.median(a)) if n else nan,
        "std": std,
        "se_mean": std / float(np.sqrt(n)) if n > 1 else nan,
    }


def compare_subsamples(alpha_a: np.ndarray, alpha_b: np.ndarray, *, n_sigma: float = 2.0) -> dict:
    """Difference of two mean spectral indices, with the sensitivity that makes a null meaningful.

    "The two are indistinguishable" is a statement about *resolution*, not about two rounded means
    agreeing. Reporting it without a sensitivity is the error CLAUDE.md records for ``frblens``:
    a null divides by what the measurement could have seen, not by how many objects went into it.

    Returns the two means and dispersions, their difference ``diff = mean_a - mean_b``, the standard
    error on that difference (Welch, i.e. no equal-variance assumption), the difference in units of
    that error (``n_sigma_observed``), and three sensitivity numbers whose meanings must not be
    conflated: ``ci95_lo``/``ci95_hi``, the 95% confidence interval on the difference -- the values
    *outside* it are what the data exclude; ``resolvable`` = ``n_sigma * se``, the offset detected
    at ``n_sigma`` with only 50% probability (a detection threshold, not an exclusion); and
    ``resolvable_90power``, the offset detected at 2 sigma with 90% power. The earlier version of
    this docstring called ``resolvable`` an exclusion, and the paper inherited the error: an offset
    of exactly 2 SE yields a >=2 sigma result half the time.

    The Welch SE assumes independent per-source flux measurements; ATNF fluxes cluster by source
    publication, so the effective N in each arm is below the row count and the SE is a floor.
    """
    a = np.asarray(alpha_a, float)
    b = np.asarray(alpha_b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    na, nb = int(a.size), int(b.size)
    nan = float("nan")
    if na < 2 or nb < 2:  # a standard error needs at least two points in each arm
        return {
            "n_a": na,
            "n_b": nb,
            "mean_a": float(np.mean(a)) if na else nan,
            "mean_b": float(np.mean(b)) if nb else nan,
            "std_a": nan,
            "std_b": nan,
            "diff": nan,
            "se_diff": nan,
            "n_sigma_observed": nan,
            "resolvable": nan,
            "ci95_lo": nan,
            "ci95_hi": nan,
            "resolvable_90power": nan,
        }
    # ddof=1: these are samples, and with n_a = 43 the difference between the two conventions is
    # not negligible in the standard error the null is quoted against.
    sa, sb = float(np.std(a, ddof=1)), float(np.std(b, ddof=1))
    se = float(np.hypot(sa / np.sqrt(na), sb / np.sqrt(nb)))
    diff = float(np.mean(a) - np.mean(b))
    return {
        "n_a": na,
        "n_b": nb,
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "std_a": sa,
        "std_b": sb,
        "diff": diff,
        "se_diff": se,
        "n_sigma_observed": abs(diff) / se if se > 0 else nan,
        "resolvable": n_sigma * se if se > 0 else nan,
        # What the data exclude is the complement of the CI on the observed difference --
        # asymmetric because the observed difference is not zero.
        "ci95_lo": diff - 1.96 * se if se > 0 else nan,
        "ci95_hi": diff + 1.96 * se if se > 0 else nan,
        # 90% power at 2 sigma: (1.96 + 1.2816) * se.
        "resolvable_90power": 3.2416 * se if se > 0 else nan,
    }


def find_spectra(psr: dict) -> dict:
    """Per-pulsar $\\alpha$ and MSP/normal class for pulsars with a tabulated S400 and S1400.

    ``psr`` carries ``p0`` (s), ``s400`` and ``s1400`` (mJy). Returns ``alpha``, ``period_s``,
    ``s400``/``s1400``, and ``is_msp`` for the subset with finite, positive fluxes at both
    frequencies (plus ``name`` when the input carries one). ``n_rejected_nonpositive_flux``
    counts rows with both fluxes tabulated but one of them non-positive -- e.g. J0540-6919,
    whose S400 is recorded as 0.0 -- which the positivity guard silently dropped before this
    count made the cut visible.
    """
    s400 = np.asarray(psr["s400"], float)
    s1400 = np.asarray(psr["s1400"], float)
    p0 = np.asarray(psr["p0"], float)
    both = np.isfinite(s400) & np.isfinite(s1400) & np.isfinite(p0)
    good = both & (s400 > 0) & (s1400 > 0)
    alpha, _ = pulsar_alpha(s400[good], s1400[good])
    out = {
        "alpha": alpha,
        "period_s": p0[good],
        "s400": s400[good],
        "s1400": s1400[good],
        "is_msp": is_millisecond(p0[good]),
        "n_rejected_nonpositive_flux": int(np.sum(both & ~((s400 > 0) & (s1400 > 0)))),
    }
    if psr.get("name") is not None:
        out["name"] = np.asarray(psr["name"], dtype=object)[good]
    return out


def completeness_limited(
    alpha: np.ndarray,
    s400: np.ndarray,
    *,
    s1400_floor_mjy: float,
    alpha_min: float = -3.0,
) -> dict:
    """The joint-detection sample above the flux where the S1400 requirement cannot lose sources.

    Requiring a tabulated flux at *both* frequencies truncates asymmetrically: near the S400 end
    of the sample a steep source has already dropped below the effective S1400 floor, so the
    joint sample's mean index is biased flat. Above
    ``s400 >= s1400_floor * (1400/400)**(-alpha_min)`` no source with ``alpha >= alpha_min``
    can have been lost to the S1400 side, so the subsample's mean is free of that truncation
    (this is ``dr20radio``'s ``alpha_complete_limit`` applied to the identical geometry).

    Note the mirror check -- raising the *S1400* floor -- is **not** a clean bias estimate:
    it conditions on the numerator of alpha and flattens the mean mechanically.
    """
    a = np.asarray(alpha, float)
    s = np.asarray(s400, float)
    ratio = NU_GHZ["s1400"] / NU_GHZ["s400"]
    s400_cut = float(s1400_floor_mjy * ratio ** (-alpha_min))
    m = np.isfinite(a) & (s >= s400_cut)
    d = spectral_distribution(a[m])
    return {
        "s1400_floor_mjy": float(s1400_floor_mjy),
        "alpha_min": float(alpha_min),
        "s400_cut_mjy": round(s400_cut, 2),
        "n": d["n"],
        "mean": round(d["mean"], 3) if d["n"] else None,
        "median": round(d["median"], 3) if d["n"] else None,
    }


def msp_cut_sweep(
    alpha: np.ndarray,
    period_s: np.ndarray,
    *,
    cuts_s: tuple[float, ...] = (0.01, 0.02, 0.03, 0.05, 0.1),
) -> dict:
    """The MSP-normal comparison as the one classification constant sweeps.

    The 30 ms boundary is the only free choice the null depends on; unlike a vacuous variant,
    moving it genuinely reassigns membership in both arms. Keys are the cuts in seconds.
    """
    a = np.asarray(alpha, float)
    p = np.asarray(period_s, float)
    out = {}
    for cut in cuts_s:
        msp = is_millisecond(p, p_max=cut)
        c = compare_subsamples(a[msp], a[~msp])
        out[str(cut)] = {
            "n_msp": c["n_a"],
            "diff": round(c["diff"], 3) if np.isfinite(c["diff"]) else None,
            "two_se": (round(2.0 * c["se_diff"], 3) if np.isfinite(c["se_diff"]) else None),
        }
    return out


def permutation_pvalue(
    alpha: np.ndarray,
    is_msp: np.ndarray,
    *,
    n_perm: int = 20000,
    seed: int = 0,
) -> float:
    """Two-sided permutation p-value for the MSP-normal mean difference.

    At n = 43 with a skewed alpha distribution the Welch z is optimistic about its own
    calibration; the permutation test needs no distributional assumption. Uses the
    (k+1)/(B+1) convention.
    """
    a = np.asarray(alpha, float)
    m = np.asarray(is_msp, bool)
    ok = np.isfinite(a)
    a, m = a[ok], m[ok]
    n_msp = int(m.sum())
    if n_msp < 2 or (a.size - n_msp) < 2:
        return float("nan")
    obs = abs(float(np.mean(a[m]) - np.mean(a[~m])))
    rng = np.random.default_rng(seed)
    k = 0
    for _ in range(n_perm):
        perm = rng.permutation(a.size) < n_msp
        d = abs(float(np.mean(a[perm]) - np.mean(a[~perm])))
        if d >= obs:
            k += 1
    return float((k + 1) / (n_perm + 1))


def synthetic_field(
    n_sources: int = 1500,
    *,
    msp_fraction: float = 0.1,
    mean_alpha: float = -1.8,
    rel_err: float = 0.15,
    seed: int = 0,
) -> dict[str, np.ndarray]:
    """Synthetic pulsar population with steep spectra (offline fixture).

    Periods are bimodal (normal $\\sim$0.5 s, MSPs $\\sim$4 ms); spectral indices are drawn around
    ``mean_alpha`` (MSPs slightly flatter), and S400/S1400 follow $\\nu^\\alpha$ with scatter, so the
    recovered mean index returns ``mean_alpha``. Returns ``p0``/``s400``/``s1400`` and ``alpha_true``.
    """
    rng = np.random.default_rng(seed)
    is_msp = rng.random(n_sources) < msp_fraction
    p0 = np.where(
        is_msp,
        10.0 ** rng.uniform(-2.6, -2.0, n_sources),
        10.0 ** rng.uniform(-1.0, 0.5, n_sources),
    )
    alpha = rng.normal(mean_alpha, 0.5, n_sources)
    alpha[is_msp] = rng.normal(mean_alpha + 0.2, 0.5, int(is_msp.sum()))  # MSPs slightly flatter
    s400 = 10.0 ** rng.uniform(-0.5, 1.5, n_sources)  # ~0.3-30 mJy at 400 MHz
    s1400 = s400 * (NU_GHZ["s1400"] / NU_GHZ["s400"]) ** alpha
    s400 = s400 * rng.normal(1.0, rel_err, n_sources)
    s1400 = s1400 * rng.normal(1.0, rel_err, n_sources)
    # a minority are undetected at one frequency
    s1400[rng.random(n_sources) < 0.2] = np.nan
    return {
        "p0": p0,
        "s400": np.clip(s400, 1e-3, None),
        "s1400": s1400,
        "alpha_true": alpha,
        "is_msp": is_msp,
        "n_catalogue": np.asarray([n_sources], float),
    }


def fetch_atnf() -> dict:  # pragma: no cover - network
    """Fetch the ATNF Pulsar Catalogue (VizieR ``B/psr``): PSRJ, P0 (s), S400, S1400 (mJy).

    ``B/psr`` is a frozen VizieR snapshot (published 2017-07-18), behind the live ATNF
    ``psrcat``; the vintage is returned (and committed) so a reader querying psrcat directly
    knows why the counts differ. An absent S400/S1400 is *no published measurement at that
    frequency* -- overwhelmingly because the source's surveys did not observe there -- not a
    non-detection.
    """
    from datetime import datetime, timezone

    import numpy as _np
    from astroquery.vizier import Vizier

    v = Vizier(columns=["PSRJ", "P0", "S400", "S1400"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs("B/psr/psr")[0]
    return {
        "p0": _np.asarray(t["P0"], float),
        "s400": _np.asarray(t["S400"], float),
        "s1400": _np.asarray(t["S1400"], float),
        "name": _np.asarray([str(x).strip() for x in t["PSRJ"]], dtype=object),
        # The parent count, so the selection has a denominator. Without it the paper said
        # "over the whole catalogue" for what is really 473 of 2536, and this slice and `ppdot`
        # disagreed about the size of the same table (2536 vs a remembered "~3500").
        "n_catalogue": _np.asarray([len(t)], float),
        "catalogue_version": "B/psr (VizieR snapshot, published 2017-07-18)",
        "fetched_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Full slice: load (or fetch) the ATNF catalogue, compute pulsar spectral indices, write artifacts."""
    from pathlib import Path

    if offline:
        psr = synthetic_field()
        source = "synthetic"
    else:  # pragma: no cover - network
        psr = fetch_atnf()
        source = "ATNF Pulsar Catalogue (B/psr)"

    res = find_spectra(psr)
    alpha, is_msp = res["alpha"], res["is_msp"]
    dist = spectral_distribution(alpha)
    msp = spectral_distribution(alpha[is_msp])
    normal = spectral_distribution(alpha[~is_msp])
    # The comparison's sensitivity, not just its two means. Storing only the means made the
    # paper's "indistinguishable" unquantifiable from the committed evidence.
    cmp = compare_subsamples(alpha[is_msp], alpha[~is_msp])
    n_cat = psr.get("n_catalogue")
    # Single-band coverage: the binding constraint on the joint sample is which surveys
    # observed at 400 MHz, not a flux limit -- these two counts are the evidence.
    s400_all = np.asarray(psr["s400"], float)
    s1400_all = np.asarray(psr["s1400"], float)
    p0_all = np.asarray(psr["p0"], float)
    valid_p = np.isfinite(p0_all) & (p0_all > 0)
    metrics = {
        "source": source,
        "n": dist["n"],
        "n_catalogue": int(n_cat[0]) if n_cat is not None and len(n_cat) else None,
        "catalogue_version": psr.get("catalogue_version"),
        "fetched_utc": psr.get("fetched_utc"),
        "n_s400_tabulated": int(np.sum(np.isfinite(s400_all))),
        "n_s1400_tabulated": int(np.sum(np.isfinite(s1400_all))),
        "n_rejected_nonpositive_flux": res["n_rejected_nonpositive_flux"],
        "n_msp_catalogue": int(np.sum(valid_p & is_millisecond(p0_all))),
        "n_normal_catalogue": int(np.sum(valid_p & ~is_millisecond(p0_all))),
        "mean_alpha": round(dist["mean"], 2),
        "mean_alpha_se": round(dist["se_mean"], 4),
        "median_alpha": round(dist["median"], 2),
        "std_alpha": round(dist["std"], 2),
        "n_msp": msp["n"],
        "n_normal": normal["n"],
        "mean_alpha_msp": round(msp["mean"], 2) if msp["n"] else 0.0,
        "mean_alpha_normal": round(normal["mean"], 2) if normal["n"] else 0.0,
        "std_alpha_msp": round(cmp["std_a"], 2),
        "std_alpha_normal": round(cmp["std_b"], 2),
        "alpha_diff": round(cmp["diff"], 3),
        "alpha_diff_se": round(cmp["se_diff"], 3),
        "alpha_diff_sigma": round(cmp["n_sigma_observed"], 2),
        "alpha_resolvable_2sigma": round(cmp["resolvable"], 2),
        "alpha_diff_ci95_lo": round(cmp["ci95_lo"], 3),
        "alpha_diff_ci95_hi": round(cmp["ci95_hi"], 3),
        "alpha_resolvable_90power": round(cmp["resolvable_90power"], 2),
        "alpha_median_diff": round(float(np.median(alpha[is_msp]) - np.median(alpha[~is_msp])), 3),
        "permutation_p_two_sided": round(permutation_pvalue(alpha, is_msp), 3),
    }
    # The completeness-limited subsample: the one number that sizes the joint-detection bias
    # the Discussion names. The S1400 floor is the sample's own 1st percentile.
    floor = float(np.percentile(res["s1400"], 1.0))
    metrics["completeness"] = completeness_limited(
        alpha, res["s400"], s1400_floor_mjy=round(floor, 3), alpha_min=-3.0
    )
    metrics["msp_cut_sweep_s"] = msp_cut_sweep(alpha, res["period_s"])

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "pulsarspec_metrics.json")
    if res.get("name") is not None:
        _write_sources(res, op / "results" / "pulsarspec_sources.csv")
    _figure(res, op / "papers" / "pulsarspec" / "figures")
    _write_macros(metrics, op / "papers" / "pulsarspec" / "generated" / "macros.tex")
    return metrics


def _write_sources(res: dict, path) -> None:
    """The per-pulsar table behind the summary scalars (real leg only: needs names)."""
    import csv
    from pathlib import Path

    pt = Path(path)
    pt.parent.mkdir(parents=True, exist_ok=True)
    with pt.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["psrj", "p0_s", "s400_mjy", "s1400_mjy", "alpha", "is_msp"])
        for i in range(len(res["alpha"])):
            w.writerow(
                [
                    res["name"][i],
                    f"{res['period_s'][i]:.10g}",
                    f"{res['s400'][i]:.6g}",
                    f"{res['s1400'][i]:.6g}",
                    f"{res['alpha'][i]:.4f}",
                    int(res["is_msp"][i]),
                ]
            )


def _figure(res: dict, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    alpha = np.asarray(res["alpha"], float)
    is_msp = np.asarray(res["is_msp"], bool)
    p = np.asarray(res["period_s"], float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.6))
    bins = np.linspace(-4, 1, 30)
    ax1.hist(alpha[~is_msp], bins=bins, color="0.6", label="normal", alpha=0.8)
    ax1.hist(alpha[is_msp], bins=bins, color="r", label="MSP", alpha=0.7)
    ax1.axvline(float(np.median(alpha)), color="k", ls="--", lw=0.8)
    ax1.set(xlabel=r"spectral index $\alpha^{1400}_{400}$", ylabel="number", title="Pulsar spectra")
    ax1.legend(fontsize=8)
    ax2.scatter(p[~is_msp], alpha[~is_msp], s=5, color="0.6")
    ax2.scatter(p[is_msp], alpha[is_msp], s=8, color="r")
    ax2.set(xscale="log", xlabel="period (s)", ylabel=r"$\alpha$", title=r"$\alpha$ vs period")
    fig.tight_layout()
    fig.savefig(out / "spectra.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        v = m.get(key)
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "--"
        return str(v)

    lines = [
        "% Auto-generated by jansky_research.pulsarspec._write_macros — do not edit by hand.",
        rf"\newcommand{{\psrSource}}{{{m['source']}}}",
        rf"\newcommand{{\psrN}}{{{m['n']}}}",
        rf"\newcommand{{\psrMeanAlpha}}{{{m['mean_alpha']}}}",
        rf"\newcommand{{\psrMedianAlpha}}{{{m['median_alpha']}}}",
        rf"\newcommand{{\psrStdAlpha}}{{{m['std_alpha']}}}",
        rf"\newcommand{{\psrNmsp}}{{{m['n_msp']}}}",
        rf"\newcommand{{\psrMeanAlphaMsp}}{{{m['mean_alpha_msp']}}}",
        rf"\newcommand{{\psrMeanAlphaNormal}}{{{m['mean_alpha_normal']}}}",
        rf"\newcommand{{\psrNnormal}}{{{_fmt('n_normal')}}}",
        rf"\newcommand{{\psrNcatalogue}}{{{_fmt('n_catalogue')}}}",
        rf"\newcommand{{\psrStdAlphaMsp}}{{{_fmt('std_alpha_msp')}}}",
        rf"\newcommand{{\psrStdAlphaNormal}}{{{_fmt('std_alpha_normal')}}}",
        # The four numbers that turn "indistinguishable" from an impression into a limit.
        rf"\newcommand{{\psrAlphaDiff}}{{{_fmt('alpha_diff')}}}",
        rf"\newcommand{{\psrAlphaDiffSE}}{{{_fmt('alpha_diff_se')}}}",
        rf"\newcommand{{\psrAlphaDiffSigma}}{{{_fmt('alpha_diff_sigma')}}}",
        rf"\newcommand{{\psrAlphaResolvable}}{{{_fmt('alpha_resolvable_2sigma')}}}",
        # The exclusion is the CI's complement; `resolvable` is a 50%-power threshold.
        rf"\newcommand{{\psrCIlo}}{{{_fmt('alpha_diff_ci95_lo')}}}",
        rf"\newcommand{{\psrCIhi}}{{{_fmt('alpha_diff_ci95_hi')}}}",
        rf"\newcommand{{\psrNinetyPower}}{{{_fmt('alpha_resolvable_90power')}}}",
        rf"\newcommand{{\psrMeanAlphaSE}}{{{_fmt('mean_alpha_se')}}}",
        rf"\newcommand{{\psrMedianDiff}}{{{_fmt('alpha_median_diff')}}}",
        rf"\newcommand{{\psrPermP}}{{{_fmt('permutation_p_two_sided')}}}",
        # Coverage and selection: the joint sample's denominators.
        rf"\newcommand{{\psrNSfourHundred}}{{{_fmt('n_s400_tabulated')}}}",
        rf"\newcommand{{\psrNSfourteenHundred}}{{{_fmt('n_s1400_tabulated')}}}",
        rf"\newcommand{{\psrNmspCat}}{{{_fmt('n_msp_catalogue')}}}",
        rf"\newcommand{{\psrNnormalCat}}{{{_fmt('n_normal_catalogue')}}}",
        rf"\newcommand{{\psrNrejected}}{{{_fmt('n_rejected_nonpositive_flux')}}}",
    ]
    comp = m.get("completeness") or {}

    def _cfmt(key: str) -> str:
        v = comp.get(key)
        return "--" if v is None else str(v)

    lines += [
        rf"\newcommand{{\psrCompleteN}}{{{_cfmt('n')}}}",
        rf"\newcommand{{\psrCompleteMean}}{{{_cfmt('mean')}}}",
        rf"\newcommand{{\psrCompleteMedian}}{{{_cfmt('median')}}}",
        rf"\newcommand{{\psrCompleteCut}}{{{_cfmt('s400_cut_mjy')}}}",
        rf"\newcommand{{\psrCompleteFloor}}{{{_cfmt('s1400_floor_mjy')}}}",
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

    p = argparse.ArgumentParser(description="Pulsar radio spectral indices (ATNF catalogue).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
