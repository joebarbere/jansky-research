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
    "fetch_atnf",
    "find_spectra",
    "is_millisecond",
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
    """Summarise an $\\alpha$ distribution: ``n``, ``mean``, ``median``, ``std`` (finite values only)."""
    a = np.asarray(alpha, float)
    a = a[np.isfinite(a)]
    n = int(a.size)
    return {
        "n": n,
        "mean": float(np.mean(a)) if n else float("nan"),
        "median": float(np.median(a)) if n else float("nan"),
        "std": float(np.std(a)) if n else float("nan"),
    }


def find_spectra(psr: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Per-pulsar $\\alpha$ and MSP/normal class for pulsars with both S400 and S1400.

    ``psr`` carries ``p0`` (s), ``s400`` and ``s1400`` (mJy). Returns ``alpha``, ``period_s``, and
    ``is_msp`` for the subset detected at both frequencies (finite, positive fluxes).
    """
    s400 = np.asarray(psr["s400"], float)
    s1400 = np.asarray(psr["s1400"], float)
    p0 = np.asarray(psr["p0"], float)
    good = np.isfinite(s400) & np.isfinite(s1400) & (s400 > 0) & (s1400 > 0) & np.isfinite(p0)
    alpha, _ = pulsar_alpha(s400[good], s1400[good])
    return {"alpha": alpha, "period_s": p0[good], "is_msp": is_millisecond(p0[good])}


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
    }


def fetch_atnf() -> dict[str, np.ndarray]:  # pragma: no cover - network
    """Fetch the ATNF Pulsar Catalogue (VizieR ``B/psr``): P0 (s), S400, S1400 (mJy)."""
    import numpy as _np
    from astroquery.vizier import Vizier

    v = Vizier(columns=["PSRJ", "P0", "S400", "S1400"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs("B/psr/psr")[0]
    return {
        "p0": _np.asarray(t["P0"], float),
        "s400": _np.asarray(t["S400"], float),
        "s1400": _np.asarray(t["S1400"], float),
    }


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Full slice: load (or fetch) the ATNF catalogue, compute pulsar spectral indices, write artifacts."""
    import json
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
    metrics = {
        "source": source,
        "n": dist["n"],
        "mean_alpha": round(dist["mean"], 2),
        "median_alpha": round(dist["median"], 2),
        "std_alpha": round(dist["std"], 2),
        "n_msp": msp["n"],
        "mean_alpha_msp": round(msp["mean"], 2) if msp["n"] else 0.0,
        "mean_alpha_normal": round(normal["mean"], 2) if normal["n"] else 0.0,
    }

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "pulsarspec_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(res, op / "papers" / "pulsarspec" / "figures")
    _write_macros(metrics, op / "papers" / "pulsarspec" / "generated" / "macros.tex")
    return metrics


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
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


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
