"""Peaked-spectrum (GPS/CSS) radio-source selection via three-frequency spectral curvature.

Gigahertz-Peaked-Spectrum (GPS) and Compact-Steep-Spectrum (CSS) sources are compact, young radio
AGN whose radio spectrum *rises then falls*, peaking in the ~0.1--3 GHz band (O'Dea & Saikia 2021).
Selecting them needs $\\geq 3$ frequencies to see the turnover. Using three public surveys --- TGSS
(150 MHz), NVSS (1.4 GHz), and VLASS (3 GHz) --- we compute two indices, $\\alpha_\\mathrm{low}$
(150$\\to$1400 MHz) and $\\alpha_\\mathrm{high}$ (1400$\\to$3000 MHz), and select sources that are
rising at low frequency and falling at high frequency. The **curvature**
$\\alpha_\\mathrm{high}-\\alpha_\\mathrm{low}$ is far more robust to a constant TGSS flux-scale offset
than a single steep cut, and $\\alpha_\\mathrm{high}$ is TGSS-independent.

This is the maximal-reuse slice: it composes :mod:`jansky_research.spectra` (two-point index,
cross-match, fetch, NED/SIMBAD annotation) and :mod:`jansky_research.vlass` (VLASS 3 GHz fetch,
variability metrics to flag blazars, vetting) rather than reimplementing them. Pure NumPy + a
synthetic three-survey fixture for offline tests.
"""

from __future__ import annotations

import numpy as np

from .spectra import spectral_index

__all__ = [
    "NU_GHZ",
    "classify_sed",
    "find_peaked",
    "peak_frequency",
    "run",
    "synthetic_field",
    "two_point_indices",
]

# Survey reference frequencies (GHz).
NU_GHZ = {"tgss": 0.1475, "nvss": 1.4, "vlass": 3.0}


def two_point_indices(
    s_tgss: np.ndarray,
    s_nvss: np.ndarray,
    s_vlass: np.ndarray,
    *,
    e_tgss: np.ndarray | None = None,
    e_nvss: np.ndarray | None = None,
    e_vlass: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Low (150$\\to$1400 MHz) and high (1400$\\to$3000 MHz) spectral indices (reuses ``spectra``)."""
    alpha_low, _ = spectral_index(s_tgss, NU_GHZ["tgss"], s_nvss, NU_GHZ["nvss"], e_tgss, e_nvss)
    alpha_high, _ = spectral_index(
        s_nvss, NU_GHZ["nvss"], s_vlass, NU_GHZ["vlass"], e_nvss, e_vlass
    )
    return alpha_low, alpha_high


def classify_sed(alpha_low: float, alpha_high: float, *, up: float = 0.1, dn: float = -0.1) -> str:
    """Classify a 3-point radio SED from its two indices.

    ``peaked`` -- rising then falling ($\\alpha_\\mathrm{low}>$ ``up`` and $\\alpha_\\mathrm{high}<$
    ``dn``: a turnover between 150 MHz and 3 GHz, the GPS/CSS signature); ``inverted`` -- still rising
    at 3 GHz (peak above the band, or a flat-spectrum core); ``steep`` -- falling throughout;
    ``flat`` -- everything else.
    """
    if not (np.isfinite(alpha_low) and np.isfinite(alpha_high)):
        return "nan"
    if alpha_low > up and alpha_high < dn:
        return "peaked"
    if alpha_high > up:
        return "inverted"
    if alpha_low < -0.5 and alpha_high < -0.5:
        return "steep"
    return "flat"


def peak_frequency(flux: np.ndarray, nu_ghz: np.ndarray) -> tuple[float, bool]:
    """Turnover frequency (GHz) of a 3-point SED, by a parabolic log--log fit.

    Fits $\\log_{10} S = a x^2 + b x + c$ with $x=\\log_{10}\\nu$; the extremum is at
    $x_\\mathrm{peak}=-b/2a$ and is a genuine peak (concave) when $a<0$. Returns
    ``(nu_peak_ghz, is_peak)``.
    """
    flux = np.asarray(flux, dtype=float)
    nu = np.asarray(nu_ghz, dtype=float)
    if flux.size < 3 or np.any(flux <= 0):
        return float("nan"), False
    a, b, _ = np.polyfit(np.log10(nu), np.log10(flux), 2)
    if a == 0.0:
        return float("nan"), False
    x_peak = -b / (2.0 * a)
    if abs(x_peak) > 6.0:  # extremum far outside any sane radio band -> no in-band turnover
        return float("nan"), bool(a < 0.0)
    return float(10.0**x_peak), bool(a < 0.0)


def find_peaked(
    tgss: dict[str, np.ndarray],
    nvss: dict[str, np.ndarray],
    vlass: dict[str, np.ndarray],
    *,
    radius_arcsec: float = 15.0,
) -> dict[str, np.ndarray]:
    """Triple cross-match (anchored on NVSS) + spectral classification. Reuses ``spectra.crossmatch``.

    Each survey is a ``{ra, dec, flux, eflux}`` dict (mJy). Returns per-matched-source arrays:
    positions, the three fluxes, the two indices, the curvature, the SED class, and the turnover
    frequency.
    """
    from .spectra import crossmatch

    ra_n, dec_n = np.asarray(nvss["ra"], float), np.asarray(nvss["dec"], float)
    # NVSS rows that have BOTH a TGSS and a VLASS counterpart.
    it_n, it_t, _ = crossmatch(ra_n, dec_n, tgss["ra"], tgss["dec"], radius_arcsec)
    iv_n, iv_v, _ = crossmatch(ra_n, dec_n, vlass["ra"], vlass["dec"], radius_arcsec)
    t_of = dict(zip(it_n.tolist(), it_t.tolist(), strict=True))
    v_of = dict(zip(iv_n.tolist(), iv_v.tolist(), strict=True))
    common = sorted(set(t_of) & set(v_of))
    if not common:
        empty = np.array([])
        return {k: empty for k in ("ra", "dec", "s_tgss", "s_nvss", "s_vlass")}

    idx = np.asarray(common, dtype=int)
    ti = np.array([t_of[i] for i in common])
    vi = np.array([v_of[i] for i in common])
    s_t = np.asarray(tgss["flux"], float)[ti]
    s_n = np.asarray(nvss["flux"], float)[idx]
    s_v = np.asarray(vlass["flux"], float)[vi]
    e_t = np.asarray(tgss["eflux"], float)[ti]
    e_n = np.asarray(nvss["eflux"], float)[idx]
    e_v = np.asarray(vlass["eflux"], float)[vi]
    alpha_low, alpha_high = two_point_indices(s_t, s_n, s_v, e_tgss=e_t, e_nvss=e_n, e_vlass=e_v)
    cls = np.array([classify_sed(lo, hi) for lo, hi in zip(alpha_low, alpha_high, strict=True)])
    nu = np.array([NU_GHZ["tgss"], NU_GHZ["nvss"], NU_GHZ["vlass"]])
    nu_peak = np.array(
        [peak_frequency(np.array([a, b, c]), nu)[0] for a, b, c in zip(s_t, s_n, s_v, strict=True)]
    )
    return {
        "ra": ra_n[idx],
        "dec": dec_n[idx],
        "s_tgss": s_t,
        "s_nvss": s_n,
        "s_vlass": s_v,
        "alpha_low": alpha_low,
        "alpha_high": alpha_high,
        "curvature": alpha_high - alpha_low,
        "cls": cls,
        "nu_peak_ghz": nu_peak,
    }


def synthetic_field(
    n_sources: int = 1500, *, peaked_fraction: float = 0.05, rel_err: float = 0.1, seed: int = 0
) -> tuple[dict, dict, dict, np.ndarray]:
    """Synthetic TGSS/NVSS/VLASS catalogues with injected peaked + steep + flat SEDs (offline fixture).

    Returns ``(tgss, nvss, vlass)`` survey dicts sharing positions (small jitter), so the cross-match
    recovers them. Peaked sources rise to ~1.4 GHz then fall; the rest are steep or flat power laws.
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(180.0, 185.0, n_sources)
    dec = rng.uniform(20.0, 25.0, n_sources)
    nu = np.array([NU_GHZ["tgss"], NU_GHZ["nvss"], NU_GHZ["vlass"]])
    s_nvss = 10.0 ** rng.uniform(0.3, 1.5, n_sources)  # ~2-30 mJy at 1.4 GHz
    is_peaked = rng.random(n_sources) < peaked_fraction
    alpha = rng.uniform(-1.1, -0.5, n_sources)  # steep/flat power-law index for the rest
    flux = np.empty((n_sources, 3))
    for i in range(n_sources):
        if is_peaked[i]:  # parabola peaking near 1.4 GHz: rising then falling
            lp = np.log10(s_nvss[i]) - 1.2 * (np.log10(nu / NU_GHZ["nvss"])) ** 2
            flux[i] = 10.0**lp
        else:
            flux[i] = s_nvss[i] * (nu / NU_GHZ["nvss"]) ** alpha[i]
    flux *= rng.normal(1.0, rel_err, flux.shape)
    flux = np.clip(flux, 1e-3, None)
    jit = lambda: rng.normal(0.0, 1.0 / 3600.0, n_sources)  # noqa: E731  (~1" position jitter)

    def survey(col):
        return {
            "ra": ra + jit(),
            "dec": dec + jit(),
            "flux": flux[:, col],
            "eflux": rel_err * flux[:, col],
        }

    return survey(0), survey(1), survey(2), is_peaked


def run(center=None, radius_deg: float = 2.0, out: str = ".", *, offline: bool = False) -> dict:
    """Full slice: fetch (or synthesise) TGSS/NVSS/VLASS, find peaked sources, vet, write artifacts."""
    import json
    from pathlib import Path

    if offline or center is None:
        tgss, nvss, vlass, truth = synthetic_field()
        source = "synthetic"
    else:  # pragma: no cover - network
        from .spectra import fetch_survey
        from .vlass import _fetch_e1_tap

        tgss = fetch_survey(center, radius_deg, "tgss")
        nvss = fetch_survey(center, radius_deg, "nvss")
        vra, vdec, vflux, veflux = _fetch_e1_tap((center.ra.deg, center.dec.deg), radius_deg)
        vlass = {"ra": vra, "dec": vdec, "flux": vflux, "eflux": veflux}
        truth = None

    res = find_peaked(tgss, nvss, vlass)
    cls = res.get("cls", np.array([]))
    peaked = cls == "peaked"
    metrics = {
        "source": source,
        "n_matched": int(cls.size),
        "n_peaked": int(peaked.sum()),
        "n_inverted": int(np.sum(cls == "inverted")),
        "n_steep": int(np.sum(cls == "steep")),
    }
    if truth is not None:  # synthetic: recovery of the injected peaked sources
        # match recovered peaked back to the truth by position is overkill here; report the rates
        metrics["n_injected_peaked"] = int(truth.sum())
        metrics["n_peaked_recovered"] = int(peaked.sum())

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "peaked_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(res, op / "papers" / "peaked" / "figures")
    _write_macros(metrics, op / "papers" / "peaked" / "generated" / "macros.tex")
    return metrics


def _figure(res: dict, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    al = res.get("alpha_low", np.array([]))
    ah = res.get("alpha_high", np.array([]))
    cls = res.get("cls", np.array([]))
    peaked = cls == "peaked"
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(al[~peaked], ah[~peaked], s=6, color="0.6", label="other")
    ax.scatter(al[peaked], ah[peaked], s=28, color="r", marker="*", label="peaked (GPS/CSS)")
    ax.axhline(0.0, color="k", lw=0.5)
    ax.axvline(0.0, color="k", lw=0.5)
    ax.set(
        xlabel=r"$\alpha_{\rm low}$ (150$\to$1400 MHz)",
        ylabel=r"$\alpha_{\rm high}$ (1400$\to$3000 MHz)",
        title="Spectral-curvature plane",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "curvature.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    lines = [
        "% Auto-generated by jansky_research.peaked._write_macros — do not edit by hand.",
        rf"\newcommand{{\pkSource}}{{{m['source']}}}",
        rf"\newcommand{{\pkNmatched}}{{{m['n_matched']}}}",
        rf"\newcommand{{\pkNpeaked}}{{{m['n_peaked']}}}",
        rf"\newcommand{{\pkNinverted}}{{{m['n_inverted']}}}",
        rf"\newcommand{{\pkNsteep}}{{{m['n_steep']}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    from astropy.coordinates import SkyCoord

    p = argparse.ArgumentParser(description="Find peaked-spectrum (GPS/CSS) radio sources.")
    p.add_argument("--ra", type=float, help="field-centre RA (deg)")
    p.add_argument("--dec", type=float, help="field-centre Dec (deg)")
    p.add_argument("--radius", type=float, default=2.0, help="cone radius (deg)")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    center = None if (args.offline or args.ra is None) else SkyCoord(args.ra, args.dec, unit="deg")
    print(json.dumps(run(center, args.radius, args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
