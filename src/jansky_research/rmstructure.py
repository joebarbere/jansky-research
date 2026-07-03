"""Galactic RM structure functions from SPICE-RACS --- the `rmsky` slice at 10x the data (plan 36).

SPICE-RACS DR2 (arXiv:2605.16917; ~2.5--3.4x10^5 RMs over 87.5% of sky, 6.7 deg^-2) was released
without a systematic Galactic structure-function analysis by latitude. This module supplies the
tooling: the second-order **RM structure function** :math:`\\mathrm{SF}(\\delta\\theta) =
\\langle[\\mathrm{RM}(\\hat n_1)-\\mathrm{RM}(\\hat n_2)]^2\\rangle` in angular-separation bins,
**noise-debiased** by subtracting :math:`\\langle\\sigma_1^2+\\sigma_2^2\\rangle` per bin (the
measurement-error term that otherwise floors the SF), with bootstrap errors --- computed per
Galactic-latitude bin so the turbulence amplitude and any coherence-scale break can be compared
across the disc--halo transition. The latitude profile / quadrant machinery reuses `rmsky`.

GATE 0 (2026-07-02): DR2's catalogue is public on the CSIRO DAP (collection csiro:64891,
`spice-racs.dr2.fits.gz`, 4.97 GB, no auth); the verified bounded first leg is DR1 on CASDA TAP ---
`AS110.spice_racs_dr1_corrected_cut_v02`, **24,758 rows** (live count) with `l, b, rm, rm_err,
snr_polint` columns. Offline, a synthetic Gaussian RM screen with a KNOWN coherence scale and a
known latitude enhancement drives the recover-a-known.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .rmsky import _ratio_bootstrap_se, enhancement_ratio

__all__ = [
    "latitude_ladder",
    "structure_function",
    "synthetic_rm_screen",
    "fetch_spice_racs_dr1",
    "run",
]

SPICE_DR1_TABLE = "AS110.spice_racs_dr1_corrected_cut_v02"
CASDA_TAP = "https://casda.csiro.au/casda_vo_tools/tap"
DR2_DAP = "https://data.csiro.au/collection/csiro:64891"


def structure_function(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    rm: np.ndarray,
    rm_err: np.ndarray,
    *,
    bins_deg: np.ndarray | None = None,
    n_boot: int = 100,
    max_pairs: int = 2_000_000,
    seed: int = 0,
) -> dict:
    r"""Noise-debiased second-order RM structure function with bootstrap errors.

    For every source pair within the separation bins:
    :math:`\mathrm{SF}(\delta\theta) = \langle(\mathrm{RM}_1-\mathrm{RM}_2)^2\rangle -
    \langle\sigma_1^2+\sigma_2^2\rangle` --- the second term removes the measurement-noise floor
    (Haverkorn et al. 2004 convention). Pairs are randomly subsampled to ``max_pairs`` (recorded)
    so the cost stays quadratic-safe; bootstrap resamples *sources* (not pairs) to respect the
    correlated pair structure. Returns bin centres, SF, its bootstrap SE, pair counts, and the
    subsample fraction.
    """
    rng = np.random.default_rng(seed)
    ra = np.radians(np.asarray(ra_deg, float))
    dec = np.radians(np.asarray(dec_deg, float))
    rm = np.asarray(rm, float)
    var = np.asarray(rm_err, float) ** 2
    n = rm.size
    if bins_deg is None:
        bins_deg = np.logspace(-1, 1.3, 12)  # 0.1 -- 20 deg
    bins = np.radians(np.asarray(bins_deg, float))

    # all pairs (i<j) when tractable; RANDOM pair draws for large n (triu at n~2.5e5 would
    # need ~3e10 index entries -- hundreds of GB). Random pairs are an unbiased SF estimator.
    n_pairs_all = n * (n - 1) // 2
    if n <= 3000:
        i_idx, j_idx = np.triu_indices(n, k=1)
        if n_pairs_all > max_pairs:
            keep = rng.choice(n_pairs_all, max_pairs, replace=False)
            i_idx, j_idx = i_idx[keep], j_idx[keep]
    else:
        n_draw = int(min(max_pairs, n_pairs_all))
        i_idx = rng.integers(0, n, n_draw)
        j_idx = rng.integers(0, n, n_draw)
        good = i_idx != j_idx
        i_idx, j_idx = i_idx[good], j_idx[good]

    # angular separation via the haversine formula (stable at small angles)
    sdlat = np.sin((dec[j_idx] - dec[i_idx]) / 2.0)
    sdlon = np.sin((ra[j_idx] - ra[i_idx]) / 2.0)
    h = sdlat**2 + np.cos(dec[i_idx]) * np.cos(dec[j_idx]) * sdlon**2
    sep = 2.0 * np.arcsin(np.sqrt(np.clip(h, 0.0, 1.0)))

    d2 = (rm[i_idx] - rm[j_idx]) ** 2
    nvar = var[i_idx] + var[j_idx]
    which = np.digitize(sep, bins) - 1
    nb = len(bins) - 1

    def sf_of(mask_sources: np.ndarray) -> np.ndarray:
        m = mask_sources[i_idx] & mask_sources[j_idx]
        out = np.full(nb, np.nan)
        for b in range(nb):
            sel = m & (which == b)
            if sel.sum() >= 20:
                out[b] = d2[sel].mean() - nvar[sel].mean()
        return out

    all_mask = np.ones(n, bool)
    sf = sf_of(all_mask)
    boots = np.full((n_boot, nb), np.nan)
    for k in range(n_boot):
        pick = rng.integers(0, n, n)
        mask = np.zeros(n, bool)
        mask[np.unique(pick)] = True  # source-level resample (unique-set approximation)
        boots[k] = sf_of(mask)
    se = np.nanstd(boots, axis=0)
    counts = np.array([(which == b).sum() for b in range(nb)])
    centres = np.degrees(np.sqrt(bins[:-1] * bins[1:]))
    return {
        "sep_deg": centres,
        "sf": sf,
        "sf_err": se,
        "n_pairs": counts,
        "pair_fraction": min(1.0, max_pairs / max(1, n_pairs_all)),
    }


def synthetic_rm_screen(
    n_sources: int = 1500,
    *,
    coherence_deg: float = 2.0,
    amp_high_b: float = 15.0,
    plane_boost: float = 5.0,
    noise: float = 2.0,
    seed: int = 0,
) -> dict:
    r"""A synthetic RM sky with a KNOWN coherence scale and a known plane enhancement.

    Sources are scattered over a patch; the RM field is a sum of Gaussian blobs of angular size
    ``coherence_deg`` (so the SF rises up to ~ the coherence scale and saturates at
    :math:`2\sigma_\mathrm{RM}^2` beyond it), with the RM amplitude boosted by ``plane_boost``
    at low latitude. Gaussian measurement noise ``noise`` (rad/m^2) is added and recorded in
    ``rm_err`` --- the debiasing target the SF must remove.
    """
    rng = np.random.default_rng(seed)
    ra = rng.uniform(0.0, 40.0, n_sources)
    dec = rng.uniform(-20.0, 20.0, n_sources)
    gb = dec  # patch geometry: treat dec as latitude for the fixture
    # Gaussian-blob random field with the requested coherence scale
    n_blob = 220
    bra = rng.uniform(-5.0, 45.0, n_blob)
    bdec = rng.uniform(-25.0, 25.0, n_blob)
    bamp = rng.normal(0.0, amp_high_b, n_blob)
    rm_true = np.zeros(n_sources)
    for k in range(n_blob):
        d2 = (ra - bra[k]) ** 2 * np.cos(np.radians(dec)) ** 2 + (dec - bdec[k]) ** 2
        rm_true += bamp[k] * np.exp(-0.5 * d2 / coherence_deg**2)
    rm_true *= 1.0 + (plane_boost - 1.0) * np.exp(-0.5 * (gb / 5.0) ** 2)
    rm_err = np.full(n_sources, noise)
    rm_obs = rm_true + rng.normal(0.0, noise, n_sources)
    return {
        "ra": ra,
        "dec": dec,
        "gal_b": gb,
        "gal_l": ra,
        "rm": rm_obs,
        "rm_err": rm_err,
        "coherence_deg": coherence_deg,
    }


def fetch_spice_racs_dr1(
    *, snr_min: float = 8.0, max_rows: int = 30000
) -> dict:  # pragma: no cover - network
    """Fetch the SPICE-RACS DR1 RM catalogue from CASDA TAP (verified 24,758 rows, no auth)."""
    from astroquery.utils.tap.core import TapPlus

    tap = TapPlus(url=CASDA_TAP)
    q = (
        f"SELECT TOP {int(max_rows)} ra, dec, l, b, rm, rm_err, snr_polint "
        f"FROM {SPICE_DR1_TABLE} WHERE snr_polint >= {snr_min} AND rm_err > 0"
    )
    t = tap.launch_job(q).get_results()
    return {
        "ra": np.asarray(t["ra"], float),
        "dec": np.asarray(t["dec"], float),
        "gal_l": np.asarray(t["l"], float),
        "gal_b": np.asarray(t["b"], float),
        "rm": np.asarray(t["rm"], float),
        "rm_err": np.asarray(t["rm_err"], float),
    }


DR2_LOCAL = Path("data/spice-racs.dr2.fits")  # gunzipped from the DAP .gz so astropy can memmap


def load_spice_racs_dr2(
    path: str | Path = DR2_LOCAL, *, snr_min: float = 8.0
) -> dict:  # pragma: no cover - needs the 5 GB DAP file
    """Load the public SPICE-RACS DR2 catalogue FITS (CSIRO DAP csiro:64891, no auth).

    Column names follow the DR1 convention (rm, rm_err, snr_polint, l, b); any variant casing is
    resolved by lookup. Applies the S/N cut and finite-error filter used throughout.
    """
    from astropy.io import fits
    from astropy.table import Table

    with fits.open(path, memmap=True) as hdul:
        t = Table(hdul[1].data)
    cols = {c.lower(): c for c in t.colnames}

    def col(*names):
        for nm in names:
            if nm in cols:
                return np.asarray(t[cols[nm]], float)
        raise KeyError(f"none of {names} in DR2 table (has: {sorted(cols)[:40]}...)")

    rm = col("rm")
    rm_err = col("rm_err", "e_rm", "rm_err_obs")
    snr = col("snr_polint", "snr_pi", "snr")
    gl = col("l", "gal_l", "glon")
    gb = col("b", "gal_b", "glat")
    ra = col("ra", "ra_deg")
    dec = col("dec", "dec_deg")
    m = (snr >= snr_min) & np.isfinite(rm) & (rm_err > 0)
    return {
        "ra": ra[m],
        "dec": dec[m],
        "gal_l": gl[m],
        "gal_b": gb[m],
        "rm": rm[m],
        "rm_err": rm_err[m],
    }


def latitude_ladder(
    s: dict,
    *,
    b_edges: tuple = (0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 90.0),
    max_pairs: int = 500_000,
    n_boot: int = 40,
) -> dict:
    """SF plateau (and its sqrt/2 = RM dispersion) per |b| bin --- the fluctuation-power profile.

    The two-bin disc--halo split showed a factor-~23 contrast; the ladder resolves HOW the RM
    fluctuation power falls with latitude. Each bin's plateau is the median of the SF's three
    largest-separation finite bins; per-bin source counts are reported (thin bins are honest
    NaNs). The intrinsic+extragalactic floor is latitude-independent, so the ladder's SHAPE is
    Galactic even though each absolute value is an upper bound.
    """
    ab = np.abs(np.asarray(s["gal_b"], float))
    out: dict[str, list] = {
        "b_lo": [],
        "b_hi": [],
        "n": [],
        "plateau": [],
        "plateau_err": [],
        "sigma_rm": [],
    }
    for lo, hi in zip(b_edges[:-1], b_edges[1:], strict=True):
        m = (ab >= lo) & (ab < hi)
        out["b_lo"].append(lo)
        out["b_hi"].append(hi)
        out["n"].append(int(m.sum()))
        if m.sum() < 200:
            out["plateau"].append(float("nan"))
            out["plateau_err"].append(float("nan"))
            out["sigma_rm"].append(float("nan"))
            continue
        sf = structure_function(
            s["ra"][m],
            s["dec"][m],
            s["rm"][m],
            s["rm_err"][m],
            max_pairs=max_pairs,
            n_boot=n_boot,
        )
        good = np.isfinite(sf["sf"])
        plat = float(np.nanmedian(sf["sf"][good][-3:])) if good.sum() >= 3 else float("nan")
        perr = float(np.nanmedian(sf["sf_err"][good][-3:])) if good.sum() >= 3 else float("nan")
        out["plateau"].append(plat)
        out["plateau_err"].append(perr)
        out["sigma_rm"].append(float(np.sqrt(plat / 2.0)) if plat > 0 else float("nan"))
    lad: dict[str, Any] = {k: np.asarray(v) for k, v in out.items()}
    # first-order floor subtraction: the highest-|b| bin estimates the latitude-independent
    # intrinsic+extragalactic floor; sigma_gal = sqrt(sigma^2 - floor^2) (NaN where <= floor)
    floor = lad["sigma_rm"][np.isfinite(lad["sigma_rm"])][-1]
    with np.errstate(invalid="ignore"):
        lad["sigma_gal"] = np.sqrt(np.clip(lad["sigma_rm"] ** 2 - floor**2, 0.0, None))
    lad["floor_sigma"] = float(floor)
    return lad


def _sf_break(sep_deg: np.ndarray, sf: np.ndarray) -> float:
    """Crude coherence-scale estimate: separation where the SF first reaches half its plateau."""
    good = np.isfinite(sf) & (sf > 0)
    if good.sum() < 4:
        return float("nan")
    plateau = np.nanmedian(sf[good][-3:])
    rising = np.where(good & (sf >= 0.5 * plateau))[0]
    return float(sep_deg[rising[0]]) if rising.size else float("nan")


def run(out: str = ".", *, offline: bool = True, dr2: bool = False) -> dict:
    """Offline: SF + latitude recover-a-known on the synthetic screen; real: SPICE-RACS DR1."""
    import json
    from pathlib import Path

    if offline:
        s = synthetic_rm_screen()
        source = "synthetic RM screen"
    elif dr2:  # pragma: no cover - needs the local DAP file
        s = load_spice_racs_dr2()
        s["coherence_deg"] = float("nan")
        source = "SPICE-RACS DR2 (CSIRO DAP csiro:64891)"
    else:  # pragma: no cover - network
        s = fetch_spice_racs_dr1()
        s["coherence_deg"] = float("nan")
        source = f"SPICE-RACS DR1 ({SPICE_DR1_TABLE})"

    lo = np.abs(s["gal_b"]) < 10.0
    hi = np.abs(s["gal_b"]) > 10.0
    sf_lo = structure_function(s["ra"][lo], s["dec"][lo], s["rm"][lo], s["rm_err"][lo])
    sf_hi = structure_function(s["ra"][hi], s["dec"][hi], s["rm"][hi], s["rm_err"][hi])
    pole = 15.0 if offline else 60.0  # the synthetic patch spans only ±20 deg
    ratio = enhancement_ratio(s["rm"], s["gal_b"], pole_deg=pole)
    ratio_se = _ratio_bootstrap_se(s["rm"], s["gal_b"], pole_deg=pole)
    break_lo = _sf_break(sf_lo["sep_deg"], sf_lo["sf"])
    break_hi = _sf_break(sf_hi["sep_deg"], sf_hi["sf"])

    ladder = latitude_ladder(s) if not offline else None  # pragma: no cover - big data only
    metrics = {
        "source": source,
        "n_sources": int(s["rm"].size),
        "enhancement_ratio": round(float(ratio), 2),
        "enhancement_ratio_se": round(float(ratio_se), 2),
        "sf_plateau_low_b": round(float(np.nanmedian(sf_lo["sf"][-3:])), 1),
        "sf_plateau_high_b": round(float(np.nanmedian(sf_hi["sf"][-3:])), 1),
        "sf_break_low_b_deg": round(break_lo, 2) if np.isfinite(break_lo) else None,
        "sf_break_high_b_deg": round(break_hi, 2) if np.isfinite(break_hi) else None,
        "true_coherence_deg": s.get("coherence_deg"),
    }
    if ladder is not None:  # pragma: no cover - big data only
        fin = np.isfinite(ladder["sigma_rm"])
        metrics.update(
            {
                "ladder_bins": [
                    {
                        "b": f"{ladder['b_lo'][i]:.0f}-{ladder['b_hi'][i]:.0f}",
                        "n": int(ladder["n"][i]),
                        "sigma_rm": round(float(ladder["sigma_rm"][i]), 1),
                        "sigma_gal": round(float(ladder["sigma_gal"][i]), 1),
                    }
                    for i in range(len(ladder["n"]))
                    if fin[i]
                ],
                "sigma_rm_plane": round(float(ladder["sigma_rm"][fin][0]), 1),
                "sigma_rm_pole": round(float(ladder["sigma_rm"][fin][-1]), 1),
                "ladder_floor_sigma": round(ladder["floor_sigma"], 1),
            }
        )
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "rmstructure_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(s, sf_lo, sf_hi, op / "papers" / "rmstructure" / "figures")
    _write_macros(metrics, op / "papers" / "rmstructure" / "generated" / "macros.tex")
    return metrics


def _figure(s, sf_lo, sf_hi, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.9))
    sc = ax1.scatter(s["ra"], s["gal_b"], c=s["rm"], s=4, cmap="RdBu_r", vmin=-40, vmax=40)
    fig.colorbar(sc, ax=ax1, label="RM (rad m$^{-2}$)")
    ax1.set(xlabel="lon (deg)", ylabel="lat (deg)", title="RM sky")
    for sf, lab, c in ((sf_lo, "|b| < 10°", "C3"), (sf_hi, "|b| > 10°", "C0")):
        g = np.isfinite(sf["sf"])
        ax2.errorbar(
            sf["sep_deg"][g],
            sf["sf"][g],
            yerr=sf["sf_err"][g],
            fmt="o-",
            ms=3,
            color=c,
            lw=1,
            label=lab,
        )
    ax2.set(
        xscale="log",
        yscale="log",
        xlabel=r"separation $\delta\theta$ (deg)",
        ylabel=r"SF(RM) (rad$^2$ m$^{-4}$)",
        title="Noise-debiased structure function",
    )
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "rmstructure.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.rmstructure._write_macros -- do not edit.",
        rf"\newcommand{{\rmsSource}}{{{m['source']}}}",
        rf"\newcommand{{\rmsN}}{{{_fmt('n_sources')}}}",
        rf"\newcommand{{\rmsRatio}}{{{_fmt('enhancement_ratio')}}}",
        rf"\newcommand{{\rmsRatioSe}}{{{_fmt('enhancement_ratio_se')}}}",
        rf"\newcommand{{\rmsPlatLo}}{{{_fmt('sf_plateau_low_b')}}}",
        rf"\newcommand{{\rmsPlatHi}}{{{_fmt('sf_plateau_high_b')}}}",
        rf"\newcommand{{\rmsBreakLo}}{{{_fmt('sf_break_low_b_deg')}}}",
        rf"\newcommand{{\rmsBreakHi}}{{{_fmt('sf_break_high_b_deg')}}}",
        rf"\newcommand{{\rmsSigPlane}}{{{_fmt('sigma_rm_plane')}}}",
        rf"\newcommand{{\rmsSigPole}}{{{_fmt('sigma_rm_pole')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="SPICE-RACS RM structure functions by latitude.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--dr2", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline, dr2=args.dr2), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
