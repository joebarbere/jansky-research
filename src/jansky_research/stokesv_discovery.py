"""RACS Stokes-V discovery: forced target-list photometry + multi-epoch V variability.

The `stokesv` slice built the methods (leakage floor, forced photometry, proper-motion vetting)
and found the honest limit: **single-epoch V is variability-limited** --- a coherent emitter caught
between bursts shows nothing. This module turns that limitation into the measurement. GATE 0
(2026-07-02, recorded in ``plans/33``) found the epoch pairs CASDA actually serves: RACS-mid1
(MJD 59233) vs RACS-mid2 (MJD 60769) is a **4.2-yr same-band pair at ~1368 MHz** (plus commensal
WALLABY V epochs between), and RACS-low2 vs low3 a near-band 2-yr pair. RACS-low2 Paper VIII
(arXiv:2606.16182) published the blind low-band V catalogue only --- forced photometry on a curated
M-dwarf list (below the blind 5-sigma threshold) and any multi-epoch V comparison are unclaimed.

Selection reuses the `stokesv` leakage floor (off-axis Stokes-I->V leakage is the killer
systematic) per epoch; association reuses `proper_motion_confirm` --- with the twist that here the
catalogue positions must be **propagated to each epoch's MJD** before the cutout is even requested
(nearby M dwarfs move arcsec/yr; 4.2 yr x 1"/yr is two RACS-mid pixels). The offline fixture
injects steady emitters, single-epoch flares, and leakage-only contaminants into a synthetic epoch
pair; the recover-a-known is that selection finds the emitters and the variability step flags
exactly the flares.
"""

from __future__ import annotations

import numpy as np

from .stokesv import (
    fractional_circular_pol,
    leakage_floor,
    select_circular_pol,
)

__all__ = [
    "epoch_position",
    "select_epoch_candidates",
    "synthetic_epoch_pair",
    "two_epoch_variability",
    "run",
]

#: GATE-0-verified epoch anchors (MJD of the RACS-mid1 / mid2 observations at the probe field).
RACS_MID1_MJD = 59233.4
RACS_MID2_MJD = 60769.2


def epoch_position(
    ra: np.ndarray,
    dec: np.ndarray,
    pmra: np.ndarray,
    pmdec: np.ndarray,
    cat_epoch_yr: float,
    obs_mjd: float,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Propagate catalogue positions to a radio epoch with Gaia proper motions.

    ``pmra`` is $\mu_{\alpha*}=\mu_\alpha\cos\delta$ in mas/yr (the Gaia convention), so the RA
    *coordinate* moves by $\mu_{\alpha*}\,\Delta t/\cos\delta$. ``cat_epoch_yr`` is the catalogue
    reference epoch (Gaia DR3: 2016.0); ``obs_mjd`` the radio observation. Nearby M dwarfs move
    up to arcsec/yr, so forced photometry MUST use the propagated position at each epoch.
    """
    dt_yr = (obs_mjd - 51544.5) / 365.25 + 2000.0 - cat_epoch_yr
    ra = np.asarray(ra, float)
    dec = np.asarray(dec, float)
    cosd = np.cos(np.radians(dec))
    ra_out = ra + (np.asarray(pmra, float) * dt_yr / 1000.0 / 3600.0) / cosd
    dec_out = dec + np.asarray(pmdec, float) * dt_yr / 1000.0 / 3600.0
    return ra_out, dec_out


def select_epoch_candidates(
    i_flux: np.ndarray,
    v_flux: np.ndarray,
    e_i: np.ndarray,
    e_v: np.ndarray,
    *,
    v_snr_min: float = 5.0,
) -> tuple[np.ndarray, float]:
    """Per-epoch candidate selection: the `stokesv` leakage floor + V-SNR cut, packaged.

    Returns ``(selected_mask, floor)`` where ``floor`` is this epoch's 7x-median |V|/I leakage
    floor --- computed per epoch because the leakage environment differs between observations.
    """
    frac, _ = fractional_circular_pol(v_flux, i_flux)
    floor = leakage_floor(frac)
    sel, _ = select_circular_pol(
        i_flux, v_flux, e_i, e_v, leakage_threshold=floor, v_snr_min=v_snr_min
    )
    return sel, float(floor)


def two_epoch_variability(
    v1: np.ndarray,
    e1: np.ndarray,
    v2: np.ndarray,
    e2: np.ndarray,
    sel1: np.ndarray,
    sel2: np.ndarray,
) -> dict[str, np.ndarray]:
    r"""Pairwise V variability between two epochs.

    Returns per-target arrays: ``dv_sig`` --- the variability significance
    $|V_2-V_1|/\sqrt{\sigma_1^2+\sigma_2^2}$ (on the signed fluxes, so a handedness flip counts as
    variability); ``appeared``/``disappeared`` --- selected in exactly one epoch (the coherent-
    emitter signature the single-epoch slice could not see); ``variable`` --- selected in at least
    one epoch with ``dv_sig`` >= 4.
    """
    v1 = np.asarray(v1, float)
    v2 = np.asarray(v2, float)
    dv_sig = np.abs(v2 - v1) / np.sqrt(np.asarray(e1, float) ** 2 + np.asarray(e2, float) ** 2)
    sel1 = np.asarray(sel1, bool)
    sel2 = np.asarray(sel2, bool)
    appeared = ~sel1 & sel2
    disappeared = sel1 & ~sel2
    variable = (sel1 | sel2) & (dv_sig >= 4.0)
    return {
        "dv_sig": dv_sig,
        "appeared": appeared,
        "disappeared": disappeared,
        "variable": variable,
    }


def synthetic_epoch_pair(
    n_stars: int = 600,
    *,
    steady_fraction: float = 0.03,
    flare_fraction: float = 0.03,
    leakage_scale: float = 0.008,
    seed: int = 0,
) -> dict[str, np.ndarray]:
    """Two synthetic forced-photometry epochs of one target list (the offline fixture).

    Three populations: **steady** coherent emitters (deep |V|/I in both epochs), **flaring**
    emitters (deep V in exactly one epoch --- the population single-epoch selection undercounts),
    and leakage-only contaminants (half-normal |V|/I ~ ``leakage_scale`` each epoch,
    uncorrelated). Truth masks ``is_steady``/``is_flare`` drive the recover-a-known.
    """
    rng = np.random.default_rng(seed)
    i_flux = 10.0 ** rng.uniform(0.7, 2.5, n_stars)  # ~5-300 mJy, both epochs
    e_i = 0.05 * i_flux + 0.25
    e_v = np.full(n_stars, 0.25)

    r = rng.random(n_stars)
    is_steady = r < steady_fraction
    is_flare = (r >= steady_fraction) & (r < steady_fraction + flare_fraction)
    sign = rng.choice([-1.0, 1.0], n_stars)
    frac = rng.uniform(0.2, 0.8, n_stars)

    def epoch(flare_on: np.ndarray, eseed: int) -> np.ndarray:
        erng = np.random.default_rng(eseed)
        leak = np.abs(erng.normal(0.0, leakage_scale, n_stars))
        v = sign * leak * i_flux + erng.normal(0.0, e_v)
        on = is_steady | (is_flare & flare_on)
        v[on] = sign[on] * frac[on] * i_flux[on] + erng.normal(0.0, e_v[on])
        return v

    flare_in_2 = rng.random(n_stars) < 0.5  # each flare is on in exactly one epoch
    v1 = epoch(~flare_in_2, seed + 1)
    v2 = epoch(flare_in_2, seed + 2)
    return {
        "i_flux": i_flux,
        "e_i": e_i,
        "e_v": e_v,
        "v1": v1,
        "v2": v2,
        "is_steady": is_steady,
        "is_flare": is_flare,
    }


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline slice: epoch-pair selection + variability recover-a-known (real leg is next PR)."""
    import json
    from pathlib import Path

    s = synthetic_epoch_pair()
    sel1, floor1 = select_epoch_candidates(s["i_flux"], s["v1"], s["e_i"], s["e_v"])
    sel2, floor2 = select_epoch_candidates(s["i_flux"], s["v2"], s["e_i"], s["e_v"])
    var = two_epoch_variability(s["v1"], s["e_v"], s["v2"], s["e_v"], sel1, sel2)

    emitters = s["is_steady"] | s["is_flare"]
    union = sel1 | sel2
    completeness = float(union[emitters].mean()) if emitters.any() else float("nan")
    purity = float(emitters[union].mean()) if union.any() else float("nan")
    flagged_flare = var["appeared"] | var["disappeared"]
    var_completeness = float(flagged_flare[s["is_flare"]].mean()) if s["is_flare"].any() else 0.0
    var_purity = float(s["is_flare"][flagged_flare].mean()) if flagged_flare.any() else 0.0
    # the headline the single-epoch slice could not measure: how many emitters one epoch misses
    single_epoch_miss = float(1.0 - sel1[emitters].mean()) if emitters.any() else float("nan")

    metrics = {
        "source": "synthetic epoch pair",
        "n_targets": int(s["i_flux"].size),
        "floor_epoch1_pct": round(100 * floor1, 2),
        "floor_epoch2_pct": round(100 * floor2, 2),
        "completeness": round(completeness, 3),
        "purity": round(purity, 3),
        "var_completeness": round(var_completeness, 3),
        "var_purity": round(var_purity, 3),
        "single_epoch_miss_frac": round(single_epoch_miss, 3),
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "stokesv_discovery_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n"
    )
    _figure(s, sel1, sel2, var, (floor1, floor2), op / "papers" / "stokesv_discovery" / "figures")
    _write_macros(metrics, op / "papers" / "stokesv_discovery" / "generated" / "macros.tex")
    return metrics


def _figure(s, sel1, sel2, var, floors, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 4.0))

    frac1 = np.abs(s["v1"]) / s["i_flux"]
    emitters = s["is_steady"] | s["is_flare"]
    ax1.loglog(s["i_flux"][~emitters], frac1[~emitters], ".", ms=3, color="0.6", label="leakage")
    ax1.loglog(s["i_flux"][emitters], frac1[emitters], "o", ms=4, color="C3", label="true emitters")
    ax1.axhline(floors[0], color="C0", ls="--", lw=1, label="epoch-1 leakage floor")
    ax1.set(xlabel="I (mJy)", ylabel="|V|/I (epoch 1)", title="Per-epoch selection")
    ax1.legend(fontsize=8)

    ax2.plot(s["v1"][~emitters], s["v2"][~emitters], ".", ms=3, color="0.6")
    ax2.plot(
        s["v1"][s["is_steady"]], s["v2"][s["is_steady"]], "o", ms=4, color="C2", label="steady"
    )
    ax2.plot(s["v1"][s["is_flare"]], s["v2"][s["is_flare"]], "s", ms=4, color="C3", label="flare")
    flagged = var["appeared"] | var["disappeared"]
    ax2.plot(
        s["v1"][flagged], s["v2"][flagged], "x", ms=8, color="C1", label="flagged appear/disappear"
    )
    lim = max(np.abs(ax2.get_xlim()).max(), np.abs(ax2.get_ylim()).max())
    ax2.plot([-lim, lim], [-lim, lim], "-", color="0.8", lw=0.8, zorder=0)
    ax2.set(xlabel="V epoch 1 (mJy)", ylabel="V epoch 2 (mJy)", title="Two-epoch V plane")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "stokesv_discovery.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.stokesv_discovery._write_macros -- do not edit.",
        rf"\newcommand{{\svdSource}}{{{m['source']}}}",
        rf"\newcommand{{\svdNtargets}}{{{_fmt('n_targets')}}}",
        rf"\newcommand{{\svdFloorA}}{{{_fmt('floor_epoch1_pct')}}}",
        rf"\newcommand{{\svdFloorB}}{{{_fmt('floor_epoch2_pct')}}}",
        rf"\newcommand{{\svdComp}}{{{_fmt('completeness')}}}",
        rf"\newcommand{{\svdPur}}{{{_fmt('purity')}}}",
        rf"\newcommand{{\svdVarComp}}{{{_fmt('var_completeness')}}}",
        rf"\newcommand{{\svdVarPur}}{{{_fmt('var_purity')}}}",
        rf"\newcommand{{\svdMissFrac}}{{{_fmt('single_epoch_miss_frac')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(
        description="RACS Stokes-V discovery: epoch-pair tooling (offline)."
    )
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=True), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
