#!/usr/bin/env python
"""vlasspm note, referee round 2: the two extra runs and the vector UV Ceti comparison.

    uv run python scripts/vlasspm_referee2.py --out .      # -> results/vlasspm_referee2.json

1. Fine-bin completeness at 1.5 mJy with the cut (same bins as the 3 mJy run of round 1).
2. Parallax-aware injections: 3 mJy, cut, fine bins, at 4, 8 and 16 pc and without parallax
   (same seed, so the runs differ only by the parallax). The domain is the smallest tested
   distance at which the worst fine-bin completeness over 1.1-5"/yr is within 5% (relative) of
   the no-parallax run; the limit is quoted from the completeness AT that distance.
3. UV Ceti: the radio proper-motion VECTOR against the SIMBAD/UCAC4 system motion and the two
   Gaia components (chi^2 with the fit covariance), and the orbital scale |UV - BL|.

Reads the round-1 evidence (results/vlasspm_referee1.json) and the work directory; no network.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import vlasspm_real as R  # noqa: E402
from vlasspm_referee1 import FINE_BINS, _cached  # noqa: E402

from jansky_research import vlasspm as v  # noqa: E402

DISTANCES_PC = (None, 4.0, 8.0, 16.0)
RATE_LO = 1.1  # the limit's rate range starts here: every fine bin above is complete to ~97%
DOMAIN_LOSS = 0.05  # relative completeness loss that ends the distance domain
KM_S_PER_AU_YR = 4.740470


def _fine(cats, model, size_noise, thr, rule, n, flux, distance_pc, seed=97) -> dict:
    c = v.completeness(
        *cats,
        n=n,
        flux_mjy=flux,
        seed=seed,
        realistic_frac=1.0,
        model=model,
        size_noise=size_noise,
        compact_max=thr,
        compact_rule=rule,
        bins=np.asarray(FINE_BINS),
        distance_pc=distance_pc,
    )
    return {k: c[k] for k in ("bin_edges", "per_bin", "n_per_bin", "overall", "n_injected")} | {
        "flux_mjy": flux,
        "distance_pc": distance_pc,
    }


def worst_above(fine: dict, lo: float = RATE_LO) -> float:
    edges = fine["bin_edges"]
    return float(min(p for p, e in zip(fine["per_bin"], edges[:-1], strict=True) if e >= lo - 1e-9))


def uv_vectors(ref1: dict) -> dict:
    uv = ref1["uvcet"]
    fit = uv["fits"]["fixed_gaia_parallax"]
    mu = np.array([fit["mu_ra"], fit["mu_dec"]])
    cov = np.asarray(fit["mu_cov"])
    g = uv["gaia"]
    ref = {
        "system (SIMBAD/UCAC4)": np.array([g["system"]["pmra_mas_yr"], g["system"]["pmdec_mas_yr"]])
        / 1000,
        "UV Cet (Gaia DR3)": np.array(
            [g["components"]["UV Cet"]["pmra"], g["components"]["UV Cet"]["pmdec"]]
        )
        / 1000,
        "BL Cet (Gaia DR3)": np.array(
            [g["components"]["BL Cet"]["pmra"], g["components"]["BL Cet"]["pmdec"]]
        )
        / 1000,
    }
    out = {"radio_mu": mu.tolist(), "radio_mu_cov": cov.tolist(), "comparisons": {}}
    inv = np.linalg.inv(cov)
    for name, r in ref.items():
        d = mu - r
        out["comparisons"][name] = {
            "mu": r.tolist(),
            "offset": d.tolist(),
            "offset_abs": float(np.hypot(*d)),
            "chi2_2dof": float(d @ inv @ d),
        }
    orb = ref["UV Cet (Gaia DR3)"] - ref["BL Cet (Gaia DR3)"]
    out["uv_minus_bl_abs"] = float(np.hypot(*orb))
    out["uv_minus_bl_parallax_arcsec"] = (
        g["components"]["UV Cet"]["parallax"] - g["components"]["BL Cet"]["parallax"]
    ) / 1000
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=str(REPO / "data/vlass/work"))
    ap.add_argument("--out", default=str(REPO))
    ap.add_argument("--n-inject", type=int, default=20000)
    args = ap.parse_args()
    work, out = Path(args.work), Path(args.out)
    floors = json.loads((work / "floors.json").read_text())
    banded = json.loads((work / "floors_banded.json").read_text())
    cats = R.catalogs(work, floors["floor_arcsec"], banded)
    calib = json.loads((work / "calibration.json").read_text())
    model = v.ErrorModel(calib["model"]["k_struct"], calib["model"]["q_beam"])
    r5 = json.loads((work / "run5_calibration.json").read_text())
    stars = json.loads((work / "run5_stars.json").read_text())
    size_noise = v.SizeNoiseModel(np.asarray(stars["samples"]))
    thr, rule = r5["threshold"], r5["rule"]
    ref1 = json.loads((out / "results" / "vlasspm_referee1.json").read_text())
    metrics = json.loads((out / "results" / "vlasspm_metrics.json").read_text())
    area = metrics["real_run5_compactness"]["area_deg2"]

    R.log("referee 2: fine bins at 1.5 mJy ...")
    f15 = _cached(
        work / "ref2_fine_1p5mJy.json",
        lambda: _fine(cats, model, size_noise, thr, rule, args.n_inject, 1.5, None),
    )
    R.log(f"  {[round(x, 3) for x in f15['per_bin']]}")
    plx = {}
    for d in DISTANCES_PC:
        tag = "none" if d is None else f"{d:g}pc"
        R.log(f"referee 2: parallax injections, 3 mJy, {tag} ...")
        plx[tag] = _cached(
            work / f"ref2_fine_3mJy_{tag}.json",
            lambda d=d: _fine(cats, model, size_noise, thr, rule, args.n_inject, 3.0, d),
        )
        R.log(f"  {[round(x, 3) for x in plx[tag]['per_bin']]}")
    ref_worst = worst_above(plx["none"])
    by_d = {}
    for tag, fine in plx.items():
        w = worst_above(fine)
        by_d[tag] = {
            "worst_1p1_5": w,
            "relative_loss": 1.0 - w / ref_worst,
            "limit_per_deg2_95": v.surface_density_limit(0, area, w),
        }
    ok = [d for d in DISTANCES_PC[1:] if by_d[f"{d:g}pc"]["relative_loss"] <= DOMAIN_LOSS]
    d_min = min(ok) if ok else None
    tag_min = f"{d_min:g}pc" if d_min is not None else "none"
    w15 = worst_above(f15)
    res = {
        "source": "real: VLASS QL epochs (vlasspm work directory); scripts/vlasspm_referee2.py",
        "is_real": True,
        "fine_completeness_1p5mJy_cut": f15,
        "parallax_injections_3mJy_cut": plx,
        "parallax_domain": {
            "rate_lo_arcsec_yr": RATE_LO,
            "criterion": f"smallest tested distance whose worst fine-bin completeness over "
            f"{RATE_LO}-5 arcsec/yr is within {DOMAIN_LOSS:.0%} (relative) of the no-parallax run",
            "by_distance": by_d,
            "distance_min_pc": d_min,
            "v_tan_min_km_s": None if d_min is None else float(KM_S_PER_AU_YR * RATE_LO * d_min),
        },
        "limits": {
            "rate_range_arcsec_yr": [RATE_LO, 5.0],
            "weighting": "none: worst fine bin over the range",
            "limit_3mJy_at_domain_distance": by_d[tag_min]["limit_per_deg2_95"],
            "completeness_3mJy_at_domain_distance": by_d[tag_min]["worst_1p1_5"],
            "limit_3mJy_no_parallax": by_d["none"]["limit_per_deg2_95"],
            "completeness_3mJy_no_parallax": by_d["none"]["worst_1p1_5"],
            "limit_1p5mJy_no_parallax": v.surface_density_limit(0, area, w15),
            "completeness_1p5mJy_no_parallax": w15,
            "area_deg2": area,
        },
        "uvcet_vectors": uv_vectors(ref1),
    }
    from jansky_research.report import write_results

    write_results(res, out / "results" / "vlasspm_referee2.json")
    R.log(f"wrote results/vlasspm_referee2.json; domain {d_min} pc; {by_d}")
    v.write_real_paper(out)
    R.log("paper macros + figure written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
