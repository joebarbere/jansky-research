#!/usr/bin/env python
"""Plan 64 / vlasspm note, referee round 1: the extra measurements, written as evidence.

    uv run python scripts/vlasspm_referee1.py --out .      # -> results/vlasspm_referee1.json

Reads the real-leg work directory (``data/vlass/work``: epoch catalogues, floors, the fitted
error model, the run-4 search and the run-5 cut) and measures:

1. UV Ceti: a three-epoch straight-line fit to the committed E2/E3/E4 detections, weighted by
   each detection's covariance, with the parallax fixed at Gaia DR3's value, without parallax,
   and free; compared with Gaia DR3 for both components (TAP) and the system motion (SIMBAD).
2. Chance coincidence: P(>= 1 Gaia DR3 / CatWISE2020 source within r) at random footprint
   positions (fixed seed), for the 5" search radius and smaller radii.
3. Completeness in fine rate bins around the 0.92"/yr edge (3 mJy, with the cut).
4. The completeness the heavy error tails could cost, bounded from the static validation table.
5. The distance below which parallax (ignored by the injections) spoils the E3 test.
6. The worst-bin surface-density limit next to the averaged one.

Network: steps 1 and 2 query Gaia/SIMBAD/VizieR. Every stage is cached in the work directory.
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

from jansky_research import vlasspm as v  # noqa: E402

UVCET, BLCET = 5140693571158946048, 5140693571158739840  # Gaia DR3 source_id (SIMBAD: GJ 65B/A)
FINE_BINS = (0.8, 0.92, 1.0, 1.1, 1.25, 1.62, 2.85, 5.0)
CHANCE_N, CHANCE_SEED = 500, 20260927
CHANCE_RADII = (0.5, 1.0, 1.5, 2.0, 3.0, 5.0)


def _cached(path: Path, fn):
    if path.exists():
        return json.loads(path.read_text())
    out = fn()
    R._save(path, out)
    return out


def gaia_system() -> dict:  # network
    from astroquery.gaia import Gaia
    from astroquery.simbad import Simbad

    q = (
        "SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, pmra_error, "
        "pmdec_error, ruwe FROM gaiadr3.gaia_source WHERE source_id IN "
        f"({UVCET}, {BLCET})"
    )
    rows = Gaia.launch_job(q).get_results()
    comp = {}
    for r in rows:
        name = "UV Cet" if int(r["source_id"]) == UVCET else "BL Cet"
        comp[name] = {
            k: float(r[k])
            for k in (
                "ra",
                "dec",
                "parallax",
                "parallax_error",
                "pmra",
                "pmdec",
                "pmra_error",
                "pmdec_error",
                "ruwe",
            )
        } | {"source_id": int(r["source_id"])}
        comp[name]["mu_arcsec_yr"] = float(np.hypot(r["pmra"], r["pmdec"])) / 1000
    s = Simbad()
    s.add_votable_fields("pmra", "pmdec", "plx_value", "pm_bibcode")
    t = s.query_object("GJ 65")
    system = {
        "simbad_main_id": str(t["main_id"][0]),
        "pmra_mas_yr": float(t["pmra"][0]),
        "pmdec_mas_yr": float(t["pmdec"][0]),
        "mu_arcsec_yr": float(np.hypot(t["pmra"][0], t["pmdec"][0])) / 1000,
        "plx_mas": float(t["plx_value"][0]),
        "pm_bibcode": str(t["pm_bibcode"][0]),
    }
    return {"components": comp, "system": system}


def uvcet_fit(cats, model, search: dict, gaia: dict) -> dict:
    c = next(
        x
        for x in search["candidates"]
        if x["triple"] == "E2-E3-E4" and abs(x["ra2"] - 24.7787) < 0.01
    )
    ra, dec, t, cov = R.candidate_detections(c, cats, model)
    plx = gaia["components"]["UV Cet"]["parallax"] / 1000
    fits = {
        "fixed_gaia_parallax": v.fit_track(ra, dec, t, cov, parallax_arcsec=plx),
        "no_parallax": v.fit_track(ra, dec, t, cov, parallax_arcsec=0.0),
        "free_parallax": v.fit_track(ra, dec, t, cov, parallax_arcsec=None),
    }
    f = v.parallax_residual_factor(ra[0], dec[0], [t[0]], [t[1]], [t[2]])[0]
    return {
        "detections": {
            "ra": ra.tolist(),
            "dec": dec.tolist(),
            "t": t.tolist(),
            "cov": cov.tolist(),
        },
        "parallax_used_arcsec": plx,
        "fits": fits,
        "triplet_rate_search": c["mu"],
        "e3_parallax_residual_arcsec": float(f * plx),
    }


def chance(cats) -> dict:  # network
    """Random footprint positions: a random E1 component displaced by 2-5 arcmin (so the point
    is in the surveyed sky but not on the radio source)."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    rng = np.random.default_rng(CHANCE_SEED)
    e1 = cats[0]
    host = rng.integers(0, len(e1), CHANCE_N)
    off = rng.uniform(120, 300, CHANCE_N) / 3600
    ang = rng.uniform(0, 2 * np.pi, CHANCE_N)
    ra = (e1.ra[host] + off * np.cos(ang) / np.cos(np.radians(e1.dec[host]))) % 360
    dec = e1.dec[host] + off * np.sin(ang)
    near = {}
    for name, cat, cols in (
        ("gaia", "I/355/gaiadr3", ["RA_ICRS", "DE_ICRS"]),
        ("catwise", "II/365/catwise", ["RAPMdeg", "DEPMdeg"]),
    ):
        d = np.full(CHANCE_N, np.inf)
        viz = Vizier(columns=cols, row_limit=-1)
        for lo in range(0, CHANCE_N, 100):
            sl = slice(lo, lo + 100)
            res = viz.query_region(
                SkyCoord(ra[sl] * u.deg, dec[sl] * u.deg), radius=5.0 * u.arcsec, catalog=cat
            )
            if not len(res):
                continue
            tab = res[0]
            q = np.asarray(tab["_q"], int) - 1 + lo
            dx, dy = v.tangent_offsets_arcsec(
                ra[q], dec[q], np.asarray(tab[cols[0]], float), np.asarray(tab[cols[1]], float)
            )
            sep = np.hypot(dx, dy)
            np.minimum.at(d, q, sep)
        near[name] = d
    both = np.minimum(near["gaia"], near["catwise"])
    return {
        "n_positions": CHANCE_N,
        "seed": CHANCE_SEED,
        "radii_arcsec": list(CHANCE_RADII),
        "p_gaia": v.chance_within(near["gaia"].tolist(), CHANCE_RADII),
        "p_catwise": v.chance_within(near["catwise"].tolist(), CHANCE_RADII),
        "p_either": v.chance_within(both.tolist(), CHANCE_RADII),
    }


def fine_completeness(cats, model, size_noise, thr, rule, n: int) -> dict:
    return v.completeness(
        *cats,
        n=n,
        flux_mjy=3.0,
        seed=97,
        realistic_frac=1.0,
        model=model,
        size_noise=size_noise,
        compact_max=thr,
        compact_rule=rule,
        bins=np.asarray(FINE_BINS),
    )


def tail_bound(cats, calib: dict, fluxes=(3.0, 1.5)) -> dict:
    """Upper bound on the completeness the heavy error tails can cost.

    A mover is lost when its E3 residual lands beyond 3 sigma; under the model that happens with
    Rayleigh probability 0.0111, and the injections reproduce it. The real residuals of compact
    static sources exceed 3 sigma more often, by f_gt3 - 0.0111, in the validation table. The
    bound takes the largest such excess over the S/N bins a detection at that flux can occupy
    (per-epoch peak 0.8-1.25 x the flux, over the 10th-90th percentile of the local rms)."""
    rms = np.concatenate([c.rms for c in cats])
    r10, r90 = np.percentile(rms, [10, 90])
    joint = calib["validation_snr_x_size"]
    out = {"rms_p10_p90_mjy": [float(r10), float(r90)]}
    for f in fluxes:
        lo, hi = 0.8 * f / r90, 1.25 * f / r10
        excess = []
        for row in joint:
            if row["snr_hi"] <= lo or row["snr_lo"] >= hi:
                continue
            compact = next(x for x in row["by_size"] if x["lo"] == 0.0)
            excess.append(compact["f_gt3"] - v.RAYLEIGH_TAIL_3SIGMA)
        out[f"{f:g}mJy"] = {"snr_range": [float(lo), float(hi)], "max_excess": float(max(excess))}
    return out


def parallax_floor(cats, model, n: int = 2000, seed: int = 5) -> dict:
    """Distance below which the parallax the injections ignore adds more than one standard error
    to the E3 residual of a 3 mJy point source (median over E1-E2-E3 sky positions and dates)."""
    rng = np.random.default_rng(seed)
    e1, e2, e3 = cats[:3]
    k = rng.integers(0, len(e1), n)
    ra, dec = e1.ra[k], e1.dec[k]
    from scipy.spatial import cKDTree

    x = v._xyz(ra, dec)
    j2 = cKDTree(v._xyz(e2.ra, e2.dec)).query(x)[1]
    j3 = cKDTree(v._xyz(e3.ra, e3.dec)).query(x)[1]
    t1, t2, t3 = e1.t_yr[k], e2.t_yr[j2], e3.t_yr[j3]
    fac = np.array(
        [v.parallax_residual_factor(ra[i], dec[i], [t1[i]], [t2[i]], [t3[i]])[0] for i in range(n)]
    )
    r = (t3 - t2) / (t2 - t1)
    s = []
    for cat, jj in ((e1, k), (e2, j2), (e3, j3)):
        cov, _, _, _, _ = v._injected_errors(cat, jj, np.full(n, 3.0), 0.0)
        s.append(0.5 * (cov[:, 0] + cov[:, 1]))
    sig = np.sqrt(r**2 * s[0] + (1 + r) ** 2 * s[1] + s[2])
    plx_max = np.median(sig / fac)  # arcsec of parallax that shifts E3 by 1 sigma
    return {
        "n_positions": n,
        "median_factor": float(np.median(fac)),
        "median_sigma_e3_arcsec": float(np.median(sig)),
        "parallax_max_arcsec": float(plx_max),
        "distance_min_pc": float(1.0 / plx_max),
    }


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
    search = json.loads((work / "search.json").read_text())
    r5 = json.loads((work / "run5_calibration.json").read_text())
    stars = json.loads((work / "run5_stars.json").read_text())
    size_noise = v.SizeNoiseModel(np.asarray(stars["samples"]))
    R.log("referee 1: Gaia / SIMBAD for the UV/BL Cet system ...")
    gaia = _cached(work / "ref1_gaia.json", gaia_system)
    R.log("referee 1: UV Ceti track fits ...")
    uv = _cached(work / "ref1_uvcet_fit.json", lambda: uvcet_fit(cats, model, search, gaia))
    R.log(f"  { ({k: round(f['mu'], 3) for k, f in uv['fits'].items()}) }")
    R.log("referee 1: chance coincidence at random positions ...")
    ch = _cached(work / "ref1_chance.json", lambda: chance(cats))
    R.log(f"  p_either {ch['p_either']}")
    R.log("referee 1: fine-bin completeness ...")
    fine = _cached(
        work / "ref1_fine_completeness.json",
        lambda: fine_completeness(
            cats, model, size_noise, r5["threshold"], r5["rule"], args.n_inject
        ),
    )
    R.log(f"  {[round(x, 3) for x in fine['per_bin']]}")
    tb = _cached(work / "ref1_tail_bound.json", lambda: tail_bound(cats, calib))
    pf = _cached(work / "ref1_parallax_floor.json", lambda: parallax_floor(cats, model))
    R.log(f"  tail bound {tb}; parallax floor {pf}")
    metrics = json.loads((out / "results" / "vlasspm_metrics.json").read_text())
    r5m = metrics["real_run5_compactness"]
    comp = r5m["completeness"]["3mJy_cut"]
    area = r5m["area_deg2"]
    edges = np.asarray(comp["bin_edges"])
    above = [p for p, lo in zip(comp["per_bin"], edges[:-1], strict=True) if lo >= 0.9]
    worst = float(min(above))
    fine_above = [p for p, lo in zip(fine["per_bin"], FINE_BINS[:-1], strict=True) if lo >= 0.9]
    res = {
        "source": "real: VLASS QL epochs (vlasspm work directory) + Gaia DR3 / SIMBAD / VizieR; "
        "scripts/vlasspm_referee1.py",
        "is_real": True,
        "uvcet": {"gaia": gaia, **uv},
        "chance_coincidence": ch,
        "fine_completeness_3mJy_cut": {
            k: fine[k] for k in ("bin_edges", "per_bin", "n_per_bin", "overall", "n_injected")
        },
        "tail_completeness_bound": tb,
        "parallax_floor": pf,
        "limits": {
            "rate_weighting": "unweighted mean over the 0.92-5 arcsec/yr log bins = log-uniform prior",
            "worst_bin_completeness_3mJy_cut": worst,
            "worst_bin_limit_per_deg2_95": v.surface_density_limit(0, area, worst),
            "worst_fine_bin_completeness_3mJy_cut": float(min(fine_above)),
            "worst_fine_bin_limit_per_deg2_95": v.surface_density_limit(0, area, min(fine_above)),
            "area_deg2": area,
        },
    }
    from jansky_research.report import write_results

    write_results(res, out / "results" / "vlasspm_referee1.json")
    R.log("wrote results/vlasspm_referee1.json")
    v.write_real_paper(out)  # macros + figure pick up the new evidence
    R.log("paper macros + figure written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
