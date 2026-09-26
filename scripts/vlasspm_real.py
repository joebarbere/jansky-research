#!/usr/bin/env python
"""Plan 64 real leg: the blind VLASS proper-motion search over all four Quick-Look epochs.

Staged and checkpointed so it can run detached (``systemd-run --user``) and resume after an
interruption: every stage writes its output under ``--work`` and is skipped when that output
exists. Nothing here writes into ``results/`` except the final stage.

Stages
  1 load        each epoch's clean components -> ``epoch<N>.npz`` (per-component decimal year)
  2 floors      per-epoch astrometric floor measured from bright static sources
  3 search      orphans -> linkage -> collinearity -> flux, for every epoch triple
  4 null        RA-scramble re-linkage per triple (expected chance candidates)
  5 complete    injected movers through the same pipeline (completeness vs rate)
  6 uvcet       the recover-a-known: is UV Ceti among the candidates?
  7 vet         Gaia DR3 / CatWISE2020 counterparts for every candidate (network)
  8 write       results/vlasspm_metrics.json + results/vlasspm_candidates.csv

Data (gitignored, fetched by hand -- see survey/vlasspm-findings.md for URLs):
  E1  data/vlass/CIRADA_VLASS1QLv3.1_table1_components.csv.gz + ..._table3_subtile.csv.gz
  E2  data/CIRADA_VLASS2QLv2_table1_components.csv.gz + data/vlass/CIRADA_VLASS2QLv1_table3_subtile.csv.gz
  E3  data/QL3.1_components.fits + data/QL3.2_components.fits      (per-component MJD)
  E4  data/vlass/QL4.1_components_v1.fits                           (per-component MJD)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research import vlasspm as v  # noqa: E402

DATA = Path("/home/joe/dev/github/joebarbere/jansky-research/data")  # shared, not the worktree's
EPOCHS = {
    1: {
        "comp": [DATA / "vlass/CIRADA_VLASS1QLv3.1_table1_components.csv.gz"],
        "subtile": DATA / "vlass/CIRADA_VLASS1QLv3_table3_subtile.csv.gz",
    },
    2: {
        "comp": [DATA / "CIRADA_VLASS2QLv2_table1_components.csv.gz"],
        "subtile": DATA / "vlass/CIRADA_VLASS2QLv1_table3_subtile.csv.gz",
    },
    3: {"comp": [DATA / "QL3.1_components.fits", DATA / "QL3.2_components.fits"]},
    4: {"comp": [DATA / "vlass/QL4.1_components_v1.fits"]},
}
# Quick-Look peak fluxes are low by an epoch-dependent factor (VLASS Memos 13/22; the same
# constants as vlass.VLASS_PEAK_CORRECTION). No published E4 factor yet: 1.0, and the flux
# cut's factor-3 tolerance absorbs it.
PEAK_CORRECTION = {1: 1.13, 2: 1.075, 3: 1.031, 4: 1.0}
UVCET_GAIA_DR3 = 5140693571158946048
AREA_CELL_DEG = 0.5
T0 = time.time()


def log(msg: str) -> None:
    print(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')} +{(time.time() - T0) / 60:6.1f} min] {msg}",
        flush=True,
    )


def _iso_to_year(s: str) -> float:
    d = datetime.fromisoformat(s.strip().replace("Z", ""))
    start = datetime(d.year, 1, 1)
    return (
        d.year + (d - start).total_seconds() / (datetime(d.year + 1, 1, 1) - start).total_seconds()
    )


def _mjd_to_year(mjd: np.ndarray) -> np.ndarray:
    return 2000.0 + (np.asarray(mjd, float) - 51544.5) / 365.25


def load_cirada(epoch: int) -> dict:
    """Stream-parse a CIRADA component CSV with the Gordon+2021 clean cuts; join subtile dates."""
    dates: dict[str, float] = {}
    with gzip.open(EPOCHS[epoch]["subtile"], "rt") as fh:
        for r in csv.DictReader(fh):
            dates[r["Subtile"]] = _iso_to_year(r["DATEOBS"])
    cols = {k: [] for k in ("ra", "dec", "t", "flux", "err", "row")}
    n_all = n_nodate = 0
    with gzip.open(EPOCHS[epoch]["comp"][0], "rt") as fh:
        rd = csv.reader(fh)
        h = next(rd)
        ix = {
            k: h.index(k)
            for k in (
                "RA",
                "DEC",
                "E_RA",
                "E_DEC",
                "Peak_flux",
                "Duplicate_flag",
                "Quality_flag",
                "S_Code",
                "Subtile",
            )
        }
        for row_i, r in enumerate(rd):
            n_all += 1
            try:
                if int(float(r[ix["Duplicate_flag"]])) >= 2:
                    continue
                if int(float(r[ix["Quality_flag"]])) not in (0, 4):
                    continue
                if r[ix["S_Code"]].strip() == "E":
                    continue
                t = dates.get(r[ix["Subtile"]])
                if t is None:
                    n_nodate += 1
                    continue
                dec = float(r[ix["DEC"]])
                era = float(r[ix["E_RA"]]) * np.cos(np.radians(dec))
                edec = float(r[ix["E_DEC"]])
                cols["ra"].append(float(r[ix["RA"]]))
                cols["dec"].append(dec)
                cols["t"].append(t)
                cols["flux"].append(float(r[ix["Peak_flux"]]) * PEAK_CORRECTION[epoch])
                cols["err"].append(3600.0 * np.sqrt(0.5 * (era**2 + edec**2)))  # per-axis, arcsec
                cols["row"].append(row_i)
            except (ValueError, IndexError):
                continue
    out = {k: np.asarray(x, np.int64 if k == "row" else float) for k, x in cols.items()}
    out["n_all"], out["n_nodate"] = n_all, n_nodate
    return out


def load_ql_fits(epoch: int) -> dict:
    """QL3.x/QL4.1 FITS: ``Flag == 0`` (Memo 22 sidelobe flag), per-component MJD."""
    from astropy.io import fits

    parts = {k: [] for k in ("ra", "dec", "t", "flux", "err", "row")}
    n_all, offset = 0, 0
    for path in EPOCHS[epoch]["comp"]:
        with fits.open(path, memmap=True) as hd:
            d = hd[1].data
            n = len(d)
            n_all += n
            keep = (np.asarray(d["Flag"]) == 0) & (np.asarray(d["S_Code"]).astype(str) != "E")
            dec = np.asarray(d["DEC"], float)[keep]
            era = np.asarray(d["E_RA"], float)[keep] * np.cos(np.radians(dec))
            edec = np.asarray(d["E_DEC"], float)[keep]
            parts["ra"].append(np.asarray(d["RA"], float)[keep])
            parts["dec"].append(dec)
            parts["t"].append(_mjd_to_year(np.asarray(d["MJD"], float)[keep]))
            parts["flux"].append(np.asarray(d["Peak_flux"], float)[keep] * PEAK_CORRECTION[epoch])
            parts["err"].append(3600.0 * np.sqrt(0.5 * (era**2 + edec**2)))
            parts["row"].append(offset + np.flatnonzero(keep))
            offset += n
    out = {k: np.concatenate(x) for k, x in parts.items()}
    out["n_all"], out["n_nodate"] = n_all, 0
    return out


def stage_load(work: Path) -> None:
    for e in EPOCHS:
        f = work / f"epoch{e}.npz"
        if f.exists():
            continue
        log(f"load E{e} ...")
        d = load_cirada(e) if e in (1, 2) else load_ql_fits(e)
        good = np.isfinite(d["ra"]) & np.isfinite(d["dec"]) & np.isfinite(d["t"]) & (d["err"] > 0)
        np.savez(
            f,
            **{k: d[k][good] for k in ("ra", "dec", "t", "flux", "err", "row")},
            n_all=d["n_all"],
            n_nodate=d["n_nodate"],
        )
        log(
            f"  E{e}: {int(good.sum()):,} clean of {d['n_all']:,} (no date: {d['n_nodate']:,}); "
            f"t {d['t'][good].min():.2f}-{d['t'][good].max():.2f}"
        )


def catalogs(
    work: Path, floors: list[float] | None, banded: dict | None = None
) -> list[v.EpochCatalog]:
    """Clean epochs with per-component errors = hypot(catalogue error, floor).

    With ``banded`` (from :func:`vlasspm.calibrate_floors_banded`) the floor depends on the
    component's declination band; bands without enough bright matches use the all-sky floor.
    """
    cats = []
    for k, e in enumerate(EPOCHS):
        d = np.load(work / f"epoch{e}.npz")
        if floors is None:
            err = d["err"]
        elif banded is None:
            err = np.hypot(d["err"], floors[k])
        else:
            err = np.hypot(d["err"], v.banded_floor(d["dec"], k, banded, fallback=floors[k]))
        cats.append(v.EpochCatalog(d["ra"], d["dec"], d["t"], d["flux"], err, d["row"]))
    return cats


def _save(path: Path, obj) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=float))
    tmp.replace(path)  # atomic: a killed job never leaves a half-written checkpoint


def _tname(t) -> str:
    return "-".join(f"E{x + 1}" for x in t)


def footprint_area(cats: list[v.EpochCatalog]) -> float:
    """deg^2 of equal-area cells occupied in EVERY given epoch (the triple's common sky)."""
    nra = int(360 / AREA_CELL_DEG)
    ds = np.radians(AREA_CELL_DEG)  # d(sin dec) step giving ~AREA_CELL_DEG^2 cells at the equator
    occupied = None
    for c in cats:
        cell = set(
            zip(
                (c.ra // AREA_CELL_DEG).astype(int) % nra,
                np.floor(np.sin(np.radians(c.dec)) / ds).astype(int),
                strict=True,
            )
        )
        occupied = cell if occupied is None else occupied & cell
    cell_sr = np.radians(AREA_CELL_DEG) * ds
    return float(len(occupied or ()) * cell_sr * (180 / np.pi) ** 2)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=str(DATA / "vlass/work"))
    ap.add_argument("--n-null", type=int, default=50)
    ap.add_argument("--n-inject", type=int, default=20000)
    ap.add_argument("--out", default=str(REPO))
    args = ap.parse_args()
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    log(f"plan 64 real leg; work={work} out={args.out}")

    stage_load(work)

    f = work / "floors.json"
    if not f.exists():
        log("floors: bright-source epoch-pair scatter ...")
        _save(f, v.calibrate_floors(catalogs(work, None), flux_min_mjy=10.0))
    floors = json.loads(f.read_text())
    log(f"  floors (arcsec): {floors['floor_arcsec']}; pairs: {floors['pair_sigma_arcsec']}")
    f = work / "floors_banded.json"
    if not f.exists():
        log("floors: per declination band ...")
        _save(f, v.calibrate_floors_banded(catalogs(work, None), flux_min_mjy=10.0))
    banded = json.loads(f.read_text())
    for b in banded["bands"]:
        fl = b["floor_arcsec"]
        log(
            f"  Dec [{b['lo']:+.0f},{b['hi']:+.0f}): "
            + (str([round(x, 3) for x in fl]) if fl else "too few")
        )
    cats = catalogs(work, floors["floor_arcsec"], banded)
    triples = v.epoch_triples(len(cats))

    f = work / "search.json"
    if not f.exists():
        log("search: every epoch triple ...")
        res = v.search_multi(cats, triples)
        out, cands = {}, []
        for t, r in res.items():
            a, b, c = r.orphans
            cd = r.candidates
            out[_tname(t)] = {
                "n_orphans": list(r.n_orphans),
                "n_pairs": len(r.pairs),
                "n_triplets": len(r.triplets),
                "n_candidates": len(cd),
                "area_deg2": footprint_area([cats[x] for x in t]),
            }
            for n in range(len(cd)):
                i, j, k = cd.i[n], cd.j[n], cd.k[n]
                cands.append(
                    {
                        "triple": _tname(t),
                        "ra1": a.ra[i],
                        "dec1": a.dec[i],
                        "t1": a.t_yr[i],
                        "row1": int(a.ident[i]),
                        "ra2": b.ra[j],
                        "dec2": b.dec[j],
                        "t2": b.t_yr[j],
                        "row2": int(b.ident[j]),
                        "ra3": c.ra[k],
                        "dec3": c.dec[k],
                        "t3": c.t_yr[k],
                        "row3": int(c.ident[k]),
                        "mu_ra": cd.mu_ra[n],
                        "mu_dec": cd.mu_dec[n],
                        "mu": cd.mu[n],
                        "resid_sigma": cd.resid_sigma[n],
                        "flux1": a.flux[i],
                        "flux2": b.flux[j],
                        "flux3": c.flux[k],
                    }
                )
            log(f"  {_tname(t)}: {out[_tname(t)]}")
        _save(f, {"per_triple": out, "candidates": cands})
    search = json.loads(f.read_text())

    for t in triples:
        f = work / f"null_{_tname(t)}.json"
        if f.exists():
            continue
        log(f"null {_tname(t)}: {args.n_null} scrambles ...")
        r = v.search(cats[t[0]], cats[t[1]], cats[t[2]])
        _save(f, v.scramble_null(r.orphans, n_reps=args.n_null, seed=hash(t) % 2**31))
        log(
            f"  {_tname(t)} chance candidates/scramble: {json.loads(f.read_text())['candidates_mean']:.2f}"
        )

    for flux in (3.0, 1.5):
        f = work / f"completeness_{flux:g}mJy.json"
        if f.exists():
            continue
        log(f"completeness: {args.n_inject} injected movers at {flux} mJy ...")
        _save(f, v.completeness(*cats, n=args.n_inject, flux_mjy=flux, seed=int(flux * 10)))
        log(f"  overall {json.loads(f.read_text())['overall']:.3f}")

    f = work / "uvcet.json"
    if not f.exists():
        log("uvcet: recover-a-known ...")
        _save(f, uvcet_check(search["candidates"], cats))
    log(f"  {json.loads(f.read_text())}")

    f = work / "vet.json"
    if not f.exists():
        log(f"vet: Gaia/CatWISE counterparts for {len(search['candidates'])} candidates ...")
        _save(f, vet_counterparts(search["candidates"]))

    write_outputs(work, Path(args.out), floors, search)
    log("done")
    return 0


def _track(c: dict, t: float) -> tuple[float, float]:
    ra, dec = v._mover_track(c["ra1"], c["dec1"], c["mu_ra"], c["mu_dec"], c["t1"], t)
    return float(ra), float(dec)


def uvcet_check(cands: list[dict], cats: list | None = None) -> dict:  # network
    """Is UV Ceti (Gaia DR3 5140693571158946048) recovered, and by which triple?"""
    try:
        from astroquery.vizier import Vizier

        g = Vizier(columns=["RA_ICRS", "DE_ICRS", "pmRA", "pmDE"], row_limit=5).query_constraints(
            catalog="I/355/gaiadr3", Source=str(UVCET_GAIA_DR3)
        )[0]
        ra, dec, pmra, pmde = (float(g[k][0]) for k in ("RA_ICRS", "DE_ICRS", "pmRA", "pmDE"))
    except Exception as exc:  # noqa: BLE001
        return {"error": f"Gaia lookup failed: {exc!r}"}
    hits = []
    for c in cands:
        # Where UV Cet was at the candidate's middle-epoch time, vs the candidate's position.
        rp, dp = v._mover_track(ra, dec, pmra / 1000, pmde / 1000, 2016.0, c["t2"])
        dx, dy = v.tangent_offsets_arcsec(rp, dp, c["ra2"], c["dec2"])
        sep = float(np.hypot(dx, dy))
        if sep < 5.0:
            hits.append(
                {
                    "triple": c["triple"],
                    "sep_arcsec": sep,
                    "mu": c["mu"],
                    "gaia_mu": float(np.hypot(pmra, pmde)) / 1000,
                }
            )
    per_epoch: dict = {}
    for e, cat in enumerate(cats or []):
        # Where UV Cet should be in this epoch: nearest component, its flux and isolation.
        rp, dp = v._mover_track(
            ra, dec, pmra / 1000, pmde / 1000, 2016.0, float(np.median(cat.t_yr))
        )
        near = np.flatnonzero((np.abs(cat.dec - dp) < 0.01) & (np.abs(cat.ra - rp) < 0.02))
        if near.size == 0:
            per_epoch[f"E{e + 1}"] = "no component within ~36 arcsec"
            continue
        rp2, dp2 = v._mover_track(ra, dec, pmra / 1000, pmde / 1000, 2016.0, cat.t_yr[near])
        dx, dy = v.tangent_offsets_arcsec(rp2, dp2, cat.ra[near], cat.dec[near])
        k = int(np.argmin(np.hypot(dx, dy)))
        dxy = v.tangent_offsets_arcsec(
            cat.ra[near[k]], cat.dec[near[k]], cat.ra[near], cat.dec[near]
        )
        per_epoch[f"E{e + 1}"] = {
            "sep_arcsec": float(np.hypot(dx[k], dy[k])),
            "flux_mjy": float(cat.flux[near[k]]),
            "n_other_within_30": int((np.hypot(*dxy) < 30.0).sum() - 1),
        }
    return {
        "gaia_mu_arcsec_yr": float(np.hypot(pmra, pmde)) / 1000,
        "recovered": bool(hits),
        "hits": hits,
        "per_epoch": per_epoch,
    }


def vet_counterparts(cands: list[dict], radius_arcsec: float = 5.0) -> list[dict]:  # network
    """Gaia DR3 / CatWISE2020 near each candidate's track at the catalogue epoch.

    For Gaia, the NEAREST match's proper motion is recorded and compared with the radio one:
    agreement (vector difference < 30% of |mu|) marks a known star recovered blind -- the
    positive control. Disagreement or no match leaves the candidate unexplained.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    gaia = Vizier(columns=["RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "Gmag"], row_limit=20)
    wise = Vizier(columns=["*"], row_limit=20)
    out = []
    for n, c in enumerate(cands):
        rec: dict = {"index": n}
        ra, dec = _track(c, 2016.0)
        try:
            r = gaia.query_region(
                SkyCoord(ra * u.deg, dec * u.deg),
                radius=radius_arcsec * u.arcsec,
                catalog="I/355/gaiadr3",
            )
            rec["gaia"] = int(len(r[0])) if len(r) else 0
            if rec["gaia"]:
                t = r[0]
                dx, dy = v.tangent_offsets_arcsec(
                    ra, dec, np.asarray(t["RA_ICRS"]), np.asarray(t["DE_ICRS"])
                )
                k = int(np.argmin(np.hypot(dx, dy)))
                pmra, pmde = float(t["pmRA"][k]) / 1000, float(t["pmDE"][k]) / 1000  # arcsec/yr
                rec |= {
                    "gaia_sep": float(np.hypot(dx[k], dy[k])),
                    "gaia_pmra": pmra,
                    "gaia_pmde": pmde,
                    "gaia_gmag": float(t["Gmag"][k]),
                }
                if np.isfinite(pmra) and np.isfinite(pmde):
                    diff = float(np.hypot(c["mu_ra"] - pmra, c["mu_dec"] - pmde))
                    rec["pm_diff"] = diff
                    rec["pm_agree"] = bool(diff < 0.3 * max(c["mu"], 1e-9))
        except Exception as exc:  # noqa: BLE001
            rec["gaia"] = f"error: {exc!r}"
        ra, dec = _track(c, 2015.4)
        try:
            r = wise.query_region(
                SkyCoord(ra * u.deg, dec * u.deg),
                radius=radius_arcsec * u.arcsec,
                catalog="II/365/catwise",
            )
            rec["catwise"] = int(len(r[0])) if len(r) else 0
        except Exception as exc:  # noqa: BLE001
            rec["catwise"] = f"error: {exc!r}"
        out.append(rec)
    return out


def write_outputs(work: Path, out: Path, floors: dict, search: dict) -> None:
    from jansky_research.report import write_results

    nulls = {p.stem[5:]: json.loads(p.read_text()) for p in sorted(work.glob("null_*.json"))}
    comp = {
        p.stem[13:]: json.loads(p.read_text()) for p in sorted(work.glob("completeness_*.json"))
    }
    uv = json.loads((work / "uvcet.json").read_text())
    vet = json.loads((work / "vet.json").read_text())
    dark = [
        c | {"vet": vt}
        for c, vt in zip(search["candidates"], vet, strict=True)
        if vt.get("gaia") == 0 and vt.get("catwise") == 0
    ]
    area = max(pt["area_deg2"] for pt in search["per_triple"].values())
    comp3 = comp.get("3mJy", {}).get("overall", 0.0)
    metrics = {
        "source": "real: VLASS Quick-Look epochs 1, 2, 3 (QL3.1+3.2) and 4.1 component catalogues",
        "is_real": True,
        "real_epochs": {
            f"E{e}": {
                "n_clean": int(len(np.load(work / f"epoch{e}.npz")["ra"])),
                "t_min": float(np.load(work / f"epoch{e}.npz")["t"].min()),
                "t_max": float(np.load(work / f"epoch{e}.npz")["t"].max()),
            }
            for e in EPOCHS
        },
        "real_floors": floors,
        "real_floors_banded": json.loads((work / "floors_banded.json").read_text()),
        "real_per_triple": search["per_triple"],
        "real_null": {
            k: {kk: vv for kk, vv in n.items() if kk != "candidates_per_rep"}
            for k, n in nulls.items()
        },
        "real_completeness": comp,
        "real_n_candidates": len(search["candidates"]),
        "real_n_optically_dark": len(dark),
        "real_n_gaia_pm_agree": sum(1 for vt in vet if vt.get("pm_agree") is True),
        "real_n_gaia_pm_disagree": sum(1 for vt in vet if vt.get("pm_agree") is False),
        "real_uvcet": uv,
        "real_area_deg2_max_triple": area,
        "real_density_limit_per_deg2_3mJy": v.surface_density_limit(len(dark), area, comp3),
    }
    (out / "results").mkdir(parents=True, exist_ok=True)
    write_results(metrics, out / "results" / "vlasspm_metrics.json")
    keys = list(search["candidates"][0]) if search["candidates"] else []
    with (out / "results" / "vlasspm_candidates.csv").open("w", newline="") as fh:
        extra = [
            "gaia",
            "catwise",
            "gaia_sep",
            "gaia_pmra",
            "gaia_pmde",
            "gaia_gmag",
            "pm_diff",
            "pm_agree",
        ]
        w = csv.DictWriter(fh, fieldnames=keys + extra)
        w.writeheader()
        for c, vt in zip(search["candidates"], vet, strict=True):
            w.writerow(c | {k: vt.get(k) for k in extra})
    log(f"wrote results: {len(search['candidates'])} candidates, {len(dark)} optically dark")


if __name__ == "__main__":
    raise SystemExit(main())
