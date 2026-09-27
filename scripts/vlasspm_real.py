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


# Columns kept per component. Every epoch's table (CIRADA E1/E2 CSVs; NRAO QL3.x/QL4.1 FITS)
# carries the same PyBDSF shape set -- fitted Maj/Min/PA, deconvolved DC_Maj/DC_Min/DC_PA, the
# restoring beam BMAJ/BMIN/BPA per component (arcsec, deg), Peak_flux and Isl_rms (mJy/beam) --
# so no epoch needs image headers for the shape model. See survey/vlasspm-findings.md
# ("Shape-aware errors") for the per-epoch check.
SHAPE_COLS = ("Maj", "Min", "PA", "DC_Maj", "DC_Min", "DC_PA", "BMAJ", "BMIN", "BPA", "Isl_rms")
SHAPE_KEYS = ("maj", "min", "pa", "dcmaj", "dcmin", "dcpa", "bmaj", "bmin", "bpa", "rms")


def _sky_errors(e_ra_deg, e_dec_deg) -> tuple:
    """Catalogue E_RA/E_DEC (deg) -> on-sky arcsec.

    E_RA is ALREADY an on-sky angle in all four catalogues: it equals the Condon (1997)
    position error rotated by PA to 1 part in 10^3 without any cos(dec). Runs 1-3 multiplied it
    by cos(dec), shrinking RA errors by up to x2.7 at high declination (x1.15 at Dec -30).
    """
    return 3600.0 * np.asarray(e_ra_deg, float), 3600.0 * np.asarray(e_dec_deg, float)


def load_cirada(epoch: int) -> dict:
    """Stream-parse a CIRADA component CSV with the Gordon+2021 clean cuts; join subtile dates."""
    dates: dict[str, float] = {}
    with gzip.open(EPOCHS[epoch]["subtile"], "rt") as fh:
        for r in csv.DictReader(fh):
            dates[r["Subtile"]] = _iso_to_year(r["DATEOBS"])
    keys = ("ra", "dec", "t", "flux", "era", "edec", *SHAPE_KEYS, "row")
    cols: dict[str, list] = {k: [] for k in keys}
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
                *SHAPE_COLS,
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
                shape = [float(r[ix[c]]) for c in SHAPE_COLS]
                cols["ra"].append(float(r[ix["RA"]]))
                cols["dec"].append(float(r[ix["DEC"]]))
                cols["t"].append(t)
                cols["flux"].append(float(r[ix["Peak_flux"]]) * PEAK_CORRECTION[epoch])
                cols["era"].append(float(r[ix["E_RA"]]))
                cols["edec"].append(float(r[ix["E_DEC"]]))
                for k, x in zip(SHAPE_KEYS, shape, strict=True):
                    cols[k].append(x)
                cols["row"].append(row_i)
            except (ValueError, IndexError):
                continue
    out = {k: np.asarray(x, np.int64 if k == "row" else float) for k, x in cols.items()}
    out["era"], out["edec"] = _sky_errors(out["era"], out["edec"])
    out["n_all"], out["n_nodate"] = n_all, n_nodate
    return out


def load_ql_fits(epoch: int) -> dict:
    """QL3.x/QL4.1 FITS: ``Flag == 0`` (Memo 22 sidelobe flag), per-component MJD."""
    from astropy.io import fits

    keys = ("ra", "dec", "t", "flux", "era", "edec", *SHAPE_KEYS, "row")
    parts: dict[str, list] = {k: [] for k in keys}
    n_all, offset = 0, 0
    for path in EPOCHS[epoch]["comp"]:
        with fits.open(path, memmap=True) as hd:
            d = hd[1].data
            n = len(d)
            n_all += n
            keep = (np.asarray(d["Flag"]) == 0) & (np.asarray(d["S_Code"]).astype(str) != "E")
            era, edec = _sky_errors(d["E_RA"][keep], d["E_DEC"][keep])
            parts["ra"].append(np.asarray(d["RA"], float)[keep])
            parts["dec"].append(np.asarray(d["DEC"], float)[keep])
            parts["t"].append(_mjd_to_year(np.asarray(d["MJD"], float)[keep]))
            parts["flux"].append(np.asarray(d["Peak_flux"], float)[keep] * PEAK_CORRECTION[epoch])
            parts["era"].append(era)
            parts["edec"].append(edec)
            for k, c in zip(SHAPE_KEYS, SHAPE_COLS, strict=True):
                parts[k].append(np.asarray(d[c], float)[keep])
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
        good = (
            np.isfinite(d["ra"])
            & np.isfinite(d["dec"])
            & np.isfinite(d["t"])
            & (d["era"] > 0)
            & (d["edec"] > 0)
            & (d["maj"] > 0)
            & (d["min"] > 0)
            & (d["rms"] > 0)
        )
        keys = ("ra", "dec", "t", "flux", "era", "edec", *SHAPE_KEYS, "row")
        np.savez(f, **{k: d[k][good] for k in keys}, n_all=d["n_all"], n_nodate=d["n_nodate"])
        log(
            f"  E{e}: {int(good.sum()):,} clean of {d['n_all']:,} (no date: {d['n_nodate']:,}); "
            f"t {d['t'][good].min():.2f}-{d['t'][good].max():.2f}"
        )


def catalogs(
    work: Path, floors: list[float] | None, banded: dict | None = None
) -> list[v.EpochCatalog]:
    """Clean epochs with the shape-aware measurement covariance.

    ``cov`` = the Condon (1997) error ellipse recomputed from each component's fitted Maj/Min/PA,
    its beam and its peak S/N (identical to the catalogue's E_RA/E_DEC, but keeping the ellipse
    and the RA-Dec correlation), plus floor^2 on both axes. The floor is the bright-static floor
    of the component's declination band (``banded``) or the all-sky one; ``floors=None`` = none.
    ``shape``/``beam`` carry the deconvolved source and restoring beam for the systematic terms.
    """
    cats = []
    for k, e in enumerate(EPOCHS):
        d = np.load(work / f"epoch{e}.npz")
        n = d["ra"].size
        if floors is None:
            floor = np.zeros(n)
        elif banded is None:
            floor = np.full(n, floors[k])
        else:
            floor = v.banded_floor(d["dec"], k, banded, fallback=floors[k])
        rms = d["rms"] * PEAK_CORRECTION[e]  # so flux / rms is the catalogue's raw peak S/N
        snr = d["flux"] / rms
        cov = v.condon_cov(d["maj"], d["min"], d["pa"], d["bmaj"], d["bmin"], snr)
        cov[:, 0] += floor**2
        cov[:, 1] += floor**2
        cats.append(
            v.EpochCatalog(
                d["ra"],
                d["dec"],
                d["t"],
                d["flux"],
                np.sqrt(0.5 * (cov[:, 0] + cov[:, 1])),
                d["row"],
                cov,
                v.ellipse_cov(d["dcmaj"], d["dcmin"], d["dcpa"]),
                v.ellipse_cov(d["bmaj"], d["bmin"], d["bpa"]),
                rms,
                floor,
            )
        )
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


DEC_EDGES_FINE = (-90.0, -35.0, -30.0, -25.0, -20.0, 0.0, 30.0, 90.0)
# Run 3's bands were (-40, -20, 0, 30, 90): components south of -40 fell back to the all-sky
# floor, and the -40..-20 band averaged a floor that is 0.26" below -35 and 0.12" at -25..-20.
# The finer southern edges flatten the statics' Dec dependence (survey/vlasspm-findings.md).
LEGACY_DEC_EDGES = (-40.0, -20.0, 0.0, 30.0, 90.0)


def _cand_row(t, a, b, c, cd, n) -> dict:
    i, j, k = cd.i[n], cd.j[n], cd.k[n]
    return {
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
        "snr1": a.flux[i] / a.rms[i],
        "snr2": b.flux[j] / b.rms[j],
        "snr3": c.flux[k] / c.rms[k],
        "size1": float(v.cov_axes(a.shape[i : i + 1])[0][0]),
        "size2": float(v.cov_axes(b.shape[j : j + 1])[0][0]),
        "size3": float(v.cov_axes(c.shape[k : k + 1])[0][0]),
    }


def run_search(cats: list[v.EpochCatalog], triples, model: v.ErrorModel) -> dict:
    res = v.search_multi(cats, triples, model=model)
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
        cands += [_cand_row(t, a, b, c, cd, n) for n in range(len(cd))]
        log(f"  {_tname(t)}: {out[_tname(t)]}")
    return {"per_triple": out, "candidates": cands}


def _legacy_cov(work: Path) -> list[np.ndarray]:
    """Runs 1-3's errors, for the comparison row: circularised catalogue error with the
    cos(dec) bug, plus the old 4-band floors (south of -40 on the all-sky floor)."""
    raw = catalogs(work, None)
    fl = v.calibrate_floors(raw, flux_min_mjy=10.0)["floor_arcsec"]
    bd = v.calibrate_floors_banded(raw, dec_edges=LEGACY_DEC_EDGES, flux_min_mjy=10.0)
    out = []
    for k, e in enumerate(EPOCHS):
        d = np.load(work / f"epoch{e}.npz")
        era = d["era"] * np.cos(np.radians(d["dec"]))
        err2 = 0.5 * (era**2 + d["edec"] ** 2) + v.banded_floor(d["dec"], k, bd, fl[k]) ** 2
        out.append(np.column_stack([err2, err2, np.zeros(err2.size)]))
    return out


def stage_calibrate(work: Path, cats: list[v.EpochCatalog]) -> dict:
    """Fit (k_struct, q_beam) on half the sky's matched statics; validate on the other half.

    Halves are alternate 1-degree RA strips (independent sky, same Dec and epoch mix). Tables
    are per bin of S/N, deconvolved size, beam change, declination and epoch pair, for three
    error models: runs 1-3 (legacy), measurement covariance only (Condon + floors), and the
    fitted model. Rayleigh: median 1.177, P(>3) = 0.0111, P(>5) = 3.7e-6.
    """
    legacy = _legacy_cov(work)
    parts, parts_legacy, pair_lab = [], [], []
    for a in range(len(cats)):
        for b in range(a + 1, len(cats)):
            m = v.match_statics(cats[a], cats[b])
            t = v.static_pair_terms(cats[a], cats[b], m)
            t["strip"] = (np.floor(cats[a].ra[m["i"]]) % 2).astype(float)
            parts.append(t)
            tl = dict(t)
            tl["cov"] = legacy[a][m["i"]] + legacy[b][m["j"]]
            tl["shape"] = np.zeros_like(t["shape"])
            tl["dbeam"] = np.zeros_like(t["dbeam"])
            parts_legacy.append(tl)
            pair_lab.append(np.full(m["i"].size, 10 * (a + 1) + (b + 1), float))
    terms = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    terms["pair"] = np.concatenate(pair_lab)
    tleg = {k: np.concatenate([p[k] for p in parts_legacy]) for k in parts_legacy[0]}
    tleg["pair"] = terms["pair"]
    fit_sel = terms["strip"] == 0
    model, surface = v.fit_error_model(v._subset_terms(terms, fit_sel))
    val = v._subset_terms(terms, ~fit_sel)
    val_leg = v._subset_terms(tleg, ~fit_sel)
    bins = {
        "snr": v.SNR_BIN_EDGES,
        "size": v.SIZE_BIN_EDGES,
        "beam_change": v.BEAM_BIN_EDGES,
        "dec": (-90.0, -35.0, -30.0, -25.0, -20.0, -10.0, 0.0, 15.0, 30.0, 50.0, 90.0),
        "pair": (12, 13, 14, 23, 24, 34, 35),
    }
    tables: dict = {}
    for name, tt, mdl in (
        ("legacy_runs1to3", val_leg, v.NO_MODEL),
        ("measurement_only", val, v.NO_MODEL),
        ("fitted_model", val, model),
    ):
        tables[name] = {
            "all": v.calibration_table(tt, mdl, by="snr", edges=(0.0, np.inf))[0],
            **{by: v.calibration_table(tt, mdl, by=by, edges=ed) for by, ed in bins.items()},
        }
        log(f"  {name}: {tables[name]['all']}")
    # Joint S/N x size table for the fitted model (validation half): where does it fail?
    joint = []
    for lo, hi in zip(v.SNR_BIN_EDGES[:-1], v.SNR_BIN_EDGES[1:], strict=True):
        sel = (val["snr"] >= lo) & (val["snr"] < hi)
        joint.append(
            {
                "snr_lo": lo,
                "snr_hi": hi,
                "by_size": v.calibration_table(
                    v._subset_terms(val, sel), model, by="size", edges=v.SIZE_BIN_EDGES
                ),
            }
        )
    # Catalogue check: the recomputed Condon ellipse vs the quoted E_RA/E_DEC.
    check = {}
    for e in EPOCHS:
        d = np.load(work / f"epoch{e}.npz")
        rms = d["rms"] * PEAK_CORRECTION[e]
        cv = v.condon_cov(d["maj"], d["min"], d["pa"], d["bmaj"], d["bmin"], d["flux"] / rms)
        check[f"E{e}"] = {
            "median_ratio_ra": float(np.median(np.sqrt(cv[:, 0]) / d["era"])),
            "median_ratio_dec": float(np.median(np.sqrt(cv[:, 1]) / d["edec"])),
            "p05_p95_ratio_ra": [
                float(x) for x in np.percentile(np.sqrt(cv[:, 0]) / d["era"], [5, 95])
            ],
            "median_catalogue_err_arcsec": float(np.median(np.hypot(d["era"], d["edec"]) / 2**0.5)),
            "frac_deconvolved_zero": float(np.mean(d["dcmaj"] == 0)),
            "median_beam_arcsec": [float(np.median(d["bmaj"])), float(np.median(d["bmin"]))],
        }
    return {
        "method": (
            'statics = isolated (no other component within 30") components matched within 5" '
            "across every epoch pair; fit on even 1-deg RA strips, validated on odd ones; "
            "loss = sum over (size x S/N) and (size x beam change) bins of log(P(z>3)/0.0111)^2"
        ),
        "n_pairs_fit": int(fit_sel.sum()),
        "n_pairs_validate": int((~fit_sel).sum()),
        "model": {"k_struct": model.k_struct, "q_beam": model.q_beam},
        "loss_surface": surface,
        "validation": tables,
        "validation_snr_x_size": joint,
        "catalogue_condon_check": check,
    }


def rescore_previous(cats: list[v.EpochCatalog], model: v.ErrorModel) -> list[dict]:
    """The run-3 candidates (results/vlasspm_candidates.csv, by catalogue row) under the old
    measurement-only errors and the new model: displacement significances and E3 residual."""
    path = REPO / "results" / "vlasspm_candidates_run3.csv"
    if not path.exists():
        return []
    rows = list(csv.DictReader(path.open()))
    idx = [{int(r): n for n, r in enumerate(c.ident)} for c in cats]
    out = []
    for n, r in enumerate(rows):
        ep = [int(x[1]) - 1 for x in r["triple"].split("-")]
        a, b, c = (cats[e] for e in ep)
        i, j, k = (idx[e][int(r[f"row{m + 1}"])] for m, e in enumerate(ep))
        rec = {"c": f"c{n}", "triple": r["triple"], "dec": float(r["dec1"])}
        rec["run3_resid_sigma"] = float(r["resid_sigma"])
        rec["measurement_only"] = v.triplet_significance(a, i, b, j, c, k, v.NO_MODEL)
        rec["fitted_model"] = v.triplet_significance(a, i, b, j, c, k, model)
        rec["snr"] = [float(x.flux[y] / x.rms[y]) for x, y in ((a, i), (b, j), (c, k))]
        rec["size_arcsec"] = [
            float(v.cov_axes(x.shape[y : y + 1])[0][0]) for x, y in ((a, i), (b, j), (c, k))
        ]
        out.append(rec)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=str(REPO / "data/vlass/work"))
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
        log("floors: per declination band (fine southern bands) ...")
        _save(
            f,
            v.calibrate_floors_banded(
                catalogs(work, None), dec_edges=DEC_EDGES_FINE, flux_min_mjy=10.0
            ),
        )
    banded = json.loads(f.read_text())
    for b in banded["bands"]:
        fl = b["floor_arcsec"]
        log(
            f"  Dec [{b['lo']:+.0f},{b['hi']:+.0f}): "
            + (str([round(x, 3) for x in fl]) if fl else "too few")
        )
    cats = catalogs(work, floors["floor_arcsec"], banded)
    triples = v.epoch_triples(len(cats))

    f = work / "calibration.json"
    if not f.exists():
        log("calibration: shape-aware error model from matched static sources ...")
        _save(f, stage_calibrate(work, cats))
    calib = json.loads(f.read_text())
    model = v.ErrorModel(calib["model"]["k_struct"], calib["model"]["q_beam"])
    log(f"  model: {model}")

    f = work / "search.json"
    if not f.exists():
        log("search: every epoch triple (shape-aware errors) ...")
        _save(f, run_search(cats, triples, model))
    search = json.loads(f.read_text())

    f = work / "search_ablation.json"
    if not f.exists():
        log("ablation: search with the measurement covariance only (no k/q systematics) ...")
        abl = run_search(cats, triples, v.NO_MODEL)
        _save(
            f,
            {
                "per_triple": abl["per_triple"],
                "n_candidates": len(abl["candidates"]),
                "candidates": abl["candidates"],
            },
        )
        log(f"  {json.loads(f.read_text())['n_candidates']} candidates without the systematics")

    for t in triples:
        f = work / f"null_{_tname(t)}.json"
        if f.exists():
            continue
        log(f"null {_tname(t)}: {args.n_null} scrambles ...")
        r = v.search(cats[t[0]], cats[t[1]], cats[t[2]], model=model)
        _save(
            f,
            v.scramble_null(r.orphans, n_reps=args.n_null, seed=hash(t) % 2**31, model=model),
        )
        log(
            f"  {_tname(t)} chance candidates/scramble: "
            f"{json.loads(f.read_text())['candidates_mean']:.2f}"
        )

    for flux in (3.0, 1.5):
        f = work / f"completeness_{flux:g}mJy.json"
        if f.exists():
            continue
        log(f"completeness: {args.n_inject} injected movers at {flux} mJy (half realistic) ...")
        _save(
            f,
            v.completeness(
                *cats,
                n=args.n_inject,
                flux_mjy=flux,
                seed=int(flux * 10),
                realistic_frac=0.5,
                model=model,
            ),
        )
        c = json.loads(f.read_text())
        log(f"  isolated {c['isolated']['overall']:.3f}  realistic {c['realistic']['overall']:.3f}")

    # A faint point source does not come out of PyBDSF with size 0: its deconvolved size is
    # noisy (UV Ceti's is 1.2", 0 and 2.9" in E2-E4; the median faint component's is ~1.7").
    # Injecting at size 0 therefore gives the injections the smallest structure term the model
    # allows; this variant brackets that choice with a 2" deconvolved size.
    f = work / "completeness_3mJy_size2.json"
    if not f.exists():
        log(f"completeness: {args.n_inject} injected movers at 3 mJy, 2 arcsec size ...")
        _save(
            f,
            v.completeness(
                *cats,
                n=args.n_inject,
                flux_mjy=3.0,
                seed=31,
                realistic_frac=0.5,
                size_arcsec=2.0,
                model=model,
            ),
        )
        c = json.loads(f.read_text())
        log(f"  isolated {c['isolated']['overall']:.3f}  realistic {c['realistic']['overall']:.3f}")

    f = work / "uvcet.json"
    if not f.exists():
        log("uvcet: recover-a-known ...")
        _save(f, uvcet_check(search["candidates"], cats))
    log(f"  {json.loads(f.read_text())}")

    f = work / "previous_candidates.json"
    if not f.exists():
        log("re-scoring the run-3 candidates under the new errors ...")
        _save(f, rescore_previous(cats, model))
    for r in json.loads(f.read_text()):
        log(f"  {r}")

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
                    "resid_sigma": c["resid_sigma"],
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
    calib = json.loads((work / "calibration.json").read_text())
    calib.pop("loss_surface", None)  # large; kept in the work directory
    abl = json.loads((work / "search_ablation.json").read_text())
    prev = json.loads((work / "previous_candidates.json").read_text())
    dark = [
        c | {"vet": vt}
        for c, vt in zip(search["candidates"], vet, strict=True)
        if vt.get("gaia") == 0 and vt.get("catwise") == 0
    ]
    area = max(pt["area_deg2"] for pt in search["per_triple"].values())
    res_path = out / "results" / "vlasspm_metrics.json"
    previous: dict = {}
    if res_path.exists():
        old = json.loads(res_path.read_text())
        previous = old.get("real_previous_run3") or {
            k: old.get(k)
            for k in (
                "real_floors_banded",
                "real_per_triple",
                "real_null",
                "real_completeness",
                "real_n_candidates",
                "real_uvcet",
            )
        }
    metrics = {
        "source": "real: VLASS Quick-Look epochs 1, 2, 3 (QL3.1+3.2) and 4.1 component catalogues",
        "is_real": True,
        "real_error_model": (
            "run 4 (shape-aware): Condon (1997) ellipse per component from fitted Maj/Min/PA, "
            "beam and peak S/N (E_RA is on-sky; runs 1-3 wrongly multiplied it by cos dec), "
            "plus bright-static floors in fine Dec bands, plus k_struct^2 x deconvolved shape "
            "and q_beam^2 x |beam change|, both fitted on matched statics"
        ),
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
        "real_calibration": calib,
        "real_per_triple": search["per_triple"],
        "real_ablation_measurement_only": {
            "per_triple": abl["per_triple"],
            "n_candidates": abl["n_candidates"],
        },
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
        "real_run3_candidates_rescored": prev,
        "real_area_deg2_max_triple": area,
        "real_previous_run3": previous,
    }
    (out / "results").mkdir(parents=True, exist_ok=True)
    written = write_results(metrics, res_path)
    # Superseded run-3 key: it treated unvetted candidates as detections (see the findings).
    if "real_density_limit_per_deg2_3mJy" in written:
        written.pop("real_density_limit_per_deg2_3mJy")
        m = written.get("_merge", {}).get("retained_from_previous_run", [])
        if "real_density_limit_per_deg2_3mJy" in m:
            m.remove("real_density_limit_per_deg2_3mJy")
        if "_merge" in written and not m:
            written.pop("_merge")
        res_path.write_text(json.dumps(written, indent=2) + "\n")
    keys = list(search["candidates"][0]) if search["candidates"] else list(_CAND_KEYS)
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


_CAND_KEYS = (
    "triple",
    "ra1",
    "dec1",
    "t1",
    "row1",
    "ra2",
    "dec2",
    "t2",
    "row2",
    "ra3",
    "dec3",
    "t3",
    "row3",
    "mu_ra",
    "mu_dec",
    "mu",
    "resid_sigma",
)


if __name__ == "__main__":
    raise SystemExit(main())
