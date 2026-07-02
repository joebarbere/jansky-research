#!/usr/bin/env python
"""RACS Stokes-V discovery — real leg driver (plan 33, jansky-research).

Forced I+V photometry of the nearest CNS5 M dwarfs in the two RACS-mid epochs
(mid1 ~2021-01, mid2 ~2025), CASDA SODA cutouts, one CSV row per (target, epoch).

Run from the jansky-research checkout:

    cd /home/joe/dev/github/joebarbere/jansky-research
    uv run python /tmp/svd_real_leg.py [--limit N]

Reuses jansky_research.stokesv (_casda_session getpass shim reading ~/.casda_pw,
measure_circular_pol, the retry-with-relogin pattern of fetch_racs_cutout) and
jansky_research.stokesv_discovery (epoch_position, RACS_MID*_MJD nominal anchors).
Resumable: (name, epoch) rows already in the CSV are skipped; flushed per row.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np

REPO = Path("/home/joe/dev/github/joebarbere/jansky-research")
sys.path.insert(0, str(REPO / "src"))

from jansky_research.stokesv import _casda_session, measure_circular_pol  # noqa: E402
from jansky_research.stokesv_discovery import (  # noqa: E402
    RACS_MID1_MJD,
    RACS_MID2_MJD,
    epoch_position,
)

CASDA_TAP = "https://casda.csiro.au/casda_vo_tools/tap"
CNS5_VIZIER = "J/A+A/670/A19"
CSV_PATH = REPO / "results" / "stokesv_discovery_realtargets.csv"
CSV_FIELDS = [
    "name",
    "gaia_id",
    "epoch",
    "ra_prop",
    "dec_prop",
    "epoch_mjd",
    "band_mhz",
    "i_mjy",
    "e_i",
    "v_mjy",
    "e_v",
    "n_products",
    "note",
]
C_M_S = 299792458.0
CUTOUT_RADIUS_DEG = 0.025  # 1.5 arcmin radius -> ~3 arcmin cutout
MIN_G_RP = 0.85  # M-dwarf red-colour cut (no SpT column in the VizieR CNS5 table)
EPOCH_SPLIT_MJD = 60000.0  # mid1 (~59200-59500) vs mid2 (~60600-60800)


def fetch_cns5_mdwarfs(n: int) -> list[dict]:
    """Nearest ~n CNS5 M dwarfs (Dec<+40, red colour, Gaia PM), parallax-descending."""
    from astroquery.vizier import Vizier

    v = Vizier(
        columns=[
            "CNS5",
            "GJ",
            "Comp",
            "GaiaDR3",
            "RAJ2000",
            "DEJ2000",
            "Epoch",
            "plx",
            "pmRA",
            "pmDE",
            "Gmag",
            "RPmag",
            "SimbadName",
        ],
        row_limit=-1,
    )
    t = v.get_catalogs(CNS5_VIZIER)[0]

    def col(name):
        return np.ma.filled(np.asarray(t[name], float), np.nan)

    ra, dec = col("RAJ2000"), col("DEJ2000")
    epoch = np.where(np.isfinite(col("Epoch")), col("Epoch"), 2016.0)
    plx, pmra, pmde = col("plx"), col("pmRA"), col("pmDE")
    g_rp = col("Gmag") - col("RPmag")
    sel = (
        np.isfinite(ra)
        & np.isfinite(dec)
        & (dec < 40.0)
        & np.isfinite(plx)
        & (plx > 0)
        & np.isfinite(pmra)
        & np.isfinite(pmde)
        & np.isfinite(g_rp)
        & (g_rp >= MIN_G_RP)
    )
    idx = np.where(sel)[0]
    idx = idx[np.argsort(-plx[idx])][:n]
    targets = []
    for i in idx:
        simbad = str(t["SimbadName"][i]).strip()
        cns5 = str(t["CNS5"][i]).strip()
        gaia = str(t["GaiaDR3"][i]).strip()
        targets.append(
            {
                "name": simbad or cns5 or f"GaiaDR3_{gaia}",
                "gaia_id": gaia if gaia and gaia != "--" else "",
                "ra": float(ra[i]),
                "dec": float(dec[i]),
                "epoch_yr": float(epoch[i]),
                "plx": float(plx[i]),
                "pmra": float(pmra[i]),
                "pmde": float(pmde[i]),
            }
        )
    return targets


def obscore_products(tap, ra: float, dec: float):
    """RACS-mid restored taylor.0 I+V products whose s_region contains (ra, dec)."""
    q = f"""
    SELECT obs_id, filename, t_min, em_min, em_max, access_url
    FROM ivoa.obscore
    WHERE obs_collection = 'The Rapid ASKAP Continuum Survey'
    AND (filename LIKE 'image.i.%.restored.%' OR filename LIKE 'image.v.%.restored.%')
    AND filename LIKE '%.taylor.0.restored.conv.fits'
    AND em_min > 0.21 AND em_max < 0.23
    AND 1 = CONTAINS(POINT('ICRS', {ra:.8f}, {dec:.8f}), s_region)
    """
    return tap.search(q).to_table()


def pick_epoch_pair(table, side: str):
    """From an obscore table, pick the (I-row, V-row, n_products) for epoch mid1/mid2.

    mid1 = earliest obs_id group with t_min < EPOCH_SPLIT_MJD; mid2 = latest group above it.
    Requires both an I and a V taylor.0 product in the same obs_id (same observation).
    """
    if len(table) == 0:
        return None, None, 0
    t_min = np.asarray(table["t_min"], float)
    in_side = t_min < EPOCH_SPLIT_MJD if side == "mid1" else t_min >= EPOCH_SPLIT_MJD
    n_products = int(in_side.sum())
    if not in_side.any():
        return None, None, 0
    sub = table[in_side]
    groups: dict[str, dict[str, int]] = {}
    for k, row in enumerate(sub):
        stokes = str(row["filename"]).split(".")[1]  # image.<stokes>.RACS_...
        groups.setdefault(str(row["obs_id"]), {})[stokes] = k
    complete = {o: g for o, g in groups.items() if "i" in g and "v" in g}
    if not complete:
        return None, None, n_products
    key = (min if side == "mid1" else max)(
        complete, key=lambda o: float(sub["t_min"][complete[o]["i"]])
    )
    return sub[complete[key]["i"]], sub[complete[key]["v"]], n_products


def fetch_cutout_image(casda, row, ra: float, dec: float):
    """SODA-stage one obscore product row, cut out at (ra, dec) -> (image_mJy, wcs).

    Mirrors stokesv.fetch_racs_cutout's download/read, but the product row comes from
    the epoch-resolved TAP query instead of Casda.query_region.
    """
    import astropy.units as u
    import requests
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.table import Table
    from astropy.wcs import WCS

    one = Table(rows=[(str(row["access_url"]),)], names=("access_url",))
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    urls = casda.cutout(one, coordinates=coord, radius=CUTOUT_RADIUS_DEG * u.deg)
    furl = next(u_ for u_ in urls if u_.endswith(".fits"))
    raw = requests.get(furl, timeout=200).content
    with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as fh:
        fh.write(raw)
        path = fh.name
    try:
        with fits.open(path) as hd:
            data = np.squeeze(np.asarray(hd[0].data, float))
            wcs = WCS(hd[0].header).celestial
    finally:
        os.unlink(path)
    return data * 1000.0, wcs  # Jy/beam -> mJy/beam


def measure_epoch(casda, username, pw_path, irow, vrow, ra, dec, retries=3):
    """I+V cutouts + forced photometry, retry-with-relogin (CASDA 401s intermittently)."""
    last = None
    for attempt in range(retries):
        try:
            if casda is None:
                casda = _casda_session(username, pw_path)
            img_i, wcs = fetch_cutout_image(casda, irow, ra, dec)
            img_v, _ = fetch_cutout_image(casda, vrow, ra, dec)
            meas = measure_circular_pol(img_i, img_v, wcs, ra, dec)
            return meas, casda
        except Exception as exc:  # noqa: BLE001
            last = exc
            print(f"    retry {attempt + 1}/{retries} after: {exc!r}", flush=True)
            casda = None  # force fresh login (handles the intermittent 401)
    raise RuntimeError(f"cutout/photometry failed after {retries} retries: {last!r}")


def load_done(path: Path) -> set[tuple[str, str]]:
    if not path.exists():
        return set()
    with path.open() as fh:
        return {(r["name"], r["epoch"]) for r in csv.DictReader(fh)}


def main() -> int:
    ap = argparse.ArgumentParser(description="RACS Stokes-V discovery real leg (plan 33)")
    ap.add_argument("--limit", type=int, default=60, help="number of CNS5 targets")
    ap.add_argument("--csv", default=str(CSV_PATH))
    args = ap.parse_args()

    import pyvo

    username = os.environ.get("CASDA_USERNAME", "joe.barbere@gmail.com")
    pw_path = "~/.casda_pw"
    out = Path(args.csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    done = load_done(out)
    print(f"resuming: {len(done)} (target, epoch) rows already in {out}", flush=True)

    print("fetching CNS5 targets from VizieR ...", flush=True)
    targets = fetch_cns5_mdwarfs(args.limit)
    print(
        f"{len(targets)} CNS5 M-dwarf targets (Dec<+40, G-RP>={MIN_G_RP}, "
        f"parallax {targets[0]['plx']:.0f} .. {targets[-1]['plx']:.0f} mas)",
        flush=True,
    )

    tap = pyvo.dal.TAPService(CASDA_TAP)
    casda = None
    new_file = not out.exists() or out.stat().st_size == 0
    fh = out.open("a", newline="")
    writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
    if new_file:
        writer.writeheader()
        fh.flush()

    def emit(rowdict):
        writer.writerow(rowdict)
        fh.flush()
        os.fsync(fh.fileno())

    nominal = {"mid1": RACS_MID1_MJD, "mid2": RACS_MID2_MJD}
    t0 = time.time()
    for it, tgt in enumerate(targets):
        pending = [e for e in nominal if (tgt["name"], e) not in done]
        if not pending:
            print(f"[{it + 1}/{len(targets)}] {tgt['name']}: already done, skipping", flush=True)
            continue
        print(
            f"[{it + 1}/{len(targets)}] {tgt['name']} (plx {tgt['plx']:.0f} mas, "
            f"pm {tgt['pmra']:.0f},{tgt['pmde']:.0f} mas/yr)",
            flush=True,
        )
        for epoch in pending:
            base = {
                "name": tgt["name"],
                "gaia_id": tgt["gaia_id"],
                "epoch": epoch,
                "ra_prop": "nan",
                "dec_prop": "nan",
                "epoch_mjd": "nan",
                "band_mhz": "nan",
                "i_mjy": "nan",
                "e_i": "nan",
                "v_mjy": "nan",
                "e_v": "nan",
                "n_products": 0,
                "note": "",
            }
            try:
                # pass 1: propagate with the nominal module MJD (fine for containment)
                ra0, dec0 = epoch_position(
                    tgt["ra"],
                    tgt["dec"],
                    tgt["pmra"],
                    tgt["pmde"],
                    tgt["epoch_yr"],
                    nominal[epoch],
                )
                prods = obscore_products(tap, float(ra0), float(dec0))
                irow, vrow, n_prod = pick_epoch_pair(prods, epoch)
                base["n_products"] = n_prod
                if irow is None:
                    base["note"] = f"no {epoch} I+V product pair in obscore"
                    print(f"  {epoch}: {base['note']}", flush=True)
                    emit(base)
                    continue
                # pass 2: re-propagate to the product's actual t_min
                mjd = float(irow["t_min"])
                ra_p, dec_p = epoch_position(
                    tgt["ra"], tgt["dec"], tgt["pmra"], tgt["pmde"], tgt["epoch_yr"], mjd
                )
                ra_p, dec_p = float(ra_p), float(dec_p)
                lam = 0.5 * (float(irow["em_min"]) + float(irow["em_max"]))
                band_mhz = C_M_S / lam / 1e6
                meas, casda = measure_epoch(casda, username, pw_path, irow, vrow, ra_p, dec_p)
                base.update(
                    ra_prop=f"{ra_p:.6f}",
                    dec_prop=f"{dec_p:.6f}",
                    epoch_mjd=f"{mjd:.5f}",
                    band_mhz=f"{band_mhz:.1f}",
                    i_mjy=f"{meas['i_peak']:.4f}",
                    e_i=f"{meas['i_rms']:.4f}",
                    v_mjy=f"{meas['v_peak']:.4f}",
                    e_v=f"{meas['v_rms']:.4f}",
                    note=f"offset_arcsec={meas['offset_arcsec']:.2f}",
                )
                print(
                    f"  {epoch}: MJD {mjd:.1f}  I={meas['i_peak']:.3f}"
                    f"+/-{meas['i_rms']:.3f} mJy  V={meas['v_peak']:.3f}"
                    f"+/-{meas['v_rms']:.3f} mJy  ({n_prod} products)",
                    flush=True,
                )
                emit(base)
            except Exception as exc:  # noqa: BLE001 - one target must not kill the run
                base["note"] = f"FAILED: {type(exc).__name__}: {exc}"[:300]
                print(f"  {epoch}: {base['note']}", flush=True)
                traceback.print_exc()
                emit(base)
                casda = None
    fh.close()
    dt = time.time() - t0
    print(f"done in {dt / 60:.1f} min ({dt:.0f} s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
