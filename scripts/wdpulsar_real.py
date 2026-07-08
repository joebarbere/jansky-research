#!/usr/bin/env python
"""WD-pulsar candidate radio survey — real leg driver (plan 41, jansky-research).

Forced I+V photometry of the 56 Pelisoli+2025 candidates (+ the AR Sco control) in EVERY
complete RACS (low or mid) I+V observation covering each position: CASDA obscore per-position
product query, SODA cutouts, one CSV row per (target, obs_id). Adapted from the plan-33 driver
(`stokesv_discovery_real.py`): same session/retry/resume mechanics, but all epochs per target
(not a fixed pair) and both RACS bands.

Run from the jansky-research checkout:

    uv run --extra vlass python scripts/wdpulsar_real.py [--limit N]

Resumable: (name, obs_id) rows already in results/wdpulsar_realtargets.csv are skipped;
flushed per row.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research.stokesv import _casda_session, measure_circular_pol  # noqa: E402
from jansky_research.wdpulsar import candidate_targets, load_candidate_table  # noqa: E402

CASDA_TAP = "https://casda.csiro.au/casda_vo_tools/tap"
CSV_PATH = REPO / "results" / "wdpulsar_realtargets.csv"
CSV_FIELDS = [
    "name",
    "type",
    "obs_id",
    "band",
    "epoch_mjd",
    "i_mjy",
    "e_i",
    "v_mjy",
    "e_v",
    "offset_arcsec",
    "note",
]
CUTOUT_RADIUS_DEG = 0.025


def obscore_products(tap, ra: float, dec: float):
    """All RACS restored taylor.0 I+V products (low + mid bands) containing (ra, dec)."""
    q = f"""
    SELECT obs_id, filename, t_min, em_min, em_max, access_url
    FROM ivoa.obscore
    WHERE obs_collection = 'The Rapid ASKAP Continuum Survey'
    AND (filename LIKE 'image.i.%.restored.%' OR filename LIKE 'image.v.%.restored.%')
    AND filename LIKE '%.taylor.0.restored.conv.fits'
    AND 1 = CONTAINS(POINT('ICRS', {ra:.8f}, {dec:.8f}), s_region)
    """
    return tap.search(q).to_table()


def complete_iv_groups(table) -> list[dict]:
    """Group products by obs_id; keep groups with BOTH an I and a V taylor.0 product."""
    groups: dict[str, dict] = {}
    for k in range(len(table)):
        filename = str(table["filename"][k])
        stokes = filename.split(".")[1]
        obs = str(table["obs_id"][k])
        g = groups.setdefault(
            obs,
            {
                "obs_id": obs,
                "t_min": float(table["t_min"][k]),
                # band from wavelength: RACS-mid ~0.21-0.23 m, RACS-low ~0.30-0.36 m
                "band": "mid" if float(table["em_min"][k]) < 0.26 else "low",
            },
        )
        g[stokes] = k
    return [g for g in groups.values() if "i" in g and "v" in g]


def fetch_cutout_image(casda, table, row_idx: int, ra: float, dec: float):
    """SODA-stage one obscore product row, cut out at (ra, dec) -> (image_mJy, wcs)."""
    import astropy.units as u
    import requests
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.table import Table
    from astropy.wcs import WCS

    one = Table(rows=[(str(table["access_url"][row_idx]),)], names=("access_url",))
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


def measure_group(casda, username, pw_path, table, group, ra, dec, retries=3):
    """I+V cutouts + forced photometry for one obs group, retry-with-relogin."""
    last = None
    for attempt in range(retries):
        try:
            if casda is None:
                casda = _casda_session(username, pw_path)
            img_i, wcs = fetch_cutout_image(casda, table, group["i"], ra, dec)
            img_v, _ = fetch_cutout_image(casda, table, group["v"], ra, dec)
            # SODA edge rounding occasionally returns I/V cutouts differing by one pixel;
            # crop both to the common overlap (the WCS stays valid for the kept corner)
            ny = min(img_i.shape[0], img_v.shape[0])
            nx = min(img_i.shape[1], img_v.shape[1])
            img_i, img_v = img_i[:ny, :nx], img_v[:ny, :nx]
            return measure_circular_pol(img_i, img_v, wcs, ra, dec), casda
        except Exception as exc:  # noqa: BLE001
            last = exc
            print(f"    retry {attempt + 1}/{retries} after: {exc!r}", flush=True)
            casda = None
    raise RuntimeError(f"cutout/photometry failed after {retries} retries: {last!r}")


def load_done(path: Path) -> set[tuple[str, str]]:
    """Rows already measured. Failed rows do NOT count: they are retried on the next run."""
    if not path.exists():
        return set()
    with path.open() as fh:
        return {
            (r["name"], r["obs_id"])
            for r in csv.DictReader(fh)
            if not r.get("note", "").startswith("failed")
        }


def main() -> int:
    ap = argparse.ArgumentParser(description="WD-pulsar candidate radio survey (plan 41)")
    ap.add_argument("--limit", type=int, default=0, help="limit targets (0 = all 57)")
    ap.add_argument("--csv", default=str(CSV_PATH))
    args = ap.parse_args()

    import pyvo

    username = os.environ.get("CASDA_USERNAME", "joe.barbere@gmail.com")
    pw_path = "~/.casda_pw"
    out = Path(args.csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    done = load_done(out)
    print(f"resuming: {len(done)} (target, obs_id) rows already in {out}", flush=True)

    targets = candidate_targets(load_candidate_table())
    if args.limit:
        targets = targets[: args.limit]
    print(f"{len(targets)} targets (56 candidates + AR Sco control)", flush=True)

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

    t0 = time.time()
    for it, tgt in enumerate(targets):
        print(f"[{it + 1}/{len(targets)}] {tgt['name']} ({tgt['type']})", flush=True)
        try:
            table = obscore_products(tap, tgt["ra_deg"], tgt["dec_deg"])
            groups = complete_iv_groups(table)
        except Exception as exc:  # noqa: BLE001
            print(f"    obscore query failed: {exc!r} (will retry next run)", flush=True)
            continue
        pending = [g for g in groups if (tgt["name"], g["obs_id"]) not in done]
        print(f"    {len(groups)} complete I+V groups, {len(pending)} pending", flush=True)
        for g in pending:
            base = {
                "name": tgt["name"],
                "type": tgt["type"],
                "obs_id": g["obs_id"],
                "band": g["band"],
                "epoch_mjd": f"{g['t_min']:.1f}",
                "i_mjy": "nan",
                "e_i": "nan",
                "v_mjy": "nan",
                "e_v": "nan",
                "offset_arcsec": "nan",
                "note": "",
            }
            try:
                meas, casda = measure_group(
                    casda, username, pw_path, table, g, tgt["ra_deg"], tgt["dec_deg"]
                )
                base.update(
                    {
                        "i_mjy": f"{meas['i_peak']:.4f}",
                        "e_i": f"{meas['i_rms']:.4f}",
                        "v_mjy": f"{meas['v_peak']:.4f}",
                        "e_v": f"{meas['v_rms']:.4f}",
                        "offset_arcsec": f"{meas['offset_arcsec']:.2f}",
                    }
                )
            except Exception as exc:  # noqa: BLE001
                base["note"] = f"failed: {type(exc).__name__}"
                print(f"    {g['obs_id']}: FAILED {exc!r}", flush=True)
            emit(base)
        elapsed = (time.time() - t0) / 60.0
        print(f"    elapsed {elapsed:.1f} min", flush=True)
    fh.close()
    print("sweep complete", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
