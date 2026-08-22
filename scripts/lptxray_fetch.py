#!/usr/bin/env python3
"""Plan 93: stage X-ray catalogue cones around every LPT and WD-pulsar candidate position.

NETWORK. One HEASARC cone per (position, catalogue) at ``--cone-arcmin``, cached to
``data/lptxray/cones.json``. Everything downstream -- association, the measured
chance-coincidence rate, the local source density -- is computed offline from these cones,
so the rigid position-shift trials cost no further queries.

    uv run python scripts/lptxray_fetch.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

CATALOGS = ("xmmssc", "rass2rxs", "erass1main", "csc")
KEEP_MAX = 400


def _rows(table) -> list[dict]:  # type: ignore[no-untyped-def]
    out = []
    for row in table[:KEEP_MAX]:
        rec = {}
        for name in table.colnames:
            val = row[name]
            try:
                if hasattr(val, "mask") and val.mask:
                    continue
                rec[name] = val.item() if hasattr(val, "item") else str(val)
            except Exception:
                rec[name] = str(val)
        out.append(rec)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cone-arcmin", type=float, default=10.0)
    ap.add_argument("--out", type=Path, default=Path("data/lptxray/cones.json"))
    args = ap.parse_args(argv)

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.heasarc import Heasarc

    targets: list[dict] = []
    for r in csv.DictReader(open("data/lpt_sample.csv")):
        targets.append(
            {
                "name": r["name"],
                "sample": "lpt",
                "ra": float(r["ra_deg"]),
                "dec": float(r["dec_deg"]),
                "type": "LPT",
                "lit_xray": r["xray"],
            }
        )
    for r in csv.DictReader(open("data/wdpulsar_candidates.csv")):
        targets.append(
            {
                "name": r["short_name"],
                "sample": "wdcand",
                "ra": float(r["ra_deg"]),
                "dec": float(r["dec_deg"]),
                "type": r["type"],
                "lit_xray": r["xray"],
                "gmag": r["gmag"],
                "class_this_work": r["class_this_work"],
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cache: dict = {}
    if args.out.exists():
        cache = json.loads(args.out.read_text())
    cones = cache.setdefault("cones", {})

    radius = args.cone_arcmin * u.arcmin
    total = len(targets) * len(CATALOGS)
    done = 0
    for t in targets:
        key = f"{t['sample']}:{t['name']}"
        entry = cones.setdefault(key, {"meta": t, "cats": {}})
        entry["meta"] = t
        coord = SkyCoord(t["ra"] * u.deg, t["dec"] * u.deg)
        for cat in CATALOGS:
            done += 1
            if cat in entry["cats"]:
                continue
            for attempt in range(3):
                try:
                    tab = Heasarc.query_region(coord, catalog=cat, radius=radius)
                    entry["cats"][cat] = {"n": len(tab), "rows": _rows(tab)}
                    break
                except Exception as exc:  # noqa: BLE001
                    if attempt == 2:
                        entry["cats"][cat] = {"error": str(exc)[:200]}
                    else:
                        time.sleep(3 * (attempt + 1))
            print(
                f"[{done:3d}/{total}] {key:40s} {cat:11s} -> {entry['cats'][cat].get('n', 'ERR')}",
                flush=True,
            )
            tmp = args.out.with_suffix(".json.part")
            tmp.write_text(json.dumps({"cone_arcmin": args.cone_arcmin, "cones": cones}))
            tmp.replace(args.out)

    print(f"\nwrote {args.out}: {len(cones)} positions x {len(CATALOGS)} catalogues")
    return 0


if __name__ == "__main__":
    sys.exit(main())
