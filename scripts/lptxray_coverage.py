#!/usr/bin/env python3
"""Plan 93: pointed X-ray observation coverage at every LPT position.

NETWORK. The serendipitous source catalogues are blind to the dedicated follow-up that
produced every published LPT X-ray detection (see ``jansky_research.lptxray``), so the
question "has this position been observed at all?" has to be asked of the observation logs
rather than the source catalogues. Caches to ``data/lptxray/coverage.json``.

    uv run python scripts/lptxray_coverage.py
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

# Field-of-view scale per mission: the radius within which a target could plausibly land on
# a detector for a pointing at the listed aimpoint. Generous on purpose -- this is a
# "was it ever looked at" census, and a near-miss is worth reporting as a near-miss.
LOGS = {"xmmmaster": 15.0, "chanmaster": 15.0, "erosmaster": 30.0}
KEEP = ("obsid", "name", "time", "exposure", "status", "public_date", "pi_lname", "ra", "dec")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("data/lptxray/coverage.json"))
    args = ap.parse_args(argv)

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.heasarc import Heasarc

    targets = [
        {
            "name": r["name"],
            "ra": float(r["ra_deg"]),
            "dec": float(r["dec_deg"]),
            "lit_xray": r["xray"],
        }
        for r in csv.DictReader(open("data/lpt_sample.csv"))
    ]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out: dict = json.loads(args.out.read_text()) if args.out.exists() else {}
    cov = out.setdefault("coverage", {})

    for t in targets:
        entry = cov.setdefault(t["name"], {"meta": t, "logs": {}})
        coord = SkyCoord(t["ra"] * u.deg, t["dec"] * u.deg)
        for log, arcmin in LOGS.items():
            if log in entry["logs"]:
                continue
            recs: list[dict] = []
            err = None
            for attempt in range(3):
                try:
                    tab = Heasarc.query_region(coord, catalog=log, radius=arcmin * u.arcmin)
                    for row in tab[:60]:
                        rec = {}
                        for k in KEEP:
                            if k in tab.colnames:
                                v = row[k]
                                try:
                                    rec[k] = v.item() if hasattr(v, "item") else str(v)
                                except Exception:
                                    rec[k] = str(v)
                        recs.append(rec)
                    err = None
                    break
                except Exception as exc:  # noqa: BLE001
                    err = str(exc)[:200]
                    time.sleep(3 * (attempt + 1))
            entry["logs"][log] = {"n": len(recs), "obs": recs} if err is None else {"error": err}
            print(
                f"{t['name'][:32]:32s} {log:11s} -> {entry['logs'][log].get('n', 'ERR')}",
                flush=True,
            )
            tmp = args.out.with_suffix(".json.part")
            tmp.write_text(json.dumps(out, indent=1))
            tmp.replace(args.out)

    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
