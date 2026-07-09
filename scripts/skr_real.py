#!/usr/bin/env python
"""Cassini SKR real-leg driver (plan 60): download KEY60S .TAB days from PDS-PPI.

Lists the PDS-PPI KEY60S day-of-year buckets, downloads the requested year's daily .TAB files
into data/skr/ (resumable: existing files skipped; the seq digit varies per day so we scrape the
directory index rather than guess). After download, run:

    uv run python -m jansky_research.skr --out .   # (no --offline) parses data/skr/ + Horizons

Usage:  uv run python scripts/skr_real.py --year 2017 [--doy-min 100 --doy-max 258]
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from jansky_research.skr import PDS_BASE  # noqa: E402

DATA_DIR = REPO / "data" / "skr"


def list_bucket(year: int, hundreds: int) -> list[str]:
    """Filenames of the KEY60S .TAB products in the T<year><hundreds>XX day-of-year bucket."""
    url = f"{PDS_BASE}/T{year}{hundreds}XX/"
    try:
        html = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "replace")
    except Exception as exc:  # noqa: BLE001
        print(f"  bucket {url} unavailable: {exc!r}", flush=True)
        return []
    # the trailing seq/version char is alphanumeric (e.g. _Z, _P, _5), not always a digit
    return sorted(set(re.findall(rf"RPWS_KEY__{year}\d{{3}}_[A-Z0-9]\.TAB", html)))


def main() -> int:
    ap = argparse.ArgumentParser(description="Download Cassini KEY60S days (plan 60)")
    ap.add_argument("--year", type=int, default=2017)
    ap.add_argument("--doy-min", type=int, default=1)
    ap.add_argument("--doy-max", type=int, default=366)
    args = ap.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    names: list[str] = []
    for hundreds in range(4):  # DOY 000-099, 100-199, 200-299, 300-366
        names += list_bucket(args.year, hundreds)
    # filter to the requested DOY window
    sel = [n for n in names if args.doy_min <= int(n[10:13]) <= args.doy_max]
    print(f"{len(sel)} KEY60S days in {args.year} DOY {args.doy_min}-{args.doy_max}", flush=True)

    got = skipped = failed = 0
    for name in sel:
        doy = int(name[10:13])
        url = f"{PDS_BASE}/T{args.year}{doy // 100}XX/{name}"
        dest = DATA_DIR / name
        if dest.exists() and dest.stat().st_size > 0:
            skipped += 1
            continue
        try:
            urllib.request.urlretrieve(url, dest)
            got += 1
            if got % 20 == 0:
                print(f"  ...{got} downloaded", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"  FAILED {name}: {exc!r}", flush=True)
            failed += 1
    print(f"done: {got} downloaded, {skipped} already present, {failed} failed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
