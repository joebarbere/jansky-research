#!/usr/bin/env python
"""OVRO-LWA type II census -- real-leg driver (plan 50, jansky-research).

STREAMING, IN MEMORY, NO DISK. The OVRO-LWA solar dynamic spectra are on AWS Open Data (public
bucket ``ovro-lwa-solar``, ``spec_fits/<YYYY>/<YYYYMMDD>.fits``; no login, no CAPTCHA -- the
Cloudflare Turnstile only gates the query UI). Each daily file is ~1.7 GB, so this driver does NOT
download them: for each requested day it opens the S3 FITS lazily (astropy ``use_fsspec``) and
range-reads only the Stokes-I plane in time-chunks, block-averaging to ~4 s bins on the fly, so a
day is held entirely in memory (peak ~one reduced chunk) and freed before the next -- nothing
touches disk.

Needs the ``typeii`` extra (fsspec + aiohttp) for the lazy cloud reads:  uv sync --extra typeii

Usage:
    # an explicit day list:
    uv run --extra typeii python scripts/typeii_real.py --dates 2024-05-14 2024-05-15
    # OR every observing day in the archive over a date range (listed from the S3 bucket):
    uv run --extra typeii python scripts/typeii_real.py --start 2024-04-01 --end 2026-07-09

The CDAW LASCO CME catalogue, the SWPC/HEK GOES flare list, and the SILSO sunspot series are all
fetched automatically for the date span. Writes results/typeii_metrics.json (is_real=True; the full
plan product set: event list + CME + GOES association + occurrence vs cycle phase) + regenerates
the paper macros.

The detector and its synthetic recover-a-known (completeness-vs-SNR curve + the Gopalswamy
fast-and-wide CME association wiring) are validated in core CI without any of this -- that is the
shippable deliverable; this driver produces the real streamed census on demand.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research.typeii import S3_SPEC_FITS, _figure, _write_macros, real_census  # noqa: E402


def list_observing_days(start: str, end: str) -> list[str]:
    """Observing days (YYYY-MM-DD) present in the S3 bucket within [start, end], via the S3 list API."""
    bucket = S3_SPEC_FITS.split("/spec_fits")[0]
    days: list[str] = []
    for year in range(int(start[:4]), int(end[:4]) + 1):
        token = ""
        while True:
            u = f"{bucket}/?list-type=2&prefix=spec_fits/{year}/&max-keys=1000{token}"
            xml = urllib.request.urlopen(u, timeout=60).read().decode("utf-8", "replace")
            for k in re.findall(r"<Key>spec_fits/\d{4}/(\d{8})\.fits</Key>", xml):
                iso = f"{k[:4]}-{k[4:6]}-{k[6:8]}"
                if start <= iso <= end:
                    days.append(iso)
            nt = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", xml)
            if not nt:
                break
            token = "&continuation-token=" + urllib.parse.quote(nt.group(1))
    return sorted(set(days))


def main() -> int:
    ap = argparse.ArgumentParser(description="OVRO-LWA type II streamed real census (plan 50)")
    ap.add_argument(
        "--dates",
        nargs="*",
        metavar="YYYY-MM-DD",
        help="explicit days to stream from AWS Open Data (in memory, no disk)",
    )
    ap.add_argument("--start", help="census start YYYY-MM-DD (lists all archive days in range)")
    ap.add_argument("--end", help="census end YYYY-MM-DD")
    args = ap.parse_args()

    if args.start and args.end:
        dates = list_observing_days(args.start, args.end)
        print(f"{len(dates)} observing days in {args.start}..{args.end}", flush=True)
    elif args.dates:
        dates = args.dates
    else:
        ap.error("give --dates ... OR --start/--end")

    metrics = real_census(dates)  # CDAW CMEs + HEK flares + SILSO fetched automatically
    (REPO / "results").mkdir(exist_ok=True)
    (REPO / "results" / "typeii_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(metrics, REPO / "papers" / "typeii" / "figures")
    _write_macros(metrics, REPO / "papers" / "typeii" / "generated" / "macros.tex")
    m = {k: v for k, v in metrics.items() if k != "event_list"}
    print(json.dumps(m, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
