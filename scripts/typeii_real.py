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
    uv run --extra typeii python scripts/typeii_real.py \
        --dates 2024-05-14 2024-05-15 2024-05-16 \
        --cme data/typeii/lasco_cme.csv

``--cme`` is a CDAW LASCO CME table with columns ``onset_hr,speed_kms,width_deg`` (hours on the
same clock as the burst days). Writes results/typeii_metrics.json (is_real=True) + regenerates the
paper macros.

The detector and its synthetic recover-a-known (completeness-vs-SNR curve + the Gopalswamy
fast-and-wide CME association wiring) are validated in core CI without any of this -- that is the
shippable deliverable; this driver produces the real streamed census on demand.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research.typeii import _figure, _write_macros, real_census  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description="OVRO-LWA type II streamed real census (plan 50)")
    ap.add_argument(
        "--dates",
        nargs="+",
        required=True,
        metavar="YYYY-MM-DD",
        help="days to stream from AWS Open Data (in memory, no disk)",
    )
    ap.add_argument(
        "--cme", required=True, help="CDAW LASCO CME CSV (onset_hr,speed_kms,width_deg)"
    )
    args = ap.parse_args()

    metrics = real_census(args.dates, args.cme)
    (REPO / "results").mkdir(exist_ok=True)
    (REPO / "results" / "typeii_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(metrics, REPO / "papers" / "typeii" / "figures")
    _write_macros(metrics, REPO / "papers" / "typeii" / "generated" / "macros.tex")
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
