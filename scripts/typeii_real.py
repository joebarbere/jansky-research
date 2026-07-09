#!/usr/bin/env python
"""OVRO-LWA type II census -- real-leg driver (plan 50, jansky-research).

DATA ACCESS (corrected 2026-07-09): the OVRO-LWA solar dynamic spectra are on **AWS Open Data** --
public bucket ``ovro-lwa-solar``, path ``spec_fits/<YYYY>/<YYYYMMDD>.fits`` -- directly
downloadable with no login and NO bot challenge. (The Cloudflare Turnstile on
`ovsa.njit.edu/lwadata-query` only gates the query *UI*, not the data.) The daily files are large
(~1.7 GB; a 4D I/V dynamic spectrum, ~15-85 MHz, ~0.26 s), so a multi-year census is a bulk
download, but it is not access-blocked.

Usage:
    # 1. download the wanted day(s) -- no auth needed:
    uv run python scripts/typeii_real.py --download 2024-05-14 2024-05-15
    #    (or directly: curl -O https://ovro-lwa-solar.s3-us-west-2.amazonaws.com/spec_fits/2024/20240514.fits)
    # 2. drop a CDAW LASCO CME table at data/typeii/lasco_cme.csv (onset_hr,speed_kms,width_deg)
    # 3. run the census:
    uv run python scripts/typeii_real.py

The detector and its synthetic recover-a-known (completeness-vs-SNR curve + the Gopalswamy
fast-and-wide CME association wiring) are validated in core CI without any of this -- that is the
shippable deliverable; this driver produces the real census once the (large) FITS are local.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research.typeii import (  # noqa: E402
    DATA_DIR,
    _figure,
    _write_macros,
    real_census,
    s3_dspec_url,
)


def main() -> int:
    ap = argparse.ArgumentParser(description="OVRO-LWA type II real census (plan 50)")
    ap.add_argument("--download", nargs="*", default=[], metavar="YYYY-MM-DD",
                    help="download these day(s) from AWS Open Data into data/typeii/")
    args = ap.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for date in args.download:
        url = s3_dspec_url(date)
        dest = DATA_DIR / f"{date.replace('-', '')}.fits"
        if dest.exists() and dest.stat().st_size > 0:
            print(f"  {dest.name} present, skipping", flush=True)
            continue
        print(f"  downloading {url} (~1.7 GB) ...", flush=True)
        urllib.request.urlretrieve(url, dest)
    if args.download:
        return 0

    if not any(DATA_DIR.glob("*.fits")):
        print(
            f"No OVRO-LWA dspec FITS in {DATA_DIR}. Download first, e.g.\n"
            "  uv run python scripts/typeii_real.py --download 2024-05-14",
            file=sys.stderr,
        )
        return 1
    metrics = real_census(DATA_DIR)
    (REPO / "results").mkdir(exist_ok=True)
    (REPO / "results" / "typeii_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(metrics, REPO / "papers" / "typeii" / "figures")
    _write_macros(metrics, REPO / "papers" / "typeii" / "generated" / "macros.tex")
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
