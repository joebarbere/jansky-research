#!/usr/bin/env python
"""OVRO-LWA type II census -- real-leg driver (plan 50, jansky-research).

DATA-ACCESS NOTE (GATE-0, 2026-07-09): the OVRO-LWA solar data portal
(`https://www.ovsa.njit.edu/lwadata-query`) is a JavaScript SPA behind a **Cloudflare Turnstile
bot challenge**, so the Level-1 beamforming dynamic-spectrum FITS cannot be fetched by a script.
The files must be downloaded INTERACTIVELY through the portal into ``data/typeii/`` (pattern
``ovro-lwa.lev1_bmf_256ms_96kHz.YYYY-MM-DD.dspec_I.fits``; ~0.6 GB/day, 13.4-86.9 MHz, 256 ms).

Once the FITS are local, also drop a CDAW LASCO CME table at ``data/typeii/lasco_cme.csv`` with
columns ``onset_hr,speed_kms,width_deg`` (hours from the same epoch as the FITS days), then run:

    uv run python scripts/typeii_real.py

It sweeps every local dspec FITS with the tested `detect_typeii`, cross-matches detections to the
LASCO CMEs, and writes results/typeii_metrics.json (is_real=True) + regenerates the paper macros.

The detector and its synthetic recover-a-known (completeness/purity + the Gopalswamy fast-and-wide
CME association) are validated in core CI without any of this --- that is the shippable deliverable;
this driver produces the real OVRO-LWA census when the FITS are in hand.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research.typeii import (  # noqa: E402
    DATA_DIR,
    _figure,
    _write_macros,
    real_census,
)


def main() -> int:
    if not any(DATA_DIR.glob("ovro-lwa.*dspec_I.fits")):
        print(
            f"No OVRO-LWA dspec FITS in {DATA_DIR}. The portal is Turnstile-gated -- download the "
            "daily FITS interactively via https://www.ovsa.njit.edu/lwadata-query first "
            "(see this file's docstring).",
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
