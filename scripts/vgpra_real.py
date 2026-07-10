#!/usr/bin/env python
"""Voyager 2 PRA rotation-period re-derivation -- real-leg driver (plan 46, F9).

Downloads the two open PDS-PPI encounter volumes (no auth; GATE-0 verified) --
`VG2-U-PRA-3-RDR-LOWBAND-6SEC-V1.0` (Uranus, 49 MB) and
`VG2-N-PRA-3-RDR-LOWBAND-6SEC-V1.0` (Neptune, 79 MB) -- parses the fixed-width ASCII major frames,
extracts the band-integrated flux, detects UKR/NKR burst episodes, and runs the `frbperiod`
Rayleigh Z^2 posterior with a bootstrap (few-cycle-honest) uncertainty. Writes
results/vgpra_metrics.json (is_real=True) with the three-way comparison against the historical radio
value and Lamy+2025, and regenerates the paper macros.

The metric, its synthetic recover-a-known, and all unit tests run offline in CI without any of this;
this driver produces the real re-derivation on demand.

Usage:
    uv run python scripts/vgpra_real.py                 # both planets, default settings
    uv run python scripts/vgpra_real.py --cache /tmp/vgpra   # reuse a download cache
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research import vgpra as vg  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Voyager 2 PRA Uranus/Neptune period re-derivation (real)."
    )
    ap.add_argument("--out", default=str(REPO))
    ap.add_argument("--cache", default=None, help="download cache dir (default: system temp)")
    args = ap.parse_args(argv)

    metrics = vg._real_analysis(cache_dir=args.cache)
    op = Path(args.out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "vgpra_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    vg._figure(metrics, op / "papers" / "vgpra" / "figures")
    vg._write_macros(metrics, op / "papers" / "vgpra" / "generated" / "macros.tex")

    print(json.dumps({k: v for k, v in metrics.items() if k != "planets"}, indent=2))
    for planet, d in metrics["planets"].items():
        pub = vg.PUBLISHED_HR[planet]
        hist = next(iter(pub.values()))
        print(
            f"  {planet:8s} P={d.get('best_period_hr')} +/- {d.get('boot_sigma_hr')} h "
            f"(band [{d.get('boot_lo_hr')},{d.get('boot_hi_hr')}]) z2={d.get('z2')} "
            f"n_bursts={d.get('n_bursts')} span={d.get('span_hr')}h | "
            f"hist {hist[0]}h -> consistent={d.get('consistent_hist')}"
            + (f" | Lamy consistent={d.get('consistent_lamy')}" if "consistent_lamy" in d else "")
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
