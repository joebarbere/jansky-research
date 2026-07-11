#!/usr/bin/env python
"""PTE-II single-pulse heavy-tail census -- real-leg driver (plan 47, F10).

Downloads the open PTE-II SQLite database (~1.5 GB, GitHub LFS, no auth), runs the per-source
giant-pulse test across all 363 pulsars, and cross-matches ATNF spin-down luminosity. Writes
results/pte2_metrics.json (is_real=True) + results/pte2_census.json (full per-pulsar table) and
regenerates the paper macros/figure. The metric + synthetic recover-a-known run offline in CI.

Usage:
    uv run python scripts/pte2_real.py --cache /tmp/pte2_cache
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research import pte2  # noqa: E402
from jansky_research.pte2_real import run_real_census  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="PTE-II heavy-tail census (real).")
    ap.add_argument("--out", default=str(REPO))
    ap.add_argument(
        "--cache", default=None, help="download/extract cache dir (default: system temp)"
    )
    ap.add_argument("--min-pulses", type=int, default=pte2.MIN_PULSES)
    args = ap.parse_args(argv)

    metrics = run_real_census(args.out, cache_dir=args.cache, min_pulses=args.min_pulses)
    op = Path(args.out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "pte2_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    pte2._figure(metrics, op / "papers" / "pte2" / "figures")
    pte2._write_macros(metrics, op / "papers" / "pte2" / "generated" / "macros.tex")

    print(json.dumps({k: v for k, v in metrics.items() if k != "top_heavy"}, indent=2))
    print("\nTop heavy-tailed sources:")
    for r in metrics.get("top_heavy", [])[:15]:
        print(
            f"  {r['jname']:14s} n={r['n']:5d} n_giant={r['n_giant']:3d} "
            f"excess={r['excess']:.2f} gamma={r['gamma']} heavy={r['heavy_tailed']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
