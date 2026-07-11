#!/usr/bin/env python
"""JBO glitch waiting-time census -- real-leg driver (plan 48, F11).

Scrapes the live Jodrell Bank glitch table, classifies every pulsar with >=5 glitches by its
inter-glitch waiting-time distribution (exponential / quasi-periodic / clustered), checks the known
quasi-periodic glitchers, and diffs against the end-2018 Basu+2022 subset. Writes
results/glitchpop_metrics.json + results/glitchpop_census.json and regenerates the paper macros. The
classifier + synthetic recover-a-known run offline in CI.

Usage:
    uv run python scripts/glitchpop_real.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research import glitchpop as gp  # noqa: E402
from jansky_research.glitchpop_real import run_real_census  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="JBO glitch waiting-time census (real).")
    ap.add_argument("--out", default=str(REPO))
    ap.add_argument("--min-glitches", type=int, default=gp.MIN_GLITCHES)
    args = ap.parse_args(argv)

    metrics = run_real_census(args.out, min_glitches=args.min_glitches)
    op = Path(args.out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "glitchpop_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    gp._figure(metrics, op / "papers" / "glitchpop" / "figures")
    gp._write_macros(metrics, op / "papers" / "glitchpop" / "generated" / "macros.tex")

    print(
        json.dumps(
            {k: v for k, v in metrics.items() if k not in ("flipped", "newly_qualified")}, indent=2
        )
    )
    print("\nClassification flips (post-2018 data changed the class):")
    for f in metrics.get("flipped", []):
        print(f"  {f['jname']:14s} {f['was']} -> {f['now']}  (n {f['n_pre']} -> {f['n_now']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
