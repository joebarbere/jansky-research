#!/usr/bin/env python
"""e-Callisto megaconstellation RFI-trend census -- real-leg driver (plan 54, F17).

The e-Callisto FITS archive (soleil.i4ds.ch, open, no auth) is an accidental 15-year RFI record.
This driver samples one representative 15-min spectrum per station-month over 2012-2026, reduces
each to a burst-immune, gain-cancelling occupancy metric, and trends the Starlink unintended-
emission (UEM) band -- attributing any post-2019 rise to the public Starlink constellation count.

PRIMARY metric = the narrowband UEM-line excess (`line_vs_adjacent`): the level at the 137/150/175
MHz Starlink UEM lines over their adjacent clean channels. It is self-normalizing within the UEM
band, so it cancels station gain drift AND survives the per-station RFI-avoidance notches that make
a fixed FM control unusable at many stations (HUMAIN notches the FM band and the 137 MHz line
entirely -- verified on real data). A station-adaptive band differential (UEM minus the best-sampled
clean control) is a cross-check. Config-stability is enforced: a station must keep the SAME sampled
UEM lines across its retained months, else the differing months are dropped.

Nothing is downloaded in bulk: each ~100 kB gzipped FITS is fetched, reduced in memory, and freed.

Usage:
    uv run python scripts/rfitrend_real.py                       # default stations, 2012-2026
    uv run python scripts/rfitrend_real.py --stations HUMAIN ALMATY --start 2012 --end 2026

Writes results/rfitrend_metrics.json (is_real=True; per-station trends + network-pooled trend +
Starlink correlation + Perez+2020 reproduction) and regenerates the paper macros. The synthetic
recover-a-known + tests run in core CI without any network; this produces the real census on demand.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from jansky_research import rfitrend as rf  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="e-Callisto megaconstellation RFI-trend real census.")
    p.add_argument("--stations", nargs="+", default=list(rf.ECALLISTO_STATIONS))
    p.add_argument("--start", type=int, default=2012)
    p.add_argument("--end", type=int, default=2026)
    p.add_argument("--out", default=str(REPO))
    args = p.parse_args(argv)

    import json

    metrics = rf._real_trend(
        args.out, stations=tuple(args.stations), start_year=args.start, end_year=args.end
    )
    op = Path(args.out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "rfitrend_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    rf._figure(metrics, op / "papers" / "rfitrend" / "figures")
    rf._write_macros(metrics, op / "papers" / "rfitrend" / "generated" / "macros.tex")
    print(json.dumps({k: v for k, v in metrics.items() if k != "per_station"}, indent=2))
    for st, d in metrics["per_station"].items():
        print(
            f"  {st:8s} n={d.get('n_months', 0):3d} lines={d.get('stable_lines')} "
            f"ctrl={d.get('control_name')} line_slope={d.get('line_excess_slope_per_yr')} "
            f"p={d.get('line_excess_p')} corr={d.get('corr_line_excess_starlink')} "
            f"perez={d.get('perez_2012_2019_change')}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
