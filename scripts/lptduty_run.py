#!/usr/bin/env python3
"""Plan 90 runner: constrain how often each LPT is on, from the committed VAST sweep.

Offline — reads only committed evidence (`results/lptv_vast_epochs.csv`,
`data/lpt_sample.csv`) and writes `results/lptduty_metrics.json`.

Every limit is a function of the assumed pulse flux, so the run reports a grid of assumed
fluxes per source rather than one number. The headline flux for each source is its own
brightest measured |V| where one exists (the pulse we know it can produce); sources with no
detection are reported across the grid only.

    uv run python scripts/lptduty_run.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jansky_research import lptduty as ld  # noqa: E402

EPOCHS = Path("results/lptv_vast_epochs.csv")
CATALOG = Path("data/lpt_sample.csv")
OUT = Path("results/lptduty_metrics.json")

# Assumed pulse fluxes (mJy) at which each source's constraint is evaluated. Spans the range
# of measured LPT pulse brightness in this sweep (3.6-164 mJy in the lptv census) down to
# near the per-epoch noise, where the constraint necessarily dies.
FLUX_GRID = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--epochs", type=Path, default=EPOCHS)
    ap.add_argument("--catalog", type=Path, default=CATALOG)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    rows = ld.load_epochs(args.epochs)
    periods = ld.read_periods(args.catalog)
    by_source: dict[str, list[ld.EpochRow]] = defaultdict(list)
    for r in rows:
        by_source[r.name].append(r)

    per_source = {}
    for name, srows in sorted(by_source.items()):
        brightest = max(abs(r.v_mjy) for r in srows)
        detections = [r for r in srows if abs(r.v_mjy) / r.e_v >= ld.DETECT_THRESHOLD_SIGMA]
        # The flux we know this source can produce, where it has produced one.
        headline_flux = max(abs(r.v_mjy) for r in detections) if detections else None
        grid = {}
        for s in FLUX_GRID:
            c = ld.duty_constraint(srows, pulse_mjy=s)
            grid[f"{s:g}"] = {
                "effective_epochs": c.effective_epochs,
                "p_point": c.p_point,
                "p_upper_95": None if math.isinf(c.p_upper_95) else c.p_upper_95,
                "unconstrained": math.isinf(c.p_upper_95),
            }
        headline = None
        if headline_flux is not None:
            c = ld.duty_constraint(srows, pulse_mjy=headline_flux)
            headline = {
                "assumed_pulse_mjy": headline_flux,
                "effective_epochs": c.effective_epochs,
                "p_point": c.p_point,
                "p_upper_95": None if math.isinf(c.p_upper_95) else c.p_upper_95,
            }
        per_source[name] = {
            "n_epochs_measured": len(srows),
            "total_exposure_s": float(sum(r.duration_s for r in srows)),
            "median_sigma_v_mjy": float(sorted(r.e_v for r in srows)[len(srows) // 2]),
            "brightest_abs_v_mjy": brightest,
            "n_detections": len(detections),
            "period_s": periods.get(name),
            "snapshot_median_s": float(sorted(r.duration_s for r in srows)[len(srows) // 2]),
            "headline": headline,
            "vs_assumed_flux_mjy": grid,
            # f_active and (w+T)/P are not separately identifiable from counts alone.
            "identifiable_factors": False,
        }

    payload = {
        "slice": "lptduty",
        "is_real": True,
        "source": ("committed lptv VAST sweep (results/lptv_vast_epochs.csv); no new data fetched"),
        "quantity": "p = f_active * (w + T) / P, the per-snapshot probability of catching a pulse",
        "caveats": [
            "Only the product p is identifiable; f_active and the in-period duty cycle "
            "separate only with an ephemeris good enough to phase each snapshot.",
            "The denominator is efficiency-weighted exposure (sum of per-epoch detection "
            "probability at the assumed flux), not the epoch count -- the frblens lesson.",
            "Efficiencies below MIN_EFFICIENCY are floored to zero so that many shallow "
            "epochs cannot sum into sensitivity that does not exist.",
            "VAST pointings are not randomly phased against any LPT period; aliasing between "
            "the survey cadence and a period is NOT yet tested (plan 90 GATE 0).",
        ],
        "threshold_sigma": ld.DETECT_THRESHOLD_SIGMA,
        "min_efficiency": ld.MIN_EFFICIENCY,
        "flux_grid_mjy": FLUX_GRID,
        "n_sources": len(per_source),
        "n_epochs_measured_total": len(rows),
        "per_source": per_source,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(args.out)
    print(f"wrote {args.out}: {len(per_source)} sources, {len(rows)} measured epochs")
    for name, d in sorted(per_source.items()):
        h = d["headline"]
        if h:
            print(
                f"  {name:32s} k={d['n_detections']} at {h['assumed_pulse_mjy']:.1f} mJy: "
                f"p={h['p_point']:.3f} (eff epochs {h['effective_epochs']:.0f})"
            )
        else:
            g = d["vs_assumed_flux_mjy"]["5"]
            lim = g["p_upper_95"]
            print(
                f"  {name:32s} k=0, p < {lim:.3f} at 5 mJy"
                if lim is not None
                else f"  {name:32s} k=0, unconstrained at 5 mJy"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
