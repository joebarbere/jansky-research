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

GATE0 = Path("results/lptduty_gate0.json")
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
    # The binomial constraint assumes snapshots are independent draws on pulse phase.
    # GATE 0 tests that per source; stamp each result with its verdict so a reader of this
    # file cannot pick up a number without the caveat that governs it.
    gate0 = {}
    if GATE0.exists():
        gate0 = json.loads(GATE0.read_text()).get("per_source", {})
    by_source: dict[str, list[ld.EpochRow]] = defaultdict(list)
    for r in rows:
        by_source[r.name].append(r)

    per_source = {}
    for name, srows in sorted(by_source.items()):
        brightest = max(abs(r.v_mjy) for r in srows)
        detections = [r for r in srows if ld._is_detection(r)]  # incl. the leakage veto
        # The flux we know this source can produce, where it has produced one.
        headline_flux = max(abs(r.v_mjy) for r in detections) if detections else None
        grid = {}
        for s in FLUX_GRID:
            c = ld.duty_constraint(srows, pulse_mjy=s)
            grid[f"{s:g}"] = {
                "effective_epochs": c.effective_epochs,
                "p_point": c.p_point,
                "p_upper_95": None if math.isinf(c.p_upper_95) else c.p_upper_95,
                # capped at 1.0 in duty_constraint: a "limit" at or above 1 excludes nothing
                "unconstrained": math.isinf(c.p_upper_95) or c.p_upper_95 >= 1.0,
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
        period = periods.get(name)
        snap_med = float(sorted(r.duration_s for r in srows)[len(srows) // 2])
        wtp = (
            (ld.PULSE_WIDTH_S.get(name, 0.0) + snap_med) / period if period and period > 0 else None
        )
        # the implied active fraction: p divided by this source's own (w+T)/P. This is the
        # physical quantity the note cares about, and quoting p without it invites reading a
        # tight-looking limit as strong when its own T/P makes it weak (or vacuous).
        implied_fa = None
        implied_fa_lim = None
        if wtp:
            if headline and headline["p_point"] is not None:
                implied_fa = round(headline["p_point"] / wtp, 3)
            lim5 = grid.get("5", {}).get("p_upper_95")
            if lim5 is not None and not detections:
                implied_fa_lim = round(min(lim5 / wtp, 9.99), 3)
        per_source[name] = {
            "n_epochs_measured": len(srows),
            "wt_over_p": round(wtp, 4) if wtp else None,
            "implied_f_active": implied_fa,
            "implied_f_active_limit_5mjy": implied_fa_lim,
            "implied_limit_vacuous": bool(implied_fa_lim is not None and implied_fa_lim >= 1.0),
            "efficiency_kept_at_headline": (
                round(headline["effective_epochs"] / len(srows), 4) if headline else None
            ),
            "efficiency_kept_at_5mjy": round(grid["5"]["effective_epochs"] / len(srows), 4),
            "total_exposure_s": float(sum(r.duration_s for r in srows)),
            "median_sigma_v_mjy": float(sorted(r.e_v for r in srows)[len(srows) // 2]),
            "brightest_abs_v_mjy": brightest,
            "n_detections": len(detections),
            "period_s": periods.get(name),
            "snapshot_median_s": float(sorted(r.duration_s for r in srows)[len(srows) // 2]),
            "headline": headline,
            "vs_assumed_flux_mjy": grid,
            "phase_sampling": {
                "verdict": gate0.get(name, {}).get("verdict", "not tested"),
                "independence_assumption_holds": (
                    None
                    if name not in gate0 or not gate0[name].get("testable")
                    else not gate0[name].get("clustered", False)
                    and gate0[name].get("period_precision_adequate", False)
                ),
            },
            # False where GATE 0 showed the snapshots are not uniform in pulse phase, or
            # could not tell: the binomial model does not apply, so p above is not a limit.
            "constraint_valid": (
                bool(
                    gate0.get(name, {}).get("period_precision_adequate")
                    and not gate0.get(name, {}).get("clustered")
                )
                if name in gate0 and gate0[name].get("testable")
                else False
            ),
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
            "VAST pointings are not randomly phased against any LPT period; GATE 0 "
            "(phase-uniformity, Rayleigh + Kuiper, results/lptduty_gate0.json) tests this "
            "per source. A uniform marginal does not prove independence: if activity "
            "persists across a cycle, correlated snapshots over-disperse the counts and the "
            "Poisson limits are too tight (by up to ~55% under a one-epoch-per-cycle "
            "collapse); the point estimates are unbiased either way.",
        ],
        "threshold_sigma": ld.DETECT_THRESHOLD_SIGMA,
        "min_efficiency": ld.MIN_EFFICIENCY,
        "flux_grid_mjy": FLUX_GRID,
        "n_sources": len(per_source),
        "n_epochs_measured_total": len(rows),
        "per_source": per_source,
    }
    # are the three detection-source rates distinguishable? (G-test vs one common rate)
    det_named = [(n, v) for n, v in per_source.items() if v["n_detections"] and v["headline"]]
    if len(det_named) >= 2:
        payload["common_rate_test"] = ld.common_rate_test(
            [v["n_detections"] for _, v in det_named],
            [v["headline"]["effective_epochs"] for _, v in det_named],
        )
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
    # regenerate the paper macros from the freshly written metrics plus the committed
    # gate0/phase JSONs, so the note cannot drift from the evidence (this call is the
    # previously-missing producer of papers/lptduty/generated/macros.tex)
    gate0_path = args.out.parent / "lptduty_gate0.json"
    phase_path = args.out.parent / "lptduty_phase.json"
    if gate0_path.exists() and phase_path.exists():
        ld.write_paper_assets(
            payload,
            json.loads(gate0_path.read_text()),
            json.loads(phase_path.read_text()),
            args.out.parent.parent / "papers" / "lptduty" / "generated" / "macros.tex",
        )
        print("wrote papers/lptduty/generated/macros.tex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
