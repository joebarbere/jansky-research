#!/usr/bin/env python3
"""Plan 90 GATE 0: is each source's snapshot set usable for a binomial constraint?

Offline. The duty-cycle constraint assumes each snapshot is an independent draw on pulse
phase. VAST observes on a roughly fortnightly cadence, which is not random with respect to
any LPT period, so this tests it rather than asserting it -- and first checks whether the
catalogued period is precise enough for the test itself to mean anything.

Writes results/lptduty_gate0.json.

    uv run python scripts/lptduty_gate0.py
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jansky_research import lptduty as ld  # noqa: E402

# Family-wise significance: a Bonferroni-corrected p below this counts as clustered.
# (Previously written as 0.01 with an unexplained "* 5" at each use, which recorded a
# significance level in the JSON that was not the one governing the verdict.)
ALPHA_FAMILYWISE = 0.05


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--epochs", type=Path, default=Path("results/lptv_vast_epochs.csv"))
    ap.add_argument("--catalog", type=Path, default=Path("data/lpt_sample.csv"))
    ap.add_argument("--out", type=Path, default=Path("results/lptduty_gate0.json"))
    args = ap.parse_args(argv)

    rows = ld.load_epochs(args.epochs)
    periods = ld.read_periods(args.catalog)
    precision = ld.read_period_precision(args.catalog)
    by: dict[str, list[ld.EpochRow]] = defaultdict(list)
    for r in rows:
        by[r.name].append(r)

    n_tested = sum(1 for n in by if n in periods)
    per_source = {}
    for name, srows in sorted(by.items()):
        period = periods.get(name)
        if period is None:
            per_source[name] = {"testable": False, "reason": "no period in the catalogue"}
            continue
        ps = ld.phase_sampling(srows, period)
        quoted = precision.get(name, float("inf"))
        # Phase stays coherent across the baseline only if the period is known better than
        # required_period_precision_s. Otherwise the computed phase smears and the test
        # cannot detect clustering that may still be there: inconclusive, not reassuring.
        adequate = quoted <= ps.required_period_precision_s
        p_corrected = min(1.0, ps.rayleigh_p * n_tested)  # Bonferroni over sources tested
        per_source[name] = {
            "testable": True,
            "n_epochs": ps.n_epochs,
            "period_s": period,
            "quoted_period_precision_s": quoted,
            "required_period_precision_s": ps.required_period_precision_s,
            "period_precision_adequate": adequate,
            "baseline_days": ps.baseline_days,
            "rayleigh_z": ps.rayleigh_z,
            "rayleigh_p": ps.rayleigh_p,
            "rayleigh_p_bonferroni": p_corrected,
            "kuiper_v": ps.kuiper_v,
            "clustered": bool(adequate and p_corrected < ALPHA_FAMILYWISE),
            "verdict": (
                "inconclusive: catalogued period too imprecise for coherent phase"
                if not adequate
                else (
                    "phase sampling NOT uniform - independence assumption fails"
                    if p_corrected < ALPHA_FAMILYWISE
                    else "phase sampling consistent with uniform"
                )
            ),
        }

    flagged = [n for n, d in per_source.items() if d.get("clustered")]
    inconclusive = [
        n for n, d in per_source.items() if d.get("testable") and not d["period_precision_adequate"]
    ]
    payload = {
        "slice": "lptduty",
        "stage": "gate0-aliasing",
        "is_real": True,
        "question": (
            "Are VAST snapshots uniformly distributed in each source's pulse phase? The "
            "binomial constraint in lptduty_metrics.json assumes they are."
        ),
        "method": (
            "Rayleigh Z and Kuiper V on phases referenced to each source's first epoch (the "
            "zero point is arbitrary; clustering is invariant under a phase shift, so no "
            "ephemeris is needed to test uniformity). Bonferroni-corrected over the sources "
            "tested."
        ),
        "caveats": [
            "quoted_period_precision_s is inferred from the catalogue's decimal places -- a "
            "proxy for the published uncertainty, not the uncertainty itself. The ephemeris "
            "audit must replace it with real values from the discovery papers.",
            "A period derivative large enough to matter over the baseline would break phase "
            "coherence even when the period is quoted precisely; pdot is not folded in here.",
            "A uniform result does not prove independence, only that this test found no "
            "departure from it.",
        ],
        "n_sources_tested": n_tested,
        "alpha_familywise": ALPHA_FAMILYWISE,
        "flagged_clustered": flagged,
        "inconclusive_period_precision": inconclusive,
        "per_source": per_source,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(args.out)

    print(f"wrote {args.out}")
    print(f"{'source':32s} {'p(Bonf)':>9} {'V':>6}  verdict")
    for name, d in sorted(per_source.items()):
        if not d.get("testable"):
            print(f"{name:32s} {'-':>9} {'-':>6}  {d['reason']}")
            continue
        print(f"{name:32s} {d['rayleigh_p_bonferroni']:9.2e} {d['kuiper_v']:6.3f}  {d['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
