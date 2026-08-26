#!/usr/bin/env python3
"""Plan 90 increment 3: split f_active from the in-period duty cycle, where an ephemeris allows.

Offline. Writes results/lptduty_phase.json.

The product `p = f_active * (w + T) / P` separates only for sources with a **published
reference epoch** and a period precise enough to keep phase coherent across the VAST
baseline. The 2026-08-21 ephemeris audit (survey/lptduty-findings.md) found that most of the
sample fails the first condition, not the second: several papers publish a precise period and
no PEPOCH at all, which makes physical phase unavailable no matter how good the period is.

Only ephemerides with a published epoch -- or, for J183950, the repo's own documented anchor
-- are encoded here. Assuming an epoch "near the campaign" would manufacture a phase and is
exactly what this file refuses to do.

    uv run python scripts/lptduty_phase.py
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

# Ephemerides with a PUBLISHED reference epoch, transcribed 2026-08-21 from the papers'
# arXiv full text. pulse_width_s is the published on-pulse width; where a range is quoted the
# WIDER end is used, because a wider window admits more snapshots and so gives the more
# conservative (larger) f_active denominator.
EPHEMERIDES = {
    "ASKAP J175534.9-252749.1": ld.Ephemeris(
        name="ASKAP J175534.9-252749.1",
        period_s=4186.3285,
        sigma_period_s=0.0002,
        pepoch_mjd=59965.03792,
        pulse_width_s=30.0,
        reference="mcsweeney2025 (arXiv:2507.14448), tabulated P and PEPOCH 59965.03792(3)",
    ),
    "ASKAP J1832-0911": ld.Ephemeris(
        name="ASKAP J1832-0911",
        period_s=2656.247,
        sigma_period_s=0.001,
        pepoch_mjd=60344.0,
        pulse_width_s=240.0,
        reference="wang2025 (arXiv:2411.16606): P=2656.247(1) s, period epoch MJD 60344, width 2-4 min",
    ),
    "ASKAP J183950.5-075635": ld.Ephemeris(
        name="ASKAP J183950.5-075635",
        period_s=23221.740,
        sigma_period_s=0.332,
        pepoch_mjd=60358.245243,
        pulse_width_s=710.0,
        reference=(
            "lee2025 (arXiv:2501.09133) P=23221.740(332); epoch is the repo's own documented "
            "anchor (lptv.J1839_T0_MJD, +/-0.013 cycle systematic), not a published PEPOCH"
        ),
    ),
}

# Sources deliberately absent, and why -- so the gap is a recorded decision, not an oversight.
NO_PUBLISHED_EPOCH = {
    "GPM J1839-10": "period known to 2e-4 s over a 34-yr baseline, but the paper tabulates no "
    "PEPOCH; phase cannot be assigned without inventing one. Also 50-70% nulling.",
    "GLEAM-X J162759.5-523504.3": "no PEPOCH; and its own best-fit Pdot (6e-10) would move "
    "phase by tens of cycles across the baseline.",
    "ASKAP J1935+2148": "no PEPOCH tabulated (campaign 2022-10 to 2023-08).",
    "GCRT J1745-3009": "no ephemeris exists; period inferred from 3 burst intervals in a "
    "single 2002 observation, never re-detected.",
    "ASKAP J142431.2-612611": "epoch published (MJD 60684.856) but the period is quoted to "
    "1e-2 s where 4e-3 s is needed, and Pdot is unconstrained; most VAST epochs predate it.",
    "ASKAP J165130.3-450520": "no PEPOCH; reference epoch ~mid-2025 postdates most VAST epochs.",
    "ASKAP J170036.6-445758": "no PEPOCH; same as above.",
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--epochs", type=Path, default=Path("results/lptv_vast_epochs.csv"))
    ap.add_argument("--out", type=Path, default=Path("results/lptduty_phase.json"))
    args = ap.parse_args(argv)

    rows = ld.load_epochs(args.epochs)
    by: dict[str, list[ld.EpochRow]] = defaultdict(list)
    for r in rows:
        by[r.name].append(r)

    per_source = {}
    for name, eph in EPHEMERIDES.items():
        srows = by.get(name)
        if not srows:
            continue
        brightest = max(abs(r.v_mjy) for r in srows)
        dets = [r for r in srows if abs(r.v_mjy) / r.e_v >= ld.DETECT_THRESHOLD_SIGMA]
        assumed = max(abs(r.v_mjy) for r in dets) if dets else 5.0
        pr = ld.phase_resolved_activity(srows, eph, pulse_mjy=assumed)
        # Record where each detection actually lands in phase. This is the diagnostic that
        # tells a single-window model it is wrong, and it is how the machinery was validated:
        # J183950's two phases reproduce the values lptv published independently.
        det_phases = []
        for r in dets:
            ph = math.fmod(
                (r.epoch_mjd - eph.pepoch_mjd) * 86400.0 / eph.period_s
                + 0.5 * r.duration_s / eph.period_s,
                1.0,
            )
            if ph < 0:
                ph += 1.0
            det_phases.append({"obs_id": r.obs_id, "v_sigma": abs(r.v_mjy) / r.e_v, "phase": ph})
        # Pdot smear across THIS source's own epochs (max |t - PEPOCH|), from the published
        # bound in the catalogue; an earlier caveat quoted ~0.22 cycles for J1832 computed
        # over half the baseline -- the honest number uses the farthest epoch
        import csv as _csv

        pdot = next(
            (
                abs(float(r["pdot_s_s"]))
                for r in _csv.DictReader(open("data/lpt_sample.csv"))
                if r["name"] == name and r.get("pdot_s_s")
            ),
            None,
        )
        max_dt_days = max(abs(r.epoch_mjd - eph.pepoch_mjd) for r in srows)
        pdot_smear = None
        if pdot is not None:
            dt_s = max_dt_days * 86400.0
            pdot_smear = round(pdot * dt_s * dt_s / (2.0 * eph.period_s**2), 3)
        per_source[name] = {
            "pdot_bound_smear_cycles_at_farthest_epoch": pdot_smear,
            "window_assignment_defensible": bool(pdot_smear is None or pdot_smear <= 0.1),
            "ephemeris": {
                "period_s": eph.period_s,
                "sigma_period_s": eph.sigma_period_s,
                "pepoch_mjd": eph.pepoch_mjd,
                "pulse_width_s": eph.pulse_width_s,
                "reference": eph.reference,
            },
            "assumed_pulse_mjy": assumed,
            "brightest_abs_v_mjy": brightest,
            "n_epochs": len(srows),
            "n_in_pulse_window": pr.n_on_window,
            "effective_in_window": pr.effective_on_window,
            "window_fraction_of_period": pr.window_fraction,
            "n_detections_in_window": pr.n_detections_in_window,
            "n_detections_outside_window": pr.n_detections_outside,
            "f_active_point": pr.f_active_point,
            # the split is CONDITIONAL on the PEPOCH being a pulse arrival; the only support
            # is the detection's own phase, whose a-priori chance of landing in the window
            # is the window fraction itself
            "apriori_window_landing_prob": round(pr.window_fraction, 3),
            "f_active_upper_95": (
                None if math.isinf(pr.f_active_upper_95) else pr.f_active_upper_95
            ),
            "max_phase_uncertainty_cycles": pr.max_phase_uncertainty,
            "usable": pr.usable and bool(pdot_smear is None or pdot_smear <= 0.1),
            "pepoch_is_published_pulse_anchor": name != "ASKAP J183950.5-075635",
            "detection_phases": det_phases,
        }

    payload = {
        "slice": "lptduty",
        "stage": "phase-resolved",
        "is_real": True,
        "question": (
            "Where a published ephemeris allows physical phase, how much of the low "
            "per-snapshot detection rate is inactivity (f_active) rather than the snapshot "
            "missing a narrow pulse?"
        ),
        "caveats": [
            "Only sources with a PUBLISHED reference epoch are included. Assuming an epoch "
            "'near the campaign' would manufacture phase; most of the sample is excluded for "
            "that reason, not for period imprecision.",
            "Detections outside the predicted window are reported, not discarded: they are "
            "evidence the ephemeris is wrong.",
            "Period derivatives ARE folded in as pdot_bound_smear_cycles_at_farthest_epoch. "
            "For ASKAP J1832-0911 the published Pdot bound, if saturated, smears 0.40 cycles "
            "at the farthest epoch (an earlier caveat said ~0.22, computed over half the "
            "baseline), so its window assignment -- and hence n_in_pulse_window and any "
            "f_active bound -- is not defensible and usable is False.",
            "J183950's epoch is the repo's own anchor with a +/-0.013 cycle systematic, not a "
            "published PEPOCH.",
        ],
        "n_with_epoch_encoded": len(per_source),
        # J183950's anchor is the repo's own documented T0, not a published PEPOCH
        "n_published_pepoch": sum(
            1 for d in per_source.values() if d["pepoch_is_published_pulse_anchor"]
        ),
        # J142431 is excluded for PERIOD PRECISION, not for lacking an epoch: count it
        # separately so \ldNNoEpoch means what it says
        "n_no_published_epoch": sum(
            1 for v in NO_PUBLISHED_EPOCH.values() if "epoch published" not in v
        ),
        "n_epoch_but_excluded": sum(
            1 for v in NO_PUBLISHED_EPOCH.values() if "epoch published" in v
        ),
        "excluded_no_published_epoch": NO_PUBLISHED_EPOCH,
        "per_source": per_source,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(args.out)

    print(f"wrote {args.out}")
    for name, d in sorted(per_source.items()):
        fa = d["f_active_point"]
        lim = d["f_active_upper_95"]
        shown = (
            f"f_active={fa:.2f}"
            if fa is not None
            else (f"f_active < {lim:.2f}" if lim is not None else "f_active unconstrained")
        )
        print(
            f"  {name:30s} window {d['window_fraction_of_period']:.3f} of P, "
            f"{d['n_in_pulse_window']:3d}/{d['n_epochs']:3d} epochs on-window, "
            f"k_in={d['n_detections_in_window']} k_out={d['n_detections_outside_window']}, "
            f"{shown}, usable={d['usable']}"
        )
    print(
        f"  ({payload['n_no_published_epoch']} sources excluded: no published epoch; "
        f"{payload['n_epoch_but_excluded']} excluded for period precision)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
