#!/usr/bin/env python3
"""Plan 91 GATE 0: can a Taylor-term spectral index be measured for any LPT pulse we have?

Offline. Reads the committed VAST epoch table for the Stokes-I signal-to-noise actually
available at each detection, runs the injection-recovery study the plan requires *before* any
alpha claim, and writes results/lptspec_gate0.json.

The plan's own instruction: "if only one or two of the seven are recoverable, that is the
result". This script decides that, without fetching a single image.

    uv run python scripts/lptspec_gate0.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jansky_research import lptduty as ld  # noqa: E402
from jansky_research import lptspec as ls  # noqa: E402

# Indices spanning the physically interesting range: steep aged synchrotron plasma through
# flat/inverted coherent emission. Recovery is tested against each.
ALPHA_GRID = [-2.0, -0.7, 0.0, 1.0]
SNR_GRID = [300.0, 150.0, 100.0, 50.0, 36.0, 25.0, 15.0, 10.0]
TARGET_SIGMA_ALPHA = 0.3


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--epochs", type=Path, default=Path("results/lptv_vast_epochs.csv"))
    ap.add_argument("--racs", type=Path, default=Path("results/lptv_realtargets.csv"))
    ap.add_argument("--out", type=Path, default=Path("results/lptspec_gate0.json"))
    args = ap.parse_args(argv)

    # The pulses we actually have, and the Stokes-I S/N each offers a spectral index.
    pulses = []
    # Both legs of the lptv census: the RACS survey epochs and the VAST sweep. The plan
    # speaks of "seven detections"; they live in two committed tables.
    sources = [(args.epochs, "VAST"), (args.racs, "RACS")]
    for path, leg in sources:
        with open(path) as fh:
            for row in csv.DictReader(fh):
                try:
                    v, ev, i_f, ei = (
                        float(row["v_mjy"]),
                        float(row["e_v"]),
                        float(row["i_mjy"]),
                        float(row["e_i"]),
                    )
                except (TypeError, ValueError):
                    continue
                # NaN must be rejected explicitly: float("nan") does not raise, and every NaN
                # comparison is False, so a bare "< threshold" test lets NaN rows through as
                # detections. This silently admitted 319 non-detections on the first run.
                if not all(math.isfinite(x) for x in (v, ev, i_f, ei)):
                    continue
                if ev <= 0 or ei <= 0 or abs(v) / ev < ld.DETECT_THRESHOLD_SIGMA:
                    continue
                if abs(v) <= ld.LEAKAGE_FRAC * abs(i_f):
                    continue
                snr = i_f / ei
                sigma_alpha_ideal = ls.alpha_uncertainty(i_f, ei, 0.0)
                sigma_alpha = ls.realistic_sigma_alpha(i_f, ei, 0.0)
                pulses.append(
                    {
                        "leg": leg,
                        "name": row["name"],
                        "obs_id": row["obs_id"],
                        "i_mjy": i_f,
                        "snr_stokes_i": snr,
                        "sigma_alpha_idealised": sigma_alpha_ideal,
                        "sigma_alpha_realistic": sigma_alpha,
                        "usable": bool(sigma_alpha <= TARGET_SIGMA_ALPHA),
                    }
                )
    pulses.sort(key=lambda p: -p["snr_stokes_i"])

    recovery = {}
    for alpha in ALPHA_GRID:
        for snr in SNR_GRID:
            r = ls.injection_recovery(snr, alpha, n_trials=40000, seed=7)
            recovery[f"alpha={alpha:g},snr={snr:g}"] = {
                "alpha_true": r.alpha_true,
                "snr": r.snr,
                "alpha_mean": r.alpha_mean,
                "alpha_std": r.alpha_std,
                "bias": r.bias,
                "frac_within_0p3": r.frac_within_0p3,
            }

    n_usable = sum(1 for p in pulses if p["usable"])
    payload = {
        "slice": "lptspec",
        "stage": "gate0-recoverability",
        "is_real": True,
        "question": (
            "Before fetching any image: at the Stokes-I S/N our LPT pulses actually have, is "
            "a Taylor-term alpha recoverable at all?"
        ),
        "method": (
            "alpha = T1/T0 from ASKAP MFS Taylor terms. The Taylor-1 image carries "
            "sigma_1 = sigma_0 * sqrt(12) * nu0 / B (~10.7x for RACS-low), because the "
            "in-band lever arm is only +/-16% in frequency, so sigma_alpha ~ 10.7/(S/N). "
            "Injection-recovery measures bias as well as scatter, since a ratio with a noisy "
            "denominator is biased and an error bar cannot show that."
        ),
        "taylor1_noise_ratio": ls.taylor1_noise_ratio(),
        "mt_mfs_penalty": ls.MT_MFS_PENALTY,
        "mt_mfs_penalty_source": (
            "Rashid et al. 2024 (arXiv:2405.18978): MT-MFS in-band indices for SNR <~100 "
            "have errors >~0.2; the idealised formula here gives 0.107 there, so real MT-MFS "
            "is ~2x noisier. A uGMRT simulation, not an ASKAP calibration -- the real-cutout "
            "injection study must replace it with a measured ASKAP value."
        ),
        "snr_for_sigma_alpha_0p3_idealised": ls.usable_snr_threshold(0.3, realistic=False),
        "snr_for_sigma_alpha_0p3_realistic": ls.usable_snr_threshold(0.3),
        "snr_for_sigma_alpha_0p1_realistic": ls.usable_snr_threshold(0.1),
        "target_sigma_alpha": TARGET_SIGMA_ALPHA,
        "n_pulses": len(pulses),
        "n_pulses_usable": n_usable,
        "pulses": pulses,
        "injection_recovery": recovery,
        "caveats": [
            "This is the idealised case: Gaussian noise, no deconvolution error, no "
            "primary-beam or bandpass systematic, and the true source spectrum assumed to be "
            "a power law across the band. Real Taylor-term alpha will be no better than this "
            "and probably worse, so a source failing here cannot be rescued by real data.",
            "Bias is negligible above S/N ~20 with the T0 > 5 sigma guard, but grows below "
            "it, and near S/N ~6 the guard itself selects upward fluctuations of T0 and "
            "shrinks |alpha| -- a selection effect, not noise.",
            "sigma_alpha here is the statistical term only. Pulse dilution (the pulse "
            "occupies a fraction of the integration) is a separate systematic and is not "
            "folded in.",
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(args.out)

    print(f"wrote {args.out}")
    print(f"Taylor-1 noise penalty for RACS-low: {ls.taylor1_noise_ratio():.1f}x")
    print(
        f"S/N for sigma_alpha=0.3: {ls.usable_snr_threshold(0.3, realistic=False):.0f} "
        f"idealised, {ls.usable_snr_threshold(0.3):.0f} with the published MT-MFS penalty"
    )
    print(f"\n{'pulse':44s} {'S/N(I)':>7} {'ideal':>8} {'realistic':>10}  usable")
    for p in pulses:
        print(
            f"  {p['name']:28s} {p['obs_id']:11s} {p['snr_stokes_i']:7.1f} "
            f"{p['sigma_alpha_idealised']:8.2f} {p['sigma_alpha_realistic']:10.2f}  {p['usable']}"
        )
    print(
        f"\n{n_usable} of {len(pulses)} pulses can support an in-band index at "
        f"sigma_alpha <= {TARGET_SIGMA_ALPHA}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
