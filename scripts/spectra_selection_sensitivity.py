"""Post-hoc sensitivity of the spectra selection-bias model to its one assumption (2026-10-11).

`selection_bias_mc` adds our formal index errors to the reference indices. This script reruns it
on the committed `results/uss_reference_check.csv` (single-component reference rows) with the
noise of the selected variable scaled to what the raw-minus-reference differences allow: two
variance bounds from the standard deviation, and two from the robust (MAD) core, since the
differences are heavy-tailed on the flat side while USS selection acts on the steep side. Writes
`results/uss_selection_sensitivity.json`; reads only committed files.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from jansky_research.spectra import selection_bias_mc

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    rows = [
        r
        for r in csv.DictReader((ROOT / "results" / "uss_reference_check.csv").open())
        if r["ref_scode"] == "S"
    ]
    ref = np.array([float(r["ref_alpha"]) for r in rows])
    e_ours = np.array([float(r["e_alpha"]) for r in rows])
    e_ref = np.array([float(r["ref_e"]) for r in rows])
    diff = np.array([float(r["alpha_raw"]) - float(r["ref_alpha"]) for r in rows])
    obs_var = float(diff.var(ddof=1))
    # Lower bound: the reference's formal errors are right, so our share of the scatter is the rest.
    k_lo = float(np.sqrt(max(obs_var - np.mean(e_ref**2), 0.0) / np.mean(e_ours**2)))
    # Upper bound: all of the observed scatter is in our (selected) index.
    k_hi = float(np.sqrt(obs_var / np.mean(e_ours**2)))
    comb = np.hypot(e_ours, e_ref)
    mad = lambda x: float(1.4826 * np.median(np.abs(x - np.median(x))))  # noqa: E731
    # Robust core: the differences' MAD in units of the combined formal error (applied to both
    # catalogues alike), and the MAD scatter attributed entirely to our index.
    k_core = mad(diff / comb)
    k_core_ours = mad(diff) / float(np.sqrt(np.mean(e_ours**2)))
    cand = [d for d, r in zip(diff, rows, strict=True) if float(r["alpha_raw"]) < -1.3]
    obs = float(np.mean(cand))
    se = float(np.std(cand, ddof=1) / np.sqrt(len(cand)))
    out = {
        "source": "real: tgss-x-nvss reference comparison (post-hoc, from committed CSV)",
        "n_rows": len(rows),
        "rms_e_ours": round(float(np.sqrt(np.mean(e_ours**2))), 4),
        "rms_e_ref": round(float(np.sqrt(np.mean(e_ref**2))), 4),
        "observed_scatter": round(float(np.sqrt(obs_var)), 4),
        "observed_scatter_mad": round(mad(diff), 4),
        "z_sd": round(float(np.std(diff / comb, ddof=1)), 3),
        "z_mad": round(k_core, 3),
        "n_abs_diff_gt_0p2": int(np.sum(np.abs(diff) > 0.2)),
        "n_diff_gt_0p2": int(np.sum(diff > 0.2)),
        "candidate_offset": round(obs, 4),
        "candidate_offset_se": round(se, 4),
        "variants": {},
    }
    variants = (
        ("core_formal_ratio", k_core),
        ("core_all_ours", k_core_ours),
        ("sd_ref_errors_trusted", k_lo),
        ("formal", 1.0),
        ("sd_all_scatter_ours", k_hi),
    )
    for name, k in variants:
        m = selection_bias_mc(ref, e_ours, noise_scale=k)
        resid = (obs - m["predicted_offset"]) / np.hypot(se, m["predicted_offset_sd"])
        out["variants"][name] = {
            "noise_scale": round(k, 3),
            **m,
            "fraction_explained": round(m["predicted_offset"] / obs, 2),
            "residual_sigma": round(float(resid), 2),
        }
    path = ROOT / "results" / "uss_selection_sensitivity.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["variants"], indent=1))


if __name__ == "__main__":
    main()
