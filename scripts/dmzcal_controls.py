"""Plan 99 step 2: C0 (null false-fail rate) and C1 (power) for the calibration rule.

Synthetic by design: both controls draw from E2's own generative model, so this file
never touches a real burst. Writes ``results/dmzcal_controls.json``; it must be committed
before any real posterior is computed (plan 99, "Statistics and controls").

* ``frozen``: exactly the plan -- generator and scorer both E2 with its comoving-volume
  prior on 0.01 <= z <= 4, N = 36 (the step-0 certification size), M = 2000.
* ``supplementary_zmax1``: POST HOC, non-gating. Added after a preliminary M = 200 run
  showed C1 host power ~0.1 under the frozen prior, whose median z ~2.4 is far above the
  certification sample's (mostly z < 0.5), where host DM matters most. Same rule
  construction with the prior truncated at z <= 1 for generator and scorer. It answers
  "how much power would an estimator whose prior sits at low z have" -- relevant to
  reading E1 (zdm, survey-selected prior) -- and does not change any frozen verdict.

Run: ``uv run python scripts/dmzcal_controls.py``  (~10 min, CPU)
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from jansky_research import dmzcal

ROOT = Path(__file__).resolve().parents[1]
N_CERT = 36  # results/dmzcal_provenance.json counts.certification_spec_z


def main() -> None:
    prov = json.loads((ROOT / "results" / "dmzcal_provenance.json").read_text())
    n = int(prov["counts"]["certification_spec_z"])
    assert n == N_CERT, f"certification size changed: {n}"
    t0 = time.time()
    frozen = dmzcal.run_controls(n=n, m=2000, seed=99)
    t1 = time.time()
    supp = dmzcal.run_controls(n=n, m=2000, seed=990, table_kw={"zmax": 1.0})
    t2 = time.time()
    out = {
        "source": "synthetic controls: plan 99 C0 (null false-fail) and C1 (power), "
        "E2 generative model only -- no real burst used",
        "plan": "plans/99-dmzcal-coverage.md step 2",
        "frozen": frozen,
        "supplementary_zmax1": {
            "post_hoc": True,
            "gating": False,
            "why": "preliminary M=200 frozen run gave C1 host power ~0.1; the frozen "
            "comoving-volume prior has median z ~2.4, unlike the certification sample",
            **supp,
        },
        "runtime_s": {"frozen": round(t1 - t0, 1), "supplementary": round(t2 - t1, 1)},
    }
    path = ROOT / "results" / "dmzcal_controls.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    for k in ("frozen", "supplementary_zmax1"):
        r = out[k]
        print(
            k,
            "C0 ff",
            r["C0"]["false_fail_rate"],
            "widened",
            r["rule"]["widened"],
            "| C1 host",
            r["C1"]["host_mean_x2"]["power"],
            "F",
            r["C1"]["F_x2"]["power"],
        )


if __name__ == "__main__":
    main()
