"""Plan 99, GATE-2 round 2 (N1): power of the post-hoc DM|z calibration check.

Synthetic, POST HOC. DM_EG is drawn at the 36 certification redshifts (the only real
input: z values, no DMs) from a mis-specified E2 world and scored with ``dmzcal.dm_pit``
under the unmodified model; power = fraction of samples KS-rejected at p < 0.01. It
says how large a host-DM / scatter error the "no DM-model failure detected at known z" statement can
exclude. Writes ``results/dmzcal_dmz_power.json``.

Run: ``uv run python scripts/dmzcal_dmz_power.py``
"""

from __future__ import annotations

import json
import zlib
from pathlib import Path

import numpy as np
from scipy import stats

from jansky_research import dmzcal

ROOT = Path(__file__).resolve().parents[1]
REPS = 1000


def main() -> None:
    prov = json.loads((ROOT / "results" / "dmzcal_provenance.json").read_text())
    z = np.array(
        [
            r["z"]
            for r in prov["bursts"]
            if r["side"] == "certification" and r["spec_z"] and r["secure_host"]
        ]
    )
    base = dmzcal.HOFFMANN_EMIN25
    cases = {
        "null": base,
        "host_mean_x1.5": dmzcal.E2Params(lmean=base.lmean + np.log10(1.5)),
        "host_mean_x2": dmzcal.E2Params(lmean=base.lmean + np.log10(2.0)),
        "host_mean_div1.5": dmzcal.E2Params(lmean=base.lmean - np.log10(1.5)),
        "F_x2": dmzcal.E2Params(F=2 * base.F),
        "MW_DM_plus30": None,  # handled by shifting the drawn DM
    }
    score = dmzcal.build_table(base)
    out = {}
    for name, gp in cases.items():
        rng = np.random.default_rng(zlib.crc32(name.encode()))
        if gp is None:
            # an unmodelled +30 pc/cc (e.g. DM_ISM or halo under-subtracted): draw from the
            # true model and add 30 before scoring
            gen = score
            i = np.abs(gen.z[None, :] - z[:, None]).argmin(axis=1)
            cum = np.cumsum(gen.prob[i], axis=1)
            cum /= cum[:, -1:]
            rej = 0
            for _ in range(REPS):
                k = np.minimum(
                    (cum < rng.random(z.size)[:, None]).sum(axis=1), gen.prob.shape[1] - 1
                )
                dm = gen.edges[k] + rng.random(z.size) * gen.dm_step + 30.0
                p = dmzcal.dm_pit(score, dm, gen.z[i])
                rej += float(stats.kstest(p, "uniform").pvalue) < 0.01
            out[name] = rej / REPS
        else:
            gen = score if name == "null" else dmzcal.build_table(gp)
            out[name] = dmzcal.dm_given_z_power(score, gen, z, REPS, rng)
        print(name, out[name])
    res = {
        "source": "synthetic power analysis (post hoc, GATE-2 round 2 N1): DM|z check at "
        "the 36 certification redshifts; no real DM used",
        "plan": "plans/99-dmzcal-coverage.md (post hoc)",
        "n": int(z.size),
        "reps": REPS,
        "alpha": 0.01,
        "rejection_rate": out,
    }
    (ROOT / "results" / "dmzcal_dmz_power.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
