"""Run only the environment-restricted-Vmax leg of fashienv and merge it into the results JSON.

The full real leg (``python -m jansky_research.fashienv --real``) computes the same ``env_vmax``
block; this script exists so the fourth-round check can be added without rerunning the 1,000-
placement nulls, which it does not change. Network: FASHI DR2, VizieR voids and groups.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from jansky_research import fashienv as fe
from jansky_research.report import write_results


def main() -> None:  # pragma: no cover - network
    path = Path("results/fashienv_metrics.json")
    m = json.loads(path.read_text())
    voids = fe.fetch_voidfinder_spheres(all_holes=True)
    grp = fe.fetch_tempel_groups()
    cat = fe._clean(fe.fetch_fashi_dr2())
    env = fe._environments(cat, voids, grp)
    comp, vcat = cat["completeness"], cat["vmax_mpc3"]
    base = np.isfinite(comp) & (comp >= fe.FASHI_DR2_C_MIN) & np.isfinite(vcat) & (vcat > 0)
    area = fe.FASHI_DR2_AREA_DEG2 * (np.pi / 180.0) ** 2
    ev = fe.env_vmax_leg(
        cat, area, env, voids, grp, base, vcat, comp, (cat["ra"], cat["dec"]),
        (grp["gal_ra"], grp["gal_dec"]), n_rand=2_000_000, n_void=200, n_group=200,
    )  # fmt: skip
    # Pairing check: the survey-wide-Vmax offsets must reproduce the committed null rows.
    old = [r["B_offset"] for r in m["random_void_null_constrained"]["rows"][:200]]
    new = [r["survey_vmax"] for r in ev["void_null_constrained"]["rows"]]
    ev["pairing_check_void_max_abs_diff"] = round(
        float(np.nanmax(np.abs(np.subtract(old, new)))), 5
    )
    oldg = [r["offset"] for r in m["random_group_null"]["rows"][:200]]
    newg = [r["survey_vmax"] for r in ev["group_null"]["rows"]]
    ev["pairing_check_group_max_abs_diff"] = round(
        float(np.nanmax(np.abs(np.subtract(oldg, newg)))), 5
    )
    m["env_vmax"] = ev
    write_results(m, path)
    print(
        json.dumps(
            {k: v for k, v in ev.items() if k not in ("void_null_constrained", "group_null")},
            indent=1,
        )
    )
    for k in ("void_null_constrained", "group_null"):
        print(k, json.dumps({a: b for a, b in ev[k].items() if a != "rows"}))


if __name__ == "__main__":
    main()
