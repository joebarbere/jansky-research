"""Run only the environment-restricted-Vmax leg of fashienv and merge it into the results JSON.

The full real leg (``python -m jansky_research.fashienv --real``) computes the same ``env_vmax``
block; this script exists so the referee-round checks can be added without rerunning the 1,000-
placement nulls, which they do not change. Network: FASHI DR2, VizieR voids and groups.

    uv run python scripts/fashienv_env_vmax.py --out <dir>   # reads/writes <dir>/results/...

``--out`` defaults to the repo root; point it at a scratch copy to experiment. The write goes
through ``report.write_results`` (``preserve_live_results``), and the previous ``env_vmax``
block's summaries are kept under ``env_vmax_previous`` so a rerun can be compared with the last.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from jansky_research import fashienv as fe
from jansky_research.report import write_results


def _summary(block: dict) -> dict:
    """An env_vmax block without its per-placement rows."""
    out = {}
    for k, v in block.items():
        if isinstance(v, dict) and "rows" in v:
            out[k] = {a: b for a, b in v.items() if a != "rows"}
        elif k != "previous":
            out[k] = v
    return out


def main(argv: list[str] | None = None) -> None:  # pragma: no cover - network
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=".", help="directory holding results/fashienv_metrics.json")
    ap.add_argument("--no-frame-correct", action="store_true", help="skip the dmax frame fix")
    args = ap.parse_args(argv)
    path = Path(args.out) / "results" / "fashienv_metrics.json"
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
        frame_correct=not args.no_frame_correct,
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
    if "env_vmax" in m:
        ev["previous"] = _summary(m["env_vmax"])
    m["env_vmax"] = ev
    write_results(m, path)
    print(json.dumps(_summary(ev), indent=1))


if __name__ == "__main__":
    main()
