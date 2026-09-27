"""Run the fashienv robustness leg (shuffle null + beam blending) and merge it into the results.

    uv run python scripts/fashienv_robustness.py --out <dir>   # reads/writes <dir>/results/...

The full real leg (``_real_leg``) computes the same ``robustness`` block. Needs the measured
offsets already in the results JSON (``void_knee_offset`` etc. and the ``env_vmax`` block).
Network: FASHI DR2, VizieR voids and groups. Writes through ``report.write_results``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from jansky_research import fashienv as fe
from jansky_research.report import write_results


def main(argv: list[str] | None = None) -> None:  # pragma: no cover - network
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=".", help="directory holding results/fashienv_metrics.json")
    ap.add_argument("--n-shuffle", type=int, default=1000)
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
    rob = fe.robustness_leg(
        cat, area, env, voids, grp, base, vcat, comp, fe.measured_offsets(m),
        n_rand=2_000_000, n_shuffle=args.n_shuffle,
    )  # fmt: skip
    m["robustness"] = rob
    write_results(m, path)
    print(json.dumps(rob, indent=1))


if __name__ == "__main__":
    main()
