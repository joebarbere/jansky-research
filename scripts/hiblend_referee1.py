"""hiblend note, referee round 1: post-hoc checks R1-R4 (survey/hiblend-findings.md step 6).

    uv run python scripts/hiblend_referee1.py --out .     # -> <out>/results/hiblend_referee1.json

Predictions were committed (9f49696) before this ran. Network: FASHI DR2 (cached) and VizieR.
--out is required, so a stray invocation cannot write into the repo root by default.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from jansky_research import fashienv as fe
from jansky_research import hiblend as h
from jansky_research.report import write_results


def main(argv: list[str] | None = None) -> None:  # pragma: no cover - network
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=101)
    args = ap.parse_args(argv)
    t0 = time.time()
    raw = fe.fetch_fashi_dr2()
    ok = (
        np.isfinite(raw["flux"]) & (raw["flux"] > 0) & np.isfinite(raw["flux_err"]) & (raw["flux_err"] > 0)
        & np.isfinite(raw["cz"]) & np.isfinite(raw["w50"]) & np.isfinite(raw["ra"]) & np.isfinite(raw["dec"])
    )  # fmt: skip
    fashi = {k: np.asarray(v)[ok] for k, v in raw.items()}
    field = h.build_field(fashi, h.fetch_alfalfa())
    res = {
        "source": "FASHI DR2 (arXiv:2606.31539) x ALFALFA alpha.100 (Haynes+2018, J/ApJ/861/49)",
        "is_real": True, "post_hoc": True, "predictions_commit": "9f49696",
        "seed": args.seed, **h.referee1_checks(field, seed=args.seed),
        "runtime_s": round(time.time() - t0, 1),
    }  # fmt: skip
    write_results(res, Path(args.out) / "results" / "hiblend_referee1.json")
    print(json.dumps({k: v for k, v in res.items() if not k.startswith("R1")}, indent=1)[:6000])


if __name__ == "__main__":
    main()
