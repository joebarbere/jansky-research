"""hiblend post-hoc diagnostics D0-D3 of the C2 failure (survey/hiblend-findings.md step 4).

    uv run python scripts/hiblend_diagnostics.py --out <dir>   # <dir>/results/hiblend_diagnostics.json

Post hoc: the predictions were committed (d35a374) before this ran. Nothing here can change the
frozen step-3 outcome; it only discriminates between the hypotheses for why C2 failed.
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
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--seed", type=int, default=97)
    args = ap.parse_args(argv)
    t0 = time.time()
    raw = fe.fetch_fashi_dr2()
    ok = (
        np.isfinite(raw["flux"]) & (raw["flux"] > 0) & np.isfinite(raw["flux_err"]) & (raw["flux_err"] > 0)
        & np.isfinite(raw["cz"]) & np.isfinite(raw["w50"]) & np.isfinite(raw["ra"]) & np.isfinite(raw["dec"])
    )  # fmt: skip
    fashi = {k: np.asarray(v)[ok] for k, v in raw.items()}
    field = h.build_field(fashi, h.fetch_alfalfa())
    diag = h.diagnostics(
        field, np.random.default_rng(args.seed), n_cat_ra=field["n_cat_ra"],
        n_cat_dec=field["n_cat_dec"], n_boot=args.n_boot,
    )  # fmt: skip
    metrics = {
        "source": "FASHI DR2 (arXiv:2606.31539) x ALFALFA alpha.100 (Haynes+2018, J/ApJ/861/49)",
        "is_real": True,
        "post_hoc": True,
        "predictions_commit": "d35a374",
        "seed": args.seed, "runtime_s": round(time.time() - t0, 1),
        **diag,
    }  # fmt: skip
    write_results(metrics, Path(args.out) / "results" / "hiblend_diagnostics.json")
    print(json.dumps(metrics, indent=1, default=str))


if __name__ == "__main__":
    main()
