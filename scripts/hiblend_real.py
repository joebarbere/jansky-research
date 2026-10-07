"""hiblend real leg (plan 97): FASHI DR2 x ALFALFA alpha.100, gated controls, then beta.

    uv run python scripts/hiblend_real.py --out <dir>      # writes <dir>/results/hiblend_metrics.json

The controls run in the order frozen in plans/97-hiblend-fashi-alfalfa.md and a failed gate
stops the run before beta is computed (see hiblend.run_gated). Network: the FASHI DR2 CSTCloud
share (cached under data/fashi_dr2/) and VizieR (CfA mirror fallback).
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
    ap.add_argument("--out", default=".")
    ap.add_argument("--n-power", type=int, default=20)
    ap.add_argument("--n-inj", type=int, default=10)
    ap.add_argument("--seed", type=int, default=97)
    args = ap.parse_args(argv)
    t0 = time.time()
    raw = fe.fetch_fashi_dr2()
    ok = (
        np.isfinite(raw["flux"]) & (raw["flux"] > 0) & np.isfinite(raw["flux_err"]) & (raw["flux_err"] > 0)
        & np.isfinite(raw["cz"]) & np.isfinite(raw["w50"]) & np.isfinite(raw["ra"]) & np.isfinite(raw["dec"])
    )  # fmt: skip
    fashi = {k: np.asarray(v)[ok] for k, v in raw.items()}
    alf = h.fetch_alfalfa()
    field = h.build_field(fashi, alf)
    shifted = h.build_field(fashi, alf, shift_dec_arcmin=h.SHIFT_ARCMIN)
    print(f"[hiblend] {len(fashi['ra'])} FASHI, {len(alf['ra'])} ALFALFA, {len(field['flux_f'])} matched "
          f"({time.time() - t0:.0f}s)", flush=True)  # fmt: skip
    res = h.run_gated(
        field, shifted, np.random.default_rng(args.seed), n_power=args.n_power, n_inj=args.n_inj
    )
    metrics = {
        "source": "FASHI DR2 (arXiv:2606.31539) x ALFALFA alpha.100 (Haynes+2018, J/ApJ/861/49)",
        "is_real": True,
        "plan": "plans/97-hiblend-fashi-alfalfa.md (controls frozen 2026-10-05)",
        "beams_arcmin": {"fast": h.FAST_FWHM_ARCMIN, "alfa": h.ALFA_FWHM_ARCMIN},
        "match": {"radius_arcmin": h.MATCH_RADIUS_ARCMIN, "dv_kms": h.MATCH_DV_KMS},
        "n_fashi": int(len(fashi["ra"])), "n_alfalfa": int(len(alf["ra"])),
        "n_alfalfa_code1": field["n_alfalfa_code1"],
        "n_alfalfa_only_neighbours": field["n_alfalfa_only_neighbours"],
        "seed": args.seed, "runtime_s": round(time.time() - t0, 1),
        **res,
    }  # fmt: skip
    path = Path(args.out) / "results" / "hiblend_metrics.json"
    write_results(metrics, path)
    print(json.dumps({k: v for k, v in metrics.items() if k != "gates"}, indent=1, default=str))
    print(json.dumps(metrics["gates"], indent=1, default=str))


if __name__ == "__main__":
    main()
