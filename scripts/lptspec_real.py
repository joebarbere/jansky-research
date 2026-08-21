#!/usr/bin/env python3
"""Plan 91: measure in-band alpha for the LPT pulses GATE 0 said could carry one.

NETWORK. Stages RACS/VAST ``taylor.0`` and ``taylor.1`` cutouts from CASDA (OPAL login) for
the three pulses `results/lptspec_gate0.json` flagged usable, forced-measures both at the
catalogued position, and writes results/lptspec_metrics.json.

Forced, not peak-searched: the position is fixed at the catalogue coordinate, as in the
`stokesv` slice, because a peak search over blank sky returns the largest of several
independent beams and biases every value high.

    CASDA_USERNAME=<opal email> uv run python scripts/lptspec_real.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np  # noqa: E402

from jansky_research import lptspec as ls  # noqa: E402

CUTOUT_RADIUS_DEG = 0.03


def forced_value(img: np.ndarray, wcs, ra: float, dec: float) -> tuple[float, float]:
    """Value at the catalogue pixel, and the local rms from an annulus around it."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    x, y = wcs.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    xi, yi = int(round(float(x))), int(round(float(y)))
    if not (0 <= yi < img.shape[0] and 0 <= xi < img.shape[1]):
        raise ValueError("catalogue position falls outside the cutout")
    value = float(img[yi, xi])
    yy, xx = np.mgrid[: img.shape[0], : img.shape[1]]
    r = np.hypot(xx - xi, yy - yi)
    ann = img[(r > 12) & (r < 40) & np.isfinite(img)]
    rms = float(np.std(ann)) if ann.size > 50 else math.nan
    return value, rms


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate0", type=Path, default=Path("results/lptspec_gate0.json"))
    ap.add_argument("--catalog", type=Path, default=Path("data/lpt_sample.csv"))
    ap.add_argument("--out", type=Path, default=Path("results/lptspec_metrics.json"))
    args = ap.parse_args(argv)

    gate0 = json.loads(args.gate0.read_text())
    usable = [p for p in gate0["pulses"] if p["usable"]]
    pos = {
        r["name"]: (float(r["ra_deg"]), float(r["dec_deg"]))
        for r in csv.DictReader(open(args.catalog))
    }

    casda = None
    out = []
    for p in usable:
        name, obs = p["name"], p["obs_id"]
        sbid = obs.split("-")[-1]
        ra, dec = pos[name]
        print(f"--- {name} {obs} (SB{sbid})", flush=True)
        got = {}
        for term in (0, 1):
            res = ls.fetch_taylor_cutout(
                ra, dec, sbid, term, radius_deg=CUTOUT_RADIUS_DEG, casda=casda
            )
            if res is None:
                print(f"    taylor.{term}: FETCH FAILED", flush=True)
                break
            img, wcs, casda = res
            val, rms = forced_value(img, wcs, ra, dec)
            got[term] = {"value_mjy": val, "local_rms_mjy": rms}
            print(f"    taylor.{term}: {val:9.3f} mJy/beam (local rms {rms:.3f})", flush=True)
        if len(got) != 2:
            out.append({**p, "measured": False, "reason": "cutout fetch failed"})
            continue

        t0, s0 = got[0]["value_mjy"], got[0]["local_rms_mjy"]
        t1 = got[1]["value_mjy"]
        try:
            alpha = ls.taylor_alpha(t0, t1)
        except ValueError as exc:
            out.append({**p, "measured": False, "reason": str(exc)})
            continue
        sigma_ideal = ls.alpha_uncertainty(abs(t0), s0, t1) if s0 > 0 else math.nan
        sigma_real = ls.MT_MFS_PENALTY * sigma_ideal if math.isfinite(sigma_ideal) else math.nan
        out.append(
            {
                **p,
                "measured": True,
                "sbid": sbid,
                "taylor0_mjy": t0,
                "taylor1_mjy": t1,
                "local_rms_mjy": s0,
                "measured_snr_taylor0": t0 / s0 if s0 > 0 else math.nan,
                "alpha": alpha,
                "sigma_alpha_idealised": sigma_ideal,
                "sigma_alpha_realistic": sigma_real,
            }
        )
        print(f"    alpha = {alpha:+.2f} +/- {sigma_real:.2f} (realistic)", flush=True)

    payload = {
        "slice": "lptspec",
        "stage": "real-taylor-alpha",
        "is_real": True,
        "source": "RACS/VAST taylor.0 and taylor.1 cutouts staged from CASDA (SODA)",
        "method": (
            "Forced value at the catalogue pixel in each Taylor term, local rms from an "
            "annulus; alpha = T1/T0. Forced rather than peak-searched, per the stokesv "
            "lesson that a peak search over blank sky is biased high."
        ),
        "caveats": [
            "sigma_alpha uses the local rms of the taylor.0 cutout scaled by the analytic "
            "Taylor-1 penalty and the published MT-MFS factor; it is NOT yet an ASKAP-"
            "measured injection result. Replacing it is the next increment.",
            "Pulse dilution is not corrected: taylor.1 is fitted over the whole synthesis "
            "while the pulse occupies part of it, so alpha here is a synthesis-averaged "
            "quantity, not the intrinsic in-pulse index.",
            "A single cutout per term; no independent re-imaging or self-calibration check.",
        ],
        "n_attempted": len(usable),
        "n_measured": sum(1 for r in out if r.get("measured")),
        "pulses": out,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(args.out)
    print(f"\nwrote {args.out}: {payload['n_measured']} of {len(usable)} measured")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
