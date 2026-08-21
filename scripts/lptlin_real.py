#!/usr/bin/env python3
"""Plan 92: linear polarization of LPT pulses, from VAST Stokes Q/U images.

NETWORK. For each detected pulse whose observation carries Q and U, stages the
``taylor.0`` cutouts in i/q/u from CASDA, forced-measures all three at the catalogued
position, and writes results/lptlin_metrics.json.

Reads ``taylor.0`` only -- the band-averaged flux -- never ``taylor.1``. Plan 91 established
that the MFS spectral term is invalid for a source that varies within the synthesis; the
band-averaged Stokes images are exactly what `lptv` already forced-measures for I and V, and
dilution cancels in any ratio taken from the same image.

    CASDA_USERNAME=<opal email> uv run python scripts/lptlin_real.py
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

from jansky_research import lptlin as ll  # noqa: E402
from jansky_research import lptspec as ls  # noqa: E402

CUTOUT_RADIUS_DEG = 0.03


def forced_value(img: np.ndarray, wcs, ra: float, dec: float) -> tuple[float, float]:
    """Value at the catalogue pixel and the local rms from a surrounding annulus."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    x, y = wcs.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    xi, yi = int(round(float(x))), int(round(float(y)))
    if not (0 <= yi < img.shape[0] and 0 <= xi < img.shape[1]):
        raise ValueError("catalogue position falls outside the cutout")
    yy, xx = np.mgrid[: img.shape[0], : img.shape[1]]
    r = np.hypot(xx - xi, yy - yi)
    ann = img[(r > 12) & (r < 40) & np.isfinite(img)]
    rms = float(np.std(ann)) if ann.size > 50 else math.nan
    return float(img[yi, xi]), rms


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate0", type=Path, default=Path("results/lptspec_gate0.json"))
    ap.add_argument("--catalog", type=Path, default=Path("data/lpt_sample.csv"))
    ap.add_argument("--out", type=Path, default=Path("results/lptlin_metrics.json"))
    args = ap.parse_args(argv)

    pulses = json.loads(args.gate0.read_text())["pulses"]
    pos = {
        r["name"]: (float(r["ra_deg"]), float(r["dec_deg"]))
        for r in csv.DictReader(open(args.catalog))
    }

    casda = None
    out = []
    for p in pulses:
        name, obs = p["name"], p["obs_id"]
        sbid = obs.split("-")[-1]
        ra, dec = pos[name]
        print(f"--- {name} {obs} (SB{sbid})", flush=True)
        vals = {}
        ok = True
        for stokes in ("i", "q", "u"):
            res = ls.fetch_taylor_cutout(
                ra, dec, sbid, 0, stokes=stokes, radius_deg=CUTOUT_RADIUS_DEG, casda=casda
            )
            if res is None:
                print(f"    stokes {stokes}: no image in this observation", flush=True)
                ok = False
                break
            img, wcs, casda = res
            v, rms = forced_value(img, wcs, ra, dec)
            vals[stokes] = (v, rms)
            print(f"    stokes {stokes}: {v:9.3f} mJy/beam (local rms {rms:.3f})", flush=True)
        if not ok:
            out.append({**p, "measured": False, "reason": "no Q/U imaging in this observation"})
            continue

        i_v = vals["i"][0]
        q_v, u_v = vals["q"][0], vals["u"][0]
        sigma_qu = float(np.mean([vals["q"][1], vals["u"][1]]))
        try:
            lin = ll.linear_fraction(q_v, u_v, i_v, sigma_qu)
        except ValueError as exc:
            out.append({**p, "measured": False, "reason": str(exc)})
            continue
        out.append(
            {
                **p,
                "measured": True,
                "sbid": sbid,
                "i_mjy": i_v,
                "q_mjy": q_v,
                "u_mjy": u_v,
                "sigma_qu_mjy": sigma_qu,
                "l_raw_mjy": lin.l_raw_mjy,
                "l_debiased_mjy": lin.l_debiased_mjy,
                "linear_fraction": lin.linear_fraction,
                "evpa_deg": lin.angle_deg,
                "significance": lin.significance,
                "above_leakage": lin.above_leakage,
                "detected": lin.detected,
            }
        )
        print(
            f"    L={lin.l_debiased_mjy:.3f} mJy ({lin.significance:.1f} sigma), "
            f"L/I={100 * lin.linear_fraction:.1f}%, EVPA={lin.angle_deg:.0f} deg, "
            f"detected={lin.detected}",
            flush=True,
        )

    payload = {
        "slice": "lptlin",
        "stage": "real-linear-polarization",
        "is_real": True,
        "source": "VAST/RACS taylor.0 Stokes i/q/u cutouts staged from CASDA (SODA)",
        "method": (
            "Forced value at the catalogue pixel in each Stokes image; L = sqrt(Q^2+U^2) "
            "debiased by subtracting the noise in quadrature; L/I with a leakage veto at "
            "0.6% of |I|, the same two-part criterion lptv applies to V. taylor.0 only -- "
            "the MFS spectral term is invalid for a transient (plan 91)."
        ),
        "caveats": [
            "The leakage floor is assumed equal to lptv's I->V figure (0.6%); a measured "
            "per-field I->Q/U leakage would be better and is not available here.",
            "Dilution cancels in L/I only while the pulse is well above the noise; the "
            "debiasing term biases the fraction slightly DOWN as the pulse approaches the "
            "noise floor (2% low at 5 sigma). Down cannot manufacture polarization.",
            "No ionospheric Faraday rotation correction: the EVPA is the observed angle, not "
            "the intrinsic one, and RM is not measured from a single band-averaged image.",
            "A single epoch per pulse; no independent re-imaging check.",
        ],
        "n_attempted": len(pulses),
        "n_measured": sum(1 for r in out if r.get("measured")),
        "n_detected": sum(1 for r in out if r.get("detected")),
        "pulses": out,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(args.out)
    print(
        f"\nwrote {args.out}: {payload['n_measured']} measured, "
        f"{payload['n_detected']} linear detections"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
