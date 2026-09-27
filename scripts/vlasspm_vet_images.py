"""Plan 64 candidate vetting in the images: does a point source walk along the fitted track?

Fetch cutouts first (radio-cutout skill, CADC SODA, 1.5 arcmin at each candidate's ra2/dec2) into
``<cutouts>/<key>/`` where ``<key> = candidate_key(row)`` (``J<ra2>_<dec2>`` to 4 decimals), so a
candidate keeps its cutouts across re-runs whatever its index in the CSV. ``--fetch`` does that
for any candidate without a directory. Writes the per-epoch measurements to ``--out-json`` and a
montage to ``--montage``. Two tests per epoch: S/N of the peak within 2.5" of the predicted track
position, and S/N within 2.5" of the FIRST detection's position -- a real mover leaves that spot
empty; a static source does not. The restoring beam of every cutout (header BMAJ/BMIN/BPA) is
recorded, because a beam that changes shape between epochs moves the centroid of a partly
resolved source (survey/vlasspm-findings.md).
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import subprocess
import sys
from pathlib import Path

import matplotlib
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from jansky_research import vlasspm as v  # noqa: E402


def candidate_key(c: dict) -> str:
    return f"J{float(c['ra2']):.4f}_{float(c['dec2']):+.4f}"


def peak_near(d, w, rms, ra, dec, rad_arcsec=2.5):
    """S/N of the brightest pixel within ``rad_arcsec`` of (ra, dec), and its offset (arcsec)."""
    good = np.isfinite(d)
    yy, xx = np.mgrid[0 : d.shape[0], 0 : d.shape[1]]
    px, py = w.world_to_pixel_values(ra, dec)
    sel = (np.hypot(xx - px, yy - py) <= rad_arcsec) & good  # 1"/pix cutouts
    if not sel.any():
        return np.nan, np.nan
    iy, ix = np.unravel_index(np.argmax(np.where(sel, d, -np.inf)), d.shape)
    pr, pd = w.pixel_to_world_values(ix, iy)
    off = np.hypot(*v.tangent_offsets_arcsec(ra, dec, pr, pd))
    return float(d[iy, ix] / rms), float(off)


def fetch(c: dict, dest: Path) -> None:  # network
    dest.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(REPO / ".claude/skills/radio-cutout/fetch_cutout.py"),
        "--survey",
        "vlass",
        "--ra",
        str(c["ra2"]),
        "--dec",
        str(c["dec2"]),
        "--size-arcmin",
        "1.5",
        "--out",
        str(dest),
    ]
    with (dest / "fetch.log").open("w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default=str(REPO / "results" / "vlasspm_candidates.csv"))
    ap.add_argument("--cutouts", default=str(REPO / "data" / "vlass" / "cutouts"))
    ap.add_argument("--out-json", default=str(REPO / "data" / "vlass" / "cutouts" / "vetting.json"))
    ap.add_argument("--montage", default=str(REPO / "data" / "vlass" / "cutouts" / "montage.png"))
    ap.add_argument("--fetch", action="store_true", help="fetch cutouts that are missing")
    args = ap.parse_args()
    cands = list(csv.DictReader(open(args.candidates)))
    root = Path(args.cutouts)
    out = {}
    fig, axes = plt.subplots(len(cands), 4, figsize=(11, 2.8 * len(cands)), squeeze=False)
    for n, c in enumerate(cands):
        key = candidate_key(c)
        cdir = root / key
        if args.fetch and not list(cdir.glob("*.fits")):
            print(f"fetching {key} ...", flush=True)
            fetch(c, cdir)
        ra1, dec1, t1 = float(c["ra1"]), float(c["dec1"]), float(c["t1"])
        mura, mudec = float(c["mu_ra"]), float(c["mu_dec"])
        files = sorted(glob.glob(str(cdir / "*.fits")), key=lambda f: fits.getheader(f)["DATE-OBS"])
        rows = []
        for k, f in enumerate(files):
            h = fits.getheader(f)
            d = np.squeeze(fits.getdata(f)).astype(float) * 1e3  # mJy/beam
            w = WCS(h).celestial
            t = Time(h["DATE-OBS"]).decimalyear
            rp, dp = v._mover_track(ra1, dec1, mura, mudec, t1, t)
            good = np.isfinite(d)
            rms = 1.4826 * np.median(np.abs(d[good] - np.median(d[good])))
            snr_track, off_track = peak_near(d, w, rms, rp, dp)
            snr_first, _ = peak_near(d, w, rms, ra1, dec1)
            beam = [h.get("BMAJ"), h.get("BMIN"), h.get("BPA")]
            rows.append(
                {
                    "date": h["DATE-OBS"][:10],
                    "rms_mjy": round(float(rms), 3),
                    "snr_at_track": round(snr_track, 1),
                    "offset_from_track_arcsec": round(off_track, 2),
                    "snr_at_first_position": round(snr_first, 1),
                    "beam_arcsec_deg": [
                        round(beam[0] * 3600, 2) if beam[0] else None,
                        round(beam[1] * 3600, 2) if beam[1] else None,
                        round(beam[2], 1) if beam[2] is not None else None,
                    ],
                }
            )
            if k < 4:
                ax = axes[n, k]
                px, py = w.world_to_pixel_values(rp, dp)
                s = 20
                sub = d[int(py) - s : int(py) + s + 1, int(px) - s : int(px) + s + 1]
                ax.imshow(
                    sub,
                    origin="lower",
                    cmap="gray_r",
                    vmin=-2 * rms,
                    vmax=max(6 * rms, np.nanmax(sub)),
                )
                ax.plot(s, s, "r+", ms=10)  # predicted track position this epoch
                fx, fy = w.world_to_pixel_values(ra1, dec1)
                ax.plot(fx - int(px) + s, fy - int(py) + s, "bx", ms=8)  # first detection
                ax.set_title(f"{key} {h['DATE-OBS'][:10]}  S/N@track {snr_track:.1f}", fontsize=7)
                ax.set_xticks([])
                ax.set_yticks([])
        for k in range(len(files), 4):
            axes[n, k].axis("off")
        out[key] = {
            "index": n,
            "triple": c["triple"],
            "ra": ra1,
            "dec": dec1,
            "mu": float(c["mu"]),
            "resid_sigma": float(c["resid_sigma"]),
            "epochs": rows,
        }
    fig.tight_layout()
    fig.savefig(args.montage, dpi=90)
    Path(args.out_json).write_text(json.dumps(out, indent=1))
    for k, r in out.items():
        print(k, r["triple"], f"({r['ra']:.3f},{r['dec']:.3f}) mu={r['mu']:.2f}")
        for e in r["epochs"]:
            print("   ", e)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
