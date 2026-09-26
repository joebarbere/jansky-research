"""Plan 64 candidate vetting in the images: does a point source walk along the fitted track?

Fetch cutouts first (radio-cutout skill, CADC SODA, 1.5 arcmin at each candidate's ra2/dec2 into
data/vlass/cutouts/c<N>/). Writes vetting.json + vetting_montage.png there. Two tests per epoch:
S/N of the peak within 2.5" of the predicted track position, and S/N within 2.5" of the FIRST
detection's position -- a real mover leaves that spot empty; a static source does not.
"""

import csv
import glob
import json
import sys

import matplotlib
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))
from jansky_research import vlasspm as v


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


cands = list(
    csv.DictReader(
        open(
            str(
                __import__("pathlib").Path(__file__).resolve().parents[1]
                / "results"
                / "vlasspm_candidates.csv"
            )
        )
    )
)
out = {}
fig, axes = plt.subplots(len(cands), 4, figsize=(11, 2.8 * len(cands)))
for n, c in enumerate(cands):
    ra1, dec1, t1 = float(c["ra1"]), float(c["dec1"]), float(c["t1"])
    mura, mudec = float(c["mu_ra"]), float(c["mu_dec"])
    files = sorted(
        glob.glob(
            f"/home/joe/dev/github/joebarbere/jansky-research/data/vlass/cutouts/c{n}/*.fits"
        ),
        key=lambda f: fits.getheader(f)["DATE-OBS"],
    )
    rows = []
    for k, f in enumerate(files):
        h = fits.getheader(f)
        d = np.squeeze(fits.getdata(f)).astype(float) * 1e3  # mJy/beam
        w = WCS(h).celestial
        t = Time(h["DATE-OBS"]).decimalyear
        rp, dp = v._mover_track(ra1, dec1, mura, mudec, t1, t)
        r0, d0 = ra1, dec1  # where it was at the first detection
        good = np.isfinite(d)
        rms = 1.4826 * np.median(np.abs(d[good] - np.median(d[good])))
        yy, xx = np.mgrid[0 : d.shape[0], 0 : d.shape[1]]

        snr_track, off_track = peak_near(d, w, rms, rp, dp)
        snr_first, _ = peak_near(d, w, rms, r0, dec1)
        rows.append(
            {
                "date": h["DATE-OBS"][:10],
                "rms_mjy": round(float(rms), 3),
                "snr_at_track": round(snr_track, 1),
                "offset_from_track_arcsec": round(off_track, 2),
                "snr_at_first_position": round(snr_first, 1),
            }
        )
        if k < 4:
            ax = axes[n, k]
            px, py = w.world_to_pixel_values(rp, dp)
            s = 20
            sub = d[int(py) - s : int(py) + s + 1, int(px) - s : int(px) + s + 1]
            ax.imshow(
                sub, origin="lower", cmap="gray_r", vmin=-2 * rms, vmax=max(6 * rms, np.nanmax(sub))
            )
            ax.plot(s, s, "r+", ms=10)  # predicted track position this epoch
            fx, fy = w.world_to_pixel_values(r0, dec1)
            ax.plot(fx - int(px) + s, fy - int(py) + s, "bx", ms=8)  # first-detection position
            ax.set_title(f"c{n} {h['DATE-OBS'][:10]}  S/N@track {snr_track:.1f}", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
    for k in range(len(files), 4):
        axes[n, k].axis("off")
    out[f"c{n}"] = {
        "triple": c["triple"],
        "ra": ra1,
        "dec": dec1,
        "mu": float(c["mu"]),
        "epochs": rows,
    }
fig.tight_layout()
fig.savefig(
    "/home/joe/dev/github/joebarbere/jansky-research/data/vlass/cutouts/vetting_montage.png", dpi=90
)
json.dump(
    out,
    open("/home/joe/dev/github/joebarbere/jansky-research/data/vlass/cutouts/vetting.json", "w"),
    indent=1,
)
for k, r in out.items():
    print(k, r["triple"], f"({r['ra']:.3f},{r['dec']:.3f}) mu={r['mu']:.2f}")
    for e in r["epochs"]:
        print("   ", e)
