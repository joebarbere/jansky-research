#!/usr/bin/env python3
"""Fetch a radio-survey FITS image cutout at a sky position, routing to the archive that serves it.

One place that knows *which archive serves which survey* and how to pull a cutout — a recurring need
across the research slices. The **CADC SODA** path (VLASS and any other CADC collection) is the
validated, default workhorse; other surveys are routed to their archive with honest notes.

  Archive map (verified during the slices):
    VLASS (2-4 GHz)            -> CADC SODA           [working — this script]
    GLEAM / GLEAM-X (72-231)   -> Data Central        [POST cutout API; see --survey gleam notes]
    RACS-I (low/mid/high)      -> Data Central / CASDA [Stokes I]
    RACS-V (Stokes V)          -> CASDA only          [use the casda-cutout-fetch skill; often down]

Exit codes: 0 ok (FITS saved + validated) · 2 bad usage · 5 no image found · 6 unsupported survey
route (not yet wired here — see notes) · 7 download/validation failure · 8 unexpected error.

See SKILL.md.
"""

from __future__ import annotations

import argparse
import io
import re
from pathlib import Path

FITS_MAGIC = b"SIMPLE  ="


def log(msg: str) -> None:
    print(f"[radio-cutout] {msg}", flush=True)


def fetch_cadc_cutout(
    collection: str, ra: float, dec: float, size_arcmin: float, *, name_filter: str | None = None
):
    """Download a SODA cutout from a CADC collection at (ra, dec). Returns (path-bytes, header) tuples.

    Generalises the VLASS cutout path: one ``query_region`` + ``get_image_list`` round-trip, then the
    SODA cutout of each matching image. ``name_filter`` is a regex to pick image variants (e.g.
    ``\\.ql\\.`` for VLASS Quick-Look). Yields ``(filename, fits_bytes)``.
    """
    import requests
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.cadc import Cadc

    c = SkyCoord(ra * u.deg, dec * u.deg)
    rad = (size_arcmin / 60.0) * u.deg
    cadc = Cadc()
    urls = cadc.get_image_list(cadc.query_region(c, radius=rad, collection=collection), c, rad)
    if name_filter:
        urls = [u_ for u_ in urls if re.search(name_filter, u_)]
    out = []
    for i, u_ in enumerate(urls):
        try:
            data = requests.get(u_, timeout=180).content
            if FITS_MAGIC in data[:80]:
                out.append((f"{collection.lower()}_{i:02d}.fits", data))
        except Exception as e:  # noqa: BLE001
            log(f"  skip {u_[:60]}: {str(e).splitlines()[0][:60]}")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--survey",
        default="vlass",
        help="vlass (CADC, working) | a raw CADC --collection | gleam/racs (routed, see notes)",
    )
    ap.add_argument(
        "--collection", help="explicit CADC collection name (overrides --survey routing)"
    )
    ap.add_argument("--ra", type=float, required=True, help="J2000 RA (deg)")
    ap.add_argument("--dec", type=float, required=True, help="J2000 Dec (deg)")
    ap.add_argument("--size-arcmin", type=float, default=1.5, help="cutout size (arcmin)")
    ap.add_argument("--out", default="cutouts", help="output directory")
    args = ap.parse_args(argv)

    out_dir = Path(args.out)
    survey = args.survey.lower()
    routes = {
        "vlass": ("VLASS", r"\.ql\."),  # Quick-Look images
    }
    if args.collection:
        collection, name_filter = args.collection, None
    elif survey in routes:
        collection, name_filter = routes[survey]
    elif survey in ("gleam", "gleam-x", "racs", "racs-i", "racs-low", "racs-mid", "racs-high"):
        log(f"survey '{survey}' is served by Data Central (Stokes I) or CASDA, not CADC.")
        log(
            "  GLEAM/GLEAM-X & RACS-I: Data Central POST cutout API "
            "(https://datacentral.org.au/api/services/cutout/; fits:true, band PK). See SKILL.md."
        )
        log("  RACS-V (Stokes V): CASDA only — use the casda-cutout-fetch skill.")
        return 6
    else:
        log(f"unknown survey '{survey}' and no --collection given")
        return 2

    try:
        log(
            f"querying CADC collection {collection} at ({args.ra}, {args.dec}) r={args.size_arcmin}'"
        )
        cutouts = fetch_cadc_cutout(
            collection, args.ra, args.dec, args.size_arcmin, name_filter=name_filter
        )
    except Exception as e:  # noqa: BLE001
        log(f"UNEXPECTED ERROR: {str(e).splitlines()[0]}")
        return 8
    if not cutouts:
        log("no image found at this position in that collection")
        return 5
    out_dir.mkdir(parents=True, exist_ok=True)
    saved = 0
    for name, data in cutouts:
        # validate it parses as a FITS image before keeping it
        try:
            from astropy.io import fits

            with fits.open(io.BytesIO(data)) as hd:
                _ = hd[0].header
        except Exception as e:  # noqa: BLE001
            log(f"  REJECT {name}: not a valid FITS ({str(e).splitlines()[0][:50]})")
            continue
        (out_dir / name).write_bytes(data)
        log(f"  saved {name} ({len(data)} bytes)")
        saved += 1
    if not saved:
        return 7
    log(f"DONE: {saved} FITS in {out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
