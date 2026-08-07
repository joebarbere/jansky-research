"""First radio-counterpart census of the SDSS-V DR20 Black Hole Mapper quasars.

DR20 (Almeida et al., arXiv:2607.26149) delivers ~500k BHM spectroscopic objects including
the first optical SDSS spectra from the southern hemisphere. No radio cross-match of any
SDSS-V catalog exists (plans/88): this module matches the DR20 quasar table against VLASS
(north, local CIRADA epoch catalogs) and RACS (south), with the false-match rate *measured*
by position-shift trials and the radio-targeted open-fiber cartons excluded from fractions
(they were selected BECAUSE they are radio sources — counting them would be circular). Their
match rate against the SELECTING survey (RACS, southern leg) is the ~100% pipeline validation;
against VLASS at 3 GHz it is a cross-frequency detection fraction (steep-spectrum sources
selected at 144/888 MHz routinely fade below VLASS depth) — measured, not assumed.

Committed-real-results pattern: real legs write force-tracked ``results/dr20radio_*.json``;
paper macros come only from committed evidence; synthetic fixtures feed tests alone.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "RADIO_CARTON_PREFIXES",
    "SPALL_URL",
    "VLASS_DEC_LIMIT_DEG",
    "crossmatch",
    "detection_fraction",
    "false_match_rate",
    "run_north",
    "select_quasars",
    "synthetic_survey",
    "wilson_interval",
]

SPALL_URL = (
    "https://data.sdss.org/sas/dr20/spectro/boss/redux/v6_2_1/summary/allepoch/"
    "spAll-lite-v6_2_1-allepoch.fits.gz"
)
# Open-fiber cartons that TARGETED radio sources — excluded from detection fractions
# (circular) and used as a positive-control validation set instead.
RADIO_CARTON_PREFIXES = (
    "openfibertargets_bhm_racsradio",
    "openfibertargets_bhm_lofarradio",
)
VLASS_DEC_LIMIT_DEG = -40.0


def select_quasars(
    cls: np.ndarray,
    zwarning: np.ndarray,
    firstcarton: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(quasar_mask, radio_carton_mask)`` from spAll columns.

    A quasar is ``CLASS == 'QSO'`` with a clean redshift (``ZWARNING == 0``). The radio-carton
    mask flags objects whose ``FIRSTCARTON`` starts with any radio-selected open-fiber prefix;
    it is reported separately, never inside the census fractions.
    """
    cls = np.char.strip(np.asarray(cls, dtype=str))
    carton = np.char.strip(np.asarray(firstcarton, dtype=str))
    quasar = (cls == "QSO") & (np.asarray(zwarning) == 0)
    radio_carton = np.zeros(carton.size, dtype=bool)
    for prefix in RADIO_CARTON_PREFIXES:
        radio_carton |= np.char.startswith(carton, prefix)
    return quasar, radio_carton


def crossmatch(
    ra1_deg: np.ndarray,
    dec1_deg: np.ndarray,
    ra2_deg: np.ndarray,
    dec2_deg: np.ndarray,
    *,
    radius_arcsec: float = 2.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-neighbour sky match of catalog 1 onto catalog 2.

    Returns ``(matched_mask, sep_arcsec)`` for catalog-1 rows (sep is to the nearest
    catalog-2 source whether or not it is within the radius).
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    c1 = SkyCoord(np.asarray(ra1_deg, float) * u.deg, np.asarray(dec1_deg, float) * u.deg)
    c2 = SkyCoord(np.asarray(ra2_deg, float) * u.deg, np.asarray(dec2_deg, float) * u.deg)
    _, sep, _ = c1.match_to_catalog_sky(c2)
    sep_as = sep.arcsec
    return sep_as <= radius_arcsec, sep_as


def false_match_rate(
    ra_q: np.ndarray,
    dec_q: np.ndarray,
    ra_r: np.ndarray,
    dec_r: np.ndarray,
    *,
    radius_arcsec: float = 2.5,
    n_trials: int = 10,
    shift_arcmin: tuple[float, float] = (5.0, 30.0),
    seed: int = 0,
) -> dict:
    """Chance-coincidence match rate, MEASURED by rigid position-shift trials.

    Each trial shifts every quasar position by a random offset (uniform in ``shift_arcmin``,
    random direction — large enough to decorrelate, small enough to stay in the same source
    density) and re-runs the match. The mean shifted match fraction is the false-match rate
    to subtract from raw detection fractions; the std across trials is its uncertainty.
    """
    rng = np.random.default_rng(seed)
    rates = []
    dec_q = np.asarray(dec_q, float)
    ra_q = np.asarray(ra_q, float)
    for _ in range(n_trials):
        amp = rng.uniform(*shift_arcmin, size=ra_q.size) / 60.0
        ang = rng.uniform(0, 2 * np.pi, size=ra_q.size)
        dec_s = np.clip(dec_q + amp * np.sin(ang), -89.9, 89.9)
        ra_s = (ra_q + amp * np.cos(ang) / np.cos(np.deg2rad(dec_s))) % 360.0
        m, _ = crossmatch(ra_s, dec_s, ra_r, dec_r, radius_arcsec=radius_arcsec)
        rates.append(float(np.mean(m)))
    return {
        "rate": float(np.mean(rates)),
        "std": float(np.std(rates)),
        "n_trials": n_trials,
        "radius_arcsec": radius_arcsec,
    }


def wilson_interval(k: int, n: int, *, z: float = 1.0) -> tuple[float, float, float]:
    """Wilson score interval for a binomial fraction: ``(fraction, lo, hi)``."""
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return p, max(center - half, 0.0), min(center + half, 1.0)


def detection_fraction(values: np.ndarray, matched: np.ndarray, *, bins: np.ndarray) -> dict:
    """Matched fraction binned over ``values`` (e.g. redshift), with Wilson errors."""
    values = np.asarray(values, float)
    matched = np.asarray(matched, bool)
    out: dict[str, list] = {
        "bin_lo": [],
        "bin_hi": [],
        "n": [],
        "k": [],
        "frac": [],
        "lo": [],
        "hi": [],
    }
    for lo, hi in zip(bins[:-1], bins[1:], strict=False):
        sel = (values >= lo) & (values < hi)
        n, k = int(sel.sum()), int((matched & sel).sum())
        p, plo, phi = wilson_interval(k, n)
        out["bin_lo"].append(float(lo))
        out["bin_hi"].append(float(hi))
        out["n"].append(n)
        out["k"].append(k)
        out["frac"].append(None if n == 0 else round(p, 5))
        out["lo"].append(None if n == 0 else round(plo, 5))
        out["hi"].append(None if n == 0 else round(phi, 5))
    return out


def synthetic_survey(
    *,
    n_quasars: int = 2000,
    n_radio: int = 5000,
    counterpart_fraction: float = 0.12,
    n_carton: int = 50,
    area_deg: float = 10.0,
    radius_arcsec: float = 2.5,
    seed: int = 0,
) -> dict:
    """Synthetic quasar + radio catalogs with known truth, for the offline round-trip.

    A ``counterpart_fraction`` of ordinary quasars gets a radio source planted at its
    position (within 1"); ``n_carton`` extra quasars are radio-TARGETED (counterpart by
    construction) and carry a radio-carton ``FIRSTCARTON`` string — the circularity the
    census must exclude. Everything lives in an ``area_deg``-sized box at Dec ~ +10.
    """
    rng = np.random.default_rng(seed)
    ra_q = rng.uniform(180.0, 180.0 + area_deg, n_quasars)
    dec_q = rng.uniform(10.0, 10.0 + area_deg, n_quasars)
    z_q = rng.uniform(0.3, 4.0, n_quasars)
    ra_r = rng.uniform(180.0, 180.0 + area_deg, n_radio)
    dec_r = rng.uniform(10.0, 10.0 + area_deg, n_radio)
    is_cp = rng.random(n_quasars) < counterpart_fraction
    jitter = 1.0 / 3600.0
    ra_r = np.concatenate([ra_r, ra_q[is_cp] + rng.normal(0, jitter, is_cp.sum())])
    dec_r = np.concatenate([dec_r, dec_q[is_cp] + rng.normal(0, jitter, is_cp.sum())])
    # radio-carton quasars: targeted because radio-detected -> counterpart by construction
    ra_c = rng.uniform(180.0, 180.0 + area_deg, n_carton)
    dec_c = rng.uniform(10.0, 10.0 + area_deg, n_carton)
    ra_r = np.concatenate([ra_r, ra_c + rng.normal(0, jitter, n_carton)])
    dec_r = np.concatenate([dec_r, dec_c + rng.normal(0, jitter, n_carton)])
    cls = np.array(["QSO"] * (n_quasars + n_carton))
    zwarning = np.zeros(n_quasars + n_carton, dtype=int)
    carton = np.array(
        ["bhm_spiders_agn"] * n_quasars + ["openfibertargets_bhm_racsradio_boss"] * n_carton
    )
    return {
        "ra_q": np.concatenate([ra_q, ra_c]),
        "dec_q": np.concatenate([dec_q, dec_c]),
        "z_q": np.concatenate([z_q, rng.uniform(0.3, 4.0, n_carton)]),
        "cls": cls,
        "zwarning": zwarning,
        "firstcarton": carton,
        "ra_r": ra_r,
        "dec_r": dec_r,
        "true_fraction": counterpart_fraction,
        "radius_arcsec": radius_arcsec,
    }


# ------------------------------------------------------------------------- real data legs


def fetch_spall(dest_dir: str = "data") -> str:  # pragma: no cover - network
    """Download the DR20 spAll-lite summary file (177 MiB, resumable); returns the path."""
    from pathlib import Path
    from urllib.request import Request, urlopen

    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / SPALL_URL.rsplit("/", 1)[1]
    if path.exists():
        return str(path)
    part = path.with_suffix(path.suffix + ".part")
    offset = part.stat().st_size if part.exists() else 0
    req = Request(SPALL_URL)
    if offset:
        req.add_header("Range", f"bytes={offset}-")
    with urlopen(req, timeout=60) as r, open(part, "ab" if offset else "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    part.rename(path)
    return str(path)


def read_spall_quasars(path: str) -> dict:
    """Read the spAll-lite columns the census needs and apply :func:`select_quasars`."""
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        d = hdul[1].data
        cls = d["CLASS"]
        zw = np.asarray(d["ZWARNING"])
        carton = d["FIRSTCARTON"]
        quasar, radio_carton = select_quasars(cls, zw, carton)
        return {
            "ra": np.asarray(d["RACAT"], float)[quasar],
            "dec": np.asarray(d["DECCAT"], float)[quasar],
            "z": np.asarray(d["Z"], float)[quasar],
            "obs": np.char.strip(np.asarray(d["OBS"], dtype=str))[quasar],
            "radio_carton": radio_carton[quasar],
            "n_total_rows": int(len(d)),
        }


def load_vlass_positions() -> dict:  # pragma: no cover - local bulk files
    """Positions of quality-cut VLASS components from the local epoch catalogs.

    Applies the same cuts as the merged `vlass` slice: ``Duplicate_flag < 2``,
    ``Quality_flag in (0, 4)``, ``S_Code != 'E'``. Returns ``{"E2": (ra, dec), "E3": ...}``.
    """
    import csv
    import gzip

    from astropy.io import fits

    out = {}
    ra, dec = [], []
    with gzip.open("data/CIRADA_VLASS2QLv2_table1_components.csv.gz", "rt") as fh:
        for r in csv.DictReader(fh):
            try:
                if int(float(r["Duplicate_flag"])) >= 2:
                    continue
                if int(float(r["Quality_flag"])) not in (0, 4):
                    continue
                if r["S_Code"].strip() == "E":
                    continue
                ra.append(float(r["RA"]))
                dec.append(float(r["DEC"]))
            except (KeyError, ValueError):
                continue
    out["E2"] = (np.array(ra), np.array(dec))
    ra3, dec3 = [], []
    for name in ("data/QL3.1_components.fits", "data/QL3.2_components.fits"):
        with fits.open(name, memmap=True) as hdul:
            d = hdul[1].data
            # The E3 interim lists (VLASS Memo 22) carry a simplified schema: a binary
            # quality `Flag` (0 = good) and no Duplicate_flag/Quality_flag columns; empty
            # islands (S_Code 'E') are already absent.
            ok = np.asarray(d["Flag"]) == 0
            ra3.append(np.asarray(d["RA"], float)[ok])
            dec3.append(np.asarray(d["DEC"], float)[ok])
    out["E3"] = (np.concatenate(ra3), np.concatenate(dec3))
    return out


def run_north(
    out: str = ".", *, radius_arcsec: float = 2.5, n_shift_trials: int = 10
) -> dict:  # pragma: no cover - network + bulk local data (pure pieces tested offline)
    """Real leg A: DR20 quasars (Dec > -40) vs the local VLASS epoch catalogs.

    Writes the committed evidence file ``results/dr20radio_north.json``.
    """
    import json
    from pathlib import Path

    spall = fetch_spall()
    q = read_spall_quasars(spall)
    north = q["dec"] > VLASS_DEC_LIMIT_DEG
    census = north & ~q["radio_carton"]
    carton = north & q["radio_carton"]
    vlass = load_vlass_positions()
    zbins = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0])
    epochs = {}
    matched_any = np.zeros(int(census.sum()), dtype=bool)
    for name, (ra_r, dec_r) in vlass.items():
        m, _ = crossmatch(
            q["ra"][census], q["dec"][census], ra_r, dec_r, radius_arcsec=radius_arcsec
        )
        fm = false_match_rate(
            q["ra"][census],
            q["dec"][census],
            ra_r,
            dec_r,
            radius_arcsec=radius_arcsec,
            n_trials=n_shift_trials,
        )
        mc, _ = crossmatch(
            q["ra"][carton], q["dec"][carton], ra_r, dec_r, radius_arcsec=radius_arcsec
        )
        matched_any |= m
        p, lo, hi = wilson_interval(int(m.sum()), int(m.size))
        epochs[name] = {
            "n_census": int(m.size),
            "n_matched": int(m.sum()),
            "raw_fraction": round(p, 5),
            "wilson_lo": round(lo, 5),
            "wilson_hi": round(hi, 5),
            "false_match": fm,
            "corrected_fraction": round(p - fm["rate"], 5),
            "carton_validation": {
                "n": int(mc.size),
                "matched": int(mc.sum()),
                "fraction": round(float(np.mean(mc)), 4) if mc.size else None,
            },
            "n_radio_sources": int(ra_r.size),
        }
    p_any, lo_any, hi_any = wilson_interval(int(matched_any.sum()), int(matched_any.size))
    metrics = {
        "source": "SDSS-V DR20 spAll-lite v6_2_1 allepoch x VLASS QL E2+E3 (local CIRADA)",
        "n_spall_rows": q["n_total_rows"],
        "n_quasars_clean": int((q["z"] > -1).sum()),
        "n_north_census": int(census.sum()),
        "n_north_radio_carton_excluded": int(carton.sum()),
        "obs_breakdown_north": {
            o: int((q["obs"][north] == o).sum()) for o in np.unique(q["obs"][north])
        },
        "radius_arcsec": radius_arcsec,
        "epochs": epochs,
        "any_epoch": {
            "n_matched": int(matched_any.sum()),
            "raw_fraction": round(p_any, 5),
            "wilson_lo": round(lo_any, 5),
            "wilson_hi": round(hi_any, 5),
        },
        "fraction_vs_z_any_epoch": detection_fraction(q["z"][census], matched_any, bins=zbins),
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "dr20radio_north.json").write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="DR20 BHM x VLASS/RACS radio census (plan 88).")
    p.add_argument("--out", default=".")
    p.add_argument("--north", action="store_true", help="run the VLASS northern leg")
    args = p.parse_args(argv)
    if args.north:
        m = run_north(args.out)
        slim = {k: v for k, v in m.items() if k != "fraction_vs_z_any_epoch"}
        print(json.dumps(slim, indent=2))
        return 0
    p.error("choose a mode: --north (RACS southern leg lands in increment 2)")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
