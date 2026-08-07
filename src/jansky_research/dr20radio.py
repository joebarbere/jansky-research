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
    "parse_racs_csv",
    "run_north",
    "run_south",
    "select_quasars",
    "synthetic_survey",
    "synthetic_two_surveys",
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


RACS_TAP_SYNC = "https://casda.csiro.au/casda_vo_tools/tap/sync"
RACS_TABLE = "AS110.racs_dr1_sources_galacticcut_v2021_08_v02"  # 2,123,638 sources (Hale+ 2021)


def synthetic_two_surveys(
    *,
    fade_fraction: float = 0.35,
    seed: int = 0,
    **kwargs,
) -> dict:
    """Two-survey variant of :func:`synthetic_survey` — the increment-1 blind spot, fixed.

    The base survey acts as the SELECTING survey: every radio-carton quasar has a counterpart
    there by construction. The second ("matching") survey keeps only ``fade_fraction`` of the
    carton counterparts (spectral fading between observing frequencies) while ordinary
    counterparts carry over unchanged — so a census pipeline must see carton-match ~100%
    against the selecting survey but only ~``fade_fraction`` against the other one.
    """
    s = synthetic_survey(seed=seed, **kwargs)
    rng = np.random.default_rng(seed + 1)
    n_carton = int((np.char.startswith(s["firstcarton"], "openfibertargets")).sum())
    # the last n_carton radio sources are the carton counterparts (see synthetic_survey)
    keep = np.ones(s["ra_r"].size, dtype=bool)
    carton_slice = np.arange(s["ra_r"].size - n_carton, s["ra_r"].size)
    keep[carton_slice] = rng.random(n_carton) < fade_fraction
    s["ra_r2"] = s["ra_r"][keep]
    s["dec_r2"] = s["dec_r"][keep]
    s["fade_fraction"] = fade_fraction
    return s


def parse_racs_csv(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse a RACS TAP CSV chunk into ``(ra, dec, peak_flux_mjy)`` arrays."""
    ra, dec, flux = [], [], []
    for line in text.splitlines()[1:]:
        parts = line.split(",")
        if len(parts) != 3:
            continue
        try:
            ra.append(float(parts[0]))
            dec.append(float(parts[1]))
            flux.append(float(parts[2]))
        except ValueError:
            continue
    return np.array(ra), np.array(dec), np.array(flux)


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
        carton_s = np.char.strip(np.asarray(carton, dtype=str))
        return {
            "ra": np.asarray(d["RACAT"], float)[quasar],
            "dec": np.asarray(d["DECCAT"], float)[quasar],
            "z": np.asarray(d["Z"], float)[quasar],
            "obs": np.char.strip(np.asarray(d["OBS"], dtype=str))[quasar],
            "radio_carton": radio_carton[quasar],
            "carton_racs": np.char.startswith(carton_s, "openfibertargets_bhm_racsradio")[quasar],
            "carton_lofar": np.char.startswith(carton_s, "openfibertargets_bhm_lofarradio")[quasar],
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


def fetch_racs_positions(
    dest_dir: str = "data/racs_dr1",
    *,
    dec_min: float = -85.0,
    dec_max: float = 30.0,
    strip_deg: float = 1.0,
) -> dict:  # pragma: no cover - network
    """Fetch the RACS-low DR1 source positions by resumable 1-degree Dec strips.

    Each strip is cached as ``<dest>/strip_<dec>.csv`` (a completed strip is never re-fetched);
    the consolidated arrays are returned and cached as ``<dest>/racs_positions.npz``.
    """
    import time
    from pathlib import Path
    from urllib.parse import urlencode
    from urllib.request import urlopen

    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    npz = dest / "racs_positions.npz"
    if npz.exists():
        d = np.load(npz)
        return {"ra": d["ra"], "dec": d["dec"], "flux": d["flux"]}
    ras, decs, fluxes = [], [], []
    lo = dec_min
    while lo < dec_max:
        hi = min(lo + strip_deg, dec_max)
        cache = dest / f"strip_{lo:+.0f}.csv"
        if not cache.exists():
            q = f"SELECT ra, dec, peak_flux FROM {RACS_TABLE} WHERE dec >= {lo} AND dec < {hi}"
            params = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q})
            for attempt in range(6):
                try:
                    with urlopen(f"{RACS_TAP_SYNC}?{params}", timeout=180) as r:
                        text = r.read().decode()
                    if not text.startswith("ra"):
                        raise OSError("unexpected TAP response")
                    cache.write_text(text)
                    break
                except (TimeoutError, OSError):
                    time.sleep(min(2**attempt, 60))
            else:
                raise OSError(f"RACS strip {lo} failed after retries")
            print(f"[dr20radio] RACS strip {lo:+.0f}: cached", flush=True)
        ra, dec, fx = parse_racs_csv(cache.read_text())
        ras.append(ra)
        decs.append(dec)
        fluxes.append(fx)
        lo = hi
    out = {
        "ra": np.concatenate(ras),
        "dec": np.concatenate(decs),
        "flux": np.concatenate(fluxes),
    }
    np.savez_compressed(npz, ra=out["ra"], dec=out["dec"], flux=out["flux"])
    return out


def run_south(
    out: str = ".", *, radius_arcsec: float = 5.0, n_shift_trials: int = 10
) -> dict:  # pragma: no cover - network + bulk data (pure pieces tested offline)
    """Real leg B: the categorical first — DR20 quasars south of -40 deg vs RACS-low DR1.

    Also computes the overlap band (-40..+30) for the VLASS cross-check, and validates the
    racsradio carton against its SELECTING survey (expected ~100%). RACS-low's 25" beam
    motivates the wider default match radius (5"), with the false-match rate measured as
    always. Writes ``results/dr20radio_south.json``.
    """
    import json
    from pathlib import Path

    spall = fetch_spall()
    q = read_spall_quasars(spall)
    racs = fetch_racs_positions()
    zbins = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0])

    def census_block(sel: np.ndarray, label: str) -> dict:
        cen = sel & ~q["radio_carton"]
        m, _ = crossmatch(
            q["ra"][cen], q["dec"][cen], racs["ra"], racs["dec"], radius_arcsec=radius_arcsec
        )
        fm = false_match_rate(
            q["ra"][cen],
            q["dec"][cen],
            racs["ra"],
            racs["dec"],
            radius_arcsec=radius_arcsec,
            n_trials=n_shift_trials,
        )
        p, lo_w, hi_w = wilson_interval(int(m.sum()), int(m.size))
        return {
            "label": label,
            "n_census": int(m.size),
            "n_matched": int(m.sum()),
            "raw_fraction": round(p, 5),
            "wilson_lo": round(lo_w, 5),
            "wilson_hi": round(hi_w, 5),
            "false_match": fm,
            "corrected_fraction": round(p - fm["rate"], 5),
            "fraction_vs_z": detection_fraction(q["z"][cen], m, bins=zbins),
            "obs_breakdown": {o: int((q["obs"][cen] == o).sum()) for o in np.unique(q["obs"][cen])},
        }

    deep_south = q["dec"] <= VLASS_DEC_LIMIT_DEG
    overlap = (q["dec"] > VLASS_DEC_LIMIT_DEG) & (q["dec"] <= 30.0)

    # Carton validation split by SELECTING survey: racsradio cartons vs RACS is the true
    # ~100% pipeline validation; lofarradio (144 MHz-selected) vs RACS is a cross-frequency
    # fraction, exactly like the VLASS case in increment 1.
    def carton_block(mask: np.ndarray) -> dict:
        mm, _ = crossmatch(
            q["ra"][mask], q["dec"][mask], racs["ra"], racs["dec"], radius_arcsec=radius_arcsec
        )
        return {
            "n": int(mm.size),
            "matched": int(mm.sum()),
            "fraction": round(float(np.mean(mm)), 4) if mm.size else None,
        }

    in_racs_sky = q["dec"] <= 30.0
    carton_racs_v = carton_block(q["carton_racs"] & in_racs_sky)
    carton_lofar_v = carton_block(q["carton_lofar"] & in_racs_sky)
    metrics = {
        "source": f"SDSS-V DR20 spAll-lite x RACS-low DR1 ({RACS_TABLE})",
        "n_racs_sources": int(racs["ra"].size),
        "radius_arcsec": radius_arcsec,
        "deep_south": census_block(deep_south, "dec <= -40 (SDSS x RACS: categorical first)"),
        "overlap_band": census_block(overlap, "-40 < dec <= +30 (VLASS cross-check band)"),
        "carton_validation": {
            "racsradio_vs_racs_selecting_survey": carton_racs_v,
            "lofarradio_vs_racs_cross_frequency": carton_lofar_v,
        },
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "dr20radio_south.json").write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="DR20 BHM x VLASS/RACS radio census (plan 88).")
    p.add_argument("--out", default=".")
    p.add_argument("--north", action="store_true", help="run the VLASS northern leg")
    p.add_argument("--south", action="store_true", help="run the RACS southern leg")
    args = p.parse_args(argv)
    if args.south:
        m = run_south(args.out)
        slim = {
            k: (
                {kk: vv for kk, vv in v.items() if kk != "fraction_vs_z"}
                if isinstance(v, dict) and "fraction_vs_z" in v
                else v
            )
            for k, v in m.items()
        }
        print(json.dumps(slim, indent=2))
        return 0
    if args.north:
        m = run_north(args.out)
        slim = {k: v for k, v in m.items() if k != "fraction_vs_z_any_epoch"}
        print(json.dumps(slim, indent=2))
        return 0
    p.error("choose a mode: --north or --south")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
