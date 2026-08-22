"""X-ray counterparts of long-period transients and white-dwarf-pulsar candidates (plan 93).

The question plan 93 set out to answer -- does accretion state predict radio loudness in
white-dwarf binaries? -- pairs Rose et al. 2026 (ASKAP J174508.9-505149 is an *accreting*
cataclysmic variable with orbitally modulated X-rays) against this repo's `wdpulsar` null
(zero persistent radio emission among 49 RACS-covered AR Sco-like candidates).

**GATE 0 cut the LPT leg of that comparison, and the reason is the slice's main methodological
result.** All three LPTs with a published X-ray detection (ASKAP J1832-0911,
ASKAP J144834-685644, ASKAP J174508.9-505149) are absent from every serendipitous X-ray source
catalogue, while each has *pointed* archival coverage: XMM 0953011101 ("VAST J1448-6856",
2024-08), XMM 0973390301 ("ASKAP J1745", 2025-10) and Chandra 26681/26682/29265/29266 (2024).
Those observations are dedicated follow-up taken *after* the catalogues were built, so an
archival catalogue cross-match is structurally blind to exactly the data that produced the
known detections: recall is 1/3 on the LPT sample.

On the white-dwarf candidates the same machinery recovers **16 of 16** of the X-ray
identifications Pelisoli et al. list independently (11 XMM-Newton, 5 ROSAT). Catalogue
cross-matching therefore has power on that sample and none on the LPTs, and this module
reports the two legs separately rather than differencing them.

Disciplines carried in from `dr20radio` and `frblens`:

- **The chance-coincidence rate is measured, not assumed** -- from the local source density in
  an annulus around each target *and* from rigid position-shift trials, both computed offline
  from one cached cone per position (`scripts/lptxray_fetch.py`).
- **A fraction is reported with both its difference and its ratio.** A difference in
  percentage points inherits the normalisation; a ratio does not.
- **A null divides by sensitivity, not by sample size.** `catalogue_recall` is the guard: a
  detection fraction from a catalogue that cannot see the known detections is not a
  measurement, and the module refuses to let one be quoted without its recall.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from .dr20radio import wilson_interval

__all__ = [
    "CATALOGS",
    "MATCH_RADIUS_ARCSEC",
    "FLUX_COLUMN",
    "ACCRETING_TYPES",
    "angular_sep_arcsec",
    "associate",
    "local_density_per_sqarcmin",
    "chance_expected",
    "shift_trial_rate",
    "catalogue_recall",
    "fraction_comparison",
    "summarize_position",
    "ON_AXIS_ARCMIN",
    "classify_coverage",
]

#: A pointing whose aimpoint lands this close to the target was almost certainly pointed *at*
#: it. Beyond it the target is somewhere in the field of view of an observation of something
#: else, which is a different kind of evidence and is counted separately.
ON_AXIS_ARCMIN = 2.0

CATALOGS = ("xmmssc", "rass2rxs", "erass1main", "csc")

#: Association radius per catalogue, set by that catalogue's positional accuracy.
#: 2RXS is the loose one (ROSAT positions carry a 10-30 arcsec error with a tail), which is
#: why the chance rate has to be measured rather than assumed.
MATCH_RADIUS_ARCSEC = {"xmmssc": 15.0, "rass2rxs": 45.0, "erass1main": 20.0, "csc": 10.0}

#: Flux-like column per catalogue. Units differ (2RXS is a count rate, the others are fluxes
#: in different bands), so these are NEVER combined into one number -- they are reported per
#: catalogue. Cross-band flux comparison would need a spectral assumption this sample cannot
#: support.
FLUX_COLUMN = {
    "xmmssc": ("ep_flux", "erg/cm2/s", "0.2-12 keV"),
    "rass2rxs": ("count_rate", "ct/s", "0.1-2.4 keV"),
    "erass1main": ("b1_flux", "erg/cm2/s", "0.2-2.3 keV"),
    "csc": ("b_flux_ap", "erg/cm2/s", "0.5-7 keV"),
}

#: Pelisoli et al. classifications that imply ongoing accretion. AR Sco itself -- the template
#: for the sample -- is NOT among these: it is a non-accreting ejector/propeller system.
ACCRETING_TYPES = ("polar", "IP", "CV")


def angular_sep_arcsec(
    ra1: float, dec1: float, ra2: np.ndarray | float, dec2: np.ndarray | float
) -> np.ndarray:
    """Great-circle separation in arcsec (haversine; exact for the small angles used here)."""
    r1, d1 = math.radians(ra1), math.radians(dec1)
    r2 = np.radians(np.asarray(ra2, dtype=float))
    d2 = np.radians(np.asarray(dec2, dtype=float))
    sin_d = np.sin((d2 - d1) / 2.0) ** 2
    sin_r = np.sin((r2 - r1) / 2.0) ** 2
    a = sin_d + math.cos(d1) * np.cos(d2) * sin_r
    return np.degrees(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))) * 3600.0


def _coords(rows: list[dict[str, Any]]) -> tuple[np.ndarray, np.ndarray]:
    if not rows:
        return np.array([]), np.array([])
    ra = np.array([float(r["ra"]) for r in rows])
    dec = np.array([float(r["dec"]) for r in rows])
    return ra, dec


def associate(
    ra: float, dec: float, rows: list[dict[str, Any]], radius_arcsec: float
) -> dict[str, Any]:
    """Nearest catalogue source to ``(ra, dec)`` and whether it falls inside ``radius_arcsec``."""
    cra, cdec = _coords(rows)
    if cra.size == 0:
        return {"matched": False, "n_within": 0, "sep_arcsec": None, "nearest": None}
    sep = angular_sep_arcsec(ra, dec, cra, cdec)
    order = int(np.argmin(sep))
    n_within = int(np.sum(sep <= radius_arcsec))
    return {
        "matched": bool(sep[order] <= radius_arcsec),
        "n_within": n_within,
        "sep_arcsec": round(float(sep[order]), 3),
        "nearest": rows[order].get("name"),
    }


def local_density_per_sqarcmin(
    ra: float,
    dec: float,
    rows: list[dict[str, Any]],
    *,
    cone_arcmin: float,
    inner_arcmin: float = 2.0,
) -> float:
    """Source surface density from an annulus around the target.

    The inner cut removes the target's own counterpart (and any physically associated
    source) so the density describes the *field*, which is what a chance rate needs.
    """
    cra, cdec = _coords(rows)
    if cra.size == 0:
        return 0.0
    sep_arcmin = angular_sep_arcsec(ra, dec, cra, cdec) / 60.0
    in_ann = np.sum((sep_arcmin > inner_arcmin) & (sep_arcmin <= cone_arcmin))
    area = math.pi * (cone_arcmin**2 - inner_arcmin**2)
    return float(in_ann / area) if area > 0 else 0.0


def chance_expected(density_per_sqarcmin: float, radius_arcsec: float) -> float:
    """Expected number of chance sources inside the association radius."""
    return float(density_per_sqarcmin * math.pi * (radius_arcsec / 60.0) ** 2)


def shift_trial_rate(
    ra: float,
    dec: float,
    rows: list[dict[str, Any]],
    *,
    radius_arcsec: float,
    cone_arcmin: float,
    shifts_arcmin: tuple[float, ...] = (3.0, 5.0, 7.0),
    n_azimuth: int = 8,
) -> dict[str, Any]:
    """Rigid position-shift trials: how often does a *displaced* position find a match?

    Every trial disc is required to lie wholly inside the cached cone, so a shifted position
    is searched against exactly the same catalogue completeness as the true one. This is the
    `dr20radio` chance-rate measurement, done offline from one cone rather than by requerying.
    """
    cra, cdec = _coords(rows)
    usable = [s for s in shifts_arcmin if s + radius_arcsec / 60.0 <= cone_arcmin]
    trials = 0
    hits = 0
    cosd = max(math.cos(math.radians(dec)), 1e-6)
    for s in usable:
        for k in range(n_azimuth):
            theta = 2.0 * math.pi * k / n_azimuth
            tra = ra + (s / 60.0) * math.cos(theta) / cosd
            tdec = dec + (s / 60.0) * math.sin(theta)
            trials += 1
            if cra.size and np.any(angular_sep_arcsec(tra, tdec, cra, cdec) <= radius_arcsec):
                hits += 1
    return {
        "n_trials": trials,
        "n_hits": hits,
        "rate": round(hits / trials, 5) if trials else float("nan"),
        "shifts_arcmin": list(usable),
    }


def catalogue_recall(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Recall of the cross-match against independently known X-ray detections.

    ``records`` carry ``lit_xray`` (the literature's own detection flag) and ``matched``.
    A detection fraction measured with a catalogue whose recall on *known* detections is low
    is a statement about the catalogue, not about the sources -- this is the number that
    decides whether a null means anything.
    """
    known = [r for r in records if r.get("lit_xray_detected")]
    recovered = [r for r in known if r.get("matched")]
    n_known = len(known)
    return {
        "n_known_detections": n_known,
        "n_recovered": len(recovered),
        "recall": round(len(recovered) / n_known, 4) if n_known else float("nan"),
        "missed": [r["name"] for r in known if not r.get("matched")],
        "usable": bool(n_known >= 3 and len(recovered) / n_known >= 0.8) if n_known else False,
    }


def fraction_comparison(k_a: int, n_a: int, k_b: int, n_b: int) -> dict[str, Any]:
    """Two detection fractions with their difference *and* their ratio.

    `dr20radio` learned this the hard way: a difference in percentage points inherits the
    normalisation of both arms, so it moves whenever a common limit moves; the ratio does
    not. Reporting only one of the two invites the wrong reading.
    """
    if n_a <= 0 or n_b <= 0:
        raise ValueError("both samples must be non-empty")
    f_a, lo_a, hi_a = wilson_interval(k_a, n_a, z=1.96)
    f_b, lo_b, hi_b = wilson_interval(k_b, n_b, z=1.96)
    return {
        "k_a": k_a,
        "n_a": n_a,
        "frac_a": round(f_a, 5),
        "ci_a": [round(lo_a, 5), round(hi_a, 5)],
        "k_b": k_b,
        "n_b": n_b,
        "frac_b": round(f_b, 5),
        "ci_b": [round(lo_b, 5), round(hi_b, 5)],
        "difference_pp": round(100.0 * (f_a - f_b), 4),
        "ratio": round(f_a / f_b, 5) if f_b > 0 else None,
    }


def classify_coverage(
    ra: float,
    dec: float,
    observations: list[dict[str, Any]],
    *,
    on_axis_arcmin: float = ON_AXIS_ARCMIN,
) -> dict[str, Any]:
    """Split archival pointings into ones aimed *at* the target and ones that merely cover it.

    Measured from the aimpoint offset rather than from the observation's target name, because
    names mislead in exactly the cases that matter: all 24 pointings near ASKAP J1935+2148 are
    observations of SGR 1935+2154, a different source a few arcmin away whose designation
    shares the same RA digits.

    A serendipitous pointing is real evidence -- the target is in the field -- but it is not
    the same as a dedicated observation, and conflating the two would overstate how much
    deliberate X-ray attention a source has had.
    """
    targeted: list[dict[str, Any]] = []
    serendipitous: list[dict[str, Any]] = []
    for obs in observations:
        if "ra" not in obs or "dec" not in obs:
            continue
        try:
            sep = float(angular_sep_arcsec(ra, dec, float(obs["ra"]), float(obs["dec"]))) / 60.0
        except (TypeError, ValueError):
            continue
        rec = {
            "obsid": obs.get("obsid"),
            "target": obs.get("name"),
            "offset_arcmin": round(sep, 3),
            "exposure": obs.get("exposure"),
        }
        (targeted if sep <= on_axis_arcmin else serendipitous).append(rec)
    return {
        "n_targeted": len(targeted),
        "n_serendipitous": len(serendipitous),
        "targeted": sorted(targeted, key=lambda r: r["offset_arcmin"]),
        "serendipitous": sorted(serendipitous, key=lambda r: r["offset_arcmin"]),
        "on_axis_arcmin": on_axis_arcmin,
    }


def summarize_position(
    meta: dict[str, Any],
    cats: dict[str, dict[str, Any]],
    *,
    cone_arcmin: float,
) -> dict[str, Any]:
    """Association, local density, expected chance count and shift trials for one position."""
    ra, dec = float(meta["ra"]), float(meta["dec"])
    out: dict[str, Any] = {
        "name": meta["name"],
        "sample": meta["sample"],
        "type": meta.get("type", ""),
        "ra": ra,
        "dec": dec,
        "lit_xray": meta.get("lit_xray", ""),
        "catalogs": {},
    }
    matched_any = False
    for cat in CATALOGS:
        entry = cats.get(cat, {})
        if "rows" not in entry:
            out["catalogs"][cat] = {"error": entry.get("error", "not queried")}
            continue
        rows = entry["rows"]
        radius = MATCH_RADIUS_ARCSEC[cat]
        assoc = associate(ra, dec, rows, radius)
        density = local_density_per_sqarcmin(ra, dec, rows, cone_arcmin=cone_arcmin)
        rec = {
            **assoc,
            "radius_arcsec": radius,
            "n_in_cone": len(rows),
            "density_per_sqarcmin": round(density, 6),
            "chance_expected": round(chance_expected(density, radius), 6),
            "shift_trials": shift_trial_rate(
                ra, dec, rows, radius_arcsec=radius, cone_arcmin=cone_arcmin
            ),
        }
        if assoc["matched"]:
            matched_any = True
            col = FLUX_COLUMN[cat][0]
            cra, cdec = _coords(rows)
            sep = angular_sep_arcsec(ra, dec, cra, cdec)
            best = rows[int(np.argmin(sep))]
            if col in best:
                rec["flux"] = best[col]
                rec["flux_units"] = FLUX_COLUMN[cat][1]
                rec["flux_band"] = FLUX_COLUMN[cat][2]
        out["catalogs"][cat] = rec
    out["matched"] = matched_any
    out["lit_xray_detected"] = str(meta.get("lit_xray", "")).strip().lower() not in ("", "no")
    return out
