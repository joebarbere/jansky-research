"""Jovian decametric (DAM) occurrence census from Juno/Waves in the (CML, Io-phase) plane (plan 37).

The classic ground-based result (Bigg 1964; Carr, Desch & Alexander 1983; Marques et al. 2017's
26-yr Nancay catalogue) is that Jovian DAM occurrence is organised by the observer's System III
central meridian longitude (CML) and Io's orbital phase into the Io-A/B/C/D source regions. The
public Juno/Waves Estimated Flux Density Dataset (doi:10.25935/6jg4-mk86; daily CDFs, 110 channels
to 40.5 MHz, per-channel ``Background``/``Sigma``) lets an occurrence census be built **from
Juno's vantage** --- the recover-a-known is Io-controlled enhancement; the vantage is the new part.

Two ephemeris subtleties (GATE 0, live-verified): the sub-Juno CML must come from JPL Horizons
(``PDObsLon``; the naive IAU :math:`W_{III}` rotation formula is wrong for Juno by up to ~40 deg),
and Io's phase uses the Lieske (1987) mean longitude
:math:`l_1 = 106.07719 + 203.488955790\\,(JD-2451545)` deg minus the CML. Units are
:math:`V^2 m^{-2} Hz^{-1}` (spectral power density, not W-flux) --- fine for an *occurrence*
census, which needs only detection above the shipped background.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "io_mean_longitude",
    "io_phase",
    "read_waves_cdf",
    "detect_active",
    "occurrence_map",
    "io_region_contrast",
    "synthetic_orbit",
    "run",
]

J2000_JD = 2451545.0
#: Lieske (1987, A&A 176, 146) E2x5 mean longitude of Io (deg, deg/day) --- <0.2 deg vs Horizons.
IO_L0, IO_RATE = 106.07719, 203.488955790
DAM_BAND_MHZ = (3.0, 40.5)
DATA_DIR = Path("data/junodam")
CDF_URL = (
    "https://maser.obspm.fr/repository/juno/waves/data/l3a/data/cdf/"
    "{y}/{m:02d}/jno_wav_cdr_lesia_{y}{m:02d}{d:02d}_v01.cdf"
)


def io_mean_longitude(jd: np.ndarray) -> np.ndarray:
    """Io's mean orbital longitude (deg, System III-adjacent frame; Lieske 1987)."""
    return (IO_L0 + IO_RATE * (np.asarray(jd, float) - J2000_JD)) % 360.0


def io_phase(jd: np.ndarray, cml_deg: np.ndarray) -> np.ndarray:
    """Io orbital phase: Io's longitude ahead of the sub-observer meridian (deg, 0-360)."""
    return (io_mean_longitude(jd) - np.asarray(cml_deg, float)) % 360.0


def fetch_cml_horizons(
    jd_start: float, jd_stop: float, *, step: str = "15m"
) -> dict:  # pragma: no cover - network
    """Sub-Juno System III CML from JPL Horizons (``PDObsLon``; observer ``500@-61``)."""
    from astroquery.jplhorizons import Horizons

    obj = Horizons(
        id="599",
        location="500@-61",
        epochs={"start": f"JD{jd_start}", "stop": f"JD{jd_stop}", "step": step},
    )
    eph = obj.ephemerides(quantities="14,20")
    return {
        "jd": np.asarray(eph["datetime_jd"], float),
        "cml": np.asarray(eph["PDObsLon"], float),
        "delta_au": np.asarray(eph["delta"], float),  # Juno--Jupiter range
    }


def read_waves_cdf(path: str | Path, *, bin_s: int = 15) -> dict:  # pragma: no cover - data file
    """Read one daily Juno/Waves L3a CDF -> DAM-band activity per ``bin_s`` time bin.

    Returns per-bin JD and the fraction of DAM-band channels whose median power in the bin
    exceeds ``Background + 5 Sigma`` (the shipped per-channel statistics).
    """
    import cdflib

    cdf = cdflib.CDF(str(path))
    epoch = cdflib.cdfepoch.to_datetime(cdf.varget("Epoch"))
    jd = np.array([e.astype("datetime64[s]").astype(float) for e in epoch]) / 86400.0 + 2440587.5
    freq_mhz = np.asarray(cdf.varget("Frequency"), float) / 1e3  # kHz -> MHz
    band = (freq_mhz >= DAM_BAND_MHZ[0]) & (freq_mhz <= DAM_BAND_MHZ[1])
    data = np.asarray(cdf.varget("Data"), float)[:, band]
    bg = np.asarray(cdf.varget("Background"), float)[band]
    sig = np.asarray(cdf.varget("Sigma"), float)[band]

    n = (data.shape[0] // bin_s) * bin_s
    d = data[:n].reshape(-1, bin_s, band.sum())
    med = np.median(d, axis=1)
    active_frac = (med > bg + 5.0 * sig).mean(axis=1)
    jd_bin = jd[:n].reshape(-1, bin_s).mean(axis=1)
    return {"jd": jd_bin, "active_frac": active_frac}


def detect_active(active_frac: np.ndarray, *, min_frac: float = 0.1) -> np.ndarray:
    """A time bin is 'DAM active' when >= ``min_frac`` of DAM-band channels exceed background."""
    return np.asarray(active_frac, float) >= min_frac


def occurrence_map(
    cml_deg: np.ndarray,
    io_phase_deg: np.ndarray,
    active: np.ndarray,
    *,
    n_bins: int = 18,
    min_exposure: int = 3,
) -> dict:
    """Occurrence probability in the (CML, Io-phase) plane, with per-cell exposure.

    Cells visited fewer than ``min_exposure`` times are masked NaN (a one-month orbit does not
    cover the plane uniformly --- exposure must be reported, not hidden).
    """
    edges = np.linspace(0.0, 360.0, n_bins + 1)
    exp_map, _, _ = np.histogram2d(cml_deg, io_phase_deg, bins=[edges, edges])
    act_map, _, _ = np.histogram2d(cml_deg[active], io_phase_deg[active], bins=[edges, edges])
    with np.errstate(invalid="ignore", divide="ignore"):
        occ = act_map / exp_map
    occ[exp_map < min_exposure] = np.nan
    return {"occ": occ, "exposure": exp_map, "edges": edges}


#: Canonical Io-controlled regions in (CML, Io-phase), after Carr et al. (1983) Fig. 7.32 /
#: Marques et al. (2017) Table 2 (deg ranges; coarse boxes for a contrast statistic, not fits).
IO_REGIONS = {
    "Io-A": ((200.0, 290.0), (195.0, 260.0)),
    "Io-B": ((90.0, 200.0), (65.0, 115.0)),
    "Io-C": ((280.0, 360.0), (200.0, 260.0)),
    "Io-D": ((0.0, 200.0), (95.0, 130.0)),
}


def io_region_contrast(m: dict) -> dict:
    """Mean occurrence inside the canonical Io boxes vs outside --- the recover-a-known statistic."""
    edges = m["edges"]
    cen = 0.5 * (edges[:-1] + edges[1:])
    cml_g, pha_g = np.meshgrid(cen, cen, indexing="ij")
    inside = np.zeros_like(cml_g, bool)
    for (c0, c1), (p0, p1) in IO_REGIONS.values():
        inside |= (cml_g >= c0) & (cml_g < c1) & (pha_g >= p0) & (pha_g < p1)
    occ = m["occ"]
    good = np.isfinite(occ)
    inside_occ = float(np.nanmean(occ[inside & good])) if (inside & good).any() else float("nan")
    outside_occ = float(np.nanmean(occ[~inside & good])) if (~inside & good).any() else float("nan")
    return {
        "occ_io_regions": inside_occ,
        "occ_elsewhere": outside_occ,
        "contrast": inside_occ / outside_occ if outside_occ > 0 else float("nan"),
        "cells_used": int(good.sum()),
    }


def synthetic_orbit(
    n_days: float = 28.0,
    *,
    bin_s: int = 15,
    p_in: float = 0.35,
    p_out: float = 0.04,
    seed: int = 0,
) -> dict:
    """Synthetic month: real CML/Io-phase rates, DAM active preferentially in the Io boxes.

    CML advances at Jupiter's System III rate (870.536 deg/day; Archinal et al. 2018) and Io's
    longitude at the Lieske rate, so the (CML, phase) plane fills exactly as it does for a real
    observer. Activity is Bernoulli: ``p_in`` inside the canonical boxes, ``p_out`` outside ---
    the census must recover a contrast ~``p_in/p_out``.
    """
    rng = np.random.default_rng(seed)
    jd = J2000_JD + np.arange(0.0, n_days, bin_s / 86400.0)
    cml = (284.95 + 870.5360000 * (jd - J2000_JD)) % 360.0
    pha = io_phase(jd, cml)
    inside = np.zeros(jd.size, bool)
    for (c0, c1), (p0, p1) in IO_REGIONS.values():
        inside |= (cml >= c0) & (cml < c1) & (pha >= p0) & (pha < p1)
    active = rng.random(jd.size) < np.where(inside, p_in, p_out)
    return {"jd": jd, "cml": cml, "io_phase": pha, "active": active, "p_in": p_in, "p_out": p_out}


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: synthetic-orbit recover-a-known. Real: one month of Juno/Waves CDFs + Horizons."""
    import json

    if offline:
        s = synthetic_orbit()
        cml, pha, active = s["cml"], s["io_phase"], s["active"]
        source = "synthetic orbit (canonical Io boxes injected)"
        expected = s["p_in"] / s["p_out"]
    else:  # pragma: no cover - data files + network
        files = sorted(DATA_DIR.glob("jno_wav_cdr_lesia_*_v01.cdf"))
        parts = [read_waves_cdf(f) for f in files]
        jd = np.concatenate([p["jd"] for p in parts])
        af = np.concatenate([p["active_frac"] for p in parts])
        eph = fetch_cml_horizons(float(jd.min()) - 0.02, float(jd.max()) + 0.02)
        # unwrap-interpolate the 15-min CML onto the data bins
        cml_unwrap = np.unwrap(np.radians(eph["cml"]))
        cml = np.degrees(np.interp(jd, eph["jd"], cml_unwrap)) % 360.0
        pha = io_phase(jd, cml)
        active = detect_active(af)
        source = f"Juno/Waves L3a v01, {len(files)} days"
        expected = float("nan")
        # the vantage dimension: Juno--Jupiter range dominates detection (proximity, not clock)
        dist = np.interp(jd, eph["jd"], eph["delta_au"])
        far = dist > np.median(dist)
        m_far = occurrence_map(cml[far], pha[far], active[far])
        con_far = io_region_contrast(m_far)
        extra = {
            "activity_near_half_pct": round(100 * float(active[~far].mean()), 2),
            "activity_far_half_pct": round(100 * float(active[far].mean()), 2),
            "io_contrast_far_half": round(con_far["contrast"], 2)
            if np.isfinite(con_far["contrast"])
            else None,
        }

    m = occurrence_map(cml, pha, active)
    con = io_region_contrast(m)
    metrics = {
        "source": source,
        "n_bins_total": int(active.size),
        "n_active": int(active.sum()),
        "occ_io_regions": round(con["occ_io_regions"], 3),
        "occ_elsewhere": round(con["occ_elsewhere"], 3),
        "io_contrast": round(con["contrast"], 2) if np.isfinite(con["contrast"]) else None,
        "expected_contrast": round(expected, 2) if np.isfinite(expected) else None,
        "cells_used": con["cells_used"],
    }
    if not offline:  # pragma: no cover - real-leg extras
        metrics.update(extra)
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "junodam_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(m, cml, pha, active, op / "papers" / "junodam" / "figures")
    _write_macros(metrics, op / "papers" / "junodam" / "generated" / "macros.tex")
    return metrics


def _figure(m: dict, cml, pha, active, out_dir) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0))
    im = ax1.pcolormesh(m["edges"], m["edges"], m["occ"].T, cmap="viridis", shading="auto")
    fig.colorbar(im, ax=ax1, label="occurrence probability")
    for name, ((c0, c1), (p0, p1)) in IO_REGIONS.items():
        ax1.add_patch(
            plt.Rectangle((c0, p0), c1 - c0, p1 - p0, fill=False, ec="w", lw=1.0, ls="--")
        )
        ax1.text(c0 + 3, p1 - 12, name, color="w", fontsize=7)
    ax1.set(xlabel="CML (System III, deg)", ylabel="Io phase (deg)", title="DAM occurrence map")
    ax2.plot(cml[~active][::7], pha[~active][::7], ".", ms=1, color="0.8")
    ax2.plot(cml[active], pha[active], ".", ms=2, color="C3")
    ax2.set(
        xlabel="CML (deg)",
        ylabel="Io phase (deg)",
        title="Active bins (red) vs coverage",
        xlim=(0, 360),
        ylim=(0, 360),
    )
    fig.tight_layout()
    fig.savefig(out / "junodam.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.junodam._write_macros -- do not edit.",
        rf"\newcommand{{\jdSource}}{{{m['source']}}}",
        rf"\newcommand{{\jdNbins}}{{{_fmt('n_bins_total')}}}",
        rf"\newcommand{{\jdNact}}{{{_fmt('n_active')}}}",
        rf"\newcommand{{\jdOccIo}}{{{_fmt('occ_io_regions')}}}",
        rf"\newcommand{{\jdOccOut}}{{{_fmt('occ_elsewhere')}}}",
        rf"\newcommand{{\jdContrast}}{{{_fmt('io_contrast')}}}",
        rf"\newcommand{{\jdExpContrast}}{{{_fmt('expected_contrast')}}}",
        rf"\newcommand{{\jdCells}}{{{_fmt('cells_used')}}}",
        rf"\newcommand{{\jdActNear}}{{{_fmt('activity_near_half_pct')}}}",
        rf"\newcommand{{\jdActFar}}{{{_fmt('activity_far_half_pct')}}}",
        rf"\newcommand{{\jdContrastFar}}{{{_fmt('io_contrast_far_half')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Juno/Waves DAM occurrence census (CML x Io phase).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
