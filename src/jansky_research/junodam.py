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
RJ_AU = 71492.0 / 1.495978707e8  # Jupiter equatorial radius in AU (for range display)
CDF_URL = (
    "https://maser.obspm.fr/repository/juno/waves/data/l3a/data/cdf/"
    "{y}/{m:02d}/jno_wav_cdr_lesia_{y}{m:02d}{d:02d}_v01.cdf"
)


def io_mean_longitude(jd: np.ndarray) -> np.ndarray:
    """Io's mean orbital longitude (deg, System III-adjacent frame; Lieske 1987)."""
    return (IO_L0 + IO_RATE * (np.asarray(jd, float) - J2000_JD)) % 360.0


def io_phase(jd: np.ndarray, cml_deg: np.ndarray) -> np.ndarray:
    r"""Conventional Io phase :math:`\Phi_{Io}` --- Io's departure from SUPERIOR conjunction.

    :math:`\Phi_{Io} = \mathrm{CML} + 180^\circ - \Lambda_{Io}` (Bigg 1964 convention as
    codified by Carr et al. 1983): zero when Io is diametrically opposite the observer as seen
    from Jupiter. GATE-2 caught an earlier sign/offset error here (:math:`\Lambda_{Io} -
    \mathrm{CML}` = :math:`180^\circ - \Phi_{Io}`), which displaced every canonical box.
    """
    return (np.asarray(cml_deg, float) + 180.0 - io_mean_longitude(jd)) % 360.0


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
    floor = bg + 5.0 * sig
    active_frac = (med > floor).mean(axis=1)
    # per-bin 90th-percentile channel SNR vs the 5-sigma floor: snr_p90 >= 1 APPROXIMATES
    # active_frac >= 0.1 (>=10% of DAM channels above floor) -- exact except for linear-interpolation
    # ties at the 10% boundary, where p90 is marginally stricter. Its value for the null is that it
    # scales LINEARLY under the 1/r^2 distance correction (percentile(c*r,90) = c*percentile(r,90)),
    # so the corrected detector is a self-consistent p90 rule applied identically at every range.
    snr_p90 = np.percentile(med / floor, 90, axis=1)
    jd_bin = jd[:n].reshape(-1, bin_s).mean(axis=1)
    return {"jd": jd_bin, "active_frac": active_frac, "snr_p90": snr_p90}


def detect_active(active_frac: np.ndarray, *, min_frac: float = 0.1) -> np.ndarray:
    """A time bin is 'DAM active' when >= ``min_frac`` of DAM-band channels exceed background."""
    return np.asarray(active_frac, float) >= min_frac


def sensitivity_corrected_active(
    snr_p90: np.ndarray, dist_au: np.ndarray, *, ref_au: float | None = None
) -> np.ndarray:
    r"""The 1/r^2 sensitivity null: distance-correct each bin's DAM SNR, then re-detect.

    DAM power falls as :math:`1/r^2` with Juno--Jupiter range, so near perijove the *same*
    intrinsic emission clears the background+5$\sigma$ floor more often --- a pure sensitivity
    effect, not intrinsic occurrence. Correcting each bin's 90th-percentile SNR to a reference
    range (:math:`S\to S\,(r/r_\mathrm{ref})^2`, default ``ref_au`` = median range) and
    re-thresholding at :math:`\ge 1` is the null: if the proximity duty-cycle trend is only
    sensitivity, the corrected occurrence is flat with range. Any residual near/far ratio after
    correction bounds the intrinsic$+$beaming part. Mirrors the `skr` slice's null model.
    """
    s = np.asarray(snr_p90, float)
    d = np.asarray(dist_au, float)
    ref = float(np.nanmedian(d)) if ref_au is None else ref_au
    return (s * (d / ref) ** 2) >= 1.0


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


#: Canonical Io-controlled regions in (CML, conventional Io phase), after Carr et al. (1983)
#: Fig. 7.32 / Marques et al. (2017) Table 2. CML ranges with c0 > c1 WRAP through 0 (Io-C).
IO_REGIONS = {
    "Io-A": ((200.0, 290.0), (205.0, 260.0)),
    "Io-B": ((90.0, 200.0), (80.0, 110.0)),
    "Io-C": ((300.0, 20.0), (225.0, 260.0)),
    "Io-D": ((0.0, 200.0), (95.0, 130.0)),
}


def _in_box(cml, pha, box):
    """Membership with CML wrap support (c0 > c1 wraps through 0 deg)."""
    (c0, c1), (p0, p1) = box
    in_c = ((cml >= c0) | (cml < c1)) if c0 > c1 else ((cml >= c0) & (cml < c1))
    return in_c & (pha >= p0) & (pha < p1)


def io_region_contrast(m: dict) -> dict:
    """Mean occurrence inside the canonical Io boxes vs outside --- the recover-a-known statistic."""
    edges = m["edges"]
    cen = 0.5 * (edges[:-1] + edges[1:])
    cml_g, pha_g = np.meshgrid(cen, cen, indexing="ij")
    inside = np.zeros_like(cml_g, bool)
    for box in IO_REGIONS.values():
        inside |= _in_box(cml_g, pha_g, box)
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
    for box in IO_REGIONS.values():
        inside |= _in_box(cml, pha, box)
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
        files = sorted(DATA_DIR.glob("jno_wav_cdr_lesia_*_v0?.cdf"))
        parts = [read_waves_cdf(f) for f in files]
        months = [f.name.split("_")[4][:6] for f in files]
        jd = np.concatenate([p["jd"] for p in parts])
        af = np.concatenate([p["active_frac"] for p in parts])
        snr_p90 = np.concatenate([p["snr_p90"] for p in parts])
        # fetch Horizons per contiguous data segment (a >2-day gap starts a new segment),
        # so multi-month runs stay inside per-query epoch limits
        order = np.argsort(jd)
        jd, af, snr_p90 = jd[order], af[order], snr_p90[order]
        starts = [0] + list(np.where(np.diff(jd) > 2.0)[0] + 1) + [jd.size]
        cml = np.empty_like(jd)
        dist = np.empty_like(jd)
        for a, b in zip(starts[:-1], starts[1:], strict=True):
            eph = fetch_cml_horizons(float(jd[a]) - 0.02, float(jd[b - 1]) + 0.02)
            cml_unwrap = np.unwrap(np.radians(eph["cml"]))
            cml[a:b] = np.degrees(np.interp(jd[a:b], eph["jd"], cml_unwrap)) % 360.0
            dist[a:b] = np.interp(jd[a:b], eph["jd"], eph["delta_au"])
        pha = io_phase(jd, cml)
        active = detect_active(af)
        active_p90 = snr_p90 >= 1.0  # the p90 surrogate detector (uncorrected baseline)
        active_corr = sensitivity_corrected_active(snr_p90, dist)  # 1/r^2 null, SAME p90 detector
        source = f"Juno/Waves L3a v01+v02, {len(files)} days"
        expected = float("nan")
        # the vantage dimension: Juno--Jupiter range dominates detection (proximity, not clock)
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
        # distance-RESOLVED Io contrast (the paper's scoped test): quartiles of Juno range
        qs = np.quantile(dist, [0.25, 0.5, 0.75])
        qmasks = [
            dist <= qs[0],
            (dist > qs[0]) & (dist <= qs[1]),
            (dist > qs[1]) & (dist <= qs[2]),
            dist > qs[2],
        ]
        for k, msk in enumerate(qmasks, start=1):
            cq = io_region_contrast(occurrence_map(cml[msk], pha[msk], active[msk]))
            extra[f"io_contrast_q{k}"] = (
                round(cq["contrast"], 2) if np.isfinite(cq["contrast"]) else None
            )
            extra[f"activity_q{k}_pct"] = round(100 * float(active[msk].mean()), 2)
            # the 1/r^2 sensitivity null: same range bins, distance-corrected detection
            extra[f"activity_q{k}_corr_pct"] = round(100 * float(active_corr[msk].mean()), 3)
        # raw vs sensitivity-corrected near/far ratio: how much of the ~180x proximity trend
        # survives once the 1/r^2 detection advantage of perijove is divided out. near_far_raw is
        # the original active_frac detector; near_far_raw_p90 is the SAME p90 detector as the
        # corrected column (so raw_p90 -> corrected isolates the distance correction, not a
        # detector swap); the two raw baselines agree closely.
        near_raw, far_raw = float(active[qmasks[0]].mean()), float(active[qmasks[3]].mean())
        near_p, far_p = float(active_p90[qmasks[0]].mean()), float(active_p90[qmasks[3]].mean())
        near_c, far_c = float(active_corr[qmasks[0]].mean()), float(active_corr[qmasks[3]].mean())
        extra["near_far_raw"] = round(near_raw / far_raw, 1) if far_raw > 0 else None
        extra["near_far_raw_p90"] = round(near_p / far_p, 1) if far_p > 0 else None
        extra["near_far_corrected"] = round(near_c / far_c, 1) if far_c > 0 else None
        extra["range_near_rj"] = round(float(np.median(dist[qmasks[0]])) / RJ_AU, 1)
        extra["range_far_rj"] = round(float(np.median(dist[qmasks[3]])) / RJ_AU, 1)
        # per-month Io contrast: the robustness spread across orbital configurations
        pm = []
        for ym in sorted(set(months)):
            sel = np.zeros(jd.size, bool)
            for f_m, p in zip(months, parts, strict=True):
                if f_m == ym:
                    sel[np.searchsorted(jd, p["jd"])] = True
            cm = io_region_contrast(occurrence_map(cml[sel], pha[sel], active[sel]))
            if np.isfinite(cm["contrast"]):
                pm.append(cm["contrast"])
        if pm:
            extra["n_months"] = len(pm)
            extra["io_contrast_month_median"] = round(float(np.median(pm)), 2)
            extra["io_contrast_month_min"] = round(float(np.min(pm)), 2)
            extra["io_contrast_month_max"] = round(float(np.max(pm)), 2)

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
        spans = [(c0, c1)] if c0 < c1 else [(c0, 360.0), (0.0, c1)]  # wrap-aware (Io-C)
        for a, b in spans:
            ax1.add_patch(
                plt.Rectangle((a, p0), b - a, p1 - p0, fill=False, ec="w", lw=1.0, ls="--")
            )
        ax1.text(spans[0][0] + 3, p1 - 12, name, color="w", fontsize=7)
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
        rf"\newcommand{{\jdCqA}}{{{_fmt('io_contrast_q1')}}}",
        rf"\newcommand{{\jdCqB}}{{{_fmt('io_contrast_q2')}}}",
        rf"\newcommand{{\jdCqC}}{{{_fmt('io_contrast_q3')}}}",
        rf"\newcommand{{\jdCqD}}{{{_fmt('io_contrast_q4')}}}",
        rf"\newcommand{{\jdAqA}}{{{_fmt('activity_q1_pct')}}}",
        rf"\newcommand{{\jdAqD}}{{{_fmt('activity_q4_pct')}}}",
        rf"\newcommand{{\jdAqAcorr}}{{{_fmt('activity_q1_corr_pct')}}}",
        rf"\newcommand{{\jdAqDcorr}}{{{_fmt('activity_q4_corr_pct')}}}",
        rf"\newcommand{{\jdNearFarRaw}}{{{_fmt('near_far_raw')}}}",
        rf"\newcommand{{\jdNearFarRawPninety}}{{{_fmt('near_far_raw_p90')}}}",
        rf"\newcommand{{\jdNearFarCorr}}{{{_fmt('near_far_corrected')}}}",
        rf"\newcommand{{\jdAqBcorr}}{{{_fmt('activity_q2_corr_pct')}}}",
        rf"\newcommand{{\jdAqCcorr}}{{{_fmt('activity_q3_corr_pct')}}}",
        rf"\newcommand{{\jdRangeNear}}{{{_fmt('range_near_rj')}}}",
        rf"\newcommand{{\jdRangeFar}}{{{_fmt('range_far_rj')}}}",
        rf"\newcommand{{\jdNmonths}}{{{_fmt('n_months')}}}",
        rf"\newcommand{{\jdCmMed}}{{{_fmt('io_contrast_month_median')}}}",
        rf"\newcommand{{\jdCmMin}}{{{_fmt('io_contrast_month_min')}}}",
        rf"\newcommand{{\jdCmMax}}{{{_fmt('io_contrast_month_max')}}}",
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
