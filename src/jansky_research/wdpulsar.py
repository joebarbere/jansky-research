"""Radio counterpart survey of the 56 white-dwarf-pulsar candidates (plan 41).

Pelisoli et al. 2025 (MNRAS 540, 821; arXiv:2505.04693) published 56 AR Sco-like binary
white-dwarf-pulsar candidates from Gaia+WISE light curves (26 newly characterised; type census
26 YSO / 16 polar / 3 CV / 2 IP / 8 unclear / 1 pulsar). GATE-0 (2026-07-07): **no systematic
radio search of the list exists** (the two citing papers are single-object work); there is NO
VizieR/CDS deposit --- the machine-readable source is Table 2 of the arXiv HTML, vendored here
as `data/wdpulsar_candidates.csv` (the `lpt` provenance pattern; positions to ~0.1 arcsec).

Two GATE-0 corrections to plan 41's anchors, both load-bearing:

- **J1912-4410 is NOT a plausible RACS re-detection.** Its MeerKAT pulses (<4 s FWHM, ~15
  mJy peaks, ~1% duty cycle; Pelisoli+2023, arXiv:2306.09272) dilute to ~0.1-0.2 mJy time-
  averaged, and it is absent from the RACS-low DR1 catalogue (verified cone search). Its
  forced limit is kept as the paper's teachable point: for duty-cycled pulsators a survey
  non-detection is NOT an absence.
- **AR Sco is the pipeline control instead** (it is the list's template, deliberately not
  among the 56): RACS-low DR1 detects it at 8.58 +- 0.92 mJy (RACS J162146.6-225311), and at
  1.5 GHz it shows |V|/I ~ 22-27% --- the forced I+V photometry must re-find it (the GJ 65
  pattern from `stokesv_discovery`).

Timing note: the RACS-low2 release (arXiv:2606.16182) ships a circular-polarization catalogue
--- simultaneously this slice's dataset and the most likely scoop vector. Reuse:
`stokesv.measure_circular_pol` / `classify_emitter` for the forced photometry and vetting,
`vlass.fetch_vlass_epoch` (local bulk catalogues) for the Dec > -40 VLASS cone checks, and the
resumable per-row-CSV driver pattern (`scripts/wdpulsar_real.py`).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .stokesv import classify_emitter, measure_circular_pol

__all__ = [
    "load_candidate_table",
    "candidate_targets",
    "synthetic_cutout",
    "injection_roundtrip",
    "summarize_sweep",
    "vlass_cone_checks",
    "run",
]

CANDIDATES_LOCAL = Path("data/wdpulsar_candidates.csv")
SWEEP_CSV = Path("results/wdpulsar_realtargets.csv")
# the list's template object (deliberately NOT one of the 56): the forced-photometry control.
# RACS-low DR1: RACS J162146.6-225311, 8.58 +- 0.92 mJy at 887.5 MHz.
AR_SCO = {
    "name": "AR_Sco",
    "ra_deg": 245.44700,
    "dec_deg": -22.88617,
    "type": "control",
    "racs_dr1_mjy": 8.58,
}
EXPECTED_TYPES = {"YSO": 26, "polar": 16, "CV": 3, "IP": 2, "unclear": 8, "pulsar": 1}


def load_candidate_table(path: str | Path = CANDIDATES_LOCAL) -> dict:
    """Load the vendored 56-row Pelisoli Table 2 mirror, with structural validation.

    Validates the row count, the per-type census against the paper's numbers, and that the
    one confirmed system (J1912-4410, type ``pulsar``) is present --- the same transcription
    discipline the `lpt` slice used (it caught a source-paper typo that way).
    """
    import csv

    with open(path) as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 56:
        raise ValueError(f"candidate table has {len(rows)} rows, expected 56")
    census: dict[str, int] = {}
    for r in rows:
        census[r["type"]] = census.get(r["type"], 0) + 1
    if census != EXPECTED_TYPES:
        raise ValueError(f"type census {census} != paper's {EXPECTED_TYPES}")
    if not any(r["short_name"] == "J1912-4410" and r["type"] == "pulsar" for r in rows):
        raise ValueError("J1912-4410 (the confirmed WD pulsar) missing from the table")
    return {
        "gaia_source_id": np.array([r["gaia_source_id"] for r in rows]),
        "simbad_id": np.array([r["simbad_id"] for r in rows]),
        "ra_deg": np.array([float(r["ra_deg"]) for r in rows]),
        "dec_deg": np.array([float(r["dec_deg"]) for r in rows]),
        "short_name": np.array([r["short_name"] for r in rows]),
        "type": np.array([r["type"] for r in rows]),
        # Table 2's asterisk = "classification determined as part of this work" (17 rows);
        # the abstract's "26 previously uncharacterised" = rows with NO Simbad ID. Both kept:
        # conflating them was a transcription bug the census validation caught.
        "class_this_work": np.array([r["class_this_work"] == "True" for r in rows]),
        "previously_uncharacterized": np.array([not r["simbad_id"].strip() for r in rows]),
    }


def candidate_targets(cat: dict, *, include_control: bool = True) -> list[dict]:
    """The sweep's target list: the 56 candidates (+ the AR Sco control, flagged)."""
    targets = [
        {
            "name": str(cat["short_name"][i]),
            "ra_deg": float(cat["ra_deg"][i]),
            "dec_deg": float(cat["dec_deg"][i]),
            "type": str(cat["type"][i]),
        }
        for i in range(cat["ra_deg"].size)
    ]
    if include_control:
        targets.append(dict(AR_SCO))
    return targets


def synthetic_cutout(
    *,
    ra0: float = 245.447,
    dec0: float = -22.886,
    i_mjy: float = 8.6,
    v_mjy: float = -2.1,
    rms_mjy: float = 0.25,
    beam_arcsec: float = 15.0,
    n_pix: int = 120,
    pix_arcsec: float = 2.5,
    seed: int = 0,
):
    """A synthetic RACS-like (I, V, wcs) cutout with a Gaussian point source at ``(ra0, dec0)``.

    The offline recover-a-known fixture: `stokesv.measure_circular_pol` must return the
    injected fluxes within the noise, and a blank field must yield an honest limit.
    """
    from astropy.wcs import WCS

    rng = np.random.default_rng(seed)
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [n_pix / 2 + 0.5, n_pix / 2 + 0.5]
    w.wcs.cdelt = [-pix_arcsec / 3600.0, pix_arcsec / 3600.0]
    yy, xx = np.mgrid[0:n_pix, 0:n_pix]
    sig_pix = beam_arcsec / 2.3548 / pix_arcsec
    psf = np.exp(
        -0.5 * ((xx - (n_pix / 2 - 0.5)) ** 2 + (yy - (n_pix / 2 - 0.5)) ** 2) / sig_pix**2
    )
    image_i = i_mjy * psf + rng.normal(0.0, rms_mjy, (n_pix, n_pix))
    image_v = v_mjy * psf + rng.normal(0.0, rms_mjy, (n_pix, n_pix))
    return image_i, image_v, w


def injection_roundtrip(*, seed: int = 0) -> dict:
    """The offline recover-a-known: injected (I, V) recovered; blank field stays a limit."""
    i_in, v_in = 8.6, -2.1
    img_i, img_v, w = synthetic_cutout(i_mjy=i_in, v_mjy=v_in, seed=seed)
    hit = measure_circular_pol(img_i, img_v, w, 245.447, -22.886)
    blank_i, blank_v, wb = synthetic_cutout(i_mjy=0.0, v_mjy=0.0, seed=seed + 1)
    miss = measure_circular_pol(blank_i, blank_v, wb, 245.447, -22.886)
    return {
        "i_in": i_in,
        "v_in": v_in,
        "i_out": hit["i_peak"],
        "v_out": hit["v_peak"],
        "frac_out": hit["frac_pol"],
        "class": classify_emitter(hit["v_peak"], hit["i_peak"]),
        "blank_i_sig": miss["i_peak"] / miss["i_rms"] if miss["i_rms"] > 0 else float("nan"),
    }


def summarize_sweep(
    csv_path: str | Path,
    targets: list[dict],
    *,
    det_sigma: float = 5.0,
    limit_sigma: float = 3.0,
) -> dict:
    """Reduce the resumable sweep CSV (one row per target x epoch) into the survey table.

    Per target: number of epochs measured, the maximum-significance I and V measurement
    across epochs, detection flags at ``det_sigma``, the DEEPEST ``limit_sigma`` upper limits
    (best epoch), and the polarization class where I is detected. A V "detection" must also
    exceed a per-measurement leakage veto (|V| > 0.6% of the coincident I --- the characteristic
    ASKAP on-axis instrumental-leakage level) to survive vetting. (This is a raw-level veto, not
    the `stokesv` slice's per-region 7-sigma-of-median-|V/I| floor; it is the weaker of the two,
    but with zero 5-sigma V candidates the choice does not affect the result.)
    """
    import csv

    with open(csv_path) as f:
        rows = [r for r in csv.DictReader(f) if r.get("i_mjy") not in ("", "nan", None)]
    by: dict[str, list[dict]] = {}
    for r in rows:
        by.setdefault(r["name"], []).append(r)
    out = []
    for t in targets:
        rs = by.get(t["name"], [])
        entry: dict = {
            "name": t["name"],
            "type": t["type"],
            "n_epochs": len(rs),
        }
        if rs:
            i_sig = [float(r["i_mjy"]) / float(r["e_i"]) for r in rs]
            v_sig = [abs(float(r["v_mjy"])) / float(r["e_v"]) for r in rs]
            k_i = int(np.argmax(i_sig))
            k_v = int(np.argmax(v_sig))
            i_peak = float(rs[k_i]["i_mjy"])
            v_peak = float(rs[k_v]["v_mjy"])
            i_det = i_sig[k_i] >= det_sigma
            # leakage vetting: |V| must beat both the noise and 0.6% of the coincident I
            v_det = bool(
                v_sig[k_v] >= det_sigma and abs(v_peak) > 0.006 * max(float(rs[k_v]["i_mjy"]), 0.0)
            )
            entry.update(
                {
                    "i_det": bool(i_det),
                    "v_det": v_det,
                    "i_mjy": round(i_peak, 3) if i_det else None,
                    "v_mjy": round(v_peak, 3) if v_det else None,
                    "i_limit_mjy": round(limit_sigma * min(float(r["e_i"]) for r in rs), 3),
                    "v_limit_mjy": round(limit_sigma * min(float(r["e_v"]) for r in rs), 3),
                    "class": classify_emitter(v_peak, i_peak) if i_det else "nan",
                }
            )
        out.append(entry)
    measured = [e for e in out if e["n_epochs"] > 0]
    i_dets = [e for e in measured if e.get("i_det")]
    v_dets = [e for e in measured if e.get("v_det")]
    control = next((e for e in out if e["name"] == AR_SCO["name"]), None)
    return {
        "per_target": out,
        "n_targets": len(targets),
        "n_measured": len(measured),
        "n_i_detections": len(i_dets),
        "n_v_detections": len(v_dets),
        "i_detection_names": sorted(e["name"] for e in i_dets),
        "v_detection_names": sorted(e["name"] for e in v_dets),
        "median_v_limit_mjy": round(float(np.median([e["v_limit_mjy"] for e in measured])), 3)
        if measured
        else None,
        "control_i_det": bool(control and control.get("i_det")),
        "control_i_mjy": control.get("i_mjy") if control else None,
    }


def vlass_cone_checks(
    targets: list[dict], *, radius_arcsec: float = 7.5, epoch: int = 3
) -> list[dict]:  # pragma: no cover - needs the local ~1 GB bulk catalogues
    """VLASS Quick-Look cone checks (3 GHz, Stokes I) for the Dec > -40 targets.

    Uses `vlass.fetch_vlass_epoch`'s cached bulk catalogues; returns matched components
    (peak flux, separation) per target. Complements RACS: higher frequency, northern reach.
    """
    from .vlass import fetch_vlass_epoch

    out = []
    for t in targets:
        if t["dec_deg"] <= -40.0:
            continue
        ra, dec, fp, efp = fetch_vlass_epoch(epoch, (t["ra_deg"], t["dec_deg"]), 0.01)
        sep = (
            np.hypot((ra - t["ra_deg"]) * np.cos(np.radians(t["dec_deg"])), dec - t["dec_deg"])
            * 3600.0
        )
        m = sep <= radius_arcsec
        out.append(
            {
                "name": t["name"],
                "n_match": int(m.sum()),
                "peak_mjy": round(float(fp[m].max()), 3) if m.any() else None,
                "sep_arcsec": round(float(sep[m].min()), 2) if m.any() else None,
            }
        )
    return out


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: table validation + injection round-trip; real: reduce the sweep CSV (+VLASS)."""
    import json

    cat = load_candidate_table()
    targets = candidate_targets(cat)
    rt = injection_roundtrip()
    metrics: dict = {
        "source": "vendored Pelisoli+2025 Table 2 (arXiv:2505.04693) + synthetic injection",
        "is_real": not offline,
        "n_candidates": int(cat["ra_deg"].size),
        "n_uncharacterized": int(cat["previously_uncharacterized"].sum()),
        "n_class_this_work": int(cat["class_this_work"].sum()),
        "injection": {k: (round(v, 3) if isinstance(v, float) else v) for k, v in rt.items()},
    }

    if not offline:  # pragma: no cover - needs the sweep CSV (CASDA) + bulk VLASS catalogues
        summary = summarize_sweep(SWEEP_CSV, targets)
        vlass = vlass_cone_checks(targets)
        metrics.update(summary)
        metrics["vlass"] = vlass
        metrics["n_vlass_matched"] = sum(1 for v in vlass if v["n_match"] > 0)
        # paper-facing derived numbers (the null is the deliverable, so name it precisely)
        by = {e["name"]: e for e in summary["per_target"]}
        cand = [e for e in summary["per_target"] if e["type"] != "control"]
        covered = [e for e in cand if e["n_epochs"] > 0]
        metrics["n_candidates_covered"] = len(covered)
        metrics["n_candidate_i_det"] = sum(1 for e in cand if e.get("i_det"))
        metrics["n_candidate_v_det"] = sum(1 for e in cand if e.get("v_det"))
        vlims = [e["v_limit_mjy"] for e in covered if e.get("v_limit_mjy") is not None]
        metrics["v_limit_min_mjy"] = round(min(vlims), 3) if vlims else None
        metrics["v_limit_max_mjy"] = round(max(vlims), 3) if vlims else None
        ctrl = by.get(AR_SCO["name"], {})
        # SNR = I / rms; the 3sigma limit encodes rms = i_limit / 3, so SNR = 3 * I / i_limit
        metrics["control_snr"] = (
            round(3.0 * ctrl["i_mjy"] / ctrl["i_limit_mjy"], 1)
            if ctrl.get("i_mjy") and ctrl.get("i_limit_mjy")
            else None
        )
        j = by.get("J1912-4410", {})
        metrics["j1912_i_limit_mjy"] = j.get("i_limit_mjy")
        metrics["j1912_v_limit_mjy"] = j.get("v_limit_mjy")
        metrics["j1912_n_epochs"] = j.get("n_epochs")
        vc = [v for v in vlass if v["n_match"] > 0 and v["name"] != AR_SCO["name"]]
        metrics["vlass_candidate_detections"] = vc
        metrics["source"] = (
            "RACS forced I+V photometry (CASDA) + VLASS QL cones, 56 Pelisoli candidates"
        )

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    per_target = metrics.pop("per_target", None)
    (op / "results" / "wdpulsar_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    if per_target is not None:
        metrics["per_target"] = per_target
        (op / "results" / "wdpulsar_table.json").write_text(json.dumps(per_target, indent=2) + "\n")
    _figure(metrics, op / "papers" / "wdpulsar" / "figures")
    _write_macros(metrics, op / "papers" / "wdpulsar" / "generated" / "macros.tex")
    _write_limits_table(metrics, op / "papers" / "wdpulsar" / "generated" / "limits_table.tex")
    return metrics


def _figure(m: dict, out_dir: str | Path) -> None:
    try:
        from .report import _agg
    except ImportError:  # pragma: no cover - minimal venvs
        return
    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    img_i, _img_v, _w = synthetic_cutout()
    ax1.imshow(img_i, origin="lower", cmap="viridis")
    ax1.set(title="injection fixture (Stokes I)", xlabel="pix", ylabel="pix")
    per = m.get("per_target") or []
    lims = [e["v_limit_mjy"] for e in per if e.get("v_limit_mjy") is not None]
    if lims:
        ax2.hist(lims, bins=15, color="C0")
        ax2.set(xlabel="3$\\sigma$ V limit (mJy)", ylabel="targets", title="V limit census")
    else:
        ax2.text(0.5, 0.5, "real sweep pending", ha="center", va="center")
        ax2.set_axis_off()
    fig.tight_layout()
    fig.savefig(out / "wdpulsar.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    def g(key: str) -> str:
        v = m.get(key)
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    pref = "wdReal" if m.get("is_real") else "wdSyn"
    lines = [
        "% Auto-generated by jansky_research.wdpulsar._write_macros -- do not edit.",
        "% Synthetic (wdSyn*) and real (wdReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders, so synthetic numbers can never masquerade",
        "% under wdReal* (an offline rebuild resets wdReal* to placeholders by design).",
        rf"\newcommand{{\wdSource}}{{{m['source']}}}",
        rf"\newcommand{{\wdNCand}}{{{m['n_candidates']}}}",
        rf"\newcommand{{\wdNUnchar}}{{{m['n_uncharacterized']}}}",
        rf"\newcommand{{\wdNClassNew}}{{{m['n_class_this_work']}}}",
    ]
    keys = (
        ("NMeasured", "n_measured"),
        ("NCandCovered", "n_candidates_covered"),
        ("NCandIDet", "n_candidate_i_det"),
        ("NCandVDet", "n_candidate_v_det"),
        ("MedianVLimit", "median_v_limit_mjy"),
        ("VLimMin", "v_limit_min_mjy"),
        ("VLimMax", "v_limit_max_mjy"),
        ("ControlI", "control_i_mjy"),
        ("ControlSnr", "control_snr"),
        ("JLimI", "j1912_i_limit_mjy"),
        ("JLimV", "j1912_v_limit_mjy"),
        ("JEpochs", "j1912_n_epochs"),
        ("NVlass", "n_vlass_matched"),
    )
    for ns in ("wdSyn", "wdReal"):
        live = ns == pref
        for macro, key in keys:
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(key) if live else '--'}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _write_limits_table(m: dict, path: str | Path, *, n_excerpt: int = 16) -> None:
    """Excerpt of the per-target census as 6-column LaTeX; JSON carries the full machine copy.

    The excerpt leads with the load-bearing rows (the AR Sco control and the confirmed pulsar
    J1912-4410), then representative covered candidates; every row is name & type & epochs &
    I (mJy or limit) & 3sigma V limit & class, matching the paper's tabular.
    """
    per = m.get("per_target") or []
    by = {e["name"]: e for e in per}
    lead = [by[n] for n in ("AR_Sco", "J1912-4410") if n in by]
    rest = [e for e in per if e["name"] not in ("AR_Sco", "J1912-4410") and e["n_epochs"] > 0]
    rows_src = lead + rest[: max(0, n_excerpt - len(lead))]
    out = [
        "% Auto-generated by jansky_research.wdpulsar._write_limits_table -- do not edit.",
        "% Columns: name & type & epochs & I (mJy, or 3sig limit) & 3sig V limit & class",
    ]
    for e in rows_src:
        name = e["name"].replace("_", r"\_")  # LaTeX-safe (AR_Sco)
        if e["n_epochs"] == 0:
            out.append(f"{name} & {e['type']} & 0 & \\multicolumn{{3}}{{c}}{{uncovered}} \\\\")
            continue
        # I column: the detection if there is one, else the 3sigma upper limit prefixed "<"
        i_s = f"{e['i_mjy']:.2f}" if e.get("i_mjy") is not None else f"$<${e['i_limit_mjy']:.2f}"
        out.append(
            f"{name} & {e['type']} & {e['n_epochs']} & {i_s} & "
            f"$<${e['v_limit_mjy']:.2f} & {e.get('class', 'nan')} \\\\"
        )
    if not rows_src:
        # Offline builds carry no per-target measurements (real data only), so emit an
        # honest placeholder row rather than an empty tabular body — matching the wdReal*
        # macro convention and keeping the paper's tabular valid (an empty body would leave
        # the header \hline immediately followed by this trailing \hline, which LaTeX rejects
        # with "Misplaced \noalign").
        out.append(r"\multicolumn{6}{c}{\emph{Per-target rows appear in the real-data build.}} \\")
    out.append(r"\hline")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(out) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="WD-pulsar candidate radio survey (plan 41).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    m = run(args.out, offline=args.offline)
    m.pop("per_target", None)
    print(json.dumps(m, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
