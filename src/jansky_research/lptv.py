"""LPT catalogue v3 + the first multi-epoch Stokes-V forced photometry at all LPT positions (plan 44).

Two coupled legs on the long-period-transient (LPT) class.

**Catalogue (v3).** `data/lpt_sample.csv` (merged `lpt` slice) is extended in place with the three
2026 discoveries verified at GATE-0 (2026-07-08): ASKAP J142431.2-612611 (arXiv:2603.07857) and
the two VASTER sources ASKAP J165130.3-450520 and J170036.6-445758 (arXiv:2606.20067) --- N goes
13 -> 16. Coordinates were checked against each source-name convention (the discipline that
caught the Rea+2026 review's own 2225-vs-3225 s typo); this module re-runs the population
statistics at N=16 and, in particular, whether the hinted ~78-min WD-binary period boundary moves
(it does not: the split stays insignificant --- reported, not spun).

**Stokes-V (first).** The merged counterpart cross-match was Stokes-I only (VLASS/LoTSS); nobody
has done forced *circular-polarization* photometry at the LPT positions. RACS-low2 Paper VIII
(arXiv:2606.16182) published a BLIND V catalogue that did not target them, and per-source
polarization exists for a handful (GLEAM-X J1627 ~90% linear; ASKAP J1935 >70% circular in its
weak state; J1424 circular ~8%; CHIME J1634 ~100% circular) --- but no systematic multi-epoch V
survey (GATE-0 PASS). A persistent circularly-polarized counterpart would discriminate
coherent-emission models; the likely outcome at RACS snapshot depth is a systematic V-limit table
(low duty cycles), which is the deliverable, framed as such.

Reuse: `lpt.load_sample`/`population_table`/`period_split_stat` for the catalogue leg;
`stokesv.measure_circular_pol`/`handedness`/`classify_emitter` and the RACS leakage floor for the
V leg; `wdpulsar.synthetic_cutout` for the offline injection recover-a-known; the resumable
per-row-CSV CASDA driver pattern (`scripts/lptv_real.py`).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .lpt import load_sample, period_split_stat, population_table
from .stokesv import classify_emitter, handedness, measure_circular_pol

__all__ = [
    "lpt_positions",
    "catalogue_stats",
    "injection_roundtrip",
    "summarize_v_sweep",
    "handedness_changes",
    "run",
]

SWEEP_CSV = Path("results/lptv_realtargets.csv")
LEAKAGE_FRAC = 0.006  # ASKAP on-axis Stokes-I->V leakage level (the V detection floor)
DET_NSIGMA = 5.0
LIMIT_NSIGMA = 3.0


def lpt_positions() -> list[dict]:
    """The v3 LPT target list (name, RA, Dec, period) from the vendored catalogue."""
    s = load_sample()
    return [
        {
            "name": str(s["name"][i]),
            "ra_deg": float(s["ra"][i]),
            "dec_deg": float(s["dec"][i]),
            "period_s": float(s["period_s"][i]),
        }
        for i in range(s["name"].size)
    ]


def catalogue_stats() -> dict:
    """Re-run the population statistics at v3 (N=16): headline table + the period-split test."""
    s = load_sample()
    pt = population_table(s)
    split = period_split_stat(s["period_s"], s["is_wd_binary"])
    # the 2026 additions, flagged by discovery year
    new = [str(n) for n, y in zip(s["name"], s["year"], strict=True) if int(y) == 2026]
    return {
        **pt,
        "period_split_delta_log": split["delta_log_median"],
        "period_split_p": split["p_perm"],
        "binary_boundary_significant": bool(split["p_perm"] < 0.05),
        "n_added_2026": len([n for n in new if "J1424" in n or "J1651" in n or "J1700" in n]),
    }


def injection_roundtrip(*, seed: int = 0) -> dict:
    """Offline recover-a-known: a synthetic V point source is measured back; a blank field is a limit."""
    from .wdpulsar import synthetic_cutout

    i_in, v_in = 6.0, -3.6  # a highly-circular (|V|/I=0.6) coherent-emitter-like injection
    img_i, img_v, w = synthetic_cutout(i_mjy=i_in, v_mjy=v_in, rms_mjy=0.2, seed=seed)
    hit = measure_circular_pol(img_i, img_v, w, 245.447, -22.886)
    blank_i, blank_v, wb = synthetic_cutout(i_mjy=0.0, v_mjy=0.0, rms_mjy=0.2, seed=seed + 1)
    miss = measure_circular_pol(blank_i, blank_v, wb, 245.447, -22.886)
    return {
        "v_in": v_in,
        "v_out": round(hit["v_peak"], 3),
        "frac_out": round(hit["frac_pol"], 3),
        "class": classify_emitter(hit["v_peak"], hit["i_peak"]),
        "handedness": handedness(hit["v_peak"]),
        "blank_v_sig": round(miss["v_peak"] / miss["v_rms"], 2)
        if miss["v_rms"] > 0
        else float("nan"),
    }


def handedness_changes(epoch_rows: list[dict]) -> str | None:
    """Inter-epoch V-sign change across a target's DETECTED epochs (LCP<->RCP), else None.

    Only epochs with a leakage-vetted V detection count; a handedness flip between epochs is a
    strong coherent-emission signature (but see `stokesv.handedness`: the absolute V sign
    convention varies between RACS pipeline versions --- flagged, not physically interpreted here).
    """
    signs = {
        handedness(float(r["v_mjy"]))
        for r in epoch_rows
        if r.get("v_det") and r.get("v_mjy") not in (None, "", "nan")
    }
    return "flip" if len(signs) > 1 else None


def summarize_v_sweep(
    csv_path: str | Path,
    targets: list[dict],
    *,
    det_sigma: float = DET_NSIGMA,
    limit_sigma: float = LIMIT_NSIGMA,
    confuse_offset: float = 4.0,
    confuse_iratio: float = 10.0,
    secure_offset: float = 2.0,
) -> dict:
    """Reduce the resumable V-sweep CSV (one row per LPT position x RACS epoch) into the V table.

    Per LPT: epochs measured, the max-|V|-significance measurement, a leakage-vetted V-detection
    flag (|V| > det_sigma*rms AND > LEAKAGE_FRAC*I), the deepest limit_sigma V upper limit, the
    circular class where detected, and any inter-epoch handedness flip.

    Two vetoes gate a *believable* detection: (1) instrumental leakage (|V| > LEAKAGE_FRAC*I);
    (2) **confusion** --- the peak is flagged suspect when it sits far from the catalogue position
    (offset > ``confuse_offset``) AND its Stokes-I is a gross outlier versus the source's own other
    epochs (> ``confuse_iratio`` x their median). Forced photometry at a fixed position can lock
    onto a nearby bright source; a detection that is simultaneously off-centre and 10-100x brighter
    than every other epoch of the same target is a confusing source, not the LPT. Among the
    surviving detections, ``secure_offset`` splits **secure** (on-centre, offset < it) from
    **candidate** (offset >= it but not confusion-vetoed): a several-arcsec offset that is not a
    gross I outlier is not confusion, but the polarized peak's association with the target is less
    certain. Thresholds (4", 10x, 2") were set with this sweep in view; they are heuristic, and the
    conjunction keeps them from being razor-tuned (the confused peak is far from the boundary in
    both offset and I-ratio). A blind snapshot survey of duty-cycled bursters is expected to yield
    mostly limits, but genuine bursts can be caught --- so detections are reported at their honest
    confidence, suspects flagged, not silently kept.
    """
    import csv

    rows = [r for r in csv.DictReader(open(csv_path)) if r.get("v_mjy") not in ("", "nan", None)]
    by: dict[str, list] = {}
    for r in rows:
        r["v_det"] = abs(float(r["v_mjy"])) >= det_sigma * float(r["e_v"]) and abs(
            float(r["v_mjy"])
        ) > LEAKAGE_FRAC * max(float(r.get("i_mjy", 0.0) or 0.0), 0.0)
        by.setdefault(r["name"], []).append(r)
    out = []
    for t in targets:
        rs = by.get(t["name"], [])
        e: dict = {
            "name": t["name"],
            "period_min": round(t["period_s"] / 60.0, 1),
            "n_epochs": len(rs),
        }
        if rs:
            v_sig = [abs(float(r["v_mjy"])) / float(r["e_v"]) for r in rs]
            k = int(np.argmax(v_sig))
            best = rs[k]
            v_det = bool(best["v_det"])
            # confusion veto on the detection epoch: off-centre AND a gross I outlier vs the source
            i_all = [float(r.get("i_mjy", 0.0) or 0.0) for r in rs]
            i_med = float(np.median(i_all)) if i_all else 0.0
            i_det = float(best.get("i_mjy", 0.0) or 0.0)
            offset = float(best.get("offset_arcsec", 0.0) or 0.0)
            i_ratio = i_det / i_med if i_med > 0 else float("inf")
            suspect = bool(v_det and offset > confuse_offset and i_ratio > confuse_iratio)
            believable = bool(v_det and not suspect)
            secure = bool(believable and offset < secure_offset)
            e.update(
                {
                    "v_det": v_det,
                    "believable": believable,
                    "secure": secure,  # on-centre; a believable-but-offset det is a candidate
                    "suspect_confusion": suspect,
                    "v_mjy": round(float(best["v_mjy"]), 3) if v_det else None,
                    "offset_arcsec": round(offset, 2) if v_det else None,
                    "i_outlier_ratio": round(i_ratio, 1) if v_det else None,
                    "v_limit_mjy": round(limit_sigma * min(float(r["e_v"]) for r in rs), 3),
                    "class": classify_emitter(float(best["v_mjy"]), i_det) if believable else "nan",
                    "handedness_change": handedness_changes(rs),
                }
            )
        out.append(e)
    measured = [e for e in out if e["n_epochs"] > 0]
    dets = [e for e in measured if e.get("believable")]
    secure_dets = [e for e in dets if e.get("secure")]
    candidate_dets = [e for e in dets if not e.get("secure")]
    suspects = [e for e in measured if e.get("suspect_confusion")]
    return {
        "per_target": out,
        "n_targets": len(targets),
        "n_measured": len(measured),
        "n_v_detections": len(dets),  # believable = secure + candidate (both confusion-vetted)
        "n_v_secure": len(secure_dets),  # on-centre, high-confidence
        "n_v_candidate": len(candidate_dets),  # believable but offset -> candidate
        "n_suspect_confusion": len(suspects),
        "v_detection_names": sorted(e["name"] for e in dets),
        "v_secure_names": sorted(e["name"] for e in secure_dets),
        "suspect_names": sorted(e["name"] for e in suspects),
        "n_handedness_flips": sum(1 for e in dets if e.get("handedness_change") == "flip"),
        "median_v_limit_mjy": round(float(np.median([e["v_limit_mjy"] for e in measured])), 3)
        if measured
        else None,
    }


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: v3 catalogue stats + injection recover-a-known; real: reduce the V sweep CSV."""
    import json

    cat = catalogue_stats()
    rt = injection_roundtrip()
    metrics: dict = {
        "source": "LPT v3 catalogue (16 LPTs) + synthetic V injection",
        "is_real": not offline,
        "n_lpt": cat["n_lpt"],
        "n_wd_binary": cat["n_wd_binary"],
        "median_period_min": cat["median_period_min"],
        "period_split_p": cat["period_split_p"],
        "binary_boundary_significant": cat["binary_boundary_significant"],
        "injection": rt,
    }
    per_target = None
    if not offline:  # pragma: no cover - needs the CASDA V-sweep CSV
        summary = summarize_v_sweep(SWEEP_CSV, lpt_positions())
        per_target = summary.pop("per_target")
        metrics.update(summary)
        metrics["source"] = "LPT v3 (16 LPTs) + RACS forced Stokes-V photometry (CASDA)"

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "lptv_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    if per_target is not None:  # pragma: no cover - real leg
        (op / "results" / "lptv_vtable.json").write_text(json.dumps(per_target, indent=2) + "\n")
        metrics["per_target"] = per_target
    _figure(metrics, op / "papers" / "lptv" / "figures")
    _write_macros(metrics, op / "papers" / "lptv" / "generated" / "macros.tex")
    _write_v_table(metrics, op / "papers" / "lptv" / "generated" / "v_table.tex")
    return metrics


def _figure(m: dict, out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    s = load_sample()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.9))
    pmin = s["period_s"] / 60.0  # the v3 period distribution
    ax1.hist(np.log10(pmin), bins=10, color="0.7")
    ax1.axvline(np.log10(78.0), color="C3", ls="--", label="~78-min hint")
    ax1.set(xlabel=r"$\log_{10}$ period (min)", ylabel="LPTs", title="v3 period distribution")
    ax1.legend(fontsize=8)
    per = m.get("per_target") or []
    lims = [e["v_limit_mjy"] for e in per if e.get("v_limit_mjy") is not None]
    if lims:
        ax2.hist(lims, bins=12, color="C0")
        ax2.set(xlabel=r"3$\sigma$ V limit (mJy)", ylabel="LPTs", title="V limit census")
    else:
        ax2.text(0.5, 0.5, "real V sweep pending", ha="center", va="center")
        ax2.set_axis_off()
    fig.tight_layout()
    fig.savefig(out / "lptv.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    def g(d: str, key: str | None = None) -> str:
        sub = m.get(d)
        v = sub.get(key) if isinstance(sub, dict) and key is not None else sub
        if v is None:
            return "--"
        return "--" if isinstance(v, float) and not np.isfinite(v) else str(v)

    pref = "lvReal" if m.get("is_real") else "lvSyn"
    lines = [
        "% Auto-generated by jansky_research.lptv._write_macros -- do not edit.",
        "% Synthetic (lvSyn*) and real (lvReal*) namespaces are BOTH always emitted; the",
        "% inactive namespace holds placeholders, so synthetic numbers can never masquerade",
        "% under lvReal* (an offline rebuild resets lvReal* to placeholders by design).",
        rf"\newcommand{{\lvSource}}{{{m['source']}}}",
        rf"\newcommand{{\lvNLpt}}{{{m['n_lpt']}}}",
        rf"\newcommand{{\lvNWdBinary}}{{{m['n_wd_binary']}}}",
        rf"\newcommand{{\lvMedianPeriodMin}}{{{m['median_period_min']}}}",
        rf"\newcommand{{\lvPeriodSplitP}}{{{m['period_split_p']}}}",
        rf"\newcommand{{\lvInjClass}}{{{m['injection']['class'].replace('_', ' ')}}}",
    ]
    for ns in ("lvSyn", "lvReal"):
        live = ns == pref
        for macro, d, key in (
            ("NMeasured", "n_measured", None),
            ("NVDet", "n_v_detections", None),
            ("NVSecure", "n_v_secure", None),
            ("NVCandidate", "n_v_candidate", None),
            ("NSuspect", "n_suspect_confusion", None),
            ("MedianVLimit", "median_v_limit_mjy", None),
            ("NFlips", "n_handedness_flips", None),
        ):
            lines.append(rf"\newcommand{{\{ns}{macro}}}{{{g(d, key) if live else '--'}}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _write_v_table(m: dict, path: str | Path) -> None:
    """Per-LPT V detection/limit rows as LaTeX (the full v3 table; JSON carries the machine copy)."""
    out = [
        "% Auto-generated by jansky_research.lptv._write_v_table -- do not edit.",
        "% Columns: name & period (min) & epochs & V (mJy or 3sig limit) & class & handedness flip",
    ]
    per_target = m.get("per_target") or []
    if not per_target:  # offline / CI: no real V sweep yet -- a buildable placeholder row
        out.append(r"(real Stokes-V sweep pending) & -- & -- & -- & -- & -- \\")
        out.append(r"\hline")
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("\n".join(out) + "\n")
        return
    for e in per_target:
        name = e["name"].replace("_", r"\_")
        if e["n_epochs"] == 0:
            out.append(
                f"{name} & {e['period_min']} & 0 & \\multicolumn{{3}}{{c}}{{uncovered}} \\\\"
            )
            continue
        flip = e.get("handedness_change") or "--"
        if e.get("believable"):  # a detection: show the V value + circular class
            cls = e.get("class", "nan").replace("_", " ")
            if not e.get("secure"):  # offset -> candidate association
                cls += f" (cand., {e['offset_arcsec']:.1f}\\arcsec)"
            v_s = f"{e['v_mjy']:.2f}"
        elif e.get("suspect_confusion"):  # off-centre I outlier: report the limit, flag confusion
            v_s, cls = f"$<${e['v_limit_mjy']:.2f}", "confused"
        else:  # non-detection: deepest 3-sigma limit
            v_s, cls = f"$<${e['v_limit_mjy']:.2f}", "--"
        out.append(f"{name} & {e['period_min']} & {e['n_epochs']} & {v_s} & {cls} & {flip} \\\\")
    out.append(r"\hline")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(out) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="LPT v3 catalogue + Stokes-V forced photometry.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    m = run(args.out, offline=args.offline)
    m.pop("per_target", None)
    print(json.dumps(m, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
