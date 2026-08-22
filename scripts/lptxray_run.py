#!/usr/bin/env python3
"""Plan 93: analyse the staged X-ray cones and write results/lptxray_metrics.json.

Offline. Consumes ``data/lptxray/cones.json`` (from ``scripts/lptxray_fetch.py``) and, when
present, ``data/lptxray/coverage.json`` (from ``scripts/lptxray_coverage.py``) plus this
repo's own radio photometry in ``results/wdpulsar_realtargets.csv``.

    uv run python scripts/lptxray_run.py
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jansky_research import lptxray as lx  # noqa: E402

RADIO_SIGMA = 5.0


def radio_leg(path: Path) -> dict:
    """Per-candidate radio limits from the `wdpulsar` forced photometry."""
    if not path.exists():
        return {"available": False}
    per: dict[str, list[dict]] = {}
    for r in csv.DictReader(open(path)):
        try:
            rec = {
                "i": float(r["i_mjy"]),
                "e_i": float(r["e_i"]),
                "v": float(r["v_mjy"]),
                "e_v": float(r["e_v"]),
            }
        except (ValueError, KeyError):
            continue
        if not all(v == v for v in rec.values()):
            continue
        per.setdefault(r["name"], []).append(rec)
    out = {}
    for name, eps in per.items():
        best_i = min(e["e_i"] for e in eps)
        best_v = min(e["e_v"] for e in eps)
        det = any(e["i"] > RADIO_SIGMA * e["e_i"] for e in eps)
        out[name] = {
            "n_epochs": len(eps),
            "i_limit_mjy": round(RADIO_SIGMA * best_i, 4),
            "v_limit_mjy": round(RADIO_SIGMA * best_v, 4),
            "detected": det,
        }
    return {"available": True, "per_source": out, "sigma": RADIO_SIGMA}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cones", type=Path, default=Path("data/lptxray/cones.json"))
    ap.add_argument("--coverage", type=Path, default=Path("data/lptxray/coverage.json"))
    ap.add_argument("--radio", type=Path, default=Path("results/wdpulsar_realtargets.csv"))
    ap.add_argument("--out", type=Path, default=Path("results/lptxray_metrics.json"))
    args = ap.parse_args(argv)

    blob = json.loads(args.cones.read_text())
    cone_arcmin = float(blob["cone_arcmin"])
    records = [
        lx.summarize_position(v["meta"], v["cats"], cone_arcmin=cone_arcmin)
        for v in blob["cones"].values()
        if v.get("cats")
    ]
    lpts = [r for r in records if r["sample"] == "lpt"]
    cands = [r for r in records if r["sample"] == "wdcand"]

    # --- the guard: what recall does this cross-match have on KNOWN detections? ---
    recall = {"lpt": lx.catalogue_recall(lpts), "wdcand": lx.catalogue_recall(cands)}

    # --- measured chance rate, two independent ways ---
    chance: dict = {}
    for cat in lx.CATALOGS:
        exp, hits, trials, dens = [], 0, 0, []
        for r in records:
            c = r["catalogs"].get(cat, {})
            if "chance_expected" not in c:
                continue
            exp.append(c["chance_expected"])
            dens.append(c["density_per_sqarcmin"])
            hits += c["shift_trials"]["n_hits"]
            trials += c["shift_trials"]["n_trials"]
        chance[cat] = {
            "radius_arcsec": lx.MATCH_RADIUS_ARCSEC[cat],
            "median_density_per_sqarcmin": round(statistics.median(dens), 6) if dens else None,
            "mean_chance_expected": round(sum(exp) / len(exp), 5) if exp else None,
            "shift_trials": trials,
            "shift_hits": hits,
            "shift_false_match_rate": round(hits / trials, 5) if trials else None,
        }

    # --- accretion split among the white-dwarf candidates ---
    acc = [r for r in cands if r["type"] in lx.ACCRETING_TYPES]
    non = [r for r in cands if r["type"] not in lx.ACCRETING_TYPES]
    k_acc = sum(r["matched"] for r in acc)
    k_non = sum(r["matched"] for r in non)
    comparison = lx.fraction_comparison(k_acc, len(acc), k_non, len(non))

    # --- the circularity check, and it is free to fail ---
    # Pelisoli et al. say the X-ray column "defined follow-up priority": sources with a
    # literature X-ray detection were preferentially followed up spectroscopically, and the
    # spectroscopy is what produced the `type` classification. So "accreting systems are
    # X-ray bright" is partly an artefact of how the types were assigned. Restricting to
    # candidates with NO literature X-ray flag removes that path -- those were never
    # X-ray-prioritised, so any surviving split is not circular.
    blind_acc = [r for r in acc if not r["lit_xray_detected"]]
    blind_non = [r for r in non if not r["lit_xray_detected"]]
    blind: dict = {
        "rationale": (
            "Candidates carrying no literature X-ray flag, i.e. those never prioritised for "
            "spectroscopic follow-up on X-ray grounds. The type classification for this "
            "subset cannot have been driven by a prior X-ray detection."
        ),
        "n_accreting": len(blind_acc),
        "n_non_accreting": len(blind_non),
    }
    if blind_acc and blind_non:
        blind["comparison"] = lx.fraction_comparison(
            sum(r["matched"] for r in blind_acc),
            len(blind_acc),
            sum(r["matched"] for r in blind_non),
            len(blind_non),
        )
    else:
        blind["comparison"] = None
        blind["note"] = "one arm is empty; the check cannot be run on this sample"

    # --- robustness: the same split under cuts that could each have broken it ---
    # The four catalogues have different footprints, so an aggregate "any match" could in
    # principle be an artefact of which group sits where. Two cuts remove that freedom: the
    # eROSITA-DE half is a common-footprint subsample, and 2RXS alone is a single all-sky
    # catalogue with no footprint asymmetry whatever.
    def _in_erosita_de(rec: dict) -> bool:
        import astropy.units as u
        from astropy.coordinates import SkyCoord

        return 180.0 <= SkyCoord(rec["ra"] * u.deg, rec["dec"] * u.deg).galactic.l.deg < 360.0

    robustness: dict = {}
    de = [r for r in cands if _in_erosita_de(r)]
    de_acc = [r for r in de if r["type"] in lx.ACCRETING_TYPES]
    de_non = [r for r in de if r["type"] not in lx.ACCRETING_TYPES]
    if de_acc and de_non:
        robustness["erosita_de_half_only"] = {
            "rationale": "common-footprint subsample: every source has eROSITA coverage",
            **lx.fraction_comparison(
                sum(r["matched"] for r in de_acc),
                len(de_acc),
                sum(r["matched"] for r in de_non),
                len(de_non),
            ),
        }
    for cat in lx.CATALOGS:
        ka = sum(r["catalogs"][cat].get("matched", False) for r in acc)
        kn = sum(r["catalogs"][cat].get("matched", False) for r in non)
        robustness[f"{cat}_only"] = lx.fraction_comparison(ka, len(acc), kn, len(non))
    robustness["rass2rxs_only"]["rationale"] = (
        "2RXS is all-sky, so this leg carries no footprint asymmetry between the groups."
    )

    radio = radio_leg(args.radio)
    radio_by_group: dict = {}
    if radio.get("available"):
        for label, group in (("accreting", acc), ("non_accreting", non)):
            names = [r["name"] for r in group if r["name"] in radio["per_source"]]
            lims = [radio["per_source"][n]["i_limit_mjy"] for n in names]
            radio_by_group[label] = {
                "n_in_group": len(group),
                "n_with_radio_limit": len(names),
                "n_radio_detected": sum(radio["per_source"][n]["detected"] for n in names),
                "median_i_limit_mjy": round(statistics.median(lims), 4) if lims else None,
            }

    coverage: dict = {"available": False}
    if args.coverage.exists():
        cblob = json.loads(args.coverage.read_text())["coverage"]
        rows = []
        for name, e in cblob.items():
            obs = [o for v in e["logs"].values() for o in v.get("obs", [])]
            cls = lx.classify_coverage(float(e["meta"]["ra"]), float(e["meta"]["dec"]), obs)
            rows.append(
                {
                    "name": name,
                    "lit_xray": e["meta"]["lit_xray"],
                    "n_obs_in_field": len(obs),
                    **cls,
                }
            )
        coverage = {
            "available": True,
            "n_total": len(rows),
            "n_any_coverage": sum(r["n_obs_in_field"] > 0 for r in rows),
            "n_targeted": sum(r["n_targeted"] > 0 for r in rows),
            "n_targeted_without_published_detection": sum(
                r["n_targeted"] > 0 and r["lit_xray"] == "no" for r in rows
            ),
            "note": (
                "A serendipitous pointing covers the position while observing something "
                "else; it is real evidence but is not dedicated attention. 'Without a "
                "published detection' does NOT mean unanalysed -- several of these "
                "observations are the source of published upper limits."
            ),
            "per_source": rows,
        }

    payload = {
        "slice": "lptxray",
        "stage": "real-xray-crossmatch",
        "is_real": True,
        "source": (
            "HEASARC cone searches: xmmssc (5XMM-DR15), rass2rxs (2RXS), erass1main "
            "(eROSITA-DE eRASS1), csc (Chandra CSC 2.1.1); one cone per position, "
            f"{cone_arcmin} arcmin."
        ),
        "method": (
            "Nearest-source association inside a per-catalogue radius set by that "
            "catalogue's positional accuracy. The chance-coincidence rate is MEASURED two "
            "ways from the same cached cone -- the field source density in a 2-10 arcmin "
            "annulus, and rigid position-shift trials at 3/5/7 arcmin in 8 azimuths, every "
            "trial disc required to lie inside the cone. Fractions are reported with both "
            "their difference and their ratio."
        ),
        "gate0": {
            "lpt_leg": "CUT",
            "reason": (
                "All three LPTs with a published X-ray detection are absent from every "
                "serendipitous source catalogue while each has dedicated pointed coverage "
                "taken after those catalogues were built. Catalogue recall on the LPT "
                "sample is 1/3, so a catalogue-derived LPT detection fraction measures the "
                "catalogues, not the sources."
            ),
        },
        "recall_guard": recall,
        "chance_rate": chance,
        "n_lpt": len(lpts),
        "n_candidates": len(cands),
        "accretion_comparison": {
            "accreting_types": list(lx.ACCRETING_TYPES),
            "xray": comparison,
            "xray_no_prior_flag": blind,
            "robustness": robustness,
            "radio": radio_by_group,
        },
        "pointed_coverage": coverage,
        "caveats": [
            "Selection is not physics: the candidates were selected photometrically from "
            "Gaia+WISE and the LPTs by being radio-detected, so the two lists are not two "
            "draws from one population and their X-ray fractions are not directly "
            "comparable. This is why the LPT leg is reported separately and not differenced "
            "against the candidates.",
            "The four catalogues have different bands, depths and footprints; eROSITA-DE "
            "covers only l = 180-360 deg (9 of 16 LPTs, 38 of 56 candidates). Fluxes are "
            "reported per catalogue and never combined, which would need a spectral "
            "assumption this sample cannot support.",
            "The LPTs lie mostly at b ~ 0 at kpc distances while the candidates are nearby "
            "and high-latitude, so absorption and distance differ by orders of magnitude "
            "between the two lists -- a further reason the two legs are not differenced.",
            "A literature X-ray flag is a heterogeneous quantity: it records that SOMEONE "
            "detected the source, at unknown depth. It is used here only to measure recall, "
            "never as a detection fraction.",
            "The accretion split is partly circular by construction: Pelisoli et al. used "
            "their X-ray column to set spectroscopic follow-up priority, and that "
            "spectroscopy produced the type classifications, so X-ray-detected candidates "
            "were likelier to be classified at all. `xray_no_prior_flag` re-runs the "
            "comparison on candidates with no literature X-ray flag, where that path is "
            "cut; read that number before the headline one.",
        ],
        "records": records,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(".json.part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(args.out)

    print(f"wrote {args.out}")
    print(
        f"  recall  LPT      : {recall['lpt']['n_recovered']}/"
        f"{recall['lpt']['n_known_detections']} usable={recall['lpt']['usable']}"
    )
    print(
        f"  recall  candidates: {recall['wdcand']['n_recovered']}/"
        f"{recall['wdcand']['n_known_detections']} usable={recall['wdcand']['usable']}"
    )
    print(
        f"  X-ray accreting   : {comparison['k_a']}/{comparison['n_a']} "
        f"({100 * comparison['frac_a']:.1f}%)"
    )
    print(
        f"  X-ray other       : {comparison['k_b']}/{comparison['n_b']} "
        f"({100 * comparison['frac_b']:.1f}%)"
    )
    print(f"  ratio {comparison['ratio']}  difference {comparison['difference_pp']} pp")
    return 0


if __name__ == "__main__":
    sys.exit(main())
