"""Plan 99 step 3: real-data coverage audit of E1-E3 (+ E4 point accuracy) and C2.

Inputs (all committed): ``results/dmzcal_provenance.json`` (sample), ``dmzcal_controls.json``
(frozen rule bands, N = 36), ``dmzcal_e1_zdm.json`` / ``dmzcal_e3.json`` / ``dmzcal_e4.json``
(per-burst outputs of the out-of-repo drivers). E2 is computed here.

Verdict wording follows the plan's post-C1 asymmetric rule: a failure is reported by name; a
pass is reported only as "consistent with calibration at a sensitivity that cannot exclude a
2x host-mean error" (``CONSISTENT_LOW_POWER``), never as CALIBRATED.

Run: ``uv run python scripts/dmzcal_real.py``
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from jansky_research import dmzcal

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
Z_OVERRIDE = {"FRB20231201A": 0.1119}
CALEB_DM_RANGE_95 = [1.628, 3.397]  # arXiv:2508.01648, "Clues" section
PASS_WORD = "CONSISTENT_LOW_POWER"


def _load(name: str) -> dict:
    return json.loads((R / name).read_text())


def _rule() -> dmzcal.CalibrationRule:
    c = _load("dmzcal_controls.json")["frozen"]["rule"]
    return dmzcal.CalibrationRule(
        bands={float(q): tuple(b) for q, b in c["bands"].items()}, ks_p=c["ks_p"]
    )


def _e2(prov: dict) -> dict:
    t = dmzcal.build_table(dmzcal.HOFFMANN_EMIN25)
    rows = list(prov["bursts"])
    rows.append(prov["recover_a_known"] | {"side": "recover_a_known", "secure_host": True})
    out = []
    for b in rows:
        z_true = Z_OVERRIDE.get(b["name"], float(b["z"]))
        dm_eg = float(b["DM"]) - float(b["DMISM"]) - dmzcal.HOFFMANN_EMIN25.dm_halo
        rec = {
            "name": b["name"],
            "telescope": b["telescope"],
            "side": b["side"],
            "spec_z": b["spec_z"],
            "secure_host": b.get("secure_host", True),
            "repeater": b.get("repeater", False),
            "z_true": z_true,
            "dm_eg": dm_eg,
        }
        if dm_eg <= 0:
            rec["excluded"] = "DM_EG <= 0"
        else:
            rec["pit"] = float(dmzcal.pit(t, np.array([dm_eg]), np.array([z_true]))[0])
            lo68, hi68 = dmzcal.central_interval(t, dm_eg, 0.68)
            lo95, hi95 = dmzcal.central_interval(t, dm_eg, 0.95)
            rec["ci68"], rec["ci95"] = [lo68, hi68], [lo95, hi95]
            rec["z_median"] = dmzcal.central_interval(t, dm_eg, 0.0)[0]
        out.append(rec)
    return {"bursts": out, "table": t}


def _verdict(rule: dmzcal.CalibrationRule, pits: np.ndarray) -> dict:
    ev = rule.evaluate(pits)
    ev["verdict"] = PASS_WORD if ev["passed"] else ev["verdict"]
    ev["n"] = int(pits.size)
    ev["bands_frozen_for_n"] = 36
    return ev


def _subset(bursts: list[dict], **want: object) -> list[dict]:
    return [b for b in bursts if all(b.get(k) == v for k, v in want.items())]


def _arm(rule: dmzcal.CalibrationRule, rows: list[dict], verdict: bool) -> dict:
    ok = [r for r in rows if "pit" in r]
    pits = np.array([r["pit"] for r in ok])
    out: dict = {
        "n_in": len(rows),
        "n_used": len(ok),
        "excluded": [r["name"] for r in rows if "pit" not in r],
    }
    if not ok:
        return out
    out["coverage"] = {str(q): dmzcal.coverage(pits, q) for q in (0.68, 0.95)}
    out["median_pit"] = float(np.median(pits))
    out["frac_pit_below_half"] = float(np.mean(pits < 0.5))
    if verdict:
        out["rule"] = _verdict(rule, pits)
    return out


def _point(rows: list[dict], key: str) -> dict:
    d = np.array([(r[key] - r["z_true"]) / (1 + r["z_true"]) for r in rows if key in r])
    if d.size == 0:
        return {"n": 0}
    return {
        "n": int(d.size),
        "median_dz1z": float(np.median(d)),
        "nmad_dz1z": float(1.4826 * np.median(np.abs(d - np.median(d)))),
    }


def _c2(
    table: dmzcal.LikelihoodTable,
    cert_e2: list[dict],
    rule: dmzcal.CalibrationRule,
    reps: int = 500,
) -> dict:
    """Planted truth through the real pipeline: each certification burst's z_true replaced
    by a draw from its own E2 posterior (cell-uniform), then PIT recomputed."""
    dm = np.array([r["dm_eg"] for r in cert_e2 if "pit" in r])
    post = dmzcal.posterior(table, dm)
    w = dmzcal._node_widths(table.z)
    rng = np.random.default_rng(0)

    def one() -> np.ndarray:
        idx = np.array([rng.choice(table.z.size, p=post[:, j]) for j in range(dm.size)])
        z = table.z[idx] + (rng.random(dm.size) - 0.5) * w[idx]
        return dmzcal.pit(table, dm, z)

    first = rule.evaluate(one())
    rate = float(np.mean([rule.evaluate(one())["passed"] for _ in range(reps)]))
    return {
        "n": int(dm.size),
        "seed": 0,
        "first_draw": first,
        "passed": first["passed"],
        "pass_rate_over_reps": rate,
        "reps": reps,
        "expected_pass_rate": 1 - _load("dmzcal_controls.json")["frozen"]["C0"]["false_fail_rate"],
    }


def main() -> None:
    prov = _load("dmzcal_provenance.json")
    rule = _rule()
    e2 = _e2(prov)
    table = e2.pop("table")
    est = {
        "E1_zdm_HoffmannEmin25": _load("dmzcal_e1_zdm.json")["bursts"],
        "E2_minimal": e2["bursts"],
        "E3_fruitbat_Batten2021": _load("dmzcal_e3.json")["bursts"],
    }
    e4 = _load("dmzcal_e4.json")["bursts"]

    def cert(rows: list[dict]) -> list[dict]:
        return _subset(rows, side="certification", spec_z=True, secure_host=True)

    audit: dict[str, dict[str, Any]] = {}
    for name, rows in est.items():
        c = cert(rows)
        audit[name] = {
            "certification": _arm(rule, c, verdict=True),
            "by_telescope": {
                tel: _arm(rule, [r for r in c if r["telescope"] == tel], verdict=False)
                for tel in sorted({r["telescope"] for r in c})
            },
            "repeaters": _arm(rule, [r for r in c if r["repeater"]], verdict=False),
            "non_repeaters": _arm(rule, [r for r in c if not r["repeater"]], verdict=False),
            "in_sample_production_LABELLED": _arm(
                rule, _subset(rows, side="production", spec_z=True, secure_host=True), verdict=False
            ),
            "unflagged_spectrum_sensitivity": _arm(
                rule, _subset(rows, side="certification", spec_z=False), verdict=False
            ),
        }
    cert_names = {r["name"] for r in cert(e2["bursts"])}
    point = {
        "E1": _point(
            [r for r in est["E1_zdm_HoffmannEmin25"] if r["name"] in cert_names], "z_median"
        ),
        "E2": _point([r for r in e2["bursts"] if r["name"] in cert_names], "z_median"),
        "E3": _point(
            [r for r in est["E3_fruitbat_Batten2021"] if r["name"] in cert_names], "z_median"
        ),
        "E4": _point([r for r in e4 if r["name"] in cert_names], "z_point"),
        "E4_excluded": [r["name"] for r in e4 if r["name"] in cert_names and "excluded" in r],
    }
    rak: dict[str, Any] = {}
    for name, rows in est.items():
        r = next(x for x in rows if x["side"] == "recover_a_known")
        rak[name] = {k: r.get(k) for k in ("pit", "z_median", "ci68", "ci95", "excluded")}
        if "ci68" in r:
            rak[name]["in68"] = r["ci68"][0] <= r["z_true"] <= r["ci68"][1]
            rak[name]["in95"] = r["ci95"][0] <= r["z_true"] <= r["ci95"][1]
    r4 = next(x for x in e4 if x["side"] == "recover_a_known")
    rak["E4_point"] = r4.get("z_point")
    rak["caleb_dm_range_95"] = CALEB_DM_RANGE_95
    rak["z_true"] = 2.148
    out: dict[str, Any] = {
        "source": "real: plan 99 coverage audit of E1-E3 + E4 point accuracy on FRBs/FRB "
        "localized hosts (certification sample per dmzcal_provenance.json)",
        "plan": "plans/99-dmzcal-coverage.md step 3",
        "verdict_rule": "asymmetric (plan 99, post-C1): failures named; a pass is "
        f"{PASS_WORD}, never CALIBRATED",
        "c1_power": {"host_mean_x2": 0.1205, "F_x2": 0.951},
        "rule_bands": {str(q): list(b) for q, b in rule.bands.items()},
        "audit": audit,
        "point_accuracy_certification": point,
        "recover_a_known_FRB20240304B": rak,
        "C2_planted_truth_E2": _c2(table, cert(e2["bursts"]), rule),
        "e2_bursts": e2["bursts"],
    }
    (R / "dmzcal_metrics.json").write_text(json.dumps(out, indent=1) + "\n")
    for name, a in audit.items():
        c = a["certification"]
        print(
            name,
            "n",
            c["n_used"],
            "cov",
            c.get("coverage"),
            "medPIT",
            round(c.get("median_pit", float("nan")), 3),
            "->",
            c["rule"]["verdict"],
        )
    print(
        "C2 passed:",
        out["C2_planted_truth_E2"]["passed"],
        "rate",
        out["C2_planted_truth_E2"]["pass_rate_over_reps"],
    )


if __name__ == "__main__":
    main()
