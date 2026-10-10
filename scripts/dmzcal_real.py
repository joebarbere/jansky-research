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
from scipy import stats

from jansky_research import dmzcal

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
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
        z_true = float(b["z"])  # provenance applies primary-source corrections
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
            rec["dm_pit"] = float(dmzcal.dm_pit(t, np.array([dm_eg]), np.array([z_true]))[0])
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


HOFFMANN_DSA_LIMIT = 183.0  # DM_obs - DM_ISM; arXiv:2408.04878v2 Sec. 2.3, main.tex l.221
DM_LIMITS = [HOFFMANN_DSA_LIMIT, 250.0, 300.0, 400.0, 500.0, float("inf")]
DM_BINS = [0.0, HOFFMANN_DSA_LIMIT, 300.0, 500.0, float("inf")]  # un-nested (GATE-2 r2 N5)
SUBSETS = {  # descriptive telescope subsets (GATE-2 r2 N6: refreshed on corrected labels)
    "CHIME+DSA": {"CHIME", "DSA"},
    "DSA+ASKAP": {"DSA", "ASKAP"},
    "CHIME": {"CHIME"},
    "non-CHIME": {"DSA", "ASKAP", "MeerKAT"},
}


def _obs_minus_ism(r: dict) -> float:
    """DM_obs - DM_ISM, from whichever form an estimator's record carries."""
    if "dm_eg" in r:
        return float(r["dm_eg"]) + dmzcal.HOFFMANN_EMIN25.dm_halo
    return float(r["dm"]) - float(r["dm_ism"])


def _describe(pits: np.ndarray) -> dict:
    if pits.size == 0:
        return {"n": 0}
    return {
        "n": int(pits.size),
        "coverage": {str(q): dmzcal.coverage(pits, q) for q in (0.68, 0.95)},
        "median": float(np.median(pits)),
        "ks_p": float(stats.kstest(pits, "uniform").pvalue),
        "tails68": dmzcal.tail_fractions(pits),
    }


def _planted(rule: dmzcal.CalibrationRule, rows: list[dict]) -> dict:
    arr = np.array([r["planted_pits"] for r in rows if "planted_pits" in r])  # (n, reps)
    if arr.size == 0:
        return {"n": 0}
    passes = [rule.evaluate(arr[:, k])["passed"] for k in range(arr.shape[1])]
    return {
        "n": int(arr.shape[0]),
        "reps": int(arr.shape[1]),
        "first_rep_passed": bool(passes[0]),
        "pass_rate": float(np.mean(passes)),
    }


def _post_hoc(
    rule: dmzcal.CalibrationRule,
    est: dict[str, list[dict]],
    cert_names: set[str],
    expected_pass: float,
) -> dict:
    """GATE-2 round-1 checks. ALL POST HOC: chosen after the real verdicts were seen."""
    out: dict[str, Any] = {"label": "POST HOC (GATE-2 round 1); no frozen verdicts"}
    for name, rows in est.items():
        c = [r for r in rows if r["name"] in cert_names and "pit" in r]
        blk: dict[str, Any] = {
            "tails68_z": dmzcal.tail_fractions(np.array([r["pit"] for r in c])),
            "sign_below_half": int(sum(r["pit"] < 0.5 for r in c)),
        }
        if all("dm_pit" in r for r in c) and c:
            blk["dm_given_z"] = _describe(np.array([r["dm_pit"] for r in c]))
        lim = {}
        for L in DM_LIMITS:
            sel = [r for r in c if _obs_minus_ism(r) < L]
            lim["all" if L == float("inf") else f"<{L:g}"] = _describe(
                np.array([r["pit"] for r in sel])
            )
        blk["by_dm_obs_minus_ism_limit_NESTED"] = lim
        bins = {}
        for lo, hi in zip(DM_BINS[:-1], DM_BINS[1:], strict=True):
            sb = np.array([r["pit"] for r in c if lo <= _obs_minus_ism(r) < hi])
            d = _describe(sb)
            if sb.size:
                d["frac_below_half"] = float(np.mean(sb < 0.5))
            bins[f"{lo:g}-{hi:g}"] = d
        blk["by_dm_obs_minus_ism_bin"] = bins
        x = np.array([_obs_minus_ism(r) for r in c])
        y = np.array([r["pit"] for r in c])
        rho = stats.spearmanr(x, y)
        blk["trend_spearman_pit_vs_dm"] = {
            "rho": float(rho.statistic),
            "p": float(rho.pvalue),
            "n": int(x.size),
        }
        blk["subsets"] = {
            k: _describe(np.array([r["pit"] for r in c if r["telescope"] in tels]))
            for k, tels in SUBSETS.items()
        }
        if any("pit_halo" in r for r in c):
            blk["halo_sensitivity"] = {
                h: _describe(np.array([r["pit_halo"][h] for r in c if h in r.get("pit_halo", {})]))
                for h in ("25", "75")
            }
        if any("planted_pits" in r for r in c):
            blk["planted_truth"] = _planted(rule, c) | {"expected_pass_rate": expected_pass}
        out[name] = blk
    pr = {r["name"]: r for r in _load("dmzcal_e1_zdm_pruned.json")["bursts"]}
    pp = np.array([pr[n]["pit"] for n in sorted(cert_names) if "pit" in pr.get(n, {})])
    out["E1_survey_files_pruned"] = _describe(pp) | {
        "verdict": _verdict(rule, pp)["verdict"],
        "max_abs_pit_shift": float(
            max(
                abs(pr[r["name"]]["pit"] - r["pit"])
                for r in est["E1_zdm_HoffmannEmin25"]
                if r["name"] in cert_names and "pit" in r
            )
        ),
    }
    e1c = [r for r in est["E1_zdm_HoffmannEmin25"] if r["name"] in cert_names and "pit" in r]
    pa = np.array([r.get("pit_alt", r["pit"]) for r in e1c])
    out["E1_unresolved_alternatives_applied"] = _describe(pa) | {
        "verdict": _verdict(rule, pa)["verdict"],
        "changed": {
            r["name"]: [r["pit"], r["pit_alt"]]
            for r in e1c
            if "pit_alt" in r and abs(r["pit_alt"] - r["pit"]) > 1e-12
        },
    }
    out["host_magnitude_stratification"] = (
        "SKIPPED: pre-stated as conditional ('if the P(O|x) column permits'); host "
        "magnitudes are available only for the 19 CHIME/Leung hosts (m_r), not for the "
        "DSA/ASKAP/MeerKAT ones, so no sample-wide stratification is possible."
    )
    return out


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
        "post_hoc": _post_hoc(
            rule,
            est,
            cert_names,
            1 - _load("dmzcal_controls.json")["frozen"]["C0"]["false_fail_rate"],
        ),
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
