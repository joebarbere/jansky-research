"""Plan 99 step 0: which localized FRBs may certify the p(z|DM) estimators.

Assembles the deduplicated localized sample from the FRBs/FRB repository at a
pinned commit, and splits it into the *production* side (bursts that entered the
fit of the zdm parameter state under audit) and the *certification* side (bursts
that could not have shaped it). Writes ``results/dmzcal_provenance.json``.

No posterior and no z-vs-DM comparison is computed here: this only counts.

Run: ``uv run python scripts/dmzcal_provenance.py``  (network on first run)
"""

from __future__ import annotations

import csv
import json
import re
import urllib.request
from collections import Counter
from pathlib import Path

from jansky_research.dmzcal import _date, _dm, _key, dedupe, provenance_side

ROOT = Path(__file__).resolve().parents[1]
FRB_COMMIT = "996fcda9b0b22431e3171208b4c8e1cf2798e823"  # FRBs/FRB HEAD, 2026-10-10
ZDM_COMMIT = "e0f985bb55ad03d8b3435bc6a334a58d3ed674ba"  # FRBs/zdm HEAD, 2026-10-10
CACHE = ROOT / "data" / "dmzcal" / f"frb_{FRB_COMMIT[:7]}"
RAW = f"https://raw.githubusercontent.com/FRBs/FRB/{FRB_COMMIT}/frb/data"
TREE = f"https://api.github.com/repos/FRBs/FRB/git/trees/{FRB_COMMIT}?recursive=1"

# --- Parameter state under audit -------------------------------------------------
# zdm states.py "HoffmannEmin25" = Hoffmann et al. 2025, PASA 42, e017
# (doi:10.1017/pasa.2024.127; arXiv:2408.04878v2). zdm's load_state() default is
# "HoffmannHalo25", which states.py marks "Unpublished (yet!)": its fit sample cannot
# be established, so it is not audited.
STATE = "HoffmannEmin25"

# Fit sample of HoffmannEmin25, from the arXiv:2408.04878v2 source (main.tex):
# - Sec. 2 / Table "DSA FRBs used in this analysis" (main.tex l.242-280): 25 DSA-110
#   bursts; DMs of all enter the likelihood, z of only the three bolded ones.
DSA_FIT = (
    "20220121B 20220204A 20220207C 20220208A 20220307B 20220310F 20220319D 20220330D "
    "20220418A 20220424E 20220506D 20220509G 20220726A 20220801A 20220825A 20220831A "
    "20220914A 20220920A 20220926A 20221002A 20221012A 20221027A 20221029A 20221101B "
    "20221101A"
).split()
DSA_FIT_Z_USED = {"20220207C", "20220319D", "20220509G"}
# - Sec. 2.4 "Updated CRAFT surveys" (main.tex l.318-319): the James+2022b / Baptista+2023 Parkes, ASKAP Fly's-Eye
#   and CRAFT/ICS samples, plus Table "Additional CRAFT FRBs" (l.281-317), "all FRBs
#   detected by CRAFT up until the end of 2023". Rows commented out of that table may
#   or may not have entered via the earlier samples; rather than guess, every ASKAP
#   burst dated <= 2023-12-31 is put on the production side (conservative: it can only
#   shrink the certification sample).
ASKAP_FIT_CUTOFF = "20231231"
# - Sec. 2.1 (main.tex l.142-153): CHIME, MeerKAT and UTMOST were explicitly excluded.
# - F is fixed to 0.32 citing Macquart+2020 and Zhang+2021 (Sec. 4, main.tex l.382),
#   not fitted, so F carries no burst-level provenance in this state.

# zdm survey model used for E1, by discovery instrument (zdm/data/Surveys at ZDM_COMMIT).
SURVEY_MODEL = {
    "ASKAP": "CRAFT_average_ICS (2024 certification bursts are Shannon+2024 ICS detections)",
    "DSA": "DSA",
    "CHIME": "CHIME (decbin files; HoffmannEmin25 did not fit CHIME - see Sec. 2.1)",
    "MeerKAT": "MeerTRAPcoherent unless SURVEY_OVERRIDE gives the beam mode",
}

# Spectroscopic-z rule (plan 99 rule 2), refined 2026-10-10 before any posterior:
# public_hosts.csv leaves `Spectrum` blank for the CHIME/KKO hosts of Leung et al. 2025
# (arXiv:2502.11217v2), whose follow-up table (redshifts_table.tex, the table cited at
# main.tex l.561) gives a spectroscopic z (Lick/Keck/Gemini or an archival spectrum)
# for each of these gold-sample hosts (P(O|x) > 0.9, main.tex l.345). Their other
# localizations ("remaining", sample_full.tex) have no follow-up and P(O|x) <= 0.75,
# so they stay out. z as printed there, used to cross-check the FRBs/FRB value:
LEUNG25_SPEC_Z = {
    "20230203A": 0.1464,
    "20230222A": 0.1223,
    "20230222B": 0.1100,
    "20230311A": 0.1918,
    "20230703A": 0.1184,
    "20230730A": 0.2115,
    "20230926A": 0.0553,
    "20231005A": 0.0713,
    "20231011A": 0.0783,
    "20231017A": 0.2450,
    "20231025B": 0.3238,
    "20231123A": 0.0729,
    "20231128A": 0.1079,
    "20231201A": 0.1119,
    "20231204A": 0.0644,
    "20231206A": 0.0659,
    "20231223C": 0.1059,
    "20231229A": 0.0190,
    "20231230A": 0.0298,
}
# Leung's 20230311A has "a secure redshift, but no secure host" (main.tex l.553): kept
# for spec-z, flagged by its P(O|x) below like any other burst.
P_OX_MIN = 0.9  # host-association floor where public_hosts.csv reports P_Ox

# Primary-source corrections to FRBs/FRB, applied here so every estimator sees one sample.
# Each carries its locator. Applied AFTER the first real run (step-0 finding 7 and the
# GATE-2 round-1 verification); survey/dmzcal-findings.md reports results before and after.
# name -> (value, locator)
Z_OVERRIDE: dict[str, tuple[float, str]] = {
    "FRB20231201A": (0.1119, "Leung+2025 arXiv:2502.11217v2 redshifts_table.tex"),
    "FRB20201124A": (
        0.0979,
        "Fong+2021 arXiv:2106.11993 abstract (MMT, 0.0979 +/- 0.0001); "
        "FRBs/FRB 0.0982 has no located source",
    ),
}
TELESCOPE_OVERRIDE: dict[str, tuple[str, str]] = {
    "FRB20201124A": (
        "CHIME",
        "Lanman+2022 arXiv:2109.09254 abstract: 'first discovered by "
        "CHIME/FRB'; MeerKAT was follow-up only",
    ),
}
DM_OVERRIDE: dict[str, tuple[float, str]] = {
    "FRB20210410D": (578.78, "Caleb+2023 arXiv:2302.09754 burst-properties table (+/- 2)"),
}
# Unresolved conflicts between primary sources: NOT applied to the headline; scored as a
# post-hoc sensitivity row (GATE-2 round 2, N4). name -> {"z"|"DM": (value, locator)}
ALTERNATIVES: dict[str, dict[str, tuple[float, str]]] = {
    "FRB20231120A": {
        "z": (
            0.0700,
            "Connor+2024 arXiv:2409.16952 final table (vs Sharma Ext. Data Table 1 0.0368)",
        )
    },
    "FRB20230124A": {"DM": (590.6, "Connor+2024 arXiv:2409.16952 final table")},
    "FRB20230307A": {"DM": (608.9, "Connor+2024 arXiv:2409.16952 final table")},
    "FRB20230501A": {"DM": (532.5, "Connor+2024 arXiv:2409.16952 final table")},
}

# zdm survey model per burst where the instrument mode is known (default: by telescope)
SURVEY_OVERRIDE: dict[str, tuple[str, str]] = {
    "FRB20210410D": (
        "MeerTRAPincoherent",
        "Caleb+2023 arXiv:2302.09754 table: 'Beam: Incoherent beam'",
    ),
}

# FRB 20240304B, added by hand: Caleb et al. 2026 (arXiv:2508.01648), Table 1 (DM,
# scattering-corrected; l, b; NE2001 DM_ISM) and Table 2 (z_spec). Chosen after its z
# was known -> excluded from every statistic by plan 99 sample rule 4.
FRB20240304B = {
    "name": "FRB20240304B",
    "DM": 2458.20,
    "DMISM": 28.1,
    "gl": 269.8676,
    "gb": 72.1035,
    "z": 2.148,
    "telescope": "MeerKAT",
    "spec_z": True,
    "refs": ["Caleb2026"],
}


def _get(url: str, dest: Path) -> Path:  # pragma: no cover - network
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(url, headers={"User-Agent": "jansky-research"})
        with urllib.request.urlopen(req, timeout=60) as r:
            dest.write_bytes(r.read())
    return dest


def fetch() -> tuple[list[dict], list[dict], list[dict]]:  # pragma: no cover - network
    base = _get(f"{RAW}/FRBs/FRBs_base.csv", CACHE / "FRBs_base.csv")
    hosts = _get(f"{RAW}/Galaxies/public_hosts.csv", CACHE / "public_hosts.csv")
    tree = json.loads(_get(TREE, CACHE / "tree.json").read_text())
    jsons = []
    for node in tree["tree"]:
        m = re.fullmatch(r"frb/data/FRBs/(FRB2\d+[A-Za-z]*)\.json", node["path"])
        if m:
            p = _get(f"{RAW}/FRBs/{m.group(1)}.json", CACHE / "json" / f"{m.group(1)}.json")
            jsons.append(json.loads(p.read_text()) | {"_file": m.group(1)})
    return (
        list(csv.DictReader(base.open())),
        list(csv.DictReader(hosts.open())),
        jsons,
    )


def side(telescope: str, name: str) -> tuple[str, str]:
    return provenance_side(telescope, name, DSA_FIT, DSA_FIT_Z_USED, ASKAP_FIT_CUTOFF)


def assemble(base: list[dict], hosts: list[dict], jsons: list[dict]) -> dict:
    tel = {_key(r["Name"]): r["telescope"] for r in base}
    tel_by_date: dict[str, set[str]] = {}
    for r in base:
        tel_by_date.setdefault(_date(r["Name"]), set()).add(r["telescope"])
    spec = {_key(r["FRB"]): bool(r.get("Spectrum", "").strip()) for r in hosts}
    pox: dict[str, float | None] = {}
    for r in hosts:
        try:
            pox[_key(r["FRB"])] = float(r["P_Ox"])
        except (TypeError, ValueError):
            pox[_key(r["FRB"])] = None
    in_hosts = set(spec)
    with_z = [j for j in jsons if j.get("z") not in (None, "")]
    kept, merges = dedupe(with_z)
    rows = []
    for j in kept:
        k = _key(j["_file"])
        t = tel.get(k) or tel.get(re.sub(r"[A-Za-z]+$", "", k), "")
        if not t and len(tel_by_date.get(_date(k), set())) == 1:
            t = next(iter(tel_by_date[_date(k)]))  # letter suffix differs between files
        t_raw = t
        if j["_file"] in TELESCOPE_OVERRIDE:
            t = TELESCOPE_OVERRIDE[j["_file"]][0]
        s, why = side(t, j["_file"])
        is_spec = spec.get(k, False) or k in LEUNG25_SPEC_Z
        p = pox.get(k)
        z_lit = LEUNG25_SPEC_Z.get(k)
        rows.append(
            {
                "name": j["_file"],
                "telescope": t,
                "z": Z_OVERRIDE[j["_file"]][0] if j["_file"] in Z_OVERRIDE else float(j["z"]),
                "z_frbs_frb": float(j["z"]),
                "telescope_frbs_frb": t_raw,
                "DM": DM_OVERRIDE[j["_file"]][0] if j["_file"] in DM_OVERRIDE else _dm(j),
                "DM_frbs_frb": _dm(j),
                "survey_override": SURVEY_OVERRIDE[j["_file"]][0]
                if j["_file"] in SURVEY_OVERRIDE
                else None,
                "DMISM": j.get("DMISM", {}).get("value")
                if isinstance(j.get("DMISM"), dict)
                else j.get("DMISM"),
                "repeater": bool(j.get("repeater")),
                "spec_z": is_spec,
                "in_public_hosts": k in in_hosts,
                "P_Ox": p,
                "secure_host": p is None or p >= P_OX_MIN,
                "z_literature_check": None
                if z_lit is None
                else {"Leung2025": z_lit, "agrees": abs(z_lit - float(j["z"])) < 5e-4},
                "refs": j.get("refs", []),
                "side": s,
                "reason": why,
            }
        )
    rows.sort(key=lambda r: r["name"])
    cert = [r for r in rows if r["side"] == "certification"]
    cert_spec = [r for r in cert if r["spec_z"] and r["secure_host"]]
    return {
        "source": "real: FRBs/FRB localized-host JSONs + HoffmannEmin25 fit sample "
        "(arXiv:2408.04878v2 source)",
        "plan": "plans/99-dmzcal-coverage.md step 0",
        "frb_commit": FRB_COMMIT,
        "zdm_commit": ZDM_COMMIT,
        "state_audited": STATE,
        "state_reference": "Hoffmann et al. 2025, PASA 42, e017, doi:10.1017/pasa.2024.127",
        "F_provenance": "F = 0.32 fixed, citing Macquart+2020 and Zhang+2021 "
        "(arXiv:2408.04878v2 Sec. 4, main.tex l.382); not fitted",
        "zdm_default_state_note": "zdm load_state() defaults to HoffmannHalo25, marked "
        "'Unpublished (yet!)' in states.py; not audited",
        "survey_models": SURVEY_MODEL,
        "counts": {
            "json_files": len(jsons),
            "with_z": len(with_z),
            "after_dedupe": len(rows),
            "by_side": dict(Counter(r["side"] for r in rows)),
            "certification_by_telescope": dict(Counter(r["telescope"] for r in cert)),
            "certification_spec_z": len(cert_spec),
            "certification_spec_z_by_telescope": dict(Counter(r["telescope"] for r in cert_spec)),
            "certification_spec_z_repeaters": sum(r["repeater"] for r in cert_spec),
            "certification_spec_z_but_insecure_host": sum(
                r["spec_z"] and not r["secure_host"] for r in cert
            ),
        },
        "corrections": {
            "z": {k: {"value": v, "locator": loc} for k, (v, loc) in Z_OVERRIDE.items()},
            "telescope": {
                k: {"value": v, "locator": loc} for k, (v, loc) in TELESCOPE_OVERRIDE.items()
            },
            "DM": {k: {"value": v, "locator": loc} for k, (v, loc) in DM_OVERRIDE.items()},
            "survey_model": {
                k: {"value": v, "locator": loc} for k, (v, loc) in SURVEY_OVERRIDE.items()
            },
            "alternatives": {
                k: {f: {"value": v, "locator": loc} for f, (v, loc) in d.items()}
                for k, d in ALTERNATIVES.items()
            },
            "unresolved_not_applied": [
                "DSA DMs 20230124A/20230307A/20230501A differ by 0.6-1.25 pc/cc from "
                "Connor+2024 arXiv:2409.16952 final table (Sharma+2024 has no DM column)",
                "Connor+2024 gives z = 0.0700 for 20231120A vs Sharma+2024 Ext. Data "
                "Table 1 0.0368 (= FRBs/FRB); Sharma kept",
            ],
        },
        "z_discrepancies": [
            {"name": r["name"], "FRBs_FRB": r["z"], **r["z_literature_check"]}
            for r in rows
            if r["z_literature_check"] and not r["z_literature_check"]["agrees"]
        ],
        "underpowered_by_rule": len(cert_spec) < 25,
        "merges": merges,
        "recover_a_known": FRB20240304B,
        "bursts": rows,
    }


def main() -> None:  # pragma: no cover - network
    out = assemble(*fetch())
    path = ROOT / "results" / "dmzcal_provenance.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["counts"], indent=1))
    print("merges:", len(out["merges"]), "underpowered:", out["underpowered_by_rule"])


if __name__ == "__main__":
    main()
