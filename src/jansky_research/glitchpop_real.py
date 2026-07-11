"""JBO glitch waiting-time census -- real leg (plan 48, F11). Scrapes the live catalogue; not in CI.

Scrapes the live Jodrell Bank glitch table, classifies every pulsar with >=5 glitches by its
inter-glitch waiting-time distribution, checks the known quasi-periodic glitchers come out
quasi-periodic (the real-data recover-a-known), and diffs the classifications against the end-2018
Basu+2022 subset. The classifier + its synthetic recover-a-known run offline in core CI.
"""

from __future__ import annotations

from pathlib import Path

from .glitchpop import (
    JBO_URL,
    KNOWN_QUASIPERIODIC,
    MIN_GLITCHES,
    classification_delta,
    group_by_pulsar,
    parse_glitch_table,
    population_census,
    population_significance,
)

# Magnetars / AXPs whose catalogued "glitches" are X-ray-outburst-driven, not rotation-powered radio
# glitches --- excluded from the population (their waiting times are even more monitoring-dominated).
MAGNETARS = ("1E_2259+586", "1RXS_J1708-4009", "1E_1841-045", "CXOU_J1714", "SGR", "PSR_J1846-0258")


def scrape_glitch_table(url: str = JBO_URL) -> str:  # pragma: no cover - network
    """Fetch the live JBO glitch table HTML (no auth; a browser UA avoids the odd server block)."""
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (jansky-research)"})
    with urllib.request.urlopen(req, timeout=60) as r:  # noqa: S310 (trusted JBO host)
        return r.read().decode("utf-8", errors="replace")


def run_real_census(out: str, *, min_glitches: int = MIN_GLITCHES) -> dict:  # pragma: no cover
    """Full real census: scrape -> classify every pulsar (>=5 glitches) -> known-QP check -> delta."""
    import json

    glitches = parse_glitch_table(scrape_glitch_table())
    by = group_by_pulsar(glitches)
    # drop magnetars/AXPs (their catalogued glitches are X-ray-outburst-driven, not rotation-powered)
    n_magnetars = sum(1 for j in by if any(m in j for m in MAGNETARS))
    by = {j: d for j, d in by.items() if not any(m in j for m in MAGNETARS)}
    rows = population_census(by, min_glitches=min_glitches)
    n_exp = sum(r["klass"] == "exponential" for r in rows)
    n_qp = sum(r["klass"] == "quasi_periodic" for r in rows)
    n_cl = sum(r["klass"] == "clustered" for r in rows)
    # real-data recover-a-known: the known quasi-periodic glitchers must NOT come out exponential
    known = {r["jname"]: r["klass"] for r in rows if r["jname"] in KNOWN_QUASIPERIODIC}
    known_ok = (
        "yes"
        if known and all(k == "quasi_periodic" for k in known.values())
        else ("no" if known else "absent")
    )
    delta = classification_delta(by, min_glitches=min_glitches)
    sigstats = population_significance(rows)
    metrics = {
        "source": "JBO glitch catalogue (jb.man.ac.uk); per-pulsar waiting-time classification + post-2018 delta",
        "is_real": True,
        "n_glitches": len(glitches),
        "n_pulsars": len(by),
        "n_magnetars_dropped": int(n_magnetars),
        "n_qualified_full": len(rows),
        "n_exponential": int(n_exp),
        "n_quasiperiodic": int(n_qp),
        "n_clustered": int(n_cl),
        **sigstats,
        "known_quasiperiodic": known,
        "known_quasiperiodic_ok": known_ok,
        "n_newly_classifiable": delta["n_newly_classifiable"],
        "n_stable_sample": delta["n_stable_sample"],
        "n_flipped": delta["n_flipped"],
        "flipped": delta["flipped"],
        "newly_classifiable": delta["newly_classifiable"][:40],
        "census": rows,
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "glitchpop_census.json").write_text(json.dumps(metrics, indent=2) + "\n")
    # drop the bulky per-pulsar table from the returned summary
    return {k: v for k, v in metrics.items() if k != "census"}
