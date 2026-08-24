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


def run_real_census(
    out: str,
    *,
    min_glitches: int = MIN_GLITCHES,
    from_csv: str | None = None,
    retrieved: str = "unrecorded",
) -> dict:  # pragma: no cover
    """Full real census: snapshot (or scrape) -> classify -> known-QP check -> delta.

    ``from_csv`` analyses a committed snapshot instead of the live table -- the live catalogue
    demonstrably mutates (three recorded scrapes returned 222, 223 and 220 pulsars; a fourth,
    2026-08-24, returned 224), so without a pinned snapshot no number in the paper is
    regenerable, by anyone, ever. ``retrieved`` is the snapshot's retrieval date, committed with
    the results and stated in the paper.
    """
    import csv as _csv

    if from_csv:
        with open(from_csv, newline="") as fh:
            glitches = [
                {
                    "jname": r["jname"],
                    "mjd": float(r["mjd"]),
                    "size": float(r["size"]) if r["size"] else float("nan"),
                    "dnudot": float(r["dnudot"]) if r["dnudot"] else float("nan"),
                    "refs": r["refs"],
                    "is_new": r["is_new"] == "True",
                }
                for r in _csv.DictReader(fh)
            ]
    else:
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
    from .glitchpop import (
        BASU_END_MJD,
        census_accounting,
        gap_factor_sweep,
        injection_surface,
        p_uniformity_check,
    )
    from .glitchpop import (
        population_significance as _popsig,
    )

    acct = census_accounting(by, min_glitches=min_glitches)
    sweep = gap_factor_sweep(by)
    ks = p_uniformity_check(rows)
    surf = injection_surface()
    sig_emp = _popsig(rows, fp_rate_by_n=surf["fp_rate_by_n"])
    # the pre-split glitch count MEASURED from the analysed table, so "the post-2018 increment"
    # is an epoch statement about this snapshot, not a difference of catalogue vintages
    n_pre_split = sum(1 for g in glitches if float(g["mjd"]) < BASU_END_MJD)
    metrics = {
        "source": (
            "JBO glitch catalogue (jb.man.ac.uk), snapshot retrieved "
            f"{retrieved}; per-pulsar waiting-time classification + post-2018 delta"
        ),
        "is_real": True,
        "retrieved": retrieved,
        # matched pairs: the abstract previously counted glitches on the raw catalogue and
        # pulsars on the magnetar-filtered sample in one sentence
        "n_glitches_raw": len(glitches),
        "n_pulsars_raw": len({g["jname"] for g in glitches}),
        "n_glitches_analysed": int(sum(len(d["mjd"]) for d in by.values())),
        "n_pulsars_analysed": len(by),
        "n_glitches_pre_split": int(n_pre_split),
        "n_glitches_post_split": int(len(glitches) - n_pre_split),
        "n_retroactive_pre_split": int(
            sum(1 for g in glitches if float(g["mjd"]) < BASU_END_MJD and g.get("is_new"))
        ),
        **acct,
        "gap_factor_sweep": sweep,
        "p_uniformity": ks,
        "injection_surface": surf["surface"],
        "exponential_outcome_rates": surf["exponential_rows"],
        "qp_poisson_binomial_p": sig_emp.get("qp_poisson_binomial_p"),
        "expected_false_qp_empirical": sig_emp.get("expected_false_qp_empirical"),
        "expected_false_clustered": sig_emp.get("expected_false_clustered"),
        "clustered_binomial_p": sig_emp.get("clustered_binomial_p"),
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
    from .report import write_results

    write_results(metrics, op / "results" / "glitchpop_census.json")
    # drop the bulky per-pulsar table from the returned summary
    return {k: v for k, v in metrics.items() if k != "census"}
