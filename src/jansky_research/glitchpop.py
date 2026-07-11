"""JBO glitch catalogue: a per-pulsar waiting-time classification + the post-2018 delta (plan 48, F11).

The Jodrell Bank glitch catalogue (jb.man.ac.uk) now holds 727 glitches from 223 pulsars; the last
population-statistics paper, Basu+2022 (MNRAS 510, 4049), stops at end-2018 (543/178), so ~184
glitches have never entered a population analysis. Basu+2022 modelled the aggregate glitch-SIZE and
spin-down-change distributions (not the per-pulsar waiting-time behaviour). Howitt+2018 (ApJ 867, 60)
DID classify per-pulsar waiting-time distributions (Poisson vs quasi-periodic) --- but only for the
handful of most-glitching pulsars, on a pre-2018 sample, and without gap-robustness. Zhu & Zheng
(2025, arXiv:2501.01862) use the same >=5-glitch sample but for glitch-CLUSTER-period-vs-age.

**The wedge this slice adds** (over Howitt+2018): the per-pulsar waiting-time classification extended
to the FULL >=5-glitch sample, made robust to MONITORING GAPS, and diffed across the post-2018 data.
For each pulsar we classify its inter-glitch waiting-time distribution --- **exponential** (memoryless
Poisson, the null), **quasi-periodic** (regular, Vela/J0537-like), or **clustered** (over-dispersed)
--- from the coefficient of variation, calibrated against the exponential null by a parametric
bootstrap. The headline is the **post-2018 delta**: which pulsars newly cross the >=5-glitch threshold
and which classifications FLIP when the post-2018 glitches are added to the end-2018 Basu subset.

**Monitoring-gap excision (the load-bearing methodological point).** The archival waiting times
conflate real inter-glitch intervals with intervals in which a pulsar kept glitching but was not
observed. J0537-6910's 2264-day RXTE->NICER gap among ~100-day intervals lifts its raw CV to ~2, so a
naive fit (and an initial gamma-BIC version of this classifier) calls the textbook quasi-periodic
glitcher EXPONENTIAL. Excising waits > 6x the median recovers it (and Vela); the CV-bootstrap replaced
the gamma-BIC after the latter also underweighted Vela's milder regularity --- so the two-object
real-data recover-a-known is partly in-sample (the statistic was reselected against it), and the
independent validation is the SYNTHETIC injection test.

**Honest scope.** (1) Glitch detection is incomplete/uneven, so regularity is a lower bound;
classification changes are data-driven, not physical claims. (2) The >=5-glitch cut selects active
glitchers. (3) Each pulsar is tested at 0.05, so ~5% false quasi-periodics are expected --- the
aggregate quasi-periodic fraction is significant (`population_significance`) but individual borderline
members are not secure. (4) Gap excision assumes large waits are monitoring gaps, so it CANNOT detect
genuine long-gap clustering (which Zhu&Zheng do report on this sample); "0-few clustered" is a method
asymmetry, not evidence against clustering. (5) Magnetars are excluded (X-ray-outburst-driven).
Reuse: the `lpt`/`ppdot` scrape discipline.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

__all__ = [
    "parse_glitch_table",
    "group_by_pulsar",
    "waiting_times",
    "waiting_time_fit",
    "size_class",
    "classify_pulsar",
    "population_census",
    "population_significance",
    "classification_delta",
    "synthetic_glitch_series",
    "inject_recover",
    "run",
]

JBO_URL = "https://www.jb.man.ac.uk/pulsar/glitches/gTable.html"
BASU_END_MJD = 58484.0  # 2019-01-01; the end-2018 cut of the last population paper (Basu+2022)
MIN_GLITCHES = (
    5  # per-pulsar minimum glitch count for a waiting-time classification (pre-registered)
)
# Known quasi-periodic glitchers, in the catalogue's own naming (it mixes J- and B-names): J0537-6910
# and Vela (listed as B0833-45). These must come out quasi-periodic (the real-data recover-a-known).
KNOWN_QUASIPERIODIC = ("J0537-6910", "B0833-45")


def parse_glitch_table(html: str) -> list[dict]:
    """Parse the JBO ``gTable.html`` into per-glitch records (JNAME carried forward, 'X' -> NaN).

    Each data row is 11 cells: index, JNAME, B-name, glitch#, MJD, MJD-err, dnu/nu (x1e-9), err,
    dnudot/nudot (x1e-3), err, references. The JNAME cell is blank on some repeat-glitch rows, so we
    carry the last non-empty name forward. ``is_new`` flags glitches marked "New" in the references
    (provenance). Rows without a numeric MJD are skipped.
    """
    rows = re.findall(r"<tr[^>]*>(.*?)</tr>", html, re.S | re.I)
    out: list[dict] = []
    last_name = ""
    for r in rows:
        cells = [
            re.sub(r"<[^>]*>", "", c).strip()
            for c in re.findall(r"<td[^>]*>(.*?)</td>", r, re.S | re.I)
        ]
        if len(cells) < 11:
            continue
        try:
            mjd = float(cells[4])
        except ValueError:
            continue
        name = cells[1] or last_name
        if not name:
            continue
        last_name = name
        out.append(
            {
                "jname": name,
                "mjd": mjd,
                "size": _to_float(cells[6]),  # dnu/nu in units of 1e-9
                "dnudot": _to_float(cells[8]),  # dnudot/nudot in units of 1e-3
                "refs": cells[10],
                "is_new": "new" in cells[10].lower(),
            }
        )
    return out


def _to_float(s: str) -> float:
    try:
        return float(s)
    except ValueError:
        return float("nan")


def group_by_pulsar(glitches: list[dict]) -> dict:
    """Group glitch records by J-name into ``{jname: {"mjd": array, "size": array, "n": int}}``."""
    by: dict = {}
    for g in glitches:
        d = by.setdefault(g["jname"], {"mjd": [], "size": []})
        d["mjd"].append(g["mjd"])
        d["size"].append(g["size"])
    for d in by.values():
        order = np.argsort(d["mjd"])
        d["mjd"] = np.asarray(d["mjd"], float)[order]
        d["size"] = np.asarray(d["size"], float)[order]
        d["n"] = int(d["mjd"].size)
    return by


def waiting_times(mjds: np.ndarray) -> np.ndarray:
    """Sorted inter-glitch waiting times in YEARS (positive successive-MJD differences)."""
    t = np.sort(np.asarray(mjds, float))
    w = np.diff(t)
    return w[w > 0] / 365.25


def waiting_time_fit(
    waits: np.ndarray,
    *,
    min_waits: int = 4,
    sig: float = 0.05,
    gap_factor: float = 6.0,
    n_boot: int = 4000,
    seed: int = 0,
) -> dict:
    """Classify a waiting-time distribution: exponential (Poisson) vs quasi-periodic vs clustered.

    **Monitoring-gap excision (the load-bearing step).** The archival catalogue conflates real
    inter-glitch intervals with MONITORING gaps --- intervals in which the pulsar kept glitching but
    was not observed. A single such gap wrecks the raw statistics: the textbook quasi-periodic
    glitcher J0537-6910 has a 2264-day gap (the RXTE->NICER hiatus) among ~100-day intervals, which
    lifts its raw coefficient of variation to ~2 so a naive fit calls it exponential. We drop waits
    longer than ``gap_factor`` x the median as monitoring gaps before classifying: for a regular
    glitcher a gap that large is almost certainly unobserved, while a genuine exponential loses only
    its far tail (P(w > 6*median) ~ 1.5%), so the excision barely biases the null.

    **The test.** The coefficient of variation ``cv = std/mean`` is the regularity measure: ~1 for an
    exponential (memoryless Poisson) process, <1 for regular/quasi-periodic, >1 for clustered. We
    calibrate it against the exponential null with a parametric bootstrap --- the CV distribution of
    ``n`` exponential waits (CV is scale-free) --- and classify quasi-periodic if the observed CV is
    significantly LOW (``p<sig``), clustered if significantly HIGH, else exponential. This is more
    powerful than a gamma-BIC (which underweights the mild-but-real regularity of, e.g., Vela).
    ``insufficient`` below ``min_waits``.
    """
    w = np.asarray(waits, float)
    w = w[np.isfinite(w) & (w > 0)]
    n_raw = int(w.size)
    if n_raw < min_waits:
        return {"n_waits": n_raw, "klass": "insufficient", "cv": float("nan")}
    med = float(np.median(w))
    keep = w <= gap_factor * med
    n_gaps = int((~keep).sum())
    w = w[keep]
    n = int(w.size)
    if n < min_waits:
        return {"n_waits": n, "n_gaps_excised": n_gaps, "klass": "insufficient", "cv": float("nan")}
    mean = float(w.mean())
    cv = float(w.std(ddof=0) / mean) if mean > 0 else float("nan")
    # parametric bootstrap null: draw n_raw exponential waits and apply the SAME gap excision before
    # taking the CV, so the null undergoes identical processing to the data (upper-truncation lowers
    # CV, so replicating it keeps the test from being biased toward "quasi-periodic"). Scale-free.
    rng = np.random.default_rng(seed)
    samp = rng.exponential(1.0, size=(n_boot, n_raw))
    med_b = np.median(samp, axis=1, keepdims=True)
    kept = np.where(samp <= gap_factor * med_b, samp, np.nan)
    with np.errstate(invalid="ignore"):
        boot_cv = np.nanstd(kept, axis=1, ddof=0) / np.nanmean(kept, axis=1)
    boot_cv = boot_cv[np.isfinite(boot_cv)]
    p_low = float(np.mean(boot_cv <= cv))  # P(CV this regular | exponential)
    p_high = float(np.mean(boot_cv >= cv))  # P(CV this clustered | exponential)
    if p_low < sig:
        klass = "quasi_periodic"
    elif p_high < sig:
        klass = "clustered"
    else:
        klass = "exponential"
    return {
        "n_waits": n,
        "n_waits_raw": n_raw,
        "n_gaps_excised": n_gaps,
        "cv": round(cv, 3),
        "p_regular": round(p_low, 4),
        "p_clustered": round(p_high, 4),
        "klass": klass,
    }


def size_class(sizes: np.ndarray) -> dict:
    """Secondary: glitch-size (dnu/nu) summary --- power-law tail index + log-normal (mu, sigma)."""
    from .frbstats import fit_power_law

    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s > 0)]
    if s.size < 3:
        return {"n_sizes": int(s.size), "gamma": float("nan"), "log_mu": float("nan")}
    ls = np.log10(s)
    out = {
        "n_sizes": int(s.size),
        "log_mu": round(float(ls.mean()), 3),
        "log_sigma": round(float(ls.std(ddof=0)), 3),
    }
    try:
        pl = fit_power_law(s, auto_xmin=True)
        out["gamma"] = round(float(pl.gamma), 3)
    except Exception:
        out["gamma"] = float("nan")
    return out


def classify_pulsar(d: dict, *, min_glitches: int = MIN_GLITCHES) -> dict:
    """Full per-pulsar classification (waiting-time class + size summary) for a grouped record."""
    mjd = np.asarray(d["mjd"], float)
    n = int(mjd.size)
    if n < min_glitches:
        return {"n": n, "klass": "insufficient"}
    wt = waiting_time_fit(waiting_times(mjd))
    sz = size_class(np.asarray(d.get("size", []), float))
    return {"n": n, "span_yr": round(float((mjd.max() - mjd.min()) / 365.25), 2), **wt, **sz}


def population_census(by_pulsar: dict, *, min_glitches: int = MIN_GLITCHES) -> list[dict]:
    """Classify every pulsar with >= ``min_glitches`` glitches; return rows sorted by glitch count."""
    rows = []
    for jname, d in by_pulsar.items():
        c = classify_pulsar(d, min_glitches=min_glitches)
        if c["klass"] == "insufficient":
            continue
        rows.append({"jname": jname, **c})
    rows.sort(key=lambda r: r["n"], reverse=True)
    return rows


def population_significance(rows: list[dict], *, sig: float = 0.05) -> dict:
    """Is the quasi-periodic FRACTION significant against the all-exponential null (multiple testing)?

    Each pulsar is tested one-sided at ``sig``, so under an all-Poisson population ~``n_fit*sig`` false
    quasi-periodics are expected by chance. This returns the expected false count and the binomial
    probability of seeing at least ``n_quasiperiodic`` --- so the paper can state honestly that the
    *aggregate* excess of quasi-periodic pulsars is real while individual borderline members are not
    secure (~1--2 are expected spurious).
    """
    from scipy.stats import binomtest

    n_fit = len(rows)
    n_qp = int(sum(r.get("klass") == "quasi_periodic" for r in rows))
    exp_false = round(n_fit * sig, 2)
    binom_p = (
        float(binomtest(n_qp, n_fit, sig, alternative="greater").pvalue) if n_fit else float("nan")
    )
    return {
        "expected_false_qp": exp_false,
        "qp_binomial_p": binom_p,
        "qp_excess_significant": bool(np.isfinite(binom_p) and binom_p < 0.05),
    }


def classification_delta(
    by_pulsar: dict, *, split_mjd: float = BASU_END_MJD, min_glitches: int = MIN_GLITCHES
) -> dict:
    """Which pulsars newly qualify, and which classifications FLIP, when post-split glitches are added.

    Reclassifies each pulsar on (a) the pre-``split_mjd`` subset (the end-2018 Basu epoch) and (b) the
    full catalogue. ``newly_qualified`` = crossed the >=``min_glitches`` threshold only with the new
    glitches; ``flipped`` = had >=``min_glitches`` in both epochs but the waiting-time class changed.
    This is the slice's headline --- the impact of the post-2018 data on the population picture.
    """
    newly, flips = [], []
    n_both = 0
    for jname, d in by_pulsar.items():
        mjd = np.asarray(d["mjd"], float)
        size = np.asarray(d["size"], float)
        pre = mjd < split_mjd
        full_c = classify_pulsar(d, min_glitches=min_glitches)
        if full_c["klass"] == "insufficient":
            continue
        pre_c = (
            classify_pulsar({"mjd": mjd[pre], "size": size[pre]}, min_glitches=min_glitches)
            if pre.sum() >= min_glitches
            else {"klass": "insufficient"}
        )
        if pre_c["klass"] == "insufficient":
            # no real class from the end-2018 subset -> newly classifiable, NOT a flip
            newly.append(
                {
                    "jname": jname,
                    "n_pre": int(pre.sum()),
                    "n_now": int(mjd.size),
                    "klass": full_c["klass"],
                }
            )
            continue
        n_both += 1
        if pre_c["klass"] != full_c["klass"]:
            flips.append(
                {
                    "jname": jname,
                    "was": pre_c["klass"],
                    "now": full_c["klass"],
                    "n_pre": int(pre.sum()),
                    "n_now": int(mjd.size),
                }
            )
    return {
        "n_qualified_full": len(newly) + n_both,
        "n_newly_classifiable": len(newly),
        "n_stable_sample": n_both,
        "n_flipped": len(flips),
        "flipped": flips,
        "newly_classifiable": newly,
    }


def synthetic_glitch_series(
    *,
    kind: str = "exponential",
    n: int = 12,
    mean_yr: float = 2.0,
    start_mjd: float = 50000.0,
    seed: int = 0,
) -> dict:
    """Synthetic glitch MJDs from a known waiting-time process, for the recover-a-known.

    ``exponential`` = memoryless Poisson (cv~1); ``quasi_periodic`` = regular intervals with small
    jitter (cv<<1, Vela/J0537-like); ``clustered`` = over-dispersed bursts (cv>1). Returns MJDs and a
    plausible size array so `classify_pulsar` can be run end-to-end.
    """
    rng = np.random.default_rng(seed)
    if kind == "exponential":
        waits = rng.exponential(mean_yr, n - 1)
    elif kind == "quasi_periodic":
        waits = np.abs(rng.normal(mean_yr, 0.12 * mean_yr, n - 1))
    elif kind == "clustered":
        # bursts of short waits separated by long gaps -> over-dispersed
        short = rng.exponential(0.15 * mean_yr, n - 1)
        gap = rng.random(n - 1) < 0.3
        waits = np.where(gap, rng.exponential(4.0 * mean_yr, n - 1), short)
    else:
        raise ValueError(f"unknown kind {kind!r}")
    mjd = start_mjd + np.concatenate([[0.0], np.cumsum(waits) * 365.25])
    sizes = np.power(10.0, rng.normal(1.5, 0.6, n))  # dnu/nu ~ lognormal, arbitrary units
    return {"mjd": mjd, "size": sizes, "kind": kind}


def inject_recover(*, counts=(6, 10, 20, 40), n_each: int = 40, seed: int = 0) -> dict:
    """Recover-a-known: classification accuracy for exponential vs quasi-periodic, vs glitch count.

    Builds ``n_each`` synthetic pulsars of each kind at each glitch count and records how often the
    waiting-time classifier labels them correctly. Accuracy rises with count (few glitches genuinely
    cannot distinguish a regular from a Poisson process) --- the honest completeness floor. Also
    returns the exponential false-positive rate (Poisson series misclassified as regular/clustered).
    """
    rng = np.random.default_rng(seed)
    exp_acc, qp_acc = {}, {}
    fp = 0
    fp_total = 0
    for c in counts:
        e_ok = q_ok = 0
        for _ in range(n_each):
            se = int(rng.integers(1 << 30))
            e = classify_pulsar(synthetic_glitch_series(kind="exponential", n=c, seed=se))
            fp_total += 1
            if e["klass"] == "exponential":
                e_ok += 1
            else:
                fp += 1
            sq = int(rng.integers(1 << 30))
            q = classify_pulsar(synthetic_glitch_series(kind="quasi_periodic", n=c, seed=sq))
            if q["klass"] == "quasi_periodic":
                q_ok += 1
        exp_acc[f"n{c}"] = round(e_ok / n_each, 3)
        qp_acc[f"n{c}"] = round(q_ok / n_each, 3)
    return {
        "exponential_accuracy_vs_count": exp_acc,
        "quasiperiodic_completeness_vs_count": qp_acc,
        "exponential_false_positive_rate": round(fp / max(fp_total, 1), 4),
    }


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Offline: synthetic recover-a-known (classify injected exponential/quasi-periodic series)."""
    import json

    if offline:
        metrics: dict = _synthetic_metrics()
    else:  # pragma: no cover - real leg scrapes the live JBO table
        metrics = _real_census(out)

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "glitchpop_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(metrics, op / "papers" / "glitchpop" / "figures")
    _write_macros(metrics, op / "papers" / "glitchpop" / "generated" / "macros.tex")
    return metrics


def _synthetic_metrics() -> dict:
    """Recover-a-known: the classifier labels injected exponential/quasi-periodic/clustered series."""
    rec = inject_recover()
    # the two gap-robust, physically-relevant classes (exponential null + the quasi-periodic signal)
    demos = {
        k: classify_pulsar(synthetic_glitch_series(kind=k, n=30, seed=3))["klass"]
        for k in ("exponential", "quasi_periodic")
    }
    # clustered is recoverable ONLY without monitoring-gap excision (its long gaps look like gaps);
    # the default excision under-detects it, so we demo it with excision off to show the class exists
    demos["clustered_no_excision"] = waiting_time_fit(
        waiting_times(synthetic_glitch_series(kind="clustered", n=30, seed=3)["mjd"]),
        gap_factor=1e9,
    )["klass"]
    return {
        "source": "synthetic glitch series (injected exponential / quasi-periodic / clustered waits)",
        "is_real": False,
        **rec,
        "demo_klass": demos,
        "recovered": bool(
            demos["exponential"] == "exponential" and demos["quasi_periodic"] == "quasi_periodic"
        ),
    }


def _real_census(out: str) -> dict:  # pragma: no cover - network (scrape the live JBO table)
    """Real leg: scrape the live JBO table, classify, and diff against the end-2018 subset."""
    from .glitchpop_real import run_real_census

    return run_real_census(out)


def _figure(m: dict, out_dir: str | Path) -> None:
    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.8))
    # left: the three synthetic waiting-time distributions
    for kind, c in (("exponential", "C0"), ("quasi_periodic", "C1"), ("clustered", "C3")):
        w = waiting_times(synthetic_glitch_series(kind=kind, n=400, seed=5)["mjd"])
        ax1.hist(w / np.median(w), bins=30, histtype="step", color=c, label=kind, density=True)
    ax1.set(
        xlabel="waiting time / median", ylabel="density", title="Injected processes", xlim=(0, 4)
    )
    ax1.legend(fontsize=8)
    # right: recover-a-known accuracy vs glitch count
    acc = m.get("quasiperiodic_completeness_vs_count", {})
    if acc:
        counts = [int(k[1:]) for k in acc]
        ax2.plot(counts, [acc[k] for k in acc], "o-", color="C1", label="quasi-periodic")
        ea = m.get("exponential_accuracy_vs_count", {})
        ax2.plot(counts, [ea[k] for k in ea], "s-", color="C0", label="exponential")
    ax2.set(
        xlabel="glitches per pulsar",
        ylabel="classification accuracy",
        title="Recover-a-known",
        ylim=(0, 1.05),
    )
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "glitchpop.pdf")
    plt.close(fig)


def _write_macros(m: dict, path: str | Path) -> None:
    """Emit both namespaces: gpSyn* (recover-a-known, always live) + gpReal* (the real census)."""

    def g(src: dict, key: str) -> str:
        v = src.get(key)
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return "--"
        return str(v)

    syn = _synthetic_metrics() if m.get("is_real") else m
    real = m if m.get("is_real") else {}
    ea = syn.get("exponential_accuracy_vs_count", {})
    qa = syn.get("quasiperiodic_completeness_vs_count", {})
    lines = [
        "% Auto-generated by jansky_research.glitchpop._write_macros -- do not edit.",
        "% gpSyn* = synthetic recover-a-known (always live); gpReal* = the real JBO census.",
        rf"\newcommand{{\gpSource}}{{{m['source']}}}",
        rf"\newcommand{{\gpSynRecovered}}{{{'yes' if syn.get('recovered') else 'no'}}}",
        rf"\newcommand{{\gpSynFP}}{{{g(syn, 'exponential_false_positive_rate')}}}",
        rf"\newcommand{{\gpSynExpAccLow}}{{{ea.get('n6', '--')}}}",
        rf"\newcommand{{\gpSynExpAccHigh}}{{{ea.get('n40', '--')}}}",
        rf"\newcommand{{\gpSynQpCompLow}}{{{qa.get('n6', '--')}}}",
        rf"\newcommand{{\gpSynQpCompHigh}}{{{qa.get('n40', '--')}}}",
    ]
    for macro, key in (
        ("NGlitches", "n_glitches"),
        ("NPulsars", "n_pulsars"),
        ("NFit", "n_qualified_full"),
        ("NExp", "n_exponential"),
        ("NQp", "n_quasiperiodic"),
        ("NClustered", "n_clustered"),
        ("NNewly", "n_newly_classifiable"),
        ("NStable", "n_stable_sample"),
        ("NFlipped", "n_flipped"),
        ("ExpectedFalseQp", "expected_false_qp"),
        ("NMagnetars", "n_magnetars_dropped"),
        ("KnownQpOK", "known_quasiperiodic_ok"),
    ):
        lines.append(rf"\newcommand{{\gpReal{macro}}}{{{g(real, key)}}}")
    from .report import _fmt_p

    bp = real.get("qp_binomial_p")
    lines.append(
        rf"\newcommand{{\gpRealQpBinomP}}{{{_fmt_p(bp) if isinstance(bp, (int, float)) and np.isfinite(bp) else '--'}}}"
    )
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    ap = argparse.ArgumentParser(description="JBO glitch waiting-time population census.")
    ap.add_argument("--out", default=".")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
