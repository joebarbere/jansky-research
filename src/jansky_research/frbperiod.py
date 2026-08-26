"""FRB repeater activity-periodicity search.

Some repeating fast radio bursts show **periodic activity windows** — bursts cluster at a preferred
phase of a multi-day cycle. The archetype is FRB 20180916B, with a 16.35-day period (CHIME/FRB
Collaboration 2020, Nature 582, 351). This module searches a repeater's burst arrival times (MJDs)
for such a period with a phase-folding **Rayleigh ($Z^2_1$) periodogram** (Buccheri et al. 1983,
A&A 128, 245) — pure NumPy, CPU-only. (CHIME report 16.35 +/- 0.15 d.)

Honest scope (and the write-up must say this): the CHIME catalogue is a **transit survey** that
sees each source roughly once per sidereal day with strongly non-uniform exposure. A catalogue-only
periodogram therefore (a) is heavily aliased near 1 day and its beats, and (b) cannot give a
rigorous false-alarm probability without the survey's exposure model. We can *recover* a known
period as a peak and set rough limits, but significance from the simple analytic FAP below is an
*upper bound on confidence* only — it ignores exposure and aliasing.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "PeriodogramResult",
    "collapse_transits",
    "detection_power",
    "false_alarm_prob",
    "mc_null",
    "measured_duty_cycle",
    "period_bootstrap",
    "period_search",
    "rayleigh_z2",
    "run",
    "search_repeaters",
    "synthetic_periodic_arrivals",
]


@dataclass(frozen=True)
class PeriodogramResult:
    """Result of a Rayleigh periodogram search over a grid of trial periods (days)."""

    periods: np.ndarray
    z2: np.ndarray
    best_period: float
    best_z2: float
    n_bursts: int
    fap: float  # analytic false-alarm probability of the peak (exposure-blind; an upper bound)


def rayleigh_z2(times: np.ndarray, period: float) -> float:
    """Rayleigh $Z^2_1$ statistic for phase concentration at a trial ``period``.

    Folds the arrival times to phase $\\phi_i = 2\\pi\\,(t_i/P \\bmod 1)$ and returns
    $Z^2_1 = (2/n)\\,[(\\sum\\cos\\phi_i)^2 + (\\sum\\sin\\phi_i)^2]$. It is large when the bursts
    cluster at one phase (a periodic activity window) and $\\sim 2$ for random phases. Under the
    no-signal hypothesis $Z^2_1$ follows a $\\chi^2$ with 2 degrees of freedom.
    """
    t = np.asarray(times, dtype=float)
    phi = 2.0 * np.pi * ((t / period) % 1.0)
    c = np.cos(phi).sum()
    s = np.sin(phi).sum()
    return float((2.0 / t.size) * (c * c + s * s))


def false_alarm_prob(z2_max: float, n_indep: int) -> float:
    """Exposure-blind false-alarm probability of a peak $Z^2_{\\max}$ over ``n_indep`` trials.

    For a single trial, $P(Z^2_1 > z) = e^{-z/2}$ ($\\chi^2_2$ survival). With ``n_indep``
    independent trial periods, $\\mathrm{FAP} = 1 - (1 - e^{-z/2})^{n_\\mathrm{indep}}$. This assumes
    the trial frequencies are independent — which a transit survey's aliased spectral window
    violates — and ignores the exposure, so it is an **approximate** number, not a rigorous
    significance. Use it only as a rough guide; trust the match to an independently-known period.
    """
    p_single = np.exp(-z2_max / 2.0)
    return float(1.0 - (1.0 - p_single) ** max(n_indep, 1))


def period_search(
    times: np.ndarray, periods: np.ndarray, *, span: float | None = None
) -> PeriodogramResult:
    """Rayleigh periodogram over ``periods`` (days); returns the peak period and its FAP.

    ``n_indep`` (for the FAP) is estimated as the number of independent frequencies across the data
    span, $\\sim (1/P_\\min - 1/P_\\max)\\,T$, capped by the grid size.
    """
    t = np.asarray(times, dtype=float)
    z2 = np.array([rayleigh_z2(t, P) for P in periods])
    k = int(np.argmax(z2))
    if span is None:
        span = float(t.max() - t.min())
    n_indep = min(int((1.0 / periods.min() - 1.0 / periods.max()) * span) + 1, periods.size)
    return PeriodogramResult(
        periods=periods,
        z2=z2,
        best_period=float(periods[k]),
        best_z2=float(z2[k]),
        n_bursts=int(t.size),
        fap=false_alarm_prob(float(z2[k]), n_indep),
    )


def collapse_transits(times: np.ndarray, *, tol_days: float = 0.05) -> np.ndarray:
    """One epoch per transit: bursts separated by less than ``tol_days`` merge to their mean.

    CHIME sees a source once per sidereal day; several catalogued bursts minutes apart in one
    transit are duplicate phase samples at a multi-day trial period, not independent draws, and
    feeding them to $Z^2_1$ as independent terms inflates the statistic (the round-10 referee
    measured $Z^2$ 33.4 $\\to$ 23.9 for FRB 20180916B: 19 bursts on 14 distinct transits). The
    same pseudo-replication is already collapsed for sub-bursts upstream; this closes it for
    distinct same-transit bursts.
    """
    t = np.sort(np.asarray(times, float))
    if t.size == 0:
        return t
    groups = [[t[0]]]
    for x in t[1:]:
        if x - groups[-1][-1] < tol_days:
            groups[-1].append(x)
        else:
            groups.append([x])
    return np.array([float(np.mean(g)) for g in groups])


def _z2_grid(times: np.ndarray, periods: np.ndarray) -> np.ndarray:
    """Vectorised $Z^2_1$ over a period grid (times x periods trig sums)."""
    t = np.asarray(times, float)
    phi = 2.0 * np.pi * ((t[:, None] / periods[None, :]) % 1.0)
    c = np.cos(phi).sum(axis=0)
    s = np.sin(phi).sum(axis=0)
    return (2.0 / t.size) * (c * c + s * s)


def measured_duty_cycle(times: np.ndarray, period: float) -> float:
    """Activity duty cycle at ``period``: the minimal phase arc containing every burst."""
    t = np.asarray(times, float)
    ph = np.sort((t / period) % 1.0)
    if ph.size < 2:
        return 1.0
    gaps = np.diff(np.concatenate([ph, [ph[0] + 1.0]]))
    return float(1.0 - gaps.max())


def mc_null(
    n_epochs: int,
    span: float,
    periods: np.ndarray,
    *,
    donor_intervals: np.ndarray | None = None,
    n_trials: int = 2000,
    seed: int = 0,
) -> dict:
    """Monte-Carlo null for the grid-max $Z^2_1$: uniform arrivals, or clustered ones.

    ``donor_intervals=None`` draws ``n_epochs`` uniform arrival times over ``span`` --- the
    random-phase null, against which the analytic FAP is roughly calibrated. Supplying another
    repeater's inter-burst intervals bootstraps a CLUSTERED, aperiodic null (resampled with
    replacement, rescaled to ``span``): burst clustering, not aliasing, is what breaks the
    analytic FAP (measured: it raises the null's median max-$Z^2$ by ~30%). Never build the null
    by shuffling the target's OWN intervals --- for a genuinely periodic source the shuffled
    intervals remain near-multiples of the period and the "null" inherits the signal. Returns the
    per-trial grid-max and its period, from which exceedance and coincidence probabilities follow.
    """
    rng = np.random.default_rng(seed)
    maxima = np.empty(n_trials)
    argp = np.empty(n_trials)
    for i in range(n_trials):
        if donor_intervals is None:
            t = np.sort(rng.uniform(0.0, span, n_epochs))
        else:
            iv = rng.choice(donor_intervals, n_epochs - 1, replace=True)
            iv = iv * (span / iv.sum())
            t = np.concatenate([[0.0], np.cumsum(iv)])
        z2 = _z2_grid(t, periods)
        k = int(np.argmax(z2))
        maxima[i] = z2[k]
        argp[i] = periods[k]
    return {"max_z2": maxima, "argmax_period": argp, "n_trials": n_trials}


def null_summary(
    null: dict, *, observed_z2: float, target_period: float, target_tol: float
) -> dict:
    """Exceedance and coincidence probabilities from an :func:`mc_null` draw.

    The coincidence probability --- how often a null realization's grid peak lands within
    ``target_tol`` of the independently published period --- is the number the recover-a-known
    argument actually rests on, and the one the paper's headline quotes.
    """
    m = null["max_z2"]
    n = m.size
    p_exceed = float((np.sum(m >= observed_z2) + 1) / (n + 1))
    p_coinc = float(
        (np.sum(np.abs(null["argmax_period"] - target_period) <= target_tol) + 1) / (n + 1)
    )
    return {
        "n_trials": n,
        "p_exceed": float(f"{p_exceed:.2e}"),
        "p_coincidence": float(f"{p_coinc:.2e}"),
        "median_max_z2": round(float(np.median(m)), 2),
        "p99_max_z2": round(float(np.percentile(m, 99)), 2),
    }


def detection_power(
    n_epochs: int,
    span: float,
    periods: np.ndarray,
    *,
    period: float = 16.33,
    duty_cycles: tuple[float, ...] = (0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
    fap_threshold: float = 0.01,
    n_inj: int = 200,
    seed: int = 0,
) -> list[dict]:
    """Power of the search against an injected periodic activity window, per duty cycle.

    "Correctly finds nothing elsewhere" is a sensitivity claim; this measures it. Bursts are
    injected at random cycles of ``period`` within ``span``, at phases uniform inside the window,
    and scored by the run's own rule (grid-max FAP below ``fap_threshold``).
    """
    rng = np.random.default_rng(seed)
    span_eff = span
    n_indep = min(int((1.0 / periods.min() - 1.0 / periods.max()) * span_eff) + 1, periods.size)
    out = []
    for duty in duty_cycles:
        n_det = 0
        n_cycles = max(int(span_eff / period), 1)
        for _ in range(n_inj):
            cyc = rng.integers(0, n_cycles, n_epochs)
            ph = rng.uniform(0.0, duty, n_epochs)
            t = np.sort((cyc + ph) * period)
            z2 = _z2_grid(t, periods)
            if false_alarm_prob(float(z2.max()), n_indep) < fap_threshold:
                n_det += 1
        out.append({"duty_cycle": duty, "power": round(n_det / n_inj, 3)})
    return out


def period_bootstrap(
    times: np.ndarray, periods: np.ndarray, *, n_boot: int = 500, seed: int = 0
) -> dict:
    """Bootstrap uncertainty on the peak period (epochs resampled with replacement).

    The grid step is not the measurement error; this is.
    """
    rng = np.random.default_rng(seed)
    t = np.asarray(times, float)
    peaks = []
    for _ in range(n_boot):
        tb = t[rng.integers(0, t.size, t.size)]
        z2 = _z2_grid(tb, periods)
        peaks.append(float(periods[int(np.argmax(z2))]))
    arr = np.asarray(peaks)
    med = float(np.median(arr))
    near = arr[np.abs(arr - med) < 2.0]  # exclude rare alias jumps from the width estimate
    return {
        "n_boot": n_boot,
        "period_sd": round(float(np.std(near, ddof=1)), 3),
        "frac_near_peak": round(near.size / arr.size, 3),
    }


def synthetic_periodic_arrivals(
    period: float = 16.0,
    n: int = 40,
    active_frac: float = 0.3,
    span: float = 350.0,
    seed: int | None = 0,
) -> np.ndarray:
    """Generate burst MJDs clustered in one phase window of ``period`` (offline test fixture).

    Each burst lands in a random cycle within ``span`` days, at a phase uniformly inside
    ``[0, active_frac]`` of the period — i.e. a periodic activity window. Used by the tests and as
    the offline fallback.
    """
    rng = np.random.default_rng(seed)
    n_cycles = max(int(span / period), 1)
    cycle = rng.integers(0, n_cycles, n)
    phase = rng.uniform(0.0, active_frac, n)
    return np.sort((cycle + phase) * period)


def _synthetic_repeaters() -> tuple[np.ndarray, np.ndarray]:
    """Two synthetic repeaters — one 16.35-day periodic, one random — for the offline run/tests."""
    a = synthetic_periodic_arrivals(16.35, n=40, active_frac=0.15, span=400, seed=0)
    rng = np.random.default_rng(1)
    b = np.sort(rng.uniform(0.0, 400.0, 15))
    mjd = 58000.0 + np.concatenate([a, b])
    names = np.array(["SYN-PER"] * a.size + ["SYN-RND"] * b.size)
    return mjd, names


def search_repeaters(
    mjd: np.ndarray, names: np.ndarray, *, min_bursts: int = 8, periods: np.ndarray | None = None
) -> tuple[list[dict], np.ndarray]:
    """Run the periodogram per repeater source; sources below ``min_bursts`` are skipped.

    Returns ``(rows, periods)`` where each row is ``{name, n, searched, best_period, z2, fap}``.
    """
    if periods is None:
        periods = np.linspace(2.0, 100.0, 12000)
    mjd = np.asarray(mjd, dtype=float)
    names = np.asarray(names)
    rows: list[dict] = []
    for nm in sorted(set(names.tolist())):
        t = np.sort(mjd[names == nm])
        if t.size < min_bursts:
            rows.append(
                {
                    "name": nm,
                    "n": int(t.size),
                    "n_epochs": int(collapse_transits(t).size),
                    "searched": False,
                    "best_period": None,
                    "z2": None,
                    "fap": None,
                }
            )
            continue
        # the searched series is ONE EPOCH PER TRANSIT (same-transit bursts are duplicate
        # phase samples at multi-day periods, not independent draws)
        te = collapse_transits(t)
        r = period_search(te, periods)
        rows.append(
            {
                "name": nm,
                "n": int(t.size),
                "n_epochs": int(te.size),
                "searched": True,
                "best_period": r.best_period,
                "z2": r.best_z2,
                "fap": r.fap,
                "span_days": round(float(te.max() - te.min()), 2),
            }
        )
    return rows, periods


def run(
    out: str = ".",
    *,
    offline: bool = False,
    min_bursts: int = 8,
    mc_trials: int | None = None,
    n_inj: int | None = None,
) -> dict:
    """Search every catalogue repeater for activity periodicity; write results + a periodogram.

    Writes ``results/period_metrics.json`` and a periodogram figure for the most significant
    detection. The per-source CSV goes to ``survey/period_results.csv`` on a real run (the committed,
    showcased real-data table) but to ``results/period_results.csv`` (git-ignored) when ``offline`` ---
    so a synthetic run never clobbers the tracked real results. Returns the metrics dict.
    """
    import csv
    from pathlib import Path

    if offline:
        mjd, names = _synthetic_repeaters()
        source = "synthetic"
    else:  # pragma: no cover - network
        from . import pipeline

        cat, source = pipeline.build_catalog(offline=False)
        rep = np.asarray(cat["repeater"], dtype=bool)
        mjd = np.asarray(cat["mjd"])[rep]
        names = np.asarray(cat["repeater_name"])[rep]

    rows, periods = search_repeaters(mjd, names, min_bursts=min_bursts)
    searched = [r for r in rows if r["searched"]]
    significant = sorted((r for r in searched if r["fap"] < 0.01), key=lambda r: r["fap"])
    fap_threshold = 0.01
    # the MC blocks are cheap enough at full depth for the real run; the offline/CI leg uses a
    # shallow draw so the smoke build stays fast (depth is recorded in the metrics)
    mc_trials = mc_trials if mc_trials is not None else (200 if offline else 2000)
    n_inj = n_inj if n_inj is not None else (50 if offline else 200)
    metrics = {
        "source": source,
        "n_sources": len(rows),
        "n_searchable": len(searched),
        "n_significant": len(significant),
        "min_bursts": min_bursts,
        "fap_threshold": fap_threshold,
        "grid": {
            "p_min_days": float(periods.min()),
            "p_max_days": float(periods.max()),
            "n_grid": int(periods.size),
        },
        # ALL searched sources, nulls included -- a metrics file that records only detections
        # leaves the null outside every automated check (the innerrc lesson)
        "searched": [
            {
                "name": r["name"],
                "n_bursts": r["n"],
                "n_epochs": r["n_epochs"],
                "best_period_days": round(r["best_period"], 3),
                "z2": round(r["z2"], 1),
                "fap": float(f"{r['fap']:.2e}"),
                "span_days": r["span_days"],
            }
            for r in searched
        ],
        "detections": [
            {
                "name": r["name"],
                "period_days": round(r["best_period"], 3),
                "z2": round(r["z2"], 1),
                "fap": r["fap"],
                "n": r["n"],
                "n_epochs": r["n_epochs"],
            }
            for r in significant
        ],
    }

    if significant:
        det = significant[0]
        t_det = collapse_transits(np.sort(mjd[names == det["name"]]))
        span = float(t_det.max() - t_det.min())
        duty = measured_duty_cycle(t_det, det["best_period"])
        metrics["duty_cycle_at_peak"] = round(duty, 3)
        metrics["period_bootstrap"] = period_bootstrap(t_det, periods)
        # Monte-Carlo nulls: uniform, and clustered (bootstrapping ANOTHER searched source's
        # intervals -- never the target's own, which would inherit the periodicity)
        nulls: dict = {}
        uni = mc_null(t_det.size, span, periods, n_trials=mc_trials)
        nulls["uniform"] = null_summary(
            uni, observed_z2=det["z2"], target_period=16.35, target_tol=0.15
        )
        donors = [r for r in searched if r["name"] != det["name"]]
        if donors:
            donor = donors[0]
            t_don = collapse_transits(np.sort(mjd[names == donor["name"]]))
            iv = np.diff(t_don)
            clu = mc_null(t_det.size, span, periods, donor_intervals=iv, n_trials=mc_trials, seed=1)
            nulls["clustered"] = null_summary(
                clu, observed_z2=det["z2"], target_period=16.35, target_tol=0.15
            )
            nulls["clustered_donor"] = donor["name"]
            # the sensitivity of the null source's own search: could it have seen a twin?
            metrics["null_source_power"] = {
                "name": donor["name"],
                "twin_duty_cycle": round(duty, 3),
                "power": detection_power(
                    donor["n_epochs"],
                    donor["span_days"],
                    periods,
                    period=det["best_period"],
                    duty_cycles=(0.05, 0.1, 0.15, 0.2, round(duty, 3), 0.3),
                    fap_threshold=fap_threshold,
                    n_inj=n_inj,
                ),
            }
        metrics["mc_nulls"] = nulls

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import _results_are_real, write_results

    json_path = op / "results" / "period_metrics.json"
    write_artifacts = True
    if source == "synthetic":
        try:
            import json as _json

            write_artifacts = not (
                json_path.is_file() and _results_are_real(_json.loads(json_path.read_text()))
            )
        except Exception:
            write_artifacts = True

    write_results(metrics, json_path)
    # real run -> the committed showcase CSV under survey/ plus the per-burst epochs table (the
    # committed evidence the null reproduction needs); offline -> a clearly-synthetic-named,
    # git-ignored file so the synthetic run cannot overwrite or masquerade as the real table
    csv_dir = op / ("results" if offline else "survey")
    csv_dir.mkdir(parents=True, exist_ok=True)
    csv_name = "period_results_synthetic.csv" if offline else "period_results.csv"
    with open(csv_dir / csv_name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["repeater", "n_bursts", "n_epochs", "searched", "best_period_days", "z2", "fap"]
        )
        for r in sorted(rows, key=lambda r: -r["n"]):
            w.writerow(
                [
                    r["name"],
                    r["n"],
                    r["n_epochs"],
                    r["searched"],
                    "" if r["best_period"] is None else f"{r['best_period']:.3f}",
                    "" if r["z2"] is None else f"{r['z2']:.1f}",
                    "" if r["fap"] is None else f"{r['fap']:.2e}",
                ]
            )
    if not offline:  # pragma: no cover - network run
        with open(op / "results" / "period_epochs.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["repeater", "mjd"])
            for nm in sorted(set(np.asarray(names).tolist())):
                for x in np.sort(np.asarray(mjd, float)[np.asarray(names) == nm]):
                    w.writerow([nm, f"{x:.6f}"])
    if significant and write_artifacts:
        _periodogram_figure(
            mjd, names, significant[0], periods, op / "papers" / "frbperiod" / "figures"
        )
    _write_macros(metrics, op / "papers" / "frbperiod" / "generated" / "macros.tex")
    return metrics


def _write_macros(m: dict, path) -> None:
    """Emit LaTeX ``\\newcommand`` macros so the paper hard-codes no number (offline + real share this)."""
    from pathlib import Path

    from .report import _fmt_p

    grid = m.get("grid") or {}
    lines = [
        "% Auto-generated by jansky_research.frbperiod._write_macros — do not edit by hand.",
        rf"\newcommand{{\fpSource}}{{{m['source']}}}",
        rf"\newcommand{{\fpNsources}}{{{m['n_sources']}}}",
        rf"\newcommand{{\fpNsearch}}{{{m['n_searchable']}}}",
        rf"\newcommand{{\fpNsig}}{{{m['n_significant']}}}",
        rf"\newcommand{{\fpNunsearchable}}{{{m['n_sources'] - m['n_searchable']}}}",
        rf"\newcommand{{\fpMinBursts}}{{{m['min_bursts']}}}",
        rf"\newcommand{{\fpFapThreshold}}{{{m.get('fap_threshold', '--')}}}",
        rf"\newcommand{{\fpPmin}}{{{grid.get('p_min_days', '--'):g}}}"
        if grid
        else r"\newcommand{\fpPmin}{--}",
        rf"\newcommand{{\fpPmax}}{{{grid.get('p_max_days', '--'):g}}}"
        if grid
        else r"\newcommand{\fpPmax}{--}",
        rf"\newcommand{{\fpNgrid}}{{{grid.get('n_grid', '--')}}}"
        if grid
        else r"\newcommand{\fpNgrid}{--}",
    ]
    if m["detections"]:
        d = m["detections"][0]
        lines += [
            rf"\newcommand{{\fpName}}{{{d['name']}}}",
            rf"\newcommand{{\fpPeriod}}{{{d['period_days']:.2f}}}",
            rf"\newcommand{{\fpZtwo}}{{{d['z2']:.1f}}}",
            rf"\newcommand{{\fpFAP}}{{{_fmt_p(d['fap'])}}}",
            rf"\newcommand{{\fpNbursts}}{{{d['n']}}}",
            rf"\newcommand{{\fpNepochs}}{{{d.get('n_epochs', '--')}}}",
        ]

    def _g(dic: dict, key: str) -> str:
        v = dic.get(key)
        return "--" if v is None else str(v)

    nulls = m.get("mc_nulls") or {}
    uni = nulls.get("uniform") or {}
    clu = nulls.get("clustered") or {}
    pb = m.get("period_bootstrap") or {}
    nsp = m.get("null_source_power") or {}
    ptab = {row["duty_cycle"]: row["power"] for row in nsp.get("power") or []}
    twin = nsp.get("twin_duty_cycle")
    lines += [
        rf"\newcommand{{\fpDuty}}{{{_g(m, 'duty_cycle_at_peak')}}}",
        rf"\newcommand{{\fpPerr}}{{{_g(pb, 'period_sd')}}}",
        rf"\newcommand{{\fpUniPexceed}}{{{_fmt_p(uni['p_exceed']) if uni else '--'}}}",
        rf"\newcommand{{\fpUniPcoinc}}{{{_fmt_p(uni['p_coincidence']) if uni else '--'}}}",
        rf"\newcommand{{\fpCluPexceed}}{{{_fmt_p(clu['p_exceed']) if clu else '--'}}}",
        rf"\newcommand{{\fpCluPcoinc}}{{{_fmt_p(clu['p_coincidence']) if clu else '--'}}}",
        rf"\newcommand{{\fpCluMedZ}}{{{_g(clu, 'median_max_z2')}}}",
        rf"\newcommand{{\fpUniMedZ}}{{{_g(uni, 'median_max_z2')}}}",
        rf"\newcommand{{\fpNullName}}{{{_g(nulls, 'clustered_donor')}}}",
        rf"\newcommand{{\fpTwinPower}}{{{ptab.get(twin, '--') if twin is not None else '--'}}}",
    ]
    # the null source's own numbers, so the paper can state what it searched and found
    null_rows = [
        r
        for r in m.get("searched") or []
        if not any(d["name"] == r["name"] for d in m.get("detections") or [])
    ]
    nr = null_rows[0] if null_rows else {}
    lines += [
        rf"\newcommand{{\fpNullZtwo}}{{{_g(nr, 'z2')}}}",
        rf"\newcommand{{\fpNullFAP}}{{{_fmt_p(nr['fap']) if nr.get('fap') is not None else '--'}}}",
        rf"\newcommand{{\fpNullPeriod}}{{{_g(nr, 'best_period_days')}}}",
        rf"\newcommand{{\fpNullNbursts}}{{{_g(nr, 'n_bursts')}}}",
        rf"\newcommand{{\fpNullNepochs}}{{{_g(nr, 'n_epochs')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and
    # would otherwise blank the other mode's macros with '--'. `make figures`
    # runs every slice offline in the repo root, so without this an offline
    # rebuild silently empties this paper. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _periodogram_figure(mjd, names, det, periods, out_dir):
    """Plot the Rayleigh periodogram of the most significant detection, marking its peak."""
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    t = collapse_transits(np.sort(np.asarray(mjd)[np.asarray(names) == det["name"]]))
    z2 = _z2_grid(t, np.asarray(periods))
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(periods, z2, lw=0.8)
    ax.axvline(det["best_period"], color="r", ls="--", label=f"peak {det['best_period']:.2f} d")
    ax.set(
        xlabel="trial period (days)",
        ylabel=r"Rayleigh $Z^2_1$",
        title=f"{det['name']} activity periodogram ({t.size} transit epochs)",
    )
    ax.legend()
    p = out / "periodogram.pdf"
    fig.tight_layout()
    fig.savefig(p)
    plt.close(fig)
    return p


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="Search CHIME repeaters for activity periodicity.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true", help="use the synthetic fixture (no network)")
    p.add_argument("--min-bursts", type=int, default=8)
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline, min_bursts=args.min_bursts), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
