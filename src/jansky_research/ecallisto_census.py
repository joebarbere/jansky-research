"""e-Callisto type III occurrence census: does the burst rate track the solar cycle?

Coincidence-confirmed type III bursts (from :mod:`jansky_research.ecallisto_catalog`) are the raw
material for an occurrence census. The first-order, well-established expectation is that the type III
rate --- electron beams from flares --- **tracks solar activity**: it rises and falls with the sunspot
number over the ~11-year cycle. This module builds that census: a coverage-corrected burst rate per
period, correlated against the SILSO sunspot number.

The honest catch is **completeness**: a burst is confirmed only if enough stations observed it, and the
active-station count varies with time, so the raw event count must be normalised by coverage before it
can be compared across epochs. We do that (rate = events / active-station coverage) and correlate the
corrected rate with the sunspot number. Pure NumPy with a synthetic offline fixture whose event stream
is generated from a synthetic solar cycle (so the census round-trips); the real run fetches the SILSO
sunspot series and samples e-Callisto days (heavy --- the Airflow ingest of ``plans/31`` is built for it).
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "audit_day_listings",
    "census_correlation",
    "coverage_corrected_rate",
    "fetch_sunspots",
    "growth_coverage",
    "parse_silso",
    "run",
    "synthetic_census",
    "synthetic_sunspots",
    "validation_suite",
]

SILSO_URL = "https://www.sidc.be/SILSO/DATA/SN_m_tot_V2.0.csv"


def coverage_corrected_rate(n_events: np.ndarray, coverage: np.ndarray) -> np.ndarray:
    r"""Coverage-corrected occurrence rate: confirmed events per unit active-station coverage.

    A burst is confirmed only where enough stations observe it, and the active-station count varies with
    time, so the raw count must be divided by the coverage before epochs can be compared:
    $\mathrm{rate}=N_\mathrm{events}/C$. Periods with zero coverage give NaN.
    """
    n = np.asarray(n_events, float)
    c = np.asarray(coverage, float)
    out = np.full(np.broadcast(n, c).shape, np.nan)
    pos = c > 0
    return np.divide(n, c, out=out, where=pos)


def census_correlation(rate: np.ndarray, sunspot: np.ndarray) -> dict:
    r"""Correlate the coverage-corrected burst rate with the sunspot number.

    Returns the Pearson $r$ and Spearman $\rho$ (rank, robust to the non-linear activity--rate relation)
    of ``rate`` versus ``sunspot`` over the periods where both are finite, plus the ordinary-least-squares
    slope of rate on sunspot and the number of periods used.
    """
    from scipy import stats as _stats

    r = np.asarray(rate, float)
    s = np.asarray(sunspot, float)
    good = np.isfinite(r) & np.isfinite(s)
    r, s = r[good], s[good]
    nan = float("nan")
    if r.size < 5 or np.ptp(r) == 0 or np.ptp(s) == 0:
        return {"n_periods": int(r.size), "pearson_r": nan, "spearman_rho": nan, "slope": nan}
    slope = float(np.polyfit(s, r, 1)[0])
    return {
        "n_periods": int(r.size),
        "pearson_r": float(np.corrcoef(r, s)[0, 1]),
        "spearman_rho": float(_stats.spearmanr(r, s).statistic),
        "slope": slope,
    }


def synthetic_sunspots(
    n_months: int = 180, *, amplitude: float = 140.0, seed: int = 0
) -> np.ndarray:
    """A synthetic monthly sunspot-number series: a ~11-year cycle with realistic asymmetry + noise.

    Models the sunspot number as a raised, fast-rise/slow-decay cycle (period 132 months) so the
    synthetic census has a realistic activity curve to track. Returns non-negative monthly values.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(n_months, dtype=float)
    phase = (t % 132.0) / 132.0
    cycle = (
        np.exp(-3.0 * phase) * (1.0 - np.cos(2.0 * np.pi * phase)) / 2.0
    )  # fast rise, slow decay
    base = amplitude * cycle / cycle.max()
    return np.clip(base + rng.normal(0.0, 5.0, n_months), 0.0, None)


def synthetic_census(
    sunspot: np.ndarray,
    *,
    k: float = 0.03,
    coverage_mean: float = 12.0,
    coverage: np.ndarray | None = None,
    c_half: float | None = None,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Synthetic monthly event counts whose true rate is proportional to the sunspot number.

    For each month with sunspot $S$ and active-station coverage $C$, the confirmed-event count is
    Poisson with mean $k\,S\,C_\mathrm{eff}$. With the default ``c_half=None`` the confirmation is
    *linear* in coverage ($C_\mathrm{eff}=C$) --- the same model the $N/C$ estimator assumes, so a
    recovery on this arm tests only arithmetic and Poisson robustness, never the correction itself
    (the round-10 referee's circularity finding). ``c_half`` switches to a **saturating**
    confirmation, $C_\mathrm{eff}=c_\mathrm{half}(1-e^{-C/c_\mathrm{half}})$ --- the shape a fixed
    two-station coincidence threshold actually produces --- under which $N/C$ over-corrects and the
    validation CAN fail; :func:`validation_suite` measures by how much. ``coverage`` supplies an
    explicit station history (e.g. :func:`growth_coverage`); default is the stationary jittered one.
    Returns ``(n_events, coverage)`` per month.
    """
    rng = np.random.default_rng(seed)
    s = np.asarray(sunspot, float)
    if coverage is None:
        coverage = np.clip(rng.normal(coverage_mean, coverage_mean * 0.25, s.size), 2.0, None)
    else:
        coverage = np.asarray(coverage, float)
    ceff = coverage if c_half is None else c_half * (1.0 - np.exp(-coverage / c_half))
    lam = k * s * ceff
    n_events = rng.poisson(np.clip(lam, 0.0, None)).astype(float)
    return n_events, coverage


def growth_coverage(
    n_months: int, *, lo: float = 2.0, hi: float = 60.0, jitter: float = 0.15, seed: int = 0
) -> np.ndarray:
    """A station history shaped like the real network's: monotone growth ``lo`` → ``hi`` plus jitter.

    The shipped stationary fixture (mean 12, ±25%) puts the coverage confound nowhere near the
    decision boundary, so the raw count already correlates with activity and the correction is
    cosmetic there. Under a growth history the raw correlation genuinely breaks and the correction
    genuinely rescues it --- the demonstration the paper's prose claims.
    """
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 1.0, n_months)
    base = lo + (hi - lo) * t
    return np.clip(base * (1.0 + rng.normal(0.0, jitter, n_months)), 1.0, None)


def validation_suite(
    *, n_months: int = 180, k: float = 0.03, seed: int = 0, n_seeds: int = 30
) -> dict:
    """Every validation arm the census statistic needs, including the ones that can fail.

    Arms: (1) the linear/stationary arm (the original recover-a-known --- correct-by-construction
    for the estimator, kept as the arithmetic check, now reported WITH the uncorrected correlation
    and the residual coverage correlation); (2) a growth-history arm where the correction genuinely
    rescues a broken raw correlation; (3) a **misspecification** arm --- saturating confirmation,
    the shape the pipeline's fixed two-station threshold produces --- where N/C over-corrects and
    the corrected rate acquires a spurious anti-correlation with coverage, measured and committed;
    (4) a seed ensemble for the slope (realization variance, not just one draw). Deterministic.
    """
    sunspot = synthetic_sunspots(n_months, seed=seed)

    def _arm(n_events: np.ndarray, coverage: np.ndarray) -> dict:
        rate = coverage_corrected_rate(n_events, coverage)
        c = census_correlation(rate, sunspot)
        raw = census_correlation(n_events, sunspot)
        good = np.isfinite(rate)
        cov_corr = float(np.corrcoef(rate[good], coverage[good])[0, 1]) if good.sum() > 4 else None
        return {
            "n_periods": c["n_periods"],
            "n_events_total": int(np.nansum(n_events)),
            "pearson_r": round(c["pearson_r"], 3),
            "spearman_rho": round(c["spearman_rho"], 3),
            "slope": round(c["slope"], 4),
            "raw_pearson_r": round(raw["pearson_r"], 3),
            "rate_coverage_corr": round(cov_corr, 3) if cov_corr is not None else None,
        }

    out: dict = {"k_true": k}
    n, c = synthetic_census(sunspot, k=k, seed=seed)
    out["linear"] = _arm(n, c)
    gcov = growth_coverage(n_months, seed=seed)
    ng, _ = synthetic_census(sunspot, k=k, coverage=gcov, seed=seed)
    out["growth"] = _arm(ng, gcov)
    ns, _ = synthetic_census(sunspot, k=k, coverage=gcov, c_half=3.0, seed=seed)
    out["growth_saturating"] = _arm(ns, gcov)
    slopes = []
    for sd in range(n_seeds):
        s_i = synthetic_sunspots(n_months, seed=sd)
        n_i, c_i = synthetic_census(s_i, k=k, seed=sd)
        slopes.append(census_correlation(coverage_corrected_rate(n_i, c_i), s_i)["slope"])
    arr = np.asarray(slopes, float)
    out["slope_ensemble"] = {
        "n_seeds": n_seeds,
        "mean": round(float(np.mean(arr)), 4),
        "sd": round(float(np.std(arr, ddof=1)), 4),
    }
    return out


def parse_silso(text: str) -> dict:
    """Parse the SILSO monthly-mean total sunspot number CSV → ``decimal_year`` and ``sunspot`` arrays.

    Format: ``year;month;decimal_year;sunspot;...`` (``-1`` = missing). Returns finite months only.
    """
    dy, sn = [], []
    for line in text.splitlines():
        parts = line.split(";")
        if len(parts) < 4:
            continue
        try:
            d, s = float(parts[2]), float(parts[3])
        except ValueError:
            continue
        if s >= 0:
            dy.append(d)
            sn.append(s)
    return {"decimal_year": np.asarray(dy, float), "sunspot": np.asarray(sn, float)}


def fetch_sunspots() -> dict:  # pragma: no cover - network
    """Fetch the SILSO monthly-mean total sunspot number (public, no auth). See :func:`parse_silso`."""
    import requests

    return parse_silso(requests.get(SILSO_URL, timeout=60).text)


def sample_real_days(
    dates: list[str], *, window_hours: tuple[int, int] = (9, 13), dt_tol_s: float = 60.0
) -> list[dict]:  # pragma: no cover - network
    """Scan a fixed UT window of each e-Callisto day → per-day confirmed events + station coverage.

    For a consistent occurrence-rate proxy, only files whose start hour is in ``window_hours`` are
    scanned (a fixed sunlit window, sampled identically every day). Each row now carries its own
    provenance --- ``n_files_listed`` / ``n_stations_listed`` (from the day index) alongside
    ``n_files_fetched`` / ``coverage`` (what actually downloaded and parsed) --- so a day where the
    ingest failed is distinguishable from a day where nobody observed. The round-10 referee found
    123 of 168 committed days with coverage 0 that the live archive lists 26--45 stations for: the
    per-file index re-download inside ``solarbursts.fetch_ecallisto`` was throttled and every
    failure vanished into a bare except. A day with files listed but none fetched now records
    coverage NaN (ingest failure), never 0. ``dt_tol_s`` defaults to the pipeline paper's published
    60 s (the earlier run used 120 s, which is recorded with the committed evidence).
    """
    from . import ecallisto_catalog as ec
    from . import solarbursts

    rows_out = []
    h0, h1 = window_hours
    for date in dates:
        listed = [
            (s, f) for (s, f) in ec.list_day_files(date) if h0 <= int(f.split("_")[2][:2]) < h1
        ]
        rows = []
        n_fetched = 0
        for station, fname in listed:
            hhmmss = fname.split("_")[2]
            try:
                spec = solarbursts.fetch_ecallisto(station, date, hhmmss[:4])
            except Exception:
                continue
            n_fetched += 1
            r = ec.scan_spectrum(spec)
            r["station"] = station
            if r.get("t_peak_s") is not None:
                start = int(hhmmss[:2]) * 3600 + int(hhmmss[2:4]) * 60 + int(hhmmss[4:6])
                r["t_peak_s"] = round(start + r["t_peak_s"], 1)
            rows.append(r)
        events = ec.coincident_events(rows, dt_tol_s=dt_tol_s)
        coverage: float = len({r["station"] for r in rows})
        if listed and n_fetched == 0:
            coverage = float("nan")  # ingest failure, not an empty sky
        rows_out.append(
            {
                "date": date,
                "n_events": len(events),
                "coverage": coverage,
                "n_files_listed": len(listed),
                "n_stations_listed": len({s for s, _ in listed}),
                "n_files_fetched": n_fetched,
            }
        )
    return rows_out


def audit_day_listings(
    dates: list[str], *, window_hours: tuple[int, int] = (9, 13), pause: float = 0.3
) -> list[dict]:  # pragma: no cover - network
    """Light audit: for each date, what the archive LISTS in the window (one index GET per day).

    No spectra are downloaded. This is the cheap check that separates "the sky was empty / nobody
    observed" from "the ingest failed": a committed census day with coverage 0 but a non-zero
    listing here was a silent fetch failure and must not enter any denominator.
    """
    import time

    from . import ecallisto_catalog as ec

    h0, h1 = window_hours
    out = []
    for date in dates:
        try:
            listed = [
                (s, f) for (s, f) in ec.list_day_files(date) if h0 <= int(f.split("_")[2][:2]) < h1
            ]
        except Exception:
            listed = None
        out.append(
            {
                "date": date,
                "n_files_listed": len(listed) if listed is not None else None,
                "n_stations_listed": len({s for s, _ in listed}) if listed is not None else None,
            }
        )
        if pause:
            time.sleep(pause)
    return out


def run(
    out: str = ".",
    *,
    offline: bool = True,
    dates: list[str] | None = None,
    audit: bool = False,
) -> dict:
    """Full slice: the deterministic validation suite, plus (optionally) a real leg or its audit.

    The validation suite is computed on EVERY invocation (it is deterministic), so the synthetic
    namespace is always populated and there is no mode-dependent macro left to clobber. ``dates``
    runs the heavy real ingest with per-day provenance; ``audit=True`` re-lists the day indexes for
    the days already committed in ``results/ecallisto_census_realdays.csv`` (one cheap GET per day,
    no spectra) and rewrites that CSV with listed-vs-fetched columns plus the real metrics file ---
    the check that turned "168 sampled days" into "45 ingested + 123 silent failures".
    """
    import csv as _csv
    from pathlib import Path

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    syn = validation_suite()
    metrics: dict = {"source": "synthetic", **syn}
    write_results(metrics, op / "results" / "ecallisto_census_metrics.json")

    real: dict | None = None
    if dates is not None and not offline:  # pragma: no cover - network
        ss = fetch_sunspots()
        samples = sample_real_days(dates)
        _write_realdays_csv(op / "results" / "ecallisto_census_realdays.csv", samples)
        real = _real_summary(samples, dt_tol_s=60.0)
        real["source"] = f"e-Callisto x SILSO ({len(dates)} days attempted)"
        real["n_silso_months"] = int(ss["sunspot"].size)
        write_results(real, op / "results" / "ecallisto_census_real_metrics.json")
    elif audit:  # pragma: no cover - network
        path = op / "results" / "ecallisto_census_realdays.csv"
        with path.open() as fh:
            committed = list(_csv.DictReader(fh))
        auditrows = audit_day_listings([r["date"] for r in committed])
        by_date = {a["date"]: a for a in auditrows}
        merged = []
        for r in committed:
            a = by_date.get(r["date"], {})
            cov = float(r["coverage"])
            listed = a.get("n_stations_listed")
            ingest_failed = bool(listed) and cov == 0
            merged.append(
                {
                    "date": r["date"],
                    "n_events": int(r["n_events"]),
                    "coverage": float("nan") if ingest_failed else cov,
                    "n_files_listed": a.get("n_files_listed"),
                    "n_stations_listed": listed,
                    "n_files_fetched": None,
                }
            )
        _write_realdays_csv(path, merged)
        real = _real_summary(merged, dt_tol_s=120.0)
        real["source"] = (
            f"e-Callisto x SILSO ({len(merged)} days attempted; ingest audited "
            "2026-08 against the live day indexes)"
        )
        write_results(real, op / "results" / "ecallisto_census_real_metrics.json")

    if real is not None:  # pragma: no cover - network
        metrics.update({f"real_{k}": v for k, v in real.items() if k != "source"})
        metrics["real_source"] = real["source"]

    # the figure and macros always describe the deterministic validation
    sunspot = synthetic_sunspots()
    n_events, coverage = synthetic_census(sunspot)
    rate = coverage_corrected_rate(n_events, coverage)
    _figure(
        np.arange(sunspot.size, dtype=float),
        rate,
        sunspot,
        op / "papers" / "ecallisto_census" / "figures",
    )
    _write_macros(metrics, op / "papers" / "ecallisto_census" / "generated" / "macros.tex")
    return metrics


def _real_summary(rows: list[dict], *, dt_tol_s: float) -> dict:  # pragma: no cover - network
    """Counts-only summary of the real leg: no correlation coefficient on a handful of events."""

    def _cov(r: dict) -> float:
        try:
            return float(r["coverage"])
        except (TypeError, ValueError):
            return float("nan")

    ing = [r for r in rows if np.isfinite(_cov(r)) and _cov(r) > 0]
    failed = [
        r
        for r in rows
        if (not np.isfinite(_cov(r)) or _cov(r) == 0) and (r.get("n_stations_listed") or 0) > 0
    ]
    dates_ing = sorted(r["date"] for r in ing)
    return {
        "n_days_attempted": len(rows),
        "n_days_ingested": len(ing),
        "n_days_ingest_failed_with_data_listed": len(failed),
        "ingested_span": [dates_ing[0], dates_ing[-1]] if dates_ing else None,
        "n_events_total": int(sum(int(r["n_events"]) for r in ing)),
        "n_days_with_events": int(sum(int(r["n_events"]) > 0 for r in ing)),
        "dt_tol_s": dt_tol_s,
    }


def _write_realdays_csv(path, rows: list[dict]) -> None:  # pragma: no cover - network
    import csv
    from pathlib import Path

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "date",
        "n_events",
        "coverage",
        "n_files_listed",
        "n_stations_listed",
        "n_files_fetched",
    ]
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            row = {c: r.get(c) for c in cols}
            cov = row["coverage"]
            if isinstance(cov, float) and not np.isfinite(cov):
                row["coverage"] = ""  # ingest failure: blank, never 0
            w.writerow(row)


def _figure(x: np.ndarray, rate: np.ndarray, sunspot: np.ndarray, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))

    # Left: the burst rate and the sunspot number over time (twin axes)
    ax1.plot(x, rate, "-", color="C3", lw=1.2, label="type III rate")
    ax1.set(xlabel="time", ylabel="burst rate (events / station)")
    axb = ax1.twinx()
    axb.plot(x, sunspot, "-", color="0.5", lw=1.0, label="sunspot number")
    axb.set_ylabel("sunspot number")
    ax1.set_title("Occurrence vs solar cycle")

    # Right: rate vs sunspot scatter with the OLS fit
    good = np.isfinite(rate) & np.isfinite(sunspot)
    ax2.plot(sunspot[good], rate[good], "o", color="C0", ms=4)
    if good.sum() >= 2 and np.ptp(sunspot[good]) > 0:
        m, b = np.polyfit(sunspot[good], rate[good], 1)
        xs = np.linspace(sunspot[good].min(), sunspot[good].max(), 20)
        ax2.plot(xs, m * xs + b, "-", color="C3", lw=1)
    ax2.set(
        xlabel="sunspot number", ylabel="burst rate (events / station)", title="Rate vs activity"
    )
    fig.tight_layout()
    fig.savefig(out / "census.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lin = m.get("linear") or {}
    grow = m.get("growth") or {}
    sat = m.get("growth_saturating") or {}
    ens = m.get("slope_ensemble") or {}
    span = m.get("real_ingested_span") or [None, None]

    def _g(dic: dict, key: str) -> str:
        val = dic.get(key)
        return "--" if val is None else str(val)

    def _mon(s: str | None) -> str:
        return "--" if not s else f"{s[:4]}-{s[4:6]}"

    # The validation suite is deterministic and recomputed on every run, so the \ecsSyn*
    # namespace can never go stale or be clobbered by a mode change; \ecsReal* carries the
    # audited real-ingest counts (no correlation coefficient -- five events cannot support one).
    lines = [
        "% Auto-generated by jansky_research.ecallisto_census._write_macros -- do not edit by hand.",
        rf"\newcommand{{\ecsSource}}{{{_fmt('source')}}}",
        rf"\newcommand{{\ecsSynKtrue}}{{{_fmt('k_true')}}}",
        rf"\newcommand{{\ecsSynNperiods}}{{{_g(lin, 'n_periods')}}}",
        rf"\newcommand{{\ecsSynNevents}}{{{_g(lin, 'n_events_total')}}}",
        rf"\newcommand{{\ecsSynPearson}}{{{_g(lin, 'pearson_r')}}}",
        rf"\newcommand{{\ecsSynSpearman}}{{{_g(lin, 'spearman_rho')}}}",
        rf"\newcommand{{\ecsSynSlope}}{{{_g(lin, 'slope')}}}",
        rf"\newcommand{{\ecsSynRawPearson}}{{{_g(lin, 'raw_pearson_r')}}}",
        rf"\newcommand{{\ecsSynCovCorr}}{{{_g(lin, 'rate_coverage_corr')}}}",
        rf"\newcommand{{\ecsSynSlopeMean}}{{{_g(ens, 'mean')}}}",
        rf"\newcommand{{\ecsSynSlopeSD}}{{{_g(ens, 'sd')}}}",
        rf"\newcommand{{\ecsSynGrowthRawPearson}}{{{_g(grow, 'raw_pearson_r')}}}",
        rf"\newcommand{{\ecsSynGrowthPearson}}{{{_g(grow, 'pearson_r')}}}",
        rf"\newcommand{{\ecsSynSatPearson}}{{{_g(sat, 'pearson_r')}}}",
        rf"\newcommand{{\ecsSynSatCovCorr}}{{{_g(sat, 'rate_coverage_corr')}}}",
        rf"\newcommand{{\ecsRealSource}}{{{_fmt('real_source')}}}",
        rf"\newcommand{{\ecsRealNdaysAttempted}}{{{_fmt('real_n_days_attempted')}}}",
        rf"\newcommand{{\ecsRealNdaysIngested}}{{{_fmt('real_n_days_ingested')}}}",
        rf"\newcommand{{\ecsRealNdaysFailed}}{{{_fmt('real_n_days_ingest_failed_with_data_listed')}}}",
        rf"\newcommand{{\ecsRealNevents}}{{{_fmt('real_n_events_total')}}}",
        rf"\newcommand{{\ecsRealNdaysWithEvents}}{{{_fmt('real_n_days_with_events')}}}",
        rf"\newcommand{{\ecsRealSpanStart}}{{{_mon(span[0])}}}",
        rf"\newcommand{{\ecsRealSpanEnd}}{{{_mon(span[1])}}}",
        rf"\newcommand{{\ecsRealTolS}}{{{_fmt('real_dt_tol_s')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and
    # would otherwise blank the other mode's macros with '--'. `make figures`
    # runs every slice offline in the repo root, so without this an offline
    # rebuild silently empties this paper. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(
        description="e-Callisto type III occurrence census vs the solar cycle."
    )
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--dates", help="comma-separated YYYYMMDD days for the real run")
    p.add_argument(
        "--audit",
        action="store_true",
        help="re-list the day indexes for the committed realdays CSV (no spectra downloads)",
    )
    args = p.parse_args(argv)
    dates = args.dates.split(",") if args.dates else None
    metrics = run(args.out, offline=args.offline or not dates, dates=dates, audit=args.audit)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
