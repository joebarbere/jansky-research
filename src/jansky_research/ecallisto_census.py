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
    "census_correlation",
    "coverage_corrected_rate",
    "fetch_sunspots",
    "parse_silso",
    "run",
    "synthetic_census",
    "synthetic_sunspots",
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


def synthetic_sunspots(n_months: int = 180, *, amplitude: float = 140.0, seed: int = 0) -> np.ndarray:
    """A synthetic monthly sunspot-number series: a ~11-year cycle with realistic asymmetry + noise.

    Models the sunspot number as a raised, fast-rise/slow-decay cycle (period 132 months) so the
    synthetic census has a realistic activity curve to track. Returns non-negative monthly values.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(n_months, dtype=float)
    phase = (t % 132.0) / 132.0
    cycle = np.exp(-3.0 * phase) * (1.0 - np.cos(2.0 * np.pi * phase)) / 2.0  # fast rise, slow decay
    base = amplitude * cycle / cycle.max()
    return np.clip(base + rng.normal(0.0, 5.0, n_months), 0.0, None)


def synthetic_census(
    sunspot: np.ndarray,
    *,
    k: float = 0.03,
    coverage_mean: float = 12.0,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Synthetic monthly event counts whose true rate is proportional to the sunspot number.

    For each month with sunspot $S$ and (varying) active-station coverage $C$, the true confirmed-event
    count is Poisson with mean $k\,S\,C$, so the coverage-corrected rate $N/C$ has expectation $k\,S$ ---
    proportional to activity. Returns ``(n_events, coverage)`` per month, from which
    :func:`coverage_corrected_rate` + :func:`census_correlation` recover $r\approx1$ and slope $\approx k$.
    """
    rng = np.random.default_rng(seed)
    s = np.asarray(sunspot, float)
    coverage = np.clip(rng.normal(coverage_mean, coverage_mean * 0.25, s.size), 2.0, None)
    lam = k * s * coverage
    n_events = rng.poisson(np.clip(lam, 0.0, None)).astype(float)
    return n_events, coverage


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
    dates: list[str], *, window_hours: tuple[int, int] = (9, 13)
) -> list[dict]:  # pragma: no cover - network
    """Scan a fixed UT window of each e-Callisto day → per-day confirmed events + station coverage.

    For a consistent occurrence-rate proxy, only files whose start hour is in ``window_hours`` are scanned
    (a fixed sunlit window, sampled identically every day). Returns one row per date with the
    coincidence-confirmed event count and the number of active stations.
    """
    from . import ecallisto_catalog as ec
    from . import solarbursts

    rows_out = []
    h0, h1 = window_hours
    for date in dates:
        files = [
            (s, f)
            for (s, f) in ec.list_day_files(date)
            if h0 <= int(f.split("_")[2][:2]) < h1
        ]
        rows = []
        for station, fname in files:
            hhmmss = fname.split("_")[2]
            try:
                spec = solarbursts.fetch_ecallisto(station, date, hhmmss[:4])
            except Exception:
                continue
            r = ec.scan_spectrum(spec)
            r["station"] = station
            if r.get("t_peak_s") is not None:
                start = int(hhmmss[:2]) * 3600 + int(hhmmss[2:4]) * 60 + int(hhmmss[4:6])
                r["t_peak_s"] = round(start + r["t_peak_s"], 1)
            rows.append(r)
        events = ec.coincident_events(rows, dt_tol_s=120.0)
        rows_out.append(
            {"date": date, "n_events": len(events), "coverage": len({r["station"] for r in rows})}
        )
    return rows_out


def run(out: str = ".", *, offline: bool = True, dates: list[str] | None = None) -> dict:
    """Full slice: build the type III occurrence census and correlate the rate with the sunspot number."""
    import json
    from pathlib import Path

    if offline or dates is None:
        sunspot = synthetic_sunspots()
        n_events, coverage = synthetic_census(sunspot)
        rate = coverage_corrected_rate(n_events, coverage)
        source = "synthetic"
        xlabel = np.arange(sunspot.size, dtype=float)
    else:  # pragma: no cover - network
        ss = fetch_sunspots()
        samples = sample_real_days(dates)
        # pair each sampled day with its month's sunspot number
        n_events = np.array([s["n_events"] for s in samples], float)
        coverage = np.array([s["coverage"] for s in samples], float)
        yrs = np.array([int(d[:4]) + (int(d[4:6]) - 0.5) / 12.0 for d in dates], float)
        sunspot = np.interp(yrs, ss["decimal_year"], ss["sunspot"])
        rate = coverage_corrected_rate(n_events, coverage)
        source = f"e-Callisto x SILSO ({len(dates)} days)"
        xlabel = yrs

    corr = census_correlation(rate, sunspot)
    metrics = {
        "source": source,
        "n_periods": corr["n_periods"],
        "n_events_total": int(np.nansum(n_events)),
        "pearson_r": round(corr["pearson_r"], 3) if np.isfinite(corr["pearson_r"]) else None,
        "spearman_rho": round(corr["spearman_rho"], 3)
        if np.isfinite(corr["spearman_rho"])
        else None,
        "slope": round(corr["slope"], 4) if np.isfinite(corr["slope"]) else None,
    }

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "ecallisto_census_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(xlabel, rate, sunspot, op / "papers" / "ecallisto_census" / "figures")
    _write_macros(metrics, op / "papers" / "ecallisto_census" / "generated" / "macros.tex")
    return metrics


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
    ax2.set(xlabel="sunspot number", ylabel="burst rate (events / station)", title="Rate vs activity")
    fig.tight_layout()
    fig.savefig(out / "census.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.ecallisto_census._write_macros -- do not edit by hand.",
        rf"\newcommand{{\ecsSource}}{{{m['source']}}}",
        rf"\newcommand{{\ecsNperiods}}{{{_fmt('n_periods')}}}",
        rf"\newcommand{{\ecsNevents}}{{{_fmt('n_events_total')}}}",
        rf"\newcommand{{\ecsPearson}}{{{_fmt('pearson_r')}}}",
        rf"\newcommand{{\ecsSpearman}}{{{_fmt('spearman_rho')}}}",
        rf"\newcommand{{\ecsSlope}}{{{_fmt('slope')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="e-Callisto type III occurrence census vs the solar cycle.")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--dates", help="comma-separated YYYYMMDD days for the real run")
    args = p.parse_args(argv)
    dates = args.dates.split(",") if args.dates else None
    metrics = run(args.out, offline=args.offline or not dates, dates=dates)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
