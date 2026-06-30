"""e-Callisto daily burst-ingest worker: scan a day of solar-radio spectra for type III candidates.

This is the in-process worker behind the e-Callisto Airflow pipeline (``airflow/dags/ecallisto_ingest``)
and the ``make ecallisto-day`` CLI --- the *same* code path, so the DAG and the command line produce
identical rows. e-Callisto is a frequently-updated archive (150+ stations, new gzipped-FITS dynamic
spectra every day, 20+ years), which is what makes a *scheduled, backfilling, per-station fan-out*
pipeline the right shape for it (unlike a one-shot static catalogue).

Each station file is scanned for a drifting type III ridge by reusing the tested ``solarbursts``
dynamic-spectrum tools (background subtraction, burst windowing, the per-channel ridge detector, the
robust drift fit). The output is a per-day table of **burst candidates** --- e-Callisto is uncalibrated
and RFI-heavy, so these are detections to be vetted, not a finished occurrence census. Pure NumPy with a
synthetic-day offline fixture; the real listing/fetch is network-gated.
"""

from __future__ import annotations

import numpy as np

from . import solarbursts

__all__ = [
    "ingest_day",
    "list_day_files",
    "run",
    "scan_day_specs",
    "scan_spectrum",
    "synthetic_day",
]

ECALLISTO_BASE = "http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto"


def scan_spectrum(
    spec: dict,
    *,
    min_channels: int = 20,
    r2_min: float = 0.5,
    pad_s: float = 10.0,
) -> dict:
    """Scan one dynamic spectrum for a type III drift ridge; return a candidate row.

    Reuses :func:`solarbursts.find_burst_window` / :func:`solarbursts.detect_burst_ridge` /
    :func:`solarbursts.fit_drift_rate`. A candidate is flagged (``is_burst``) when the ridge spans at
    least ``min_channels`` channels, drifts the right way (frequency falling with time, ``drift < 0``),
    and the ridge fit is coherent (``r2 >= r2_min``) --- the discriminant the ``solarbursts`` slice
    found trustworthy. Returns the channel count, frequency span, drift rate (MHz/s), fit R^2, peak
    time, and the boolean flag.
    """
    data, freqs, times = spec["data"], spec["freqs"], spec["times"]
    window = solarbursts.find_burst_window(data, times, pad_s=pad_s)
    rf, rt = solarbursts.detect_burst_ridge(data, freqs, times, window=window)
    nan = float("nan")
    if rf.size < 2:
        return {"n_channels": int(rf.size), "is_burst": False, "drift_mhz_s": nan, "r2": nan}
    drift = solarbursts.fit_drift_rate(rf, rt)
    slope, icpt, keep = solarbursts._robust_linfit(rt, rf)
    fk, tk = rf[keep], rt[keep]
    ss_res = float(np.sum((fk - (slope * tk + icpt)) ** 2))
    ss_tot = float(np.sum((fk - np.mean(fk)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else nan
    is_burst = bool(rf.size >= min_channels and drift < 0 and np.isfinite(r2) and r2 >= r2_min)
    return {
        "n_channels": int(rf.size),
        "f_lo_mhz": round(float(np.min(rf)), 2),
        "f_hi_mhz": round(float(np.max(rf)), 2),
        "drift_mhz_s": round(drift, 3) if np.isfinite(drift) else None,
        "r2": round(r2, 3) if np.isfinite(r2) else None,
        "t_peak_s": round(float(np.median(rt)), 1),
        "is_burst": is_burst,
    }


def scan_day_specs(specs: list[tuple[str, dict]], **kw) -> list[dict]:
    """Scan a day's (station, spectrum) pairs; return one candidate row per station."""
    rows = []
    for station, spec in specs:
        row = scan_spectrum(spec, **kw)
        row["station"] = station
        rows.append(row)
    return rows


def synthetic_day(
    *, n_stations: int = 8, n_bursts: int = 3, seed: int = 0
) -> list[tuple[str, dict]]:
    """A synthetic observing day: ``n_bursts`` stations with an injected type III, the rest quiet.

    Each "burst" station gets a :func:`solarbursts.synthetic_burst` spectrum; the quiet stations get a
    same-shaped pure-noise spectrum. Returns (station, spectrum) pairs so :func:`scan_day_specs`
    recovers exactly the injected bursts offline (no network, no large files).
    """
    rng = np.random.default_rng(seed)
    template = solarbursts.synthetic_burst(seed=seed)
    shape = template["data"].shape
    specs: list[tuple[str, dict]] = []
    for i in range(n_stations):
        name = f"STATION{i:02d}"
        if i < n_bursts:
            spec = solarbursts.synthetic_burst(seed=seed + i + 1)
        else:
            spec = {
                "data": rng.normal(0.0, 1.0, shape),
                "freqs": template["freqs"],
                "times": template["times"],
            }
        specs.append((name, spec))
    return specs


def list_day_files(date_yyyymmdd: str) -> list[tuple[str, str]]:  # pragma: no cover - network
    """List the e-Callisto archive files for one day → ``(station, filename)`` pairs.

    Parses the public day-directory HTML index for ``<station>_<date>_<hhmmss>_NN.fit.gz`` files.
    """
    import re

    import requests

    yyyy, mm, dd = date_yyyymmdd[:4], date_yyyymmdd[4:6], date_yyyymmdd[6:8]
    day_url = f"{ECALLISTO_BASE}/{yyyy}/{mm}/{dd}/"
    idx = requests.get(day_url, timeout=60).text
    out = []
    for m in re.finditer(rf"([A-Za-z0-9\-]+)_{date_yyyymmdd}_([0-9]{{6}})_[0-9]+\.fit\.gz", idx):
        out.append((m.group(1), m.group(0)))
    return out


def ingest_day(
    date_yyyymmdd: str, *, stations: list[str] | None = None, max_files: int | None = None, **kw
) -> list[dict]:  # pragma: no cover - network
    """Fetch and scan a day's e-Callisto spectra (optionally restricted to ``stations``).

    Lists the day, fetches each gzipped-FITS spectrum (reusing ``solarbursts``' parser via
    :func:`solarbursts.fetch_ecallisto`), scans it, and returns the candidate rows. ``max_files`` caps
    the fan-out for a quick run.
    """
    files = list_day_files(date_yyyymmdd)
    if stations is not None:
        files = [(s, f) for (s, f) in files if s in stations]
    if max_files is not None:
        files = files[:max_files]
    rows = []
    for station, fname in files:
        hhmm = fname.split("_")[2][:4]
        try:
            spec = solarbursts.fetch_ecallisto(station, date_yyyymmdd, hhmm)
        except Exception:
            continue
        row = scan_spectrum(spec, **kw)
        row["station"] = station
        row["file"] = fname
        rows.append(row)
    return rows


def _metrics(rows: list[dict], source: str) -> dict:
    n = len(rows)
    bursts = [r for r in rows if r.get("is_burst")]
    drifts = [r["drift_mhz_s"] for r in bursts if r.get("drift_mhz_s") is not None]
    return {
        "source": source,
        "n_scanned": n,
        "n_bursts": len(bursts),
        "burst_fraction": round(len(bursts) / n, 3) if n else None,
        "median_drift_mhz_s": round(float(np.median(drifts)), 3) if drifts else None,
    }


def run(out: str = ".", *, offline: bool = True, date: str | None = None, **kw) -> dict:
    """Full worker: scan a day (synthetic offline, or the real archive) → catalogue + metrics + figure."""
    import csv
    import json
    from pathlib import Path

    example: dict | None = (
        None  # a representative detected-burst spectrum, for the illustration panel
    )
    if offline or date is None:
        specs = synthetic_day()
        rows = scan_day_specs(specs, **kw)
        source = "synthetic-day"
        by_station = dict(specs)
        example = next((by_station[r["station"]] for r in rows if r.get("is_burst")), None)
    else:  # pragma: no cover - network
        rows = ingest_day(date, **kw)
        source = f"e-Callisto {date}"
        hit = next((r for r in rows if r.get("is_burst")), None)
        if hit is not None:
            example = solarbursts.fetch_ecallisto(
                hit["station"], date, hit["file"].split("_")[2][:4]
            )

    metrics = _metrics(rows, source)
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "ecallisto_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    if rows:
        cols = ["station", "is_burst", "n_channels", "f_lo_mhz", "f_hi_mhz", "drift_mhz_s", "r2"]
        with (op / "results" / "ecallisto_catalog.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
    _figure(rows, example, op / "papers" / "ecallisto_pipeline" / "figures")
    _write_macros(metrics, op / "papers" / "ecallisto_pipeline" / "generated" / "macros.tex")
    return metrics


def _figure(rows: list[dict], example: dict | None, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.6))

    # Left: a representative detected burst -- the dynamic spectrum with the ridge the worker fits
    if example is not None:
        data, freqs, times = example["data"], example["freqs"], example["times"]
        clean = solarbursts.background_subtract(data)
        window = solarbursts.find_burst_window(data, times)
        rf, rt = solarbursts.detect_burst_ridge(data, freqs, times, window=window)
        ax1.pcolormesh(times, freqs, clean, cmap="inferno", shading="auto")
        ax1.plot(rt, rf, ".", color="cyan", ms=2, label="detected ridge")
        ax1.set(xlabel="time (s)", ylabel="frequency (MHz)", title="Example type III detection")
        ax1.legend(loc="upper right", fontsize=8)
    else:
        ax1.set_axis_off()

    # Right: the day's scan summary -- candidates vs quiet stations
    n_burst = sum(1 for r in rows if r.get("is_burst"))
    n_quiet = len(rows) - n_burst
    ax2.bar(["burst\ncandidate", "quiet"], [n_burst, n_quiet], color=["C3", "0.6"])
    ax2.set(ylabel="stations", title="Day scan summary")
    fig.tight_layout()
    fig.savefig(out / "ecallisto.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.ecallisto_catalog._write_macros -- do not edit by hand.",
        rf"\newcommand{{\ecSource}}{{{m['source']}}}",
        rf"\newcommand{{\ecNscanned}}{{{_fmt('n_scanned')}}}",
        rf"\newcommand{{\ecNbursts}}{{{_fmt('n_bursts')}}}",
        rf"\newcommand{{\ecBurstFrac}}{{{_fmt('burst_fraction')}}}",
        rf"\newcommand{{\ecMedDrift}}{{{_fmt('median_drift_mhz_s')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(
        description="Scan a day of e-Callisto spectra for type III candidates."
    )
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--date", help="YYYYMMDD")
    p.add_argument("--max-files", type=int, default=None)
    args = p.parse_args(argv)
    kw = {} if args.max_files is None else {"max_files": args.max_files}
    metrics = run(args.out, offline=args.offline or not args.date, date=args.date, **kw)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
