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
    "coincident_events",
    "ingest_day",
    "list_day_files",
    "run",
    "scan_day_specs",
    "scan_spectrum",
    "synthetic_coincident_day",
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


def coincident_events(
    rows: list[dict], *, dt_tol_s: float = 60.0, min_stations: int = 2
) -> list[dict]:
    """Cluster per-station burst candidates into cross-station-**coincident** events.

    A real solar radio burst is seen at (near) the same universal time by every station on the sunlit
    side; RFI and local artefacts are single-station. So a candidate confirmed at ``>= min_stations``
    distinct stations within ``dt_tol_s`` of each other is a real burst, while an isolated single-station
    candidate is rejected --- the coincidence QC that turns raw candidates into a trustworthy catalogue.
    Groups the ``is_burst`` rows by peak time (single-linkage in time) and returns the confirmed events:
    mean peak time, the number and list of stations, and the median drift rate.
    """
    bursts = sorted(
        (r for r in rows if r.get("is_burst") and r.get("t_peak_s") is not None),
        key=lambda r: r["t_peak_s"],
    )
    clusters: list[list[dict]] = []
    for r in bursts:
        if clusters and r["t_peak_s"] - clusters[-1][-1]["t_peak_s"] <= dt_tol_s:
            clusters[-1].append(r)
        else:
            clusters.append([r])
    events = []
    for c in clusters:
        stations = sorted({r["station"] for r in c})
        if len(stations) < min_stations:
            continue
        drifts = [r["drift_mhz_s"] for r in c if r.get("drift_mhz_s") is not None]
        events.append(
            {
                "t_peak_s": round(float(np.mean([r["t_peak_s"] for r in c])), 1),
                # single-linkage has no span cap: adjacent members are within dt_tol_s but the
                # CLUSTER can span longer, so the span is reported per event rather than
                # implied to be <= dt_tol_s (round-11 referee: 8 stations 50 s apart -- true
                # span 350 s -- confirmed as one "60 s" event)
                "span_s": round(float(c[-1]["t_peak_s"] - c[0]["t_peak_s"]), 1),
                "n_stations": len(stations),
                "stations": stations,
                "median_drift_mhz_s": round(float(np.median(drifts)), 3) if drifts else None,
            }
        )
    return events


def chance_coincidence_rate(
    candidate_counts: list[int], *, dt_tol_s: float = 60.0, day_s: float = 86400.0
) -> float:
    """Expected number of CHANCE cross-station coincidences per day, in closed form.

    With ``c_i`` independent, uniformly-timed candidates at station ``i``, each cross-station
    pair falls within ``dt_tol_s`` with probability ~``2 dt / T``; summing over pairs gives the
    expected spurious event count. This is the number that sets the coincidence QC's
    reliability -- two chance-coincident interference candidates ARE promoted by construction
    -- and it was previously never quoted.
    """
    total = 0.0
    for i, ci in enumerate(candidate_counts):
        for cj in candidate_counts[i + 1 :]:
            total += ci * cj * (2.0 * dt_tol_s / day_s)
    return total


def synthetic_coincident_day(
    *,
    n_coincident: int = 4,
    n_rfi: int = 3,
    n_quiet: int = 3,
    t_burst_s: float = 300.0,
    include_morph_contaminants: bool = True,
    seed: int = 0,
) -> list[tuple[str, dict]]:
    """A synthetic observing day with one real (multi-station) burst plus single-station RFI.

    ``n_coincident`` stations see the **same** type III at the common time ``t_burst_s`` (real burst;
    independent noise); ``n_rfi`` stations each carry a spurious burst at a *distinct* time well outside
    the coincidence window (local interference); ``n_quiet`` stations are pure noise. The burst UT is set
    by shifting each spectrum's time axis. So :func:`coincident_events` recovers exactly one confirmed
    event (the real burst) and rejects the single-station RFI --- the offline recover-a-known.
    """
    rng = np.random.default_rng(seed)
    template = solarbursts.synthetic_burst(seed=seed)
    specs: list[tuple[str, dict]] = []
    idx = 0
    for i in range(n_coincident):
        b = solarbursts.synthetic_burst(seed=seed + i)
        specs.append((f"STATION{idx:02d}", {**b, "times": b["times"] + t_burst_s}))
        idx += 1
    for i in range(n_rfi):
        b = solarbursts.synthetic_burst(seed=seed + 100 + i)
        t_rfi = (
            t_burst_s + 200.0 + i * 150.0
        )  # distinct, > dt_tol from the burst and from each other
        specs.append((f"STATION{idx:02d}", {**b, "times": b["times"] + t_rfi}))
        idx += 1
    for _ in range(n_quiet):
        noise = rng.normal(0.0, 1.0, template["data"].shape)
        specs.append(
            (
                f"STATION{idx:02d}",
                {"data": noise, "freqs": template["freqs"], "times": template["times"]},
            )
        )
        idx += 1
    if include_morph_contaminants:
        # contaminants the detector must reject on MORPHOLOGY, not timing -- and two of them
        # sit within the 60 s coincidence tolerance of each other, so if morphology failed
        # they WOULD be promoted (the earlier fixture's "RFI" was the same synthetic_burst
        # function as the real burst, so the validation could not fail on this axis)
        shape = template["data"].shape
        # (a) constant narrowband carrier: no drift, few channels
        car = rng.normal(0.0, 1.0, shape)
        car[shape[0] // 3 : shape[0] // 3 + 2, :] += 9.0
        specs.append(
            (
                f"STATION{idx:02d}",
                {"data": car, "freqs": template["freqs"], "times": template["times"] + 430.0},
            )
        )
        idx += 1
        # (b) broadband zero-drift impulse, 30 s from the carrier (INSIDE the tolerance)
        imp = rng.normal(0.0, 1.0, shape)
        imp[:, shape[1] // 2 : shape[1] // 2 + 4] += 10.0
        specs.append(
            (
                f"STATION{idx:02d}",
                {"data": imp, "freqs": template["freqs"], "times": template["times"] + 460.0},
            )
        )
        idx += 1
        # (c) reverse-drift ridge: a type III flipped in frequency drifts the wrong way
        rev = solarbursts.synthetic_burst(seed=seed + 500)
        specs.append(
            (
                f"STATION{idx:02d}",
                {
                    "data": rev["data"][::-1],
                    "freqs": rev["freqs"],
                    "times": rev["times"] + 650.0,
                },
            )
        )
        idx += 1
    return specs


def synthetic_ensemble(n_seeds: int = 60, **kw) -> dict:
    """Seed-ensemble statistics for the synthetic-day validation (one seed is one realization).

    "Recovers exactly the injected bursts" on seed 0 alone is the single-realization pattern
    this repo has been bitten by (rmstructure); this measures the RATE: over ``n_seeds``
    independent days, how often the burst is confirmed as exactly one event at the injected
    multiplicity, how many contaminant/quiet stations are ever false-flagged, and the drift
    spread. Deterministic given ``n_seeds``.
    """
    n_ok = 0
    n_false_flag = 0
    n_contaminant_flag = 0
    drifts = []
    for s in range(n_seeds):
        specs = synthetic_coincident_day(seed=s, **kw)
        rows = scan_day_specs(specs)
        events = coincident_events(rows)
        n_burst_st = sum(1 for r in rows[:4] if r.get("is_burst"))
        quiet_flags = sum(1 for r in rows[7:10] if r.get("is_burst"))
        morph_flags = sum(1 for r in rows[10:] if r.get("is_burst"))
        n_false_flag += quiet_flags
        n_contaminant_flag += morph_flags
        if len(events) == 1 and events[0]["n_stations"] == 4 and n_burst_st == 4:
            n_ok += 1
        for e in events:
            if e.get("median_drift_mhz_s") is not None:
                drifts.append(e["median_drift_mhz_s"])
    return {
        "n_seeds": n_seeds,
        "n_recovered_exactly": n_ok,
        "n_quiet_false_flags": n_false_flag,
        "n_contaminant_false_flags": n_contaminant_flag,
        "drift_mean": round(float(np.mean(drifts)), 3) if drifts else None,
        "drift_sd": round(float(np.std(drifts, ddof=1)), 3) if len(drifts) > 1 else None,
    }


def list_day_files(date_yyyymmdd: str) -> list[tuple[str, str]]:  # pragma: no cover - network
    """List the e-Callisto archive files for one day → ``(station, filename)`` pairs.

    Parses the public day-directory HTML index for ``<station>_<date>_<hhmmss>_NN.fit.gz`` files.
    """
    import re

    import requests

    yyyy, mm, dd = date_yyyymmdd[:4], date_yyyymmdd[4:6], date_yyyymmdd[6:8]
    day_url = f"{ECALLISTO_BASE}/{yyyy}/{mm}/{dd}/"
    idx = requests.get(day_url, timeout=60).text
    seen = set()
    out = []
    # the index HTML renders each filename TWICE (href attribute + link text); without the
    # de-dup every file was double-counted and n_scanned was 2x the spectra (round-11 referee)
    for m in re.finditer(rf"([A-Za-z0-9\-]+)_{date_yyyymmdd}_([0-9]{{6}})_[0-9]+\.fit\.gz", idx):
        if m.group(0) not in seen:
            seen.add(m.group(0))
            out.append((m.group(1), m.group(0)))
    return out


def scan_file(
    station: str, date_yyyymmdd: str, fname: str, **kw
) -> dict:  # pragma: no cover - network
    """Fetch EXACTLY one listed file, scan it, and return the row with a universal-time peak.

    The single shared worker for both the CLI (:func:`ingest_day`) and the Airflow DAG's
    ``scan_station`` task, so the two paths cannot drift (the round-11 referee found the DAG
    reimplementing this inline with different error semantics). Fetching by FILENAME --- not by
    re-resolving the HHMM --- is what guarantees the scanned file is the listed file, and hence
    that the universal-time conversion below uses the right file start.
    """
    spec = solarbursts.fetch_ecallisto(station, date_yyyymmdd, "", filename=fname)
    row = scan_spectrum(spec, **kw)
    row["station"] = station
    row["file"] = fname
    # convert the local (from-file-start) peak time to universal time-of-day so coincidence
    # compares the same clock across stations whose 15-min files begin at different UTs
    if row.get("t_peak_s") is not None:
        hhmmss = fname.split("_")[2]
        start_sod = int(hhmmss[:2]) * 3600 + int(hhmmss[2:4]) * 60 + int(hhmmss[4:6])
        row["t_peak_s"] = round(start_sod + row["t_peak_s"], 1)
    return row


def ingest_day(
    date_yyyymmdd: str, *, stations: list[str] | None = None, max_files: int | None = None, **kw
) -> tuple[list[dict], dict]:  # pragma: no cover - network
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
    n_failed = 0
    for station, fname in files:
        try:
            row = scan_file(station, date_yyyymmdd, fname, **kw)
        except Exception:
            # counted, never silent: a throttled day must be distinguishable from an empty one
            n_failed += 1
            continue
        rows.append(row)
    ingest_stats = {
        "n_files_listed": len(files),
        "n_stations_listed": len({s for s, _ in files}),
        "n_fetch_failed": n_failed,
        "stations_scope": sorted({s for s, _ in files}) if stations is not None else "all",
    }
    return rows, ingest_stats


def _metrics(
    rows: list[dict], events: list[dict], source: str, ingest_stats: dict | None = None
) -> dict:
    n = len(rows)
    bursts = [r for r in rows if r.get("is_burst")]
    drifts = [r["drift_mhz_s"] for r in bursts if r.get("drift_mhz_s") is not None]
    n_confirmed_det = sum(e["n_stations"] for e in events)
    per_station_counts: dict[str, int] = {}
    for r in bursts:
        per_station_counts[r["station"]] = per_station_counts.get(r["station"], 0) + 1
    out = {
        "source": source,
        "n_scanned": n,
        "n_stations_scanned": len({r["station"] for r in rows}),
        "n_bursts": len(bursts),
        "n_candidate_stations": len(per_station_counts),
        "burst_fraction": round(len(bursts) / n, 3) if n else None,
        "median_drift_mhz_s": round(float(np.median(drifts)), 3) if drifts else None,
        # cross-station coincidence QC
        "n_events": len(events),
        "event_spans_s": [e["span_s"] for e in events],
        "max_event_stations": max((e["n_stations"] for e in events), default=0),
        "n_rfi_rejected": len(bursts) - n_confirmed_det,
        # the reliability the QC actually has: the expected CHANCE coincidences per day given
        # the observed per-station candidate counts (two chance-coincident candidates ARE
        # promoted by construction)
        "expected_chance_events_per_day": round(
            chance_coincidence_rate(list(per_station_counts.values())), 4
        ),
    }
    if ingest_stats:
        out.update(ingest_stats)
    return out


def run(out: str = ".", *, offline: bool = True, date: str | None = None, **kw) -> dict:
    """Full worker: scan a day (synthetic offline, or the real archive) → catalogue + metrics + figure."""
    import csv
    from pathlib import Path

    example: dict | None = (
        None  # a representative detected-burst spectrum, for the illustration panel
    )
    if offline or date is None:
        specs = synthetic_coincident_day()
        rows = scan_day_specs(specs, **kw)
        source = "synthetic-day"
        by_station = dict(specs)
        example = next((by_station[r["station"]] for r in rows if r.get("is_burst")), None)
        ingest_stats: dict | None = None
    else:  # pragma: no cover - network
        rows, ingest_stats = ingest_day(date, **kw)
        source = f"e-Callisto {date}"
        hit = next((r for r in rows if r.get("is_burst")), None)
        if hit is not None:
            example = solarbursts.fetch_ecallisto(hit["station"], date, "", filename=hit["file"])

    events = coincident_events(rows)  # cross-station coincidence QC -> confirmed events
    metrics = _metrics(rows, events, source, ingest_stats)
    if offline or date is None:
        # the seed-ensemble rate the single-realization sentence needs (deterministic)
        metrics["ensemble"] = synthetic_ensemble(n_seeds=60)
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    # Each leg keeps its own evidence file. Sharing one name made the synthetic macros
    # unauditable: the paper cites the synthetic day, and only the real day's JSON survived.
    stem = "ecallisto_synthetic" if source.startswith("synthetic") else "ecallisto"
    write_results(metrics, op / "results" / f"{stem}_metrics.json")
    if rows:
        # t_peak_s and file are the two fields coincident_events consumes: a catalogue that
        # omits them cannot audit n_events (the innerrc lesson, found live here in round 11)
        cols = [
            "station",
            "file",
            "t_peak_s",
            "is_burst",
            "n_channels",
            "f_lo_mhz",
            "f_hi_mhz",
            "drift_mhz_s",
            "r2",
        ]
        with (op / "results" / f"{stem}_catalog.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
    # one figure per leg: the committed real leg's output once sat under the synthetic
    # caption because both legs wrote one path (round-11 blocker)
    figname = "ecallisto_syn" if source.startswith("synthetic") else "ecallisto_real"
    _figure(rows, example, events, op / "papers" / "ecallisto_pipeline" / "figures", name=figname)
    _write_macros(metrics, op / "papers" / "ecallisto_pipeline" / "generated" / "macros.tex")
    return metrics


def _figure(
    rows: list[dict], example: dict | None, events: list[dict], out_dir, *, name: str = "ecallisto"
) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))

    # Left: a representative detected burst -- the dynamic spectrum with the ridge the worker fits
    if example is not None:
        data, freqs, times = example["data"], example["freqs"], example["times"]
        clean = solarbursts.background_subtract(data)
        window = solarbursts.find_burst_window(data, times)
        rf, rt = solarbursts.detect_burst_ridge(data, freqs, times, window=window)
        ax1.pcolormesh(times - times[0], freqs, clean, cmap="inferno", shading="auto")
        ax1.plot(rt - times[0], rf, ".", color="cyan", ms=2, label="detected ridge")
        ax1.set(xlabel="time (s)", ylabel="frequency (MHz)", title="Example type III detection")
        ax1.legend(loc="upper right", fontsize=8)
    else:
        ax1.set_axis_off()

    # Right: the coincidence timeline -- each candidate's peak time, confirmed (multi-station) vs single
    confirmed_t = {round(e["t_peak_s"], 1) for e in events}
    bursts = [r for r in rows if r.get("is_burst") and r.get("t_peak_s") is not None]
    stations = sorted({r["station"] for r in bursts})
    ymap = {s: k for k, s in enumerate(stations)}
    for r in bursts:
        near = any(abs(r["t_peak_s"] - t) <= 60.0 for t in confirmed_t)
        ax2.scatter(
            r["t_peak_s"],
            ymap[r["station"]],
            color="C3" if near else "0.6",
            s=40,
            marker="o" if near else "x",
        )
    for e in events:
        ax2.axvline(e["t_peak_s"], color="C3", ls="--", lw=0.6, alpha=0.6)
    ax2.set(
        xlabel="burst peak time (s)",
        ylabel="station",
        yticks=range(len(stations)),
        title=f"Coincidence QC: {len(events)} confirmed",
    )
    ax2.set_yticklabels(stations, fontsize=6)
    fig.tight_layout()
    fig.savefig(out / f"{name}.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    # Every value here is mode-dependent, so every name is namespaced. Un-namespaced, the two
    # legs wrote the same seven names and the last run won: the committed file held the real
    # 2011-09-14 archive day (n_events 0, because no burst was confirmed that day) while all six
    # macro uses in the paper describe the *synthetic* day, so the abstract's only result
    # typeset as "a real burst at 0 stations ... confirms exactly 0 event". See CLAUDE.md,
    # "Merging is not enough -- mode-dependent macros must be NAMESPACED".
    keys = [
        ("Nscanned", "n_scanned"),
        ("NstationsScanned", "n_stations_scanned"),
        ("Nbursts", "n_bursts"),
        ("NcandStations", "n_candidate_stations"),
        ("BurstFrac", "burst_fraction"),
        ("MedDrift", "median_drift_mhz_s"),
        ("Nevents", "n_events"),
        ("MaxEventStations", "max_event_stations"),
        ("NrfiRejected", "n_rfi_rejected"),
        ("NfilesListed", "n_files_listed"),
        ("NfetchFailed", "n_fetch_failed"),
        ("ChanceEvents", "expected_chance_events_per_day"),
    ]
    # per-leg SOURCE macros: a single \ecSource made the downgrade guard see-saw (a synthetic
    # rerun looked like a downgrade of the last real write and could not update its OWN
    # namespace); with namespaced sources every write is a plain placeholder merge
    mode = "Syn" if "synthetic" in str(m["source"]).lower() else "Real"
    lines = [
        "% Auto-generated by jansky_research.ecallisto_catalog._write_macros -- do not edit by hand.",
        rf"\newcommand{{\ecSynSource}}{{{m['source'] if mode == 'Syn' else '--'}}}",
        rf"\newcommand{{\ecRealSource}}{{{m['source'] if mode == 'Real' else '--'}}}",
        r"\newcommand{\ecTolS}{60}",
        r"\newcommand{\ecMinStations}{2}",
    ]
    ens = m.get("ensemble") or {}
    if ens:
        lines += [
            rf"\newcommand{{\ecSynEnsN}}{{{ens['n_seeds']}}}",
            rf"\newcommand{{\ecSynEnsOk}}{{{ens['n_recovered_exactly']}}}",
            rf"\newcommand{{\ecSynEnsQuietFP}}{{{ens['n_quiet_false_flags']}}}",
            rf"\newcommand{{\ecSynEnsContFP}}{{{ens['n_contaminant_false_flags']}}}",
            rf"\newcommand{{\ecSynEnsDriftSd}}{{{ens['drift_sd']}}}",
        ]
    else:
        lines += [
            rf"\newcommand{{\ecSynEns{k}}}{{--}}"
            for k in ("N", "Ok", "QuietFP", "ContFP", "DriftSd")
        ]
    # Both namespaces are emitted on every run: this run fills its own and leaves the other
    # mode's as placeholders, which preserve_live_macros then restores from the committed file.
    # Emitting only one namespace would DELETE the other -- the merge rewrites the lines the new
    # run emits and cannot carry over a name the new text never mentions.
    for suffix, key in keys:
        for ns in ("Syn", "Real"):
            value = _fmt(key) if ns == mode else "--"
            lines.append(rf"\newcommand{{\ec{ns}{suffix}}}{{{value}}}")
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
