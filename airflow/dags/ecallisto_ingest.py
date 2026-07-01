"""Airflow DAG: daily ingestion of the e-Callisto archive into a rolling type III burst catalogue.

This is the automation layer that *uses Airflow for what Airflow is for* --- a **frequently-updated**
archive (e-Callisto adds new dynamic spectra from 150+ stations every day, 20+ years deep), ingested on
a **daily schedule** with **catchup/backfill** over any date range and a **per-station fan-out** via
dynamic task mapping. (Contrast the static-catalogue slices, which a server-less file-DAG runner drives
instead --- see ``workflow/`` / ``plans/32``.)

Shape:  wait_for_day (sensor) -> list_stations -> scan_station.expand(...) -> reduce_day

Every task is an **in-process** call into ``jansky_research.ecallisto_catalog`` --- the same code path as
``make ecallisto-day`` --- not a nested container. Outputs land under ``$JR_OUTPUT`` (a bind-mounted
repo checkout) as one CSV per day, written idempotently so re-running a date overwrites that day's rows.
CI never runs this stack (too heavy/flaky); it is exercised by ``make dag-test``, and the worker it
calls is unit-tested offline.
"""

from __future__ import annotations

import os
from pathlib import Path

import pendulum
from airflow.decorators import dag, task
from airflow.sensors.base import PokeReturnValue

OUT = Path(os.environ.get("JR_OUTPUT", "/opt/airflow/project"))
# cap the per-day fan-out so a backfill demo stays bounded; remove for a full ingest
MAX_FILES = int(os.environ.get("EC_MAX_FILES", "12"))


@dag(
    dag_id="ecallisto_ingest",
    schedule="@daily",
    catchup=True,
    start_date=pendulum.datetime(2011, 9, 14, tz="UTC"),  # the solarbursts recover-a-known day
    max_active_runs=2,
    tags=["jansky-research", "solar", "streaming"],
    doc_md=__doc__,
)
def ecallisto_ingest():
    @task.sensor(poke_interval=120, timeout=3600, mode="reschedule")
    def wait_for_day(ds: str) -> PokeReturnValue:
        """Poke the archive until the day's directory has files --- the frequently-updated hook."""
        from jansky_research import ecallisto_catalog

        files = ecallisto_catalog.list_day_files(ds.replace("-", ""))
        return PokeReturnValue(is_done=len(files) > 0, xcom_value=ds.replace("-", ""))

    @task()
    def list_stations(date: str) -> list[dict]:
        """List the day's station files (bounded by MAX_FILES for the backfill demo)."""
        from jansky_research import ecallisto_catalog

        files = ecallisto_catalog.list_day_files(date)[:MAX_FILES]
        return [{"date": date, "station": s, "file": f} for s, f in files]

    @task()
    def scan_station(item: dict) -> dict:
        """Fetch + scan one station's spectrum (one mapped task per station --- the fan-out)."""
        from jansky_research import ecallisto_catalog, solarbursts

        hhmmss = item["file"].split("_")[2]
        spec = solarbursts.fetch_ecallisto(item["station"], item["date"], hhmmss[:4])
        row = ecallisto_catalog.scan_spectrum(spec)
        row["station"] = item["station"]
        row["date"] = item["date"]
        # peak time as universal time-of-day (file start + local peak) so coincidence compares one clock
        if row.get("t_peak_s") is not None:
            start = int(hhmmss[:2]) * 3600 + int(hhmmss[2:4]) * 60 + int(hhmmss[4:6])
            row["t_peak_s"] = round(start + row["t_peak_s"], 1)
        return row

    @task()
    def reduce_day(rows: list[dict], date: str) -> str:
        """Idempotently write this day's candidate rows + cross-station-coincident events (one CSV each)."""
        import csv

        from jansky_research import ecallisto_catalog

        out = OUT / "results" / "ecallisto_catalog"
        out.mkdir(parents=True, exist_ok=True)
        path = out / f"{date}.csv"
        cols = ["date", "station", "is_burst", "n_channels", "f_lo_mhz", "f_hi_mhz", "drift_mhz_s", "r2"]
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        # coincidence QC: promote multi-station candidates to confirmed events (idempotent per day)
        events = ecallisto_catalog.coincident_events(rows)
        epath = out / f"{date}_events.csv"
        with epath.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["t_peak_s", "n_stations", "median_drift_mhz_s"])
            w.writeheader()
            w.writerows({k: e[k] for k in ("t_peak_s", "n_stations", "median_drift_mhz_s")} for e in events)
        n_burst = sum(1 for r in rows if r.get("is_burst"))
        return f"{date}: {len(events)} confirmed events from {n_burst}/{len(rows)} candidates -> {path}"

    date = wait_for_day()
    items = list_stations(date)
    rows = scan_station.expand(item=items)  # dynamic task mapping: one task per station
    reduce_day(rows, date)


ecallisto_ingest()
