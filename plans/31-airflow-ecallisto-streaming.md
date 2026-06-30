# 31 — Repoint Airflow at a frequently-updated archive: an e-Callisto daily burst-ingest pipeline

Status: ✅ done (worker + streaming DAG + paper + frbstats Airflow claim dropped; CI tests the worker, not the stack)

## Context

The frbstats paper's Airflow-on-Podman section is the project's weakest publishable claim: the dataset
is a **single ~216 KB static CSV** (CHIME/FRB Catalog 1, published once in 2021), and the DAG is
`schedule=None`, `catchup=False`, four in-process `PythonOperator`s in a straight line — a same-code-path
convenience wrapper that exercises **none** of Airflow's strengths (scheduling, backfill, sensors,
parallel fan-out, retries). A referee reads it as over-engineering a CSV, which hurts the submission.

Airflow earns its keep on **large and/or frequently-updated** archives. e-Callisto (the `solarbursts`
data source) is exactly that: 150+ ground stations, **new FITS every day**, 20+ years of archive
(`http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/YYYY/MM/DD/`,
`STATION_YYYYMMDD_HHMMSS_NN.fit.gz`). A scheduled DAG that ingests each day's spectra across stations,
detects type II/III bursts with the tested `solarbursts` tooling, and maintains a rolling burst
catalogue is a genuine Airflow use case — and an honest **astro-ph.IM** automation contribution.

This plan replaces the frbstats Airflow demo with a real streaming pipeline (and pairs with plan #32,
which gives the *static* slices a right-sized orchestrator so Airflow is reserved for streaming/large).

## Deliverables

- `airflow/dags/ecallisto_ingest.py` — a DAG that genuinely uses Airflow:
  - **daily `schedule`** with **`catchup=True`** so a date range backfills (the 20-year archive is the
    backfill story); `max_active_runs` bounded.
  - a **`@task.expand` (dynamic task mapping) fan-out over stations** for a given day — each station's
    file fetched, RFI-flagged, and burst-scanned in parallel (`solarbursts.find_burst_window` /
    `detect_burst_ridge` / drift fit), the per-day results reduced into the rolling catalogue.
  - an **`HttpSensor`/poke** (or a listing check) that waits for a day's directory to populate before
    the fan-out — the "frequently-updated" hook.
  - idempotent writes (re-running a day overwrites that day's rows, not appends), and Podman/Fedora
    gotchas already solved in the repo (fully-qualified images, SELinux `:z`, rootless `AIRFLOW_UID`).
- `src/jansky_research/ecallisto_catalog.py` — the in-process worker the DAG calls: list a day's
  station files, fetch+parse (reuse `solarbursts.fetch_ecallisto`/parser), flag bursts, return rows;
  a `synthetic_day` fixture for offline tests; a `make ecallisto-day DATE=...` CLI that runs one day
  *without* Airflow (the shared code path).
- `tests/test_ecallisto_catalog.py` — offline synthetic-day tests; 85% floor. **CI never runs the
  Airflow stack** (too heavy/flaky) — it tests the worker; the DAG is exercised by `make dag-test`.
- `papers/ecallisto_pipeline/` *(or fold into a revised automation paper)* — astro-ph.IM: "a
  rootless-Podman Airflow pipeline for continuous ingestion of the e-Callisto archive into a rolling
  solar-radio-burst catalogue," with backfill, per-station fan-out, and idempotency as the contribution.
- Update `papers/frbstats/` to **drop the Airflow over-claim** (frbstats becomes a clean tool +
  validation for JOSS/RNAAS; the automation story moves here). Update `README`/`TODO` accordingly.

## Approach

1. **Worker first (CI-testable).** `ecallisto_catalog` over a synthetic day, then a real single day
   (`make ecallisto-day DATE=20110914` reusing the `solarbursts` 2011-09-14 event) — the same code the
   DAG calls. 85% floor on the synthetic.
2. **DAG that uses Airflow's strengths.** Daily schedule + `catchup` backfill over a bounded date range;
   dynamic task mapping fan-out over stations; a sensor on the day's directory; idempotent reduce into
   the catalogue. Demonstrate a **multi-day backfill** and a **parallel station fan-out** under Podman.
3. **Honest scope.** The bursts are *detected candidates* (e-Callisto is uncalibrated, RFI-heavy,
   station-heterogeneous), not a vetted occurrence census; the contribution is the **automation pattern
   on a frequently-updated archive**, not a new burst catalogue. State this plainly.
4. **Write-up / paper surgery.** Either a short standalone astro-ph.IM paper, or revise the frbstats
   paper to remove Airflow and add a new automation paper here. Decide at GATE-2.

## Verification

- `make test` / `cov` green on the synthetic day (no network, no Airflow), 85% floor; `ruff` + `mypy`.
- `make ecallisto-day DATE=...` reproduces a day's rows without Airflow (shared code path).
- `make dag-test` runs `airflow dags test ecallisto_ingest <date>` green under rootless Podman; a
  **multi-day backfill** and the **per-station fan-out** both execute; artifacts land idempotently.
- The frbstats paper no longer claims an Airflow contribution.

## Risks & mitigations

- **Airflow-under-Podman friction (highest) →** keep the existing LocalExecutor+Postgres, in-process
  operators, `:z` relabels, rootless UID; the `make ecallisto-day` path ships the science even if the
  stack misbehaves; CI tests the worker, not Airflow.
- **e-Callisto heterogeneity/RFI →** detection candidates only, honest framing; per-station flags;
  start with a few reliable stations before the full fan-out.
- **Backfill volume →** bound the date range and `max_active_runs`; the *point* is to demonstrate
  catchup/fan-out, not to ingest 20 years in CI.
- **Scope overlap with `solarbursts` →** none in physics (reuses its tooling); the new thing is the
  streaming orchestration.
