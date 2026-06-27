# Plan 05 — Airflow on Podman automation 📋

> Context: the automation layer. Depends on 03 (needs pipeline.py); parallel with 04. Scope: medium.

## Context

Automate the analysis with Apache Airflow on local Podman/Docker infrastructure. The design
favours the simplest robust option for a single-machine research pipeline — and crucially the
science must remain reproducible **without** Airflow, so the DAG is a demonstrable automation
layer, not a critical-path dependency.

## Deliverables

- `airflow/compose.yaml` — trimmed official Airflow stack: **postgres + airflow-init + scheduler +
  webserver**. **LocalExecutor** (no Redis/Celery).
- `airflow/Dockerfile` — extends the base image; installs `jansky_research` (and thus `jansky`).
- `airflow/.env.example` — `AIRFLOW_UID` etc., with the rootless-podman note.
- `airflow/dags/research_pipeline.py` — `fetch_dataset → run_analysis → make_figures →
  make_tables → assemble_paper_inputs` (`schedule=None`, `catchup=False`), all **in-process
  PythonOperators** calling `jansky_research.pipeline` — not DockerOperator/PodmanOperator.
  `assemble_paper_inputs` emits `paper/generated/macros.tex`.
- `make airflow-up` / `airflow-down` / `dag-test`.

## Approach — Podman/Fedora specifics (must-haves)

- **SELinux `:z` relabel suffix on every bind mount** (the #1 docker→podman trap on this Fedora
  host).
- **`AIRFLOW_UID=0` for rootless podman** (or `userns_mode: keep-id`) to avoid root-owned `logs/`.
- Idempotent `airflow-init` + DB-connection retry in scheduler/webserver (podman-compose
  `depends_on: service_healthy` lags docker).
- In-process operators import `jansky_research.pipeline` — same code path as `make pipeline` and
  the notebooks (no socket, no nested containers).

## Verification

- `make dag-test` (`airflow dags test research_pipeline 2026-01-01`) runs green under podman;
  artifacts land in `paper/figures/` + `results/metrics.json`.
- `make pipeline` reproduces **identical** artifacts without Airflow.
