# Reproducing this work

Everything regenerates from a clean checkout. Three ways to run the analysis (all the same code
path — `jansky_research.pipeline`), plus the one-command full chain.

## Prerequisites

- [`uv`](https://docs.astral.sh/uv/) and the **`jansky`** repo checked out next to this one
  (`../jansky`) for local development — `pyproject.toml` uses a path source. (CI and the containers
  install `jansky` from its pinned git tag `v0.1.0` instead.)
- For the Airflow and paper steps: **Podman** (rootless is fine). On Fedora there is often no
  `podman compose` provider, so use `podman-compose`; no install needed if you have uv:
  `make ... COMPOSE="uvx podman-compose"`.

## One command

```bash
make reproduce        # fetch data -> run the pipeline -> build the PDF
```

This fetches the CHIME/FRB Catalog 1 CSV, runs the analysis, writes `results/metrics.json` +
`paper/figures/*.pdf` + `paper/generated/macros.tex`, and compiles `paper/main.pdf` in the tectonic
container. Every number in the paper is `\input` from `macros.tex`, so the PDF always matches the
pipeline.

## Step by step

```bash
uv sync                                   # env + jansky
make test                                 # 24 tests, offline, on synthetic fixtures
make pipeline                             # real-data run -> results/ + paper/ artifacts
#   (add ARGS="--offline" to use the synthetic fixture instead)
make paper COMPOSE="uvx podman-compose"   # tectonic -> paper/main.pdf
```

## Via the Airflow automation layer

The same pipeline, orchestrated. Identical artifacts to `make pipeline`.

```bash
cp airflow/.env.example airflow/.env      # AIRFLOW_UID=0 (rootless), JR_OFFLINE
export COMPOSE="uvx podman-compose"
make airflow-up                           # build + start postgres/scheduler/webserver
make dag-test                             # run research_pipeline once -> success
make airflow-down
```

Web UI at http://localhost:8080 (admin/admin) once `airflow-up` is running.

## Notes

- `data/`, `results/`, `paper/figures/`, `paper/generated/`, and `paper/main.pdf` are generated and
  git-ignored — they are products of the pipeline, never committed by hand.
- Offline: with no network, the pipeline falls back to a synthetic catalogue with known ground
  truth (so tests and CI never depend on the network); pass `--offline` to force it.
- The non-chosen survey candidates are preserved in `survey/candidate-gaps.md` as future work.
