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

This fetches the public datasets, runs every slice, writes `results/*.json` +
`papers/<slice>/figures/*.pdf` + `papers/<slice>/generated/macros.tex`, compiles each
`papers/<slice>/main.pdf` in the tectonic container, and assembles an arXiv package per paper. Every
number in each paper is `\input` from its `macros.tex`, so the PDFs always match the pipeline.

## Step by step

```bash
uv sync                                   # env + jansky
make test                                 # offline tests, on synthetic fixtures
make figures                              # run every slice -> results/ + papers/<slice>/ artifacts
#   (the `figures` target uses the synthetic fixtures; run a slice directly without --offline for
#    real public data, e.g. `uv run python -m jansky_research.hi`)
make paper                                # tectonic -> every papers/<slice>/main.pdf
make arxiv                                # assemble + validate an arXiv package per paper
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

- `data/`, `results/`, and each `papers/<slice>/`'s `figures/`, `generated/`, `arxiv-submission/`,
  and `main.pdf` are generated and git-ignored — products of the pipeline, never committed by hand
  (only each paper's `main.tex` + `refs.bib` are tracked).
- Offline: with no network, the pipeline falls back to a synthetic catalogue with known ground
  truth (so tests and CI never depend on the network); pass `--offline` to force it.
- The non-chosen survey candidates are preserved in `survey/candidate-gaps.md` as future work.
