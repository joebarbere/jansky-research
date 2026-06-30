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
make figures                              # build every static slice's artifacts via the Snakemake DAG
#   (a server-less file-target DAG, workflow/Snakefile; uses the synthetic fixtures and rebuilds only
#    what changed. Needs the `workflow` extra: `uv sync --extra workflow`, Python >=3.11. Run a slice
#    directly without --offline for real public data, e.g. `uv run python -m jansky_research.hi`)
make figures-dry                          # show the static-slice DAG without running it
make paper                                # tectonic -> every papers/<slice>/main.pdf
make arxiv                                # assemble + validate an arXiv package per paper
```

### VLASS multi-epoch variability (real data)

The VLASS slice queries real survey catalogues over a sky region (Epoch 1 via VizieR TAP; Epochs
2–3 are bulk NRAO files of ~0.3–1 GB each, cached under `data/`):

```bash
uv sync --extra vlass                                         # adds pyvo
uv run python -m jansky_research.vlass --ra 180 --dec 30 --radius 1.0 --epochs 1,2,3
```

It applies the per-epoch Quick-Look flux-scale corrections (VLASS Memos 13/22) before computing
variability, so the epoch-to-epoch scale offset is not mistaken for variability.

## Right-sized orchestration: Snakemake for static slices, Airflow for streaming

The orchestration is matched to the data. The **static** slices (one catalogue / CDF / CSV each, run
once) are a file-target dependency graph, driven by **Snakemake** (`workflow/Snakefile`, server-less,
parallel, rebuilds only what changed) — this is what `make figures` runs. A **frequently-updated**
archive — e-Callisto, new spectra every day across 150+ stations — is an ingestion problem, and that is
**Airflow**'s job: a daily-scheduled, backfilling, per-station fan-out DAG.

```bash
# static slices (Snakemake)
uv sync --extra workflow                  # adds snakemake (Python >=3.11)
make figures                              # = snakemake -j over workflow/Snakefile

# streaming ingest (Airflow on rootless Podman)
cp airflow/.env.example airflow/.env      # AIRFLOW_UID=0 (rootless)
export COMPOSE="uvx podman-compose"
make airflow-up                           # build + start postgres/scheduler/webserver
make dag-test DATE=2011-09-14             # backfill one day of ecallisto_ingest -> success
make ecallisto-day DATE=20110914          # the same day's scan WITHOUT Airflow (the shared worker)
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
