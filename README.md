# jansky-research

**Amateur radio-astronomy research, end to end.** A sibling of the
[`jansky`](https://github.com/joebarbere/jansky) teaching course: where jansky *teaches* radio
astronomy, this repo *does* a small slice of original amateur research and writes it up.

The arc is a single **vertical slice** — depth over breadth:

1. **Survey** the open-source and scientific landscape of radio astronomy (GitHub tooling +
   current literature) and find a real **gap**.
2. Build **new, tested Python tooling** that fills the gap, reusing jansky's helpers.
3. Use **Claude Code agents** to drive that tooling over an **existing public dataset**.
4. **Automate** the analysis with **Apache Airflow** on **Podman/Docker** local infrastructure.
5. Produce an **amateur research paper** (PDF, LaTeX / AASTeX) presenting the survey, the new
   code, the analysis, and any findings — authored by **Joe Barbere** and **Claude**.

> Status: **scaffolding + survey (P0–P1).** The science focus and the gap-filling tool are
> chosen by the survey and committed at **GATE 1** (see `plans/02-gap-analysis.md`).

## How it relates to `jansky`

This repo **depends on `jansky` as a library** and reuses its tested helpers
(`jansky.transients`, `jansky.rfi`, `jansky.timing`, `jansky.seti`, `jansky.sourcecounts`,
`jansky.formats`, `jansky.data`, …) rather than reimplementing them. It also mirrors jansky's
conventions: `uv`-managed, ruff + mypy + pytest with an 85% coverage floor, Podman containers,
and a `plans/NN-slug.md` workflow.

The `jansky` dependency is wired as a **local path source** (`../jansky`) for cross-repo
development; for a clean-checkout release it switches to a pinned git tag (`jansky@v0.1.0`). See
`pyproject.toml`.

## Quickstart

```bash
# Requires the jansky repo checked out next to this one (../jansky) for local dev.
uv sync          # creates the env; installs jansky from ../jansky
make test        # unit tests (offline, on synthetic fixtures)
make cov         # tests + 85% coverage floor
make fetch-data  # list the research datasets
```

## Layout

```
jansky-research/
  src/jansky_research/   # the tooling package (tested-helper pattern)
    data.py              # dataset registry + offline synthetic fallback
    pipeline.py          # fetch -> analyze -> metrics (shared by Make/notebook/Airflow) [P3]
    report.py            # figure/table emitters -> paper inputs [P3]
  survey/                # the deep-research survey outputs (markdown)
  airflow/               # Airflow-on-Podman stack + the research DAG [P5]
  paper/                 # LaTeX (AASTeX) sources; figures/ + generated/ are produced [P7]
  containers/            # paper (tectonic) build image
  notebooks/             # survey synthesis + analysis walkthrough
  plans/                 # the numbered project plans (deleted after merge)
  .claude/               # data-analysis agents [P4]
```

## The plan

The work is broken into numbered plans under `plans/` (`00-scaffold` … `08-reproducibility`),
each in the jansky house format (Context / Deliverables / Approach / Verification). Two human
**gates** guard quality: **GATE 1** commits the domain + dataset before any tooling is built, and
**GATE 2** confirms the results are real, honest, and non-trivial before the paper is drafted.

## License

MIT — see [LICENSE](LICENSE).
