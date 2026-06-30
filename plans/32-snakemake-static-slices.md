# 32 — Right-size the orchestration: Snakemake for the static slices, Airflow only for streaming

Status: ✅ done (workflow/Snakefile drives all static slices; make figures + paper.yml rewired; byte-identical verified)

## Context

Airflow is a *scheduler* for large, frequently-updated, dependency-heavy pipelines (its strengths:
cron schedules, backfill/catchup, sensors, parallel fan-out, retries, a running webserver+metadata DB).
Seventeen of the nineteen slices are the opposite: a single small **static** input (a VizieR catalogue,
one SPDF CDF, a 216 KB CSV) run once through `fetch → analyze → figure → macros → paper`. For those,
Airflow is the wrong tool — a daemon and a Postgres DB to run a four-step DAG over a CSV. The right tool
is a **file-target workflow runner**: declarative input→output rules, a dependency DAG built from the
files, parallel execution, dry-run/provenance, container/conda integration, and **no server**.

**Snakemake** (Mölder et al. 2021) is the standard such tool in reproducible computational science and
is the appropriate orchestrator for the static slices: each slice is a rule
(`results/<slice>_metrics.json`, `papers/<slice>/figures/*.pdf`, `generated/macros.tex`,
`papers/<slice>/main.pdf`) with explicit inputs/outputs, so `snakemake` rebuilds exactly what changed,
in parallel, reproducibly, and is itself citable. Make stays as the thin human entry point
(`make pipeline`/`make reproduce` shell out to Snakemake); Airflow is reserved for the genuinely
streaming/large datasets (plan #31, e-Callisto; and a future VLASS-epoch pipeline).

Net result: an honest **right-sized orchestration** story — Snakemake for the 17 static reproducible
slices, Airflow for streaming ingest — instead of one over-powered scheduler stretched over a CSV.

## Deliverables

- `workflow/Snakefile` (+ `workflow/rules/*.smk`) — a rule per slice expressing the real DAG:
  `fetch/synthetic → run → figures+macros → paper`, with `results/…`, `papers/<slice>/figures/…`,
  `papers/<slice>/generated/macros.tex`, `papers/<slice>/main.pdf` as targets and the slice module as
  the command. An `offline` config flag selects synthetic vs real inputs (mirrors `run(offline=...)`).
  Containerised paper builds via the existing tectonic image. `all`, `papers`, and per-slice targets.
- `snakemake` added as a dev/optional dependency (not a core runtime dep); `[project.optional-dependencies]`.
- `Makefile` rewired: `make pipeline` / `make figures` / `make paper` / `make reproduce` shell out to
  `snakemake -j` (so one DAG, one provenance graph), keeping the same human-facing target names.
- `.github/workflows/` — the offline-synthetic build runs `snakemake --use-... -j` (replacing the long
  hand-listed slice loop in `paper.yml` with the DAG; CI still never runs Airflow).
- Docs: `REPRODUCING.md` + `README` updated to the two-tier story (Snakemake static / Airflow streaming);
  the frbstats paper's automation claim retired here and re-homed in #31.

## Approach

1. **One slice as the pattern.** Write the Snakefile rules for `sourcecounts` (fetch→run→figure+macros
   →paper), confirm `snakemake -j` reproduces byte-identical artifacts to `python -m
   jansky_research.sourcecounts`. Establish the offline/real config switch.
2. **Fan out to all static slices.** Parameterise one rule template over the `SLICES` list (Snakemake
   wildcards), so adding a slice is a one-line registry edit, not N Makefile lines. Keep `triangulate`,
   `swaves`, etc. real-data commands as the non-offline branch.
3. **Rewire Make + CI** to drive Snakemake; verify `make reproduce` still runs fetch→pipeline→papers
   →arXiv end to end, now as a single dependency DAG with parallelism and correct incremental rebuilds.
4. **Retire the Airflow-on-static claim.** Coordinated with #31: Airflow → streaming only; the frbstats
   paper drops Airflow; the automation narrative becomes "right-sized orchestration."

## Verification

- `snakemake -n` (dry run) shows the correct per-slice DAG; `snakemake -j` reproduces **byte-identical**
  `results/*.json`, figures, and macros vs the direct `python -m jansky_research.<slice>` path.
- Incremental correctness: touching one slice's module rebuilds only its targets (and its paper), not
  the whole repo.
- `make test`/`cov`/`ruff`/`mypy` unaffected (the orchestrator wraps the same tested code, 85% floor).
- `make reproduce` runs end to end via Snakemake; CI's offline-synthetic build is green driven by the DAG.
- Airflow is no longer invoked for any static slice; `make airflow-*` / `dag-test` target only the
  streaming pipeline (#31).

## Risks & mitigations

- **Two orchestrators looks like more, not less →** frame and document it as *right-sizing* (server-less
  file-DAG for static, scheduler for streaming); Make remains the single human entry point so day-to-day
  usage is unchanged.
- **Snakemake/tectonic-container interaction →** reuse the existing paper image; run the LaTeX build as a
  rule action (same `podman run` the Makefile uses), not Snakemake's conda integration.
- **CI churn →** swap `paper.yml`'s hand-listed loop for `snakemake -j` only after the dry-run DAG is
  verified to cover every slice; keep the offline-synthetic guarantee.
- **Scope creep into a rewrite →** the slice *modules* are untouched; this is an orchestration layer
  over them, added one slice at a time behind the unchanged `make` targets.
