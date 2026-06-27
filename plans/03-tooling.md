# Plan 03 — New tooling for the chosen gap 📋

> Context: the novelty of the project. Depends on 02 (GATE 1). Pairs with 04, 05. Scope: medium.

## Context

With the domain + dataset committed at GATE 1, build the new, tested Python tool that fills the
gap — reusing `jansky`'s helpers rather than reimplementing them. This is the original
contribution the paper reports.

## Deliverables

- `src/jansky_research/<domain>.py` — the gap-filling algorithm (pure NumPy/astropy, `__all__`,
  docstring maths), **composing** `jansky` helpers (e.g. `from jansky import transients, rfi,
  timing, sourcecounts`).
- `src/jansky_research/pipeline.py` — `fetch → analyze → metrics` returning plain
  dicts/dataclasses; the **single entry point** shared by the Makefile, the notebooks, and the
  Airflow DAG; thin `python -m jansky_research.pipeline` CLI.
- `src/jansky_research/report.py` — matplotlib-`Agg` figure/table emitters writing
  `paper/figures/*.pdf` and `paper/generated/macros.tex`.
- `data.py` registry entry for the chosen dataset + a synthetic fallback generator.
- `tests/` — every `pipeline`/`report`/`<domain>` function tested against the synthetic fixture
  (deterministic seeds); ≥85% coverage; network `fetch` tested with a monkeypatched downloader.

## Approach

Mirror `jansky`'s tested-helper style exactly. Keep all I/O behind explicit out-paths so the same
functions run identically from Make, a notebook, and an Airflow task. Add the chosen `jansky`
extra to `pyproject.toml`. Optional **science-reviewer** pass on the algorithm.

## Verification

- `make cov` ≥85% on synthetic fixtures (no network); `mypy` + `ruff` clean.
- `python -m jansky_research.pipeline` runs end-to-end on the synthetic fixture and writes a
  `results/metrics.json` + at least one figure.
