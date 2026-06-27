# Plan 00 — Scaffold the jansky-research repo ✅

> Context: foundation for the whole vertical slice. Pairs with every later plan. Scope: small.

## Context

`jansky-research` is a new sibling of the `jansky` course that turns *learning* radio astronomy
into a small slice of *original amateur research*. Everything downstream (survey → tooling →
automation → paper) needs a repo skeleton that depends on `jansky` as a library and mirrors its
conventions so the quality bar (ruff/mypy/pytest, 85% coverage, Podman, plan-driven work) carries
over.

## Deliverables

- `pyproject.toml` — depends on `jansky` (local path source `../jansky` for dev; pinned git tag
  `v0.1.0` for release/CI). Mirrors jansky's ruff (line 100) / mypy / pytest / coverage(85) config.
  Airflow kept out of core install as an extra.
- `src/jansky_research/` — `__init__.py` + `data.py` (dataset registry + offline synthetic
  fallback, `JANSKY_RESEARCH_DATA_DIR`), modelled on `jansky.data`.
- `tests/` — offline `test_data.py` + `test_smoke.py` (imports `jansky` to prove the dependency
  resolves).
- `Makefile` (superset of jansky's), `.github/workflows/ci.yml` (checks out jansky as a sibling),
  `README.md`, `LICENSE` (MIT — Joe Barbere & Claude), `.gitignore`.
- Skeleton dirs: `survey/ airflow/ paper/ notebooks/ plans/ .claude/ containers/ results/`.
- Tag `jansky` `v0.1.0` (the library dependency).

## Approach

Mechanical. Copy jansky's `data.py` registry + offline-fallback pattern. No agents.

## Verification

- `uv sync` resolves `jansky`; `make test` / `cov` pass on synthetic fixtures (≥85%); `make lint`
  + `mypy` clean. CI green on an empty-ish package.

**Status: done** — commit `c3682ec` (uv sync ok, ruff+mypy clean, 10 tests, 93% coverage).
