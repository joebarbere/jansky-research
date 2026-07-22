# Contributing to jansky-research

Thanks for your interest in this project! `jansky-research` does original,
reproducible amateur radio-astronomy research as self-contained "slices" — one gap →
one tested tool → real public data → an honest write-up. This guide covers how to get
support, report a problem, and contribute a change, and the checks CI enforces.

## Getting support and reporting issues

- **Questions, support, or ideas** → open a [GitHub Discussion or Issue](https://github.com/joebarbere/jansky-research/issues)
  with the `question` label, or start with the [README](README.md), [`docs/usage.md`](docs/usage.md),
  and [`docs/faq.md`](docs/faq.md) (how others use the toolkit; how the papers in this repo work).
- **Bugs** → open an [issue](https://github.com/joebarbere/jansky-research/issues/new)
  with: what you ran (the exact `python -m jansky_research.<slice> …` command or `make`
  target), what you expected, what happened (full traceback), and your OS + Python
  version (`python --version`, `uv --version`). A failing case that reproduces with
  `--offline` (no network) is easiest to act on.
- **A result looks wrong** → that is the most valuable report of all. This project's
  first rule is honesty; if a number, figure, or claim looks off, open an issue and cite
  the slice's `survey/<slice>-findings.md` and `papers/<slice>/`.
- **Security / sensitive** → email the author (see `CITATION.cff`) rather than filing a
  public issue.

By participating you agree to the [Code of Conduct](CODE_OF_CONDUCT.md).

## Set up the environment

The project is managed with [uv](https://docs.astral.sh/uv/) and depends on the
[`jansky`](https://github.com/joebarbere/jansky) course library. For local development
`jansky` is resolved as a **sibling checkout** (`../jansky`), so clone both repos next
to each other:

```bash
git clone https://github.com/joebarbere/jansky.git
git clone https://github.com/joebarbere/jansky-research.git
cd jansky-research
uv sync                                   # core + dev toolchain (+ ../jansky, editable)
uv run python -c "import jansky_research"  # sanity check
```

Everything runs through `uv run`, so you never activate a venv by hand. `make help`
lists the shortcuts. See [`docs/usage.md`](docs/usage.md) for running a slice.

## The checks (run these before opening a PR)

CI runs all of these on every PR, across Python 3.10 and 3.12. Reproduce locally:

| Command | What it checks |
|---|---|
| `make lint`      | `ruff check src/ tests/` — lint |
| `make fmt`       | `ruff format` — auto-format |
| `make typecheck` | `mypy` over `src/jansky_research` |
| `make test`      | `pytest` unit tests (offline, on synthetic fixtures) |
| `make cov`       | unit tests with the **85% coverage floor** |

A one-liner before pushing: `make lint typecheck cov`. Papers additionally build with
`make paper` (containerised tectonic; needs `podman`).

## The slice pattern

Every result is built the same way, and a new analysis should follow it:

1. **A tested tool** in `src/jansky_research/<slice>.py` — pure NumPy/SciPy/Astropy (or
   pure PyTorch where device portability pays), with `__all__`, a documented `run()`,
   and a **mandatory synthetic offline fixture** so tests never touch the network. Keep
   it above the 85% coverage floor.
2. **A real-data run** on a public dataset (network code is marked `# pragma: no cover`).
3. **Honest framing** — validations, limits, and negatives reported plainly, never
   dressed-up discoveries. Every headline number must regenerate from the pipeline
   (figures as PDFs, numbers as `generated/macros.tex`), not be typed by hand.
4. **The write-up** — `survey/<slice>-findings.md` (the honest assessment) and an AASTeX
   paper in `papers/<slice>/` (`main.tex` + `refs.bib` tracked; figures/macros generated).

Plans live in `plans/NN-*.md`. `CLAUDE.md` documents the conventions in more depth.

## Pull requests

- **Branch before committing — never commit on `main`.** Open a PR; it is squash-merged
  and the branch deleted.
- Keep PRs focused (ideally one slice or one concern).
- **Versioning:** this repo follows SemVer per [`VERSIONING.md`](VERSIONING.md). Every PR
  adds an entry to the `## [Unreleased]` section of [`CHANGELOG.md`](CHANGELOG.md) under
  the right heading (Added / Changed / Fixed / …). `python scripts/next_version.py`
  recommends the next release number from those entries.
- New physics/constants belong in code with a unit test, not hard-coded in a paper.

## Use of AI assistance

This project is developed collaboratively with an AI coding assistant, disclosed in each
paper and in `joss/paper.md`. An AI/LLM is **not** an eligible author; contributors
remain responsible for reviewing, verifying, and taking accountability for all code,
results, and citations. If you use AI assistance in a contribution, disclose it in the
PR and confirm you have verified the output.
