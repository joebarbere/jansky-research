# jansky-research

**Amateur radio-astronomy research, end to end.** A sibling of the
[`jansky`](https://github.com/joebarbere/jansky) teaching course: where jansky *teaches* radio
astronomy, this repo *does* original amateur research — surveys the landscape, builds tested
tooling, analyzes public data, automates it, and writes it up — with honesty as the first rule.

It began as a single vertical slice (survey → gap → tool → Airflow pipeline → reproducible paper)
and has grown into a **deep-research survey plus a set of self-contained research slices**, each:
a tested CPU-only tool reusing jansky's helpers, run on real public data, put through an
adversarial science-review gate, and written up — wins *and* negatives reported plainly.

## Status — five slices, honestly tallied

| Slice | Tool | Outcome |
|-------|------|---------|
| FRB burst-statistics | `jansky_research.frbstats` | ✅ reproduced the CHIME repeater **width** result |
| Ultra-steep-spectrum hunt | `jansky_research.spectra` | ➖ USS candidates **failed** the de Gasperin cross-check (negative) |
| FRB repeater periodicity | `jansky_research.frbperiod` | ✅ recovered FRB 20180916B's **16.35-day** period |
| SETI drift-search benchmark | `jansky_research.driftsearch` | ➖ benchmark built; the "Voyager detection" was a **DC-spike artifact** (retracted) |
| HI rotation curve | `jansky_research.hi` | ✅ recovered the **flat** (non-Keplerian) inner Milky Way curve |

Three clean validations, two honest negatives, **zero overclaims that survived review** — the
science-reviewer caught the USS candidates evaporating against the authoritative catalog and the
SETI "detection" being an instrument artifact. The negatives are arguably the most instructive part.
Each slice's honest assessment is in `survey/*-findings.md`.

### Papers

Every slice is written up as its own **AASTeX paper** under `papers/<slice>/` (authored by Joseph
Barbere, with Claude credited via an AI-use disclosure + a `\software{}` citation — an AI/LLM is not
an eligible author). Each is built reproducibly with containerized tectonic, takes every headline
number from a pipeline-generated `generated/macros.tex` (no figure typed by hand), and is honest
about what it is — three validations and two negatives:

| Paper | `papers/…` | Framing |
|-------|-----------|---------|
| FRB burst statistics, validated on CHIME/FRB Cat 1 | `frbstats/` | validation + Airflow-on-Podman automation |
| Recovering FRB 20180916B's 16.35-day period | `frbperiod/` | validation |
| The flat inner Milky Way rotation curve from LAB HI | `hi/` | validation |
| A CPU-only SETI drift-search benchmark + Voyager-1 null | `driftsearch/` | benchmark + honest negative |
| TGSS×NVSS USS selection is dominated by the flux scale | `spectra/` | cautionary negative |

`make paper` builds all five PDFs; `make arxiv` runs the bundled **`arxiv-submit` skill**
(`.claude/skills/arxiv-submit/`) to assemble and validate an upload package per paper
(`papers/<slice>/arxiv-submission/`: the LaTeX-source tarball with its `.bbl`, plus a `metadata.yaml`
capturing every arXiv submission property and a `CHECKLIST.md`). The first slice is also automated by
an **Apache Airflow pipeline on rootless Podman** (`airflow/`).

## How it relates to `jansky`

This repo **depends on `jansky` as a library** and reuses its tested helpers (`jansky.transients`,
`jansky.rfi`, `jansky.timing`, `jansky.seti`, `jansky.sourcecounts`, `jansky.formats`,
`jansky.data`, …) rather than reimplementing them. It mirrors jansky's conventions: `uv`-managed,
ruff + mypy + pytest with an 85% coverage floor, Podman containers, and a `plans/NN-slug.md`
workflow. The `jansky` dependency is a local path source (`../jansky`) for cross-repo dev, switching
to the pinned git tag `jansky@v0.1.0` for clean-checkout CI. See `pyproject.toml`.

## Quickstart

```bash
# Requires the jansky repo checked out next to this one (../jansky) for local dev.
uv sync                                   # env + jansky (from ../jansky)
make test                                 # unit tests (offline, on synthetic fixtures)
make cov                                  # tests + 85% coverage floor

# Run a slice on real public data (each writes results/ + figures + macros into papers/<slice>/):
uv run python -m jansky_research.pipeline     # FRB burst statistics (CHIME catalog)
uv run python -m jansky_research.frbperiod    # FRB repeater periodicity
uv run python -m jansky_research.spectra --ra 180 --dec 30 --radius 3   # USS hunt
uv run python -m jansky_research.driftsearch  # SETI injection-recovery benchmark
uv run python -m jansky_research.hi           # Milky Way HI rotation curve
# (append --offline to run any slice on its synthetic fixture, no network)

# The papers + automation:
make figures                              # regenerate every slice's figures + macros (offline)
make paper                                # tectonic -> all papers/<slice>/main.pdf (in a container)
make arxiv                                # assemble + validate an arXiv package per paper
make airflow-up COMPOSE="uvx podman-compose" && make dag-test   # run the DAG under Podman
make reproduce                            # fetch -> figures -> papers -> arXiv packages, end to end
```

See `REPRODUCING.md` for the full reproduction, the Airflow-on-Podman notes, and offline mode.

## Layout

```
jansky-research/
  src/jansky_research/   # the tooling package (tested-helper pattern, 85% floor)
    data.py              # dataset registry + offline synthetic fallback
    frbstats.py          # FRB burst statistics (Weibull / power-law / KS)
    spectra.py           # radio spectral index + ultra-steep-spectrum hunt
    frbperiod.py         # FRB repeater activity-periodicity (Rayleigh periodogram)
    driftsearch.py       # SETI Doppler-drift injection-recovery benchmark
    hi.py                # Milky Way HI tangent-point rotation curve
    pipeline.py          # the FRB pipeline (shared by Make / notebook / Airflow)
    report.py            # figure/macro emitters -> paper inputs
  survey/                # PERMANENT: literature.md, github-landscape.md, gap-analysis.md,
                         #   candidate-gaps.md (backlog), and each slice's *-findings.md
  airflow/               # Airflow-on-Podman stack + the research DAG
  papers/<slice>/        # one AASTeX paper per slice (main.tex + refs.bib tracked;
                         #   figures/, generated/, arxiv-submission/ are produced by make)
  containers/            # tectonic paper-build image
  .claude/skills/        # arxiv-submit (assemble + validate an arXiv upload package)
  .claude/agents/        # dataset-analyst, pipeline-runner, results-interpreter (+ reused jansky)
  plans/                 # numbered project plans (deleted after merge; survey/ is the keep-file)
```

## Method & gates

Each slice follows the same discipline, recorded as a `plans/NN-slug.md` spec: build a tested tool →
run it on real public data → **science-reviewer gate** (an adversarial pass that has caught real
blockers every time) → honest write-up. Findings are framed as *candidates / validations / limits*,
never dressed-up discoveries; every reported number is reproducible. The non-chosen survey gaps are
preserved as a backlog in `survey/candidate-gaps.md`.

## License

MIT — see [LICENSE](LICENSE).
