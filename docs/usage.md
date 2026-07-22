# Usage

A short guide to installing `jansky-research` and running an analysis slice. For the
full reproduction (every paper, both orchestrators), see [`REPRODUCING.md`](../REPRODUCING.md);
for the per-slice inventory and outcomes, see the [README](../README.md).

## Install

The project is managed with [uv](https://docs.astral.sh/uv/) and depends on the
[`jansky`](https://github.com/joebarbere/jansky) course library, resolved as a **sibling
checkout**. Clone both repositories next to each other:

```bash
git clone https://github.com/joebarbere/jansky.git
git clone https://github.com/joebarbere/jansky-research.git
cd jansky-research
uv sync                                    # builds the env, incl. ../jansky (editable)
uv run python -c "import jansky_research"   # sanity check
```

Requires Python ≥3.10, <3.13. GPU pieces are opt-in: the core install and tests are
CPU-only.

## Run a slice offline (no network)

Every slice ships a synthetic fixture and runs fully offline with `--offline` — the
fastest way to see a tool work and to run its tests:

```bash
uv run python -m jansky_research.frbstats --offline     # FRB burst statistics
uv run python -m jansky_research.hi --offline           # Milky Way HI rotation curve
uv run python -m jansky_research.vgpra --offline        # Voyager 2 PRA reanalysis
make test                                               # the full offline test suite (85% floor)
```

Each slice is a Python module with a documented `run()`; `python -m jansky_research.<slice> --help`
lists its flags.

## Run on real public data

Drop `--offline` to fetch the real public dataset and run the analysis. Each slice
writes its results, figures (vector PDFs), and headline numbers (`generated/macros.tex`)
into `papers/<slice>/`:

```bash
uv run python -m jansky_research.frbperiod                              # FRB repeater periodicity (CHIME)
uv run python -m jansky_research.spectra --ra 180 --dec 30 --radius 3   # ultra-steep-spectrum hunt
uv run python -m jansky_research.vlass --ra 190 --dec 20 --radius 15    # VLASS variability (needs --extra vlass)
```

Some slices need an optional extra (e.g. `uv sync --extra vlass`, `--extra windwaves`);
the module prints a clear message if one is missing. GPU slices (`--extra fdmt`,
`--extra sbi`) run their accelerated leg from a ROCm/CUDA venv but fall back to CPU.

## Reproduce the papers

```bash
make figures     # build every static slice's inputs via the Snakemake DAG (offline)
make paper       # tectonic → papers/<slice>/main.pdf (in a container; needs podman)
make reproduce   # fetch → figures → papers → arXiv packages, end to end
```

The streaming e-Callisto ingest runs separately on Airflow + Podman — see
[`REPRODUCING.md`](../REPRODUCING.md).

## Cite

If you use this software or its analyses, cite the archived release — see
[`CITATION.cff`](../CITATION.cff) and the DOI badge in the [README](../README.md).
