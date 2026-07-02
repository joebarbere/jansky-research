# Convenience targets for jansky-research. Run `make help` for the list.
# Everything goes through uv (pinned Python 3.12). Mirrors the jansky course's
# conventions and supersets them with survey/airflow/paper targets.

.DEFAULT_GOAL := help
.PHONY: help setup test cov typecheck lint fmt fetch-data pipeline figures figures-dry airflow-up airflow-down dag-test ecallisto-day paper-image paper arxiv reproduce clean

# The research slices, each with a paper under papers/<slice>/.
SLICES ?= frbstats frbperiod driftsearch spectra hi vlass peaked southern offsets pulsarspec stacking vlbi solarbursts rmsky ppdot windwaves swaves triangulate sourcecounts type3synthesis ecallisto_pipeline ecallisto_census torchfdmt stokesv stokesv_discovery

# Compose command. Fedora/podman often has no `podman compose` provider; `podman-compose`
# is the reliable driver. No install needed if you have uv:  COMPOSE="uvx podman-compose"
COMPOSE ?= podman-compose

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-16s\033[0m %s\n", $$1, $$2}'

setup: ## Create the environment (Python 3.12 + deps, incl. jansky)
	uv sync

test: ## Run the unit tests
	uv run pytest

cov: ## Run the unit tests with a coverage report (enforces the 85% floor)
	uv run pytest --cov=jansky_research --cov-report=term-missing

typecheck: ## Type-check the package with mypy
	uv run mypy

lint: ## Lint with ruff
	uv run ruff check src/ tests/

fmt: ## Auto-format with ruff
	uv run ruff format src/ tests/

fetch-data: ## List research datasets (use ARGS="--fetch NAME" to download)
	uv run python -m jansky_research.data $(ARGS) --list

pipeline: ## Run the FRB burst-statistics analysis WITHOUT Airflow (the shared code path)
	uv run python -m jansky_research.pipeline $(ARGS)

figures: ## Regenerate every static slice's figures + macros via the Snakemake DAG (offline synthetic)
	uv run --extra workflow snakemake -s workflow/Snakefile -j$(or $(J),4)

figures-dry: ## Show the static-slice DAG without running it (snakemake dry-run)
	uv run --extra workflow snakemake -s workflow/Snakefile -n

airflow-up: ## Stand up the local Airflow stack (podman compose)
	$(COMPOSE) -f airflow/compose.yaml up -d

airflow-down: ## Tear down the local Airflow stack
	$(COMPOSE) -f airflow/compose.yaml down

dag-test: ## Backfill one day of the e-Callisto ingest DAG under Podman (DATE=YYYY-MM-DD)
	$(COMPOSE) -f airflow/compose.yaml run --rm airflow-scheduler \
		airflow dags test ecallisto_ingest $(or $(DATE),2011-09-14)

ecallisto-day: ## Scan one day of e-Callisto for type III candidates WITHOUT Airflow (DATE=YYYYMMDD)
	uv run python -m jansky_research.ecallisto_catalog --date $(or $(DATE),20110914) --out .

paper-image: ## Build the tectonic image used to compile the papers
	podman build -t jansky-research-paper:latest -f containers/paper.Dockerfile .

paper: paper-image ## Build every papers/<slice>/*.tex (main + e.g. rnaas) in the tectonic container
	@for s in $(SLICES); do \
		for tex in papers/$$s/*.tex; do \
			grep -ql '\\documentclass' "$$tex" || continue; \
			echo "==> building $$tex"; \
			podman run --rm -v "$(CURDIR)/papers/$$s":/paper:z -w /paper jansky-research-paper:latest \
				tectonic --keep-intermediates --keep-logs "$$(basename $$tex)" || exit 1; \
		done; \
	done

arxiv: ## Assemble + validate an arXiv package for every paper (papers/<slice>/arxiv-submission/)
	@for s in $(SLICES); do \
		echo "==> packaging papers/$$s"; \
		uv run python .claude/skills/arxiv-submit/assemble_arxiv.py \
			--paper papers/$$s --out papers/$$s/arxiv-submission || exit 1; \
	done

reproduce: ## Full reproduction on REAL public data -> figures+macros -> papers -> arXiv packages
	uv run python -m jansky_research.pipeline --out .
	uv run python -m jansky_research.frbperiod --out .
	uv run python -m jansky_research.driftsearch --out .
	uv run python -m jansky_research.spectra --ra 180 --dec 30 --radius 3 --out .
	uv run python -m jansky_research.hi --out .
	uv run python -m jansky_research.peaked --ra 180 --dec 30 --radius 2 --validate --out .
	uv run python -m jansky_research.southern --ra 30 --dec -30 --radius 3 --validate --out .
	uv run python -m jansky_research.offsets --out .
	uv run python -m jansky_research.pulsarspec --out .
	uv run python -m jansky_research.stacking --ra 180 --dec 25 --radius 2.5 --max-sources 250 --out .
	uv run python -m jansky_research.vlbi --online --out .
	uv run python -m jansky_research.solarbursts --recover --out .
	uv run python -m jansky_research.rmsky --out .
	uv run python -m jansky_research.ppdot --out .
	uv run python -m jansky_research.windwaves --date 20031028 --receiver rad2 --out .
	uv run python -m jansky_research.swaves --date 20130515 --spacecraft a --out .
	uv run python -m jansky_research.triangulate --date 20130515 --out .
	uv run python -m jansky_research.sourcecounts --ra 180 --dec 30 --radius 8 --out .
	uv run python -m jansky_research.type3synthesis --out .
	uv run python -m jansky_research.ecallisto_catalog --date 20110914 --out .
	uv run python -m jansky_research.ecallisto_census --offline --out .  # synthetic census; real multi-cycle e-Callisto ingest is future work (coverage-limited, see survey/)
	uv run python -m jansky_research.singlepulse --benchmark --out .  # real Crab recover-a-known + CPU benchmark (GPU: rerun with --device cuda in a ROCm venv)
	CASDA_USERNAME=$${CASDA_USERNAME:?set CASDA_USERNAME for the stokesv Stokes-V leg} uv run python -m jansky_research.stokesv --out .
	uv run python -m jansky_research.stokesv_discovery --out .  # summarises results/stokesv_discovery_realtargets.csv (regenerate: uv run python scripts/stokesv_discovery_real.py)
	$(MAKE) paper && $(MAKE) arxiv

clean: ## Remove caches and build artefacts
	rm -rf .pytest_cache/ .ruff_cache/ .mypy_cache/ .snakemake/ site/
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
