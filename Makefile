# Convenience targets for jansky-research. Run `make help` for the list.
# Everything goes through uv (pinned Python 3.12). Mirrors the jansky course's
# conventions and supersets them with survey/airflow/paper targets.

.DEFAULT_GOAL := help
.PHONY: help setup test cov typecheck lint fmt fetch-data pipeline figures airflow-up airflow-down dag-test paper-image paper arxiv reproduce clean

# The research slices, each with a paper under papers/<slice>/.
SLICES ?= frbstats frbperiod driftsearch spectra hi vlass peaked southern offsets pulsarspec stacking vlbi solarbursts rmsky

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

figures: ## Regenerate every slice's figures + macros into papers/<slice>/ (offline synthetic; add ARGS= for real data per slice)
	uv run python -m jansky_research.pipeline --out . --offline
	uv run python -m jansky_research.frbperiod --out . --offline
	uv run python -m jansky_research.driftsearch --out .
	uv run python -m jansky_research.spectra --out . --offline
	uv run python -m jansky_research.hi --out . --offline
	uv run python -m jansky_research.vlass --out . --offline
	uv run python -m jansky_research.peaked --out . --offline
	uv run python -m jansky_research.stokesv --out . --offline
	uv run python -m jansky_research.southern --out . --offline
	uv run python -m jansky_research.offsets --out . --offline
	uv run python -m jansky_research.pulsarspec --out . --offline
	uv run python -m jansky_research.stacking --out . --offline
	uv run python -m jansky_research.vlbi --out . --offline
	uv run python -m jansky_research.solarbursts --out . --offline
	uv run python -m jansky_research.rmsky --out . --offline

airflow-up: ## Stand up the local Airflow stack (podman compose)
	$(COMPOSE) -f airflow/compose.yaml up -d

airflow-down: ## Tear down the local Airflow stack
	$(COMPOSE) -f airflow/compose.yaml down

dag-test: ## Run the research DAG once without the scheduler loop
	$(COMPOSE) -f airflow/compose.yaml run --rm airflow-scheduler \
		airflow dags test research_pipeline 2026-01-01

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
	$(MAKE) paper && $(MAKE) arxiv

clean: ## Remove caches and build artefacts
	rm -rf .pytest_cache/ .ruff_cache/ .mypy_cache/ site/
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
