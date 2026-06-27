# Convenience targets for jansky-research. Run `make help` for the list.
# Everything goes through uv (pinned Python 3.12). Mirrors the jansky course's
# conventions and supersets them with survey/airflow/paper targets.

.DEFAULT_GOAL := help
.PHONY: help setup test cov typecheck lint fmt fetch-data pipeline airflow-up airflow-down dag-test paper reproduce clean

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

pipeline: ## Run the analysis end-to-end WITHOUT Airflow (the shared code path)
	uv run python -m jansky_research.pipeline $(ARGS)

airflow-up: ## Stand up the local Airflow stack (podman compose)
	$(COMPOSE) -f airflow/compose.yaml up -d

airflow-down: ## Tear down the local Airflow stack
	$(COMPOSE) -f airflow/compose.yaml down

dag-test: ## Run the research DAG once without the scheduler loop
	$(COMPOSE) -f airflow/compose.yaml run --rm airflow-scheduler \
		airflow dags test research_pipeline 2026-01-01

paper: ## Build paper/main.pdf in the tectonic container
	$(COMPOSE) -f containers/compose.yaml run --rm paper

reproduce: ## Full reproduction: fetch data -> pipeline -> paper
	$(MAKE) fetch-data && $(MAKE) pipeline && $(MAKE) paper

clean: ## Remove caches and build artefacts
	rm -rf .pytest_cache/ .ruff_cache/ .mypy_cache/ site/
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
