# Plan 04 — Data-analysis agents 📋

> Context: Claude agents that drive the tooling over real data. Depends on 03. Pairs with 05, 06.
> Scope: small.

## Context

The project calls for Claude Code agents that *use* the new tooling to analyze existing public
datasets. Define a small set of focused agents (tight tool allowlists, like jansky's agents) and
wire in the reused jansky agents/skills.

## Deliverables

- `.claude/agents/dataset-analyst.md` (`Bash, Read, Glob, Grep`) — drives
  `jansky_research.pipeline`/CLI over a fetched public dataset, sanity-checks outputs (shapes,
  units, NaNs), reports metrics + which figures/tables were produced.
- `.claude/agents/pipeline-runner.md` (`Bash, Read`) — operates the local Airflow-on-Podman stack
  (`make airflow-up`, `airflow dags test`), diagnoses podman/SELinux mount failures, confirms
  artifacts landed.
- `.claude/agents/results-interpreter.md` (`Read, WebFetch, WebSearch, Glob, Grep`) — writes the
  honest results interpretation for the paper, cross-checks the literature, **flags overclaiming**.
- `notebooks/02_analysis_walkthrough.ipynb` — a human-readable walkthrough of the analysis.
- Reuse wiring: copy (or path-include) jansky's `archive-scout`, `radio-research-assistant`,
  `science-reviewer` into this repo's `.claude/` if the cross-repo indirection proves fragile.

## Approach

Model the agent files on `jansky/.claude/agents/archive-scout.md` (frontmatter: name,
description, tight `tools`). Keep the science agents clean and the ops agent separate.

## Verification

- Each agent file is valid and invokable; a dry run of `dataset-analyst` over the synthetic
  fixture produces a correct metrics summary.
