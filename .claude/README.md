# Claude Code automation for jansky-research

This directory holds the project's data-analysis **agents** (added in P4, see
`plans/04-agents.md`):

- `dataset-analyst` — drives `jansky_research.pipeline` over a fetched public dataset.
- `pipeline-runner` — operates the local Airflow-on-Podman stack.
- `results-interpreter` — writes the honest results interpretation, flags overclaiming.

It **reuses** the `jansky` course's agents and skills (resolved when you open this repo
alongside `../jansky`, or copied here in P4 if the cross-repo indirection proves fragile):

- Agents: `archive-scout` (data discovery), `radio-research-assistant` (literature),
  `science-reviewer` (the GATE 1 / GATE 2 / paper correctness gates).
- Skills: `find-radio-papers`, `dataset-watch`, `radio-source-lookup`.

Until P4, only this note lives here.
