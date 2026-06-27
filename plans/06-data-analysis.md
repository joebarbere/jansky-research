# Plan 06 — Run the analysis on real data (GATE 2) 📋

> Context: the actual research result. Depends on 03, 04, 05. Pairs with 07. Scope: medium.

## Context

Run the new tooling end-to-end on the chosen **real** public dataset (not synthetic), capture the
results, and decide honestly whether they are worth writing up. This is where the "interesting
findings" either materialise or don't — and the gate exists to keep the paper honest.

## Deliverables

- Generated `paper/figures/*.pdf` (vector) and `paper/generated/macros.tex` (`\newcommand` per
  headline number — nothing typed by hand).
- `results/metrics.json` — the machine-readable result manifest the paper reads.
- `notebooks/02_analysis_walkthrough.ipynb` — the analysis narrative on the real data.
- A short `survey/findings.md` — what was found, with uncertainties and caveats.

## Approach

Run via `make pipeline` (and confirm the Airflow DAG produces byte-identical artifacts). Drive the
real-data run with the **dataset-analyst** agent; **results-interpreter** drafts the honest
interpretation and cross-checks the literature.

## Verification

- **GATE 2 (hard stop):** **science-reviewer** + the human confirm the findings are real, honest,
  and non-trivial (no overclaiming; uncertainties quantified; nulls reported as nulls) before any
  paper drafting.
- `results/metrics.json` + figures regenerate deterministically from `make pipeline`.
