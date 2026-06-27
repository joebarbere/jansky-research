# Plan 08 — Reproducibility & release (capstone) 📋

> Context: make the whole slice reproducible in one command. Depends on 07. Scope: small.

## Context

A research artifact is only as good as its reproducibility. Close the loop so anyone can
regenerate the analysis and the paper from a clean checkout, and document how.

## Deliverables

- `make reproduce` — `fetch-data → make pipeline → make paper` end to end.
- `REPRODUCING.md` — exact steps (local path-source vs git-tag jansky; the Airflow-on-Podman
  alternative; offline/synthetic mode).
- README polish; final pass over CI (`ci`, `notebooks`, `links`, `paper`).
- Cleanup: delete merged **plan** files per the jansky convention. **Never delete `survey/`** —
  it is a permanent research artifact, and `survey/candidate-gaps.md` is the future-work backlog
  of the gaps GATE 1 did *not* pick.

## Approach

Confirm the DAG and `make pipeline` produce byte-comparable paper inputs. Pre-warm the tectonic
package cache in the paper image for offline builds.

## Verification

- `make reproduce` runs clean from a fresh checkout and yields `paper/main.pdf`.
- The DAG and `make pipeline` produce identical `paper/figures/` + `paper/generated/macros.tex`.
- All CI workflows green.
