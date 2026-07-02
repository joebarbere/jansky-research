# 35 — The first long-period-transient (LPT) population catalogue and P–Ṗ diagram

Status: 📋 planned — the top new-findings pick from `survey/opportunity-scan-2026-07.md` Tier 1 #3; data compilation running

## Context

Long-period radio transients (periods minutes–hours, far beyond the pulsar death line) grew from 2
objects (2022) to ~15 confirmed by mid-2026, spanning magnetar and white-dwarf-binary
interpretations. The authoritative review (Rea, Hurley-Walker & Caleb 2026, arXiv:2601.10393)
**explicitly notes no population synthesis exists** for the confirmed sample. This slice builds the
first one: a verified machine-readable table of every confirmed LPT (period, Ṗ or limit, activity
fraction, fluxes, distance, binary/X-ray status, provenance per value) and the class's first
**P–Ṗ diagram** — placing LPTs against the pulsar population, death line, and magnetar/WD tracks
(direct `ppdot.py` reuse) — plus a VLASS-QL2/LoTSS-DR3-cutout counterpart cross-check per source.
A compilation is honest by construction: every number carries its arXiv provenance; the
contribution is the systematic table + diagram, not a discovery.

## Deliverables

- `data/lpt_sample.csv` (vendored, committed — the verified table with provenance columns) +
  `src/jansky_research/lpt.py`: loaders, unit normalisation, P–Ṗ placement (reuse
  `ppdot`/`jansky.transients.surface_bfield/characteristic_age/death_line_pdot` with the honest
  caveat that dipole formulae assume NS values), period-distribution + activity stats,
  `crossmatch_counterparts` (VLASS QL2 cone + LoTSS DR3 cutout API, `# pragma`), synthetic
  round-trip (inject a fake population → recover distribution), `run/_figure/_write_macros/_main`.
- `tests/test_lpt.py` (85% floor, offline); `papers/lpt/` (AASTeX); `survey/lpt-findings.md`;
  Makefile/Snakefile wiring.

## Approach

0. **GATE 0 (running):** the compilation agent returns the verified per-object table; spot-check
   3 objects against their discovery papers before vendoring. STOP if values can't be verified.
1. Tooling + synthetic round-trip; diagram machinery validated on pulsar anchors (Crab).
2. Real leg: the vendored table → P–Ṗ diagram + stats; VLASS/LoTSS counterpart fluxes per source
   (bounded: ~15 cone queries + cutouts).
3. GATE-2: every table value re-checked against provenance; NS-formula caveat for WD systems;
   "first population catalogue" claim hedged against the review's own tables.
4. Paper: value-added catalogue + diagram, "not a discovery".

## Verification

Checks green (85% floor, ruff+mypy); offline round-trip recovers the injected population; every
CSV row carries an arXiv provenance; GATE-2 sign-off; `make reproduce` regenerates figure+macros.

## Risks & mitigations

- **The review team publishes the same diagram** → ship the bounded version fast; an independent
  reproducible compilation retains cross-check value.
- **Heterogeneous published values** (bands, definitions) → provenance + flag columns, never
  averaging across definitions; disputed values flagged, not chosen.
- **NS-formula misuse on WD systems** → B/τ columns computed only where the NS interpretation is
  viable; WD binaries plotted but not assigned dipole values.
