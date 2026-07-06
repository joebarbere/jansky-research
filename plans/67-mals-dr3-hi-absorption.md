# 67 — MALS DR3 Galactic HI-absorption demographics: covering fraction and optical-depth stats

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

MALS DR3 (arXiv:2504.00097; own portal `mals.iucaa.in`) released 3,640 Galactic HI-absorption
features across 19,130 sightlines — the largest homogeneous 21-cm absorption set yet. The release
paper is catalogue-first; the plain demographic questions are open (fable-ideas F30): covering
fraction of absorption vs Galactic |b| and ℓ, and the optical-depth distribution compared against
21-SPONGE. This is deliberately a low-ceiling, zero-drama slice — a clean statistical companion
to the release, not a discovery hunt — and the plan says so plainly. Spectral-line muscle memory
comes from the merged `hi` slice (`hi.py` readers and profile statistics extend naturally).

## Deliverables

- `src/jansky_research/malshi.py`: `fetch_mals_dr3` (portal/VizieR table, `# pragma`),
  `covering_fraction` (per-|b|/ℓ bins with binomial errors), `tau_distribution` (optical-depth
  histogram + survival statistics for non-detections), `compare_21sponge` (published
  distribution as the external anchor), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/malshi/`; `survey/malshi-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text read of arXiv:2504.00097 — confirm the demographic analysis is not in
   the release paper or an announced companion; verify the machine-readable catalogue is
   fetchable from `mals.iucaa.in` (or VizieR) and record columns/row count.
1. Tooling + synthetic recover-a-known: a mock sightline population with an injected |b|-dependent
   covering fraction and a known τ distribution round-trips through the binning/survival code.
2. Real leg: covering fraction vs |b| and ℓ, τ distribution with sensitivity-aware upper limits,
   21-SPONGE comparison on the overlapping τ range.
3. GATE-2 (heterogeneous per-sightline optical-depth sensitivity; feature-blending in the
   catalogue; Galactic vs intervening classification cuts) → paper (note-scale; RNAAS-adjacent).

## Verification

Synthetic round-trip recovers the injected covering-fraction profile; the τ distribution is
benchmarked against the published 21-SPONGE distribution before any new claim; checks green;
GATE-2 sign-off.

## Risks & mitigations

- **MALS team companion paper in flight** → GATE-0 checks their publication list; scope stays
  strictly demographic so overlap with physics papers is minimal.
- **Sensitivity heterogeneity across 19,130 sightlines** → censored (survival) statistics, not
  detection-only histograms; report the completeness cut explicitly.
- **Low ceiling by design** → framed as a statistics companion from the outset; no overclaiming.
