# 70 — LoTSS DR3 faint-count injection spot-check: an independent completeness harness

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

LoTSS DR3 ships its own faint source counts and completeness corrections; an *independent*
injection-recovery spot-check of those counts, run with this repo's own harness, is a
reproducibility note, not discovery — and is framed as exactly that (fable-ideas F33). The
injection-recovery completeness machinery derives directly from the merged `stacking` slice
(inject known sources into real mosaics, measure the recovered fraction vs flux). This slice is
also the natural migration of `plans/29` (the DR2-era count work) to DR3 — that lineage is
stated here so the older plan is superseded, not duplicated. Data: LoTSS DR3 mosaics/cutouts and
catalogue via the public lofar-surveys.org access points (no auth).

## Deliverables

- `src/jansky_research/dr3counts.py`: `fetch_dr3_field` / `fetch_dr3_catalog` (network,
  `# pragma`), `inject_sources` (harness adapted from `stacking`), `completeness_curve`
  (recovered fraction vs flux per field), `corrected_counts` (Euclidean-normalized counts with
  this harness's corrections vs the DR3 paper's own), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/dr3counts/`; `survey/dr3counts-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text read of the DR3 release paper's counts/completeness section — scope the
   spot-check against what they publish (no duplication of their headline analysis); verify
   mosaic/cutout download for 2–3 representative fields; formally mark plans/29 as migrated here.
1. Tooling + synthetic recover-a-known: injections into a noise-only fixture recover the input
   completeness curve exactly.
2. Real leg: injection sweeps in a small set of representative DR3 fields (rms-diverse), derive
   independent completeness curves, apply to the public catalogue, compare corrected faint
   counts against the DR3 paper's own published counts.
3. GATE-2 (field-to-field rms variation vs the small field sample; source-morphology assumptions
   in the injected population; resolution bias) → paper (RNAAS/note-scale reproducibility check).

## Verification

The recover-a-known anchor is the DR3 paper's own counts: agreement within stated errors is the
expected result and the deliverable; disagreement gets bounded, not spun; checks green; GATE-2
sign-off.

## Risks & mitigations

- **Reproducibility note, not discovery** → stated up front; the value is an independent harness
  and a citable consistency check for every DR3-counts user.
- **Few fields ≠ the full survey** → choose rms-stratified fields and report per-field spread as
  the systematic; no all-survey extrapolation beyond that.
- **Injected-morphology mismatch** → inject both point sources and modestly-resolved sources;
  report the difference as part of the error budget.
