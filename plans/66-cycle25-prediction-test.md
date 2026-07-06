# 66 — Cycle-25 prediction test: Compagnino & Zuccarello 2021 vs the real cycle past maximum

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

Compagnino & Zuccarello 2021 (arXiv:2103.13699) published a falsifiable cycle-25 prediction —
radio-loud and halo-CME rates extrapolated from cycle-24 correlations. Cycle 25 is now past
maximum and nobody has published a check against the realized cycle (fable-ideas F29: "fastest
honest note"). Data are three public lists: LASCO CME catalogue v2, SILSO sunspot numbers, and
RSTN radio-burst event lists — all no-auth, all already familiar from the merged solar slices.
Near-zero new code: the occurrence-rate/correlation statistics reuse `ecallisto_census` almost
wholesale. Whether the prediction verifies or fails, the comparison note is citable either way.

## Deliverables

- `src/jansky_research/cycle25.py`: `fetch_lasco_cme` / `fetch_silso` / `fetch_rstn` (network,
  `# pragma`), `refit_cycle24` (reproduce their published correlation coefficients),
  `cycle25_rates` (observed radio-loud/halo-CME rates per cycle phase), `compare_prediction`
  (predicted vs observed with honest few-cycle uncertainties), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/cycle25/`; `survey/cycle25-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text read of arXiv:2103.13699 to extract the exact predicted quantities and
   fit coefficients; confirm no published cycle-25 verification exists (ADS citation sweep);
   verify LASCO v2 / SILSO / RSTN list URLs are live.
1. Tooling + synthetic recover-a-known: a fixture cycle with known injected correlation
   coefficients round-trips through the refit machinery.
2. Cycle-24 leg: refit their coefficients from the same public lists — the load-bearing anchor.
3. Cycle-25 leg: observed rates through the current epoch vs their prediction, with a stated
   partial-cycle caveat (declining phase not yet complete).
4. GATE-2 (partial-cycle coverage; catalogue-definition drift between LASCO versions; RSTN
   station-availability gaps) → paper. This is an RNAAS/note-scale deliverable — say so.

## Verification

Refit reproduces the published cycle-24 coefficients within stated errors before any cycle-25
claim; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Cycle 25 not finished** → frame as a to-date test with explicit phase coverage, not a final
  verdict; the note can be updated at cycle end.
- **Catalogue definition drift** (halo-CME criteria, RSTN event classes) → apply their stated
  cuts verbatim; document any ambiguity in the findings file.
- **Someone publishes the check first** → GATE-0 citation sweep; the note is one-day-scale, so
  ship fast or drop.
