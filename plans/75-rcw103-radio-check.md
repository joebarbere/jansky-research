# 75 — RCW 103 (6.67-hr ultra-slow magnetar): a one-target radio counterpart check

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

1E 161348−5055 in RCW 103 is the 6.67-hr ultra-slow magnetar — the period sits squarely in
long-period-transient territory, yet no systematic radio counterpart check at its position has
been *stated* in the literature (fable-ideas F38: a never-stated gap). This is a one-target
appendix in the pattern of the merged `stokesv` and `lpt` slices: forced Stokes I+V photometry
at the known X-ray position across every available RACS epoch (CASDA, verified access) plus
VLASS quick-look cone checks, with leakage vetting. The expected outcome is a clean null — and
the expected null, with quantified per-epoch limits at a physically interesting period, is the
deliverable. A detection would be major; the plan does not assume one.

## Deliverables

- `src/jansky_research/rcw103.py`: `fetch_epochs` (CASDA cutouts + VLASS QL, network,
  `# pragma`), forced I/V photometry at the X-ray position (`measure_circular_pol` reuse from
  `stokesv`), leakage vetting (`classify_emitter` reuse), `limit_table` (per-epoch 3σ I and V
  limits + luminosity limits at the SNR distance), a 6.67-hr phase-coverage summary of the
  epochs, `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/rcw103/`; `survey/rcw103-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text pass on the RCW 103 literature and the LPT reviews to confirm no
   published targeted radio-continuum limit table exists; pin the best X-ray position and the
   adopted distance; verify epoch availability at Dec ≈ −51° (RACS yes; VLASS marginal — check).
1. Tooling + synthetic recover-a-known: injected point source at the target position in a real
   cutout recovered at the correct flux; the SRSC/GJ 65-style positive-control pattern from
   `stokesv_discovery` re-run before touching the target.
2. Real leg: forced photometry across all epochs; limit table; epoch times folded at 6.67 hr to
   state which rotation phases the limits actually cover.
3. GATE-2 (SNR-field confusion at the position; leakage floor in V; phase-coverage honesty) →
   paper (RNAAS/appendix-scale; the expected null framed as the result).

## Verification

Positive controls (known V emitter + injected source) recovered before the target is measured;
a clean expected-null with quantified limits is the point, not a failure mode; checks green;
GATE-2 sign-off.

## Risks & mitigations

- **SNR diffuse emission contaminates forced photometry** → local-background annulus sized
  against the SNR shell; contamination quoted in the limit budget.
- **Low duty cycle means epochs can all miss** → the 6.67-hr phase-coverage table makes the
  miss probability explicit rather than hidden.
- **Modest ceiling** → one-target appendix framing from the outset; near-total `stokesv`/`lpt`
  reuse keeps the cost a few days.
