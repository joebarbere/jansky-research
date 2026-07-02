# 36 — Galactic RM structure vs latitude from SPICE-RACS (rmsky at 10× the data)

Status: 📋 planned — Tier 1 #2 of the opportunity scan; GATE-0 table hunt running

## Context

`rmsky` (#22) recovered the Taylor+2009 Faraday-sky results from 37,543 NVSS RMs. SPICE-RACS DR2
(arXiv:2605.16917, 2026-05) provides ~2.5–3.4×10⁵ RMs over 87.5% of the sky — 5× every previous
catalogue combined — released **without** a systematic Galactic structure-function analysis by
latitude/longitude sector. This slice extends `rmsky` to that catalogue: |RM| latitude profile +
enhancement ratio + quadrant sign structure (direct reuse), plus the new piece — the **RM
structure function** SF(δθ) per |b| bin, the coherence-scale measurement the 6.7 deg⁻² density
newly enables. DR1 (`AS110.spice_racs_dr1_corrected_cut_v02`, verified on CASDA TAP) is the
bounded fallback/first leg if the DR2 table is not yet public (GATE 0 decides).

## Deliverables

- `src/jansky_research/rmstructure.py`: `fetch_spice_racs` (TAP, `# pragma`), `structure_function`
  (pair-count SF with error bars via bootstrap; reuse `rmsky.latitude_profile`/`sign_asymmetry`
  wholesale), `synthetic_rm_screen` (injected coherence scale → SF recovers it), Taylor-2009
  cross-anchor on overlap sky, `run/_figure/_write_macros/_main`.
- Tests to the floor; `papers/rmstructure/`; `survey/rmstructure-findings.md`; wiring.

## Approach

0. **GATE 0 (running):** locate the DR2 machine-readable table (CASDA TAP / DAP / VizieR / paper
   DAS); record row count + columns; else scope to DR1 (~thousands of RMs, 1300 deg²) with the
   DR2 upgrade noted as the follow-on.
1. Tooling + synthetic recover-a-known (inject an RM screen with a known coherence scale and a
   known |b| enhancement; SF + profile recover both).
2. Real leg: latitude profile + quadrants (Taylor cross-anchor on overlap) + SF per |b| bin;
   pair-count workload is CPU-friendly (binned pair sums, no GPU).
3. GATE-2 (nπ ambiguity handling, extragalactic RM scatter floor, survey-edge selection) → paper.

## Verification

Round-trip recovers injected coherence scale/enhancement; Taylor+2009 overlap anchor within
errors; checks green; GATE-2 sign-off; `make reproduce` regenerates from the live catalogue.

## Risks & mitigations

- **POSSUM team publishes the SF analysis** → DEFROST (arXiv:2605.13605) exists for foregrounds
  but the |b|-binned SF is unclaimed; ship bounded and fast.
- **DR2 not yet public** → DR1 leg first (still a new SF measurement on that sky), DR2 swap later.
- **Structure-function pitfalls** (uneven sampling, error-bias) → subtract the noise term
  2σ_err², bootstrap errors, verify on the synthetic screen.
