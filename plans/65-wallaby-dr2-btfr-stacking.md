# 65 — WALLABY DR2 pair: BTFR/angular momentum vs environment + cube stacking at DESI positions

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — fetch the kinematic
products via CSIRO DAP DOI 10.25919/7w8n-9h19 (team page, not obscore) + check DESI footprint
overlap before committing to leg (b)

## Context

WALLABY DR2's 126 kinematic models remain the ceiling — no DR3 until ~2029, and the team's
post-release papers are morphometrics and dark-source work, not BTFR-vs-environment (fable-ideas
F28, open as of the 2026-07 scan). Two legs: (a) baryonic Tully–Fisher and Fall-relation
(angular-momentum) residuals split by field/local-density environment at N=126 — small enough
that the null framing is pre-registered up front; (b) sub-threshold spectral stacking of WALLABY
cubes (~500 MB/field) at DESI-redshift positions, velocity-aligned, adapting the merged
`stacking` injection-recovery harness from continuum images to spectral cubes. Data: kinematic
model products via the CSIRO DAP team page (DOI 10.25919/7w8n-9h19 — obscore does not expose
them), WALLABY DR2 cubes via CASDA (verified auth), DESI DR1 redshifts + group/density catalogues
for the environment axis. CPU throughout; the cube I/O is the only heavy part.

## Deliverables

- `src/jansky_research/wallaby.py`: `fetch_kinematic_models` (DAP team page, `# pragma`),
  `fetch_cubes` (CASDA, `# pragma`), `btfr_fit` (baryonic mass vs rotation velocity with
  orthogonal-scatter errors), `fall_relation` (specific angular momentum vs mass),
  `environment_split` (DESI local density / group membership terciles + residual tests),
  `cube_stack_at_positions` (velocity-aligned sub-threshold spectral stack at DESI-z positions),
  `injection_recovery_cubes` (planted HI profiles in real cubes → stack recovery, `stacking`
  harness adapted), `synthetic_btfr` (mock relation + planted environment offset → recovery),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/wallaby/`; `survey/wallaby-findings.md`; wiring.

## Approach

0. GATE 0: fetch the kinematic products from the DAP team page (DOI 10.25919/7w8n-9h19) and
   confirm the 126-model contents; check DESI DR1 footprint overlap with the DR2 fields (this
   decides whether leg (b) proceeds); full-text pass on the WALLABY DR2 + follow-up papers to
   confirm BTFR-vs-environment and DESI-position stacking remain unclaimed.
1. Tooling + synthetic recover-a-known: mock BTFR with a planted environment offset →
   `environment_split` recovers it at N=126 (this also calibrates the detectable effect size);
   planted HI profiles in real cubes → stack recovery at stated completeness.
2. Real leg (a): BTFR + Fall-relation fits, environment split, residual tests; report the
   minimum detectable offset alongside any result.
3. Real leg (b): velocity-aligned stacking at DESI positions below the WALLABY detection
   threshold, per environment bin; recover catalogued WALLABY sources in the stack first.
4. GATE-2 science review: N=126 statistical honesty (pre-registered null framing), distance and
   inclination systematics in the BTFR, stacking continuum-residual + velocity-error caveats.
5. Paper: the environment tests + the stacked detections/limits, nulls framed as deliverables.

## Verification

Cube stacking must recover catalogued WALLABY sources before any sub-threshold claim; mock BTFR
environment offset recovered at the calibrated effect size; injection-recovery completeness
stated; checks green; GATE-2 sign-off.

## Risks & mitigations

- **N=126 ceiling** → pre-register the null framing: the paper reports the minimum detectable
  environment offset and where the data land, not a forced detection; no DR3 relief until ~2029,
  so this is the definitive WALLABY-DR2-era statement either way.
- **DESI overlap may be thin** → GATE-0 checks footprint overlap before leg (b) is committed;
  leg (a) stands alone if stacking positions are too few.
- **Stacking systematics** (continuum residuals, DESI velocity errors) → the adapted
  injection-recovery harness quantifies both before any stacked signal is claimed.
