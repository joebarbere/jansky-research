# 42 — Lensed-repeater test: recurring time-delay patterns in CHIME Catalog 2

Status: ✅ done — GATE 0 passed 2026-07-06 (data via CANFAR DOI 10.11570/25.0066, shared mirror
with plan 39; arXiv:2605.19653 checked: ms intra-burst microlensing on one-offs — different
regime, cited as fence). Three design deltas the real data + GATE-2 forced: (1) the plan's day-scramble null
false-positived the three most clustered repeaters → phase-permutation null + 5-s tolerance;
(2) mjd_400 is topocentric → all TOAs barycentred (Roemer drift ±150 s at Δ=26 d ≫ tolerance);
(3) the 1 pc/cc DM floor would reject genuine image pairs → per-source 3√2·σ_DM cut, costed in
the injections. All documented as methods lessons. Result: 0 detections / 33 searched (cleanest
p=0.81) → first empirical lensed-repeater fraction limit < 0.091 per searched repeater (95%),
injection-mapped sensitivity (transit comb selection made explicit). See
survey/frblens-findings.md.

## Context

Theory fully specifies how a strongly-lensed one-off FRB masquerades as a "repeater" with a
fixed pattern of mutual delays (Dai & Lu 2017, 10.3847/1538-4357/aa8873; Li+2018, Nat. Comm.
9:3833; review arXiv:2412.01536) — but no observational catalogue-level search has been
published. CHIME's own lensing searches are intra-burst/microsecond baseband work (Leung+2022 /
Kader+2022 — a different regime, listed as a closed door in fable-ideas.md's Corrections for
sub-second echoes; this slice is the burst-to-burst regime they did not cover). Cat 2
(arXiv:2601.09399) is the first dataset big enough. Data is the same as plan 39 (per-burst TOAs,
DMs, morphologies, exposure functions) — bundle the mirror. Expected yield ≈ 0 (lensing optical
depth ~10⁻⁴): the paper *is* the first empirical upper limit on the lensed fraction, framed that
way from day one. Reuses `frbperiod`'s scramble/FAP machinery (~60%).

## Deliverables

- `src/jansky_research/frblens.py`: `all_pairs_delays` (per-repeater TOA delay histogram,
  GPU-trivial torch reduction), `recurring_delay_search`, `consistency_cuts` (DM agreement at
  the measurement floor + morphology consistency), `exposure_aware_scramble_fap` (CHIME transit
  cadence aliases — mandatory), `inject_lensed_pairs` (maps the detectable
  delay/magnification-ratio region), `lensed_fraction_limit`, `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/frblens/`; `survey/frblens-findings.md`; wiring.

## Approach

0. GATE 0: full-text pass on Dai & Lu 2017, Li+2018, arXiv:2412.01536, and an ADS search
   confirming no catalogue-level lensed-repeater search has appeared; data gate shared with
   plan 39 — confirm chime-frb.ca recovery and mirror Cat 2 tables + exposure functions into
   `data/` in the same pass.
1. Tooling + synthetic recover-a-known: inject lensed pairs (fixed delay, DM-matched,
   magnification-ratio scaled) into realistic burst trains under the real exposure function;
   recover at the stated FAP; map the sensitivity region in (delay, magnification ratio).
2. Real leg (CPU/GPU, days): per repeater, all-pairs delay histograms → recurring-delay
   candidates → DM + morphology consistency cuts → exposure-aware TOA-scramble FAPs; convert
   the (expected) null into the first lensed-fraction upper limit via the injection map.
3. GATE-2 science review: transit-cadence aliasing (the scramble is load-bearing), DM
   measurement-floor choice, morphology-consistency subjectivity, limit-statement scope.
4. Paper: the first empirical catalogue-level upper limit on the lensed-repeater fraction.

## Verification

Injected lensed pairs recovered at stated FAP; known non-lensed high-count repeaters yield clean
nulls under the exposure-aware scramble; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Expected yield ≈ 0 (optical depth ~10⁻⁴)** → the upper limit is the paper; framed as such
  from day one, never as a failed detection.
- **CHIME transit cadence aliases delay histograms** → exposure-aware scrambles are mandatory
  and validated on synthetics before any real FAP is quoted.
- **Data gate shared with plan 39** → if chime-frb.ca stays down, both slices wait; mirror once,
  run both.
