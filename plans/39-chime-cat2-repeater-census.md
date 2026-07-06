# 39 — CHIME/FRB Catalog 2: uniform repeater wait-time & duty-cycle census

Status: ✅ done — GATE 0 passed 2026-07-06 (chime-frb.ca still 503, but the official CANFAR DOI
10.11570/25.0066 serves the table + exposure maps — mirrored to `data/`; gap confirmed open:
Cook et al. arXiv:2605.08410 does rates only). Plan corrections: 20240114A is NOT in Cat 2
(post-cutoff; its "112.9 d" is chromatic anyway) — 20180916B is the single in-catalogue anchor,
recovered at 16.33 d (p=0.001, 107 cycles, duty 0.21); the public exposure product is
time-integrated → the null is a transit-comb-preserving sidereal scramble, not per-epoch
correction. Census: 83 sources, 15 above the cut, median k=0.83, 3 clustered; the two non-anchor
p≤0.01 peaks have ≤5 cycles and are labelled activity-epoch degeneracies, not periods. See
survey/frbwait-findings.md.

## Context

CHIME/FRB Catalog 2 (arXiv:2601.09399; 4,539 bursts, 83 repeaters) has had its obvious
population angles closed — injection debiasing (arXiv:2606.26334), 80-repeater uniform rates +
DM drift (arXiv:2605.08410), spectral split (arXiv:2601.16048), ML repeater classification
(arXiv:2509.02645, 2512.08308), polarization dichotomy (arXiv:2401.17378 + companion); the
Corrections section of fable-ideas.md fences all of these off. What survives: no paper computes
one uniform statistic — Weibull clustering k + exposure-corrected activity-window/duty-cycle —
across all ~83 repeaters from the public table. All periodicity results are single-source
campaigns (20240209A arXiv:2502.11215; 20240114A ~112.9 d; 20220912A arXiv:2604.09098). Extends
the merged `frbstats` (Weibull) and `frbperiod` (Rayleigh Z², LS, scramble FAP) slices, ~80%
reuse; the Cat-1 slice's honest lesson (`survey/period-findings.md`: exposure-blind FAPs are not
rigorous) is the design driver here.

## Deliverables

- `src/jansky_research/frbwait.py`: `fetch_catalog2` + `fetch_exposure_maps` (# pragma),
  `per_repeater_weibull_k` (extends `frbstats`), `exposure_corrected_periodogram` (extends
  `frbperiod` with the catalogue's own exposure maps), `population_k_distribution`
  (hierarchical pass with a stated burst-count completeness cut), `synthetic_repeater_set`
  (injected k + injected period under a real exposure function), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/frbwait/`; `survey/frbwait-findings.md`; wiring.

## Approach

0. GATE 0: full-text pass on arXiv:2601.09399, 2605.08410, 2601.16048, 2606.26334 to confirm no
   uniform wait-time/duty-cycle census exists; confirm `chime-frb.ca/catalog2` has recovered
   from its 503s and mirror the CSV/FITS tables + exposure/sensitivity maps into `data/`
   immediately (no auth). If the site stays down, the slice waits — do not substitute.
1. Tooling + synthetic recover-a-known: inject Weibull-clustered burst trains with known k and a
   known period through a real exposure function; recover both.
2. Real leg (CPU, days–weeks): per-repeater k, exposure-corrected periodogram, population
   k-distribution above the completeness cut; anchor by reproducing 20180916B's 16.35 d and
   20240114A's ~112.9 d before trusting the other ~81 repeaters.
3. GATE-2 science review: exposure-map fidelity, burst-count completeness cut, low-N repeater
   k posteriors, transit-cadence aliasing in periodograms.
4. Paper: the first uniform ~83-repeater wait-time/duty-cycle census.

## Verification

Pipeline reproduces both published anchors from the same dataset — 20180916B's 16.35 d AND
20240114A's ~112.9 d — plus the synthetic k/period round-trip under a real exposure function;
checks green; GATE-2 sign-off.

## Risks & mitigations

- **Highest scoop risk in fable-ideas.md** — CHIME/FRB is the natural author. Move fast or not
  at all: mirror data the day the site recovers, ship bounded.
- **Exposure-blind FAPs are not rigorous** (Cat-1 lesson) → the catalogue's own exposure maps
  are mandatory inputs, not optional refinements.
- **Low-burst-count repeaters give unconstrained k** → stated completeness cut; report
  population results only above it, per-source posteriors honestly wide below it.
