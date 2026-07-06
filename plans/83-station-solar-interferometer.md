# 83 — Two-dish additive interferometer: the Sun's angular diameter at 1420 MHz

Status: 📋 planned (hardware-gated) — GATE-0 is ANALYTIC and precedes any hardware purchase:
compute the baseline required to resolve the ~0.5° solar disc on this roof before buying the
second dish/LNA/combiner; sequenced last of the station slices per `station/interferometry.md`

## Context

The additive (phase-switched) interferometer in `station/interferometry.md` is the plan of
record for dish B: two filtered front ends into a combiner and one SDR; fringes as the Sun
drifts through the overlapped beams give the angular size — 1950s radio astronomy for ~$200.
The station sweep found amateur two-element write-ups stop at "we saw fringes"; a solar-diameter
number with an honest fringe-spacing error budget is the rigor axis that makes it citable. The
physical caveat shapes everything: a ~0.5° disc nulls the visibility at modest baselines, so the
required baseline must be shown to fit the roof *analytically before purchase* — a gate that can
kill the slice cheaply. Sequence last: dish B may first serve plan 79 as a Dicke-switched load.

## Deliverables

- `src/jansky_research/solarfringe.py`: `disc_visibility` (uniform/limb-darkened disc visibility
  vs baseline at 21 cm), `required_baseline` (target fringe contrast → baseline range + roof
  feasibility check), `simulate_drift_fringes` (beam envelope × fringes for a given geometry),
  `fit_fringe_spacing` (fringe rate → spacing with covariance), `fit_diameter` (spacing +
  visibility model → angular diameter + error budget), `run`/`_figure`/macros.
- Tests to the 85% floor (simulated fringes — no sky data needed for tests);
  `papers/solarfringe/`; `survey/solarfringe-findings.md`; wiring.

## Approach

0. GATE 0 (analytic, BEFORE hardware purchase): with `disc_visibility`/`required_baseline`,
   derive the east–west baseline window that keeps fringe contrast measurable for a ~0.5° disc
   at 1420 MHz, versus the roof's usable extent; include combiner/cable-mismatch loss terms.
   If the window does not fit the roof, the slice dies here for the cost of a notebook.
1. Tooling + synthetic recover-a-known (still pre-hardware): `simulate_drift_fringes` with a
   known diameter/baseline → `fit_fringe_spacing` + `fit_diameter` round-trip within stated
   errors, including baseline-length and pointing uncertainties in the budget.
2. Hardware leg (only if GATE 0 passes): second dish + LNA + combiner + matched-batch filters
   (already purchased as a pair) and matched cables; commission per `station/test-equipment.md`.
3. Observation leg (`# pragma: no cover`): solar drift scans across several days/baselines;
   fit diameter per scan; compare to the known 21 cm radio Sun (larger than optical).
4. GATE-2 science review → paper: the diameter with its full fringe-spacing error budget.

## Verification

Simulated fringes round-trip the injected diameter within the stated budget; the real-data
anchor is an honest fringe-spacing error budget (baseline, timing, beam envelope, RFI) on the
fitted diameter; checks green; GATE-2 sign-off.

## Risks & mitigations

- **The required baseline may be impractical on this roof** → that is what the analytic GATE-0
  exists for; it kills the slice before a dollar is spent, and the feasibility note still lands
  in `survey/solarfringe-findings.md`.
- **Hardware delay / dish-B contention with plan 79's Dicke switch** → sequenced last by design;
  steps 0–1 are hardware-free and ship regardless.
- **Gain drift between the two arms mimics fringe-contrast loss** → matched-batch filters,
  matched cables, and per-scan contrast checks against the simulated envelope.
