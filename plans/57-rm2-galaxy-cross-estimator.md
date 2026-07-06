# 57 — First real-data application of the ⟨RM²×g⟩ cross-correlation estimator (SPICE-RACS DR2)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — specifically: read
Zhang & Lidz (arXiv:2512.06584) full-text — it is an SKA forecast only — and confirm no
real-data application of the estimator has appeared since

## Context

Zhang & Lidz (arXiv:2512.06584) define an ⟨RM²×g⟩ estimator — cross-correlating rotation-measure
variance with galaxy overdensity to probe magnetized large-scale-structure gas — but published
it as an **SKA forecast only**; no real-data application exists (fable-ideas F20). SPICE-RACS
DR2 (arXiv:2605.16917; CSIRO DAP `csiro:64891`, on disk from the merged `rmstructure` slice)
has ~30× fewer RMs than the forecast assumes, so the honest headline is likely a first **upper
limit** — said up front. Hardest of the three plan-38 siblings (with plans/55, 56); a stretch
goal — sequence last, after plan 38, reusing the same data and the `rmstructure`
Galactic-floor/residual machinery. Fence from Corrections: the total-intensity dipole space is
crowded — RM-based statistics are the open differentiators this cluster owns.

## Deliverables

- `src/jansky_research/rmcross.py`: `rm2_residual_map` (extragalactic RM² residuals after the
  `rmstructure` per-|b| Galactic-floor subtraction, HEALPix), `galaxy_overdensity_map` (public
  galaxy catalogue fetch + masking, `# pragma`), `cross_estimator` (the ⟨RM²×g⟩ statistic per
  arXiv:2512.06584, angular bins), `mock_injection` (correlated mock RM²/galaxy fields on the
  real footprint), `scramble_covariance` (RA scrambles + jackknife),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/rmcross/`; `survey/rmcross-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text read of arXiv:2512.06584 (exact estimator definition, forecast
   assumptions, the 30× RM-count gap quantified against DR2); confirm no real-data application
   landed; pick + verify the galaxy catalogue (DESI/WISE-class, DR2-footprint overlap);
   standing full-text novelty pass.
1. Tooling + recover-a-known on mocks: correlated mock RM²/galaxy-overdensity fields with a
   known cross-amplitude on the real DR2 footprint/noise; the estimator must recover the
   injection with calibrated errors — the load-bearing validation, since no real-data anchor
   exists for a forecast-only estimator.
2. Real leg: DR2 RM² residuals × galaxy overdensity; covariance from footprint-preserving RA
   scrambles + jackknife; report a measurement or a first upper limit with the injection-based
   sensitivity curve alongside it.
3. GATE-2 science review naming the fable-ideas caveats: likely an upper limit at DR2 size
   (30× fewer RMs than the forecast) — framed as such from the abstract down; Galactic-residual
   leakage into RM² is the dominant systematic.
4. Paper: first real-data application of the estimator — method port + limit + the honest
   statement of what catalogue size makes it a detection.

## Verification

Injected-signal recovery on mocks with the real footprint/noise (the estimator is forecast-only,
arXiv:2512.06584, so mocks are the anchor); scramble nulls consistent with zero; checks green;
GATE-2 sign-off.

## Risks & mitigations

- **Likely an upper limit at DR2 size (30× fewer RMs than the forecast)** → say so up front;
  the deliverable is the first real-data application + sensitivity curve, not a detection.
- **Hardest of the three siblings; stretch goal** → sequence after plans 38/55/56, timebox the
  real-data leg; the mock-validated tooling stands alone if the limit is uninformative.
- **Galactic-RM residuals correlate with masks/extinction and fake a signal** → `rmstructure`
  floor subtraction + extinction-map deprojection + scramble nulls.
