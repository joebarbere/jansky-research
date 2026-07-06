# 38 — The first rotation-measure dipole/anisotropy test (SPICE-RACS DR2)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + same-week ADS dipole
re-search (the fable-ideas scan ran egress-blocked; see the standing caveat there); data risk ≈ 0
(DR2 table already on disk from `rmstructure`)

## Context

The cosmic-dipole-anomaly literature is entirely source-count/flux based: the joint
NVSS+RACS+LoTSS-DR2 5.4σ result (arXiv:2509.16732, PRL), arXiv:2509.18689, the RMP colloquium
(arXiv:2505.23526), and Mittal & Lewis (arXiv:2605.27520) — none uses Faraday rotation measures
as the tracer. fable-ideas.md flags total-intensity source counts as a crowded/closed axis
(Corrections section); the RM dipole is the open differentiator. SPICE-RACS DR2
(arXiv:2605.16917; ~2.5–3.4×10⁵ RMs over 87.5% of sky) is the first catalogue large enough, and
the merged `rmstructure` slice already fetched it (CSIRO DAP `csiro:64891`, ~5 GB, no auth) and
built the per-|b| Galactic-floor machinery this slice reuses (~70%). Likely outcome is an
isotropy test with honest nulls — no clean kinematic-expectation amplitude exists for RM, so a
positive is hard to interpret and the plan does not pretend otherwise.

## Deliverables

- `src/jansky_research/rmdipole.py`: `extragalactic_residuals` (per-source RM residuals at
  |b|>30–45° via `rmstructure` Galactic-floor subtraction, DEFROST regime arXiv:2605.13605),
  `healpix_binned_map`, `fit_dipole` (likelihood/Bayesian dipole in mean |RM| and RM²),
  `footprint_scramble_null` (RA scrambles preserving the Dec≤+49° footprint),
  `synthetic_dipole_catalogue` (injected dipole on the real footprint/noise),
  `compare_directions` (CMB dipole + Böhme+ source-count dipole), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/rmdipole/`; `survey/rmdipole-findings.md`; wiring.

## Approach

0. GATE 0: full-text pass on arXiv:2509.16732, 2509.18689, 2505.23526, 2605.27520 and a
   same-week ADS re-search (the dipole field moves fast) to confirm no RM-dipole paper landed;
   re-verify the DR2 table on disk matches DAP `csiro:64891`.
1. Tooling + synthetic recover-a-known: inject a known dipole (amplitude + direction) into a
   mock all-sky RM catalogue with the real DR2 footprint and noise; confirm recovery (mirrors
   `rmstructure`'s synthetic-screen validation). Pure NumPy/healpy, CPU, days.
2. Real leg: residuals → HEALPix bins → dipole fit in mean |RM| and RM²; null distribution from
   footprint-preserving RA scrambles (the footprint is the dominant systematic — the scramble
   test is load-bearing); compare direction to CMB and source-count dipoles.
3. GATE-2 science review: interpretation limits of a positive (no kinematic-expectation
   amplitude), footprint/Galactic-floor leakage, nπ-ambiguity inheritance from DR2.
4. Paper: framed as the first RM isotropy/dipole test, nulls reported plainly.

## Verification

Synthetic catalogue round-trip recovers injected dipole amplitude and direction on the real
footprint; footprint-scramble null distribution is consistent with zero when no dipole injected;
checks green; GATE-2 sign-off.

## Risks & mitigations

- **No kinematic-expectation amplitude for RM** → frame as an isotropy test with honest nulls;
  any positive is reported as an anisotropy, not a velocity measurement.
- **Dipole field moves fast** → GATE-0 includes a same-week ADS re-search; ship bounded and fast.
- **Footprint systematics (Dec≤+49°)** → the footprint-preserving scramble null is mandatory and
  validated on the synthetic catalogue before the real fit is trusted.
