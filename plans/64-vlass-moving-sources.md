# 64 — Blind moving-source search across VLASS E1/E2/E3: radio proper motions without Gaia

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — confirm the VLASS
epoch catalogue versions/astrometric docs at CIRADA and re-check for any blind radio-PM paper

## Context

Every published radio proper motion starts from a known object — the field is Gaia-anchored
(arXiv:2409.18466); no blind, radio-only moving-source search across the three VLASS epochs
exists (fable-ideas F27). Method: GPU all-pairs linkage between the E1 and E2 component
catalogues in the 0.3–5″/yr annulus, E3 collinearity confirmation (three epochs on a line, with
consistent rate), flux/morphology consistency cuts — and only *then* exclusion of Gaia/WISE
counterparts. The survivors are optically-dark movers (Y dwarfs, high-proper-motion pulsars) or,
at zero yield, the first blind surface-density limit on optically-dark radio movers. Data: VLASS
epoch component catalogues via CIRADA/CADC (no auth). Fit: GPU all-pairs (pure torch, the
`dr3ghosts` chunking pattern), weeks-long runs explicitly OK on this workstation. The published
VLASS astrometric floor (0.25–0.5″) sets the honest lower rate bound and the false-pair budget.

## Deliverables

- `src/jansky_research/vlasspm.py`: `fetch_vlass_epochs` (CIRADA component catalogues,
  `# pragma`), `epoch_pair_linkage` (GPU-chunked E1×E2 all-pairs in the 0.3–5″/yr annulus),
  `collinearity_test` (E3 three-epoch line + rate consistency), `flux_morphology_cuts`
  (compactness + epoch-to-epoch flux ratio), `counterpart_exclusion` (Gaia DR3 + CatWISE
  cross-match, applied last, `# pragma`), `scramble_false_pair_rate` (position-scrambled epoch
  null → expected chance-pair count), `surface_density_limit` (yield → limit with completeness
  from injections), `inject_movers` (planted PM sources → end-to-end recovery),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/vlasspm/`; `survey/vlasspm-findings.md`; wiring.

## Approach

0. GATE 0: full-text pass on arXiv:2409.18466 and the VLASS epoch/astrometry papers (pin the
   per-epoch astrometric floor numbers); ADS check that no blind radio-PM search has landed;
   verify CIRADA catalogue URLs and epoch coverage overlap.
1. Tooling + synthetic recover-a-known: inject movers at 0.3–5″/yr into mock epoch catalogues
   with the real astrometric floor; linkage + collinearity must recover them and the scrambled
   null must predict the chance-pair count.
2. Real leg: E1×E2 linkage (multi-day GPU job, checkpointed), E3 collinearity, flux/morphology
   cuts, then Gaia/WISE exclusion; per-candidate vetting sheet; completeness from injections.
3. GATE-2 science review: false-pair budget honesty (variable AGN pairs), astrometric-floor
   propagation into the rate annulus, candidate-vs-limit framing discipline.
4. Paper: candidates (if any survive vetting) or the first blind surface-density limit on
   optically-dark radio movers — the null is pre-framed as a deliverable.

## Verification

UV Ceti (3.4″/yr) must fall out of the blind pipeline before the counterpart-exclusion step;
injected movers recovered at stated completeness; scramble null matches observed chance-pair
count; checks green; GATE-2 sign-off.

## Risks & mitigations

- **VLASS astrometric floor (0.25–0.5″)** → the 0.3″/yr lower rate bound sits at the floor;
  propagate per-epoch astrometric errors into the annulus and report completeness vs rate
  honestly rather than claiming the nominal range.
- **Variable AGN creating false E1–E2 pairs** → the E3 collinearity test and the
  position-scramble false-pair rate are load-bearing; no candidate survives on two epochs alone.
- **Likely zero survivors** → the surface-density limit is the paper; pre-register that framing.
