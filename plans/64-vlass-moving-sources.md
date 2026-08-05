# 64 — Blind moving-source search across VLASS E1/E2/E3: radio proper motions without Gaia

Status: 📋 planned — **GATE 0 discharged 2026-08-04** (full-text novelty pass + data
verification): the blind radio-PM niche is **still open** — the field remains entirely
Gaia/optical-anchored (De et al. 2024 arXiv:2409.18466 and Driessen et al. 2023
arXiv:2306.08059 both anchor to known optical PMs; the E1–E2 variability census
arXiv:2508.00976 is flux-only; the July 2026 VLASS-completion press names transients but no
PM program, with catalog processing "continuing over the next several years"). Data pinned:
E1 = CIRADA VLASS1QLv3.1 (2.42M components, VizieR J/ApJS/255/30), E2 = VLASS2QLv2 CSV
(2.37M), E3 = QL3.1/QL3.2 interim lists per **VLASS Memo 22** (Lacy & Dong 2025-06-03;
combined 2.38M + a 3.38M median-stack list); no SE/"final" catalog supersedes QL yet; epoch
4.1 observed (Feb 2026) but uncatalogued. **One correction to carry:** replace the flat
"0.25–0.5″ astrometric floor" with per-epoch, declination-dependent floors — E1 ~0.5″
(Dec>−20°) to ~1″ (south), E3 ≈0.1″ MAD vs the VLBI RFC per Memo 22 — E1 is the limiting
epoch and drives the false-pair budget and the annulus's honest lower rate bound (cite
Gordon et al. 2021 + Memo 22, not a blanket number). Baselines confirmed: E1–E2 ≈2.8–3 yr,
E1–E3 ≈5.5–6 yr. Remaining standing item: same-week ADS re-search before first commit.

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
