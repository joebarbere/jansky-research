# 56 — CGM/halo magnetization: RM² stacking behind DESI group/cluster halos (SPICE-RACS DR2)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — specifically: read
arXiv:2605.16924 full-text (the Mg II angle is CLOSED there), confirm the group/cluster-halo
tracer variant remains unclaimed, and verify DESI DR1 group-catalogue access

## Context

CGM magnetization via RM stacking on **Mg II absorbers is closed** — Van Eck/Malik et al., A&A
2026, arXiv:2605.16924 (SPICE-RACS DR2, 612 sightlines; the Corrections fence). What survives
(fable-ideas F19) is the **halo-mass/group tracer variant**: stacking RM² behind DESI DR1
group/cluster halos by richness, asking whether magnetized halo gas shows up when selected by
halo mass rather than by cool-gas absorption. Data: the SPICE-RACS DR2 RM table (CSIRO DAP
`csiro:64891`, already on disk from the merged `rmstructure` slice) + a public DESI DR1
group/cluster catalogue. This is a sibling of plan 38 (F1, RM dipole) — sequence it after that
slice; same catalogue, and the `rmstructure` Galactic-floor machinery is the foreground layer.

## Deliverables

- `src/jansky_research/cgmrm.py`: `load_groups` (DESI DR1 group/cluster table, `# pragma`),
  `annulus_grm_subtract` (re-implementation of the 2605.16924 annulus Galactic-RM subtraction),
  `stack_rm2` (excess RM² vs impact parameter behind halos, bootstrap errors),
  `richness_bins`, `mgii_anchor` (reproduce their Mg II excess from their sightline selection),
  `synthetic_halo_injection`, `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/cgmrm/`; `survey/cgmrm-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text read of arXiv:2605.16924 (method to re-implement; confirm the
   group-tracer variant is not in it or in a companion); locate the machine-readable DESI DR1
   group/cluster catalogue and check footprint overlap with DR2; standing novelty pass.
1. Tooling + synthetic recover-a-known: inject a halo RM² excess of known amplitude and profile
   into mock sightlines with realistic Galactic foreground + measurement noise; the annulus
   subtraction + stack must recover the amplitude vs impact parameter.
2. Recover-a-known on real data: reproduce the Mg II excess of arXiv:2605.16924 using their
   annulus-GRM subtraction on the same DR2 table — the pipeline is not trusted until it does.
3. Real leg: stack RM² behind DESI DR1 groups/clusters in richness × impact-parameter bins;
   quote per-bin sightline counts; bootstrap + control-stack (random-position) errors.
4. GATE-2 science review naming the fable-ideas caveat: small N behind massive groups — a null
   must not be spun as contradicting the Mg II detection (different tracer, different gas
   phase); pre-register the null framing.
5. Paper: first halo-mass-selected RM² stacking measurement (or honest bounded null).

## Verification

The Mg II excess of arXiv:2605.16924 reproduced via their annulus-GRM subtraction; injected
synthetic halo signal recovered in the offline round-trip; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Small N behind massive groups** → quote per-richness-bin counts up front; pre-registered
  null framing; do not spin a null as contradiction of the Mg II result.
- **Galactic-RM residuals leak into the stack** → their annulus subtraction re-implemented and
  anchored, plus |b| cuts and the random-position control stacks.
- **Group-catalogue impurity/miscentring** → report richness-cut sensitivity; treat miscentring
  as a stated systematic, not a fit parameter, at this sample size.
