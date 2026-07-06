# 53 — Broadband technosignature EIRP limits from survey Stokes V (empty haystack cell)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — specifically: confirm
no survey-V EIRP-limit paper has landed since the Lenc+2018 lineage, and re-read the narrowband
haystack anchors (arXiv:2103.16250, 2606.04304) full-text

## Context

All published technosignature EIRP-limit work is narrowband (BL GBT, MeerKAT commensal
arXiv:2103.16250, OVRO-LWA narrowband arXiv:2606.04304). Lenc+2018 showed Stokes V isolates
artificial emitters, but nobody has converted survey-V non-detections into broadband
technosignature limits — an empty haystack cell (fable-ideas F16). Data: RACS-low/mid V via
CASDA (verified access pattern, OPAL + `~/.casda_pw`); targets from the Gaia CNS 100 pc sample.
Near-total reuse of the merged `stokesv`/`stokesv_discovery` slices (`measure_circular_pol`,
per-field leakage floor, exclusion list). This merges naturally with F22's RACS-V UCD cutout
sweep (plans/59) — same cutout I/O, so run both target lists in one CASDA pass.

## Deliverables

- `src/jansky_research/eirp.py`: `fetch_target_cutouts` (CASDA, `# pragma`, shared with the
  plans/59 pass), `measure_v_limits` (wraps `measure_circular_pol` at ~10⁴-target scale,
  resumable), `eirp_from_limit` (EIRP = 4πd²·S·Δν from 3σ V limits + Gaia distances),
  `haystack_axes` (place limits on the standard haystack axes vs the narrowband literature),
  `leakage_floor_check`, `synthetic_v_injection`, `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/eirp/`; `survey/eirp-findings.md`; wiring.

## Approach

0. **GATE 0:** standing full-text novelty pass on the closest papers (Lenc+2018 lineage,
   arXiv:2103.16250, 2606.04304); confirm no survey-V EIRP conversion exists; re-verify the
   CASDA cutout path and Gaia CNS table access.
1. Tooling + synthetic recover-a-known: inject V point sources of known flux into real RACS
   cutouts; `measure_v_limits` + `eirp_from_limit` must round-trip flux and EIRP.
2. Real leg: forced V at ~10⁴ Gaia CNS (100 pc) positions in RACS-low/mid; 3σ V limits → EIRP;
   CASDA-I/O-bound, resumable, ~a week wall-clock — executed as one pass with the plans/59
   UCD sweep (same cutout I/O). Cross-check detections against the `stokesv_discovery`
   exclusion list (known V stars are astrophysical, not candidates).
3. GATE-2 science review naming the fable-ideas caveats: broadband-transmitter priors are
   contestable — frame strictly as parameter-space cartography; the per-beam leakage floor
   sets the depth and must be quoted as such.
4. Paper: haystack-cell limits table + figure; RNAAS/arXiv-able.

## Verification

Known V stars and satellite artifacts must appear in the sweep (positive controls, as GJ 65 did
for `stokesv_discovery`); injected-V round-trip recovers flux/EIRP; checks green; GATE-2
sign-off.

## Risks & mitigations

- **Referee pushback on broadband-transmitter priors** → frame strictly as parameter-space
  cartography (which haystack cell is now bounded, at what depth), never as an occurrence claim.
- **Leakage sets the floor** → carry the `stokesv` per-field leakage model through to the EIRP
  axis; quote leakage-limited vs noise-limited targets separately.
- **Likely all-limits outcome** → the limit table *is* the paper; say so from day one.
