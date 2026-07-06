# 41 — Radio counterpart survey of the 56 optically-selected white-dwarf-pulsar candidates

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + locating the
machine-readable 56-row candidate table (the fable-ideas scan ran egress-blocked; see the
standing caveat there)

## Context

Pelisoli et al. 2025/26 (MNRAS 540, 821; arXiv:2505.04693) published 56 AR Sco-like white-dwarf
-pulsar candidates from Gaia+WISE light curves, 26 previously uncharacterized. No systematic
radio search of the list exists — the one confirmed WD pulsar found via this route, J1912−4410
(arXiv:2306.09272), was followed up object-by-object. A detection is a new white-dwarf pulsar;
the non-detection table is honest value regardless, and the plan treats the limit table as a
first-class deliverable. Data: the candidate table from the paper's supplementary/VizieR
(GATE-0), RACS-low/mid I+V via CASDA (verified auth pattern), and VLASS QL2 via CADC/CIRADA
(no auth) for Dec>−40 targets. Near-total reuse of the merged `stokesv`/`stokesv_discovery`
machinery (`measure_circular_pol`, `classify_emitter`, the resumable M-dwarf-run pattern).

## Deliverables

- `src/jansky_research/wdpulsar.py`: `load_candidate_table` (56 rows, provenance-flagged),
  `fetch_racs_epochs` + `fetch_vlass_ql2` (# pragma), `forced_iv_photometry` (per epoch per
  candidate via `measure_circular_pol`), `vet_leakage` (via `classify_emitter`),
  `vlass_cone_check`, `limits_table` (detections, |V|/I, upper limits), `synthetic_candidate_field`
  (injected V source into a real cutout → recovered), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/wdpulsar/`; `survey/wdpulsar-findings.md`; wiring.

## Approach

0. GATE 0: full-text pass on arXiv:2505.04693 and arXiv:2306.09272 (plus an ADS search for any
   radio follow-up of the list) to confirm no systematic radio survey exists; locate and mirror
   the machine-readable 56-row table (supplementary/VizieR) — access is unverified.
1. Tooling + synthetic recover-a-known: inject a V-bright point source at a candidate-like
   position into a real RACS cutout; forced photometry + leakage vetting recover it.
2. Real leg (CPU + CASDA I/O, ~a week wall-clock, resumable like the M-dwarf run): forced I/V
   photometry per epoch per candidate across RACS-low/mid, leakage vetting, VLASS QL2 cone
   checks for Dec>−40; assemble the detection/limit table.
3. GATE-2 science review: leakage-floor caveats on faint |V|/I claims, epoch-dependent duty
   cycles (a non-detection is not an absence), VLASS/RACS frequency mismatch in limit framing.
4. Paper: detections (if any) + the systematic 56-target radio limit table.

## Verification

J1912−4410 is in the candidate list and is radio-detected — the pipeline must re-find it (same
pattern as the GJ 65 control in `stokesv_discovery`); synthetic injection round-trip passes;
checks green; GATE-2 sign-off.

## Risks & mitigations

- **Pelisoli's team is the natural competitor** (may hold MeerKAT time) → move promptly once the
  table is in hand; the RACS/VLASS archival angle is complementary to targeted MeerKAT work.
- **Table access unverified** → GATE-0 blocks the slice until the 56-row table is mirrored.
- **Likely mostly non-detections at survey depth** → the limit table is framed as the deliverable
  from day one, not spun as a near-miss.
