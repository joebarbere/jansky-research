# 59 — RACS Stokes-V two-epoch census of the Gaia UCD sample

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — confirm the Gaia UCD
table used by arXiv:2506.21169 is machine-readable and that no RACS-V UCD paper has landed

## Context

The blind VLASS Stokes-I UCD search is CLOSED: arXiv:2506.21169 checked 14,915 Gaia
ultracool dwarfs against 3 VLASS epochs and found zero brown-dwarf counterparts. What survives
(fable-ideas F22) is the V-selected, lower-frequency two-epoch companion: UCD radio emission is
coherent and highly circularly polarized (new VLITE detection arXiv:2512.11120; the Driessen
review arXiv:2606.27706 flags exactly this lever), so forced Stokes-V photometry in RACS-low1 and
RACS-low2 probes a different emission mechanism at a different frequency than the closed search.
Data: RACS V cutouts via CASDA (verified auth + SODA pattern). Tooling: `stokesv.py`'s
`measure_circular_pol` + leakage vetting, near-total reuse. This slice merges naturally with the
F16 EIRP sweep (plans/53) — same cutout I/O, so the two target lists run in one CASDA pass.
Honest framing: likely near-zero detections — the two-epoch V limit table is the paper.

## Deliverables

- `src/jansky_research/ucdv.py`: `fetch_gaia_ucds` (the 2506.21169 parent sample, `# pragma`),
  `forced_v_photometry` (per epoch per position via `stokesv.measure_circular_pol`),
  `vet_leakage` (per-field leakage floor via `classify_emitter` reuse), `two_epoch_compare`
  (detections, handedness, variability), `limit_table` (3σ V limits + fractional-pol limits),
  `synthetic_injection` (V sources injected into real cutouts → recovery),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/ucdv/`; `survey/ucdv-findings.md`; wiring.

## Approach

0. GATE 0: full-text pass on arXiv:2506.21169 (confirm the VLASS closure and grab the UCD table),
   2512.11120, and 2606.27706; ADS check that no RACS-V UCD census exists; verify CASDA V-cutout
   access for both low1 and low2 epochs; coordinate the shared cutout pass with plans/53.
1. Tooling + synthetic recover-a-known: inject V point sources into real RACS cutouts at UCD-like
   flux levels; `measure_circular_pol` recovery + leakage-floor false-positive rate.
2. Real leg: forced V photometry at all UCD positions in both epochs (CASDA-I/O-bound, resumable,
   ~a week wall-clock, bundled with the F16 sweep); leakage vetting; detection/limit table.
3. GATE-2 science review: leakage-floor honesty, duty-cycle caveat (snapshot epochs miss bursty
   emitters), non-detection framing against the coherent-emission luminosity expectations.
4. Paper: the two-epoch V census + limit table, framed as the companion to the VLASS-I null.

## Verification

Pipeline must re-detect WISE J0623 (arXiv:2306.15219) and the Pritchard RACS UCD detections
before any new claim; injected-source recovery at stated completeness; checks green; GATE-2
sign-off.

## Risks & mitigations

- **Likely near-zero detections** → the limit table is the paper; pre-commit to that framing —
  a uniform two-epoch V census of 14,915 UCDs is citable regardless of yield.
- **Leakage masquerading as V detections** → per-field leakage floor + `classify_emitter` vetting
  (the merged `stokesv_discovery` discipline) before any candidate is named.
