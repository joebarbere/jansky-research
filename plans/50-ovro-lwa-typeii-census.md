# 50 — OVRO-LWA metric type II census × LASCO CME catalogue

Status: ✅ done — an HONEST NULL. GATE 0 2026-07-09: novelty PASS (OVRO-LWA detector arXiv:2603.25446
is type-III-only). Built the slow-drift+harmonic detector (synthetic SNR-completeness curve
1.0→0.33), the in-memory streaming pipeline (data on AWS Open Data, `ovro-lwa-solar` S3; ~30 s/day,
nothing on disk), and the full cross-match (CDAW CMEs + HEK GOES flares + SILSO occurrence).
**Ran the full census**: all 765 observing days 2024-04→2026-07, 0 failures → 331 candidates that
are **false-positive dominated** (matched-CME median 478≈background 379 km/s; observed match rate
0.55 < chance 0.64; drift⊥CME-speed r=0.09; 83% window-saturated; harmonic cut worsens it; the
flare-gated "signal" is a flare↔fast-CME confound). A blind spectral type II census fails in this
RFI-heavy band — why the archive detector is type-III-only. GATE-2 PASS w/ fixes (verified the null
is correct + not a missed sub-population; required the confound + window-saturation evidence be
pipeline-generated, now in `purity_diagnostics`). No census/rate/detection claimed. See
survey/typeii-findings.md.

## Context

Every type-II-adjacent census is cycle-24, ascending-phase-only, or N=10 case studies
(arXiv:2512.21846; the RSTN study, Solar Phys. 2024), and the published OVRO-LWA real-time
detector (ApJ 1003, 57, arXiv:2603.25446) is **type-III/IIIb-only** — that fence from the
Corrections section is simultaneously the gap evidence: the detector lineage skips type II.
The repo's own TODO already wants a slower-drift template. Data: OVRO-LWA Level-1 spectrograms,
2024-04→now, via the `ovsa.njit.edu/lwadata-query` portal (GATE 0), and/or e-Callisto as a
cross-check stream. Method: a slow-drift ridge + band-split heuristic type II detector,
cross-matched against the LASCO CME catalogue v2 and GOES flare lists, with occurrence vs cycle
phase from SILSO. Reuses the `solarbursts` detector pattern and the `ecallisto_census`
coverage-correction machinery. Sibling on the same data (type I noise-storm census) is noted but
explicitly *not* bundled — one first.

## Deliverables

- `src/jansky_research/typeii.py`: `fetch_lwa_spectrograms` (portal query + download,
  `# pragma`), `detect_typeii` (slow-drift ridge tracker + fundamental/harmonic band-split
  heuristic), `coverage_correction` (`ecallisto_census` reuse — per-day observing fraction),
  `crossmatch_cme_flare` (LASCO CME v2 + GOES windows), `occurrence_vs_phase` (SILSO cycle
  phase), `synthetic_typeii` (injected slow-drift band-split ridges + type III contaminants →
  detector completeness/purity), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/typeii/`; `survey/typeii-findings.md`; wiring.

## Approach

0. **GATE 0:** verify `ovsa.njit.edu/lwadata-query` access, file format, cadence, and per-day
   volume for Level-1 spectrograms; full-text pass on arXiv:2603.25446 (confirm type-III-only),
   arXiv:2512.21846, and the RSTN paper to confirm no type II census on this archive exists;
   confirm LASCO CME v2 + GOES + SILSO access.
1. Tooling + synthetic recover-a-known: injected type II ridges (with band-split) among type III
   contaminants and RFI; detector must separate them at stated completeness/purity.
2. Real leg: run the detector over 2024-04→now (CPU/GPU, batched); coverage-corrected event
   list; LASCO/GOES cross-match; occurrence vs cycle phase.
3. GATE-2 science review: the ~2 yr baseline caveat (no strong cycle-phase claims), 13–40 MHz
   RFI contamination handling, detector-threshold selection effects, e-Callisto cross-check.
4. Paper: the census + CME-association statistics.

## Verification

The census must reproduce the established type II–fast-CME association fraction before any new
occurrence claim; synthetic completeness/purity at stated thresholds; checks green; GATE-2
sign-off.

## Risks & mitigations

- **~2 yr baseline limits cycle-phase claims** → frame as a max-phase census with the phase
  trend as indicative only; the coverage-corrected event list is the durable product.
- **RFI at 13–40 MHz** → RFI masking before ridge tracking; e-Callisto cross-check on a burst
  subsample; report the masked-band fraction.
- **OVRO/NJIT team extends their detector to type II** → they are the natural authors; ship
  bounded and fast, and keep the e-Callisto-only fallback in scope.
