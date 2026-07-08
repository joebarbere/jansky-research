# 44 — LPT catalogue v3 + first multi-epoch Stokes-V forced photometry at all LPT positions

Status: ✅ done — GATE 0 2026-07-08: 3 new rows verified (ASKAP J1424-6126 arXiv:2603.07857;
J1651-4505 + J1700-4457 arXiv:2606.20067), coordinates decoded from source names (<1"); Stokes-V
novelty PASS (RACS-low2 blind V catalogue is a fence, per-source pol on a handful cited). Catalogue
leg: N 13→16; **the ~78-min WD-binary period boundary is still not significant (p=0.52)** — the
plan's headline question, answered. V leg: first systematic multi-epoch forced Stokes-V at all LPT
positions via CASDA (reuses stokesv.measure_circular_pol + the wdpulsar CASDA driver; 191 rows,
15/16 covered). Not all-limits: **1 secure single-epoch circular detection** (ASKAP J1745-5051,
the accreting CV, 15% circular at 21.6σ, on-centre) + **1 candidate** (VASTER ASKAP J1651-4505,
59% circular at 12.5σ but 3.2″ off) + 1 confusion peak vetoed (240 mJy/5.3″); median 3σ V limit
0.474 mJy. Both are single-epoch burst states, not persistent counterparts. GATE-2 PASS w/ fixes
(J1651 downgraded to candidate; confusion + secure/candidate vetoes disclosed as heuristics). See
survey/lptv-findings.md.

## Context

Two coupled gaps. (a) `data/lpt_sample.csv` (13 objects, merged `lpt` slice) predates ≥3 2026
discoveries: ASKAP J142431.2−612611 (arXiv:2603.07857), ASKAP J165130.3−450520 and
J170036.6−445758 (VASTER, arXiv:2606.20067) — a ~30% undercount — and the Rea+2026 review's "no
population synthesis" flag still stands. (b) The merged counterpart cross-match was Stokes-I
only (VLASS/LoTSS); nobody has done forced V photometry at the LPT positions across RACS
low1/low2/mid — RACS-low2 Paper VIII's blind V catalogue (arXiv:2606.16182) did not target them.
A persistent V counterpart discriminates coherent-emission models; systematic V limits sharpen
the burst-only character. Likely all-limits outcome at RACS depth (low duty cycles) — the limit
table is the paper. Near-total reuse of two merged slices (`lpt`, `stokesv`); data via
discovery-paper tables + CASDA RACS V epochs at ~16–17 positions. CPU + CASDA I/O, days.

## Deliverables

- `src/jansky_research/lpt.py` (extend v2 in place) + `src/jansky_research/lptv.py`:
  `add_2026_rows` (schema-checked, provenance-flagged), re-run `crossmatch_counterparts` and the
  population statistics at N=16–17 (`binary_boundary_test` — does the ~78-min p move?);
  `fetch_racs_v_epochs` (# pragma), `forced_v_photometry` (per epoch per position via
  `measure_circular_pol` + leakage vetting), `handedness_changes` (inter-epoch V sign),
  `v_limits_table`, `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/lptv/`; `survey/lptv-findings.md`; wiring.

## Approach

0. GATE 0: full-text pass on arXiv:2603.07857, 2606.20067, 2606.16182, and the Rea+2026 review
   to confirm the new rows, verify no targeted LPT V photometry has appeared, and spot-check
   coordinates/periods for the v3 rows; re-verify CASDA V-epoch availability at each position.
1. Catalogue leg: add the 2026 rows with the `lpt` flag discipline (the process that caught the
   review's own data-file typo); re-run the counterpart cross-match + population statistics and
   report whether N=16–17 moves the ~78-min binary-boundary test's p-value.
2. Stokes-V tooling + recover-a-known: SRSC positive controls (176/176) + GJ 65 re-detection
   before touching any LPT position (the `stokesv_discovery` control pattern).
3. Real leg: `measure_circular_pol` per epoch per position across RACS low1/low2/mid with
   leakage vetting; report V detections/limits + any inter-epoch handedness/sign changes.
4. GATE-2 science review: leakage floor vs faint-V claims, duty-cycle blindness of snapshot
   epochs, small-N caveats on the population statistics. 5. Paper: v3 catalogue + V table.

## Verification

SRSC 176/176 positive controls and GJ 65 re-detected before LPT positions are measured; new CSV
rows spot-checked against the discovery papers; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Likely all-limits outcome at RACS depth (low duty cycles)** → the V limit table is the
  paper; framed as the deliverable, not spin.
- **Transcription errors in new rows** → keep the `lpt` flag/provenance discipline; every value
  traced to a discovery-paper table.
- **Discovery teams may publish targeted polarimetry** → GATE-0 full-text check; the uniform
  multi-epoch survey V table remains complementary to single-object campaigns.
