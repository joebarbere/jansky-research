# 46 — Voyager 2 PRA: modern re-derivation of the Uranus & Neptune radio rotation periods

Status: ✅ done 2026-07-10 — a CONTROLLED NULL: blind Lomb-Scargle of the PRA total-power flux
recovers a clean injected rotation in synthetic tests (same wide 14–20 h window) but **recovers
NEITHER real ice-giant period** — Uranus peak 18.44 h (wanders across sub-bands, +1.2 h off 17.24 h),
Neptune rails to the 20 h search bound (+3.9 h off 16.11 h). The failure is a data limitation (the
auroral total-power flux isn't a clean rotational sinusoid over the short flyby), not a pipeline bug;
the historical beaming/magnetic-longitude modelling was essential, and the ~2 h flyby precision is
hundreds× coarser than the 28-s HST shift regardless. GATE-2 PASS (caught: the earlier narrow-window
"Neptune 16.04≈16.11" was a window artifact; band-stability is a coherence check, NOT a right/wrong
gate — Uranus's wrong peak is more band-stable than the truth; "independent" → same dataset,
re-analysis; citation fixes). Module `vgpra.py`, driver `scripts/vgpra_real.py`, paper
`papers/vgpra/`, findings `survey/vgpra-findings.md`. — GATE 0 done
2026-07-10: novelty PASS + data pinned. **Dataset IDs**:
`VG2-U-PRA-3-RDR-LOWBAND-6SEC-V1.0` (`DATA/VG2_URN_PRA_6SEC.TAB`, 49 MB) and
`VG2-N-PRA-3-RDR-LOWBAND-6SEC-V1.0` (`DATA/VG2_NEP_PRA_6SEC.TAB`, 79 MB) on PDS-PPI
(`pds-ppi.igpp.ucla.edu/data/<id>/DATA/`), direct HTTP, egress verified (curl 200).
**Format** (pinned from the .LBL + real bytes): fixed-width ASCII, line-based (2284 data bytes +
LF; NOT the 2286 RECORD_BYTES stride). Each row = one 48-s major frame: `DATE` (YYMMDD, chars
0:6) + `SECOND` (sec-of-day, 6:12) + 8 sweeps × (1 status word + 70 channels), each field 4 chars
ASCII int in **MILLIBELL** (=0.01 dB). Sweep k (0..7) starts at `SECOND + 6*k` s. **70 low-band
channels**: `f_i = 1326.0 − 19.2*i` kHz, i=0..69 (1326.0 → 1.2 kHz). Neptune 36028 rows ≈ 20 d ≈
30 rotations; Uranus ≈ 12 d ≈ 17 rotations (tens of cycles, not just a few — period recoverable,
but 28-s-level shifts are below the achievable precision → the honest deliverable). Real sample
(first 8 Neptune rows) vendored at `tests/data/vg2_nep_pra_6sec_sample.tab` for the parser test.
**Novelty PASS**: Lamy+2025 (Nat.Astron. 9, 658; 17.247864±0.000010 h) is HST-UV aurora, NOT
radio; the recent arXiv:2604.19863 is radio *occultations* (atmospheric geometry), not a PRA
burst re-derivation; Cecconi+2017 (arXiv:1710.10471) refurbished only Jupiter/Saturn. No modern
radio re-analysis of the Uranus/Neptune PRA rotation periods exists. **Published anchors**: Uranus
17.24±0.01 h (Warwick+1986), Neptune 16.11±0.05 h (Warwick+1989); 16.108±0.006 h (Lecacheux+1993).
Reuse: `frbperiod` (`rayleigh_z2`, `period_search`, scramble FAPs), `report._agg/_fmt_p`.

## Context

Lamy+2025 (Nat. Astron. 9, 658) moved Uranus's rotation period by 28 s using 11 yr of HST UV
aurora — yet the 1986 *radio* value (17.24±0.01 h) that underlies System III was never reanalysed
with modern statistics, and Neptune's 16.11 h has had no independent check in 28 yr. The
Cecconi+2017 PRA refurbishing (arXiv:1710.10471) covers only Jupiter/Saturn, so the Uranus and
Neptune encounter volumes are an unworked niche. Data: the PDS-PPI VG2-PRA Uranus/Neptune
encounter volumes (small; exact dataset IDs to be pinned at GATE 0). This slice reuses the merged
`frbperiod` machinery (Lomb-Scargle + Rayleigh Z², scramble-based FAPs) on planetary radio burst
time series, and delivers a three-way comparison: the 1986 radio value, Lamy+2025, this work.
Honest framing: with days-long flybys the posteriors may be wide — "modern uncertainties say the
radio data cannot distinguish" is itself the citable result.

## Deliverables

- `src/jansky_research/vgpra.py`: `fetch_pra_volumes` (PDS-PPI download, `# pragma`),
  `read_pra_series` (encounter-volume parser → burst/flux time series), `detect_bursts`
  (background + kσ per channel/bin), `period_posterior` (LS/Rayleigh Z² via `frbperiod` reuse,
  few-cycle-honest uncertainties from scrambles + bootstrap), `synthetic_flyby` (injected period
  into a flyby-length noise series → recovery), `compare_periods` (three-way table),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/vgpra/`; `survey/vgpra-findings.md`; wiring.

## Approach

0. **GATE 0:** pin the PDS-PPI VG2-PRA Uranus and Neptune encounter dataset IDs (volume names,
   format docs, sizes) and confirm direct download; full-text pass on Lamy+2025, Cecconi+2017
   (arXiv:1710.10471), and the original 1986/1989 radio-period papers to confirm no modern radio
   reanalysis has appeared.
1. Tooling + synthetic recover-a-known: inject a known period into a flyby-length synthetic
   series; posterior must recover it with honestly wide few-cycle error bars.
2. Real leg, Uranus: parse the encounter volume, extract the burst/flux series, run the period
   posterior; repeat for Neptune. Small data, CPU, hours.
3. GATE-2 science review: few-cycle uncertainty honesty (no overclaiming precision the flyby
   cannot support), burst-selection sensitivity, three-way-comparison framing.
4. Paper: modern posteriors for both periods + the comparison table.

## Verification

Synthetic flyby round-trip recovers the injected period; consistency with ~17.24 h (Uranus) and
~16.11 h (Neptune) is itself the pipeline validation before any discrepancy claim; checks green;
GATE-2 sign-off.

## Risks & mitigations

- **Days-long flybys → wide posteriors** — the honest paper may be "the radio data cannot
  distinguish the 1986 value from Lamy+2025"; frame that as the deliverable from day one, still
  citable.
- **Encounter-volume format archaeology** (1980s PDS layouts) → pin format docs at GATE 0; keep
  the reader tested against a vendored sample block.
