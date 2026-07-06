# 46 — Voyager 2 PRA: modern re-derivation of the Uranus & Neptune radio rotation periods

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — pin the PDS-PPI
VG2-PRA Uranus/Neptune encounter dataset IDs before any code

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
