# 47 — Parkes Transient Events II: giant-pulse/heavy-tail census across 363 pulsars

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — confirm the PTE-II
sqlite DB download path and that raw segments are included

## Context

PTE-II (arXiv:2508.14403) is a brand-new ~1.5 GB sqlite database of 165,592 single pulses from
363 pulsars (Parkes 1997–2001 archival reprocessing, raw segments preserved), and its own paper
is infrastructure-only — no science analysis. Only a handful of pulsars have published
single-pulse energy-tail fits, so a uniform giant-pulse/heavy-tail census across all 363 is an
open gap. Method: per-pulsar energy histograms with lognormal vs power-law-tail model comparison,
tail index ranked against Ė from the ATNF catalogue, exposure-normalized via injection-recovery
(the merged `stacking` slice's discipline). The merged `fdmt`/`singlepulse` tooling handles any
re-extraction from the preserved raw segments, and the GPU makes that cheap. Tiny footprint —
well inside the 275 GB disk budget.

## Deliverables

- `src/jansky_research/pte2.py`: `fetch_pte2_db` (download + SHA check, `# pragma`),
  `load_pulse_energies` (sqlite → per-pulsar energy tables), `fit_energy_tail` (lognormal vs
  power-law-tail, BIC/likelihood-ratio comparison), `tail_vs_edot` (ATNF cross-match + rank
  statistics), `inject_recover_tail` (injected tail slope into synthetic pulse sets → recovery,
  exposure normalization), optional `reextract_segment` (`fdmt`/`singlepulse` reuse on raw
  segments, GPU), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/pte2/`; `survey/pte2-findings.md`; wiring.

## Approach

0. **GATE 0:** verify the PTE-II data release URL (sqlite DB, ~1.5 GB) and schema; full-text
   pass on arXiv:2508.14403 and its citing papers to confirm no energy-tail census has appeared
   since the scan; confirm which pulsars have published tail fits (the recover-a-known set).
1. Tooling + synthetic recover-a-known: inject lognormal and power-law-tailed pulse-energy
   populations at PTE-II-like counts and sensitivity; the model comparison must classify each
   correctly and recover the injected tail index.
2. Real leg: per-pulsar histograms and fits across all 363; ATNF Ė cross-match; heavy-tail
   ranking; spot re-extraction of a few segments as a consistency check (GPU).
3. GATE-2 science review: heterogeneous 1997–2001 sensitivity (normalize or state lower bounds),
   small-N pulsars excluded by a stated count cut, selection effects of the original PTE search.
4. Paper: the census table + tail-index-vs-Ė result (or its honest absence).

## Verification

Synthetic tails round-trip at stated confidence; Vela and B1937+21-class known heavy-tail
emitters must come out heavy-tailed before any new classification is trusted; checks green;
GATE-2 sign-off.

## Risks & mitigations

- **Heterogeneous 1997–2001 sensitivity** across observing setups → injection-recovery
  normalization per epoch/backend where metadata allows; otherwise report tail fractions as
  lower bounds, plainly labelled.
- **Low-count pulsars** → pre-registered minimum-pulse-count cut; report only above it.
