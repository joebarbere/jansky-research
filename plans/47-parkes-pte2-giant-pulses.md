# 47 — Parkes Transient Events II: giant-pulse/heavy-tail census across 363 pulsars

Status: ✅ done 2026-07-10 — an HONEST NULL. Uniform per-source giant-pulse (heavy-tail) test over
all 363 PTE-II pulsars, 136 fitted (>=50 pulses): 26 (19%) flagged heavy-tailed, but (i) NO
correlation with Edot (Spearman -0.03, p=0.76; MW p=0.13), (ii) the classification is
**detection-power limited** (`count_limited`: heavy sources 2.4x more pulses, MW p=0.02; flag rate
0.088 low-count -> 0.294 high-count half), (iii) the tails are steep (median index 11.6, not the
~2-3 of true giant pulses), and the classic GP archetypes (Crab/B1937+21/Vela) don't survive the
count cut (only B0950+08 present, correctly flagged). Synthetic recover-a-known confirms the test
works (completeness 0.33->1.0 vs count, FP ~0.05) so the null is a real data limit, not tooling.
No giant-pulse population and no Edot trend claimed; confound is pipeline-generated. Module
`pte2.py` (+`pte2_real.py`), paper `papers/pte2/`, findings `survey/pte2-findings.md`. — GATE 0
done 2026-07-10: **CONDITIONAL PASS, novelty NARROWED**. Paper =
Yang+2025, *Parkes transient events II* (arXiv:2508.14403, **ApJS** 10.3847/1538-4365/adfe67):
165,592 single pulses / 363 pulsars / Parkes MB 1997-2001; SQLite3 (~1.5 GB zip on GitHub LFS
`Astroyx/Pulsar_collection` -> `Pulsar_fits_database_v1.zip`, sha256 775dffa0...; mirror CSIRO DAP
10.25919/34am-zx04; **open, no auth -> real leg runnable**). **Correction: the paper is NOT
infrastructure-only** -- it already fits log-normal/Gaussian fluence distributions for 98 pulsars
(>100 detections) + a population power-law (alpha=-1.7+/-0.1) + a 25.3% giant-pulse-fraction stat. And
**HTRU-V (Burke-Spolaor+2012, MNRAS 423, 1351)** is a prior 315-pulsar log-normal single-pulse
energy census. So a generic "363-pulsar energy-distribution census" is **anticipated and NOT novel**.
**The defensible wedge (this slice)**: uniform PER-SOURCE heavy-tail **model selection** -- log-normal
vs power-law-**tail** via a principled likelihood ratio (Vuong/Clauset) -- across ALL 363, correlated
with Edot; HTRU-V tested log-normal/Gaussian (not heavy-tail) and PTE-II fit only 98 (population
power-law, not per-source tails). **Schema**: `pulsar`(jname, p0, s1400, spNumber, smax/smin) <-
`file`(pfLinkID, timeStartMJD, gain_factor) <- `fileSegment`/`seg_file`(pfLinkID, snr_max/snr_min,
data BLOB). Per-pulse FLUENCE is NOT stored -- derived on-the-fly from blobs by the paper's tool; we
use per-pulse **S/N** (snr_max) as the energy proxy (proportional to fluence within a source; standard
for tail-shape work). Join to ATNF Edot by `jname`=PSRJ. Recover-a-known heavy-tail set: B1937+21,
Vela (J0835-4510), + PTE-II giant-pulse-fraction pulsars. Reuse: `frbstats.fit_power_law`/`select_xmin`,
`ppdot.spindown_luminosity` + `fetch_atnf_ppdot`, `stacking.injection_recovery` discipline.

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
