# Findings — a verified, provenance-carrying LPT population catalogue (plan 35)

`jansky_research.lpt` + `data/lpt_sample.csv` (every value carrying its arXiv provenance — compiled
2026-07 from the discovery papers, cross-checked against the Rea+2026 review arXiv:2601.10393 and
its GitHub data file). GATE-2 note: the review's own Fig. 3 already plots the class on P–Ṗ — the
contribution here is the provenance-typed table + regenerable statistics, not the first diagram.
Companion of GLEAM-X J0704-37 corrected to M3 (not M5).

**v3 update (plan 44, 2026-07-08):** the catalogue was extended in place from 13 → **16** with the
three 2026 discoveries (ASKAP J142431.2−612611; VASTER's J165130.3−450520 and J170036.6−445758),
each coordinate verified against its source-name convention. The stats below are the v3 (N=16)
values; the period-split test moved to Δlog P = 0.176, **p_perm = 0.52 — still not significant**
(the 2 long-period VASTER sources have unreported binary status, weakening the WD-binary/rest
period contrast). The Stokes-V forced-photometry leg lives in the sibling `lptv` slice.

## What the compilation itself found

- **A transcription error in the review's own machine-readable table**: ASKAP J1935+2148 period
  2225.309 s vs the discovery paper's 3225.313±0.002 s (dropped leading digit). Flagged; the
  discovery value is used.
- **Only 2 Ṗ measurements exist in the whole class**: CHIME J0630+25 (+5.2e-12, through a glitch —
  the review conservatively downgrades it to a limit; both framings carried) and CHIME
  J1634+44 (−9.03e-12 — the class's only firm spin-UP, natural for a binary).
- Post-review member included: ASKAP J1745−5051 (accreting WD binary, Rose+2026).

## Population statistics (all regenerate from the CSV)

| quantity | value |
|---|---|
| N (confirmed, v3 2026-07-08) | 16 |
| WD binaries / candidates | 7 |
| X-ray detected | 3 |
| period range | 7.0 min – 6.48 hr (median 73.4 min) |
| Ṗ-constrained objects below the pulsar death line | **9 / 9** |
| WD-binary vs rest period offset | Δlog P = 0.176, **p_perm = 0.52 — NOT significant at N=16** |
| same test on a synthetic real split | p ≈ 0.022 (the test has power; the non-detection is informative) |

The 9/9 below-death-line fact is the class's central puzzle made quantitative; the hinted ~78-min
binary boundary is genuinely open, not yet established.

## Counterpart cross-match (live, 2026-07-02; `lpt.crossmatch_counterparts`)

Per-object VLASS QL2 cone (20″, epoch 1) + LoTSS DR3 forced cutout peak at each LPT position:

- **VLASS (10/13 in Dec>−40 coverage): no persistent 2–4 GHz counterpart above ~0.7 mJy (5σ QL)
  for ANY object.**
- **LoTSS DR3 (3/13 in returned coverage: ASKAP J1935+2148, ILT J1101+5521, CHIME J1634+44): all
  undetected in the mosaic** (peak/rms mJy: 0.30/0.15, 0.01/0.008, −0.02/0.019). ILT J1101's
  non-detection in the *mosaic* is consistent with its burst-only conjunction emission.
- Remaining objects: outside VLASS dec range (3) / LoTSS returned HTTP 500 = no footprint (10).

Reading: the class is burst-only at current survey depths — persistent emission ≲0.7 mJy at
2–4 GHz and ≲0.04–0.76 mJy (5σ) at 144 MHz where covered. Kept in findings (not the paper) since
these numbers don't flow through the macro pipeline yet.

## Honesty rails

- Compilation, not discovery: every number is someone else's measurement with provenance.
- NS dipole quantities (B, τ) NOT assigned to binary members (orbital periods).
- GCRT J1745−3009 flagged as the weakest member; disputed values carry flags in the CSV.
- Selection effects on the period distribution are not modelled (stated).
