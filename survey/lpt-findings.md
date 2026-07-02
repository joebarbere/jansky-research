# Findings — the first LPT population catalogue + P–Ṗ diagram (plan 35)

`jansky_research.lpt` + `data/lpt_sample.csv` (13 confirmed LPTs, every value carrying its arXiv
provenance — compiled 2026-07 from the discovery papers, cross-checked against the Rea+2026 review
arXiv:2601.10393 and its GitHub data file).

## What the compilation itself found

- **A transcription error in the review's own machine-readable table**: ASKAP J1935+2148 period
  2225.309 s vs the discovery paper's 3225.313±0.002 s (dropped leading digit). Flagged; the
  discovery value is used.
- **Only 2 Ṗ measurements exist in the whole class**: CHIME J0630+25 (+5.2e-12, through a glitch —
  the review conservatively downgrades it to a limit; both framings carried) and CHIME/ILT
  J1634+44 (−9.03e-12 — the class's only firm spin-UP, natural for a binary).
- Post-review member included: ASKAP J1745−5051 (accreting WD binary, Rose+2026).

## Population statistics (all regenerate from the CSV)

| quantity | value |
|---|---|
| N (confirmed, mid-2026) | 13 |
| WD binaries / candidates | 6 |
| X-ray detected | 3 |
| period range | 7.0 min – 6.45 hr (median 69.8 min) |
| Ṗ-constrained objects below the pulsar death line | **9 / 9** |
| WD-binary vs rest period offset | Δlog P = 0.294, **p_perm = 0.27 — NOT significant at N=13** |
| same test on a synthetic real split | p ≈ 0.022 (the test has power; the non-detection is informative) |

The 9/9 below-death-line fact is the class's central puzzle made quantitative; the hinted ~78-min
binary boundary is genuinely open, not yet established.

## Honesty rails

- Compilation, not discovery: every number is someone else's measurement with provenance.
- NS dipole quantities (B, τ) NOT assigned to binary members (orbital periods).
- GCRT J1745−3009 flagged as the weakest member; disputed values carry flags in the CSV.
- Selection effects on the period distribution are not modelled (stated).
