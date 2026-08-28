# Plan 94: driftsearch RNAAS note — the replication note nobody has written

## Gate 0 (novelty, checked 2026-08-27)

ADS: RNAAS contains four notes mentioning Voyager (flybys, quasi-thermal noise — none about
carrier detection) and ~30 SETI/technosignature/Breakthrough-Listen notes, none of which is an
independent replication of the Voyager-1 validation or a drift-search benchmark. RNAAS is BL's
own house venue for this scale (Enriquez 2018 'Oumuamua; Perez 2022 WOW!; Sheikh 2021 blc1
non-redetection; Price 2021 data recorder — an instrumentation/software note; Jacobson-Bell
2025 3I/ATLAS), so software-and-validation content is demonstrably in scope.

## The deliberate venue-policy exception

The repo policy sends recover-a-knowns to repo+Zenodo. The exception is justified by what is
NOT in the literature: (1) a citable, fully-seeded injection–recovery benchmark cell any
teaching-grade detector can regression against; (2) the DC-spike trap measured in the exact
public file every newcomer downloads (brightest channel = artifact at ~10^3× the carrier);
(3) the locate-don't-assert failure mode with its number (asserted carrier 0.92 MHz off →
false null at S/N ~4.9 where the located carrier sits at S/N 997). Frame as an independent
minimal-code replication of the standard validation (Enriquez et al. 2017; Lebofsky et al.
2019 data), led by the benchmark and the pitfalls. No discovery claimed.

## Work

1. `papers/driftsearch/rnaas.tex` via the condensation workflow: ≤~1000 countable words, every
   number from `generated/macros.tex` (no hand-typed values), one figure
   (`figures/drift_recovery.pdf`), refs subset of the existing verified `refs.bib`.
2. Gates: `make paper SLICE=driftsearch` (builds any `*.tex` with a documentclass),
   `triage_papers.py --paper driftsearch` (note gets its own row, zero HIGH/MED),
   `prose_lint --file rnaas.tex --genre rnaas` (+ diff-guard N/A for a new file), word count.
3. Record in `survey/drift-findings.md`; README short-notes list gains the fourth note;
   CHANGELOG entry; PR.

## Open gate before submission

The condensation referee round (the wdpulsar lesson: shortening introduces errors — its note
inflated an efficiency while condensing). This session's subagent budget is exhausted, so the
note ships to the repo with an inline self-check only; a `paper-referee` round on the note is
REQUIRED in a future session before it goes to journals.aas.org. Tracked here and in findings.
