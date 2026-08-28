# Plan 95: paper titles enter traditional-style scope, and the repo is remediated

## Why

The lptduty note's title ("...Are Caught in a Few Per Cent of ASKAP Snapshots, and the
Published Ephemerides Rarely Say Why") survived the entire style campaign because `\title{}`
was never in scope: not fingerprinted, not linted, and the batch instructions explicitly
carved headings out. Joe's call (2026-08-28): titles are in scope, for conversion and for
new papers.

## The empirical rule (measured 2026-08-28, corpus metadata)

- 113,829 refereed journal titles: **0.9%** sentence-like (case-insensitive finite-verb
  match, "?"-titles exempt — question titles are corpus-attested, including the
  question-plus-subtitle form); p50 11 words, p90 17.
- 1,035 RNAAS titles: **6.0%** sentence-like — the "verdict title" is genuinely attested
  in that genre, but noun phrases dominate 15:1.
- Rule (guide §S4b): a title is a noun phrase ≤~18 words; RNAAS may carry a short
  declarative verdict; editorializing clauses have no analogue anywhere.

## The machinery

`stylecorpus.latex_paper_title` + `title_is_sentence_like` + two fingerprint metrics
(`title_words`, `title_sentence_like`); `lint_paper` gains genre-aware title rules
(sentence-shaped: MED for papers, LOW for RNAAS; length > 18: LOW), thresholds embedded as
measured constants since the committed corpus percentiles predate the metrics. SKILL.md,
style-guide §S4b, and the style-editor agent all updated; the agent proposes titles rather
than silently rewording (a title is the author's most visible wording).

## The remediation (audit of all 46 papers + notes: 21 offenders)

Retitled (19), each verified to lose no claim absent from its abstract:
sentence-shaped — junodam, spectra (main), stokesv, typeii, vgpra (main; the note keeps its
RNAAS-attested verdict form but is shortened); length — driftsearch, fashienv, innerrc,
offsets, pte2, solarbursts, stacking, svsbi, triangulate, type3synthesis, vlass, vlbi,
wdpulsar (main).

Deliberate exceptions (2):
- **dr20radio** (19 words, noun phrase): its "First Southern-Hemisphere SDSS Quasar
  Spectra" priority claim appears NOWHERE in the abstract — trimming the title would delete
  the claim from the paper's front matter entirely. Kept; LOW accepted.
- **spectra rnaas** (10 words, verdict sentence): the RNAAS verdict form is corpus-attested
  (6.0%), and "Is Neither Pure Nor Complete" is the note's one result. Kept; LOW accepted.

lptduty's note was retitled separately (PR #301) and is the pattern case recorded in §S4b:
its sentence title had also fossilized a superseded framing — a noun-phrase title naming
the quantity survives revision; a sentence title states a conclusion, and conclusions rot.
